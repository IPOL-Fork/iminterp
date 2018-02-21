/**
 * @file sset.c
 * @brief Contour stencil set construction
 * @author Pascal Getreuer <getreuer@gmail.com>
 *
 *
 * Copyright (c) 2010-2011, Pascal Getreuer
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD License. You
 * should have received a copy of this license along this program. If
 * not, see <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <stdio.h>
#include "sset.h"
#include <ipol/imageio.h>

/** @brief The angular resolution for stencil quantization */
#define NUMANGLES               64
/** @brief Number of quadrature panels to use for integrals */
#define QUAD_RES                8
/** @brief Initial capacity of a stencil set */
#define INITIAL_CAPACITY        16

/** @brief Stability threshold */
#define STABILITY_THRESH        1e-4


/**
 * @brief A single contour stencil \f$\mathcal{S}\f$
 *
 * This data structure represents a single contour stencil \f$\mathcal{S}\f$.
 * It comprises a list of vectors (Alpha[k],Beta[k]) =
 *     \f$(\alpha_k, \beta_k),\f$
 * indicating the value of \f$\mathcal{S}\f$ in the kth cell \f$\Omega_k\f$, a
 * function pointer Phi (and supporting parameters) for the distance function
 * \f$\varphi\f$ from which the stencil was discretized, and "glyph" data used
 * by DrawContours() for a graphical representation of the stencil.
 *
 * The cell locations are stored in ssetstruct::Cell.
 */
typedef struct
{
    /** @brief Alpha[k] is the x-component of the stencil in the kth cell */
    double *Alpha;
    /** @brief Beta[k] is the y-component of the stencil in the kth cell */
    double *Beta;
    
    /**
     * @brief Function pointer to the stencil's distance function \f$\varphi\f$
     *
     * Phi is a function with the syntax Phi(x,y,PhiParam), where
     * stencilentry::PhiParam may be used to pass additional parameters to Phi.
     * The function is rotated by stencilentry::PhiTrans.
     */
    double (*Phi)(double, double, const void*);
    /**
     * @brief Cosine and sine values of the stencil rotation
     *
     * The values PhiTrans[0] = cos(Rotation), PhiTrans[1] = sin(Rotation) are
     * computed once and reused for all rotation calculations.  The rotation is
     * performed as
     * \f[ \begin{pmatrix} x' \\ y' \end{pmatrix}
     *     = \begin{pmatrix} \mathtt{PhiTrans}[0] & \mathtt{PhiTrans}[1] \\
     *     -\mathtt{PhiTrans}[1] & \mathtt{PhiTrans}[0] \end{pmatrix}
     *     \begin{pmatrix} x \\ y \end{pmatrix} \f]
     */
    double PhiTrans[2];
    /**
     * @brief User-specified parameters to pass to stencilentry::Phi
     *
     * PhiParam is passed as the third argument to Phi.  PhiParam should be
     * used to point to a data structure containing whatever additional
     * parameters are needed for evaluating Phi.  For example, CurveDist()
     * uses PhiParam to specify the curviness parameter "a" of the parabola.
     */
    const void *PhiParam;
    
    /**
     * @brief Function to draw the stencil
     *
     * This function should draw a representation of the stencil using the
     * graphical commands in pen.h.  The given pen object is transformed such
     * that the origin is the center of the stencil and rotated according to
     * stencilentry::PhiTrans.
     */
    int (*DrawFun)(pen *Pen, const void *Param);
    /** @brief Parameter for DrawFun */
    const void *DrawParam;
    /**
     * @brief Line drawing color
     *
     * (DrawColor[0], DrawColor[1], DrawColor[2]) is the RGB color to use
     * for drawing the stencil.
     */
    const float *DrawColor;
} stencilentry;


/** @brief A cell in a contour stencil, \f$\Omega_k \subset \mathbb{R}^2\f$
 *
 * The cell is a square with unit side lengths and (x,y) specifying the
 * bottom-left corner.  The cell is
 *    \f[\Omega_k = [x,x+1] \times [y,y+1].\f]
 * The stencilsetcell type is used by ssetstruct::Cell to represent all the
 * cells in the stencil.
 */
typedef struct
{
    /** @brief x-component of the left edge of the cell */
    int x;
    /** @brief y-component of the bottom edge of the cell */
    int y;
} stencilsetcell;


/**
 * @brief Quantized vector
 *
 * For computational efficiency in FitStencils(), the stencil vectors
 * (stencilentry::Alpha[k], stencilentry::Beta[k]) are converted to polar form
 * and quantized to #NUMANGLES possible angles.  stencilqvec::AngleIndex is a
 * number between 0 and #NUMANGLES-1, and stencilqvec::Weight is the magnitude
 * of the vector.
 */
typedef struct
{
    /** @brief Quantized angle index, a number between 0 and #NUMANGLES-1 */
    int32_t AngleIndex;
    /** @brief Vector magnitude */
    float Weight;
} stencilqvec;


/**
 * @brief A contour stencil set \f$\Sigma\f$
 *
 * A stencil set is a collection of contour stencils.  This is represented as
 * an array ssetstruct::Stencil of type #stencilentry.  For computational
 * efficiency, an array of quantized stencil vectors #stencilqvec is also
 * stored in ssetstruct::StencilTable.  The locations of the cells are stored
 * in stencilsetcell::Cell.  The point spread function (PSF) to model the
 * sampling of the input image is specified by ssetstruct::Psf.
 *
 * #ssetstruct is the definition of #sset, which is encapsulated by forward
 * declaration in sset.h.
 *
 * Stencil sets are managed with the following functions:
 *    - NewStencilSet():   create a new sset
 *    - FreeStencilSet():  free an sset
 *    - AddStencil():      add a stencil to an sset
 *    - GetNumStencils():  get the number of stencils in the set
 *    - GetPsf():          get the PSF function pointer
 *    - GetPsfParam():     get the PSF parameters PsfParam
 *    - GetPsfRadius():    get the PSF radius
 */
typedef struct ssetstruct
{
    /** @brief Number of stencils in the set */
    int NumStencils;
    /** @brief Allocated capacity of Stencil and StencilTable */
    int Capacity;
    /** @brief Array of stencils \f$(\mathcal{S}_0,\mathcal{S}_1,\ldots)\f$ */
    stencilentry *Stencil;
    /** @brief Quantized stencil vectors for FitStencils() */
    stencilqvec *StencilTable;
    
    /** @brief Number of cells in the stencil */
    int NumCells;
    /** @brief The cells, \f$\Omega = \bigcup(\Omega_0,\Omega_1,\ldots)\f$ */
    stencilsetcell *Cell;
    /** @brief Quadrature weights for computing integrals over the cell */
    double *QuadDx;
    /** @brief Quadrature weights for computing integrals over the cell */
    double *QuadDy;
    
    /** @brief Point spread function (PSF) used to sample the input image */
    double (*Psf)(double, double, const void*);
    /** @brief User-defined parameters to pass to Psf */
    const void *PsfParam;
    /** @brief Support radius of the PSF */
    double PsfRadius;
} ssettype;


/**
 * @brief Allocates more memory for the stencil set
 * @param Set the stencil set
 * @return 1 on success, 0 on failure
 */
static int ReallocStencilSet(sset *Set)
{
    stencilentry *Stencil;
    int S, NewCapacity;
    
    
    if(!Set)
        return 0;
    
    /* Either allocate INITIAL_CAPACITY or increase by 10% */
    if(Set->Stencil == NULL)
        NewCapacity = INITIAL_CAPACITY;
    else
        NewCapacity = Set->Capacity + (int)(Set->Capacity/10) + 1;
    
    /* Reallocate the Stencil and StencilTable arrays */
    if(!(Stencil = Set->Stencil = (stencilentry *)
        Realloc(Set->Stencil, sizeof(stencilentry)*NewCapacity))
        || !(Set->StencilTable = (stencilqvec *)
        Realloc(Set->StencilTable, sizeof(stencilqvec)*
        NewCapacity*Set->NumCells)))
        return 0;
    
    /* Set the newly allocated pointers to NULL */
    for(S = Set->Capacity; S < NewCapacity; S++)
    {
        Stencil[S].Phi = NULL;
        Stencil[S].Alpha = NULL;
        Stencil[S].Beta = NULL;
        Stencil[S].DrawFun = NULL;
    }
    
    /* Allocate Alpha and Beta arrays for each new stencil */
    for(S = Set->Capacity; S < NewCapacity; S++)
        if(!(Stencil[S].Alpha = (double *)
            Malloc(sizeof(double)*Set->NumCells))
            || !(Stencil[S].Beta = (double *)
            Malloc(sizeof(double)*Set->NumCells)))
            return 0;
    
    Set->Capacity = NewCapacity;
    return 1;
}


/**
 * @brief Free a stencil set
 * @param Set the stencil set to free
 *
 * Frees all memory associated with a stencil set.
 */
void FreeStencilSet(sset *Set)
{
    if(Set)
    {
        if(Set->Stencil)
        {
            stencilentry *Stencil = Set->Stencil;
            int S;
            
            for(S = 0; S < Set->Capacity; S++)
            {
                Free(Stencil[S].Beta);
                Free(Stencil[S].Alpha);
            }
            
            Free(Stencil);
        }
        
        Free(Set->StencilTable);
        Free(Set->Cell);
        Free(Set->QuadDy);
        Free(Set->QuadDx);
        Free(Set);
    }
}


/**
 * @brief Compute samples of the PSF
 * @param Set the stencil set
 * @return a pointer to an array of PSF samples, NULL on failure
 */
static double *SamplePsf(sset *Set)
{
    const int R = (int)ceil(QUAD_RES * Set->PsfRadius);
    const int W = 2*R + 1;
    double *Samples = NULL;
    int i, j, Offset;
    
    
    if(!(Samples = (double *)Malloc(sizeof(double)*W*W)))
        return NULL;
    
    if(!Set->Psf)
        Samples[0] = QUAD_RES*QUAD_RES;
    else
        for(j = -R, Offset = 0; j <= R; j++)
            for(i = -R; i <= R; i++, Offset++)
                Samples[Offset] = Set->Psf(i/((double)QUAD_RES),
                    j/((double)QUAD_RES), Set->PsfParam);
    
    return Samples;
}


/**
 * @brief Compute quadrature weights for integrating over cells
 * @param Set the stencil set
 * @return 1 on success, 0 on failure
 *
 * Quadrature weights are computed for evaluating integrals of the form
 * \f[ \int_{\Omega_k} \nabla (h * \varphi)^\perp \, dx. \f]
 * For efficiency, these weights are stored in ssetstruct::QuadDx and
 * ssetstruct::QuadDy so that they are computed only once and then reused.
 */
static int ComputeQuadWeights(sset *Set)
{
    const int R = (int)ceil(QUAD_RES * Set->PsfRadius);
    const int W = 2*R + 1;
    double *QuadDx = NULL, *QuadDy = NULL, *PsfSamples = NULL;
    double Sum;
    int x, y, i, j, Offset;
    
    
    if(!(QuadDx = (double *)Malloc(sizeof(double)*
        (QUAD_RES + W)*(QUAD_RES + W)))
        || !(QuadDy = (double *)Malloc(sizeof(double)*
        (QUAD_RES + W)*(QUAD_RES + W)))
        || !(PsfSamples = SamplePsf(Set)))
    {
        Free(PsfSamples);
        Free(QuadDy);
        Free(QuadDx);
        return 0;
    }
    
    /* Compute quadrature weights that incorporate convolution with the PSF */
    for(y = -R, Offset = 0; y <= QUAD_RES + R; y++)
        for(x = -R; x <= QUAD_RES + R; x++, Offset++)
        {
            for(j = -y, Sum = 0; j <= QUAD_RES - y; j++)
            {
                if(abs(j) <= R)
                {
                    i = QUAD_RES - x;
                    
                    if(abs(i) <= R)
                        Sum += PsfSamples[(i + R) + W*(j + R)];
                    
                    i = -x;
                    
                    if(abs(i) <= R)
                        Sum -= PsfSamples[(i + R) + W*(j + R)];
                }
            }
            
            QuadDx[Offset] = Sum;
            
            for(i = -x, Sum = 0; i <= QUAD_RES - x; i++)
            {
                if(abs(i) <= R)
                {
                    j = QUAD_RES - y;
                    
                    if(abs(j) <= R)
                        Sum += PsfSamples[(i + R) + W*(j + R)];
                    
                    j = -y;
                    
                    if(abs(j) <= R)
                        Sum -= PsfSamples[(i + R) + W*(j + R)];
                }
            }
            
            QuadDy[Offset] = Sum;
        }
    
    Set->QuadDx = QuadDx;
    Set->QuadDy = QuadDy;
    Free(PsfSamples);
    return 1;
}


/** @brief phi function for the circle stencil */
static double Circle(double x, double y, ATTRIBUTE_UNUSED const void *Param)
{
    return sqrt(x*x + y*y);
}


/**
 * @brief Create a new stencil set
 * @param StencilRadius the circular radius of \f$\Omega\f$
 * @param Psf point spread function modeling the sampling of the input image
 * @param PsfParam parameters to pass to Psf
 * @param PsfRadius support radius of Psf
 * @return an sset* on success, NULL on failure
 */
sset *NewStencilSet(double StencilRadius,
    double (*Psf)(double, double, const void*),
    const void *PsfParam, double PsfRadius)
{
    sset *Set = NULL;
    int x, y, n, R = (int)floor(StencilRadius);
            
    
    if(R <= 0 || PsfRadius < 0
        || !(Set = (sset *)Malloc(sizeof(struct ssetstruct))))
        return NULL;
    
    /* Initializes an empty stencil set */
    Set->NumStencils = Set->Capacity = Set->NumCells = 0;
    Set->QuadDx = NULL;
    Set->QuadDy = NULL;
    Set->Cell = NULL;
    Set->Stencil = NULL;
    Set->StencilTable = NULL;
    Set->Psf = Psf;
    Set->PsfParam = PsfParam;
    Set->PsfRadius = (Psf) ? PsfRadius : 0;
    
    /* Allocate and precompute cell quadrature weights */
    if(!ComputeQuadWeights(Set)
        || !(Set->Cell = (stencilsetcell *)
        Malloc(sizeof(stencilsetcell)*4*R*R)))
    {
        FreeStencilSet(Set);
        return NULL;
    }
    
    /* Construct the circular neighborhood */
    for(y = -R, n = 0; y < R; y++)
        for(x = -R; x < R; x++)
            if((x + 0.5)*(x + 0.5) + (y + 0.5)*(y + 0.5)
                <= StencilRadius*StencilRadius)
            {
                Set->Cell[n].x = x;
                Set->Cell[n].y = y;
                n++;
            }
            
    Set->NumCells = n;
    
    if(!AddStencil(Set, Circle, NULL, 0, NULL, NULL, NULL))
    {
        FreeStencilSet(Set);
        return NULL;
    }
    
    return Set;
}


/**
 * @brief Get the number of stencils in a stencil set
 * @param Set the stencil set
 * @return the number of stencils in Set
 */
int GetNumStencils(const sset *Set)
{
    return (Set) ? Set->NumStencils : 0;
}


/**
 * @brief Get the PSF
 * @param Set the stencil set
 * @return function pointer to the PSF
 */
double (*GetPsf(const sset *Set))(double, double, const void*)
{
    return (Set) ? Set->Psf : NULL;
}


/**
 * @brief Get the PSF parameter pointer
 * @param Set the stencil set
 * @return the PSF parameters
 */
const void *GetPsfParam(const sset *Set)
{
    return (Set) ? Set->PsfParam : NULL;
}


/**
 * @brief Get the support radius of the PSF
 * @param Set the stencil set
 * @return the PSF support radius
 */
double GetPsfRadius(const sset *Set)
{
    return (Set) ? Set->PsfRadius : 0;
}


/**
 * @brief Evaluate the \f$\varphi_\mathcal{S}\f$ associated with a stencil
 * @param Set the stencil set
 * @param S index of the stencil
 * @param x,y location to evaluate
 *
 * Evaluate the distance function \f$\varphi_\mathcal{S}\f$ associated with
 * stencil S at (x,y).  The stencil rotation from stencilentry::PhiTrans is
 * applied.
 */
double EvalPhi(const sset *Set, int S, double x, double y)
{
    const stencilentry *Entry = &Set->Stencil[S];
    
    return Entry->Phi(Entry->PhiTrans[0]*x + Entry->PhiTrans[1]*y,
        -Entry->PhiTrans[1]*x + Entry->PhiTrans[0]*y, Entry->PhiParam);
}


/**
 * @brief Add a stencil to a stencil set
 * @param Set the stencil set
 * @param Phi the distance function from which to derive the stencil
 * @param PhiParam parameters for the distance function
 * @param Rotation rotation of the stencil
 * @param DrawFun function that draws the stencil
 * @param DrawParam parameter for DrawFun
 * @param DrawColor RGB color for drawing the stencil
 * @return 1 on success, 0 on failure
 */
int AddStencil(sset *Set, double (*Phi)(double, double, const void*),
    const void *PhiParam, double Rotation,
    int (*DrawFun)(pen *Pen, const void *Param), const void *DrawParam,
    const float *DrawColor)
{
    stencilentry *Entry;
    const double *QuadDx, *QuadDy;
    double Sum, Alpha, Beta, Temp;
    int n, x, y, i, R, S, Offset;
    
    
    if(!Set || !Set->Stencil || Set->Capacity <= Set->NumStencils)
        if(!ReallocStencilSet(Set))  /* Allocate more memory if needed */
            return 0;
    
    S = Set->NumStencils;
    Set->NumStencils++;
    Entry = &Set->Stencil[S];
    Entry->Phi = Phi;
    Entry->PhiParam = PhiParam;
    Entry->PhiTrans[0] = cos(Rotation);
    Entry->PhiTrans[1] = sin(Rotation);
    Entry->DrawFun = DrawFun;
    Entry->DrawParam = DrawParam;
    Entry->DrawColor = DrawColor;
    
    R = (int)ceil(QUAD_RES * Set->PsfRadius);
    QuadDx = Set->QuadDx;
    QuadDy = Set->QuadDy;
    
    for(n = 0, Sum = 0; n < Set->NumCells; n++)
    {
        Alpha = Beta = 0;
        
        /* Compute the average of (grad (Psf * Phi))^perp over the cell */
        for(y = -R, Offset = 0; y <= QUAD_RES + R; y++)
            for(x = -R; x <= QUAD_RES + R; x++, Offset++)
            {
                Temp = EvalPhi(Set, S, Set->Cell[n].x + ((double)x)/QUAD_RES,
                               Set->Cell[n].y + ((double)y)/QUAD_RES);
                Beta -= QuadDx[Offset] * Temp;
                Alpha += QuadDy[Offset] * Temp;
            }
        
        Entry->Alpha[n] = Alpha;
        Entry->Beta[n] = Beta;
        Sum += sqrt(Alpha*Alpha + Beta*Beta);
    }
    
    /* Normalize the stencil */
    for(n = 0; n < Set->NumCells; n++)
    {
        Entry->Alpha[n] /= Sum;
        Entry->Beta[n] /= Sum;
    }
    
    /* Quantize the stencil */
    for(n = 0; n < Set->NumCells; n++)
    {
        i = (int)ROUND(atan2(Entry->Beta[n], Entry->Alpha[n])*NUMANGLES/M_PI);
        
        while(i < 0)
            i += NUMANGLES;
        while(i >= NUMANGLES)
            i -= NUMANGLES;
        
        Set->StencilTable[n + Set->NumCells*S].AngleIndex = i;
        Set->StencilTable[n + Set->NumCells*S].Weight =
            (float)(sqrt(Entry->Alpha[n]*Entry->Alpha[n]
            + Entry->Beta[n]*Entry->Beta[n]));
    }
    
    return 1;
}


/**
 * @brief Compute the confusion matrix of a stencil set
 * @param ConfusionMatrix pointer to an array to hold the result
 * @param Set the stencil set
 * @return 1 on success, 0 on failure
 *
 * For each distance function \f$\varphi\f$ used in the stencil set, all
 * stencil TVs are computed.  The results are stored in column-major order in
 * ConfusionMatrix as
 *
 * \f[ \text{\tt ConfusionMatrix[i + j*GetNumStencils(Set)]} =
 *      V(\mathcal{S}_j, h * \varphi_i). \f]
 *
 * The stencil associated with a distance function should have the lowest TV
 * value, that is, the diagonal entries should be the smallest in each column.
 * Small off-diagonal entries indicate stencils that are likely to be confused.
 */
int StencilConfusion(double *ConfusionMatrix, sset *Set)
{
    const int R = (int)ceil(QUAD_RES * Set->PsfRadius);
    const int W = 2*R + 1;
    double *PsfSamples = NULL, *Psi = NULL;
    double Sum, Alpha, Beta, a, b, c, d;
    int MinX, MaxX, MinY, MaxY, PsiWidth, PsiHeight;
    int i, j, n, nx, ny, x, y, Success = 0;

    MinX = Set->Cell[0].x;
    MaxX = MinX + 1;
    MinY = Set->Cell[0].y;
    MaxY = MinY + 1;
    
    for(n = 1; n < Set->NumCells; n++)
    {
        if(Set->Cell[n].x >= MaxX)
            MaxX = Set->Cell[n].x + 1;
        else if(Set->Cell[n].x < MinX)
            MinX = Set->Cell[n].x;
        
        if(Set->Cell[n].y >= MaxY)
            MaxY = Set->Cell[n].y + 1;
        else if(Set->Cell[n].y < MinY)
            MinY = Set->Cell[n].y;
    }
    
    PsiWidth = MaxX - MinX + 1;
    PsiHeight = MaxY - MinY + 1;
    
    if(!(PsfSamples = SamplePsf(Set))
        || !(Psi = (double *)Malloc(sizeof(double)*PsiWidth*PsiHeight)))
        goto Catch;
    
    for(i = 0; i < Set->NumStencils; i++)
    {
        /* Compute psi_i = (h * phi_i) */
        for(ny = MinY; ny <= MaxY; ny++)
            for(nx = MinX; nx <= MaxX; nx++)
            {
                for(y = -R, Sum = 0; y <= R; y++)
                    for(x = -R; x <= R; x++)
                    {
                        Sum += PsfSamples[(R - x) + W*(R - y)]
                            * EvalPhi(Set, i, nx + ((double)x)/QUAD_RES,
                            ny + ((double)y)/QUAD_RES);
                    }
                
                Psi[(nx - MinX) + PsiWidth*(ny - MinY)] = Sum;
            }
        
        /* Compute the stencil TVs, V(S_j, psi_i) */
        for(j = 0; j < Set->NumStencils; j++)
        {
            for(n = 0, Sum = 0; n < Set->NumCells; n++)
            {
                nx = Set->Cell[n].x;
                ny = Set->Cell[n].y;
                Alpha = Set->Stencil[j].Alpha[n];
                Beta = Set->Stencil[j].Beta[n];
                
                /* Get the four corners of the cell
                 *
                 *    a---b
                 *    |   |
                 *    c---d
                 */
                a = Psi[(nx - MinX) + PsiWidth*(ny + 1 - MinY)];
                b = Psi[(nx + 1 - MinX) + PsiWidth*(ny + 1 - MinY)];
                c = Psi[(nx - MinX) + PsiWidth*(ny - MinY)];
                d = Psi[(nx + 1 - MinX) + PsiWidth*(ny - MinY)];
                                
                if(Alpha*Beta >= 0)
                    Sum += fabs(Alpha*b - (Alpha - Beta)*a - Beta*c)
                        + fabs(Beta*b + (Alpha - Beta)*d - Alpha*c);
                else
                    Sum += fabs(Alpha*a - (Alpha + Beta)*b + Beta*d)
                        + fabs(Beta*a - (Alpha + Beta)*c + Alpha*d);
            }
            
            ConfusionMatrix[i + Set->NumStencils*j] = Sum;
        }
    }
    
    Success = 1;
Catch:
    Free(Psi);
    Free(PsfSamples);
    return Success;
}


/**
 * @brief Print stencil confusion information about a stencil set
 * @param Set the stencil set
 * @return 1 on success, 0 on failure
 *
 * Prints information from the stencil confusion matrix computed by
 * StencilConfusion().
 */
int PrintStencilConfusion(sset *Set)
{
    double *ConfusionMatrix = NULL;
    double TV, MinTV, Min2TV;
    int i, j, MinS, Success = 0;
    
    
    if(!(ConfusionMatrix = (double *)Malloc(sizeof(double)*
        GetNumStencils(Set)*GetNumStencils(Set)))
        || !StencilConfusion(ConfusionMatrix, Set))
        goto Catch;
       
    /* Print table header */
    printf("Stencil confusion for the ith distance function phi_i\n"
           "  i     V(S_i)   S^*   V(S^*)   V(S^*2)   S^*2\n"
           "------------------------------------------------------------\n");
    
    for(i = 0; i < GetNumStencils(Set); i++)
    {
        MinS = -1;
        MinTV = Min2TV = 1e30;
        
        /* Find the best and second best stencil TVs */
        for(j = 0; j < GetNumStencils(Set); j++)
        {
            TV = ConfusionMatrix[i + GetNumStencils(Set)*j];
                
            if(TV < Min2TV)
            {
                if(TV < MinTV)
                {
                    Min2TV = MinTV;
                    MinS = j;
                    MinTV = TV;
                }
                else
                    Min2TV = TV;
            }
        }
        
        printf("%3d   %8.4f   %3d %8.4f  %8.4f   {",
               i, ConfusionMatrix[i + GetNumStencils(Set)*i],
               MinS, MinTV, Min2TV);
        
        /* Print all the indices of stencil who have TV equal to Min2TV */
        for(j = 0; j < GetNumStencils(Set); j++)
            if(fabs(ConfusionMatrix[i + GetNumStencils(Set)*j] - Min2TV)
                < 1e-12)
                printf("%3d,", j);
        
        printf("\b}\n");
    }
    
    Success = 1;
Catch:
    Free(ConfusionMatrix);
    return Success;
}


/**
 * @brief Convert an RGB color to YPbPr
 * @param Src the input RGB color
 * @param Dest the resulting YPbPr color
 */
static void Rgb2YPbPr(float *Dest, const float *Src)
{
    float Red = ((float *)Src)[0];
    float Green = ((float *)Src)[1];
    float Blue = ((float *)Src)[2];
        
    Dest[0] =  0.299f*Red    + 0.587f*Green    + 0.114f*Blue;
    Dest[1] = -0.168736f*Red - 0.331264f*Green + 0.5f*Blue;
    Dest[2] =  0.5f*Red      - 0.418688f*Green - 0.081312f*Blue;
}


/**
 * @brief Determine the best-fitting stencil at each point in an image
 * @param Stencil the selected best-fitting stencils \f$\mathcal{S}^\star\f$
 * @param Set stencil set to use
 * @param Image the input image
 * @param ImageWidth, ImageHeight image dimensions
 *
 * This is the main routine for determining the best-fitting stencils on an
 * image.  The result is stored in Stencil, which should be an array of size
 * ImageWidth by ImageHeight.  For each pixel, the index of the best-fitting
 * stencil is stored in Stencil.
 *
 * For computational efficiency, the stencil fitting uses the quantized stencil
 * vectors ssetstruct::StencilTable to approximate the cell total variations
 * \f[ V_{i,j}(\mathcal{S}, v) \approx |\mathcal{S}_{i,j}| \cdot
 *     V_{i,j}\bigl( (\!\begin{smallmatrix} \cos \theta \\
 *     \sin \theta \end{smallmatrix}\!), v\bigr), \quad
 *     |\theta - \angle \mathcal{S}_{i,j}| \le \tfrac{\pi}{128}. \f]
 * First, the cell total variation are computed for each pixel at each angle
 * 0, ..., #NUMANGLES-1.  Then then stencils are computed by summing the cell
 * total variations.  Finally, the stencil with the smallest total variation is
 * identified as the best-fitting stencil,
 * \f[ \mathcal{S}^\star(k) = \operatorname*{arg\,min}_{\mathcal{S}\in\Sigma}
 *                            \, V\bigl(\mathcal{S},v(\cdot-k)\bigr). \f]
 * If the difference between the smallest and second smallest total variations
 * is below a threshold, the circle stencil is selected.
 */
void FitStencils(int *Stencil, const sset *Set,
    const float *Image, int ImageWidth, int ImageHeight)
{
    const float TvThresh = (float)(STABILITY_THRESH*(4*sqrt(2)));
    float *TVCell = NULL;
    float *TVCellPtr;
    stencilqvec *TablePtr;
    const int CellsPerRow = ImageWidth - 1;
    const float *ImagePtr;
    const stencilsetcell *Cell;
    int *NeighborOffset = NULL;
    float TV, MinTV, Min2TV, Alpha[NUMANGLES], Beta[NUMANGLES];
    float a[3], b[3], c[3], d[4], ab[3], ac[3], bd[3], cd[3];
    int x, y, n, nx, ny, S, BestS, NumCells, NumStencils;
    
    
    if(!Set || !(TVCell = (float *)Malloc(sizeof(float)*
        NUMANGLES*CellsPerRow*(ImageHeight - 1)))
        || !(NeighborOffset = (int *)Malloc(sizeof(int)*Set->NumCells)))
        goto Catch;
    
    NumCells = Set->NumCells;
    NumStencils = Set->NumStencils;
    Cell = Set->Cell;
    
    for(n = 0; n < NUMANGLES; n++)
    {
        Alpha[n] = (float)cos((M_PI*n)/NUMANGLES);
        Beta[n] = (float)sin((M_PI*n)/NUMANGLES);
    }
    
    for(n = 0; n < NumCells; n++)
        NeighborOffset[n] = NUMANGLES*(Cell[n].x - CellsPerRow*Cell[n].y);
    
    for(y = 0, TVCellPtr = TVCell; y < ImageHeight - 1; y++)
    {
        for(x = 0; x < CellsPerRow; x++, TVCellPtr += NUMANGLES)
        {
            ImagePtr = Image + 3*(x + ImageWidth*y);
            /* Get the four corners of the cell in YPbPr representation
             *
             *    a---b
             *    |   |
             *    c---d
             */
            Rgb2YPbPr(a, ImagePtr);
            Rgb2YPbPr(b, ImagePtr + 3);
            Rgb2YPbPr(c, ImagePtr + 3*ImageWidth);
            Rgb2YPbPr(d, ImagePtr + 3 + 3*ImageWidth);
            
            ab[0] = a[0] - b[0];
            ab[1] = a[1] - b[1];
            ab[2] = a[2] - b[2];
            ac[0] = a[0] - c[0];
            ac[1] = a[1] - c[1];
            ac[2] = a[2] - c[2];
            bd[0] = b[0] - d[0];
            bd[1] = b[1] - d[1];
            bd[2] = b[2] - d[2];
            cd[0] = c[0] - d[0];
            cd[1] = c[1] - d[1];
            cd[2] = c[2] - d[2];
            
            /* Horizontal TV */
            TVCellPtr[0] = (float)(fabs(ab[0]) + fabs(ab[1]) + fabs(ab[2])
                + fabs(cd[0]) + fabs(cd[1]) + fabs(cd[2]));
            
            /* Angles between 0 and pi/2 */
            for(n = 1; n < (NUMANGLES/2); n++)
                TVCellPtr[n] = (float)(
                    fabs(Alpha[n]*ab[0] - Beta[n]*ac[0]) +
                    fabs(Alpha[n]*ab[1] - Beta[n]*ac[1]) +
                    fabs(Alpha[n]*ab[2] - Beta[n]*ac[2]) +
                    fabs(Beta[n]*bd[0] - Alpha[n]*cd[0]) +
                    fabs(Beta[n]*bd[1] - Alpha[n]*cd[1]) +
                    fabs(Beta[n]*bd[2] - Alpha[n]*cd[2]));

            /* Vertical TV */
            TVCellPtr[n++] = (float)(fabs(bd[0]) + fabs(bd[1]) + fabs(bd[2])
                + fabs(ac[0]) + fabs(ac[1]) + fabs(ac[2]));
                    
            /* Angles between pi/2 and pi */
            for(; n < NUMANGLES; n++)
                TVCellPtr[n] = (float)(
                    fabs(Alpha[n]*ab[0] - Beta[n]*bd[0]) +
                    fabs(Alpha[n]*ab[1] - Beta[n]*bd[1]) +
                    fabs(Alpha[n]*ab[2] - Beta[n]*bd[2]) +
                    fabs(Beta[n]*ac[0] - Alpha[n]*cd[0]) +
                    fabs(Beta[n]*ac[1] - Alpha[n]*cd[1]) +
                    fabs(Beta[n]*ac[2] - Alpha[n]*cd[2]));
        }
    }
    
    for(y = 0; y < ImageHeight; y++)
    {
        for(x = 0; x < ImageWidth; x++)
        {
            TVCellPtr = TVCell + NUMANGLES*(x + CellsPerRow*(y - 1));
            TablePtr = Set->StencilTable;
            BestS = 0;
            MinTV = Min2TV = 1e30f;
            
            for(S = 0; S < NumStencils; S++, TablePtr += NumCells)
            {
                TV = 0;
                
                /* Compute the stencil TV for the Sth stencil */
                for(n = 0; n < NumCells; n++)
                {
                    nx = x + Cell[n].x;
                                        
                    if(nx >= 0 && nx < CellsPerRow)
                    {
                        ny = y - Cell[n].y - 1;
                        
                        if(ny >= 0 && ny < ImageHeight - 1)
                            TV += TablePtr[n].Weight*TVCellPtr[
                                NeighborOffset[n] + TablePtr[n].AngleIndex ];
                    }
                }
                
                /* Determine the best-fitting stencil */
                if(TV < Min2TV)
                {
                    if(TV < MinTV)
                    {
                        BestS = S;
                        Min2TV = MinTV;
                        MinTV = TV;
                    }
                    else
                        Min2TV = TV;
                }
                
                /*if(TV <= MinTV)
                {
                    BestS = S;
                    Min2TV = MinTV;
                    MinTV = TV;
                }*/
            }
            
            /* If the separation between smallest and second smallest TVs
               is below the threshold, select the circle stencil. */
            if(Min2TV - MinTV < TvThresh)
                BestS = 0;
            
            *Stencil = BestS;
            Stencil++;
        }
    }
    
Catch:
    Free(NeighborOffset);
    Free(TVCell);
}


/***************************************************************************/


static float *MakeContoursBgImage(int OutputWidth, int OutputHeight,
    double ScaleFactor, const float *Input, int InputWidth, int InputHeight,
    int IsColor)
{
    const int NumChannels = (IsColor) ? 3 : 1;
    const int OutputStride = NumChannels*OutputWidth;
    float *Output, *OutputPtr;
    float Gray;
    int x, y, bx, bx1, bx2, by, by1, by2;
    
    if(!(Output = (float *)Malloc(sizeof(float)*OutputStride*OutputHeight)))
        return NULL;
    
    for(y = 0, by1 = 0; y < InputHeight; y++)
    {
        by2 = (y < InputHeight - 1) ?
            (int)ROUND(ScaleFactor*(y + 1)) : OutputHeight;
        
        for(x = 0, bx1 = 0; x < InputWidth; x++)
        {
            /* Darken and desaturate the pixel */
            Gray = 0.05f + 0.7f*(
                0.299f*Input[3*(x + InputWidth*y) + 0] +
                0.587f*Input[3*(x + InputWidth*y) + 1] +
                0.114f*Input[3*(x + InputWidth*y) + 2]);
            
            bx2 = (x < InputWidth - 1) ?
                ((int)ROUND(ScaleFactor*(x + 1))) : OutputWidth;
            
            if(IsColor)
                for(by = by1, OutputPtr = Output + OutputStride*by; by < by2;
                    by++, OutputPtr += OutputStride)
                    for(bx = bx1; bx < bx2; bx++)
                    {
                        OutputPtr[3*bx + 0] = Gray;
                        OutputPtr[3*bx + 1] = Gray;
                        OutputPtr[3*bx + 2] = Gray;
                    }
            else
                for(by = by1, OutputPtr = Output + OutputStride*by; by < by2;
                    by++, OutputPtr += OutputStride)
                    for(bx = bx1; bx < bx2; bx++)
                        OutputPtr[bx] = Gray;
            
            bx1 = bx2;
        }
        
        by1 = by2;
    }
        
    return Output;
}


/**
 * @brief Visualize the estimated contours
 * @param Output output image
 * @param OutputWidth, OutputHeight dimensions of the output image
 * @param ScaleFactor scale factor between input and output images
 * @param Stencil best-fitting stencils \f$\mathcal{S}^\star\f$
 * @param Set the stencil set
 * @param Input input image
 * @param InputWidth, InputHeight dimensions of the input image
 *
 * Renders a visualization of the estimated contours to the output image.  A
 * nearest neighbor interpolation of the input image is created where the colors
 * are darkened and converted to grayscale.  For each pixel, the best-fitting
 * stencil is drawn by calling the stencil's stencilentry::DrawFun.
 */
int DrawContours(const char *FileName, int OutputWidth, int OutputHeight,
    double ScaleFactor, const int *Stencil, const sset *Set,
    const float *Input, int InputWidth, int InputHeight)
{
    pen *Pen = NewPen();
    pentrans Trans;
    float *Output = NULL;
    stencilentry *S;
    int i, x, y, Success = 0;
    
    if(!Pen
        || !(Output = MakeContoursBgImage(OutputWidth, OutputHeight,
            ScaleFactor, Input, InputWidth, InputHeight, 1))
        || !PenRenderToImage(Pen, Output, OutputWidth, OutputHeight))
        goto Catch;
    
    for(i = 1; i < GetNumStencils(Set); i++)
    {
        S = &Set->Stencil[i];
        PenSetColor(Pen, S->DrawColor);
        
        for(y = 0; y < InputHeight; y++)
            for(x = 0; x < InputWidth; x++)
                if(i == Stencil[x + InputWidth*y])
                {
                    Trans = PenGetTrans(Pen);
                    PenTransformCanvas(Pen,
                        ScaleFactor*S->PhiTrans[0],
                        -ScaleFactor*S->PhiTrans[1],
                        -ScaleFactor*S->PhiTrans[1],
                        -ScaleFactor*S->PhiTrans[0],
                        ScaleFactor*(x + 0.5),
                        ScaleFactor*(y + 0.5));
                    
                    if(S->DrawFun && !(S->DrawFun(Pen, S->DrawParam)))
                    {
                        ErrorMessage("Error drawing stencil %d at (%d,%d).\n",
                            i, x, y);
                        goto Catch;
                    }
                    
                    PenSetTrans(Pen, Trans);
                }
    }
    
    /* Write the output image */
    if(!WriteImage(Output, OutputWidth, OutputHeight, FileName,
        IMAGEIO_FLOAT | IMAGEIO_RGB, 85))
        goto Catch;
    
    Success = 1;
Catch:
    Free(Output);
    FreePen(Pen);
    return Success;
}


/**
 * @brief Visualize the estimated contours as an EPS image
 * @param EpsFile output EPS file name
 * @param Stencil best-fitting stencils \f$\mathcal{S}^\star\f$
 * @param Set the stencil set
 * @param Input input image
 * @param InputWidth, InputHeight dimensions of the input image
 * @return 1 on success, 0 on failure
 *
 * Renders a visualization of the estimated contours as an EPS image.  For each
 * pixel, the best-fitting stencil is drawn by calling the stencil's
 * stencilentry::DrawFun.
 */
int DrawContoursEps(const char *EpsFile,
    const int *Stencil, const sset *Set,
    const float *Input, int InputWidth, int InputHeight)
{
    const int DimScale = 6;
    pen *Pen = NewPen();
    pentrans Trans;
    stencilentry *S;
    uint8_t *BgImage = NULL;
    int i, x, y, Success = 0;
    
    if(!Pen
        || !(BgImage = (uint8_t *)Malloc(
            sizeof(uint8_t)*InputWidth*InputHeight))
        || !EpsOpen(Pen, EpsFile,
            DimScale*InputWidth, DimScale*InputHeight))
        goto Catch;
    
    fprintf(PenGetFile(Pen), "%d %d scale\n", DimScale, DimScale);
    
    /* Darken and desaturate image */
    for(y = 0; y < InputHeight; y++)
        for(x = 0; x < InputWidth; x++)
            BgImage[x + InputWidth*y] =
                (uint8_t)floor(255*(0.05f + 0.7f*(
                    0.299f*Input[3*(x + InputWidth*y) + 0]
                    + 0.587f*Input[3*(x + InputWidth*y) + 1]
                    + 0.114f*Input[3*(x + InputWidth*y) + 2])) + 0.5f);
    
    if(!EpsWriteGrayImage(PenGetFile(Pen),
        BgImage, InputWidth, InputHeight))
        goto Catch;
    
    PenTranslateCanvas(Pen, 0, InputHeight);
    PenScaleCanvas(Pen, 1, -1);
                
    for(i = 1; i < GetNumStencils(Set); i++)
    {
        S = &Set->Stencil[i];
        PenSetColor(Pen, S->DrawColor);
        
        for(y = 0; y < InputHeight; y++)
            for(x = 0; x < InputWidth; x++)
                if(i == Stencil[x + InputWidth*y])
                {
                    Trans = PenGetTrans(Pen);
                    PenTransformCanvas(Pen,
                        S->PhiTrans[0], -S->PhiTrans[1],
                        -S->PhiTrans[1], -S->PhiTrans[0],
                        (x + 0.5), (y + 0.5));
                    
                    if(S->DrawFun && !(S->DrawFun(Pen, S->DrawParam)))
                    {
                        ErrorMessage("Error writing \"%s\".\n",
                            EpsFile);
                        goto Catch;
                    }
                    
                    PenSetTrans(Pen, Trans);
                }
    }
    
    if(!EpsClose(Pen))
    {
        ErrorMessage("Error writing \"%s\".\n", EpsFile);
        goto Catch;
    }
    
    Success = 1;
Catch:
    Free(BgImage);
    FreePen(Pen);
    return Success;
}


/**
 * @brief Visualize the estimated contours as an SVG image
 * @param SvgFile output SVG file name
 * @param BgFile output background image file name (PNG or JPEG)
 * @param Stencil best-fitting stencils \f$\mathcal{S}^\star\f$
 * @param Set the stencil set
 * @param Input input image
 * @param InputWidth, InputHeight dimensions of the input image
 * @return 1 on success, 0 on failure
 *
 * Renders a visualization of the estimated contours as an SVG image.  For each
 * pixel, the best-fitting stencil is drawn by calling the stencil's
 * stencilentry::DrawFun.
 */
int DrawContoursSvg(const char *SvgFile, const char *BgFile,
    const int *Stencil, const sset *Set,
    const float *Input, int InputWidth, int InputHeight)
{
    const int DimScale = 9;
    pen *Pen = NewPen();
    pentrans Trans;
    stencilentry *S;
    float *BgImage = NULL;
    int i, x, y, Success = 0;
    
    if(!Pen)
        goto Catch;
    
    PenSetNumDigits(Pen, 0);
    
    for(i = 1; i < GetNumStencils(Set); i++)
        if(!PenDefineColor(Pen, Set->Stencil[i].DrawColor))
            goto Catch;
    
    if(!SvgOpen(Pen, SvgFile, DimScale*InputWidth, DimScale*InputHeight))
        goto Catch;
    
    /* Make the background image */
    if(BgFile)
    {
        if(!(BgImage = MakeContoursBgImage(InputWidth, InputHeight,
            1, Input, InputWidth, InputHeight, 0))
            || !WriteImage(BgImage, InputWidth, InputHeight, BgFile,
                IMAGEIO_FLOAT | IMAGEIO_GRAYSCALE, 85))
            goto Catch;
                
        fprintf(PenGetFile(Pen),
            "<image x=\"0\" y=\"0\" width=\"%d\" height=\"%d\"\n"
            "  xlink:href=\"%s\"\n"
            "  style=\"image-rendering:-moz-crisp-edges;"
            " -ms-interpolation-mode: nearest-neighbor;\" />\n",
            DimScale*InputWidth, DimScale*InputHeight, BgFile);
    }
    
    for(i = 1; i < GetNumStencils(Set); i++)
    {
        S = &Set->Stencil[i];
        PenSetColor(Pen, S->DrawColor);
        
        for(y = 0; y < InputHeight; y++)
            for(x = 0; x < InputWidth; x++)
                if(i == Stencil[x + InputWidth*y])
                {
                    Trans = PenGetTrans(Pen);
                    PenTransformCanvas(Pen,
                        DimScale*S->PhiTrans[0], -DimScale*S->PhiTrans[1],
                        -DimScale*S->PhiTrans[1], -DimScale*S->PhiTrans[0],
                        DimScale*(x + 0.5), DimScale*(y + 0.5));
                    
                    if(S->DrawFun && !(S->DrawFun(Pen, S->DrawParam)))
                    {
                        ErrorMessage("Error writing \"%s\".\n",
                            SvgFile);
                        goto Catch;
                    }
                    
                    PenSetTrans(Pen, Trans);
                }
    }
    
    if(!SvgClose(Pen))
    {
        ErrorMessage("Error writing \"%s\".\n", SvgFile);
        goto Catch;
    }
    
    Success = 1;
Catch:
    Free(BgImage);
    FreePen(Pen);
    return Success;
}
