/**
 * @file sinterp.c
 * @brief Contour stencil interpolation
 * @author Pascal Getreuer <getreuer@gmail.com>
 *
 * This file implements contour stencil interpolation for color images by
 * integer and non-integer scale factors.  The main routines are
 *  - NewSInterp() to compute data to use a stencil set for interpolation,
 *  - Prefilter() to perform iterative refinement,
 *  - IntegerScalePass() for interpolation by an integer scale factor,
 *  - ArbitraryScale() for interpolation by a non-integer scale factor.
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

#include <string.h>
#include "invmat.h"
#include "sinterp.h"


/** @brief The size of the synthesis neighborhood \f$\mathcal{N}\f$ */
#define NEIGHRADIUS         1
/** @brief Quadrature resolution to use for integrals */
#define QUAD_RES            8
/** @brief Size of the trignometry tables used in RhoFast */
#define TRIG_TABLE_SIZE     256
/** @brief Size of the exponential table used in RhoFast */
#define EXP_TABLE_SIZE      1024

/* Quantities derived from NEIGHRADIUS */
/** @brief Support radius of the window function \f$w(x)\f$ */
#define WINDOWRADIUS            (1 + NEIGHRADIUS)
/** @brief Diameter of the neighborhood */
#define NEIGHDIAMETER           (1 + 2*NEIGHRADIUS)
/** @brief Index of the neighborhood center */
#define NEIGHCENTER             (NEIGHRADIUS + NEIGHRADIUS*NEIGHDIAMETER)
/** @brief Number of neighbors in the neighborhood, \f$|\mathcal{N}|\f$ */
#define NUMNEIGH                (NEIGHDIAMETER*NEIGHDIAMETER)


/** @brief Data for interpolating with a contour stencil */
typedef struct
{
    /** @brief Synthesis function orientations \f$\theta_\mathcal{S}^n\f$ */
    float Orientation[NUMNEIGH];
    /** @brief Interpolation matrix \f$A_\mathcal{S}^{-1}\f$ */
    float Matrix[NUMNEIGH*NUMNEIGH];
    /** @brief Anisotropy \f$\mu_\mathcal{S}\f$ */
    float Anisotropy;
} sinterpentry;


/** @brief Data for interpolating with a contour stencil set */
typedef struct sinterpstruct
{
    /** @brief The stencil set \f$\Sigma\f$ */
    const sset *StencilSet;
    /** @brief Interpolation data for each \f$\mathcal{S}\in\Sigma\f$ */
    sinterpentry *StencilInterp;
    /** @brief Spread of \f$\rho\f$ in the tagential direction */
    float RhoSigmaTangent;
    /** @brief Spread of \f$\rho\f$ in the normal direction */
    float RhoSigmaNormal;
} sinterptype;


/** @brief The oriented function rho used in local reconstructions */
static double Rho(double x, double y, double theta,
    float SigmaTangent, float SigmaNormal, float Anisotropy)
{
    double t, n;

    SigmaNormal = Anisotropy*SigmaNormal + (1 - Anisotropy)*SigmaTangent;
    t = (cos(theta)*x + sin(theta)*y) / SigmaTangent;
    n = (-sin(theta)*x + cos(theta)*y) / SigmaNormal;
    return exp(-0.5*(t*t + n*n));
}


/** @brief The oriented function rho used in local reconstructions */
static float RhoFast(float x, float y, float theta,
    float SigmaTangent, float SigmaNormal, float Anisotropy)
{
    static float ExpTable[EXP_TABLE_SIZE];
    static float CosTable[TRIG_TABLE_SIZE];
    static float SinTable[TRIG_TABLE_SIZE];
    static int TablesInitialized = 0;
    static const float ExpArgScale = 37.0236f;
    static const float TrigScale = (float)((TRIG_TABLE_SIZE - 1)/M_PI);
    float t, n;
    int i;
            
    
    if(!TablesInitialized)
    {
        /* Precompute exp and trig tables */
        for(i = 0; i < EXP_TABLE_SIZE; i++)
            ExpTable[i] = (float)exp(-(i - 0.5)/ExpArgScale);
        
        for(i = 0; i < TRIG_TABLE_SIZE; i++)
        {
            CosTable[i] = (float)cos(i/TrigScale);
            SinTable[i] = (float)sin(i/TrigScale);
        }
        
        TablesInitialized = 1;
    }
    
    SigmaNormal = Anisotropy*SigmaNormal + (1 - Anisotropy)*SigmaTangent;
    
    i = (int)(theta*TrigScale + 0.5f);
    t = (CosTable[i]*x + SinTable[i]*y) / SigmaTangent;
    n = (-SinTable[i]*x + CosTable[i]*y) / SigmaNormal;
    t = (t*t + n*n)*ExpArgScale/2;
    
    if(t >= EXP_TABLE_SIZE)
        return 0;
    else
        return ExpTable[(int)t];
    
    /*return exp(-0.5*(t*t + n*n));*/
}


/** @brief The window function, a cubic B-spline */
static float Window(float x)
{
    x = (float)fabs(x)/(((float)WINDOWRADIUS)/2);
    
    if(x < 1)
        return (x/2 - 1)*x*x + 0.66666666666666667f;
    else if(x < 2)
    {
        x = 2 - x;
        return x*x*x/6;
    }
    else
        return 0;
}


/**
 * @brief Compute interpolation information for a stencil
 * @param SInterp stencil set interpolation data
 * @param S index of the stencil to interpolate
 * @param PsfSamples precomputed samples of the PSF
 * @return 1 on success, 0 on failure
 *
 * The orientation and anisotropy of the \f$\rho_\mathcal{S}^n\f$ synthesis
 * functions are computed,
 * \f[\begin{aligned}
 *    \theta_\mathcal{S}^n &= \angle \int_{[-\frac{1}{2},+\frac{1}{2}]^2}
 *    \nabla\varphi_\mathcal{S}^\perp \bigl(x - n) \, dx, \\
 *    \mu_\mathcal{S} &= \min_{n\in \mathcal{N}}\,
 *    \Bigl\lvert\int_{[-\frac{1}{2},+\frac{1}{2}]^2} \nabla\varphi_\mathcal{S}
 *    \bigl(x - n) \, dx\Bigr\rvert.
 *    \end{aligned}\f]
 * The interpolation matrix
 * \f[ (A_\mathcal{S})_{m,n} = (h * \rho^n_\mathcal{S})(m - n),
 *     \quad m,n\in\mathcal{N} \f]
 * is computed and inverted with InvertMatrix().
 */
int InterpolateStencil(sinterp *SInterp, int S,
    const double *PsfSamples)
{
    const sset *StencilSet;
    sinterpentry *Entry;
    double A[NUMNEIGH*NUMNEIGH], InverseA[NUMNEIGH*NUMNEIGH];
    double X, Y, Alpha, Beta, Sum, Mag, MinMag;
    float RhoSigmaTangent, RhoSigmaNormal;
    int x, y, m, mx, my, n, nx, ny, Weight, R, Offset;
    
    
    if(!SInterp || !(StencilSet = SInterp->StencilSet)
        || S < 0 || S >= GetNumStencils(StencilSet))
        return 0;
    
    RhoSigmaTangent = SInterp->RhoSigmaTangent;
    RhoSigmaNormal = SInterp->RhoSigmaNormal;
    Entry = &SInterp->StencilInterp[S];
    R = (int)ceil(QUAD_RES*GetPsfRadius(StencilSet));
    MinMag = 1e10;
    
    /* Compute orientations for the synthesis functions */
    for(ny = -NEIGHRADIUS, n = 0; ny <= NEIGHRADIUS; ny++)
        for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n++)
        {
            for(y = 0, Beta = 0; y <= QUAD_RES; y++)
            {
                Weight = (y == 0 || y == QUAD_RES) ? 1 : 2;
                X = nx + 0.5;
                Y = ny - 0.5 + ((double)y)/QUAD_RES;
                Beta += Weight*(EvalPhi(StencilSet, S, X, Y)
                    - EvalPhi(StencilSet, S, X - 1, Y));
            }
            
            for(x = 0, Alpha = 0; x <= QUAD_RES; x++)
            {
                Weight = (x == 0 || x == QUAD_RES) ? 1 : 2;
                X = nx - 0.5 + ((double)x)/QUAD_RES;
                Y = ny + 0.5;
                Alpha -= Weight*(EvalPhi(StencilSet, S, X, Y)
                    - EvalPhi(StencilSet, S, X, Y - 1));
            }
            
            Alpha /= (2*QUAD_RES);
            Beta /= (2*QUAD_RES);
            Mag = Alpha*Alpha + Beta*Beta;
            
            if(Mag < MinMag)
                MinMag = Mag;
                                   
            Entry->Orientation[n] = (float)atan2(Beta, Alpha);
                
            while(Entry->Orientation[n] < 0)
                Entry->Orientation[n] += (float)M_PI;
            while(Entry->Orientation[n] >= (float)M_PI)
                Entry->Orientation[n] -= (float)M_PI;
        }
    
    MinMag = sqrt(MinMag);
    Entry->Anisotropy = (float)pow(MinMag,4);
    /* printf("S=%2d: %g\n", S, Entry->Anisotropy);*/
        
     /* Compute interpolation matrix A_{m,n} = (Psf * Rho^n)(m - n) */
    for(my = -NEIGHRADIUS, m = 0; my <= NEIGHRADIUS; my++)
        for(mx = -NEIGHRADIUS; mx <= NEIGHRADIUS; mx++, m++)
        {
            for(ny = -NEIGHRADIUS, n = 0; ny <= NEIGHRADIUS; ny++)
                for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n++)
                {
                    Sum = 0;
                    
                    for(y = -R, Offset = 0; y <= R; y++)
                        for(x = -R; x <= R; x++, Offset++)
                            Sum += PsfSamples[Offset]
                                * Rho(mx - nx - ((double)x)/QUAD_RES,
                                my - ny - ((double)y)/QUAD_RES,
                                SInterp->StencilInterp[S].Orientation[n],
                                RhoSigmaTangent, RhoSigmaNormal,
                                SInterp->StencilInterp[S].Anisotropy);
                                        
                    A[m + NUMNEIGH*n] = Sum/(QUAD_RES*QUAD_RES);
                }
        }
    
    if(!InvertMatrix(InverseA, A, NUMNEIGH))
        return 0;
    
    for(n = 0; n < NUMNEIGH*NUMNEIGH; n++)
        Entry->Matrix[n] = (float)InverseA[n];
    
    return 1;
}


/** @brief Compute samples of the PSF */
static double *SamplePsf(double (*Psf)(double, double, const void*),
    const void *PsfParam, int R, double Step)
{
    const int W = 2*R + 1;
    double *Samples = NULL;
    int i, j, Offset;
    
    
    if(!(Samples = (double *)Malloc(sizeof(double)*W*W)))
        return NULL;
    
    if(!Psf)
        Samples[0] = 1/(Step*Step);
    else
        for(j = -R, Offset = 0; j <= R; j++)
            for(i = -R; i <= R; i++, Offset++)
                Samples[Offset] = Psf(Step*i, Step*j, PsfParam);
    
    return Samples;
}


/**
 * @brief Create a new sinterp object
 * @param StencilSet stencil set \f$\Sigma\f$ to perform interpolation with
 * @param RhoSigmaTangent, RhoSigmaNormal parameters for Rho
 * @return an sinterp*, or NULL on failure
 *
 * Given a stencil set \f$\Sigma\f$, this routine performs the necessary
 * precomputations for interpolation and stores it in an sinterp data
 * structure.  Each stencil \f$\mathcal{S}\in\Sigma\f$ is processed with
 * InterpolateStencil() to create an #sinterpentry.
 */
sinterp *NewSInterp(const sset *StencilSet,
    float RhoSigmaTangent, float RhoSigmaNormal)
{
    sinterp *SInterp = NULL;
    double *PsfSamples = NULL;
    int S;
    
    
    if(!StencilSet
        || !(SInterp = (sinterp *)Malloc(sizeof(struct sinterpstruct))))
        return NULL;
    
    SInterp->StencilSet = StencilSet;
    SInterp->StencilInterp = NULL;
    SInterp->RhoSigmaTangent = RhoSigmaTangent;
    SInterp->RhoSigmaNormal = RhoSigmaNormal;
    
    if(!(SInterp->StencilInterp = (sinterpentry *)Malloc(sizeof(sinterpentry)*
        GetNumStencils(StencilSet)))
        || !(PsfSamples = SamplePsf(GetPsf(StencilSet), GetPsfParam(StencilSet),
        (int)ceil(QUAD_RES*GetPsfRadius(StencilSet)), 1.0/QUAD_RES)))
        goto Catch;
    
    for(S = 0; S < GetNumStencils(StencilSet); S++)
        if(!InterpolateStencil(SInterp, S, PsfSamples))
            goto Catch;
    
    Free(PsfSamples);
    return SInterp;
Catch:
    Free(PsfSamples);
    FreeSInterp(SInterp);
    return NULL;
}


/**
 * @brief Free an sinterp object
 * @param SInterp the sinterp to free
 *
 * Frees all memory associated with SInterp.
 */
void FreeSInterp(sinterp *SInterp)
{
    if(SInterp)
    {
        Free(SInterp->StencilInterp);
        Free(SInterp);
    }
}


/**
* @brief Boundary handling function for constant extension
* @param N is the data length
* @param i is an index into the data
* @return an index that is always between 0 and N - 1
*/
static int ConstExtension(int N, int i)
{
    if(i < 0)
        return 0;
    else if(i >= N)
        return N - 1;
    else
        return i;
}


/**
 * @brief Scale a color image by a possibly non-integer scale factor
 * @param Output the output image
 * @param OutputWidth, OutputHeight dimensions of the output image
 * @param Stencil the best-fitting stencils \f$\mathcal{S}^\star\f$
 * @param SInterp stencil interpolation data
 * @param Input the input image in row-major interleaved order
 * @param InputWidth, InputHeight dimensions of the input image
 * @param ScaleFactor scale factor between input and output images
 * @param CenteredGrid use centered grid if nonzero or top-left otherwise
 * @return 1 on success, 0 on failure
 *
 * This routine implements for a possibly non-integer scale factor the
 * interpolation formula
 * \f[ u(x) = \sum_{k \in \mathbb{Z}^2} w(x-k) \Bigl[ v_k +
 *     \sum_{m\in\mathcal{N}} c_m^k \,
 *     \rho_{\mathcal{S}^\star(k)}^m(x - m - k) \Bigr], \f]
 * where the coefficients \f$c_m^k\f$ are computed as
 * \f[ c_m^k = \sum_{n\in\mathcal{N}} (A_{\mathcal{S}^\star(k)}^{-1})_{m,n}
 *     (v_{k+n} - v_k), \quad m \in \mathcal{N}. \f]
 * For more accurate deconvolution, Prefilter() should be called first to
 * prefilter the input image.
 *
 * In the code, the mathematical meanings of the variables are as follows:
 *  - Output: \f$u\f$
 *  - Input: \f$v\f$
 *  - Coeff[3*(NUMNEIGH + 1)*k + m]: \f$c_m^k\f$
 *  - SInterp->StencilInterp[S].Matrix: \f$A_{\mathcal{S}}^{-1}\f$
 *  - WindowWeight: \f$w(x-k)\f$
 *
 * For computational efficiency, RhoFast() is used to evaluate \f$\rho\f$
 * instead of Rho(), as \f$\rho\f$ needs to be evaluated here many times.
 */
int ArbitraryScale(float *Output, int OutputWidth, int OutputHeight,
    const int *Stencil, const sinterp *SInterp,
    const float *Input, int InputWidth, int InputHeight,
    float ScaleFactor, int CenteredGrid)
{
    int (*Extension)(int, int) = ConstExtension;
    const int InputNumEl = 3*InputWidth*InputHeight;
    float *Coeff = NULL;
    const float *Matrix, *CoeffPtr;
    float u[3], uk[3], v[3], c[3*(NUMNEIGH + 1)], X, Y, Weight, DenomSum;
    float WindowWeightX[2*WINDOWRADIUS], WindowWeightY[2*WINDOWRADIUS];
    float XStart, YStart, Xp, Yp, WindowWeight;
    float RhoSigmaTangent, RhoSigmaNormal;
    int i, ix, iy, k, x, y, m, n, mx, my, nx, ny, S, Success = 0;
    
    
    if(!(Coeff = (float *)Malloc(sizeof(float)*(NUMNEIGH + 1)*InputNumEl)))
        goto Catch;
    
    if(CenteredGrid)
    {
        XStart = (1/ScaleFactor - 1
            + (InputWidth - OutputWidth/ScaleFactor))/2;
        YStart = (1/ScaleFactor - 1
            + (InputHeight - OutputHeight/ScaleFactor))/2;
    }
    else
        XStart = YStart = 0;
    
    RhoSigmaTangent = SInterp->RhoSigmaTangent;
    RhoSigmaNormal = SInterp->RhoSigmaNormal;
    
    /* Compute the coefficients c_m^k for every pixel */
    for(y = 0, k = 0; y < InputHeight; y++)
        for(x = 0; x < InputWidth; x++, k++)
        {
            S = Stencil[x + InputWidth*y];
            Matrix = SInterp->StencilInterp[S].Matrix;
            c[0] = Input[3*k + 0];
            c[1] = Input[3*k + 1];
            c[2] = Input[3*k + 2];
            
            for(m = 3; m < 3*(NUMNEIGH + 1); m++)
                c[m] = 0;
                                    
            for(ny = -NEIGHRADIUS, n = 0; ny <= NEIGHRADIUS; ny++)
                for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n++)
                {
                    if(nx != 0 || ny != 0)
                    {
                        i = 3*(Extension(InputWidth, x + nx)
                            + InputWidth*Extension(InputHeight, y - ny));
                        v[0] = Input[i + 0] - c[0];
                        v[1] = Input[i + 1] - c[1];
                        v[2] = Input[i + 2] - c[2];
                        
                        /* Compute c_m^k */
                        for(m = 0; m < NUMNEIGH; m++)
                        {
                            c[3*(m + 1) + 0] += Matrix[m + NUMNEIGH*n]*v[0];
                            c[3*(m + 1) + 1] += Matrix[m + NUMNEIGH*n]*v[1];
                            c[3*(m + 1) + 2] += Matrix[m + NUMNEIGH*n]*v[2];
                        }
                    }
                }
             
            /* Save the coefficients */
            for(m = 0; m < 3*(NUMNEIGH + 1); m++)
                Coeff[3*(NUMNEIGH + 1)*k + m] = c[m];
        }
    
    /* Compute the interpolation */
    for(y = 0, k = 0; y < OutputHeight; y++)
    {
        Y = YStart + y/ScaleFactor;
        iy = (int)ceil(Y - WINDOWRADIUS);
        
        for(n = 0; n < 2*WINDOWRADIUS; n++)
                WindowWeightY[n] = Window(Y - iy - n);
        
        for(x = 0; x < OutputWidth; x++, k += 3)
        {
            X = XStart + x/ScaleFactor;
            ix = (int)ceil(X - WINDOWRADIUS);
            
            for(n = 0; n < 2*WINDOWRADIUS; n++)
                WindowWeightX[n] = Window(X - ix - n);
            
            DenomSum = 0;
            u[0] = u[1] = u[2] = 0;
            
            for(my = 0; my < 2*WINDOWRADIUS; my++)
            {
                if((iy + my) < 0 || (iy + my) >= InputHeight)
                    continue;
                
                Yp = Y - (iy + my);
                i = ix + InputWidth*(iy + my);
                
                for(mx = 0; mx < 2*WINDOWRADIUS; mx++, i++)
                {
                    if((ix + mx) < 0 || (ix + mx) >= InputWidth)
                        continue;
                    
                    Xp = X - (ix + mx);
                    CoeffPtr = Coeff + i*3*(NUMNEIGH + 1);
                    S = Stencil[i];
                    
                    uk[0] = CoeffPtr[0];
                    uk[1] = CoeffPtr[1];
                    uk[2] = CoeffPtr[2];
                    
                    for(ny = -NEIGHRADIUS, n = 0; ny <= NEIGHRADIUS; ny++)
                        for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n++)
                        {
                            Weight = RhoFast(Xp - nx, -(Yp + ny),
                                SInterp->StencilInterp[S].Orientation[n],
                                RhoSigmaTangent, RhoSigmaNormal,
                                SInterp->StencilInterp[S].Anisotropy);
                            uk[0] += CoeffPtr[3*(n + 1) + 0] * Weight;
                            uk[1] += CoeffPtr[3*(n + 1) + 1] * Weight;
                            uk[2] += CoeffPtr[3*(n + 1) + 2] * Weight;
                        }
                    
                    WindowWeight = WindowWeightX[mx] * WindowWeightY[my];
                    DenomSum += WindowWeight;
                    u[0] += WindowWeight * uk[0];
                    u[1] += WindowWeight * uk[1];
                    u[2] += WindowWeight * uk[2];
                }
            }
                
            Output[k + 0] = u[0] / DenomSum;
            Output[k + 1] = u[1] / DenomSum;
            Output[k + 2] = u[2] / DenomSum;
        }
    }
      
    Success = 1;
Catch:
    Free(Coeff);
    return Success;
}


/**
 * @brief Precompute samples of \f$\Tilde\rho\f$ for an integer scale factor
 * @param SInterp stencil interpolation data
 * @param ScaleFactor the integer scale factor
 * @param CenteredGrid use centered grid if nonzero or top-left otherwise
 * @return pointer to samples, or NULL on failure
 *
 * This routine evaluates samples of
 * \f[ \Tilde{\rho}_{\mathcal{S}}^n(x) = w(x) \sum_{m\in\mathcal{N}}
 *     (A_{\mathcal{S}}^{-1})_{m,n} \, \rho_{\mathcal{S}}^m(x-m). \f]
 * at locations on a fine grid with step size 1/ScaleFactor, where ScaleFactor
 * is an integer.  This precomputation significantly accelerates image scaling
 * by an integer factor.
 *
 * The samples are computed for every stencil \f$\mathcal{S}\in\Sigma\f$ for
 * every neighbor \f$m\in\mathcal{N}\f$ over a square window of size
 * SupportWidth by SupportWidth, where
 *    SupportWidth = 2*WINDOWRADIUS*ScaleFactor
 * or this value minus 1, depending on the choice of grid and the oddness of
 * ScaleFactor.
 */
float *ComputeRhoSamples(const sinterp *SInterp,
    int ScaleFactor, int CenteredGrid)
{
    float *Matrix;
    float *RhoSamples = NULL;
    double Weight, RhoS[NUMNEIGH];
    double x, y, SampleOffset;
    float RhoSigmaTangent, RhoSigmaNormal;
    int SupportWidth, SupportSize;
    int S, sx, sy, i, m, mx, my, n, nx, ny;

    
    if(CenteredGrid && ScaleFactor % 2 == 0)
    {   /* If using the centered grid and the scale factor is even */
        SupportWidth = 2*WINDOWRADIUS*ScaleFactor;
        SampleOffset = 0.5/ScaleFactor - WINDOWRADIUS;
    }
    else
    {   /* If using the top-left grid or scale factor is odd */
        SupportWidth = 2*WINDOWRADIUS*ScaleFactor - 1;
        SampleOffset = 1.0/ScaleFactor - WINDOWRADIUS;
    }
    
    SupportSize = SupportWidth*SupportWidth;
        
    if(!SInterp || !(RhoSamples = (float *)Malloc(sizeof(float)*
        NUMNEIGH*GetNumStencils(SInterp->StencilSet)*SupportSize)))
        return NULL;
        
    RhoSigmaTangent = SInterp->RhoSigmaTangent;
    RhoSigmaNormal = SInterp->RhoSigmaNormal;
    
    /* Compute the samples of the Rho functions */
    for(S = 0; S < GetNumStencils(SInterp->StencilSet); S++)
    {
        Matrix = SInterp->StencilInterp[S].Matrix;
        
        for(sy = SupportWidth - 1, i = 0; sy >= 0; sy--)
            for(sx = 0; sx < SupportWidth; sx++, i++)
            {
                x = SampleOffset + ((double)sx)/ScaleFactor;
                y = SampleOffset + ((double)sy)/ScaleFactor;
                
                for(n = 0; n < NUMNEIGH; n++)
                    RhoS[n] = 0;
                
                RhoS[NEIGHCENTER] = 1;
                
                for(my = -NEIGHRADIUS, m = 0; my <= NEIGHRADIUS; my++)
                    for(mx = -NEIGHRADIUS; mx <= NEIGHRADIUS; mx++, m++)
                    {
                        Weight = Rho(x - mx, y - my,
                            SInterp->StencilInterp[S].Orientation[m],
                            RhoSigmaTangent, RhoSigmaNormal,
                            SInterp->StencilInterp[S].Anisotropy);
                        
                        for(ny = -NEIGHRADIUS, n = 0; ny <= NEIGHRADIUS; ny++)
                            for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n++)
                                if(n != NEIGHCENTER)
                                    RhoS[n] += Matrix[m + NUMNEIGH*n]*Weight;
                    }
                    
                for(n = 0; n < NUMNEIGH; n++)
                    if(n != NEIGHCENTER)
                        RhoS[NEIGHCENTER] -= RhoS[n];
                
                Weight = Window((float)x)*Window((float)y);

                for(n = 0; n < NUMNEIGH; n++)
                    RhoSamples[i + SupportSize*(n + NUMNEIGH*S)]
                        = (float)(RhoS[n] * Weight);
            }
    }
    
    return RhoSamples;
}


/**
 * @brief Scale a color image by an integer factor
 * @param Output the output image
 * @param Stencil the best-fitting stencils \f$\mathcal{S}^\star\f$
 * @param RhoSamples precomputed samples of \f$\Tilde\rho\f$
 * @param Input the input image in row-major interleaved order
 * @param InputWidth, InputHeight dimensions of the input image
 * @param ScaleFactor scale factor between input and output images
 * @param CenteredGrid use centered grid if nonzero or top-left otherwise
 *
 * This routine implements for a integer scale factor the interpolation
 * formula
 * \f[ u(x) = \sum_{k \in \mathbb{Z}^2} \biggl[ w(x-k) v_k +
 *     \sum_{n \in\mathcal{N}\backslash\{0\}} (v_{k+n} - v_k) \,
 *     \Tilde{\rho}_{\mathcal{S}^\star(k)}^n(x-k) \biggr]. \f]
 */
void IntegerScalePass(float *Output,
    const int *Stencil, const float *RhoSamples,
    const float *Input, int InputWidth, int InputHeight,
    int ScaleFactor, int CenteredGrid)
{
    const int OutputWidth = ScaleFactor*InputWidth;
    const int OutputHeight = ScaleFactor*InputHeight;
    int (*Extension)(int, int) = ConstExtension;
    const int InputStride = 3*InputWidth;
    const int OutputStride = 3*OutputWidth;
    const float *RhoSamplePtr, *InputPtr;
    float v[3], Weight;
    int sx, sy, SampleOffset, SupportWidth, SupportSize, ni;
    int ix, iy, InteriorLower, InteriorUpperX, InteriorUpperY;
    int i, i0, i00, x, y, nx, ny, OutputOffset;
    
    
    if(CenteredGrid && ScaleFactor % 2 == 0)
    {   /* If using the centered grid and the scale factor is even */
        SupportWidth = 2*WINDOWRADIUS*ScaleFactor;
        SampleOffset = WINDOWRADIUS*ScaleFactor - ScaleFactor/2;
    }
    else
    {   /* If using the top-left grid or scale factor is odd */
        SupportWidth = 2*WINDOWRADIUS*ScaleFactor - 1;
        SampleOffset = WINDOWRADIUS*ScaleFactor - 1;
        
        if(CenteredGrid)
            SampleOffset -= (ScaleFactor - 1)/2;
    }
    
    SupportSize = SupportWidth*SupportWidth;
    OutputOffset = 3*SampleOffset*(1 + OutputWidth);
    
    /* Determine the interior of the image */
    InteriorLower = (SampleOffset + ScaleFactor - 1)/ScaleFactor;
    InteriorUpperX = (OutputWidth + SampleOffset
        - SupportWidth + 1)/ScaleFactor;
    InteriorUpperY = (OutputHeight + SampleOffset
        - SupportWidth + 1)/ScaleFactor;

    if(InteriorLower < NEIGHRADIUS)
        InteriorLower = NEIGHRADIUS;
    if(InteriorUpperX > InputWidth - NEIGHRADIUS - 1)
        InteriorUpperX = InputWidth - NEIGHRADIUS - 1;
    if(InteriorUpperY > OutputHeight - NEIGHRADIUS - 1)
        InteriorUpperY = OutputHeight - NEIGHRADIUS - 1;
    
    for(iy = -WINDOWRADIUS; iy < InputHeight + WINDOWRADIUS; iy++)
        for(ix = -WINDOWRADIUS; ix < InputWidth + WINDOWRADIUS; ix++)
        {
            if(!(InteriorLower <= ix && ix < InteriorUpperX
                && InteriorLower <= iy && iy < InteriorUpperY))
            {   /* General code for near the boundaries of the image */
                RhoSamplePtr = RhoSamples + SupportSize*NUMNEIGH*
                    Stencil[(Extension(InputWidth, ix)
                    + InputWidth*Extension(InputHeight, iy))];
                
                for(ny = -NEIGHRADIUS; ny <= NEIGHRADIUS; ny++)
                    for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++)
                    {
                        ni = 3*(Extension(InputWidth, ix + nx)
                            + InputWidth*Extension(InputHeight, iy - ny));
                        v[0] = Input[ni + 0];
                        v[1] = Input[ni + 1];
                        v[2] = Input[ni + 2];
                        y = ScaleFactor*iy - SampleOffset;
                        
                        for(sy = 0; sy < SupportWidth;
                            sy++, y++, RhoSamplePtr += SupportWidth)
                            if(0 <= y && y < OutputHeight)
                            {
                                x = ScaleFactor*ix - SampleOffset;
                                i = 3*(x + OutputWidth*y);
                            
                                for(sx = 0; sx < SupportWidth;
                                    sx++, x++, i += 3)
                                    if(0 <= x && x < OutputWidth)
                                    {
                                        Weight = RhoSamplePtr[sx];
                                        Output[i + 0] += Weight*v[0];
                                        Output[i + 1] += Weight*v[1];
                                        Output[i + 2] += Weight*v[2];
                                    }
                            }
                    }
            }
            else
            {   /* Faster code for the interior of the image */
                RhoSamplePtr = RhoSamples +
                    SupportSize*NUMNEIGH*Stencil[ix + InputWidth*iy];
                i00 = 3*(ScaleFactor*ix + OutputWidth*ScaleFactor*iy) - OutputOffset;
                InputPtr = Input + 3*(ix - NEIGHRADIUS + InputWidth*(iy + NEIGHRADIUS));
                
                for(ny = 0; ny < NEIGHDIAMETER; ny++, InputPtr -= InputStride)
                {
                    for(nx = 0; nx < 3*NEIGHDIAMETER; nx += 3)
                    {
                        v[0] = InputPtr[nx + 0];
                        v[1] = InputPtr[nx + 1];
                        v[2] = InputPtr[nx + 2];
                        
                        for(sy = 0, i0 = i00; sy < SupportWidth; sy++,
                            i0 += OutputStride, RhoSamplePtr += SupportWidth)
                            for(sx = 0, i = i0; sx < SupportWidth; sx++, i += 3)
                            {
                                Weight = RhoSamplePtr[sx];
                                Output[i + 0] += Weight*v[0];
                                Output[i + 1] += Weight*v[1];
                                Output[i + 2] += Weight*v[2];
                            }
                    }
                }
            }
        }
}


/**
 * @brief Prefilter a color image to improve deconvolution accuracy
 * @param Filtered the prefiltered image
 * @param Stencil the best-fitting stencils \f$\mathcal{S}^\star\f$
 * @param RhoSamples precomputed samples of \f$\Tilde\rho\f$
 * @param Input the input image in row-major interleaved order
 * @param InputWidth, InputHeight dimensions of the input image
 * @param ScaleFactor scale factor between input and output images
 * @param CenteredGrid use centered grid if nonzero or top-left otherwise
 * @param PsfSigma Gaussian PSF standard deviation in units of input pixels
 * @param NumPasses number of passes of iterative refinement to perform
 *
 * This routine performs the iterative refinement prefiltering
 * \f[ v^0 = v, \quad v^{i+1} = v^i
 *     + \bigl(v - \operatorname{sample}(h*\mathcal{R}v^i)\bigr). \f]
 * Interpolation with a prefiltered image improves the resulting
 * deconvolution accuracy.  Interpolation is performed using the prefiltered
 * image as input with either IntegerScalePass() for an integer scale factor or
 * by ArbitraryScale() for a non-integer scale factor.
 *
 * Prefiltering is implemented with IntegerScalePass() to perform the
 * interpolations \f$\mathcal{R}v^i\f$ and DeconvResidual() to compute the
 * deconvolution residuals \f$v - \operatorname{sample}(h*\mathcal{R}v^i)\f$.
 */
int Prefilter(float *Filtered,
    const int *Stencil, const float *RhoSamples,
    const float *Input, int InputWidth, int InputHeight,
    int ScaleFactor, int CenteredGrid, float PsfSigma, int NumPasses)
{
    const int InputNumEl = 3*InputWidth*InputHeight;
    const int InterpWidth = ScaleFactor*InputWidth;
    const int InterpHeight = ScaleFactor*InputHeight;
    const int InterpNumEl = 3*InterpWidth*InterpHeight;
    float *Interp = NULL, *Residual = NULL;
    float MaxResidual;
    int i, Pass, Success = 0;
    
    
    if(!Filtered || !Stencil || !RhoSamples || !Input
        || !(Interp = (float *)Malloc(sizeof(float)*InterpNumEl))
        || !(Residual = (float *)Malloc(sizeof(float)*InputNumEl)))
        goto Catch;
    
    memcpy(Filtered, Input, sizeof(float)*InputNumEl);
    
    for(i = 0; i < InterpNumEl; i++)
        Interp[i] = 0;
    
    printf("\n  Iteration   Residual norm\n  -------------------------\n");
    
    for(Pass = 0; Pass < NumPasses; Pass++)
    {
        if(Pass == 0)
            IntegerScalePass(Interp, Stencil, RhoSamples,
                Input, InputWidth, InputHeight, ScaleFactor, CenteredGrid);
        else
            IntegerScalePass(Interp, Stencil, RhoSamples,
                Residual, InputWidth, InputHeight, ScaleFactor, CenteredGrid);
        
        if((MaxResidual = DeconvResidual(Residual,
            Input, InputWidth, InputHeight,
            Interp, InterpWidth, InterpHeight,
            (float)ScaleFactor, CenteredGrid, PsfSigma)) < 0)
            goto Catch;
        
        printf("  %8d %15.8f\n", Pass + 1, MaxResidual);
        
        for(i = 0; i < InputNumEl; i++)
            Filtered[i] += Residual[i];
    }
    
    Success = 1;
Catch:
    Free(Residual);
    Free(Interp);
    return Success;
}


/** @brief Evaluate x^2 */
static float Sqr(float x)
{
    return x*x;
}


/**
 * @brief Compute deconvolution residual for prefiltering
 * @param Residual the residual image
 * @param Coarse the coarse image
 * @param CoarseWidth, CoarseHeight dimensions of the coarse image
 * @param Interp the interpolated image
 * @param InterpWidth, InterpHeight dimensions of the interpolated image
 * @param ScaleFactor scale factor between coarse and interpolated images
 * @param CenteredGrid use centered grid if nonzero or top-left otherwise
 * @param PsfSigma Gaussian PSF standard deviation in units of input pixels
 *
 * This routine computes the deconvolution residual
 * \f[ r = v - \operatorname{sample}(h * u), \f]
 * where v is the coarse image and u is the interpolation.  This residual is
 * used in the iterative refinement procedure implemented in Prefilter().
 */
float DeconvResidual(float *Residual,
    const float *Coarse, int CoarseWidth, int CoarseHeight,
    const float *Interp, int InterpWidth, int InterpHeight,
    float ScaleFactor, int CenteredGrid, float PsfSigma)
{
    /*int (*Extension)(int, int) = ExtensionMethod[Param.Boundary];*/
    int (*Extension)(int, int) = ConstExtension;
    const int CoarseStride = 3*CoarseWidth;
    const float PsfRadius = 4*PsfSigma*ScaleFactor;
    const int PsfWidth = (int)ceil(2*PsfRadius);
    float *Temp = NULL, *PsfBuf = NULL;
    float ExpDenom, Weight, Sum[3], DenomSum, MaxResidual = 0;
    float XStart, YStart, X, Y;
    int IndexX0, IndexY0, SrcOffset, DestOffset;
    int x, y, i, n, c;
    
    
    if(!(Temp = (float *)Malloc(sizeof(float)*3*CoarseWidth*InterpHeight))
        || !(PsfBuf = (float *)Malloc(sizeof(float)*PsfWidth)))
    {
        free(Temp);
        free(PsfBuf);
        return -1;
    }
    
    if(CenteredGrid)
    {
        XStart = (1/ScaleFactor - 1
            + (CoarseWidth - InterpWidth/ScaleFactor))/2;
        YStart = (1/ScaleFactor - 1
            + (CoarseHeight - InterpHeight/ScaleFactor))/2;
    }
    else
        XStart = YStart = 0;
    
    ExpDenom = 2 * Sqr(PsfSigma*ScaleFactor);
    
    for(x = 0; x < CoarseWidth; x++)
    {
        X = (-XStart + x)*ScaleFactor;
        IndexX0 = (int)ceil(X - PsfRadius);
        
        /* Evaluate the PSF */
        for(n = 0; n < PsfWidth; n++)
            PsfBuf[n] = (float)exp(-Sqr(X - (IndexX0 + n)) / ExpDenom);
        
        for(y = 0, SrcOffset = 0, DestOffset = 3*x; y < InterpHeight;
            y++, SrcOffset += InterpWidth, DestOffset += CoarseStride)
        {
            Sum[0] = Sum[1] = Sum[2] = DenomSum = 0;
            
            for(n = 0; n < PsfWidth; n++)
            {
                Weight = PsfBuf[n];
                DenomSum += Weight;
                i = 3*(Extension(InterpWidth, IndexX0 + n) + SrcOffset);
                
                for(c = 0; c < 3; c++)
                    Sum[c] += Weight * Interp[i + c];
            }
            
            for(c = 0; c < 3; c++)
                Temp[DestOffset + c] = Sum[c] / DenomSum;
        }
    }
    
    for(y = 0; y < CoarseHeight; y++,
        Residual += CoarseStride, Coarse += CoarseStride)
    {
        Y = (-YStart + y)*ScaleFactor;
        IndexY0 = (int)ceil(Y - PsfRadius);
            
        /* Evaluate the PSF */
        for(n = 0; n < PsfWidth; n++)
            PsfBuf[n] = (float)exp(-Sqr(Y - (IndexY0 + n)) / ExpDenom);
        
        for(x = 0; x < CoarseStride; x += 3)
        {
            Sum[0] = Sum[1] = Sum[2] = DenomSum = 0;
            
            for(n = 0; n < PsfWidth; n++)
            {
                SrcOffset = x + CoarseStride*Extension(InterpHeight, IndexY0 + n);
                Weight = PsfBuf[n];
                DenomSum += Weight;
                
                for(c = 0; c < 3; c++)
                    Sum[c] += Weight * Temp[SrcOffset + c];
            }
            
            for(c = 0; c < 3; c++)
            {
                Residual[x + c] = Coarse[x + c] - Sum[c] / DenomSum;
            
                if(fabs(Residual[x + c]) > MaxResidual)
                    MaxResidual = (float)fabs(Residual[x + c]);
            }
        }
    }
    
    Free(PsfBuf);
    Free(Temp);
    return MaxResidual;
}

