/**
 * @file cwinterp.c
 * @brief Contour stencil windowed interpolation
 * @author Pascal Getreuer <getreuer@gmail.com>
 *
 * This file implements contour stencil windowed interpolation for integer
 * scale factors.  The interpolation model supposes that the input image
 * \f$v\f$ was created by convolution followed by downsampling
 * \f[ v = {\downarrow_r} (h * u) \f]
 * where \f$u\f$ is the underlying high resolution image, \f$h\f$ is the point-
 * spread function (PSF), and \f$\downarrow_r\f$ denotes downsampling by factor
 * \f$r\f$.  For simplicity, this code requires that \f$r\f$ is integer and
 * \f$h\f$ is a Gaussian.  The standard deviation \f$h\f$ is controlled by
 * \c PsfSigma.
 *
 * The image is interpolated by the following steps.  First, contour stencils
 * are applied to estimate the local contour orientations (in routine
 * \c FitStencils).  The image is then interpolated by the formula
 * \f[ u(x) = \sum_{k\in\mathbb{Z}^2} w(x - k) \Bigl[ v_k +
 *     \sum_{n\in\mathcal{N}\backslash\{0\}} (v_{k+n} - v_k)
 *     \psi^n_{\mathcal{S}^\star(k)}(x - k) \Bigr] \f]
 * (routine \c CWFirstPass).  This initial interpolation approximately
 *  satisfies the degradation model, \f$v \approx {\downarrow_r} (h * u)\f$.  To
 * improve the accuracy, several correction passes are applied.  In each pass,
 * the residual is computed (routine \c CWResidual) and then applied to refine
 * the interpolation (routine \c CWRefinementPass).  Under typical settings,
 * the degradation model is accurately satisfied after two or three correction
 * passes.
 *
 * The main computations \c CWFirstPass and \c CWRefinementPass are done in
 * fixed-point integer arithmetic.  The downsides of using fixed-point
 * arithmetic compared to floating-point arithmetic are more involved code for
 * multiplication and division (it is harder to read) and the range and
 * precision of fixed-point integers must be explicitly managed for accurate
 * results (need to be careful).  The upside is that when it does work, fixed-
 * point arithmetic is significantly faster.  Fortunately, in this application,
 * the only needed operations are fixed-point additions and fixed-point
 * multiplies where both factors are in a predictable range of values.  These
 * conditions are good for fixed-point arithmetic.
 *
 * The number of fractional bits used to represent \f$v\f$ and the residual is
 * controlled by \c INPUT_FRACBITS.  The number of fractional bits in
 * representing \f$\psi\f$ is \c PSI_FRACBITS.  Since \f$v\f$ and \f$\psi\f$
 * are multiplied, the output has
 *    \c OUTPUT_FRACBITS = (\c INPUT_FRACBITS + \c PSI_FRACBITS)
 * fractional bits.  These constants must be balanced between precision and
 * avoiding overflow.  Fortunately this balance is easy to find for this
 * application (neither extreme range nor extreme precision are needed).
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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ipol/basic.h>
#include "cwinterp.h"
#include "fitsten.h"
#include "invmat.h"
#include "nninterp.h"
#include "drawline.h"


/** @brief The number of contour stencils, cardinality of \f$\Sigma\f$ */
#define NUMSTENCILS         8

/** @brief Cardinality of the neighborhood \f$\mathcal{N}\f$ */
#define NEIGHRADIUS     1
#define NEIGHDIAMETER   (2*NEIGHRADIUS+1)
#define NUMNEIGH        (NEIGHDIAMETER*NEIGHDIAMETER)

/**
* @brief \c CWRefinementPass residual tolerance.
*
* In the correction passes, pixels will be skipped if they have magnitude
* less than 2^(-INPUT_FRACBITS + CORRECTION_IGNOREBITS), which improves the
* speed.  A larger value of CORRECTION_IGNOREBITS makes the correction passes
* faster, but less accurate.
*/
#define CORRECTION_IGNOREBITS	3

/** @brief Number of fractional bits in the input array */
#define INPUT_FRACBITS		8
/** @brief Number of fractional bits in the psi sample arrays */
#define PSI_FRACBITS		12
/** @brief Number of fractional bits in the output array */
#define OUTPUT_FRACBITS		(INPUT_FRACBITS + PSI_FRACBITS)

/** @brief The number 1.0 in fixed-point with N fractional bits */
#define FIXED_ONE(N)	(1 << (N))
/** @brief The number 0.5 in fixed-point with N fractional bits */
#define FIXED_HALF(N)	(1 << ((N) - 1))


/**
* @brief Number of elements between successive RGB fixed-point pixels
*
* \c PIXEL_STRIDE is the number of elements between successive pixels in the
* RGB fixed-point representation used internally by the interpolation.
*
* @note Unlike the other defines here, this one cannot be changed without
* also rewriting much of the code.  Its purpose is more for clarity rather
* than working as a parameter.
*/
#define PIXEL_STRIDE    3


/* Generic macros */

/** @brief Clamp X to [A, B] */
#define CLAMP(X,A,B)    (((X) < (A)) ? (A) : (((X) > (B)) ? (B) : (X)))

/** @brief Round and clamp double X to integer */
#define ROUNDCLAMP(X,A,B) (((X) < (A)) ? (A) : (((X) > (B)) ? (B) : ROUND(X)))

#ifndef M_2PI
/** @brief The constant 2*pi */
#define M_2PI       6.283185307179586476925286766559
#endif


/* Orientation in radians of each stencil */
static const double StencilOrientation[NUMSTENCILS] = {-1.178097245096172,
    -0.785398163397448, -0.392699081698724, 0.0, 0.392699081698724,
    0.785398163397448, 1.178097245096172, 1.570796326794897};
    
    
/** @brief The point spread function (PSF) */
static double Psf(double x, double y, cwparams Param)
{
    double SigmaSqr = Param.PsfSigma*Param.PsfSigma;
    return exp(-(x*x + y*y)/(2.0*SigmaSqr))/(M_2PI*SigmaSqr);
}


/** @brief The oriented functions phi used in local reconstructions */
static double Phi(double x, double y,
    double theta, double PhiSigmaTangent, double PhiSigmaNormal)
{
    double t, n;
    
    /* Oriented Gaussian */
    t = (cos(theta)*x + sin(theta)*y) / PhiSigmaTangent;
    n = (-sin(theta)*x + cos(theta)*y) / PhiSigmaNormal;
    
    return exp(-0.5*(t*t + n*n));
}


/** @brief Cubic B-spline */
static float CubicBSpline(float x)
{
    x = (float)fabs(x);

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


/** @brief The window used to sum the global solution */
static double Window(double x, double y)
{
    double Temp;
    
    x *= 2.0/(NEIGHRADIUS + 1.0);
    y *= 2.0/(NEIGHRADIUS + 1.0);
    
    /* Cubic B-spline */
    if(-2.0 < x && x < 2.0 && -2.0 < y && y < 2.0)
    {
        x = fabs(x);
        Temp = fabs(1.0 - x);
        x = 1.0 - x + (x*x*x - 2.0*Temp*Temp*Temp)/6.0;
        
        y = fabs(y);
        Temp = fabs(1.0 - y);
        y = 1.0 - y + (y*y*y - 2.0*Temp*Temp*Temp)/6.0;
        return (x * y) * 4.0/((NEIGHRADIUS + 1.0)*(NEIGHRADIUS + 1.0));
    }
    else
        return 0.0;
}


/** @brief Quadrature weights and abscissas for composite Gauss-Lobatto */
static void QuadraturePoint(double *Weight, double *Abscissa, int Index, int NumPanels)
{
    switch(Index % 3)
    {
    case 0:
        *Weight = (Index == 0 || Index == NumPanels)? 0.25 : 0.5;
        *Abscissa = Index;
        break;
    case 1:
        *Weight = 1.25;
        /* Abscissa location is Index - (3.0/sqrt(5.0) - 1.0)/2.0 */
        *Abscissa = Index - 1.70820393249936919e-1;
        break;
    case 2:
        *Weight = 1.25;
        /* Abscissa location is Index + (3.0/sqrt(5.0) - 1.0)/2.0 */
        *Abscissa = Index + 1.70820393249936919e-1;
        break;
    }
}


/** @brief Compute the convolution of the PSF and phi_theta at the point (x,y) */
static double PsfPhiConvolution(int x, int y, double Theta, cwparams Param)
{
    /* Integrate over the square [-R,R]x[-R,R] */
    const double R = 4.0*Param.PsfSigma;
    /* Number of panels along each dimension, must be divisible by 3 */
    const int NumPanels = 3*16;
    const double PanelSize = 2.0*R/NumPanels;
    double Integral = 0.0, Slice = 0.0;
    double u, v, wu, wv;
    int IndexX, IndexY;
    
    
    /* Specially handle the case where PSF is Dirac delta */
    if(Param.PsfSigma == 0.0)
        return Phi(x, y, Theta, Param.PhiSigmaTangent, Param.PhiSigmaNormal);
    
    /* Approximate 2D integral */
    for(IndexY = 0; IndexY <= NumPanels; IndexY++)
    {
        QuadraturePoint(&wv, &v, IndexY, NumPanels);
        v = PanelSize*v - R;
        
        for(Slice = 0.0, IndexX = 0; IndexX <= NumPanels; IndexX++)
        {
            QuadraturePoint(&wu, &u, IndexX, NumPanels);
            u = PanelSize*u - R;
            Slice += wu*( Psf(u, v, Param) *
                Phi(x - u, y - v, Theta, Param.PhiSigmaTangent,
                Param.PhiSigmaNormal) );
        }
        
        Integral += wv*Slice;
    }
    
    Integral *= PanelSize*PanelSize;
    return Integral;
}


/** @brief Compute the deconvolution matrices (A_S)^-1 */
static int ComputeMatrices(double *InverseA, cwparams Param)
{
    double *A = NULL;
    int m, mx, my, n, nx, ny, S;
    int Status = 0;
    

    if(!(A = (double *)Malloc(sizeof(double)*NUMNEIGH*NUMNEIGH)))
        goto Catch;

    for(S = 0; S < NUMSTENCILS; S++)
    {
        for(ny = -NEIGHRADIUS, n = 0; ny <= NEIGHRADIUS; ny++)
        for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n++)
            for(my = -NEIGHRADIUS, m = 0; my <= NEIGHRADIUS; my++)
            for(mx = -NEIGHRADIUS; mx <= NEIGHRADIUS; mx++, m++)
            {
                A[m + NUMNEIGH*n] = PsfPhiConvolution(mx - nx, my - ny,
                    StencilOrientation[S], Param);
            }

        /* Compute the inverse of A */
        if(!InvertMatrix(InverseA + S*(NUMNEIGH*NUMNEIGH), A, NUMNEIGH))
            goto Catch;
    }
    
    Status = 1;
Catch:  /* This label is used for error handling.  If something went wrong
        above (which may be out of memory or a computation error), then
        execution jumps to this point to clean up and exit. */
    Free(A);
    return Status;
}


#include <ipol/imageio.h>

/**
* @brief Precomputations before windowed interpolation \c CWInterp
*
* @param Param cwparams struct of interpolation parameters
*
* @return Pointer to \f$\psi\f$ samples array, or null on failure
*
* \c PreCWInterp precomputes samples of the \f$\psi\f$ functions,
* \f[ \psi^n_\mathcal{S}(x) = \sum_{m\in\mathcal{N}}
*        (A_\mathcal{S})^{-1}_{m,n} \varphi^m_\mathcal{S}(x - m). \f]
* The routine allocates memory to store the samples and returns a pointer to
* this memory.  It is the responsibility of the caller to call \c free
* on this pointer when done to release the memory.
*
* A non-null pointer indicates success.  On failure, the returned pointer
* is null.
*/
int32_t *PreCWInterp(cwparams Param)
{
    int32_t *Psi = NULL;
    const int ScaleFactor = (int)ceil(Param.ScaleFactor);
    const int SupportRadius = (NEIGHRADIUS+1)*ScaleFactor - 1;
    const int SupportWidth = 2*SupportRadius + 1;
    const int SupportSize = SupportWidth*SupportWidth;
    double *InverseA = NULL;
    double x, y, Wxy, Psi0, Psim, XStart, YStart, WSum;
    int S, sx, sy, i, m0, m, mx, my, n, nx, ny, Success = 0;
    
    
    if(!(Psi = (int32_t *)Malloc(sizeof(int32_t)*SupportSize*NUMNEIGH*NUMSTENCILS))
        || !(InverseA = (double *)Malloc(sizeof(double)*NUMNEIGH*NUMNEIGH*NUMSTENCILS)))
        goto Catch;
    
    /* Compute the matrices, the results are stored in InverseA. */
    if(!ComputeMatrices(InverseA, Param))
        goto Catch;
    
    if(Param.CenteredGrid)
    {
        XStart = (1/Param.ScaleFactor - 1)/2;
        YStart = (1/Param.ScaleFactor - 1)/2;
    }
    else
        XStart = YStart = 0;
    
    m0 = NEIGHRADIUS + NEIGHRADIUS*NEIGHDIAMETER;
    
    /* Precompute the samples of the Psi functions */
    for(S = 0; S < NUMSTENCILS; S++)
        for(i = 0, sy = -SupportRadius; sy <= SupportRadius; sy++)
            for(sx = -SupportRadius; sx <= SupportRadius; sx++, i++)
            {
                /* Compute the sum of window translates.  This sum should be
                   exactly constant, but there can be small variations.  We
                   divide the Psi samples computed below by WSum to compensate.
                 */
                for(ny = -(int)floor((sy + SupportRadius)/ScaleFactor), WSum = 0;
                    ny <= (2*NEIGHRADIUS + 1) && sy + ny*ScaleFactor <= SupportRadius; ny++)
                    for(nx = -(int)floor((sx + SupportRadius)/ScaleFactor);
                        nx <= (2*NEIGHRADIUS + 1) && sx + nx*ScaleFactor <= SupportRadius; nx++)
                    {
                        x = XStart + nx + ((double)sx)/((double)ScaleFactor);
                        y = YStart + ny + ((double)sy)/((double)ScaleFactor);
                        WSum += ROUND(Window(x,y)*FIXED_ONE(PSI_FRACBITS))
                            / (double)FIXED_ONE(PSI_FRACBITS);
                    }
                                
                x = XStart + ((double)sx)/((double)ScaleFactor);
                y = YStart + ((double)sy)/((double)ScaleFactor);
                Psi0 = Wxy = Window(x, y);
                
                for(my = -NEIGHRADIUS, m = 0; my <= NEIGHRADIUS; my++)
                for(mx = -NEIGHRADIUS; mx <= NEIGHRADIUS; mx++, m++)
                {
                    if(m != m0)
                    {
                        Psim = 0.0;
                        
                        for(ny = -NEIGHRADIUS, n = 0; ny <= NEIGHRADIUS; ny++)
                        for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n++)
                            Psim += InverseA[m + NUMNEIGH*(n + NUMNEIGH*S)]
                                * Phi(x - nx, y - ny,
                                StencilOrientation[S], Param.PhiSigmaTangent,
                                Param.PhiSigmaNormal);

                        Psim *= Wxy;
                        Psi0 -= Psim;
                        
                        Psi[i + SupportSize*(m + NUMNEIGH*S)] =
                            (int32_t)ROUND((Psim/WSum)*FIXED_ONE(PSI_FRACBITS));
                    }
                }
                
                Psi[i + SupportSize*(m0 + NUMNEIGH*S)] =
                    (int32_t)ROUND((Psi0/WSum)*FIXED_ONE(PSI_FRACBITS));
            }
    
    Success = 1;
Catch:  /* This label is used for error handling.  If something went wrong
        above (which may be out of memory or a computation error), then
        execution jumps to this point to clean up and exit. */
    Free(InverseA);
    if(!Success && Psi)
    {
        Free(Psi);
        Psi = NULL;
    }
    return Psi;
}


/**
* @brief Main interpolation computation for the first pass
* @param Interpolation pointer where to store the result
* @param ScaleFactor the interpolation scale factor
* @param Input pointer to the input image
* @param InputWidth, InputHeight dimensions of the input image
* @param Stencil pointer to the selected stencils
* @param Sample array of precomputed w*psi samples
*
* This is the main computation for refinement passes: it adds to the
* interpolation
* \f[ u(x) = \sum_{k\in\mathbb{Z}^2} w(x - k) \Bigl[ v_k +
*     \sum_{n\in\mathcal{N}\backslash\{0\}} (v_{k+n} - v_k)
*     \psi^n_{\mathcal{S}^\star(k)}(x - k) \Bigr]. \f]
*/
static void CWFirstPass(int32_t *Interpolation, int ScaleFactor, const int32_t *Input,
    int InputWidth, int InputHeight, const int *Stencil, const int32_t *Psi)
{
    const int32_t *PsiPtr, *SrcWindow;
    int32_t *DestWindow;
    
    const int Pad = 2;
    const int SampleRange = (NEIGHRADIUS+1)*ScaleFactor - 1;
    const int SampleWidth = 2*SampleRange + 1;
    
    const int OutputWidth = ScaleFactor*InputWidth;
    const int DestWindowJump = PIXEL_STRIDE*(OutputWidth - SampleWidth);
    const int DestStep = PIXEL_STRIDE*ScaleFactor;
    const int DestJump = PIXEL_STRIDE*(ScaleFactor-1)*OutputWidth + 2*Pad*DestStep;
    const int SrcWindowJump = PIXEL_STRIDE*(InputWidth - NEIGHDIAMETER);
    const int SrcJump = 2*PIXEL_STRIDE*Pad;
    const int StencilJump = 2*Pad;
    
    int x, y, NeighX, NeighY, SampleX, SampleY;
    int32_t cr, cg, cb;

    
    Interpolation += PIXEL_STRIDE*(Pad*ScaleFactor - SampleRange)*(1 + OutputWidth);
    Input += PIXEL_STRIDE*(Pad - NEIGHRADIUS)*(1 + InputWidth);
    Stencil += Pad*(1 + InputWidth);
    
    for(y = InputHeight - 2*Pad; y; y--,
        Stencil += StencilJump, Input += SrcJump, Interpolation += DestJump)
    for(x = InputWidth - 2*Pad; x; x--,
        Stencil++, Input += PIXEL_STRIDE, Interpolation += DestStep)
    {
        PsiPtr = Psi + *Stencil;
        SrcWindow = Input;
        
        for(NeighY = NEIGHDIAMETER; NeighY; NeighY--, SrcWindow += SrcWindowJump)
        for(NeighX = NEIGHDIAMETER; NeighX; NeighX--, SrcWindow += PIXEL_STRIDE)
        {
            cr = SrcWindow[0];
            cg = SrcWindow[1];
            cb = SrcWindow[2];
            DestWindow = Interpolation;
            SampleY = SampleWidth;
            
            for(SampleY = SampleWidth; SampleY;
                SampleY--, DestWindow += DestWindowJump, PsiPtr += SampleWidth)
            for(SampleX = 0; SampleX < SampleWidth;
                SampleX++, DestWindow += PIXEL_STRIDE)
            {
                int32_t Temp = PsiPtr[SampleX];
                DestWindow[0] += cr * Temp;
                DestWindow[1] += cg * Temp;
                DestWindow[2] += cb * Temp;
            }
        }
    }
}


/**
* @brief Main interpolation computation for refinement passes
* @param Interpolation pointer to where interpolation is stored
* @param ScaleFactor the interpolation scale factor
* @param Residual pointer to the residual
* @param InputWidth, InputHeight dimensions of the input image
* @param Stencil pointer to the selected stencils
* @param Sample array of precomputed w*psi samples
*
* This is the main computation for refinement passes: it adds to the
* interpolation
* \f[ u(x) = u(x) + \sum_{k\in\mathbb{Z}^2} w(x - k) \Bigl[ r_k +
*     \sum_{n\in\mathcal{N}\backslash\{0\}} (r_{k+n} - r_k)
*     \psi^n_{\mathcal{S}^\star(k)}(x - k) \Bigr]. \f]
*/
static void CWRefinementPass(int32_t *Interpolation, int ScaleFactor,
    const int32_t *Residual, int InputWidth, int InputHeight,
    const int *Stencil, const int32_t *Sample)
{
    const int Pad = 4;

    const int32_t *SamplePtr, *SrcWindow;
    int32_t *DestWindow;
    const int SampleRange = (NEIGHRADIUS+1)*ScaleFactor - 1;
    const int SampleWidth = 2*SampleRange + 1;
    const int SampleSize = SampleWidth*SampleWidth;
    
    const int OutputWidth = ScaleFactor*InputWidth;
    const int DestWindowJump = PIXEL_STRIDE*(OutputWidth - SampleWidth);
    const int DestStep = PIXEL_STRIDE*ScaleFactor;
    const int DestJump = PIXEL_STRIDE*(ScaleFactor-1)*OutputWidth + 2*Pad*DestStep;
    const int SrcWindowJump = PIXEL_STRIDE*(InputWidth - NEIGHDIAMETER);
    const int SrcJump = 2*PIXEL_STRIDE*Pad;
    const int StencilJump = 2*Pad;
    
    int x, y, NeighX, NeighY, SampleX, SampleY;
    int32_t cr, cg, cb;

    
    Interpolation += PIXEL_STRIDE*(Pad*ScaleFactor - SampleRange)*(1 + OutputWidth);
    Residual += PIXEL_STRIDE*(Pad - NEIGHRADIUS)*(1 + InputWidth);
    Stencil += Pad*(1 + InputWidth);

    for(y = InputHeight - 2*Pad; y; y--,
        Stencil += StencilJump, Residual += SrcJump, Interpolation += DestJump)
    for(x = InputWidth - 2*Pad; x; x--,
        Stencil++, Residual += PIXEL_STRIDE, Interpolation += DestStep)
    {
        SamplePtr = Sample + *Stencil;
        SrcWindow = Residual;
        
        for(NeighY = NEIGHDIAMETER; NeighY; NeighY--, SrcWindow += SrcWindowJump)
        for(NeighX = NEIGHDIAMETER; NeighX; NeighX--, SrcWindow += PIXEL_STRIDE)
        {
            cr = SrcWindow[0];
            cg = SrcWindow[1];
            cb = SrcWindow[2];
            
            if( ((cr >> CORRECTION_IGNOREBITS) && (-cr >> CORRECTION_IGNOREBITS))
                || ((cg >> CORRECTION_IGNOREBITS) && (-cg >> CORRECTION_IGNOREBITS))
                || ((cb >> CORRECTION_IGNOREBITS) && (-cb >> CORRECTION_IGNOREBITS)) )
            {
                DestWindow = Interpolation;
                SampleY = SampleWidth;
                
                for(SampleY = SampleWidth; SampleY;
                    SampleY--, DestWindow += DestWindowJump, SamplePtr += SampleWidth)
                for(SampleX = 0; SampleX < SampleWidth;
                    SampleX++, DestWindow += PIXEL_STRIDE)
                {
                    int32_t Temp = SamplePtr[SampleX];
                    DestWindow[0] += cr * Temp;
                    DestWindow[1] += cg * Temp;
                    DestWindow[2] += cb * Temp;
                }
            }
            else
                SamplePtr += SampleSize;
        }
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


static float Sqr(float x)
{
    return x*x;
}


/** @brief Computes the residual, Residual = Input - sample(PSF * Interpolation) */
static int32_t CWResidual(int32_t *Residual, const int32_t *Interpolation,
        const int32_t *Input, int CoarseWidth, int CoarseHeight, cwparams Param)
{
    int Pad = 4;
    const int ScaleFactor = (int)ceil(Param.ScaleFactor);
    const int InterpWidth = ScaleFactor*CoarseWidth;
    const int InterpHeight = ScaleFactor*CoarseHeight;
    const int CoarseStride = 3*CoarseWidth;
    const float PsfRadius = (float)(4*Param.PsfSigma*ScaleFactor);
    const int PsfWidth = (int)ceil(2*PsfRadius);
    float *Temp = NULL, *PsfBuf = NULL;
    float ExpDenom, Weight, Sum[3], DenomSum;
    float XStart, YStart, X, Y;
    int IndexX0, IndexY0, SrcOffset, DestOffset;
    int x, y, i, n, c, Success = 0;
    int x1, x2;
    int32_t ResNorm = 0;
    
    
    if(!(Temp = (float *)Malloc(sizeof(float)*3*CoarseWidth*InterpHeight))
        || !(PsfBuf = (float *)Malloc(sizeof(float)*PsfWidth)))
        goto Catch;
    
    if(Param.CenteredGrid)
    {
        XStart = (1.0f/ScaleFactor - 1)/2;
        YStart = (1.0f/ScaleFactor - 1)/2;
    }
    else
        XStart = YStart = 0;
    
    if(Param.PsfSigma)
        ExpDenom = 2 * Sqr((float)(Param.PsfSigma*ScaleFactor));
    else
        ExpDenom = 2 * Sqr(1e-2f*ScaleFactor);

    if(Pad < ScaleFactor)
        Pad = ScaleFactor;
    
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
                i = 3*(ConstExtension(InterpWidth, IndexX0 + n) + SrcOffset);
                
                for(c = 0; c < 3; c++)
                    Sum[c] += Weight * Interpolation[i + c];
            }
            
            for(c = 0; c < 3; c++)
                Temp[DestOffset + c] = Sum[c] / DenomSum;
        }
    }
    
    x1 = 3*Pad;
    x2 = CoarseStride - 3*Pad;
    
    for(y = 0; y < CoarseHeight; y++,
        Residual += CoarseStride, Input += CoarseStride)
    {
        if(!(y >= Pad && y < CoarseHeight-Pad))
            continue;
        
        Y = (-YStart + y)*ScaleFactor;
        IndexY0 = (int)ceil(Y - PsfRadius);
            
        /* Evaluate the PSF */
        for(n = 0; n < PsfWidth; n++)
            PsfBuf[n] = (float)exp(-Sqr(Y - (IndexY0 + n)) / ExpDenom);
                
        for(x = x1; x < x2; x += 3)
        {
            Sum[0] = Sum[1] = Sum[2] = DenomSum = 0;
            
            for(n = 0; n < PsfWidth; n++)
            {
                SrcOffset = x + CoarseStride*ConstExtension(InterpHeight, IndexY0 + n);
                Weight = PsfBuf[n];
                DenomSum += Weight;
                
                for(c = 0; c < 3; c++)
                    Sum[c] += Weight * Temp[SrcOffset + c];
            }
            
            DenomSum *= FIXED_ONE(PSI_FRACBITS);
            
            for(c = 0; c < 3; c++)
            {
                Sum[c] = Input[x + c] - Sum[c] / DenomSum;
                Residual[x + c] = (int32_t)ROUND(Sum[c]);
                                
                if(abs(Residual[x + c]) > ResNorm)
                    ResNorm = abs(Residual[x + c]);
            }
        }
    }
    
    Success = 1;
Catch:
    Free(PsfBuf);
    Free(Temp);
    return (Success) ? ResNorm : -1;
}


/**
* @brief Simultaneously pad and convert to RGB fixed-point image
*
* @param FixedRgb pointer to hold the padded and converted image data
* @param Input the input 32-bit RGBA image
* @param InputWidth, InputHeight input image dimensions
* @param Padding number of padding pixels
*
* \c ConvertInput is used by \c CWInterp to prepare the input image in a
* format convenient for computations.
*
* \c ConvertInput converts the input RGBA image with 8-bits per component to
* an RGB image with 32-bit fixed-point components, where the number of
* fractional bits is INPUT_FRACBITS.  At the same time, the function pads the
* image so that the result has size
*      (InputWidth + 2*Padding) by (InputHeight + 2*Padding).
* The padding is constant extension (pixel replication).
*/
static void ConvertInput(int32_t *FixedRgb, const uint32_t *Input, int InputWidth,
    int InputHeight, int Padding)
{
    const uint8_t *InputPtr = (uint8_t *)Input;
    const int InputStride = 4*InputWidth;
    const int Stride = PIXEL_STRIDE*(InputWidth + 2*Padding);
    int32_t r, g, b;
    int i, Row;
    
    
    FixedRgb += Padding*Stride;
    
    for(Row = InputHeight; Row; Row--, InputPtr += InputStride)
    {
        r = ((int32_t)InputPtr[0]) << INPUT_FRACBITS;
        g = ((int32_t)InputPtr[1]) << INPUT_FRACBITS;
        b = ((int32_t)InputPtr[2]) << INPUT_FRACBITS;
        
        /* Pad left side by copying pixel */
        for(i = Padding; i; i--)
        {
            *(FixedRgb++) = r;
            *(FixedRgb++) = g;
            *(FixedRgb++) = b;
        }
        
        /* Convert the interior of the image */
        for(i = 0; i < InputStride; i += 4)
        {
            *(FixedRgb++) = ((int32_t)InputPtr[i+0]) << INPUT_FRACBITS;
            *(FixedRgb++) = ((int32_t)InputPtr[i+1]) << INPUT_FRACBITS;
            *(FixedRgb++) = ((int32_t)InputPtr[i+2]) << INPUT_FRACBITS;
        }
        
        r = ((int32_t)InputPtr[i-4]) << INPUT_FRACBITS;
        g = ((int32_t)InputPtr[i-3]) << INPUT_FRACBITS;
        b = ((int32_t)InputPtr[i-2]) << INPUT_FRACBITS;
        
        /* Pad right side by copying pixel */
        for(i = Padding; i; i--)
        {
            *(FixedRgb++) = r;
            *(FixedRgb++) = g;
            *(FixedRgb++) = b;
        }
    }
    
    /* Pad bottom by copying rows */
    for(Row = Padding; Row; Row--, FixedRgb += Stride)
        memcpy(FixedRgb, FixedRgb - Stride, sizeof(int32_t)*Stride);
    
    FixedRgb -= Stride*(InputHeight + Padding);
    
    /* Pad top by coping rows */
    for(Row = Padding; Row; Row--, FixedRgb -= Stride)
        memcpy(FixedRgb - Stride, FixedRgb, sizeof(int32_t)*Stride);
}


/**
* @brief Crop and convert RGB fixed-point image to 32-bit RGBA
*
* @param Output pointer to hold the output converted image data
* @param OutputWidth, OutputHeight cropped image dimensions
* @param FixedRgb the input RGB fixed-point image
* @param Width width of the fixed-point image
*
* \c ConvertOutput is used by \c CWInterp to convert the final interpolation
* from the computation format back to RGBA.
*
* \c ConvertOutput converts from an RGB image with 32-bit fixed-point
* components, where the number of fractional bits is OUPUT_FRACBITS, to RGBA
* with 8-bits per component.  The function also crops the image to have
* dimensions OutputWidth by OutputHeight.  (Width is also needed to know the
* stride length in memory between successive rows of the input image.)  The
* upper-left corner of the cropped image can be specified by adjusted
* \c FixedRgb:
@code
    ConvertOutput(Output, OutputWidth, OutputHeight,
        FixedRgb + x0 + y0*Width, Width);
@endcode
*/
static void ConvertOutput(uint32_t *Output, int OutputWidth, int OutputHeight,
    const int32_t *FixedRgb, int Width)
{
    uint8_t *OutputPtr = (uint8_t *)Output;
    const int CroppedStride = PIXEL_STRIDE*OutputWidth, Stride = PIXEL_STRIDE*Width;
    int32_t r, g, b;
    int i, Row;
    
    
    for(Row = OutputHeight; Row; Row--, FixedRgb += Stride)
        for(i = 0; i < CroppedStride; i += PIXEL_STRIDE)
        {
            /* Convert fixed-point values to integer */
            r = (FixedRgb[i+0] + FIXED_HALF(OUTPUT_FRACBITS)) >> OUTPUT_FRACBITS;
            g = (FixedRgb[i+1] + FIXED_HALF(OUTPUT_FRACBITS)) >> OUTPUT_FRACBITS;
            b = (FixedRgb[i+2] + FIXED_HALF(OUTPUT_FRACBITS)) >> OUTPUT_FRACBITS;
            
            /* Clamp range to [0, 255] and store in Output */
            *(OutputPtr++) = CLAMP(r, 0, 255);
            *(OutputPtr++) = CLAMP(g, 0, 255);
            *(OutputPtr++) = CLAMP(b, 0, 255);
            *(OutputPtr++) = 0xFF;
        }
}


/**
* @brief Contour stencil windowed interpolation
*
* @param Output pointer to memory for holding the interpolated image
* @param Input the input image
* @param InputWidth, InputHeight input image dimensions
* @param Psi \f$\psi\f$ samples computed by \c PreCWInterp
* @param Param cwparams struct of interpolation parameters
*/
int CWInterp(uint32_t *Output, const uint32_t *Input,
    int InputWidth, int InputHeight, const int32_t *Psi, cwparams Param)
{
    const int ScaleFactor = (int)ceil(Param.ScaleFactor);
    const int SupportRadius = (NEIGHRADIUS+1)*ScaleFactor - 1;
    const int SupportWidth = 2*SupportRadius + 1;
    const int SupportSize = SupportWidth*SupportWidth;
    const int StencilMul = NUMNEIGH*SupportSize;
    int *Stencil = NULL;
    int32_t *InputFixed = NULL, *OutputFixed = NULL, *Residual = NULL;
    unsigned long StartTime, StopTime;
    int i, PadInput, pw, ph, Success = 0;
    int32_t ResNorm;
    
    
    /* Iterative refinement is unnecessary if PSF is the Dirac delta */
    if(Param.PsfSigma == 0.0)
        Param.RefinementSteps = 0;
    
    PadInput = 4 + (ScaleFactor + 1)/2;
    pw = InputWidth + 2*PadInput;
    ph = InputHeight + 2*PadInput;
    
    if( !(OutputFixed = (int32_t *)Malloc(sizeof(int32_t)*
            PIXEL_STRIDE*pw*ScaleFactor*ph*ScaleFactor))
        || !(InputFixed = (int32_t *)Malloc(sizeof(int32_t)*PIXEL_STRIDE*pw*ph))
        || !(Residual = (int32_t *)Malloc(sizeof(int32_t)*PIXEL_STRIDE*pw*ph))
        || !(Stencil = (int *)Malloc(sizeof(int)*pw*ph)) )
        goto Catch;

    /* Start timing */
    StartTime = Clock();

    /* Convert 32-bit RGBA pixels to integer array */
    ConvertInput(InputFixed, Input, InputWidth, InputHeight, PadInput);
    
    /* Select the best-fitting contour stencils */
    if(!FitStencils(Stencil, InputFixed, pw, ph, StencilMul))
        goto Catch;
    
    memset(OutputFixed, 0, sizeof(int32_t)*
        3*pw*ScaleFactor*ph*ScaleFactor);
    memset(Residual, 0, sizeof(int32_t)*3*pw*ph);
    
    printf("\n  Iteration   Residual norm\n  -------------------------\n");
    
    /* First interpolation pass */
    CWFirstPass(OutputFixed, ScaleFactor, InputFixed, pw, ph, Stencil, Psi);
    
    /* Iterative refinement */
    for(i = 1; i <= Param.RefinementSteps; i++)
    {
        /* Compute the residual */
        if((ResNorm = CWResidual(Residual, OutputFixed, InputFixed,
            pw, ph, Param)) < 0.0)
            goto Catch;
        
        printf("  %8d %15.8f\n", i, ResNorm/(255.0*256.0));
        
        /* Interpolation refinement pass */
        CWRefinementPass(OutputFixed, ScaleFactor, Residual, pw, ph, Stencil, Psi);
    }
    
    /* Convert output integer array to 32-bit RGBA */
    ConvertOutput(Output, InputWidth*ScaleFactor, InputHeight*ScaleFactor,
        OutputFixed + PIXEL_STRIDE*PadInput*ScaleFactor*(1 + pw*ScaleFactor),
        pw*ScaleFactor);
    
    /* The final interpolation is now complete, stop timing. */
    StopTime = Clock();

    /* Compute the residual norm of the final interpolation.  This
    computation is not included in the CPU timing since it is for
    information purposes only. */
    ResNorm = CWResidual(Residual, OutputFixed, InputFixed, pw, ph, Param);
    printf("  %8d %15.8f\n\n", Param.RefinementSteps + 1, ResNorm/(255.0*256.0));
    
    /* Display the CPU time spent performing the interpolation. */
    printf("  CPU time: %.3f s\n\n", 0.001*(StopTime - StartTime));

    Success = 1;
    
Catch:  /* This label is used for error handling.  If something went wrong
        above (which may be out of memory or a computation error), then
        execution jumps to this point to clean up and exit. */
    Free(Stencil);
    Free(Residual);
    Free(InputFixed);
    Free(OutputFixed);
    return Success;
}


/** @brief Round and clamp double X to integer */
#define ROUNDCLAMP(X,A,B) (((X) < (A)) ? (A) : (((X) > (B)) ? (B) : ROUND(X)))
    

/* The following parameters define the number of fractional bits used for
   signed 32-bit fixedpoint arithmetic in CWSynth2Fixed.  These parameters
   should be large enough for reasonable precision but small enough to
   avoid overflow.  Additionally, the implementation constraints on
   choosing these parameters are
 
       WINDOW_FRACBITS >= UK_FRACBITS,
       TRIG_FRACBITS + XY_FRACBITS - PHITN_FRACBITS >= 1,
       COEFF_FRACBITS + PHI_FRACBITS - UK_FRACBITS >= 1.
 */
#define XY_FRACBITS         15
#define TRIG_FRACBITS       10
#define PHITN_FRACBITS      9
#define PHI_FRACBITS        10
#define COEFF_FRACBITS      10
#define WINDOW_FRACBITS     10
#define UK_FRACBITS         10

#define ROUND_FIXED(X,N)    (((X) + FIXED_HALF(N)) >> (N))
#define FLOAT_TO_FIXED(X,N) ((int32_t)ROUND((X) * FIXED_ONE(N)))
    
/** @brief Arbitrary scale factor interpolation */
int CWArbitraryInterp(uint32_t *Output, int OutputWidth, int OutputHeight,
    const int32_t *Input, int InputWidth, int InputHeight,
    const int *Stencil, const double *InverseA, cwparams Param)
{
    /*int (*Extension)(int, int) = ExtensionMethod[Param.Boundary];*/
    int (*Extension)(int, int) = ConstExtension;
    const int InputNumEl = 3*InputWidth*InputHeight;
    const int ExpTableSize = 1024;
    const double ExpArgScale = 37.0236;
    const double PhiTScale = sqrt(ExpArgScale/2)/Param.PhiSigmaTangent;
    const double PhiNScale = sqrt(ExpArgScale/2)/Param.PhiSigmaNormal;
    
    int32_t *Coeff = NULL, *CoeffPtr, *ExpTable = NULL;
    
    float X, Y, XStart, YStart;
    float Temp, cr[NUMNEIGH], cg[NUMNEIGH], cb[NUMNEIGH];
    int32_t v0[3], v[3], WindowWeight, WindowWeightX[4], WindowWeightY[4];
    int32_t Xpf, Ypf, Weight, uk[3], u[3], DenomSum;
    int32_t CosTableTf[NUMSTENCILS], SinTableTf[NUMSTENCILS];
    int32_t CosTableNf[NUMSTENCILS], SinTableNf[NUMSTENCILS];
    int32_t Pixel;
    int i, k, x, y, m, n, mx, my, nx, ny, S, Success = 0;
    int ix, iy, Cur, Offset;

    
    if(!(Coeff = (int32_t *)Malloc(sizeof(int32_t)*(NUMNEIGH + 1)*InputNumEl))
        || !(ExpTable = (int32_t *)Malloc(sizeof(int32_t)*ExpTableSize)))
        goto Catch;
    
    if(Param.CenteredGrid)
    {
        XStart = (float)(1/Param.ScaleFactor - 1)/2;
        YStart = (float)(1/Param.ScaleFactor - 1)/2;
    }
    else
        XStart = YStart = 0;
    
    for(S = 0; S < NUMSTENCILS; S++)
    {
        CosTableTf[S] = FLOAT_TO_FIXED(
            PhiTScale * cos(StencilOrientation[S]), TRIG_FRACBITS);
        SinTableTf[S] = FLOAT_TO_FIXED(
            PhiTScale * sin(StencilOrientation[S]), TRIG_FRACBITS);
        CosTableNf[S] = FLOAT_TO_FIXED(
            PhiNScale * cos(StencilOrientation[S]), TRIG_FRACBITS);
        SinTableNf[S] = FLOAT_TO_FIXED(
            PhiNScale * sin(StencilOrientation[S]), TRIG_FRACBITS);
    }
    
    for(i = 0; i < ExpTableSize; i++)
        ExpTable[i] = FLOAT_TO_FIXED(exp( -(double)(i + 0.5f)/ExpArgScale), PHI_FRACBITS);

    for(y = 0, k = 0, CoeffPtr = Coeff; y < InputHeight; y++)
        for(x = 0; x < InputWidth; x++, k++, CoeffPtr += 3*(NUMNEIGH + 1))
        {
            S = NUMNEIGH*Stencil[k];
            
            v0[0] = Input[3*k + 0];
            v0[1] = Input[3*k + 1];
            v0[2] = Input[3*k + 2];
             
            for(m = 0; m < NUMNEIGH; m++)
                cr[m] = 0;
            for(m = 0; m < NUMNEIGH; m++)
                cg[m] = 0;
            for(m = 0; m < NUMNEIGH; m++)
                cb[m] = 0;
            
            for(ny = -NEIGHRADIUS, n = 0; ny <= NEIGHRADIUS; ny++)
            {
                Offset = InputWidth*Extension(InputHeight, y + ny);
                
                for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n++)
                {
                    if(n != 1 + (2*NEIGHRADIUS + 1))
                    {
                        i = 3*(Extension(InputWidth, x + nx) + Offset);
                        v[0] = Input[i + 0];
                        v[1] = Input[i + 1];
                        v[2] = Input[i + 2];
                        
                        for(m = 0; m < NUMNEIGH; m++)
                        {
                            Temp = (float)InverseA[m + NUMNEIGH*(n + S)];
                            cr[m] += Temp * (v[0] - v0[0]);
                            cg[m] += Temp * (v[1] - v0[1]);
                            cb[m] += Temp * (v[2] - v0[2]);
                        }
                    }
                }
            }
            
            /* The first three coeff values have UK_FRACBITS fractional bits */
            CoeffPtr[0] = v0[0] << (UK_FRACBITS - INPUT_FRACBITS);
            CoeffPtr[1] = v0[1] << (UK_FRACBITS - INPUT_FRACBITS);
            CoeffPtr[2] = v0[2] << (UK_FRACBITS - INPUT_FRACBITS);
            
            /* The other values have COEFF_FRACBITS fractional bits */
            for(m = 0; m < NUMNEIGH; m++)
            {
                CoeffPtr[3*(m + 1) + 0] = FLOAT_TO_FIXED(cr[m]
                    / FIXED_ONE(INPUT_FRACBITS), COEFF_FRACBITS);
                CoeffPtr[3*(m + 1) + 1] = FLOAT_TO_FIXED(cg[m]
                    / FIXED_ONE(INPUT_FRACBITS), COEFF_FRACBITS);
                CoeffPtr[3*(m + 1) + 2] = FLOAT_TO_FIXED(cb[m]
                    / FIXED_ONE(INPUT_FRACBITS), COEFF_FRACBITS);
            }
            
        }
    
    for(y = 0, k = 0; y < OutputHeight; y++)
    {
        Y = YStart + (float)(y/Param.ScaleFactor);
        iy = (int)ceil(Y - 2*NEIGHRADIUS);

        /* Precompute y-factor of the window weights */
        for(my = 0; my < 4*NEIGHRADIUS; my++)
            WindowWeightY[my] = FLOAT_TO_FIXED(
                CubicBSpline(Y - (iy + my)), WINDOW_FRACBITS);

        for(x = 0; x < OutputWidth; x++, k++)
        {
            X = XStart + (float)(x/Param.ScaleFactor);
            ix = (int)ceil(X - 2*NEIGHRADIUS);

            /* Precompute x-factor of the window weights */
            for(mx = 0; mx < 4*NEIGHRADIUS; mx++)
                WindowWeightX[mx] = FLOAT_TO_FIXED(
                    CubicBSpline(X - (ix + mx)), WINDOW_FRACBITS);
            
            DenomSum = 0;
            u[0] = u[1] = u[2] = 0;
            
            for(my = 0, Ypf = (int32_t)ROUND((Y - iy) * FIXED_ONE(XY_FRACBITS)); my < 4*NEIGHRADIUS;
                my++, Ypf -= FIXED_ONE(XY_FRACBITS))
            if((iy + my) >= 0 && (iy + my) < InputHeight)
            {
                for(mx = 0, Xpf = (int32_t)ROUND((X - ix) * FIXED_ONE(XY_FRACBITS)),
                    i = ix + InputWidth*(iy + my), CoeffPtr = Coeff + i*3*(NUMNEIGH + 1);
                    mx < 4*NEIGHRADIUS; mx++, i++, CoeffPtr += 3*(NUMNEIGH + 1), Xpf -= FIXED_ONE(XY_FRACBITS))
                {
                    if((ix + mx) < 0 || (ix + mx) >= InputWidth)
                        continue;
                    
                    /* WindowWeight has 2*WINDOW_FRACBITS fractional bits. */
                    WindowWeight = (WindowWeightX[mx] * WindowWeightY[my]);
                    /* DenomSum is computed using 2*WINDOW_FRACBITS. */
                    DenomSum += WindowWeight;
                    /* Now reduce to WindowWeight to WINDOW_FRACBITS. */
                    WindowWeight = (WindowWeight + FIXED_HALF(WINDOW_FRACBITS))
                        >> WINDOW_FRACBITS;
                    
                    if(!WindowWeight)
                        continue;
                                        
                    S = Stencil[i];
                                        
                    uk[0] = CoeffPtr[0];
                    uk[1] = CoeffPtr[1];
                    uk[2] = CoeffPtr[2];
                    
                    for(ny = -NEIGHRADIUS, n = 3; ny <= NEIGHRADIUS; ny++)
                        for(nx = -NEIGHRADIUS; nx <= NEIGHRADIUS; nx++, n += 3)
                        {
                            int32_t phit, phin;
                            
                            /* The tables use TRIG_FRACBITS fractional bits,
                               and X and Y use XY_FRACBITS, so the products
                               have (TRIG_FRACBITS + XY_FRACBITS) fractional
                               bits.  The shift reduces the result to
                               PHITN_FRACBITS. */
                            phit = ( CosTableTf[S]*(Xpf - nx*FIXED_ONE(XY_FRACBITS))
                                + SinTableTf[S]*(Ypf - ny*FIXED_ONE(XY_FRACBITS))
                                + FIXED_HALF(TRIG_FRACBITS + XY_FRACBITS - PHITN_FRACBITS))
                                >> (TRIG_FRACBITS + XY_FRACBITS - PHITN_FRACBITS);
                            phin = ( -SinTableNf[S]*(Xpf - nx*FIXED_ONE(XY_FRACBITS))
                                + CosTableNf[S]*(Ypf - ny*FIXED_ONE(XY_FRACBITS))
                                + FIXED_HALF(TRIG_FRACBITS + XY_FRACBITS - PHITN_FRACBITS))
                                >> (TRIG_FRACBITS + XY_FRACBITS - PHITN_FRACBITS);
                             
                            /* phit and phin have PHITN_FRACBITS, so the
                               products have 2*PHITN_FRACBITS.  The result is
                               shifted by 2*PHITN_FRACBITS to convert to
                               quantity (by floor rounding) to an integer. */
                            Cur = (phit*phit + phin*phin) >> (2*PHITN_FRACBITS);
                                
                            if(Cur >= ExpTableSize)
                                continue;
                            
                            /* Compute exp(-Cur) via table look up.  The result
                               has PHI_FRACBITS fractional bits. */
                            Weight = ExpTable[Cur];
  
                            /* The Coeff values have COEFF_FRACBITS fractional
                               bits and Weight has PHI_FRACBITS.  The products
                               are shifted so that the result has
                               WINDOW_FRACBITS. */
                            uk[0] += (CoeffPtr[n + 0] * Weight
                                + FIXED_HALF(COEFF_FRACBITS + PHI_FRACBITS - WINDOW_FRACBITS))
                                >> (COEFF_FRACBITS + PHI_FRACBITS - UK_FRACBITS);
                            uk[1] += (CoeffPtr[n + 1] * Weight
                                + FIXED_HALF(COEFF_FRACBITS + PHI_FRACBITS - WINDOW_FRACBITS))
                                >> (COEFF_FRACBITS + PHI_FRACBITS - UK_FRACBITS);
                            uk[2] += (CoeffPtr[n + 2] * Weight
                                + FIXED_HALF(COEFF_FRACBITS + PHI_FRACBITS - WINDOW_FRACBITS))
                                >> (COEFF_FRACBITS + PHI_FRACBITS - UK_FRACBITS);
                        }
                    
                    /* u is computed using WINDOW_FRACBITS + UK_FRACBITS. */
                    u[0] += WindowWeight * uk[0];
                    u[1] += WindowWeight * uk[1];
                    u[2] += WindowWeight * uk[2];
                }
            }

            if(DenomSum >= FIXED_ONE(2*WINDOW_FRACBITS) - 1)
            {
                u[0] = ROUND_FIXED(u[0], WINDOW_FRACBITS + UK_FRACBITS);
                u[1] = ROUND_FIXED(u[1], WINDOW_FRACBITS + UK_FRACBITS);
                u[2] = ROUND_FIXED(u[2], WINDOW_FRACBITS + UK_FRACBITS);
            }
            else
            {
                /* Reduce DenomSum from 2*WINDOW_FRACBITS fractional bits to
                (WINDOW_FRACBITS + UK_FRACBITS) fractional bits. */
#if WINDOW_FRACBITS > UK_FRACBITS
                DenomSum = (DenomSum + FIXED_HALF(WINDOW_FRACBITS - UK_FRACBITS))
                    >> (WINDOW_FRACBITS - UK_FRACBITS);
#endif
                /* u and DenomSum both have (WINDOW_FRACBITS + UK_FRACBITS)
                fractional bits, so the quotient is integer. */
                u[0] = (u[0] + DenomSum/2) / DenomSum;
                u[1] = (u[1] + DenomSum/2) / DenomSum;
                u[2] = (u[2] + DenomSum/2) / DenomSum;
            }
            
            Pixel = 0xFFFFFFFF;
            ((uint8_t *)&Pixel)[0] = CLAMP(u[0],0,255);
            ((uint8_t *)&Pixel)[1] = CLAMP(u[1],0,255);
            ((uint8_t *)&Pixel)[2] = CLAMP(u[2],0,255);
            Output[k] = Pixel;
        }
    }
              
    Success = 1;
Catch:
    Free(ExpTable);
    Free(Coeff);
    return Success;
}


/** @brief Adds residual back to input for refinement passes */
static void AddResidual(int32_t *InputAdjusted, const int32_t *Residual,
    int InputWidth, int InputHeight, int PadInput)
{
    const int PadWidth = InputWidth + 2*PadInput;
    const int RowEl = 3*InputWidth;
    const int PadRowEl = 3*PadWidth;
    int i, Row;
    
    Residual += 3*(PadInput + PadInput*PadWidth);
        
    for(Row = InputHeight; Row; Row--)
    {
        for(i = 0; i < RowEl; i++)
            InputAdjusted[i] += Residual[i];
        
        InputAdjusted += RowEl;
        Residual += PadRowEl;
    }
}


/** @brief Removes padding from Stencil array */
static void StencilStripPad(int *Stencil,
    int InputWidth, int InputHeight, int PadInput, int StencilMul)
{
    const int PadWidth = InputWidth + 2*PadInput;
    int *Src;
    int i, Row;
    
    Src = Stencil + (PadInput + PadInput*PadWidth);
        
    for(Row = InputHeight; Row; Row--)
    {
        for(i = 0; i < InputWidth; i++)
            Stencil[i] = Src[i] / StencilMul;
            
        Stencil += InputWidth;
        Src += PadWidth;
    }
}

/**
 * @brief Contour stencil windowed interpolation for arbitrary scale factors
 *
 * @param Output pointer to memory for holding the interpolated image
 * @param OutputWidth, OutputHeight output image dimensions
 * @param Input the input image
 * @param InputWidth, InputHeight input image dimensions
 * @param Psi \f$\psi\f$ samples computed by \c PreCWInterp
 * @param Param cwparams struct of interpolation parameters
 */
int CWInterpEx(uint32_t *Output, int OutputWidth, int OutputHeight,
    const uint32_t *Input, int InputWidth, int InputHeight,
    const int32_t *Psi, cwparams Param)
{
    const int ScaleFactor = (int)ceil(Param.ScaleFactor);
    const int SupportRadius = (NEIGHRADIUS+1)*ScaleFactor - 1;
    const int SupportWidth = 2*SupportRadius + 1;
    const int SupportSize = SupportWidth*SupportWidth;
    const int StencilMul = NUMNEIGH*SupportSize;
    double *InverseA = NULL;
    int *Stencil = NULL;
    int32_t *InputFixed = NULL, *InputAdjusted = NULL, *OutputFixed = NULL, *Residual = NULL;
    unsigned long StartTime, StopTime;
    int i, PadInput, pw, ph, Success = 0;
    int32_t ResNorm;
    
    
    /* Iterative refinement is unnecessary if PSF is the Dirac delta */
    if(Param.PsfSigma == 0.0)
        Param.RefinementSteps = 0;
    
    PadInput = 4 + (ScaleFactor + 1)/2;
    pw = InputWidth + 2*PadInput;
    ph = InputHeight + 2*PadInput;
    
    if( !(InverseA = (double *)Malloc(sizeof(double)*NUMNEIGH*NUMNEIGH*NUMSTENCILS))
        || !(OutputFixed = (int32_t *)Malloc(sizeof(int32_t)*
            PIXEL_STRIDE*pw*ScaleFactor*ph*ScaleFactor))
        || !(InputFixed = (int32_t *)Malloc(sizeof(int32_t)*PIXEL_STRIDE*pw*ph))
        || !(InputAdjusted = (int32_t *)Malloc(sizeof(int32_t)*PIXEL_STRIDE*InputWidth*InputHeight))
        || !(Residual = (int32_t *)Malloc(sizeof(int32_t)*PIXEL_STRIDE*pw*ph))
        || !(Stencil = (int *)Malloc(sizeof(int)*pw*ph)) )
        goto Catch;
    
    if(!ComputeMatrices(InverseA, Param))
        goto Catch;
    
    /* Start timing */
    StartTime = Clock();

    if(Param.RefinementSteps > 0)
    {
        /* Convert 32-bit RGBA pixels to integer array */
        ConvertInput(InputFixed, Input, InputWidth, InputHeight, PadInput);
        
        memset(InputAdjusted, 0, sizeof(int32_t)*3*InputWidth*InputHeight);
        AddResidual(InputAdjusted, InputFixed, InputWidth, InputHeight, PadInput);
        
        /* Select the best-fitting contour stencils */
        if(!FitStencils(Stencil, InputFixed, pw, ph, StencilMul))
            goto Catch;
        
        memset(OutputFixed, 0, sizeof(int32_t)*
            3*pw*ScaleFactor*ph*ScaleFactor);
        memset(Residual, 0, sizeof(int32_t)*3*pw*ph);
        
        printf("\n  Iteration   Residual norm\n  -------------------------\n");
        
        /* First interpolation pass */
        CWFirstPass(OutputFixed, ScaleFactor, InputFixed, pw, ph, Stencil, Psi);
        
        /* Iterative refinement */
        for(i = 1; i <= Param.RefinementSteps; i++)
        {
            /* Compute the residual */
            if((ResNorm = CWResidual(Residual, OutputFixed, InputFixed,
                pw, ph, Param)) < 0.0)
                goto Catch;
            
            printf("  %8d %15.8f\n", i, ResNorm/(255.0*256.0));
        
            AddResidual(InputAdjusted, Residual, InputWidth, InputHeight, PadInput);
            
            if(i < Param.RefinementSteps)
            {
                /* Interpolation refinement pass */
                CWRefinementPass(OutputFixed, ScaleFactor, Residual, pw, ph, Stencil, Psi);
            }
        }
        
        StencilStripPad(Stencil, InputWidth, InputHeight, PadInput, StencilMul);
    }
    else
    {
        /* Convert 32-bit RGBA pixels to integer array */
        ConvertInput(InputAdjusted, Input, InputWidth, InputHeight, 0);
        
        /* Select the best-fitting contour stencils */
        if(!FitStencils(Stencil, InputAdjusted, InputWidth, InputHeight, 1))
            goto Catch;
    }
    
    if(!CWArbitraryInterp(Output, OutputWidth, OutputHeight,
        InputAdjusted, InputWidth, InputHeight, Stencil, InverseA, Param))
        goto Catch;
        
    /* The final interpolation is now complete, stop timing. */
    StopTime = Clock();

    if(Param.RefinementSteps > 1)
        printf("  %8d   (not computed)\n\n", Param.RefinementSteps + 1);
    
    /* Display the CPU time spent performing the interpolation. */
    printf("  CPU time: %.3f s\n\n", 0.001*(StopTime - StartTime));

    Success = 1;
    
Catch:  /* This label is used for error handling.  If something went wrong
        above (which may be out of memory or a computation error), then
        execution jumps to this point to clean up and exit. */
    Free(Stencil);
    Free(Residual);
    Free(InputAdjusted);
    Free(InputFixed);
    Free(OutputFixed);
    Free(InverseA);
    return Success;
}


/** @brief Display the estimated contour orientations */
int DisplayContours(uint32_t *Output, int OutputWidth, int OutputHeight,
    uint32_t *Input, int InputWidth, int InputHeight, cwparams Param)
{
    const int Pad = 2;
    const float LineColor[3] = {0, 0, 0};
    int *Stencil = 0;
    float dx, dy;
    int32_t *InputInt = NULL;
    uint32_t Pixel;
    int x, y, S, pw, ph, Success = 0;
    
    
    pw = InputWidth + 2*Pad;
    ph = InputHeight + 2*Pad;
    
    if( !(InputInt = (int32_t *)Malloc(sizeof(int32_t)*3*pw*ph))
        || !(Stencil = (int *)Malloc(sizeof(int)*pw*ph)) )
        goto Catch;
    
    /* Convert 32-bit RGBA pixels to integer array */
    ConvertInput(InputInt, Input, InputWidth, InputHeight, Pad);
    
    /* Select the best-fitting contour stencils */
    if(!FitStencils(Stencil, InputInt, pw, ph, 1))
        goto Catch;
    
    /* Lighten the image */
    for(y = 0; y < InputHeight; y++)
        for(x = 0; x < InputWidth; x++)
        {
            Pixel = Input[x + InputWidth*y];
            ((uint8_t*)&Pixel)[0] = (uint8_t)(((uint8_t*)&Pixel)[0]/2 + 128);
            ((uint8_t*)&Pixel)[1] = (uint8_t)(((uint8_t*)&Pixel)[1]/2 + 128);
            ((uint8_t*)&Pixel)[2] = (uint8_t)(((uint8_t*)&Pixel)[2]/2 + 128);
            Input[x + InputWidth*y] = Pixel;
        }
        
    /* Nearest neighbor interpolation */
    NearestInterp(Output, OutputWidth, OutputHeight,
                  Input, InputWidth, InputHeight,
                  (float)Param.ScaleFactor, Param.CenteredGrid);
    
    /* Draw contour orientation lines */
    for(y = 0; y < InputHeight; y++)
        for(x = 0; x < InputWidth; x++)
        {
            S = Stencil[(x + Pad) + pw*(y + Pad)];
            dx = (float)cos(StencilOrientation[S])*0.6f;
            dy = (float)sin(StencilOrientation[S])*0.6f;
            DrawLine(Output, OutputWidth, OutputHeight,
                (float)Param.ScaleFactor*(x - dx + 0.5f) - 0.5f,
                (float)Param.ScaleFactor*(y - dy + 0.5f) - 0.5f,
                (float)Param.ScaleFactor*(x + dx + 0.5f) - 0.5f,
                (float)Param.ScaleFactor*(y + dy + 0.5f) - 0.5f, LineColor);
        }

    Success = 1;
Catch:  /* This label is used for error handling.  If something went wrong
        above (which may be out of memory or a computation error), then
        execution jumps to this point to clean up and exit. */
    Free(Stencil);
    Free(InputInt);
    return Success;
}

    