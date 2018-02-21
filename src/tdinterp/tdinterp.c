/**
 * @file tdinterp.c
 * @brief Roussos-Maragos interpolation by tensor-driven diffusion
 * @author Pascal Getreuer <getreuer@gmail.com>
 *
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
#include <fftw3.h>

#include "basic.h"
#include "finterp.h"
#include "tdinterp.h"
#include "conv.h"


/** @brief Compute the X-derivative for Weicker-Scharr scheme */
static void XDerivative(float *Dest, float *ConvTemp, const float *Src,
    int Width, int Height)
{
    int i, il, ir, y;
    int iEnd;
    
    
    for(y = 0, i = 0, iEnd = -1; y < Height; y++)
    {
        ConvTemp[i] = Src[i + 1] - Src[i];
        i++;
        
        for(iEnd += Width; i < iEnd; i++)
            ConvTemp[i] = Src[i + 1] - Src[i - 1];
        
        ConvTemp[i] = Src[i] - Src[i - 1];
        i++;
    }
    
    for(y = 0, i = iEnd = 0; y < Height; y++)
    {
        il = (y > 0) ? -Width : 0;
        ir = (y < Height - 1) ? Width : 0;
        
        for(iEnd += Width; i < iEnd; i++)
            Dest[i] = 0.3125f*ConvTemp[i]
                + 0.09375f*(ConvTemp[i + ir] + ConvTemp[i + il]);
    }
}


/** @brief Compute the Y-derivative for Weicker-Scharr scheme */
static void YDerivative(float *Dest, float *ConvTemp, const float *Src,
    int Width, int Height)
{
    int i, il, ir, iEnd, y;
    
    
    for(y = 0, i = iEnd = 0; y < Height; y++)
    {
        il = (y > 0) ? -Width : 0;
        ir = (y < Height - 1) ? Width : 0;
        
        for(iEnd += Width; i < iEnd; i++)
            ConvTemp[i] = Src[i + ir] - Src[i + il];
    }
    
    for(y = 0, i = 0, iEnd = -1; y < Height; y++)
    {
        Dest[i] = 0.40625f*ConvTemp[i] + 0.09375f*ConvTemp[i + 1];
        i++;
        
        for(iEnd += Width; i < iEnd; i++)
            Dest[i] =  0.3125f*ConvTemp[i]
                + 0.09375f*(ConvTemp[i + 1] + ConvTemp[i - 1]);
        
        Dest[i] =  0.40625f*ConvTemp[i] + 0.09375f*ConvTemp[i - 1];
        i++;
    }
}


/** @brief Construct the tensor T */
static void ComputeTensor(float *Txx, float *Txy, float *Tyy,
    float *uSmooth, float *ConvTemp, const float *u, int Width, int Height,
    filter PreSmooth, filter PostSmooth, double K)
{
    const float dt = 2;
    const double KSquared = K*K;
    const int NumPixels = Width*Height;
    const int NumEl = 3*NumPixels;
    boundaryext Boundary = GetBoundaryExt("sym");
    double Tm, SqrtTm, Trace, Lambda1, Lambda2, EigVecX, EigVecY, Temp;
    float ux, uy;
    int i, iu, id, il, ir, x, y, Channel;
    

    /* Set the tensor to zero and perform pre-smoothing on u.  Note that it is
    not safely portable to use memset for this purpose.
    http://c-faq.com/malloc/calloc.html  */
#ifdef _OPENMP
    #pragma omp parallel sections private(i)
    {	/* The following six operations can run in parallel */
        #pragma omp section
        {
            for(i = 0; i < NumPixels; i++)
                Txx[i] = 0.0f;
        }
        #pragma omp section
        {
            for(i = 0; i < NumPixels; i++)
                Txy[i] = 0.0f;
        }
        #pragma omp section
        {
            for(i = 0; i < NumPixels; i++)
                Tyy[i] = 0.0f;
        }
        #pragma omp section
        {
            SeparableConv2D(uSmooth,
                ConvTemp, u,
                PreSmooth, PreSmooth, Boundary, Width, Height, 1);
        }
        #pragma omp section
        {
            SeparableConv2D(uSmooth + NumPixels,
                ConvTemp + NumPixels, u + NumPixels,
                PreSmooth, PreSmooth, Boundary, Width, Height, 1);
        }
        #pragma omp section
        {
            SeparableConv2D(uSmooth + 2*NumPixels,
                ConvTemp + 2*NumPixels, u + 2*NumPixels,
                PreSmooth, PreSmooth, Boundary, Width, Height, 1);
        }
    }
#else
    for(i = 0; i < NumPixels; i++)
        Txx[i] = 0.0f;
    
    for(i = 0; i < NumPixels; i++)
        Txy[i] = 0.0f;
    
    for(i = 0; i < NumPixels; i++)
        Tyy[i] = 0.0f;

    SeparableConv2D(uSmooth, ConvTemp, u,
        PreSmooth, PreSmooth, Boundary, Width, Height, 3);
#endif
    
    /* Compute the structure tensor */
    for(Channel = 0; Channel < NumEl; Channel += NumPixels)
    {
        for(y = 0, i = 0; y < Height; y++)
        {
            iu = (y > 0) ? -Width : 0;
            id = (y < Height - 1) ? Width : 0;
            
            for(x = 0; x < Width; x++, i++)
            {
                il = (x > 0) ? -1 : 0;
                ir = (x < Width - 1) ? 1 : 0;
                
                ux = (uSmooth[i + ir + Channel]
                    - uSmooth[i + il + Channel]) / 2;
                uy = (uSmooth[i + id + Channel]
                    - uSmooth[i + iu + Channel]) / 2;
                Txx[i] += ux * ux;
                Txy[i] += ux * uy;
                Tyy[i] += uy * uy;
            }
        }
    }
    
    /* Perform the post-smoothing convolution with PostSmooth */
#ifdef _OPENMP
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            SeparableConv2D(Txx, ConvTemp, Txx,
                PostSmooth, PostSmooth, Boundary, Width, Height, 1);
        }
        #pragma omp section
        {
            SeparableConv2D(Txy, ConvTemp + NumPixels, Txy,
                PostSmooth, PostSmooth, Boundary, Width, Height, 1);
        }
        #pragma omp section
        {
            SeparableConv2D(Tyy, ConvTemp + 2*NumPixels, Tyy,
                PostSmooth, PostSmooth, Boundary, Width, Height, 1);
        }
    }
#else
    SeparableConv2D(Txx, ConvTemp, Txx,
        PostSmooth, PostSmooth, Boundary, Width, Height, 1);
    SeparableConv2D(Txy, ConvTemp, Txy,
        PostSmooth, PostSmooth, Boundary, Width, Height, 1);
    SeparableConv2D(Tyy, ConvTemp, Tyy,
        PostSmooth, PostSmooth, Boundary, Width, Height, 1);
#endif
    
    /* Refactor the structure tensor */
    for(i = 0; i < NumPixels; i++)
    {
        /* Compute the eigenspectra */
        Trace = 0.5*(Txx[i] + Tyy[i]);
        Temp = sqrt(Trace*Trace - Txx[i]*Tyy[i] + Txy[i]*Txy[i]);
        Lambda1 = Trace - Temp;
        Lambda2 = Trace + Temp;
        EigVecX = Txy[i];
        EigVecY = Lambda1 - Txx[i];
        Temp = sqrt(EigVecX*EigVecX + EigVecY*EigVecY);
        
        if(Temp >= 1e-9)
        {
            EigVecX /= Temp;
            EigVecY /= Temp;
            Tm = KSquared/(KSquared + (Lambda1 + Lambda2));
            SqrtTm = sqrt(Tm);
            
            /* Construct new tensor from the spectra */
            Txx[i] = (float)(dt*(SqrtTm*EigVecX*EigVecX + Tm*EigVecY*EigVecY));
            Txy[i] = (float)(dt*((SqrtTm - Tm)*EigVecX*EigVecY));
            Tyy[i] = (float)(dt*(SqrtTm*EigVecY*EigVecY + Tm*EigVecX*EigVecX));
        }
        else
        {
            Txx[i] = dt;
            Txy[i] = 0.0f;
            Tyy[i] = dt;
        }
    }
}


/** @brief Perform 5x5 Weicker-Scharr explicit diffusion steps */
static void DiffuseWithTensor(float *u, float *vx, float *vy,
    float *SumX, float *SumY,
    float *ConvTemp, const float *Txx, const float *Txy, const float *Tyy,
    int Width, int Height, int DiffIter)
{
    const int NumPixels = Width*Height;
    int i, Channel, Step;
    
    
    for(Channel = 0; Channel < 3; Channel++, u += NumPixels)
    {
        for(Step = 0; Step < DiffIter; Step++)
        {
#ifdef _OPENMP
            #pragma omp parallel
            {
                #pragma omp sections
                {	/* XDerivative and YDerivative can run in parallel */
                    #pragma omp section
                    {
                        XDerivative(vx, ConvTemp, u, Width, Height);
                    }
                    #pragma omp section
                    {
                        YDerivative(vy, ConvTemp + NumPixels, u, Width, Height);
                    }
                }
                
                #pragma omp for schedule(static)
                for(i = 0; i < NumPixels; i++)
                {
                    SumX[i] = Txx[i]*vx[i] + Txy[i]*vy[i];
                    SumY[i] = Txy[i]*vx[i] + Tyy[i]*vy[i];
                }

                #pragma omp sections
                {	/* XDerivative and YDerivative can run in parallel */
                    #pragma omp section
                    {
                        XDerivative(SumX, ConvTemp,
                                    SumX, Width, Height);
                    }
                    #pragma omp section
                    {
                        YDerivative(SumY, ConvTemp + NumPixels,
                                    SumY, Width, Height);
                    }
                }
                
                #pragma omp for schedule(static)
                for(i = 0; i < NumPixels; i++)
                    u[i] += SumX[i] + SumY[i];
            }
#else
            XDerivative(vx, ConvTemp, u, Width, Height);
            YDerivative(vy, ConvTemp, u, Width, Height);
            
            for(i = 0; i < NumPixels; i++)
            {
                SumX[i] = Txx[i]*vx[i] + Txy[i]*vy[i];
                SumY[i] = Txy[i]*vx[i] + Tyy[i]*vy[i];
            }
            
            XDerivative(SumX, ConvTemp, SumX, Width, Height);
            YDerivative(SumY, ConvTemp, SumY, Width, Height);
            
            for(i = 0; i < NumPixels; i++)
                u[i] += SumX[i] + SumY[i];
#endif
        }
    }
}


/** @brief Sum aliases (used by MakePhi) */
static void SumAliases(float *u, int d, int Width, int Height)
{
    const int BlockWidth = Width / d;
    const int BlockHeight = Height / d;
    const int StripSize = Width*BlockHeight;
    float *Src, *Dest;
    int k, x, y;
    
    
    /* Sum along x for every y */
    for(y = 0; y < Height; y++)
    {
        Dest = u + y*Width;
        
        for(k = 1; k < d; k++)
        {
            Src = Dest + k*BlockWidth;
            
            for(x = 0; x < BlockWidth; x++)
                Dest[x] += Src[x];
        }
    }
    
    /* Sum along y for 0 <= x < BlockWidth */
    for(k = 1; k < d; k++)
    {
        Dest = u;
        Src = u + Width*(k*BlockHeight);
        
        for(y = 0; y < BlockHeight; y++)
        {
            for(x = 0; x < BlockWidth; x++)
                Dest[x] += Src[x];
            
            Dest += Width;
            Src += Width;
        }
    }
    
    /* Copy block [0, BlockWidth) x [0, BlockHeight) in the x direction*/
    for(y = 0, Src = Dest = u; y < Height; y++, Src += Width)
    {
        for(k = 1, Dest += BlockWidth; k < d; k++, Dest += BlockWidth)
        {
            memcpy(Dest, Src, sizeof(float)*BlockWidth);
        }
    }
    
    /* Copy strip [0, OutputWidth) x [0, InputHeight) in the y direction*/
    for(k = 1, Dest = u + StripSize; k < d; k++, Dest += StripSize)
        memcpy(Dest, u, sizeof(float)*StripSize);
}


/** @brief Precompute the function phi used in projections */
static void MakePhi(float *Phi, double PsfSigma, int d, int Width, int Height)
{
    /* Number of terms to use in truncated sum, larger is more accurate */
    const int NumOverlaps = 2;
    const int NumPixels = Width*Height;
    float *PhiLastRow;
    float Sum, Sigma, Denom;
    int i, x, y, t, k;
    
    float *Temp;

    
    Temp = (float *)Malloc(sizeof(float)*NumPixels);

    /* Construct the Fourier transform of a Gaussian with spatial standard
     * deviation d*PsfSigma.  Using the Gaussian and the transform's
     * separability, the result is formed as the tensor product of the 1D
     * transforms.  This significantly reduces the number of exp calls.
     */
    if(PsfSigma == 0.0)
    {
        for(i = 0; i < NumPixels; i++)
            Phi[i] = 1.0f;
    }
    else
    {
        Sigma = (float)(Width/(2*M_PI*d*PsfSigma));
        Denom = 2*Sigma*Sigma;
        PhiLastRow = Phi + Width*(Height - 1);

        /* Construct the transform along x */
        for(x = 0; x < Width; x++)
        {
            for(k = -NumOverlaps, Sum = 0; k < NumOverlaps; k++)
            {
                t = x + k*Width;
                Sum += (float)exp(-(t*t)/Denom);
            }

            PhiLastRow[x] = Sum;
        }

        Sigma = (float)(Height/(2*M_PI*d*PsfSigma));
        Denom = 2*Sigma*Sigma;

        /* Construct the transform along y */
        for(y = 0, i = 0; y < Height; y++, i += Width)
        {
            for(k = -NumOverlaps, Sum = 0; k < NumOverlaps; k++)
            {
                t = y + k*Height;
                Sum += (float)exp(-(t*t)/Denom);
            }

            /* Tensor product */
            for(x = 0; x < Width; x++)
                Phi[x + i] = Sum*PhiLastRow[x];
        }
    }
    
    /* Square all the elements so that Temp = abs(fft2(psf)).^2 */
    for(i = 0; i < NumPixels; i++)
        Temp[i] = Phi[i]*Phi[i];
    
    SumAliases(Temp, d, Width, Height);
    
    /* Obtain the projection kernel Phi in the Fourier domain */
    for(i = 0; i < NumPixels; i++)
        Phi[i] /= (float)sqrt(Temp[i] * NumPixels);
    
    Free(Temp);
}


/** @brief Sum the Fourier aliases */
static void ImageSumAliases(float *u, int d, int Width, int Height,
    int NumChannels)
{
    const int BlockWidth = Width / d;
    const int BlockWidth2 = 2*BlockWidth;
    const int BlockHeight = Height / d;
    const int StripSize = 2*Width*BlockHeight;
    float *uChannel, *Src, *Dest;
    int k, i, y, Channel;
    
    
    for(Channel = 0; Channel < NumChannels; Channel++)
    {
        uChannel = u + 2*Width*Height*Channel;
        
        /* Sum along x for every y */
        for(y = 0; y < Height; y++)
        {
            Dest = uChannel + 2*Width*y;

            for(k = 1; k < d; k++)
            {
                Src = Dest + 2*BlockWidth*k;

                for(i = 0; i < BlockWidth2; i++)
                    Dest[i] += Src[i];
            }
        }

        /* Sum along y for 0 <= x < BlockWidth */
        for(k = 1; k < d; k++)
        {
            Dest = uChannel;
            Src = uChannel + 2*Width*(BlockHeight*k);

            for(y = 0; y < BlockHeight; y++)
            {
                for(i = 0; i < BlockWidth2; i++)
                    Dest[i] += Src[i];

                Dest += 2*Width;
                Src += 2*Width;
            }
        }

        /* Copy block [0, BlockWidth) x [0, BlockHeight) in the x direction */
        for(y = 0, Src = Dest = uChannel; y < Height; y++, Src += 2*Width)
        {
            for(k = 1, Dest += BlockWidth2; k < d; k++, Dest += BlockWidth2)
            {
                memcpy(Dest, Src, sizeof(float)*BlockWidth2);
            }
        }

        /* Copy strip [0, OutputWidth) x [0, InputHeight) in the y direction */
        for(k = 1, Dest = uChannel + StripSize; k < d; k++, Dest += StripSize)
            memcpy(Dest, uChannel, sizeof(float)*StripSize);
    }
}


/** @brief Complex multiply of Dest by Phi */
static void MultiplyByPhiInterleaved(float *Dest, const float *Phi,
    int NumPixels, int NumChannels)
{
    int i, Channel;

    
    for(Channel = 0; Channel < NumChannels; Channel++, Dest += 2*NumPixels)
    {
        for(i = 0; i < NumPixels; i++)
        {
            Dest[2*i] *= Phi[i];
            Dest[2*i + 1] *= Phi[i];
        }
    }
}


/** @brief Fill the redundant spectra in a Fourier transform */
static void MirrorSpectra(float *Dest, int Width, int Height, int NumChannels)
{
    const int H = Width/2 + 1;
    int i, sx, sy, x, y, Channel;
    
    
    /* Fill in the redundant spectra */
    for(Channel = 0, i = 0; Channel < NumChannels; Channel++)
    {
        for(y = 0; y < Height; y++)
        {
            sy = (y > 0) ? Height - y : 0;
            sy = 2*Width*(sy + Height*Channel);
            
            for(x = H, i += 2*H; x < Width; x++, i += 2)
            {
                sx = 2*(Width - x);
                /* Set Dest(x,y,Channel) = conj(Dest(sx,sy,Channel)) */
                Dest[i] = Dest[sx + sy];
                Dest[i + 1] = -Dest[sx + sy + 1];
            }
        }
    }
}


/** @brief Project u onto the solution set W_v */
static void Project(float *u, float *Temp, float *ConvTemp,
    fftwf_plan ForwardPlan, fftwf_plan InversePlan, const float *u0,
    const float *Phi, int ScaleFactor,
    int OutputWidth, int OutputHeight, int Padding)
{
    const int OutputNumPixels = OutputWidth*OutputHeight;
    const int OutputNumEl = 3*OutputNumPixels;
    int TransWidth, TransHeight, TransNumPixels;
    int i, x, y, sx, sy, OffsetX, OffsetY, Channel;
    
    
    TransWidth = OutputWidth + 2*Padding*ScaleFactor;
    TransHeight = OutputHeight + 2*Padding*ScaleFactor;
        
    if(TransWidth > 2*OutputWidth)
        TransWidth = 2*OutputWidth;
    if(TransHeight > 2*OutputHeight)
        TransHeight = 2*OutputHeight;
    
    TransNumPixels = TransWidth*TransHeight;
    OffsetX = (TransWidth - OutputWidth)/2;
    OffsetY = (TransHeight - OutputHeight)/2;
    OffsetX -= OffsetX % ScaleFactor;
    OffsetY -= OffsetY % ScaleFactor;
    
    for(Channel = 0, i = 0; Channel < OutputNumEl; Channel += OutputNumPixels)
    {
        for(y = 0; y < TransHeight; y++)
        {
            sy = y - OffsetY;
            
            while(1)
            {
                if(sy < 0)
                    sy = -1 - sy;
                else if(sy >= OutputHeight)
                    sy = 2*OutputHeight - 1 - sy;
                else
                    break;
            }
            
            for(x = 0; x < TransWidth; x++, i++)
            {
                sx = x - OffsetX;
                
                while(1)
                {
                    if(sx < 0)
                        sx = -1 - sx;
                    else if(sx >= OutputWidth)
                        sx = 2*OutputWidth - 1 - sx;
                    else
                        break;
                }
                
                Temp[i] = u[sx + OutputWidth*sy + Channel]
                    - u0[sx + OutputWidth*sy + Channel];
            }
        }
    }

    fftwf_execute(ForwardPlan);
    MirrorSpectra(ConvTemp, TransWidth, TransHeight, 3);
    
    MultiplyByPhiInterleaved(ConvTemp, Phi, TransNumPixels, 3);
    ImageSumAliases(ConvTemp, ScaleFactor, TransWidth, TransHeight, 3);
    MultiplyByPhiInterleaved(ConvTemp, Phi, TransNumPixels, 3);
    
    fftwf_execute(InversePlan);
    
    /* Subtract a halved version of Temp from u */
    Temp += OffsetX + TransWidth*OffsetY;
    
    for(Channel = 0, i = 0; Channel < 3; Channel++)
    {
        for(y = 0; y < OutputHeight; y++)
        {
            for(x = 0; x < OutputWidth; x++, i++)
            {
                u[i] -= Temp[x + TransWidth*y];
            }
        }
        
        Temp += TransNumPixels;
    }
}


static float ComputeDiff(const float *u, const float *uLast,
    int OutputNumEl)
{
    float Temp, Diff = 0;
    int i;
    
    for(i = 0; i < OutputNumEl; i++)
    {
        Temp = u[i] - uLast[i];
        Diff += Temp*Temp;
    }
    
    return sqrt(Diff/OutputNumEl);
}


/**
 * @brief Roussos-Maragos interpolation by tensor-driven diffusion
 *
 * @param u pointer to memory for holding the interpolated image
 * @param OutputWidth, OutputHeight output image dimensions
 * @param Input the input image
 * @param InputWidth, InputHeight input image dimensions
 * @param PsfSigma Gaussian PSF standard deviation
 * @param K parameter for constructing the tensor
 * @param Tol convergence tolerance
 * @param MaxMethodIter maximum number of iterations
 * @param DiffIter number of diffusion iterations per method iteration
 *
 * @return 1 on success, 0 on failure
 *
 * This routine implements the projected tensor-driven diffusion method of
 * Roussos and Maragos.  Following Tschumperle, a tensor T is formed from the
 * eigenspectra of the image structure tensor, and the image is diffused as
 * \f[ \partial_t u = \operatorname{div}(T \nabla u). \f]
 * In this implementation, the fast 5x5 explicit filtering scheme proposed by
 * Weickert and Scharr is used to evolve the diffusion.  The diffusion is
 * orthogonally projected onto the feasible set (the set of images for which
 * convolution with the PSF followed by downsampling recovers the input image).
 * Projection is done in the Fourier domain, this is the main computational
 * bottleneck of the method (approximately 60% of the run time is spent in DFT
 * transforms).
 *
 * Beware that this routine is relatively computationally intense, requiring
 * around 2 to 20 seconds for outputs of typical sizes.  Multithreading is
 * applied in some computations if compiling with OpenMP.  Multithreading
 * appears to have little effect for smaller images but reduces run time by
 * about 25% for larger images.
 */
int RoussosInterp(float *u, int OutputWidth, int OutputHeight,
    const float *Input, int InputWidth, int InputHeight, double PsfSigma,
    double K, float Tol, int MaxMethodIter, int DiffIter)
{
    const int Padding = 5;
    const int OutputNumPixels = OutputWidth*OutputHeight;
    const int OutputNumEl = 3*OutputNumPixels;
    float *u0 = NULL, *Txx = NULL, *Txy = NULL, *Tyy = NULL, *Temp = NULL,
        *ConvTemp = NULL, *vx = NULL, *vy = NULL, *Phi = NULL, *uLast = NULL;
    fftwf_plan ForwardPlan = 0, InversePlan = 0;
    fftw_iodim Dims[2];
    fftw_iodim HowManyDims[1];
    filter PreSmooth = {NULL, 0, 0}, PostSmooth = {NULL, 0, 0};
    float PreSmoothSigma, PostSmoothSigma, Diff;
    unsigned long StartTime, StopTime;
    int TransWidth, TransHeight, TransNumPixels;
    int Iter, ScaleFactor, Success = 0;

    
    ScaleFactor = OutputWidth / InputWidth;
    PreSmoothSigma = 0.3f * ScaleFactor;
    PostSmoothSigma = 0.4f * ScaleFactor;
    
    TransWidth = OutputWidth + 2*Padding*ScaleFactor;
    TransHeight = OutputHeight + 2*Padding*ScaleFactor;
    
    if(TransWidth > 2*OutputWidth)
        TransWidth = 2*OutputWidth;
    if(TransHeight > 2*OutputHeight)
        TransHeight = 2*OutputHeight;
    
    TransNumPixels = TransWidth*TransHeight;
    
    printf("Initial interpolation\n");
    
    if(!FourierScale2d(u, OutputWidth, 0, OutputHeight, 0,
        Input, InputWidth, InputHeight, 3, PsfSigma, BOUNDARY_HSYMMETRIC))
        goto Catch;
    else if(MaxMethodIter <= 0)
    {
        Success = 1;
        goto Catch;
    }
    
    /* Allocate a fantastic amount of memory */
    if(ScaleFactor <= 1
        || OutputWidth != ScaleFactor*InputWidth
        || OutputHeight != ScaleFactor*InputHeight
        || !(ConvTemp = (float *)fftwf_malloc(sizeof(float)*24*OutputNumPixels))
        || !(Temp = (float *)fftwf_malloc(sizeof(float)*12*OutputNumPixels))
        || !(Phi = (float *)Malloc(sizeof(float)*4*OutputNumPixels))
        || !(u0 = (float *)Malloc(sizeof(float)*3*OutputNumPixels))
        || !(uLast = (float *)Malloc(sizeof(float)*3*OutputNumPixels))
        || !(Txx = (float *)Malloc(sizeof(float)*OutputNumPixels))
        || !(Txy = (float *)Malloc(sizeof(float)*OutputNumPixels))
        || !(Tyy = (float *)Malloc(sizeof(float)*OutputNumPixels))
        || !(vx = (float *)Malloc(sizeof(float)*OutputNumPixels))
        || !(vy = (float *)Malloc(sizeof(float)*OutputNumPixels))
        || IsNullFilter(PreSmooth = GaussianFilter(PreSmoothSigma,
            (int)ceil(2.5*PreSmoothSigma)))
        || IsNullFilter(PostSmooth = GaussianFilter(PostSmoothSigma,
            (int)ceil(2.5*PostSmoothSigma))))
        goto Catch;

    /* All arrays in the main computation are in planar order so that data
    access in convolutions and DFTs are more localized. */

    HowManyDims[0].n = 3;
    HowManyDims[0].is = TransNumPixels;
    HowManyDims[0].os = TransNumPixels;
    
    Dims[0].n = TransHeight;
    Dims[0].is = TransWidth;
    Dims[0].os = TransWidth;
    Dims[1].n = TransWidth;
    Dims[1].is = 1;
    Dims[1].os = 1;
    
    /* Create plans for the 2D DFT of Src (vectorized over channels).
     * After applying the forward transform,
     * Dest[2*(x + Width*(y + Height*k))] = real component of (x,y,k)th
     * Dest[2*(x + Width*(y + Height*k)) + 1] = imag component of (x,y,k)th
     * where for 0 <= x < Width/2 + 1, 0 <= y < Height.
     */
    if(!(ForwardPlan = fftwf_plan_guru_dft_r2c(2, Dims, 1, HowManyDims,
        Temp, (fftwf_complex *)ConvTemp, FFTW_ESTIMATE | FFTW_DESTROY_INPUT))
        ||  !(InversePlan = fftwf_plan_guru_dft_c2r(2, Dims, 1, HowManyDims,
        (fftwf_complex *)ConvTemp, Temp, FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        goto Catch;
    
    printf("Roussos-Maragos interpolation\n");
    StartTime = Clock();
        
    MakePhi(Phi, PsfSigma, ScaleFactor, TransWidth, TransHeight);
    memcpy(u0, u, sizeof(float)*3*OutputNumPixels);

    /* Projected tensor-driven diffusion main loop */
    for(Iter = 1; Iter <= MaxMethodIter; Iter++)
    {
        memcpy(uLast, u, sizeof(float)*OutputNumEl);
        
        ComputeTensor(Txx, Txy, Tyy, Temp, ConvTemp, u,
            OutputWidth, OutputHeight, PreSmooth, PostSmooth, K);
        
        DiffuseWithTensor(u, vx, vy, Temp, Temp + OutputNumPixels, ConvTemp,
            Txx, Txy, Tyy, OutputWidth, OutputHeight, DiffIter);
    
        Project(u, Temp, ConvTemp, ForwardPlan, InversePlan, u0, Phi,
            ScaleFactor, OutputWidth, OutputHeight, Padding);
        
        Diff = ComputeDiff(u, uLast, OutputNumEl);
        
        if(Iter >= 2 && Diff <= Tol)
        {
            printf("Converged in %d iterations.\n", Iter);
            break;
        }
    }
    
    StopTime = Clock();
    
    if(Diff > Tol)
        printf("Maximum number of iterations exceeded.\n");
        
    printf("CPU Time: %.3f s\n\n", 0.001*(StopTime - StartTime));
    Success = 1;
Catch:
    fftwf_destroy_plan(InversePlan);
    fftwf_destroy_plan(ForwardPlan);
    fftwf_cleanup();
    Free(PostSmooth.Coeff);
    Free(PreSmooth.Coeff);
    Free(vy);
    Free(vx);
    Free(Tyy);
    Free(Txy);
    Free(Txx);
    Free(uLast);
    Free(u0);
    Free(Phi);
    
    if(Temp)
        fftwf_free(Temp);
    if(ConvTemp)
        fftwf_free(ConvTemp);
    return Success;
}
