/**
 * @file fitsten.c
 * @brief Select the best-fitting contour stencils
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

#include <stdlib.h>
#include "fitsten.h"


#define PIXEL_STRIDE	3

/* Constants related to stencil edge weights */
/* mu = 3/4 * sqrt(2) */
#define W_MU        (1.060660171779821)
/* Divisor for axial stencils = 6 */
#define W_PI_2      (6.0)
/* Divisor for diagonal stencils = 4 sqrt(2) */
#define W_PI_4      (5.656854249492381)
/* Divisor for pi/8 stencils = 6(1 + sqrt(2))sqrt(2 - sqrt(2)) + 12 */
#define W_PI_8      (23.086554390135440)



/** @brief Compute the distance between to color pixels in YCbCr space */
static int Dist(int32_t *A, int32_t *B)
{
    int iDistR = A[0] - B[0];
    int iDistG = A[1] - B[1];
    int iDistB = A[2] - B[2];
    register int iDistY = 2*iDistR + 4*iDistG + iDistB;
    register int iDistU = -iDistR - 2*iDistG + 3*iDistB;
    register int iDistV = 4*iDistR - 3*iDistG - iDistB;
    return (abs(iDistY) + abs(iDistU) + abs(iDistV)) >> 4;
}


/**
* @brief Select the best-fitting contour stencils
*
* @param Stencil pointer to array of size Width by Height to be filled with
*        the best-fitting stencil for each pixels
* @param Image the input RGB image
* @param Width, Height image dimensions
* @param StencilMul multiply the stencil index by this factor
*
* \c FitStencils finds the best-fitting stencil at each pixel of the input
* image
* \f[ \mathcal{S}^\star(k) = \arg\min_{\mathcal{S}\in\Sigma}
*     \tfrac{1}{|\mathcal{S}|} (\mathcal{S} \star [v])(k) \f]
* to estimate the local contour orientation.  The result \c Stencil is a
* two-dimensional array with the index (0, 1, ..., or 7) of the best-fitting
* stencil for each pixel.  The stencil TV estimates
* \f$(\mathcal{S} \star [v])(k)\f$ are filtered by [1,2,1; 2,4,2; 1,2,1]/16
* to improve reliability.
*
* For efficiency in \c CWFirstPass and \c CWCorrectionPass, the stencil
* indices are multiplied by \c StencilMul (so that the index may be used
* directly as an offset into the samples table).
*/
int FitStencils(int *Stencil, int32_t *Image, int Width, int Height, int StencilMul)
{
    const int Stride = PIXEL_STRIDE*Width;
    const int TVStride = 8*Width;
    int *StencilTv; /* Contour stencil total variations estimates TV[S]            */
    int Dh[3][2];   /* Horizontal differences computed over the current 3x3 window */
    int Dv[2][3];   /* Vertical differences                                        */
    int Da[2][2];   /* Diagonal differences |Image(m,n) - Image(m+1,n+1)|          */
    int Db[2][2];   /* Diagonal differences |Image(m+1,n) - Image(m,n+1)|          */
    int *TvPtr, TvMin, Temp;
    int x, y, k, S;

    
    if(!(StencilTv = (int *)malloc(sizeof(int)*8*Width*Height)))
        return 0;
    
    TvPtr = StencilTv;
    
    for(y = 0; y < Height; y++)
    {
        Dh[0][1] = 0;
        Dh[1][1] = 0;
        Dh[2][1] = 0;
        
        Dv[0][1] = 0;
        Dv[1][1] = 0;
        Dv[0][2] = (y > 0) ? Dist(Image - Stride, Image) : 0;
        Dv[1][2] = (y < Height-1) ? Dist(Image + Stride, Image) : 0;
        
        Da[0][1] = 0;
        Da[1][1] = 0;
        Db[0][1] = 0;
        Db[1][1] = 0;
    
        for(x = 0; x < Width; x++, Image += PIXEL_STRIDE, TvPtr += 8)
        {
            Dh[0][0] = Dh[0][1];
            Dh[1][0] = Dh[1][1];
            Dh[2][0] = Dh[2][1];
            
            Dv[0][0] = Dv[0][1];
            Dv[1][0] = Dv[1][1];
            Dv[0][1] = Dv[0][2];
            Dv[1][1] = Dv[1][2];
            
            Da[0][0] = Da[0][1];
            Da[1][0] = Da[1][1];
            
            Db[0][0] = Db[0][1];
            Db[1][0] = Db[1][1];
            
            if(x < Width-1)
            {
                Dh[1][1] = Dist(Image, Image+PIXEL_STRIDE);
            
                if(y > 0)
                {
                    Dh[0][1] = Dist(Image-Stride, Image-Stride+PIXEL_STRIDE);
                    Dv[0][2] = Dist(Image-Stride+PIXEL_STRIDE, Image+PIXEL_STRIDE);
                    Da[0][1] = Dist(Image-Stride, Image+PIXEL_STRIDE);
                    Db[0][1] = Dist(Image-Stride+PIXEL_STRIDE, Image);
                }
                
                if(y < Height-1)
                {
                    Dh[2][1] = Dist(Image+Stride, Image+Stride+PIXEL_STRIDE);
                    Dv[1][2] = Dist(Image+PIXEL_STRIDE, Image+Stride+PIXEL_STRIDE);
                    Da[1][1] = Dist(Image, Image+Stride+PIXEL_STRIDE);
                    Db[1][1] = Dist(Image+PIXEL_STRIDE, Image+Stride);
                }
            }
            else
            {
                Dh[0][1] = 0;
                Dh[1][1] = 0;
                Dh[2][1] = 0;
                
                Dv[0][2] = 0;
                Dv[1][2] = 0;
                
                Da[0][1] = 0;
                Da[1][1] = 0;
                
                Db[0][1] = 0;
                Db[1][1] = 0;
            }

            TvPtr[1] = (Db[1][0] + Db[0][1]) + Db[0][0] + Db[1][1];
            TvPtr[3] = (Dh[1][0] + Dh[1][1]) + Dh[0][0] + Dh[0][1] + Dh[2][0] + Dh[2][1];
            TvPtr[5] = (Da[0][0] + Da[1][1]) + Da[1][0] + Da[0][1];
            TvPtr[7] = (Dv[0][1] + Dv[1][1]) + Dv[0][0] + Dv[1][0] + Dv[0][2] + Dv[1][2];
        
            TvPtr[0] = (int)((TvPtr[7] + TvPtr[1]*W_MU) / W_PI_8 + 0.5);
            TvPtr[2] = (int)((TvPtr[1]*W_MU + TvPtr[3]) / W_PI_8 + 0.5);
            TvPtr[4] = (int)((TvPtr[3] + TvPtr[5]*W_MU) / W_PI_8 + 0.5);
            TvPtr[6] = (int)((TvPtr[5]*W_MU + TvPtr[7]) / W_PI_8 + 0.5);
            
            TvPtr[1] = (int)(TvPtr[1] / W_PI_4 + 0.5);
            TvPtr[3] = (int)(TvPtr[3] / W_PI_2 + 0.5);
            TvPtr[5] = (int)(TvPtr[5] / W_PI_4 + 0.5);
            TvPtr[7] = (int)(TvPtr[7] / W_PI_2 + 0.5);
        }
    }
    
    TvPtr = StencilTv;
    
    for(y = 0; y < Height; y++)
    {
        for(x = 0; x < Width; x++, TvPtr += 8, Stencil++)
        {
            if(1 <= x && x < Width - 1 && 1 <= y && y < Height - 1)
            {
                /* Filter the TV estimates by [1,2,1; 2,4,2; 1,2,1] */
                TvMin = TvPtr[-TVStride - 8] + TvPtr[-TVStride + 8]
                        + TvPtr[TVStride - 8] + TvPtr[TVStride + 8]
                        + ((TvPtr[-TVStride] + TvPtr[-8]
                        + TvPtr[8] + TvPtr[TVStride]) << 1)
                        + (TvPtr[0] << 2);
                S = 0;
            
                for(k = 1; k < 8; k++)
                {
                    Temp = TvPtr[k - TVStride - 8] + TvPtr[k - TVStride + 8]
                        + TvPtr[k + TVStride - 8] + TvPtr[k + TVStride + 8]
                        + ((TvPtr[k - TVStride] + TvPtr[k - 8]
                        + TvPtr[k + 8] + TvPtr[k + TVStride]) << 1)
                        + (TvPtr[k] << 2);
                    
                    if(Temp < TvMin)
                    {
                        TvMin = Temp;
                        S = k;
                    }
                        
                }
            }
            else
            {
                TvMin = TvPtr[0];
                S = 0;
                
                for(k = 1; k < 8; k++)
                    if(TvPtr[k] < TvMin)
                        TvMin = TvPtr[S = k];
            }
            
            *Stencil = StencilMul * S;
        }
    }
    
    free(StencilTv);
    return 1;
}
