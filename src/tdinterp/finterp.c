/**
 * @file finterp.c
 * @brief Fourier zero-padding interpolation
 * @author Pascal Getreuer <getreuer@gmail.com>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

#include "basic.h"
#include "finterp.h"


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
* @brief Boundary handling function for half-sample symmetric extension
* @param N is the data length
* @param i is an index into the data
* @return an index that is always between 0 and N - 1
*/
static int HSymExtension(int N, int i)
{
    while(1)
    {
        if(i < 0)
            i = -1 - i;
        else if(i >= N)
            i = (2*N - 1) - i;
        else
            return i;
    }
}


/**
* @brief Boundary handling function for whole-sample symmetric extension
* @param N is the data length
* @param i is an index into the data
* @return an index that is always between 0 and N - 1
*/
static int WSymExtension(int N, int i)
{
    while(1)
    {
        if(i < 0)
            i = -i;
        else if(i >= N)
            i = (2*N - 2) - i;
        else
            return i;
    }
}


int (*ExtensionMethod[3])(int, int) =
    {ConstExtension, HSymExtension, WSymExtension};

    
static int FourierScaleScan(float *Dest,
    int DestStride, int DestScanStride, int DestChannelStride, int DestScanSize,
    const float *Src,
    int SrcStride, int SrcScanStride, int SrcChannelStride, int SrcScanSize,
    int NumScans, int NumChannels, float XStart, double PsfSigma,
    boundaryhandling Boundary)
{
    const int SrcPadScanSize = 2*(SrcScanSize
        - ((Boundary == BOUNDARY_WSYMMETRIC) ? 1:0));
    const int DestPadScanSize = (Boundary == BOUNDARY_HSYMMETRIC) ?
        (2*DestScanSize)
        : (2*DestScanSize - (int)floor((2.0f*DestScanSize)/SrcScanSize + 0.5f));
    const int ReflectOffset = SrcPadScanSize
        - ((Boundary == BOUNDARY_HSYMMETRIC) ? 1:0);
    const int SrcDftSize = SrcPadScanSize/2 + 1;
    const int DestDftSize = DestPadScanSize/2 + 1;
    const int BufSpatialNumEl = DestPadScanSize*NumScans*NumChannels;
    const int BufDftNumEl = 2*DestDftSize*NumScans*NumChannels;
    float *BufSpatial = NULL, *BufDft = NULL, *Modulation = NULL, *Ptr;
    fftwf_plan Plan = 0;
    fftwf_iodim Dims[1], HowManyDims[1];
    float Temp, Denom;
    int i, Scan, Channel, Success = 0;
    
    
    if((Boundary != BOUNDARY_HSYMMETRIC && Boundary != BOUNDARY_WSYMMETRIC)
        || !(BufSpatial = (float *)fftwf_malloc(sizeof(float)*BufSpatialNumEl))
        || !(BufDft = (float *)fftwf_malloc(sizeof(float)*BufDftNumEl)))
        goto Catch;
    
    if(XStart != 0)
    {
        if(!(Modulation = (float *)Malloc(sizeof(float)*2*DestDftSize)))
            goto Catch;
        
        for(i = 0; i < DestDftSize; i++)
        {
            Temp = (float)(M_2PI*XStart*i/SrcPadScanSize);
            Modulation[2*i + 0] = (float)cos(Temp);
            Modulation[2*i + 1] = (float)sin(Temp);
        }
    }

    /* Fill BufSpatial with the input and symmetrize */
    for(Channel = 0; Channel < NumChannels; Channel++)
    {
        for(Scan = 0; Scan < NumScans; Scan++)
        {
            for(i = 0; i < SrcScanSize; i++)
                BufSpatial[i + SrcPadScanSize*(Scan + NumScans*Channel)]
                    = Src[SrcStride*i + SrcScanStride*Scan + SrcChannelStride*Channel];

            for(; i < SrcPadScanSize; i++)
                BufSpatial[i + SrcPadScanSize*(Scan + NumScans*Channel)]
                    = Src[SrcStride*(ReflectOffset - i)
                    + SrcScanStride*Scan + SrcChannelStride*Channel];
        }
    }
    
    /* Initialize DFT buffer to zeros (there is no "fftwf_calloc").  Note that
    it is not safely portable to use memset for this purpose.
    http://c-faq.com/malloc/calloc.html  */
    for(i = 0; i < BufDftNumEl; i++)
        BufDft[i] = 0.0f;

    /* Perform DFT real-to-complex transform */
    Dims[0].n = SrcPadScanSize;
    Dims[0].is = 1;
    Dims[0].os = 1;
    HowManyDims[0].n = NumScans*NumChannels;
    HowManyDims[0].is = SrcPadScanSize;
    HowManyDims[0].os = DestDftSize;

    if(!(Plan = fftwf_plan_guru_dft_r2c(1, Dims, 1, HowManyDims, BufSpatial,
        (fftwf_complex *)BufDft, FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        goto Catch;

    fftwf_execute(Plan);
    fftwf_destroy_plan(Plan);
    
    if(PsfSigma == 0)
        for(Channel = 0, Ptr = BufDft; Channel < NumChannels; Channel++)
            for(Scan = 0; Scan < NumScans; Scan++, Ptr += 2*DestDftSize)
                for(i = 0; i < SrcDftSize; i++)
                {
                    Ptr[2*i + 0] /= SrcPadScanSize;
                    Ptr[2*i + 1] /= SrcPadScanSize;
                }
    else
    {
        /* Also divide by the Gaussian point spread function in this case */
        Temp = (float)(SrcPadScanSize / (M_2PI * PsfSigma));
        Temp = 2*Temp*Temp;

        for(i = 0; i < SrcDftSize; i++)
        {
            if(i <= DestScanSize)
                Denom = (float)exp(-(i*i)/Temp);
            else
                Denom = (float)exp(
					-((DestPadScanSize - i)*(DestPadScanSize - i))/Temp);

            Denom *= SrcPadScanSize;

            for(Channel = 0; Channel < NumChannels; Channel++)
                for(Scan = 0; Scan < NumScans; Scan++)
                {
                    BufDft[2*(i + DestDftSize*(Scan + NumScans*Channel)) + 0] /= Denom;
                    BufDft[2*(i + DestDftSize*(Scan + NumScans*Channel)) + 1] /= Denom;
                }
        }
    }
    
    /* If XStart is nonzero, modulate the DFT to translate the result */
    if(XStart != 0)
        for(Channel = 0, Ptr = BufDft; Channel < NumChannels; Channel++)
            for(Scan = 0; Scan < NumScans; Scan++, Ptr += 2*DestDftSize)
                for(i = 0; i < SrcDftSize; i++)
                {
                    /* Complex multiply */
                    Temp = Ptr[2*i + 0]*Modulation[2*i + 0]
                        - Ptr[2*i + 1]*Modulation[2*i + 1];
                    Ptr[2*i + 1] = Ptr[2*i + 0]*Modulation[2*i + 1]
                        + Ptr[2*i + 1]*Modulation[2*i + 0];
                    Ptr[2*i + 0] = Temp;
                }
    
    /* Perform inverse DFT complex-to-real transform */
    Dims[0].n = DestPadScanSize;
    Dims[0].is = 1;
    Dims[0].os = 1;
    HowManyDims[0].n = NumScans*NumChannels;
    HowManyDims[0].is = DestDftSize;
    HowManyDims[0].os = DestPadScanSize;
    
    if(!(Plan = fftwf_plan_guru_dft_c2r(1, Dims, 1, HowManyDims,
        (fftwf_complex *)BufDft, BufSpatial,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT)))
        goto Catch;
    
    fftwf_execute(Plan);
    fftwf_destroy_plan(Plan);

    /* Fill Dest with the result (and trim padding) */
    for(Channel = 0; Channel < NumChannels; Channel++)
    {
        for(Scan = 0; Scan < NumScans; Scan++)
        {
            for(i = 0; i < DestScanSize; i++)
                Dest[DestStride*i + DestScanStride*Scan + DestChannelStride*Channel]
                    = BufSpatial[i + DestPadScanSize*(Scan + NumScans*Channel)];
        }
    }

    Success = 1;
Catch:
    Free(Modulation);
    if(BufDft)
        fftwf_free(BufDft);
    if(BufSpatial)
        fftwf_free(BufSpatial);
    fftwf_cleanup();
    return Success;
}


/**
 * @brief Scale image with Fourier zero padding
 *
 * @param Dest pointer to memory for holding the interpolated image
 * @param DestWidth output image width
 * @param XStart leftmost sample location (in input coordinates)
 * @param DestHeight output image height
 * @param YStart uppermost sample location (in input coordinates)
 * @param Src the input image
 * @param SrcWidth, SrcHeight, NumChannels input image dimensions
 * @param PsfSigma Gaussian PSF standard deviation
 * @param Boundary boundary handling
 * @return 1 on success, 0 on failure.
 *
 * The image is first mirror folded with half-sample even symmetry to avoid
 * boundary artifacts, then transformed with a real-to-complex DFT.
 *
 * The interpolation is computed so that Dest[m + DestWidth*n] is the
 * interpolation of Input at sampling location
 *    (XStart + m*SrcWidth/DestWidth, YStart + n*SrcHeight/DestHeight)
 * for m = 0, ..., DestWidth - 1, n = 0, ..., DestHeight - 1, where the
 * pixels of Src are located at the integers.
 */
int FourierScale2d(float *Dest, int DestWidth, float XStart,
    int DestHeight, float YStart,
    const float *Src, int SrcWidth, int SrcHeight, int NumChannels,
    double PsfSigma, boundaryhandling Boundary)
{
    float *Buf = NULL;
    unsigned long StartTime, StopTime;
    int Success = 0;
        
    
    if(!Dest || DestWidth < SrcWidth || DestHeight < SrcHeight || !Src
        || SrcWidth <= 0 || SrcHeight <= 0 || NumChannels <= 0 || PsfSigma < 0
        || !(Buf = (float *)Malloc(sizeof(float)*SrcWidth*DestHeight*3)))
        return 0;

    StartTime = Clock();
    
    /* Scale the image vertically */
    if(!FourierScaleScan(Buf, SrcWidth, 1, SrcWidth*DestHeight, DestHeight,
        Src, SrcWidth, 1, SrcWidth*SrcHeight, SrcHeight,
        SrcWidth, 3, YStart, PsfSigma, Boundary))
        goto Catch;
    
    /* Scale the image horizontally */
    if(!FourierScaleScan(Dest, 1, DestWidth, DestWidth*DestHeight, DestWidth,
        Buf, 1, SrcWidth, SrcWidth*DestHeight, SrcWidth,
        DestHeight, 3, XStart, PsfSigma, Boundary))
        goto Catch;

    StopTime = Clock();
    printf("CPU Time: %.3f s\n\n", 0.001*(StopTime - StartTime));
    
    Success = 1;
Catch:
    Free(Buf);
    return Success;
}
