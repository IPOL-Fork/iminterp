/**
 * @file imnoise.c
 * @brief Corrupt an image with Gaussian, Laplace, or Poisson noise
 * @author Pascal Getreuer <getreuer@gmail.com>
 * 
 * Copyright (c) 2012, Pascal Getreuer
 * All rights reserved.
 * 
 * This program is free software: you can use, modify and/or 
 * redistribute it under the terms of the simplified BSD License. You 
 * should have received a copy of this license along this program. If 
 * not, see <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <stdio.h>
#include <string.h>
#include <ipol/imageio.h>
#include "randmt.h"

/** @brief Display intensities in the range [0,DISPLAY_SCALING] */
#define DISPLAY_SCALING             255

/** @brief Quality for writing JPEG images */
#define JPEGQUALITY                 95


void PrintHelpMessage()
{
    puts("Image noise simulator, P. Getreuer, 2012\n\n"
        "Syntax: imintnoise <model>:<sigma> <input> <output>\n");
        puts("where <input> and <output> are " 
    READIMAGE_FORMATS_SUPPORTED " images.\n");        
        puts(
    "The <model> argument denotes the noise model.  The parameter <sigma>\n"
    "is the noise level, which is defined to be the square root of the\n"
    "expected mean squared error.\n");
        puts(
    "The pixel intensities are denoted below by X[n] and Y[n], and they\n"
    "are scaled as values between 0 and 255.  Values of Y[n] outside of\n"
    "this range are saturated.\n");
        puts(
    "  gaussian:<sigma>  Additive white Gaussian noise\n"
    "                    Y[n] ~ Normal(X[n], sigma^2)\n"
    "                    p(Y[n]|X[n]) = exp( -|Y[n] - X[n]|^2/(2 sigma^2) )\n");
        puts(
    "  laplace:<sigma>   Laplace noise\n"
    "                    Y[n] ~ Laplace(X[n], sigma/sqrt(2))\n"
    "                    p(Y[n]|X[n]) = exp( -|Y[n] - X[n]| sqrt(2)/sigma )\n");
        puts(
    "  poisson:<sigma>   Poisson noise\n"
    "                    Y[n] ~ Poisson(X[n]/a) a\n"
    "                    where a = 255 sigma^2 / (mean value of X)\n");
        puts("Example:\n"
            "  imintnoise laplace:10 input.bmp noisy.bmp\n");
}

void GaussianNoise(float *Image, long NumEl, float Sigma);
void LaplaceNoise(float *Image, long NumEl, float Sigma);
void PoissonNoise(float *Image, long NumEl, float Sigma);
int IsGrayscale(const float *Image, long NumPixels);


int main(int argc, char **argv)
{
    const char *Model, *InputFile, *OutputFile;
    char *ParamString;
    float *Image;
    float Param;
    long NumPixels, NumEl;
    int Width, Height, NumChannels, Status = 1;
    
    if(argc != 4 || !(ParamString = strchr(argv[1], ':')))
    {
        PrintHelpMessage();
        return 0;
    }
    
    *ParamString = '\0';
    
    /* Read command line arguments */
    Model = argv[1];
    Param = (float)(atof(ParamString + 1) / DISPLAY_SCALING);
    InputFile = argv[2];
    OutputFile = argv[3];
    
    /* Read the input image */
    if(!(Image = (float *)ReadImage(&Width, &Height, InputFile, 
        IMAGEIO_RGB | IMAGEIO_PLANAR | IMAGEIO_FLOAT)))
        goto Catch;
    
    NumPixels = ((long)Width) * ((long)Height);
    NumChannels = (IsGrayscale(Image, NumPixels)) ? 1 : 3;
    NumEl = NumChannels * NumPixels;
    
    /* Initialize random number generator */
    init_randmt_auto();
    
    if(!strcmp(Model, "gaussian"))
        GaussianNoise(Image, NumEl, Param);
    else if(!strcmp(Model, "laplace"))
        LaplaceNoise(Image, NumEl, Param);
    else if(!strcmp(Model, "poisson"))
        PoissonNoise(Image, NumEl, Param);    
    else
    {
        fprintf(stderr, "Unknown noise model, \"%s\".\n", Model);
        goto Catch;
    }
    
    /* Write noisy and denoised images */
    if(!WriteImage(Image, Width, Height, OutputFile, 
        ((NumChannels == 1) ? IMAGEIO_GRAYSCALE : IMAGEIO_RGB)
        | IMAGEIO_PLANAR | IMAGEIO_FLOAT, JPEGQUALITY))
    {
        fprintf(stderr, "Error writing to \"%s\".\n", OutputFile);
        goto Catch;
    }
    
    Status = 0;
Catch:
    Free(Image);
    return Status;
}


void GaussianNoise(float *Image, long NumEl, float Sigma)
{
    long n;
    
    for(n = 0; n < NumEl; n++)
        Image[n] += (float)(Sigma*rand_normal());
}


void LaplaceNoise(float *Image, long NumEl, float Sigma)
{
    const float Mu = (float)(M_1_SQRT2 * Sigma);
    long n;
    
    for(n = 0; n < NumEl; n++)
        Image[n] += (float)(rand_exp(Mu) * ((rand_unif() < 0.5) ? -1 : 1));
}
 
 
void PoissonNoise(float *Image, long NumEl, float Sigma)
{
    double a, Mean = 0;
    long n;
    
    for(n = 0; n < NumEl; n++)
        Mean += Image[n];
    
    Mean /= NumEl;
    a = Sigma * Sigma / ((Mean > 0) ? Mean : (0.5/255));
    
    for(n = 0; n < NumEl; n++)
        Image[n] = (float)(rand_poisson(Image[n] / a) * a);
}


int IsGrayscale(const float *Image, long NumPixels)
{
    const float *Red = Image;
    const float *Green = Image + NumPixels;
    const float *Blue = Image + 2*NumPixels;
    long n;
    
    for(n = 0; n < NumPixels; n++)
        if(Red[n] != Green[n] || Red[n] != Blue[n])
            return 0;
    
    return 1;
}
