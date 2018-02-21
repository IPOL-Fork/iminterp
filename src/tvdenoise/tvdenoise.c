/**
 * @file tvdenoise.c
 * @brief Total variation regularized denoising demo for IPOL
 * @author Pascal Getreuer <getreuer@gmail.com>
 * 
 * 
 * Copyright (c) 2011-2012, Pascal Getreuer
 * All rights reserved.
 * 
 * This program is free software: you can use, modify and/or 
 * redistribute it under the terms of the simplified BSD License. You 
 * should have received a copy of this license along this program. If 
 * not, see <http://www.opensource.org/licenses/bsd-license.html>.
 */

/**
 * @mainpage
 * @verbinclude readme.txt
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "num.h"
#include "tvreg.h"
#include <ipol/imageio.h>

/** @brief Display intensities in the range [0,DISPLAY_SCALING] */
#define DISPLAY_SCALING             255

/** @brief Number of iterations for tuning lambda */
#define LAMBDA_TUNE_ITERATIONS      5

/** @brief Quality for writing JPEG images */
#define JPEGQUALITY                 95

#ifdef NUM_SINGLE
#define IMAGEIO_NUM           (IMAGEIO_SINGLE)
#else
#define IMAGEIO_NUM           (IMAGEIO_DOUBLE)
#endif

/** @brief struct representing an image */
typedef struct
{
    /** @brief Float image data */
    num *Data;
    /** @brief Image width */
    int Width;
    /** @brief Image height */
    int Height;
    /** @brief Number of channels */
    int NumChannels;
} image;


/** @brief struct of program parameters */
typedef struct
{
    /** @brief Input file (noisy) */
    char *InputFile;
    /** @brief Output file (denoised) */
    char *OutputFile;    
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;
    
    /** @brief Noise model */
    char *Model;
    /** @brief Noise standard deviation */
    num Sigma;
    /** @brief Fidelity strength */
    num Lambda;
} programparams;  


/** @brief Print program explanation and usage */
void PrintHelpMessage()
{
    puts(
    "Total variation regularized denoising IPOL demo, P. Getreuer, 2012\n\n"
    "Syntax: iminttvdenoise [options] <noisy> <denoised>\n");
        puts("where <noisy> and <denoised> are " 
    READIMAGE_FORMATS_SUPPORTED " images.\n");       
    puts(
    "Either lambda (the fidelty strength) or sigma (the noise standard\n"
    "deviation) should be specified.  If sigma is specified, then lambda is\n"
    "selected automatically using Chambolle's algorithm.\n");        
    puts("Options:");
    puts("  -n <model>          Specify noise model, where <model> is\n"
    "                      gaussian  Additive white Gaussian noise\n"
    "                                Y[n] ~ Normal(X[n], sigma^2)");
    puts("                      laplace   Laplace noise\n"    
    "                                Y[n] ~ Laplace(X[n], sigma/sqrt(2))\n"
    "                      poisson   Poisson noise\n"
    "                                Y[n] ~ Poisson(X[n]/a) a\n"
    "                                where a = 255 sigma^2 / (mean X)");
    puts("  -n <model>:<sigma>  Specify sigma, the noise standard deviation");
    puts("  -l <number>         Specify lambda, the fidelity stength\n");
#ifdef USE_LIBJPEG
    puts("  -q <number>         Quality for saving JPEG images (0 to 100)\n");
#endif    
    puts("Example:\n"
        "  iminttvdenoise -n laplace:10 noisy.bmp denoised.bmp\n");
}

int Denoise(image u, image f, const char *Model, num Sigma, num Lambda);
num ComputeRmse(image f, image u);
int IsGrayscale(image f);
static int ParseParams(programparams *Param, int argc, char *argv[]);


int main(int argc, char **argv)
{
    programparams Param;
    image f = {NULL, 0, 0, 0}, u = {NULL, 0, 0, 0};
    int Status = 1;
    
    /* Read command line arguments */
    if(!ParseParams(&Param, argc, argv))
        return 0;
    
    /* Read the input image */
    if(!(f.Data = (num *)ReadImage(&f.Width, &f.Height, Param.InputFile,
        IMAGEIO_RGB | IMAGEIO_PLANAR | IMAGEIO_NUM)))
        goto Catch;
    
    f.NumChannels = IsGrayscale(f) ? 1 : 3;
    u = f;
    
    /* Allocate space for the noisy and denoised image */
    if(!(u.Data = (num *)Malloc(sizeof(num) * ((size_t)f.Width)
        * ((size_t)f.Height) * f.NumChannels)))
        goto Catch;
    
    /* Denoise the image */
    if(!Denoise(u, f, Param.Model, Param.Sigma, Param.Lambda))    
        goto Catch;
    
    /* Write the denoised image */
    if(!WriteImage(u.Data, u.Width, u.Height, Param.OutputFile, 
        ((u.NumChannels == 1) ? IMAGEIO_GRAYSCALE : IMAGEIO_RGB)
        | IMAGEIO_PLANAR | IMAGEIO_NUM, JPEGQUALITY))    
        goto Catch;
    
    Status = 0;
Catch:
    if(u.Data)
        Free(u.Data);
    if(f.Data)
        Free(f.Data);
    return Status;
}


/** 
 * @brief Tune lambda according to the discrepancy principle 
 * @param Opt tvreopt options variable
 * @param u the denoised image
 * @param f the given noisy image
 * @param Model string specifying either "gaussian", "laplace", or "poisson"
 * @param Sigma the standard deviation of the noise
 * @return 1 on success, 0 on failure
 * 
 * This routine tunes lambda according to the discrepancy principle.  
 * Empirical estimates of the optimal lambda value are used as the 
 * initialization.  
 * 
 * The TV denoising computation itself is performed by TvRestore().
 */
int LambdaTune(tvregopt *Opt, image u, image f, const char *Model, num Sigma)
{
    num Lambda, Rmse;
    int k;
    
    /* Empirical estimates of the optimal lambda value 
       (here, Sigma is scaled relative to intensities in [0,1]) */
    if(!strcmp(Model, "gaussian"))
        Lambda = 0.7079 / Sigma + 0.002686 / (Sigma * Sigma);
    else if(!strcmp(Model, "laplace"))
        Lambda = (-0.00416*Sigma + 0.001301) 
            / (((Sigma - 0.2042)*Sigma + 0.01635)*Sigma + 5.836e-4);
    else if(!strcmp(Model, "poisson"))
        Lambda = 0.2839 / Sigma + 0.001502 / (Sigma * Sigma);
    else
    {
        fprintf(stderr, "Unrecognized noise model \"%s\"\n", Model);
        return 0;
    }
    
    if(Lambda < 1e-4)   /* Prevent nonpositive Lambda */
        Lambda = 1e-4;
    
    TvRegSetLambda(Opt, Lambda);
    
    printf("Tuning lambda...\n\n");
    printf("  lambda    distance (target = %.5f)\n", DISPLAY_SCALING * Sigma);
    printf(" --------------------\n");
    printf("  %-9.4f", Lambda);
    
    for(k = 0; k < LAMBDA_TUNE_ITERATIONS; k++)
    {
        /* Each TvRestore uses the current u as the initial guess and 
           overwrites it with the denoising result.  This speeds up the 
           computation because the result using the previous Lambda value
           is a good estimate for the next Lambda value. */
        if(!TvRestore(u.Data, f.Data, f.Width, f.Height, f.NumChannels, Opt))
        {
            fprintf(stderr, "Error in computation.\n");
            return 0;
        }
        
        Rmse = ComputeRmse(f, u);
        
        if(!strcmp(Model, "laplace"))
            Lambda *= sqrt(Rmse/Sigma);
        else
            Lambda *= Rmse/Sigma;
        
        TvRegSetLambda(Opt, Lambda);
        printf(" %.5f\n  %-9.4f", DISPLAY_SCALING * Rmse, Lambda);
    }
    
    return 1;
}


/** 
 * @brief TV regularized denoising 
 * @param u denoised image
 * @param f given noisy image
 * @param Model string specifying either "gaussian", "laplace", or "poisson"
 * @param Sigma the standard deviation of the noise
 * @param Lambda the fidelity stength
 * @return 1 on success, 0 on failure
 */
int Denoise(image u, image f, const char *Model, num Sigma, num Lambda)
{
    tvregopt *Opt = NULL;
    int Success = 0;
    
    if(!(Opt = TvRegNewOpt()))
    {
        fprintf(stderr, "Memory allocation failed\n");
        return 0;
    }
    else if(!(TvRegSetNoiseModel(Opt, Model)))
    {
        fprintf(stderr, "Unrecognized noise model \"%s\"\n", Model);
        return 0;
    }
    
    printf("TV regularized denoising with %c%s noise model\n", 
        toupper(*Model), Model + 1);    
    
    /* Set initial guess as u = f */
    memcpy(u.Data, f.Data, sizeof(num)*f.Width*f.Height*f.NumChannels);
    TvRegSetPlotFun(Opt, NULL, NULL);
    TvRegSetTol(Opt, (num)1e-2);
    TvRegSetMaxIter(Opt, 40);
    
    if(Sigma <= 0)
        TvRegSetLambda(Opt, Lambda);
    else if(!LambdaTune(Opt, u, f, Model, Sigma))
        goto Catch;
    
    /* Final denoising */
    TvRegSetTol(Opt, (num)5e-4);
    TvRegSetMaxIter(Opt, 100);
    
    if(!TvRestore(u.Data, f.Data, f.Width, f.Height, f.NumChannels, Opt))
    {
        fprintf(stderr, "Error in computation.\n");
        goto Catch;
    }
    
    if(Sigma > 0)
        printf(" %.5f\n\n", DISPLAY_SCALING * ComputeRmse(f, u));
    
    Success = 1;
Catch:
    TvRegFreeOpt(Opt);
    return Success;
}


/** @brief Compute the root mean-square-error between two images */
num ComputeRmse(image f, image u)
{
    num Temp, Rmse = 0;
    long n, N = ((long)u.Width) * ((long)u.Height) * u.NumChannels;
    
    for(n = 0; n < N; n++) 
    {
        Temp = f.Data[n] - u.Data[n];
        Rmse += Temp*Temp;
    }
    
    return sqrt(Rmse/N);
}


/** @brief Test whether image is grayscale */
int IsGrayscale(image f)
{
    const long NumPixels = ((long)f.Width) * ((long)f.Height);
    const num *Red = f.Data;
    const num *Green = f.Data + NumPixels;
    const num *Blue = f.Data + 2*NumPixels;
    long n;
    
    for(n = 0; n < NumPixels; n++)
        if(Red[n] != Green[n] || Red[n] != Blue[n])
            return 0;
    
    return 1;
}


static int ParseParams(programparams *Param, int argc, char *argv[])
{
    char *OptionString, *Sigma;
    char OptionChar;
    int i;
    
    if(argc < 2)
    {
        PrintHelpMessage();
        return 0;
    }
    
    /* Set parameter defaults */
    Param->InputFile = NULL;
    Param->OutputFile = NULL;    
    Param->Model = "gaussian";
    Param->Sigma = -1;
    Param->Lambda = -1;    
    Param->JpegQuality = 95;
        
    for(i = 1; i < argc;)
    {
        if(argv[i] && argv[i][0] == '-')
        {
            if((OptionChar = argv[i][1]) == 0)
            {
                ErrorMessage("Invalid parameter format.\n");
                return 0;
            }
            
            if(argv[i][2])
                OptionString = &argv[i][2];
            else if(++i < argc)
                OptionString = argv[i];
            else
            {
                ErrorMessage("Invalid parameter format.\n");
                return 0;
            }
            
            switch(OptionChar)
            {
            case 'n':
                Sigma = strchr(OptionString, ':');
                Param->Model = OptionString;
                
                if((Sigma = strchr(OptionString, ':')))
                {
                    *Sigma = '\0';
                    Param->Sigma = (num)(atof(Sigma + 1) 
                        / DISPLAY_SCALING);
                    
                    if(Param->Sigma <= 0)
                    {
                        ErrorMessage("sigma must be positive.\n");
                        return 0;
                    }
                }
                break;
            case 'l':
                Param->Lambda = (num)atof(OptionString);

                if(Param->Lambda <= 0)
                {
                    ErrorMessage("lambda must be positive.\n");
                    return 0;
                }
                break;
#ifdef LIBJPEG_SUPPORT
            case 'q':
                Param->JpegQuality = atoi(OptionString);

                if(Param->JpegQuality <= 0 || Param->JpegQuality > 100)
                {
                    ErrorMessage("JPEG quality must be between 0 and 100.\n");
                    return 0;
                }
                break;
#endif
            case '-':
                PrintHelpMessage();
                return 0;
            default:
                if(isprint(OptionChar))
                    ErrorMessage("Unknown option \"-%c\".\n", OptionChar);
                else
                    ErrorMessage("Unknown option.\n");

                return 0;
            }
            
            i++;
        }
        else
        {
            if(!Param->InputFile)
                Param->InputFile = argv[i];
            else
                Param->OutputFile = argv[i];
            
            i++;
        }
    }
    
    if(!Param->InputFile)
    {
        PrintHelpMessage();
        return 0;
    }
    
    if(Param->Sigma < 0 && Param->Lambda < 0)
    {
        ErrorMessage("Either sigma or lambda must be specified.\n");
        return 0;
    }
    
    return 1;
}

