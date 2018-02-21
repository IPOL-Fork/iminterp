/**
 * @file sinterpcli.c
 * @brief Contour stencil interpolation command line program
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

/**
 * @mainpage
 * @htmlinclude readme.html
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <ipol/imageio.h>
#include "sinterp.h"


#ifndef DOXYGEN_SHOULD_SKIP_THIS

/* Default program parameters */
#define DEFAULT_PSF_SIGMA           0.35
#define DEFAULT_RHO_SIGMA_TANGENT   1.2
#define DEFAULT_RHO_SIGMA_NORMAL    0.6
#define DEFAULT_REFINEMENT_PASSES   2


/* Colors for drawing contour lines, adapted from the
 * Tango Desktop Project
 * http://tango.freedesktop.org/Tango_Desktop_Project
 */
const float ColorLine[3]   = {0.35f, 0.60f, 0.86f};
const float ColorCurve1[3] = {0.54f, 0.88f, 0.20f};
const float ColorCurve2[3] = {0.98f, 0.91f, 0.22f};
const float ColorCorner[3] = {1.00f, 0.17f, 0.16f};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** @brief struct of program parameters */
typedef struct
{
    /** @brief Input file name */
    const char *InputFile;
    /** @brief Output file name */
    const char *OutputFile;
    /** @brief Background image file name (used only for -s with SVG output) */
    const char *BgFile;
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;
    
    /** @brief If true, show the estimated contours instead of interpolating */
    int ShowContours;
        
    /** @brief Interpolation scale factor */
    double ScaleFactor;
    /** @brief If true, sample on the centered grid */
    int CenteredGrid;
    /** @brief Point-spread function (PSF) */
    double (*Psf)(double, double, const void *Param);
    /** @brief PSF parameter (e.g., Gaussian standard deviation) */
    double PsfParam;
    /** @brief PSF support radius */
    double PsfRadius;
    
    /* Method parameters */
    /** @brief \f$\sigma_\tau\f$, tagential spread of \f$\rho\f$ */
    double RhoSigmaNormal;
    /** @brief \f$\sigma_\nu\f$, normal spread of \f$\rho\f$ */
    double RhoSigmaTangent;
    /** @brief Method number of refinement passes */
    int RefinementPasses;
} programparams;


static int ParseParams(programparams *Param, int argc, char *argv[]);


/** @brief Print helpful program usage instructions */
static void PrintHelpMessage()
{
    printf("Contour stencil interpolation demo, P. Getreuer 2011\n\n");
    printf("Usage: iminterps [options] <input file> <output file>\n\n"
           "Only " READIMAGE_FORMATS_SUPPORTED " images are supported.\n\n");
    printf("Options:\n");
    printf("   -x <number>  the scale factor\n");
    printf("   -p <number>  sigma_h, the blur size of the point spread function\n");
    printf("   -g <grid>    grid to use for resampling, choices for <grid> are\n"
           "                centered     grid with centered alignment (default)\n"
           "                topleft      the top-left anchored grid\n\n");
    printf("   -s           show the estimated contours instead of interpolating\n\n");
    printf("   -t <number>  sigma_tau, spread of rho in the tagential direction\n");
    printf("   -n <number>  sigma_nu, spread of rho in the normal direction\n");
    printf("   -r <number>  the number of refinement passes\n\n");
#ifdef USE_LIBJPEG
    printf("   -q <number>  quality for saving JPEG images (0 to 100)\n\n");
#endif
    printf("When the -s option is given, the output file may be EPS, SVG, or a\n");
    printf("supported raster format.  Only for SVG, the background image must be\n");
    printf("written as a separate (raster) image, with the syntax\n");
    printf("   sinterp [options] -s <input file> <output.svg> <bg image file>\n\n");
    printf("Example: 4x scaling, sigma_h = 0.5, 2 refinement passes\n"
        "   iminterps -x 4 -p 0.5 -r 2 frog.bmp frog-4x.bmp\n");
}


/* The point spread function (PSF) */
static double GaussianPsf(double x, double y, const void *Param)
{
    const double Sigma = *((double *)Param);
    const double SigmaSqr = Sigma*Sigma;
    return exp(-(x*x + y*y)/(2.0*SigmaSqr))/(M_2PI*SigmaSqr);
}


/** @brief Distance function for a line-shaped stencil */
double LinearPatch(ATTRIBUTE_UNUSED double x, double y,
    ATTRIBUTE_UNUSED const void *Param)
{
    return y;
}

/** @brief Drawing commands to draw a line-shaped stencil */
int DrawLineStencil(pen *Pen, ATTRIBUTE_UNUSED const void *Param)
{
    return PenDrawLine(Pen, -1, 0, 1, 0);
}


/** @brief Distance function for a corner stencil */
double Corner(double x, double y, ATTRIBUTE_UNUSED const void *Param)
{
    return (x < y) ? x : y;
}

/** @brief Drawing commands to draw a corner stencil */
int DrawCornerStencil(pen *Pen, ATTRIBUTE_UNUSED const void *Param)
{
    return PenDrawLine(Pen, 1, 0, 0, 0) && PenDrawLine(Pen, 0, 0, 0, 1);
}


/** @brief Distance function for a T-junction stencil */
double TJunction(double x, double y, ATTRIBUTE_UNUSED const void *Param)
{
    y += 1;
    return (fabs(x) < y) ? fabs(x) : y;
}


/** @brief Distance function for a T-junction stencil */
double CurveDist(double x, double y, const void *Param)
{
    const double a = *((const double *)Param);
    const double a2 = a*a;
    const double xa = fabs(x)/a2;
    const double z = (1 - a*y)*(2/(3*a2));
    double x0, Temp;
    
    if((Temp = xa*xa + z*z*z) >= 0)
    {
        Temp = sqrt(Temp);
        x0 = pow(xa + Temp, 1.0/3.0)
            + ((a*y < 1) ? -1:1)*pow(fabs(-xa + Temp), 1.0/3.0);
    }
    else
    {
        Temp = sqrt(-z);
        x0 = 2*Temp*cos(acos(xa/(-z*Temp))/3);
    }
    
    Temp = x0*x0;
    return sqrt(a2*Temp + 1)*((a/2)*Temp - y);
}

/** @brief Drawing commands to draw a parabola stencil */
int DrawCurveStencil(pen *Pen, const void *Param)
{
    const double a = *((const double *)Param);
    const double x = 0.985;
    
    return PenDrawQBezier(Pen, -x, 0.5*a*x*x, 0, -0.5*a*x*x, x, 0.5*a*x*x);
}


/**
 * @brief Perform contour stencil image interpolation
 * @param Param program parameters
 * @param StencilSet the stencil set \f$\Sigma\f$
 * @return 1 on success, 0 on failure
 */
int Interpolate(programparams Param, const sset *StencilSet)
{
    const int FilterScaleFactor = 2;
    sinterp *SInterp;
    float *Input = NULL, *Filtered = NULL, *Output = NULL;
    float *FilterRhoSamples = NULL, *RhoSamples = NULL;
    int *BestStencil = NULL;
    unsigned long StartTime0, StartTime;
    int IntegerScaleFactor = (int)(Param.ScaleFactor + 0.5);
    int InputWidth, InputHeight, OutputWidth, OutputHeight, OutputNumEl;
    int i, Success = 1;

    
    printf("Interp stencils... \t");
    StartTime = Clock();
    
    if(!(SInterp = NewSInterp(StencilSet,
        (float)Param.RhoSigmaTangent, (float)Param.RhoSigmaNormal)))
        goto Catch;
    
    printf("%7.3f s\n", (Clock() - StartTime)*0.001f);
    
    if(Param.RefinementPasses || IntegerScaleFactor == Param.ScaleFactor)
    {
        printf("Precomputing rho... \t");
        StartTime = Clock();
        
        if(Param.RefinementPasses)
            if(!(FilterRhoSamples = ComputeRhoSamples(SInterp,
                FilterScaleFactor, Param.CenteredGrid)))
                goto Catch;
        
        if(IntegerScaleFactor == Param.ScaleFactor)
            if(!(RhoSamples = ComputeRhoSamples(SInterp,
                (int)ceil(Param.ScaleFactor), Param.CenteredGrid)))
                goto Catch;
        
        printf("%7.3f s\n", (Clock() - StartTime)*0.001f);
    }
    
    if(!(Input = (float *)ReadImage(&InputWidth, &InputHeight,
        Param.InputFile, IMAGEIO_FLOAT | IMAGEIO_RGB)))
        goto Catch;
    
    OutputWidth = (int)ceil(Param.ScaleFactor*InputWidth);
    OutputHeight = (int)ceil(Param.ScaleFactor*InputHeight);
    OutputNumEl = 3*OutputWidth*OutputHeight;
    printf("\n%gx interpolation, %dx%d -> %dx%d\n", Param.ScaleFactor,
        InputWidth, InputHeight, OutputWidth, OutputHeight);

    if(Param.CenteredGrid)
        printf("centered grid, ");
    else
        printf("top-left grid, ");

    if(Param.RefinementPasses == 1)
        printf("1 refinement pass");
    else if(Param.RefinementPasses > 1)
        printf("%d refinement passes", Param.RefinementPasses);
        
    printf("\nsigma_h = %g, sigma_tau = %g, sigma_nu = %g\n",
        Param.PsfParam, Param.RhoSigmaTangent,
        Param.RhoSigmaNormal);
        
    if(!(Output = (float *)Malloc(sizeof(float)*3*OutputWidth*OutputHeight))
        || !(Filtered = (float *)Malloc(sizeof(float)*3*InputWidth*InputHeight))
        || !(BestStencil = (int *)Malloc(sizeof(int)*InputWidth*InputHeight)))
        goto Catch;
        
    StartTime0 = Clock();
        
    printf("Fitting stencils... \t");
    StartTime = Clock();
    FitStencils(BestStencil, StencilSet, Input, InputWidth, InputHeight);
    printf("%7.3f s\n", (Clock() - StartTime)*0.001f);
        
    if(Param.RefinementPasses)
    {
        printf("Prefiltering... \t");
        StartTime = Clock();
        
        if(!(Prefilter(Filtered, BestStencil, FilterRhoSamples,
            Input, InputWidth, InputHeight, FilterScaleFactor,
            Param.CenteredGrid, (float)Param.PsfParam,
            Param.RefinementPasses)))
            goto Catch;
        
        printf("%7.3f s\n", (Clock() - StartTime)*0.001f);
    }
    
    if(IntegerScaleFactor == Param.ScaleFactor)
    {
        printf("Interpolating int... \t");
        StartTime = Clock();
        
        for(i = 0; i < OutputNumEl; i++)
            Output[i] = 0;
        
        IntegerScalePass(Output, BestStencil, RhoSamples,
            (Param.RefinementPasses) ? Filtered : Input, InputWidth,
            InputHeight, IntegerScaleFactor, Param.CenteredGrid);
    }
    else
    {
        printf("Interpolating arb... \t");
        StartTime = Clock();

        ArbitraryScale(Output, OutputWidth, OutputHeight, BestStencil,
            SInterp, (Param.RefinementPasses) ? Filtered : Input,
            InputWidth, InputHeight, (float)Param.ScaleFactor,
            Param.CenteredGrid);
    }
    
    printf("%7.3f s\n", (Clock() - StartTime)*0.001f);
    
    printf("Total                   %7.3f s\n", (Clock() - StartTime0)*0.001f);
    
    /* Write the output image */
    if(!WriteImage(Output, OutputWidth, OutputHeight, Param.OutputFile,
        IMAGEIO_FLOAT | IMAGEIO_RGB, Param.JpegQuality))
        goto Catch;
    
    Success = 1;
Catch:
    Free(Filtered);
    Free(Output);
    Free(BestStencil);
    Free(Input);
    Free(RhoSamples);
    Free(FilterRhoSamples);
    return Success;
}


/** @brief Case-insensitive test to see if String ends with Suffix */
static int StringEndsWith(const char *String, const char *Suffix)
{
    unsigned i, StringLength = strlen(String), SuffixLength = strlen(Suffix);
    
    if(StringLength < SuffixLength)
        return 0;
    
    String += StringLength - SuffixLength;
    
    for(i = 0; i < SuffixLength; i++)
        if(tolower(String[i]) != tolower(Suffix[i]))
            return 0;
    
    return 1;
}


/**
 * @brief Show contours estimated by contour stencils
 * @param Param program parameters
 * @param StencilSet the stencil set \f$\Sigma\f$
 * @return 1 on success, 0 on failure
 */
int ShowContours(programparams Param, const sset *StencilSet)
{
    float *Input = NULL, *Output = NULL;
    int *BestStencil = NULL;
    unsigned long StartTime;
    int InputWidth, InputHeight, OutputWidth, OutputHeight;
    int Success = 0;

    
    if(!(Input = (float *)ReadImage(&InputWidth, &InputHeight,
        Param.InputFile, IMAGEIO_FLOAT | IMAGEIO_RGB)))
        goto Catch;
    
    OutputWidth = (int)ceil(Param.ScaleFactor*InputWidth);
    OutputHeight = (int)ceil(Param.ScaleFactor*InputHeight);
        
    if(!(BestStencil = (int *)Malloc(sizeof(int)*InputWidth*InputHeight)))
        goto Catch;
        
    printf("Fitting stencils... \t");
    StartTime = Clock();
    FitStencils(BestStencil, StencilSet, Input, InputWidth, InputHeight);
    printf("%7.3f s\n", (Clock() - StartTime)*0.001f);
    
    if(!(Output = (float *)Malloc(sizeof(float)*3*OutputWidth*OutputHeight)))
        goto Catch;
    
    printf("Drawing contours... \t");
    StartTime = Clock();
    
    if(StringEndsWith(Param.OutputFile, ".eps"))
    {   /* Draw the image as an EPS vector graphic.  The background image data
           is embedded in the EPS. */
        if(!DrawContoursEps(Param.OutputFile,
        BestStencil, StencilSet, Input, InputWidth, InputHeight))
            goto Catch;
    }
    else if(StringEndsWith(Param.OutputFile, ".svg"))
    {   /* Draw the image as an SVG vector graphic.  The background image is
           drawn as a raster graphic to Param.BgFile. */
        if(!DrawContoursSvg(Param.OutputFile, Param.BgFile,
            BestStencil, StencilSet,
            Input, InputWidth, InputHeight))
            goto Catch;
    }
    else
    {   /* Draw the image as a raster graphic. */
        if(!DrawContours(Param.OutputFile, OutputWidth, OutputHeight,
            Param.ScaleFactor, BestStencil, StencilSet,
            Input, InputWidth, InputHeight))
            goto Catch;
    }
        
    printf("%7.3f s\n", (Clock() - StartTime)*0.001f);
    Success = 1;
Catch:
    Free(Output);
    Free(BestStencil);
    Free(Input);
    return Success;
}


int main(int argc, char *argv[])
{
    programparams Param;
    static const double a1 = M_1_SQRT2;
    static const double a2 = 1;
    sset *StencilSet;
    unsigned long StartTime;
    int S, Status = 1;
    
    
    /* Parse command line parameters */
    if(!ParseParams(&Param, argc, argv))
        return 1;
    
    printf("Building stencils... \t");
    StartTime = Clock();
    
    if(!(StencilSet = NewStencilSet(2, Param.Psf,
        (void *)&Param.PsfParam, Param.PsfRadius)))
        goto Catch;
    
    for(S = 0; S < 32; S++)
        AddStencil(StencilSet, LinearPatch, NULL, S*M_PI_8/4,
            DrawLineStencil, NULL, ColorLine);

    for(S = 0; S < 8; S++)
        AddStencil(StencilSet, CurveDist, &a1, S*M_PI_4,
            DrawCurveStencil, (const void *)&a1, ColorCurve1);
    
    for(S = 0; S < 8; S++)
        AddStencil(StencilSet, CurveDist, &a2, S*M_PI_4,
            DrawCurveStencil, (const void *)&a2, ColorCurve2);
        
    for(S = 0; S < 8; S++)
        AddStencil(StencilSet, Corner, NULL, S*M_PI_4,
            DrawCornerStencil, NULL, ColorCorner);

    printf("%7.3f s\n", (Clock() - StartTime)*0.001f);
    
    /* PrintStencilConfusion(StencilSet); */
    
    if(Param.ShowContours)
    {
        if(!ShowContours(Param, StencilSet))
            goto Catch;
    }
    else
    {
        if(!Interpolate(Param, StencilSet))
            goto Catch;
    }
    
    Status = 0;
Catch:
    FreeStencilSet(StencilSet);
    return Status;
}


/**
 * @brief Parse command line parameters
 * @param Param struct to hold the program parameters
 * @param argc, argv argument list received by main()
 * @return 1 on success, 0 on failure
 */
static int ParseParams(programparams *Param, int argc, char *argv[])
{
    static const char *DefaultOutputFile = (const char *)"out.bmp";
    char *OptionString;
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
    Param->BgFile = NULL;
    Param->JpegQuality = 80;
    
    Param->ShowContours = 0;
    
    Param->ScaleFactor = 4;
    Param->CenteredGrid = 1;
    Param->Psf = GaussianPsf;
    Param->PsfParam = DEFAULT_PSF_SIGMA;
            
    Param->RhoSigmaTangent = DEFAULT_RHO_SIGMA_TANGENT;
    Param->RhoSigmaNormal = DEFAULT_RHO_SIGMA_NORMAL;
    Param->RefinementPasses = DEFAULT_REFINEMENT_PASSES;

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
            case 'x':
                Param->ScaleFactor = atof(OptionString);

                if(Param->ScaleFactor < 1)
                {
                    ErrorMessage("Scale factor cannot be less than 1.0.\n");
                    return 0;
                }
                break;
            case 'g':
                if(!strcmp(OptionString, "centered")
                    || !strcmp(OptionString, "center"))
                    Param->CenteredGrid = 1;
                else if(!strcmp(OptionString, "topleft")
                    || !strcmp(OptionString, "top-left"))
                    Param->CenteredGrid = 0;
                else
                {
                    ErrorMessage("Grid must be either \"centered\" or \"topleft\".\n");
                    return 0;
                }
                break;
            case 'r':
                Param->RefinementPasses = atoi(OptionString);

                if(Param->RefinementPasses < 0 || Param->RefinementPasses > 10000)
                {
                    ErrorMessage("Invalid number of refinement steps.\n");
                    return 0;
                }
                break;
            case 'p':
                Param->PsfParam = atof(OptionString);

                if(Param->PsfParam < 0.0)
                {
                    ErrorMessage("Point spread blur size must be nonnegative.\n");
                    return 0;
                }
                else if(Param->PsfParam > 2.0)
                {
                    ErrorMessage("Point spread blur size is too large.\n");
                    return 0;
                }
                break;
            case 't':
                Param->RhoSigmaTangent = atof(OptionString);

                if(Param->RhoSigmaTangent <= 0.0)
                {
                    ErrorMessage("sigma_tau must be positive.\n");
                    return 0;
                }
                break;
            case 'n':
                Param->RhoSigmaNormal = atof(OptionString);

                if(Param->RhoSigmaNormal <= 0.0)
                {
                    ErrorMessage("sigma_nu must be positive.\n");
                    return 0;
                }
                break;
            case 's':
                Param->ShowContours = 1;
                i--;
                break;
                      
#ifdef USE_LIBJPEG
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
            else if(!Param->OutputFile)
                Param->OutputFile = argv[i];
            else
                Param->BgFile = argv[i];

            i++;
        }
    }
    
    if(!Param->InputFile)
    {
        PrintHelpMessage();
        return 0;
    }
    
    if(!Param->OutputFile)
        Param->OutputFile = DefaultOutputFile;
        
    if(Param->PsfParam == 0)
    {
        Param->Psf = NULL;
        Param->RefinementPasses = 0;
    }
    
    Param->PsfRadius = 4*Param->PsfParam;
    return 1;
}
