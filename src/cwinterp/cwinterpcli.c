/**
 * @file cwinterpcli.c
 * @brief Image Interpolation with Contour Stencils command line program
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

#include <math.h>
#include <string.h>
#include <ctype.h>

#include <ipol/imageio.h>
#include "cwinterp.h"

/** @brief Set to 1 for informative program output, 0 for quiet */
#define VERBOSE     0

/** @brief struct representing an image */
typedef struct
{
    /** @brief 32-bit RGBA image data */
    uint32_t *Data;
    /** @brief Image width */
    int Width;
    /** @brief Image height */
    int Height;
} image;

/** @brief struct of program parameters */
typedef struct
{
    /** @brief Input file name */
    char *InputFile;
    /** @brief Output file name */
    char *OutputFile;
    /** @brief If true, only show the estimated contour orientations */
    int OnlyShowContours;
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;
    /** @brief interpolation parameters */
    cwparams Cw;
    
    int TestFlag;
} programparams;


static int ParseParams(programparams *Param, int argc, char *argv[]);


static void PrintHelpMessage()
{
    puts("Image Interpolation with Contour Stencils, P. Getreuer 2010-2011\n");
    puts("Usage: iminterpcw [options] <input file> <output file>\n\n"
        "Only " READIMAGE_FORMATS_SUPPORTED " images are supported.\n");
    puts("Options:");
    puts("  -x <number>  the scale factor (may be non-integer)");
    puts("  -p <number>  sigma_h, the blur size of the point spread function");
    puts("  -g <grid>    grid to use for resampling, choices for <grid> are\n"
    "               centered     grid with centered alignment (default)\n"
    "               topleft      the top-left anchored grid\n");
    puts("  -s           "
    "show the estimated orientations instead of interpolating\n");
    puts("  -t <number>  sigma_tau, spread of phi in the tagential direction");
    puts("  -n <number>  sigma_nu, spread of phi in the normal direction");
    puts("  -r <number>  the number of refinement passes\n");
#ifdef USE_LIBJPEG
    puts("  -q <number>  quality for saving JPEG images (0 to 100)\n");
#endif
    puts("Example: 4x scaling, sigma_h = 0.35, 2 refinement passes\n"
    "  iminterpcw -x 4 -p 0.35 -r 2 frog.bmp frog-4x.bmp");
}

int main(int argc, char *argv[])
{
    programparams Param;
    image v = {NULL, 0, 0}, u = {NULL, 0, 0};
    int32_t *Psi = NULL;
    int Status = 1;
    
    
    /* Parse command line parameters */
    if(!ParseParams(&Param, argc, argv))
        return 0;

    /* Perform precomputations (but not when only showing contours) */
    if(!Param.OnlyShowContours)
        if(!(Psi = PreCWInterp(Param.Cw)))
            goto Catch;
    
    /* Read the input image */
    if(!(v.Data = (uint32_t *)ReadImage(&v.Width, &v.Height, Param.InputFile,
        IMAGEIO_U8 | IMAGEIO_RGBA)))
        goto Catch;
    
    /* Allocate the output image */
    u.Width = (int)ceil(Param.Cw.ScaleFactor * v.Width);
    u.Height = (int)ceil(Param.Cw.ScaleFactor * v.Height);
    
    if(!(u.Data = (uint32_t *)Malloc(sizeof(uint32_t)*
        ((long int)u.Width)*((long int)u.Height))))
        goto Catch;
    
    if(!Param.OnlyShowContours)
    {
        if(!Param.TestFlag && Param.Cw.ScaleFactor == ceil(Param.Cw.ScaleFactor))
        {
            printf("Integer scale factor %dx%d input -> %dx%d output\n",
                v.Width, v.Height, u.Width, u.Height);
            
            /* Perform interpolation by an integer scale factor. */
            if(!CWInterp(u.Data, v.Data, v.Width, v.Height, Psi, Param.Cw))
                goto Catch;
        }
        else
        {
            printf("Arbitrary scale factor %dx%d input -> %dx%d output\n",
                v.Width, v.Height, u.Width, u.Height);
            
            /* Perform interpolation by an arbitrary scale factor. */
            if(!CWInterpEx(u.Data, u.Width, u.Height,
                v.Data, v.Width, v.Height, Psi, Param.Cw))
                goto Catch;
        }
    }
    else
    {
        /* Program has been called with the -s flag.  Show the contour
        orientation estimates using on the best-fitting contour stencils. */
        if(!(DisplayContours(u.Data, u.Width, u.Height,
            v.Data, v.Width, v.Height, Param.Cw)))
            goto Catch;
    }
    
    /* Write the output image */
    if(!WriteImage(u.Data, u.Width, u.Height, Param.OutputFile,
        IMAGEIO_U8 | IMAGEIO_RGBA, Param.JpegQuality))
        goto Catch;
#if VERBOSE > 0
    else
        printf("Output written to \"%s\".\n", Param.OutputFile);
#endif
    
    Status = 0;	/* Finished successfully, set exit status to zero. */
    
Catch:  /* This label is used for error handling.  If something went wrong
        above (which may be out of memory, file not found, or a computation
        error), then execution jumps to this point to clean up and exit. */
    Free(u.Data);
    Free(v.Data);
    Free(Psi);
    return Status;
}


static int ParseParams(programparams *Param, int argc, char *argv[])
{
    static char *DefaultOutputFile = (char *)"out.bmp";
    char *OptionString;
    char OptionChar;
    int i;

    
    if(argc < 2)
    {
        PrintHelpMessage();
        return 0;
    }

    /* Set parameter defaults */
    Param->InputFile = 0;
    Param->OutputFile = DefaultOutputFile;
    Param->OnlyShowContours = 0;
    Param->JpegQuality = 70;
    
    Param->Cw.ScaleFactor = 4;
    Param->Cw.CenteredGrid = 1;
    Param->Cw.RefinementSteps = 2;
    Param->Cw.PsfSigma = 0.35;
    Param->Cw.PhiSigmaTangent = 1.2;
    Param->Cw.PhiSigmaNormal = 0.6;
    
    Param->TestFlag = 0;

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
                Param->Cw.ScaleFactor = atof(OptionString);

                if(Param->Cw.ScaleFactor < 1)
                {
                    ErrorMessage("Scale factor cannot be less than 1.0.\n");
                    return 0;
                }
                break;
            case 'g':
                if(!strcmp(OptionString, "centered")
                    || !strcmp(OptionString, "center"))
                    Param->Cw.CenteredGrid = 1;
                else if(!strcmp(OptionString, "topleft")
                    || !strcmp(OptionString, "top-left"))
                    Param->Cw.CenteredGrid = 0;
                else
                {
                    ErrorMessage("Grid must be either \"centered\" or \"topleft\".\n");
                    return 0;
                }
                break;
            case 'r':
                Param->Cw.RefinementSteps = atoi(OptionString);

                if(Param->Cw.RefinementSteps < 0)
                {
                    ErrorMessage("Invalid number of refinement passes.\n");
                    return 0;
                }
                break;
            case 'p':
                Param->Cw.PsfSigma = atof(OptionString);

                if(Param->Cw.PsfSigma < 0.0)
                {
                    ErrorMessage("Point spread blur size must be nonnegative.\n");
                    return 0;
                }
                else if(Param->Cw.PsfSigma > 2.0)
                {
                    ErrorMessage("Point spread blur size is too large.\n");
                    return 0;
                }
                break;
            case 't':
                Param->Cw.PhiSigmaTangent = atof(OptionString);

                if(Param->Cw.PhiSigmaTangent <= 0.0)
                {
                    ErrorMessage("sigma_tau must be positive.\n");
                    return 0;
                }
                break;
            case 'n':
                Param->Cw.PhiSigmaNormal = atof(OptionString);

                if(Param->Cw.PhiSigmaNormal <= 0.0)
                {
                    ErrorMessage("sigma_nu must be positive.\n");
                    return 0;
                }
                break;
            case 's':
                Param->OnlyShowContours = 1;
                i--;
                break;
                
            case 'T':
                Param->TestFlag = atoi(OptionString);
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
    
    if(Param->Cw.PsfSigma == 0)
        Param->Cw.RefinementSteps = 0;
    
    if(!Param->OnlyShowContours)
    {
        /* Display the chosen parameters */
        printf("Factor-%g interpolation, ", Param->Cw.ScaleFactor);
        
        if(Param->Cw.CenteredGrid)
            printf("centered grid, ");
        else
            printf("top-left grid, ");
                
        if(Param->Cw.RefinementSteps == 1)
            printf("1 refinement pass,");
        else if(Param->Cw.RefinementSteps > 1)
            printf("%d refinement passes,", Param->Cw.RefinementSteps);
        
        printf("\nsigma_h = %g, sigma_tau = %g, sigma_nu = %g\n",
            Param->Cw.PsfSigma, Param->Cw.PhiSigmaTangent,
            Param->Cw.PhiSigmaNormal);
    }

    return 1;
}
