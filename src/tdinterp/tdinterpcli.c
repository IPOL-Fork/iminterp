/**
 * @file tdinterpcli.c
 * @brief Roussos-Maragos interpolation command line program
 * @author Pascal Getreuer <getreuer@gmail.com>
 *
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "imageio.h"
#include "tdinterp.h"


#define DEFAULT_PSF_SIGMA           0.35
#define DEFAULT_K                   (1.0/255)
#define DEFAULT_TOL                 (3e-4)
#define DEFAULT_MAXMETHODITER          50
#define DEFAULT_DIFFITER            5

#define VERBOSE     0


/** @brief struct representing an image */
typedef struct
{
    /** @brief 32-bit RGBA image data */
    float *Data;
    /** @brief Image width */
    int Width;
    /** @brief Image height */
    int Height;
} image;


/** @brief struct of program parameters */
typedef struct
{
    /** @brief Input file name */
    const char *InputFile;
    /** @brief Output file name */
    const char *OutputFile;
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;
    
    /** @brief Interpolation scale factor */
    int ScaleFactor;
    /** @brief PSF standard deviation */
    double PsfSigma;
    
    /** @brief Parameter K in constructing the tensor */
    double K;
    /** @brief Convergence tolerance */
    float Tol;
    /** @brief Maximum number of method iterations */
    int MaxMethodIter;
    /** @brief Number of diffusion iterations per method iteration */
    int DiffIter;
} programparams;


static int ParseParams(programparams *Param, int argc, char *argv[]);


static void PrintHelpMessage()
{
    puts("Tensor-driven diffusion interpolation demo, P. Getreuer 2010-2011\n");
    puts("Usage: iminterptd [options] <input file> <output file>\n\n"
        "Only " READIMAGE_FORMATS_SUPPORTED " images are supported.\n");
    puts("Options:");
    puts("  -x <number>  the scale factor (must be integer)");
    puts("  -p <number>  sigma_h, the blur size of the point spread function");
    puts("  -K <number>  K, parameter in constructing the tensor (default 1/255)");
    puts("  -t <number>  tol, convergence tolerance (default 3e-4)");
    puts("  -N <number>  N, maximum number of method iterations (default 50)");
    puts("  -n <number>  n, number of diffusion steps per method iteration (default 5)\n");
#ifdef USE_LIBJPEG
    puts("  -q <number>  quality for saving JPEG images (0 to 100)\n");
#endif
    puts("Example: 4x scaling, sigma_h = 0.5\n"
         "  iminterptd -x 4 -p 0.5 frog.bmp frog-4x.bmp");
}


int main(int argc, char *argv[])
{
    programparams Param;
    image v = {NULL, 0, 0}, u = {NULL, 0, 0};
    int Status = 1;
    
    
    if(!ParseParams(&Param, argc, argv))
        return 0;

    /* Read the input image */
    if(!(v.Data = (float *)ReadImage(&v.Width, &v.Height, Param.InputFile,
        IMAGEIO_FLOAT | IMAGEIO_RGB | IMAGEIO_PLANAR)))
        goto Catch;
    
    if(v.Width < 3 || v.Height < 3)
    {
        ErrorMessage("Image is too small (%dx%d).\n", v.Width, v.Height);
        goto Catch;
    }
    
    /* Allocate the output image */
    u.Width = Param.ScaleFactor * v.Width;
    u.Height = Param.ScaleFactor * v.Height;
    
    if(!(u.Data = (float *)Malloc(sizeof(float)*3*
        ((long int)u.Width)*((long int)u.Height))))
        goto Catch;
    
    /* Call the interpolation routine */
    if(!(RoussosInterp(u.Data, u.Width, u.Height,
        v.Data, v.Width, v.Height, Param.PsfSigma, Param.K,
        Param.Tol, Param.MaxMethodIter, Param.DiffIter)))
        goto Catch;
    
    /* Write the output image */
    if(!WriteImage(u.Data, u.Width, u.Height, Param.OutputFile,
        IMAGEIO_FLOAT | IMAGEIO_RGB | IMAGEIO_PLANAR, Param.JpegQuality))
        goto Catch;
    
    Status = 0; /* Finished successfully, set exit status to zero. */
    
Catch:
    Free(u.Data);
    Free(v.Data);
    return Status;
}


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
    Param->InputFile = 0;
    Param->OutputFile = DefaultOutputFile;
    Param->JpegQuality = 80;
    
    Param->ScaleFactor = 4;
    Param->PsfSigma = DEFAULT_PSF_SIGMA;
            
    Param->K = DEFAULT_K;
    Param->Tol = (float)DEFAULT_TOL;
    Param->MaxMethodIter = DEFAULT_MAXMETHODITER;
    Param->DiffIter = DEFAULT_DIFFITER;

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
                Param->ScaleFactor = atoi(OptionString);

                if(Param->ScaleFactor < 1)
                {
                    ErrorMessage("Scale factor cannot be less than 1.0.\n");
                    return 0;
                }
                break;
            case 'p':
                Param->PsfSigma = atof(OptionString);

                if(Param->PsfSigma < 0.0)
                {
                    ErrorMessage("Point spread blur size must be nonnegative.\n");
                    return 0;
                }
                else if(Param->PsfSigma > 1.0)
                {
                    ErrorMessage("Point spread blur size is too large.\n");
                    return 0;
                }
                break;
            case 'K':
                Param->K = atof(OptionString);

                if(Param->K <= 0.0)
                {
                    ErrorMessage("K must be positive.\n");
                    return 0;
                }
                break;
            case 't':
                Param->Tol = atof(OptionString);

                if(Param->Tol <= 0.0)
                {
                    ErrorMessage("Convergence tolerance must be positive.\n");
                    return 0;
                }
                break;
            case 'N':
                Param->MaxMethodIter = atoi(OptionString);

                if(Param->MaxMethodIter < 0)
                {
                    ErrorMessage("Number of method iterations must be nonnegative.\n");
                    return 0;
                }
                break;
            case 'n':
                Param->DiffIter = atoi(OptionString);

                if(Param->DiffIter <= 0)
                {
                    ErrorMessage("Number of diffusion iterations must be positive.\n");
                    return 0;
                }
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
    
    return 1;
}
