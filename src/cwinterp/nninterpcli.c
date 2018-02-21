/**
 * @file nninterpcli.c
 * @brief Nearest neighbor image interpolation utility
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

#include <math.h>
#include <string.h>
#include <ctype.h>
#include "nninterp.h"

#define VERBOSE 0


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
    /** @brief Quality for saving JPEG images (0 to 100) */
    int JpegQuality;
    /** @brief If true, sample on the centered grid */
    int CenteredGrid;
    /** @brief Scale factor */
    float ScaleFactor;
} programparams;


int ParseParams(programparams *Param, int argc, char *argv[]);


/** @brief Print program usage help message */
void PrintHelpMessage()
{
    puts("Nearest neighbor interpolation utility, P. Getreuer 2011\n");
    puts("Usage: iminterpnn [options] <input file> <output file>\n"
           "Only " READIMAGE_FORMATS_SUPPORTED " images are supported.\n");
    puts("Options:\n");
    puts("   -x <number>  the interpolation factor\n");
    puts(
    "   -g <grid>    grid to use for resampling, choices for <grid> are\n"
    "                centered     grid with centered alignment (default)\n"
    "                topleft      the top-left anchored grid\n");
#ifdef USE_LIBJPEG
    puts("   -q <number>  quality for saving JPEG images (0 to 100)\n");
#endif
    puts("Example: interpolate by factor 2.5\n"
         "   iminterpnn -x 2.5 frog.bmp coarse.bmp");
}


int main(int argc, char *argv[])
{
    programparams Param;
    image u = {NULL, 0, 0}, v = {NULL, 0, 0};
    int Status = 1;
    
    
    if(!ParseParams(&Param, argc, argv))
        return 0;
            
    /* Read the input image */
    if(!(v.Data = (uint32_t *)ReadImage(&v.Width, &v.Height, Param.InputFile,
        IMAGEIO_U8 | IMAGEIO_RGBA)))
        goto Catch;
    
    /* Allocate the output image */
    u.Width = (int)ceil(v.Width * Param.ScaleFactor);
    u.Height = (int)ceil(v.Height * Param.ScaleFactor);
#if VERBOSE > 0
    printf("%dx%d input --> %dx%d output\n", v.Width, v.Height, u.Width, u.Height);
#endif
    
    if(!(u.Data = (uint32_t *)Malloc(sizeof(uint32_t)*
        ((long int)u.Width)*((long int)u.Height))))
    {
        fprintf(stderr, "Memory allocation failed.\n");
        goto Catch;
    }
        
    NearestInterp(u.Data, u.Width, u.Height,
                  v.Data, v.Width, v.Height,
                  Param.ScaleFactor, Param.CenteredGrid);
    
    /* Write the output image */
    if(!WriteImage(u.Data, u.Width, u.Height, Param.OutputFile,
        IMAGEIO_U8 | IMAGEIO_RGBA, Param.JpegQuality))
        goto Catch;
#if VERBOSE > 0
    else
        printf("Output written to \"%s\".\n", Param.OutputFile);
#endif
    
    Status = 0; /* Finished successfully, set exit status to zero. */
Catch:
    Free(u.Data);
    Free(v.Data);
    return Status;
}


int ParseParams(programparams *Param, int argc, char *argv[])
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
    Param->JpegQuality = 99;
    Param->ScaleFactor = 1;
    Param->CenteredGrid = 1;
    
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
                Param->ScaleFactor = (float)atof(OptionString);

                if(Param->ScaleFactor < 0)
                {
                    ErrorMessage("Invalid scale factor.\n");
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
