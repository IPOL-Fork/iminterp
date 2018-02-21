/**
 * @file imcoarsen.c
 * @brief Image coarsening utility program
 * @author Pascal Getreuer <getreuer@gmail.com>
 *
 * This file implements the imcoarsen program, a command line tool for
 * coarsening an image by Gaussain smoothing followed by downsampling.
 * See the routine PrintHelpMessage for usage details.
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

#include <ipol/imageio.h>

#define VERBOSE 0

/** @brief Approximate Gaussian with this number of standard deviations */
#define NUMSTDS 4


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

typedef enum
{
    BOUNDARY_CONSTANT = 0,
    BOUNDARY_HSYMMETRIC = 1,
    BOUNDARY_WSYMMETRIC = 2,
    BOUNDARY_PERIODIC = 3
} boundaryhandling;

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
    /** @brief Type of boundary handling */
    boundaryhandling Boundary;
    /** @brief Scale factor */
    float ScaleFactor;
    /** @brief Gaussian point spread function standard deviation */
    float PsfSigma;
} programparams;


int ParseParams(programparams *Param, int argc, char *argv[]);
int Coarsen(image v, image u, programparams Param);

/** @brief Print program usage help message */
void PrintHelpMessage()
{
    puts("Image coarsening utility, P. Getreuer 2010-2011\n");
    puts("Usage: imintcoarsen [options] <input file> <output file>\n\n"
        "Only " READIMAGE_FORMATS_SUPPORTED " images are supported.\n");
    puts("Options:");
    puts("  -x <number>  the coarsening factor (>= 1.0, may be non-integer)");
    puts("  -p <number>  sigma_h, blur size of the Gaussian point spread");
    puts("               function in units of output pixels.");
    puts("  -b <ext>     boundary handling extension, choices for <ext> are");
    puts("               const       constant extension");
    puts("               hsym        half-sample symmetric (default)");
    puts("               wsym        whole-sample symmetric");
    puts("               per         periodic");
    puts("  -g <grid>    grid to use for resampling, choices for <grid> are");
    puts("               centered    grid with centered alignment (default)");
    puts("               topleft     the top-left anchored grid\n");
#ifdef USE_LIBJPEG
    puts("  -q <number>  quality for saving JPEG images (0 to 100)\n");
#endif
    puts("Example: coarsen by factor 2.5");
    puts("   imintcoarsen -x 2.5 -p 0.35 frog.bmp coarse.bmp");
}


int main(int argc, char *argv[])
{
    programparams Param;
    image u = {NULL, 0, 0}, v = {NULL, 0, 0};
    int Status = 1;
    
    
    if(!ParseParams(&Param, argc, argv))
        return 0;
            
    /* Read the input image */
    if(!(u.Data = (uint32_t *)ReadImage(&u.Width, &u.Height, Param.InputFile,
        IMAGEIO_U8 | IMAGEIO_RGBA)))
        goto Catch;

    if(Param.ScaleFactor >= u.Width || Param.ScaleFactor >= u.Height)
    {
        ErrorMessage("Image is too small for scale factor.\n");
        goto Catch;
    }
    
    /* Allocate the output image */
    v.Width = (int)ceil(u.Width / Param.ScaleFactor);
    v.Height = (int)ceil(u.Height / Param.ScaleFactor);
#if VERBOSE > 0
    printf("%dx%d input -> %dx%d output\n", u.Width, u.Height, v.Width, v.Height);
#endif
    
    if(!(v.Data = (uint32_t *)Malloc(sizeof(uint32_t)*
        ((long int)v.Width)*((long int)v.Height))))
        goto Catch;
       
    /* Convolution followed by downsampling */
    if(!Coarsen(v, u, Param))
        goto Catch;
    
    /* Write the output image */
    if(!WriteImage(v.Data, v.Width, v.Height, Param.OutputFile,
        IMAGEIO_U8 | IMAGEIO_RGBA, Param.JpegQuality))
        goto Catch;
#if VERBOSE > 0
    else
        printf("Output written to \"%s\".\n", Param.OutputFile);
#endif
    
    Status = 0;	/* Finished successfully, set exit status to zero. */
    
Catch:
    Free(v.Data);
    Free(u.Data);
    return Status;
}


float Sqr(float x)
{
    return x*x;
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


/**
 * @brief Boundary handling function for periodic extension
 * @param N is the data length
 * @param i is an index into the data
 * @return an index that is always between 0 and N - 1
 */
static int PerExtension(int N, int i)
{
    while(1)
    {
        if(i < 0)
            i += N;
        else if(i >= N)
            i -= N;
        else
            return i;
    }
}


int (*ExtensionMethod[4])(int, int) =
    {ConstExtension, HSymExtension, WSymExtension, PerExtension};


int Coarsen(image v, image u, programparams Param)
{
    int (*Extension)(int, int) = ExtensionMethod[Param.Boundary];
    const float PsfRadius = NUMSTDS*Param.PsfSigma*Param.ScaleFactor;
    const int PsfWidth = 1 + ((PsfRadius == 0) ? 1 : (int)ceil(2*PsfRadius));
    float *Temp = NULL, *PsfBuf = NULL;
    float ExpDenom, Weight, Sum[4], DenomSum;
    float XStart, YStart, X, Y;
    uint32_t Pixel;
    int IndexX0, IndexY0, SrcOffset, DestOffset;
    int x, y, n, c, Success = 0;
    
    
    if(!(Temp = (float *)Malloc(sizeof(float)*4*v.Width*u.Height))
        || !(PsfBuf = (float *)Malloc(sizeof(float)*PsfWidth)))
        goto Catch;
    
    ExpDenom = 2 * Sqr(Param.PsfSigma*Param.ScaleFactor);
    
    if(Param.CenteredGrid)
    {
        XStart = (1/Param.ScaleFactor - 1)/2;
        YStart = (1/Param.ScaleFactor - 1)/2;
    }
    else
        XStart = YStart = 0;
    
    for(x = 0; x < v.Width; x++)
    {
        X = (-XStart + x)*Param.ScaleFactor;
        IndexX0 = (int)ceil(X - PsfRadius - 0.5f);
        DenomSum = 0;
        
        /* Evaluate the PSF */
        for(n = 0; n < PsfWidth; n++)
        {
            PsfBuf[n] = Sqr(X - (IndexX0 + n));
            
            if(!n || PsfBuf[n] < DenomSum)
                DenomSum = PsfBuf[n];
        }
        
        if(ExpDenom > 0)
            for(n = 0; n < PsfWidth; n++)
                PsfBuf[n] = (float)exp((DenomSum - PsfBuf[n]) / ExpDenom);
        else if(IndexX0 == (int)floor(X - PsfRadius + 0.5f))
        {   /* If PsfSigma = 0, sample the nearest neighbor */
            PsfBuf[0] = 1;
            PsfBuf[1] = 0;
        }
        else /* At a half integer, average the two nearest neighbors */
            PsfBuf[0] = PsfBuf[1] = 0.5f;
        
        DenomSum = 0;
        
        for(n = 0; n < PsfWidth; n++)
            DenomSum += PsfBuf[n];
        
        for(y = 0, SrcOffset = 0, DestOffset = x; y < u.Height;
            y++, SrcOffset += u.Width, DestOffset += v.Width)
        {
            Sum[0] = Sum[1] = Sum[2] = Sum[3] = 0;
            
            for(n = 0; n < PsfWidth; n++)
            {
                Weight = PsfBuf[n];
                Pixel = u.Data[Extension(u.Width, IndexX0 + n) + SrcOffset];
                
                for(c = 0; c < 4; c++)
                    Sum[c] += (float)((uint8_t *)&Pixel)[c] * Weight;
            }
            
            for(c = 0; c < 4; c++)
                Temp[4*DestOffset + c] = Sum[c] / DenomSum;
        }
    }
    
    for(y = 0; y < v.Height; y++, v.Data += v.Width)
    {
        Y = (-YStart + y)*Param.ScaleFactor;
        IndexY0 = (int)ceil(Y - PsfRadius - 0.5f);
        DenomSum = 0;
        
        /* Evaluate the PSF */
        for(n = 0; n < PsfWidth; n++)
        {
            PsfBuf[n] = Sqr(Y - (IndexY0 + n));
            
            if(!n || PsfBuf[n] < DenomSum)
                DenomSum = PsfBuf[n];
        }
        
        if(ExpDenom > 0)
            for(n = 0; n < PsfWidth; n++)
                PsfBuf[n] = (float)exp((DenomSum - PsfBuf[n]) / ExpDenom);
        else if(IndexY0 == (int)floor(Y - PsfRadius + 0.5f))
        {   /* If PsfSigma = 0, sample the nearest neighbor */
            PsfBuf[0] = 1;
            PsfBuf[1] = 0;
        }
        else /* At a half integer, average the two nearest neighbors */
            PsfBuf[0] = PsfBuf[1] = 0.5f;
        
        DenomSum = 0;
        
        for(n = 0; n < PsfWidth; n++)
            DenomSum += PsfBuf[n];
        
        for(x = 0; x < v.Width; x++)
        {
            Sum[0] = Sum[1] = Sum[2] = Sum[3] = 0;
            
            for(n = 0; n < PsfWidth; n++)
            {
                SrcOffset = x + v.Width*Extension(u.Height, IndexY0 + n);
                Weight = PsfBuf[n];
                
                for(c = 0; c < 4; c++)
                    Sum[c] += Temp[4*SrcOffset + c] * Weight;
            }
            
            for(c = 0; c < 4; c++)
                ((uint8_t *)&Pixel)[c] = (int)(Sum[c] / DenomSum + 0.5f);
            
            v.Data[x] = Pixel;
        }
    }
    
    Success = 1;
Catch:
    Free(PsfBuf);
    Free(Temp);
    return Success;
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
    Param->PsfSigma = 0.35f;
    Param->CenteredGrid = 1;
    Param->Boundary = BOUNDARY_HSYMMETRIC;
    
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

                if(Param->ScaleFactor < 1)
                {
                    ErrorMessage("Invalid scale factor.\n");
                    return 0;
                }
                break;
            case 'p':
                Param->PsfSigma = (float)atof(OptionString);

                if(Param->PsfSigma < 0.0)
                {
                    ErrorMessage("Point spread blur size must be nonnegative.\n");
                    return 0;
                }
                break;
            case 'b':
                if(!strcmp(OptionString, "const"))
                    Param->Boundary = BOUNDARY_CONSTANT;
                else if(!strcmp(OptionString, "hsym"))
                    Param->Boundary = BOUNDARY_HSYMMETRIC;
                else if(!strcmp(OptionString, "wsym"))
                    Param->Boundary = BOUNDARY_WSYMMETRIC;
                else if(!strcmp(OptionString, "per"))
                    Param->Boundary = BOUNDARY_PERIODIC;
                else
                {
                    ErrorMessage("Boundary extension must be either \"const\", \"hsym\", or \"wsym\".\n");
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
