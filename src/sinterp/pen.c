/**
 * @file pen.c
 * @brief Vector drawing abstraction
 * @author Pascal Getreuer <getreuer@gmail.com>
 *
 * Copyright (c) 2010-2011, Pascal Getreuer
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD License. You
 * should have received a copy of this license along this program. If
 * not, see <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include <stdio.h>
#include "pen.h"
#include "svd2x2.h"


/** @brief Buffer size to use for file output */
#define FILE_BUFFER_CAPACITY    (1024*4)


/** @brief Identity transform */
const pentrans IdentityPenTrans = {1, 0, 0, 1, 0, 0};

/** @brief Powers of 10, must have up to 10^#PEN_MAXDIGITS */
static const unsigned long Pow10[] = {1,10,100,1000,10000,100000,
                                      1000000,10000000,100000000};


/** @brief Stub function for pen::DrawLine callback */
static int StubDrawLine(ATTRIBUTE_UNUSED pen *Pen,
    ATTRIBUTE_UNUSED float x1, ATTRIBUTE_UNUSED float y1,
    ATTRIBUTE_UNUSED float x2, ATTRIBUTE_UNUSED float y2)
{
    ErrorMessage("Pen::DrawLine callback has not been defined.\n");
    return 0;
}


/** @brief Stub function for pen::DrawQBezier callback */
static int StubDrawQBezier(pen *Pen, float x1, float y1,
    float x2, float y2, float x3, float y3)
{
    float L, t, xlast, ylast, xcur, ycur;
    int i, NumLines;
    
    /* Check whether DrawLine is available */
    if(!Pen || !(Pen->DrawLine) || Pen->DrawLine == StubDrawLine)
    {
        ErrorMessage("Pen::DrawQBezier callback has not been defined.\n");
        return 0;
    }
    
    /* Approximate the quadratic Bezier curve using lines */
    L = EstArcLenQBezier(x1, y1, x2, y2, x3, y3, Pen->Trans);
    NumLines = (int)(L*0.04);
    NumLines = (NumLines < 4) ? 4 : NumLines;
    xlast = x1;
    ylast = y1;
    
    for(i = 1; i <= NumLines; i++)
    {
        t = ((float)i)/NumLines;
        xcur = (1 - t)*((1 - t)*x1 + 2*t*x2) + t*t*x3;
        ycur = (1 - t)*((1 - t)*y1 + 2*t*y2) + t*t*y3;
        Pen->DrawLine(Pen, xlast, ylast, xcur, ycur);
        xlast = xcur;
        ylast = ycur;
    }
    
    return 1;
}


/** @brief Stub function for pen::DrawEllipse callback */
static int StubDrawEllipse(pen *Pen,
    float x, float y, float rx, float ry, float Theta)
{
    const float CosTheta = (float)cos(Theta), SinTheta = (float)sin(Theta);
    float C, t, xlast, ylast, xcur, ycur;
    int i, NumLines;
    
    /* Check whether DrawLine is available */
    if(!Pen || !(Pen->DrawLine) || Pen->DrawLine == StubDrawLine)
    {
        ErrorMessage("Pen::DrawEllipse callback has not been defined.\n");
        return 0;
    }
    
    /* Approximate the ellipse using lines */
    C = EstArcLenEllipse(x, y, rx, ry, Theta, Pen->Trans);
    NumLines = (int)(C*0.07);
    NumLines = (NumLines < 8) ? 8 : NumLines;
    
    xlast = x + CosTheta*rx;
    ylast = y - SinTheta*rx;
    
    for(i = 1; i <= NumLines; i++)
    {
        t = (i < NumLines) ? ((float)i)*M_2PI/NumLines : 0;
        xcur = x + CosTheta*rx*cos(t) + SinTheta*ry*sin(t);
        ycur = y - SinTheta*rx*cos(t) + CosTheta*ry*sin(t);
        Pen->DrawLine(Pen, xlast, ylast, xcur, ycur);
        xlast = xcur;
        ylast = ycur;
    }
    
    return 0;
}


/**
 * @brief Create a new pen
 * @return a pen*, or NULL on allocation failure
 *
 * Allocates a new pen object.  The caller must call FreePen() to free the pen
 * when done.
 */
pen *NewPen()
{
    static const pen InitPen = {NULL, NULL, 0, 0, PEN_DEFAULT_NUMDIGITS,
        {1, 0, 0, 1, 0, 0}, NULL, 0, 0, {0, 0, 0}, -1,
        NULL, StubDrawLine, StubDrawQBezier, StubDrawEllipse};
    pen *Pen;
        
    if((Pen = (pen *)Malloc(sizeof(pen))))
        *Pen = InitPen;
    
    return Pen;
}
    
    
/**
 * @brief Free all memory associated with a pen
 * @param Pen the pen object to free
 */
void FreePen(pen *Pen)
{
    if(Pen)
    {
        Free(Pen->Palette);
        Free(Pen);
    }
}


/**
 * @brief Get the stdio FILE* pointer associated with a pen
 * @param Pen the pen object
 */
FILE *PenGetFile(pen *Pen)
{
    return (Pen) ? Pen->File : NULL;
}


/**
 * @brief Draw a line
 * @param Pen the pen object
 * @param x1, y1, x2, y2 endpoints of the line
 * @return 1 on success, 0 on failure
 *
 * This function is a wrapper for pen::DrawLine.
 */
int PenDrawLine(pen *Pen, float x1, float y1, float x2, float y2)
{
    return (Pen && Pen->DrawLine(Pen, x1, y1, x2, y2));
}


/**
 * @brief Draw a rectangle
 * @param Pen the pen object
 * @param x1, y1, x2, y2 corners of the rectangle
 * @return 1 on success, 0 on failure
 */
int PenDrawRectangle(pen *Pen, float x1, float y1, float x2, float y2)
{
    return (Pen
        && Pen->DrawLine(Pen, x1, y1, x2, y1)
        && Pen->DrawLine(Pen, x2, y1, x2, y2)
        && Pen->DrawLine(Pen, x2, y2, x1, y2)
        && Pen->DrawLine(Pen, x1, y2, x1, y1));
}


/**
 * @brief Draw a quadratic Bezier curve
 * @param Pen the pen object
 * @param x1, y1, x2, y2, x3, y3 starting point, control point, and endpoint
 * @return 1 on success, 0 on failure
 *
 * This function is a wrapper for pen::DrawQBezier.
 */
int PenDrawQBezier(pen *Pen, float x1, float y1,
    float x2, float y2, float x3, float y3)
{
    return (Pen && Pen->DrawQBezier(Pen, x1, y1, x2, y2, x3, y3));
}


/**
 * @brief Draw a circle
 * @param Pen the pen object
 * @param x, y circle center
 * @param r circle radius
 * @return 1 on success, 0 on failure
 *
 * This function is a wrapper for pen::DrawEllipse.
 */
int PenDrawCircle(pen *Pen, float x, float y, float r)
{
    return (Pen && Pen->DrawEllipse(Pen, x, y, r, r, 0));
}


/**
 * @brief Draw an ellipse
 * @param Pen the pen object
 * @param x, y center of the ellipse
 * @param rx, ry radii of the ellipse
 * @param Theta ellipse rotation
 * @return 1 on success, 0 on failure
 *
 * This function is a wrapper for pen::DrawEllipse.
 */
int PenDrawEllipse(pen *Pen,
    float x, float y, float rx, float ry, float Theta)
{
    return (Pen && Pen->DrawEllipse(Pen, x, y, rx, ry, Theta));
}


/**
 * @brief Write a floating-point value without trailing zeros
 * @param File the file or stream
 * @param Value floating-point value to write
 * @param NumDigits
 *
 * Value is written in decimal notation using up to NumDigits digits after the
 * decimal point.  Trailing zeros are omitted.
 */
void PenWriteDouble(FILE *File, double Value, int NumDigits)
{
    int i;
    unsigned long Digit, Whole, Frac;

    /* Write sign and negate if Value is negative */
    if(Value < 0)
    {
        putc('-', File);
        Value = -Value;
    }
        
    if(NumDigits > 0)
    {
        if(NumDigits > PEN_MAXDIGITS)
            NumDigits = PEN_MAXDIGITS;
        
        /* Write the whole part */
        Whole = (unsigned long)Value;
        /* Represent fractional part as an integer */
        Frac = (unsigned long)(Pow10[NumDigits]*(Value - Whole) + 0.5f);
        
        /* Carry the rounding +0.5/Pow10[NumDigits] if necessary */
        if(Frac >= Pow10[NumDigits])
        {
            Whole++;
            Frac -= Pow10[NumDigits];
        }
        
        fprintf(File, "%lu", Whole);
        
        /* Write fractional part without trailing zeros */
        if(Frac > 0)
        {
            putc('.', File);
            
            for(i = NumDigits - 1; Frac > 0 && i >= 0; i--)
            {
                Digit = Frac / Pow10[i];
                Frac %= Pow10[i];
                putc('0' + (char)Digit, File);
            }
        }
    }
    else
        fprintf(File, "%lu", (unsigned long)(Value + 0.5f));
}


/**
 * @brief Get the number of digits used for writing floats
 * @param Pen the pen object
 * @return number of digits
 *
 * Gets the number of digits used in PenWriteDouble().
 */
int PenGetNumDigits(pen *Pen)
{
    return (Pen) ? Pen->NumDigits : PEN_DEFAULT_NUMDIGITS;
}


/**
 * @brief Set the number of digits used for writing floats
 * @param Pen the pen object
 * @param NumDigits the number of digits
 *
 * Sets the number of digits used in PenWriteDouble().
 */
void PenSetNumDigits(pen *Pen, int NumDigits)
{
    if(Pen)
        Pen->NumDigits = (NumDigits < PEN_MAXDIGITS)
            ? NumDigits : PEN_MAXDIGITS;
}

                            
/**
 * @brief Get the pen canvas transformation
 * @param Pen the pen object
 * @return pentrans canvas transformation
 *
 * Gets pen::Trans.  This is useful in combination with PenSetTrans() to save
 * and later restore the canvas transformation.
 */
pentrans PenGetTrans(pen *Pen)
{
    return (Pen) ? Pen->Trans : IdentityPenTrans;
}


/**
 * @brief (Re)set the pen canvas transformation
 * @param Pen the pen object
 * @param Trans the pentrans canvas transformation
 *
 * Sets pen::Trans to Trans.  This is useful in combination with PenGetTrans()
 * to save and later restore the canvas transformation.
 */
void PenSetTrans(pen *Pen, pentrans Trans)
{
    if(Pen)
        Pen->Trans = Trans;
}


/**
 * @brief Transform the pen's canvas
 * @param Pen the pen object
 * @param Trans the pentrans canvas transformation
 */
void PenTransformCanvas(pen *Pen, double a, double b,
    double c, double d, double e, double f)
{
    if(Pen)
    {
        pentrans NewTrans;
        
        NewTrans.a = Pen->Trans.a*a
            + Pen->Trans.c*b;
        NewTrans.b = Pen->Trans.b*a
            + Pen->Trans.d*b;
        NewTrans.c = Pen->Trans.a*c
            + Pen->Trans.c*d;
        NewTrans.d = Pen->Trans.b*c
            + Pen->Trans.d*d;
        
        NewTrans.e = Pen->Trans.e
            + Pen->Trans.a*e
            + Pen->Trans.c*f;
        NewTrans.f = Pen->Trans.f
            + Pen->Trans.b*e
            + Pen->Trans.d*f;
            
        Pen->Trans = NewTrans;
    }
}


/**
 * @brief Translate the pen's canvas
 * @param Pen the pen object
 * @param tx, ty translation
 */
void PenTranslateCanvas(pen *Pen, double tx, double ty)
{
    if(Pen)
        PenTransformCanvas(Pen, 1, 0, 0, 1, tx, ty);
}


/**
 * @brief Scale the pen's canvas
 * @param Pen the pen object
 * @param XScale, YScale scale factors
 */
void PenScaleCanvas(pen *Pen, double XScale, double YScale)
{
    if(Pen)
        PenTransformCanvas(Pen, XScale, 0, 0, YScale, 0, 0);
}


/**
 * @brief Rotate the pen's canvas
 * @param Pen the pen object
 * @param Theta counter-clockwise rotation in radians
 */
void PenRotateCanvas(pen *Pen, double Theta)
{
    if(Pen)
    {
        const double CosTheta = cos(Theta);
        const double SinTheta = sin(Theta);
        PenTransformCanvas(Pen, CosTheta, -SinTheta, SinTheta, CosTheta, 0, 0);
    }
}


/**
 * @brief Horizontally skew the pen's canvas
 * @param Pen the pen object
 * @param Skew skew amount
 */
void PenXSkewCanvas(pen *Pen, double Skew)
{
    if(Pen)
        PenTransformCanvas(Pen, 1, 0, Skew, 1, 0, 0);
}


/**
 * @brief Vertically skew the pen's canvas
 * @param Pen the pen object
 * @param Skew skew amount
 */
void PenYSkewCanvas(pen *Pen, double Skew)
{
    if(Pen)
        PenTransformCanvas(Pen, 1, Skew, 0, 1, 0, 0);
}


/**
 * @brief Find the pen palette index of a color
 * @param Pen the pen object
 * @param Color color to find
 * @return palette index, or -1 if not in palette
 */
int PenColorToIndex(pen *Pen, const float *Color)
{
    int i;
    
    if(!Pen || !(Pen->Palette))
        return -1;
    
    for(i = 0; i < Pen->NumPalette; i++)
        if(Pen->Palette[i][0] == Color[0]
            && Pen->Palette[i][1] == Color[1]
            && Pen->Palette[i][2] == Color[2])
            return i;
    
    return -1;
}


/**
 * @brief Add a color to a pen's palette
 * @param Pen the pen object
 * @param Color color to add
 * @return 1 on success, 0 on failure
 *
 * Color is added to Pen->Palette if it is not already in the palette.  The
 * palette is reallocated if more space is needed.
 */
int PenDefineColor(pen *Pen, const float *Color)
{
    const int NumPalette = Pen->NumPalette;
    int i;
    
    if(!Pen || !Color)
        return 0;
    
    /* Check if Color is already in Pen->Palette */
    if(PenColorToIndex(Pen, Color) >= 0)
        return 1;
    
    /* Allocate space for the new color if necessary */
    if(NumPalette >= Pen->CapacityPalette)
    {
        float (*NewPalette)[3];
        int NewCapacity = NumPalette + 1 + (NumPalette + 5)/10;
        
        if(!(NewPalette = (float (*)[3])Realloc(Pen->Palette,
            sizeof(float)*3*NewCapacity)))
            return 0;
        
        Pen->Palette = NewPalette;
        Pen->CapacityPalette = NewCapacity;
    }
    
    for(i = 0; i < 3; i++)
        Pen->Palette[NumPalette][i] = Color[i];
    
    (Pen->NumPalette)++;
    return 1;
}


/**
 * @brief Set the current pen drawing color
 * @param Pen the pen object
 * @param Color the drawing color
 * @return the color's palette index, or -1 if not in palette
 */
int PenSetColor(pen *Pen, const float *Color)
{
    if(Pen && Color)
    {
        Pen->Color[0] = Color[0];
        Pen->Color[1] = Color[1];
        Pen->Color[2] = Color[2];
        
        if(Pen->SetColor)
            Pen->SetColor(Pen, Color);
    }
    
    return (Pen->ColorIndex = PenColorToIndex(Pen, Color));
}


/**
 * @brief Compute canvas transformation for an ellipse
 * @param x, y, rx, ry, Theta parameters for pen::DrawEllipse
 * @param Trans the canvas transformation
 *
 * The input values of the parameters x, y, rx, ry, Theta are the untransformed
 * parameters as given to pen::DrawEllipse.  They are replaced with transformed
 * values.
 */
void TransformEllipse(float *x, float *y, float *rx, float *ry, float *Theta,
    pentrans Trans)
{
    const double CosTheta = cos(*Theta), SinTheta = sin(*Theta);
    double XTrans, YTrans, A[2][2], ThetaTrans, RxTrans, RyTrans, Phi;
    int Sign1, Sign2;
    
    A[0][0] = (Trans.a*CosTheta - Trans.c*SinTheta)*(*rx);
    A[1][0] = (Trans.b*CosTheta - Trans.d*SinTheta)*(*rx);
    A[0][1] = (Trans.a*SinTheta + Trans.c*CosTheta)*(*ry);
    A[1][1] = (Trans.b*SinTheta + Trans.d*CosTheta)*(*ry);
    
    /* Compute singular value decomposition of A */
    Svd2x2(&ThetaTrans, &RxTrans, &RyTrans, &Phi, &Sign1, &Sign2, A);
    
    /* If ThetaTrans = +/- pi/2, swap x and y so that the rotation is zero. */
    if(fabs(fabs(ThetaTrans) - M_PI/2) < 1e-12)
    {
        *rx = (float)RyTrans;
        *ry = (float)RxTrans;
        *Theta = 0;
    }
    else
    {
        *rx = (float)RxTrans;
        *ry = (float)RyTrans;
        *Theta = (float)ThetaTrans;
    }
    
    XTrans = Trans.a*(*x) + Trans.c*(*y) + Trans.e;
    YTrans = Trans.b*(*x) + Trans.d*(*y) + Trans.f;
    *x = (float)XTrans;
    *y = (float)YTrans;
}


/**
 * @brief Estimate arc length of a quadratic Bezier curve
 * @param Pen the pen object
 * @param x1, y1, x2, y2, x3, y3 starting point, control point, and endpoint
 * @return estimated arc length
 *
 * Computes a crude overestimate of the arc length.
 */
float EstArcLenQBezier(float x1, float y1, float x2, float y2,
    float x3, float y3, pentrans Trans)
{
    float Temp;
    
    Temp = Trans.a*x1 + Trans.c*y1 + Trans.e;
    y1 = Trans.b*x1 + Trans.d*y1 + Trans.f;
    x1 = Temp;
    Temp = Trans.a*x2 + Trans.c*y2 + Trans.e;
    y2 = Trans.b*x2 + Trans.d*y2 + Trans.f;
    x2 = Temp;
    Temp = Trans.a*x3 + Trans.c*y3 + Trans.e;
    y3 = Trans.b*x3 + Trans.d*y3 + Trans.f;
    x3 = Temp;
    
    /* Crude overestimate of the arc length */
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2))
        + sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3));
}


/**
 * @brief Estimate the arc length of an ellipse
 * @param x, y, rx, ry, Theta parameters for pen::DrawEllipse
 * @param Trans the canvas transformation
 * @return estimated arc length
 */
float EstArcLenEllipse(float x, float y, float rx, float ry, float Theta,
    pentrans Trans)
{
    float Temp;
    
    TransformEllipse(&x, &y, &rx, &ry, &Theta, Trans);
    Temp = (rx - ry)/(rx + ry);
    Temp *= Temp;
    return M_PI*(rx + ry)*(1 + 3*Temp/(10 + sqrt(4 - 3*Temp)));
}



/*****************************************************************************
 * Raster-specific code                                                      *
 *****************************************************************************/


/** @brief Swap A and B */
#define SWAP(A,B)   \
    Temp = A; \
    A = B; \
    B = Temp;


static float fpart(float x)
{
    return (float)(x - floor(x));
}


static void PlotPixel(float *Image, int Width, int Height, int x, int y,
    float Alpha, const float *Color)
{
    if(0 <= x && x < Width && 0 <= y && y < Height)
    {
        const float CAlpha = 1 - Alpha;
        
        Image += 3*(x + Width*y);
        Image[0] = CAlpha*Image[0] + Alpha*Color[0];
        Image[1] = CAlpha*Image[1] + Alpha*Color[1];
        Image[2] = CAlpha*Image[2] + Alpha*Color[2];
    }
}


/**
 * @brief Draw an anti-aliased line in a raster image
 * @param Pen the pen object
 * @param x1, y1, x2, y2 endpoints of the line
 * @return 1 on success, 0 on failure
 */
static int RasterDrawLine(pen *Pen, float x1, float y1, float x2, float y2)
{
    float *Image, *Color;
    float yend, Gap, dx, dy, Gradient, y, Temp;
    int Width, Height, xend, ix, iy, ix1, iy1, ix2, iy2, Swapped = 0;
    
            
    if(!Pen || !(Image = Pen->Image))
        return 0;
    
    Width = Pen->ImageWidth;
    Height = Pen->ImageHeight;
    Color = Pen->Color;
    
    /* Transform endpoint coordinates */
    Temp = Pen->Trans.a*x1 + Pen->Trans.c*y1 + Pen->Trans.e;
    y1 = Pen->Trans.b*x1 + Pen->Trans.d*y1 + Pen->Trans.f;
    x1 = Temp;
    Temp = Pen->Trans.a*x2 + Pen->Trans.c*y2 + Pen->Trans.e;
    y2 = Pen->Trans.b*x2 + Pen->Trans.d*y2 + Pen->Trans.f;
    x2 = Temp;
    
    dx = x2 - x1;
    dy = y2 - y1;
    
    if(fabs(dx) < fabs(dy))
    {
        SWAP(x1,y1)
        SWAP(x2,y2)
        SWAP(dx,dy)
        Swapped = 1;
    }
        
    if(x2 < x1)
    {
        SWAP(x1,x2)
        SWAP(y1,y2)
    }
    
    if(dx == 0)
        return 1;
    
    Gradient = dy/dx;
    
    xend = (int)floor(x1 + 0.5f);
    yend = y1 + Gradient * (xend - x1);
    y = yend + Gradient;
    Gap = 1 - fpart(x1 + 0.5f);
    ix1 = xend;
    iy1 = (int)floor(yend);
    
    if(!Swapped)
    {
        PlotPixel(Image, Width, Height, ix1, iy1, (1 - fpart(yend)) * Gap, Color);
        PlotPixel(Image, Width, Height, ix1, iy1 + 1, fpart(yend) * Gap, Color);
    }
    else
    {
        PlotPixel(Image, Width, Height, iy1, ix1, (1 - fpart(yend)) * Gap, Color);
        PlotPixel(Image, Width, Height, iy1 + 1, ix1, fpart(yend) * Gap, Color);
    }
        
    xend = (int)floor(x2 + 0.5f);
    yend = y2 + Gradient * (xend - x2);
    Gap = fpart(x2 + 0.5f);
    ix2 = xend;
    iy2 = (int)floor(yend);
    
    if(!Swapped)
    {
        PlotPixel(Image, Width, Height, ix2, iy2, (1 - fpart(yend)) * Gap, Color);
        PlotPixel(Image, Width, Height, ix2, iy2 + 1, fpart(yend) * Gap, Color);
        
        for(ix = ix1 + 1; ix < ix2; ix++)
        {
            iy = (int)y;
            PlotPixel(Image, Width, Height, ix, iy, 1 - fpart(y), Color);
            PlotPixel(Image, Width, Height, ix, iy + 1, fpart(y), Color);
            y += Gradient;
        }
    }
    else
    {
        PlotPixel(Image, Width, Height, iy2, ix2, (1 - fpart(yend)) * Gap, Color);
        PlotPixel(Image, Width, Height, iy2 + 1, ix2, fpart(yend) * Gap, Color);
        
        for(ix = ix1 + 1; ix < ix2; ix++)
        {
            iy = (int)y;
            PlotPixel(Image, Width, Height, iy, ix, 1 - fpart(y), Color);
            PlotPixel(Image, Width, Height, iy + 1, ix, fpart(y), Color);
            y += Gradient;
        }
    }
    
    return 1;
}


/**
 * @brief Set pen to render to a raster image
 * @param Pen the pen object
 * @param Image the row-major interleaved RGB image to render to
 * @param Width, Height dimensions of Image
 * @return 1 on success, 0 on failure
 */
int PenRenderToImage(pen *Pen, float *Image, int Width, int Height)
{
    if(!Pen || !Image || Width <= 0 || Height <= 0)
        return 0;
    
    Pen->Image = Image;
    Pen->SetColor = NULL;
    Pen->ImageWidth = Width;
    Pen->ImageHeight = Height;
    Pen->DrawLine = RasterDrawLine;
    return 1;
}



/*****************************************************************************
 * EPS-specific code                                                         *
 *****************************************************************************/

/**
 * @brief Set drawing color in an EPS
 * @param Pen the pen object
 * @param Color drawing color
 */
void EpsSetColor(pen *Pen, const float *Color)
{
    FILE *File;
    
    if(Pen && (File = Pen->File))
    {
        PenWriteDouble(File, Color[0], Pen->NumDigits);
        putc(' ', File);
        PenWriteDouble(File, Color[1], Pen->NumDigits);
        putc(' ', File);
        PenWriteDouble(File, Color[2], Pen->NumDigits);
        fputs(" setrgbcolor\n", File);
    }
}


/**
 * @brief Draw a line in an EPS
 * @param Pen the pen object
 * @param x1, y1, x2, y2 endpoints of the line
 * @return 1 on success, 0 on failure
 */
static int EpsDrawLine(pen *Pen, float x1, float y1, float x2, float y2)
{
    FILE *File;
    double xt1, yt1, xt2, yt2;
    int NumDigits;
        
    if(!Pen || !(File = Pen->File))
        return 0;
    
    NumDigits = Pen->NumDigits;
    
    xt1 = Pen->Trans.a*x1 + Pen->Trans.c*y1 + Pen->Trans.e;
    yt1 = Pen->Trans.b*x1 + Pen->Trans.d*y1 + Pen->Trans.f;
    xt2 = Pen->Trans.a*x2 + Pen->Trans.c*y2 + Pen->Trans.e;
    yt2 = Pen->Trans.b*x2 + Pen->Trans.d*y2 + Pen->Trans.f;
    
    PenWriteDouble(File, xt2 - xt1, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, yt2 - yt1, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, xt1, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, yt1, NumDigits);
    fputs(" m1\n", File);
    return 1;
}


/**
 * @brief Draw a quadratic Bezier curve in an EPS
 * @param Pen the pen object
 * @param x1, y1, x2, y2, x3, y3 starting point, control point, and endpoint
 * @return 1 on success, 0 on failure
 */
static int EpsDrawQBezier(pen *Pen, float x1, float y1,
    float x2, float y2, float x3, float y3)
{
    FILE *File;
    double xt1, yt1, xt2, yt2, xt3, yt3;
    int NumDigits;
        
    if(!Pen || !(File = Pen->File))
        return 0;
    
    NumDigits = Pen->NumDigits;
    xt1 = Pen->Trans.a*x1 + Pen->Trans.c*y1 + Pen->Trans.e;
    yt1 = Pen->Trans.b*x1 + Pen->Trans.d*y1 + Pen->Trans.f;
    xt2 = Pen->Trans.a*x2 + Pen->Trans.c*y2 + Pen->Trans.e;
    yt2 = Pen->Trans.b*x2 + Pen->Trans.d*y2 + Pen->Trans.f;
    xt3 = Pen->Trans.a*x3 + Pen->Trans.c*y3 + Pen->Trans.e;
    yt3 = Pen->Trans.b*x3 + Pen->Trans.d*y3 + Pen->Trans.f;
    
    PenWriteDouble(File, (2.0/3.0)*(xt2 - xt1), NumDigits);
    putc(' ', File);
    PenWriteDouble(File, (2.0/3.0)*(yt2 - yt1), NumDigits);
    putc(' ', File);
    PenWriteDouble(File, (xt3 - xt1) + (2.0/3.0)*(xt2 - xt3), NumDigits);
    putc(' ', File);
    PenWriteDouble(File, (yt3 - yt1) + (2.0/3.0)*(yt2 - yt3), NumDigits);
    putc(' ', File);
    PenWriteDouble(File, xt3 - xt1, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, yt3 - yt1, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, xt1, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, yt1, NumDigits);
    fputs(" m2\n", File);
    return 1;
}


/**
 * @brief Draw a quadratic Bezier curve in an EPS
 * @param Pen the pen object
 * @param x, y center of the ellipse
 * @param rx, ry radii of the ellipse
 * @param Theta ellipse rotation
 * @return 1 on success, 0 on failure
 */
static int EpsDrawEllipse(pen *Pen, float x, float y,
    float rx, float ry, float Theta)
{
    FILE *File;
    int NumDigits;
        
    if(!Pen || !(File = Pen->File))
        return 0;
    
    NumDigits = Pen->NumDigits;
    TransformEllipse(&x, &y, &rx, &ry, &Theta, Pen->Trans);
        
    PenWriteDouble(File, rx, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, ry, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, Theta*(180/M_PI), NumDigits);
    putc(' ', File);
    PenWriteDouble(File, x, NumDigits);
    putc(' ', File);
    PenWriteDouble(File, y, NumDigits);
    fputs(" m3\n", File);
    return 1;
}


/**
 * @brief Start writing a new EPS file
 * @param Pen the pen object
 * @param FileName the EPS file name
 * @param Width,Height dimensions of the SVG canvas
 * @return 1 on success, 0 on failure
 */
int EpsOpen(pen *Pen, const char *FileName, float Width, float Height)
{
    FILE *File = NULL;
    
    if(!(File = fopen(FileName, "wb")))
    {
        ErrorMessage("Failed to open for writing \"%s\".\n", FileName);
        return 0;
    }
    
    /* Tell File to use buffering */
    setvbuf(File, 0, _IOFBF, FILE_BUFFER_CAPACITY);
    
    fputs("%%!PS-Adobe-2.0\n%%%%BoundingBox: 0 0 ", File);
    PenWriteDouble(File, Width, Pen->NumDigits);
    putc(' ', File);
    PenWriteDouble(File, Height, Pen->NumDigits);
    putc('\n', File);
    
    /* Set PDF distiller info */
    fputs("<< /PageSize [", File);
    PenWriteDouble(File, Width, Pen->NumDigits);
    putc(' ', File);
    PenWriteDouble(File, Height, Pen->NumDigits);
    fputs("] >> setpagedevice\n", File);
    fputs("systemdict /setdistillerparams known {\n"
        "<< /AutoFilterGrayImages false /GrayImageFilter /FlateEncode "
        "/AutoFilterColorImages false /ColorImageFilter /FlateEncode "
        ">> setdistillerparams\n"
        "} if\n", File);
    /* Set PS macros and line width */
    fputs("/bdef {bind def} bind def\n"
        "/m1 {moveto rlineto stroke} bdef\n"
        "/m2 {moveto rcurveto stroke} bdef\n"
        "/m3 {/savematrix matrix currentmatrix def\n"
        "translate rotate scale 0 0 1 0 360 arc\n"
        "savematrix setmatrix stroke} bdef\n"
        "0.2 setlinewidth\n", File);
        
    if(ferror(File))
    {
        ErrorMessage("Failed to write \"%s\".\n", FileName);
        goto Catch;
    }
    
    Pen->File = File;
    Pen->SetColor = EpsSetColor;
    Pen->DrawLine = EpsDrawLine;
    Pen->DrawQBezier = EpsDrawQBezier;
    Pen->DrawEllipse = EpsDrawEllipse;
    return 1;
Catch:
    if(File)
        fclose(File);
    return 0;
}


/**
 * @brief Finish writing an EPS file
 * @param Pen the pen object
 * @return 1 on success, 0 on failure
 */
int EpsClose(pen *Pen)
{
        if(!Pen || !(Pen->File)
        || fputs("showpage\n", Pen->File) == EOF
        || ferror(Pen->File)
        || fclose(Pen->File))
        return 0;
    
    return 1;
}


/**
 * @brief Write ASCII85 encoded data
 * @param File the destination stream or File
 * @param Data binary data to encode
 * @param NumBytes number of bytes in Data
 * @return 1 on success, 0 on failure
 */
int WriteASCII85(FILE *File, const uint8_t *Data, int NumBytes)
{
    unsigned long Tuple, Plain[4];
    unsigned int Encoded[5];
    int i, k, LineCount, Padding;
    
    
    /* Write ASCII85-encoded data */
    for(i = 0, LineCount = 0; i < NumBytes; i += 4)
    {
        for(k = 0; k < 4; k++)  /* Get four bytes */
            Plain[k] = Data[i + k];
        
        Tuple = (Plain[0] << 24) | (Plain[1] << 16)
            | (Plain[2] << 8) | Plain[3];
        
        for(k = 4; k >= 0; k--) /* Convert to radix 85 */
        {
            Encoded[k] = Tuple % 85;
            Tuple /= 85;
        }
        
        for(k = 0; k < 5; k++)  /* Write ASCII85 tuple */
            fputc(33 + Encoded[k], File);
        
        /* Periodically emit newlines */
        if(++LineCount >= 15)
        {
            LineCount = 0;
            
            if(fprintf(File, "\n") < 0)
                return 0;
        }
    }
    
    /* Write final tuple */
    if(i < NumBytes)
    {
        for(k = 0; i + k < NumBytes; k++)
            Plain[k] = Data[i + k];
        
        for(Padding = 0; k < 4; k++, Padding++)
            Plain[k] = 0;
        
        Tuple = (Plain[0] << 24) | (Plain[1] << 16)
            | (Plain[2] << 8) | Plain[3];
        
        for(k = 4; k >= 0; k--) /* Convert to radix 85 */
        {
            Encoded[k] = Tuple % 85;
            Tuple /= 85;
        }
        
        for(k = 0; k < 5 - Padding; k++)  /* Write ASCII85 tuple */
            fputc(33 + Encoded[k], File);
        
        if(++LineCount >= 15
            && fprintf(File, "\n") < 0)
            return 0;
    }
    
    if(fprintf(File, "~>\n") < 0 || ferror(File))
        return 0;
    
    return 1;
}


/**
 * @brief Writes a grayscale image to an EPS file
 * @param File the stdio file handle to write to
 * @param Image grayscale image data
 * @param Width, Height the image dimensions
 *
 * The image is plotted on the canvas in the rectangle
 *     [0,Width] x [0,Height]
 * where the lower-left corner is at the origin.
 *
 * This routine only writes the image data.  EpsOpen should be called first,
 * then this routine and other drawing commands, and finally EpsClose.
 *
 * For relative simplicity, the image data is written uncompressed in
 * ASCII85 encoding.  The file size is approximately 25% larger than in the PGM
 * file format.
 */
int EpsWriteGrayImage(FILE *File,
    const uint8_t *Image, int Width, int Height)
{
    /* Specify ASCII85 sRGB 8-bit color image data */
    if(!File || !Image || fprintf(File,
        "gsave\n"
        "/DeviceGray setcolorspace\n"
        "0 %d translate %d %d scale\n"
        "<< /ImageType 1\n"
        "   /Width %d\n"
        "   /Height %d\n"
        "   /ImageMatrix [%d 0 0 -%d 0 0]\n"
        "   /BitsPerComponent 8\n"
        "   /Decode [0 1]\n"
        "   /DataSource currentfile /ASCII85Decode filter\n"
        "   /Interpolate false\n"
        ">> image\n",
        Height, Width, Height,
        Width, Height, Width, Height) < 0
        || !WriteASCII85(File, Image, Width*Height)
        || fprintf(File, "grestore\n") < 0)
        return 0;
    else
        return 1;
}


/**
 * @brief Writes an RGB color image to an  EPS file
 * @param File the stdio file handle to write to
 * @param Image interleaved RGB color image data
 * @param Width, Height the image dimensions
 *
 * The image is plotted on the canvas in the rectangle
 *     [0,Width] x [0,Height]
 * where the lower-left corner is at the origin.
 *
 * This routine only writes the image data.  EpsOpen should be called first,
 * then this routine and other drawing commands, and finally EpsClose.
 *
 * For relative simplicity, the image data is written uncompressed in
 * ASCII85 encoding.  The file size is approximately 25% larger than in the PGM
 * file format.
 */
int EpsWriteColorImage(FILE *File,
    const uint8_t *Image, int Width, int Height)
{
    /* Specify ASCII85 sRGB 8-bit color image data */
    if(!File || !Image || fprintf(File,
        "gsave\n"
        "/DeviceRGB setcolorspace\n"
        "0 %d translate %d %d scale\n"
        "<< /ImageType 1\n"
        "   /Width %d\n"
        "   /Height %d\n"
        "   /ImageMatrix [%d 0 0 -%d 0 0]\n"
        "   /BitsPerComponent 8\n"
        "   /Decode [0 1 0 1 0 1]\n"
        "   /DataSource currentfile /ASCII85Decode filter\n"
        "   /Interpolate false\n"
        ">> image\n",
        Height, Width, Height,
        Width, Height, Width, Height) < 0
        || !WriteASCII85(File, Image, 3*Width*Height)
        || fprintf(File, "grestore\n") < 0)
        return 0;
    else
        return 1;
}



/*****************************************************************************
 * SVG-specific code                                                         *
 *****************************************************************************/

/** @brief Write a color in hex web notation */
static void WriteHexColor(FILE *File, const float *Color)
{
    fprintf(File, "#%02X%02X%02X",
        (int)(255*Color[0] + 0.5f),
        (int)(255*Color[1] + 0.5f),
        (int)(255*Color[2] + 0.5f));
}


/**
 * @brief Draw a line in an SVG
 * @param Pen the pen object
 * @param x1, y1, x2, y2 endpoints of the line
 * @return 1 on success, 0 on failure
 */
static int SvgDrawLine(pen *Pen, float x1, float y1, float x2, float y2)
{
    FILE *File;
    double xt, yt;
    int NumDigits;
        
    if(!Pen || !(File = Pen->File))
        return 0;
    
    NumDigits = Pen->NumDigits;
    fputs("<line ", File);
    
    if(Pen->ColorIndex >= 0)
        fprintf(File, "class=\"%c\" ", 'a' + Pen->ColorIndex);
    else
    {
        fputs("style=\"stroke:", File);
        WriteHexColor(File, Pen->Color);
        fputs("\" ", File);
    }
    
    xt = Pen->Trans.a*x1 + Pen->Trans.c*y1 + Pen->Trans.e;
    yt = Pen->Trans.b*x1 + Pen->Trans.d*y1 + Pen->Trans.f;
    fputs("x1=\"", File);
    PenWriteDouble(File, xt, NumDigits);
    fputs("\" y1=\"", File);
    PenWriteDouble(File, yt, NumDigits);
    
    xt = Pen->Trans.a*x2 + Pen->Trans.c*y2 + Pen->Trans.e;
    yt = Pen->Trans.b*x2 + Pen->Trans.d*y2 + Pen->Trans.f;
    fputs("\" x2=\"", File);
    PenWriteDouble(File, xt, NumDigits);
    fputs("\" y2=\"", File);
    PenWriteDouble(File, yt, NumDigits);
    fputs("\" />\n", File);
    return 1;
}


/**
 * @brief Draw a quadratic Bezier curve in an SVG
 * @param Pen the pen object
 * @param x1, y1, x2, y2, x3, y3 starting point, control point, and endpoint
 * @return 1 on success, 0 on failure
 */
static int SvgDrawQBezier(pen *Pen,
    float x1, float y1, float x2, float y2, float x3, float y3)
{
    FILE *File;
    double xt0, yt0, xt, yt;
    int NumDigits;
        
    if(!Pen || !(File = Pen->File))
        return 0;
    
    NumDigits = Pen->NumDigits;
    fputs("<path ", File);
    
    if(Pen->ColorIndex >= 0)
        fprintf(File, "class=\"%c\" ", 'a' + Pen->ColorIndex);
    else
    {
        fputs("style=\"stroke:", File);
        WriteHexColor(File, Pen->Color);
        fputs("; fill:none\" ", File);
    }
        
    xt0 = Pen->Trans.a*x1 + Pen->Trans.c*y1 + Pen->Trans.e;
    yt0 = Pen->Trans.b*x1 + Pen->Trans.d*y1 + Pen->Trans.f;
    fputs("d=\"M", File);
    PenWriteDouble(File, xt0, NumDigits);
    putc(',', File);
    PenWriteDouble(File, yt0, NumDigits);
    fputs(" q", File);
    
    xt = Pen->Trans.a*x2 + Pen->Trans.c*y2 + Pen->Trans.e;
    yt = Pen->Trans.b*x2 + Pen->Trans.d*y2 + Pen->Trans.f;
    PenWriteDouble(File, xt - xt0, NumDigits);
    putc(',', File);
    PenWriteDouble(File, yt - yt0, NumDigits);
    putc(' ', File);
    
    xt = Pen->Trans.a*x3 + Pen->Trans.c*y3 + Pen->Trans.e;
    yt = Pen->Trans.b*x3 + Pen->Trans.d*y3 + Pen->Trans.f;
    PenWriteDouble(File, xt - xt0, NumDigits);
    putc(',', File);
    PenWriteDouble(File, yt - yt0, NumDigits);
    fputs("\" />\n", File);
    return 1;
}


/**
 * @brief Draw an ellipse in an SVG
 * @param Pen the pen object
 * @param x, y center of the ellipse
 * @param rx, ry radii of the ellipse
 * @param Theta ellipse rotation
 * @return 1 on success, 0 on failure
 */
static int SvgDrawEllipse(pen *Pen,
    float x, float y, float rx, float ry, float Theta)
{
    FILE *File;
    int NumDigits, IsCircle;
        
    if(!Pen || !(File = Pen->File))
        return 0;
    
    NumDigits = Pen->NumDigits;
    TransformEllipse(&x, &y, &rx, &ry, &Theta, Pen->Trans);
    IsCircle = (fabs(rx - ry) < 1e-2);
    
    if(IsCircle)
        fputs("<circle ", File);
    else
        fputs("<ellipse ", File);
    
    if(Pen->ColorIndex >= 0)
        fprintf(File, "class=\"%c\" ", 'a' + Pen->ColorIndex);
    else
    {
        fputs("style=\"stroke:", File);
        WriteHexColor(File, Pen->Color);
        fputs("\" ", File);
    }
    
    Theta *= 180/M_PI;
    
    if(!IsCircle && fabs(Theta) >= 1e-2 && abs(Theta - 180) >= 1e-2)
    {
        fputs("transform=\"rotate(", File);
        PenWriteDouble(File, Theta, NumDigits);
        putc(' ', File);
        PenWriteDouble(File, x, NumDigits);
        putc(' ', File);
        PenWriteDouble(File, y, NumDigits);
        fputs(")\" ", File);
    }
    
    fputs("cx=\"", File);
    PenWriteDouble(File, x, NumDigits);
    fputs("\" cy=\"", File);
    PenWriteDouble(File, y, NumDigits);
    
    if(IsCircle)
        fputs("\" r=\"", File);
    else
    {
        fputs("\" rx=\"", File);
        PenWriteDouble(File, rx, NumDigits);
        fputs("\" ry=\"", File);
    }
    
    PenWriteDouble(File, ry, NumDigits);
    fputs("\" />\n", File);
    return 1;
}


/**
 * @brief Start writing a new SVG file
 * @param Pen the pen object
 * @param FileName the SVG file name
 * @param Width,Height dimensions of the SVG canvas
 * @return 1 on success, 0 on failure
 */
int SvgOpen(pen *Pen, const char *FileName, float Width, float Height)
{
    FILE *File = NULL;
    int i;
        
    if(!FileName || !(File = fopen(FileName, "wb")))
    {
        ErrorMessage("Failed to open for writing \"%s\".\n", FileName);
        goto Catch;
    }
    
    /* Tell File to use buffering */
    setvbuf(File, 0, _IOFBF, FILE_BUFFER_CAPACITY);
        
    /* Write the XML and DOCTYPE header */
    fputs( "<?xml version=\"1.0\" standalone=\"no\"?>\n"
           "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n"
           "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n", File);
    
    /* Write SVG header */
    fputs("<svg width=\"", File);
    PenWriteDouble(File, Width, Pen->NumDigits);
    fputs("\" height=\"", File);
    PenWriteDouble(File, Height, Pen->NumDigits);
    fputs("\" version=\"1.1\"\n"
          "  xmlns=\"http://www.w3.org/2000/svg\""
          " xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n", File);
    
    /* Write colors in Pen->Palette in an internal stylesheet */
    if(Pen->NumPalette > 0)
    {
        fputs("<defs>\n<style type=\"text/css\"><![CDATA[\n", File);
        
        for(i = 0; i < Pen->NumPalette; i++)
        {
            fprintf(File, ".%c {stroke:", 'a' + i);
            WriteHexColor(File, Pen->Palette[i]);
            fputs("; stroke-width:2; fill:none}\n", File);
        }
        
        fputs("]]></style>\n</defs>\n", File);
    }
    
    if(ferror(File))
    {
        ErrorMessage("Failed to write \"%s\".\n", FileName);
        goto Catch;
    }
    
    Pen->File = File;
    Pen->SetColor = NULL;
    Pen->DrawLine = SvgDrawLine;
    Pen->DrawQBezier = SvgDrawQBezier;
    Pen->DrawEllipse = SvgDrawEllipse;
    return 1;
Catch:
    if(File)
        fclose(File);
    return 0;
}


/**
 * @brief Finish writing an SVG file
 * @param Pen the pen object
 * @return 1 on success, 0 on failure
 */
int SvgClose(pen *Pen)
{
    if(!Pen || !(Pen->File)
        || fputs("</svg>\n", Pen->File) == EOF
        || ferror(Pen->File)
        || fclose(Pen->File))
        return 0;
    
    return 1;
}




