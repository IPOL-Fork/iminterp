/**
 * @file pen.h
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

#ifndef _PEN_H_
#define _PEN_H_

#include <ipol/basic.h>


/** @brief Default number of digits for writing floats */
#define PEN_DEFAULT_NUMDIGITS       2
/** @brief Maximum number of digits supported by PenWriteFloat() */
#define PEN_MAXDIGITS               8


/** @brief Coordinate transformation matrix */
typedef struct
{
    double a, b, c, d, e, f;
} pentrans;


/** @brief Pen drawing object */
typedef struct penstruct
{
    FILE *File;
    float *Image;
    int ImageWidth;
    int ImageHeight;
    int NumDigits;
    pentrans Trans;
    float (*Palette)[3];
    int NumPalette;
    int CapacityPalette;
    float Color[3];
    int ColorIndex;
    void (*SetColor)(struct penstruct*, const float*);
    int (*DrawLine)(struct penstruct*,
        float, float, float, float);
    int (*DrawQBezier)(struct penstruct*,
        float, float, float, float, float, float);
    int (*DrawEllipse)(struct penstruct*,
        float, float, float, float, float);
} pen;


pen *NewPen();
void FreePen(pen *Pen);

FILE *PenGetFile(pen *Pen);

/* Drawing functions */
int PenDrawLine(pen *Pen, float x1, float y1, float x2, float y2);
int PenDrawRectangle(pen *Pen, float x1, float y1, float x2, float y2);
int PenDrawQBezier(pen *Pen, float x1, float y1,
    float x2, float y2, float x3, float y3);
int PenDrawCircle(pen *Pen, float x, float y, float r);
int PenDrawEllipse(pen *Pen,
    float x, float y, float rx, float ry, float Theta);

/* Canvas transformation */
pentrans PenGetTrans(pen *Pen);
void PenSetTrans(pen *Pen, pentrans Trans);
void PenTransformCanvas(pen *Pen, double a, double b,
    double c, double d, double e, double f);
void PenTranslateCanvas(pen *Pen, double tx, double ty);
void PenScaleCanvas(pen *Pen, double XScale, double YScale);
void PenRotateCanvas(pen *Pen, double Theta);
void PenXSkewCanvas(pen *Pen, double Skew);
void PenYSkewCanvas(pen *Pen, double Skew);

/* Managing drawing colors */
int PenColorToIndex(pen *Pen, const float *Color);
int PenDefineColor(pen *Pen, const float *Color);
int PenSetColor(pen *Pen, const float *Color);

/* Conversion of floating-point values to strings */
void PenWriteDouble(FILE *File, double Value, int NumDigits);
int PenGetNumDigits(pen *Pen);
void PenSetNumDigits(pen *Pen, int NumDigits);

/* Geometric calcuations */
void TransformEllipse(float *x, float *y, float *rx, float *ry, float *Theta,
    pentrans Trans);
float EstArcLenQBezier(float x1, float y1, float x2, float y2,
    float x3, float y3, pentrans Trans);
float EstArcLenEllipse(float x, float y, float rx, float ry, float Theta,
    pentrans Trans);


/* Drawing to a raster image */
int PenRenderToImage(pen *Pen, float *Image, int Width, int Height);

/* Drawing as an EPS */
int EpsOpen(pen *Pen, const char *FileName, float Width, float Height);
int EpsClose(pen *Pen);
int WriteASCII85(FILE *File, const uint8_t *Data, int NumBytes);
int EpsWriteGrayImage(FILE *File,
    const uint8_t *Image, int Width, int Height);
int EpsWriteColorImage(FILE *File,
    const uint8_t *Image, int Width, int Height);

/* Drawing as an SVG */
int SvgOpen(pen *Pen, const char *FileName, float Width, float Height);
int SvgClose(pen *Pen);


extern const pentrans IdentityPenTrans;

#endif /* _PEN_H_ */
