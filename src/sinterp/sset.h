/**
 * @file sset.h
 * @brief Contour stencil set construction
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

#ifndef _SSET_H_
#define _SSET_H_

#include <ipol/basic.h>
#include "pen.h"

/* sset is encapsulated by forward declaration */
typedef struct ssetstruct sset;


sset *NewStencilSet(double StencilRadius,
    double (*Psf)(double, double, const void*),
    const void *PsfParam, double PsfRadius);

void FreeStencilSet(sset *Set);

int AddStencil(sset *Set, double (*Phi)(double, double, const void*),
    const void *PhiParam, double Rotation,
    int (*DrawFun)(pen *Pen, const void *Param), const void *DrawParam,
    const float *DrawColor);

int GetNumStencils(const sset *Set);
double (*GetPsf(const sset *Set))(double, double, const void*);
const void *GetPsfParam(const sset *Set);
double GetPsfRadius(const sset *Set);
double EvalPhi(const sset *Set, int S, double x, double y);

int StencilConfusion(double *ConfusionMatrix, sset *Set);

void FitStencils(int *Stencil, const sset *Set,
    const float *Image, int ImageWidth, int ImageHeight);

int DrawContours(const char *FileName, int OutputWidth, int OutputHeight,
    double ScaleFactor, const int *Stencil, const sset *Set,
    const float *Input, int InputWidth, int InputHeight);

int DrawContoursEps(const char *EpsFile,
    const int *Stencil, const sset *Set,
    const float *Input, int InputWidth, int InputHeight);

int DrawContoursSvg(const char *SvgFile, const char *BgFile,
    const int *Stencil, const sset *Set,
    const float *Input, int InputWidth, int InputHeight);

#endif  /* _SSET_H_ */
