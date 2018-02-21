/**
 * @file sinterp.h
 * @brief Contour stencil interpolation
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

#ifndef _SINTERP_H_
#define _SINTERP_H_

#include "sset.h"


/* sinterp is encapsulated by forward declaration */
typedef struct sinterpstruct sinterp;


sinterp *NewSInterp(const sset *StencilSet,
    float RhoSigmaTangent, float RhoSigmaNormal);

void FreeSInterp(sinterp *SInterp);

float *ComputeRhoSamples(const sinterp *SInterp,
    int ScaleFactor, int CenteredGrid);

int Prefilter(float *Filtered,
    const int *Stencil, const float *RhoSamples,
    const float *Input, int InputWidth, int InputHeight,
    int ScaleFactor, int CenteredGrid, float PsfSigma, int NumPasses);

float DeconvResidual(float *Residual,
    const float *Coarse, int CoarseWidth, int CoarseHeight,
    const float *Interp, int InterpWidth, int InterpHeight,
    float ScaleFactor, int CenteredGrid, float PsfSigma);

void IntegerScalePass(float *Output,
    const int *Stencil, const float *RhoSamples,
    const float *Input, int InputWidth, int InputHeight,
    int ScaleFactor, int CenteredGrid);

int ArbitraryScale(float *Output, int OutputWidth, int OutputHeight,
    const int *Stencil, const sinterp *SInterp,
    const float *Input, int InputWidth, int InputHeight,
    float ScaleFactor, int CenteredGrid);

#endif /* _SINTERP_H_ */
