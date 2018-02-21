/**
 * @file cwinterp.h
 * @brief Contour stencil windowed interpolation
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

#ifndef _CWINTERP_H_
#define _CWINTERP_H_

#include <ipol/basic.h>

/** @brief struct of interpolation parameters */
typedef struct
{
    /** The interpolation scale factor */
    double ScaleFactor;
    int CenteredGrid;
    /** The number of steps of iterative refinement */
    int RefinementSteps;
    /** \f$\sigma_h\f$, the blur size of the point spread function */
    double PsfSigma;
    /** \f$\sigma_\tau\f$, tagential spread of \f$\varphi\f$ */
    double PhiSigmaTangent;
    /** \f$\sigma_\nu\f$, normal spread of \f$\varphi\f$ */
    double PhiSigmaNormal;
} cwparams;


int32_t *PreCWInterp(cwparams Param);

int CWInterp(uint32_t *Output, const uint32_t *Input,
    int InputWidth, int InputHeight, const int32_t *Psi, cwparams Param);

int CWInterpEx(uint32_t *Output, int OutputWidth, int OutputHeight,
    const uint32_t *Input, int InputWidth, int InputHeight,
    const int32_t *Psi, cwparams Param);

int DisplayContours(uint32_t *Output, int OutputWidth, int OutputHeight,
    uint32_t *Input, int InputWidth, int InputHeight, cwparams Param);
    
#endif /* _CWINTERP_H_ */
