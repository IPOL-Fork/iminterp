/**
 * @file nninterp.h
 * @brief Nearest neighbor image interpolation
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

#ifndef _NNINTERP_H_
#define _NNINTERP_H_
#include "imageio.h"

void NearestInterp(uint32_t *Output, int OutputWidth, int OutputHeight,
    uint32_t *Input, int InputWidth, int InputHeight,
    float ScaleFactor, int CenteredGrid);

#endif /* _NNINTERP_H_ */
