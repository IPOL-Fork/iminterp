/**
 * @file fitsten.h
 * @brief Select the best-fitting contour stencils
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

#ifndef _FITSTEN_H_
#define _FITSTEN_H_

#include <ipol/basic.h>

int FitStencils(int *Stencil, int32_t *Image, int Width, int Height, int StencilMul);

#endif /* _FITSTEN_H_ */
