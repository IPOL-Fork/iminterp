/**
 * @file drawline.h
 * @brief Draw a black antialiased line
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

#ifndef _DRAWLINE_H_
#define _DRAWLINE_H_

#include <ipol/basic.h>

void DrawLine(uint32_t *Image, int Width, int Height,
    float x1, float y1, float x2, float y2, const float *Color);

#endif /* _DRAWLINE_H_ */
