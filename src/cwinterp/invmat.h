/**
 * @file invmat.h
 * @brief Invert matrix through QR decomposition
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

#ifndef _INVMAT_H_
#define _INVMAT_H_

int InvertMatrix(double *InverseData, double *AData, int N);

#endif /* _INVMAT_H_ */
