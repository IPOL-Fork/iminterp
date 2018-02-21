/**
 * @file svd2x2.h
 * @brief Compute the singular value decomposition of a 2x2 matrix
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

#ifndef _SVD2X2_H_
#define _SVD2X2_H_

void Svd2x2(double *Theta, double *Sigma1, double *Sigma2,
    double *Phi, int *Sign1, int *Sign2, double A[2][2]);

#endif /* _SVD2X2_H_ */
