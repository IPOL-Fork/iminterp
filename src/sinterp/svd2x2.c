/**
 * @file svd2x2.c
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

#include <math.h>


/**
 * @brief Compute the singular value decomposition of a 2x2 matrix
 * @param Theta rotation of the U matrix
 * @param Sigma1, Sigma2 singular values
 * @param Phi rotation of the V matrix
 * @param Sign1, Sign2 correction signs
 * @param A input 2x2 matrix
 *
 * Computes the singular value decomposition of A,
 *    \f[ A = U \Sigma V^T. \f]
 * U is a rotation by Theta,
 *    \f[ U = \begin{pmatrix}
 *             \cos \theta & -\sin \theta \\
 *             \sin \theta &  \cos \theta
 *            \end{pmatrix}. \f]
 * Sigma is the diagonal matrix with entries Sigma1 and Sigma2, which satisfy
 * Sigma1 >= Sigma2 >= 0.  V is in general not representable by a rotation
 * alone and is further decomposed as a rotation matrix W and a diagonal sign
 * correction matrix C so that
 *    \f[ A = U \Sigma C W^T. \f]
 * W is a rotation by Phi,
 *    \f[ W = \begin{pmatrix}
 *             \cos \phi & -\sin \phi \\
 *             \sin \phi &  \cos \phi
 *            \end{pmatrix}, \f]
 * and C is
 *    \f[ C = \begin{pmatrix}
 *             \text{Sign1} & 0 \\
 *             0            & \text{Sign2}
 *            \end{pmatrix}, \f]
 * where Sign1 and Sign2 are signs, values of either +1 or -1.
 *
 * The matrices U and V can be computed from Theta, Phi, Sign1, Sign2 as
@code
    U[0][0] = cos(Theta);
    U[1][0] = sin(Theta);
    U[0][1] = -sin(Theta);
    U[1][1] = code(Theta);
    
    V[0][0] = Sign1*cos(Phi);
    V[1][0] = Sign1*sin(Phi);
    V[0][1] = -Sign2*sin(Phi);
    V[1][1] = Sign2*cos(Phi);
@endcode
 *
 * The algorithm used here is based on the notes
 *    <http://www.ualberta.ca/~mlipsett/ENGM541/Readings/svd_ellis.pdf>.
 * The derivation is, first, by expanding the following in terms of the
 * components
 *    \f[ \begin{pmatrix}
 *         a & b \\
 *         b & c
 *        \end{pmatrix} = A A^T = U \Sigma^2 U^T, \f]
 * one can show
 *    \f[ \begin{aligned}
 *        2b &= (\sin 2\theta)(\sigma_1^2 - \sigma_2^2), \\
 *     a - c &= (\cos 2\theta)(\sigma_1^2 - \sigma_2^2).
 *        \end{aligned} \f]
 * Therefore,
 *    \f[ \theta = \tfrac{1}{2} \mathrm{atan2}(2b, a - c), \f]
 * and similarly for \f$ \phi \f$ with \f$ A^T A = W \Sigma^2 W^T \f$.
 *
 * Second, now that U and W are known, the singular values and signs are
 * obtained as
 *    \f[ \Sigma C = U^T A W \f]
 */
void Svd2x2(double *Theta, double *Sigma1, double *Sigma2,
    double *Phi, int *Sign1, int *Sign2, double A[2][2])
{
    double Temp, U[2][2], W[2][2];
    
    *Theta = 0.5*atan2(2*(A[0][0]*A[1][0] + A[0][1]*A[1][1]),
        (A[0][0]*A[0][0] + A[0][1]*A[0][1])
        - (A[1][0]*A[1][0] + A[1][1]*A[1][1]));
    *Phi = 0.5*atan2(2*(A[0][0]*A[0][1] + A[1][0]*A[1][1]),
        (A[0][0]*A[0][0] + A[1][0]*A[1][0])
        - (A[0][1]*A[0][1] + A[1][1]*A[1][1]));
    
    U[0][0] = cos(*Theta);
    U[1][0] = sin(*Theta);
    U[0][1] = -U[1][0];
    U[1][1] = U[0][0];
    
    W[0][0] = cos(*Phi);
    W[1][0] = sin(*Phi);
    W[0][1] = -W[1][0];
    W[1][1] = W[0][0];
    
    Temp = U[0][0]*(A[0][0]*W[0][0] + A[0][1]*W[1][0])
            + U[1][0]*(A[1][0]*W[0][0] + A[1][1]*W[1][0]);
    *Sigma1 = fabs(Temp);
    *Sign1 = (Temp < 0) ? -1 : 1;
    
    Temp = U[0][1]*(A[0][0]*W[0][1] + A[0][1]*W[1][1])
            + U[1][1]*(A[1][0]*W[0][1] + A[1][1]*W[1][1]);
    *Sigma2 = fabs(Temp);
    *Sign2 = (Temp < 0) ? -1 : 1;
}
