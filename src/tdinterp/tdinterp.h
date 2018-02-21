/**
 * @file tdinterp.h
 * @brief Roussos-Maragos interpolation by tensor-driven diffusion
 * @author Pascal Getreuer <getreuer@gmail.com>
 */
#ifndef _TDINTERP_H_
#define _TDINTERP_H_

int RoussosInterp(float *u, int OutputWidth, int OutputHeight,
    const float *Input, int InputWidth, int InputHeight, double PsfSigma,
    double K, float Tol, int MaxMethodIter, int DiffIter);

#endif /* _TDINTERP_H_ */
