/**
 * @file finterp.h
 * @brief Fourier zero-padding interpolation
 * @author Pascal Getreuer <getreuer@gmail.com>
 */
#ifndef _FINTERP_H_
#define _FINTERP_H_

typedef enum
{
	BOUNDARY_CONSTANT = 0,
	BOUNDARY_HSYMMETRIC = 1,
	BOUNDARY_WSYMMETRIC = 2
} boundaryhandling;

int FourierScale2d(float *Dest, int DestWidth, float XStart,
    int DestHeight, float YStart,
    const float *Src, int SrcWidth, int SrcHeight, int NumChannels,
    double PsfSigma, boundaryhandling Boundary);

#endif /* _FINTERP_H_ */
