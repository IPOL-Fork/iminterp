/**
 * @file invmat.c
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

#include <stdio.h>
#include <ipol/basic.h>

/**
 * @brief Invert matrix through QR decomposition
 *
 * @param InverseData pointer to memory for holding the result
 * @param AData pointer to column-major matrix data
 * @param N the number of dimensions
 *
 * @return 1 on success, 0 on failure
 *
 * The input data is overwritten during the computation.  \c InverseData
 * should be allocated before calling this function with space for at least
 * N^2 doubles.  Matrices are represented in column-major format, meaning
 *    A(i,j) = AData[i + N*j], 0 <= i, j < N.
 */
int InvertMatrix(double *InverseData, double *AData, int N)
{
    double *c = 0, *d = 0, *Aj, *Ak, *Inversek;
    double Temp, Scale, Sum;
    int i, j, k, Success = 0;
    

    if(!(c = (double *)Malloc(sizeof(double)*N))
        || !(d = (double *)Malloc(sizeof(double)*N)))
        goto Catch;
    
    for(k = 0, Ak = AData; k < N - 1; k++, Ak += N)
    {
        Scale = 0.0;
        
        for(i = k; i < N; i++)
            if((Temp = fabs(Ak[i])) > Scale)
                Scale = Temp;
        
        if(Scale == 0.0)
        {
            ErrorMessage("Singular matrix.\n");
            goto Catch; /* Singular matrix */
        }
        
        for(i = k; i < N; i++)
            Ak[i] /= Scale;
        
        for(Sum = 0.0, i = k; i < N; i++)
            Sum += Ak[i]*Ak[i];
        
        Temp = (Ak[k] >= 0.0)? sqrt(Sum) : -sqrt(Sum);
        Ak[k] += Temp;
        c[k] = Temp*Ak[k];
        d[k] = -Scale*Temp;
        
        for(j = k + 1, Aj = Ak + N; j < N; j++, Aj += N)
        {
            for(Scale = 0.0, i = k; i < N; i++)
                Scale += Ak[i] * Aj[i];
                
            Scale /= c[k];
            
            for(i = k; i < N; i++)
                Aj[i] -= Scale*Ak[i];
        }
    }

    d[N-1] = Ak[k];
    
    if(d[N-1] == 0.0)
    {
        ErrorMessage("Singular matrix.\n");
        goto Catch; /* Singular matrix */
    }
    
    for(k = 0, Inversek = InverseData; k < N; k++, Inversek += N)
    {
        for(i = 0; i < N; i++)
            Inversek[i] = -AData[k]*AData[i]/c[0];
            
        Inversek[k] += 1.0;
            
        for(j = 1, Aj = AData + N; j < N-1; j++, Aj += N)
        {
            for(Scale = 0.0, i = j; i < N; i++)
                Scale += Aj[i]*Inversek[i];
            
            Scale /= c[j];
            
            for(i = j; i < N; i++)
                Inversek[i] -= Scale*Aj[i];
        }
        
        Inversek[j] /= d[N-1];
        
        for(i = N - 2; i >= 0; i--)
        {
            for(Sum = 0.0, j = i + 1, Aj = AData + N*(i + 1); j < N; j++, Aj += N)
                Sum += Aj[i]*Inversek[j];
            
            Inversek[i] = (Inversek[i] - Sum)/d[i];
        }
    }
    
    Success = 1; /* Finished successfully */
Catch: /* Clean up */
    Free(d);
    Free(c);
    return Success;
}
