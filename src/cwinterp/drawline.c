/**
 * @file drawline.c
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

#include <math.h>
#include "drawline.h"

#define PIXEL_STRIDE	3

#define SWAP(A,B)	\
    Temp = A; \
    A = B; \
    B = Temp;

    
static float fpart(float x)
{
    return (x - (float)floor(x));
}


static void PlotPixel(uint32_t *Image, int Width, int Height, int x, int y,
    float Alpha, const float *Color)
{
    if(0 <= x && x < Width && 0 <= y && y < Height)
    {
        const float CAlpha = 1 - Alpha;
        
        Image += x + Width*y;
        ((uint8_t *)Image)[0] = (uint8_t)(
            CAlpha*((uint8_t *)Image)[0] + Alpha*Color[0] + 0.5f);
        ((uint8_t *)Image)[1] = (uint8_t)(
            CAlpha*((uint8_t *)Image)[1] + Alpha*Color[1] + 0.5f);
        ((uint8_t *)Image)[2] = (uint8_t)(
            CAlpha*((uint8_t *)Image)[2] + Alpha*Color[2] + 0.5f);
    }
}


/**
 * @brief Draw an antialiased line
 *
 * @param Image image to draw on
 * @param Width, Height image dimensions
 * @param x1, y1, x2, y2 line coordinates
 * @param Color RGB color of the line
 */
void DrawLine(uint32_t *Image, int Width, int Height,
    float x1, float y1, float x2, float y2, const float *Color)
{
    float yend, Gap, dx, dy, Gradient, y, Temp;
    int xend, ix, iy, ix1, iy1, ix2, iy2, Swapped = 0;
       
    
    dx = x2 - x1;
    dy = y2 - y1;
    
    if(fabs(dx) < fabs(dy))
    {
        SWAP(x1,y1)
        SWAP(x2,y2)
        SWAP(dx,dy)
        Swapped = 1;
    }
        
    if(x2 < x1)
    {
        SWAP(x1,x2)
        SWAP(y1,y2)
    }
    
    if(dx == 0)
        return;
    
    Gradient = dy/dx;
    
    xend = (int)floor(x1 + 0.5f);
    yend = y1 + Gradient * (xend - x1);
    y = yend + Gradient;
    Gap = 1 - fpart(x1 + 0.5f);
    ix1 = xend;
    iy1 = (int)floor(yend);
    
    if(!Swapped)
    {
        PlotPixel(Image, Width, Height, ix1, iy1, (1 - fpart(yend)) * Gap, Color);
        PlotPixel(Image, Width, Height, ix1, iy1 + 1, fpart(yend) * Gap, Color);
    }
    else
    {
        PlotPixel(Image, Width, Height, iy1, ix1, (1 - fpart(yend)) * Gap, Color);
        PlotPixel(Image, Width, Height, iy1 + 1, ix1, fpart(yend) * Gap, Color);
    }
        
    xend = (int)floor(x2 + 0.5f);
    yend = y2 + Gradient * (xend - x2);
    Gap = fpart(x2 + 0.5f);
    ix2 = xend;
    iy2 = (int)floor(yend);
    
    if(!Swapped)
    {
        PlotPixel(Image, Width, Height, ix2, iy2, (1 - fpart(yend)) * Gap, Color);
        PlotPixel(Image, Width, Height, ix2, iy2 + 1, fpart(yend) * Gap, Color);
        
        for(ix = ix1 + 1; ix < ix2; ix++)
        {
            iy = (int)y;
            PlotPixel(Image, Width, Height, ix, iy, 1 - fpart(y), Color);
            PlotPixel(Image, Width, Height, ix, iy + 1, fpart(y), Color);
            y += Gradient;
        }
    }
    else
    {
        PlotPixel(Image, Width, Height, iy2, ix2, (1 - fpart(yend)) * Gap, Color);
        PlotPixel(Image, Width, Height, iy2 + 1, ix2, fpart(yend) * Gap, Color);
        
        for(ix = ix1 + 1; ix < ix2; ix++)
        {
            iy = (int)y;
            PlotPixel(Image, Width, Height, iy, ix, 1 - fpart(y), Color);
            PlotPixel(Image, Width, Height, iy + 1, ix, fpart(y), Color);
            y += Gradient;
        }
    }
}
