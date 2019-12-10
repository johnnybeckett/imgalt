/*
imgalt - Image Alignment Tool
Author: GreatAttractor

version 0.5
2014/05/22

This code can be freely distributed and used for any purpose.

File description:
    Common types.

*/
#ifndef IMGALT_TYPES_H_
#define IMGALT_TYPES_H_

template <typename T>
struct strPoint
{
    T x, y;
    strPoint(T x, T y): x(x), y(y) { }
    strPoint(): x(0), y(0) { }

    strPoint & operator +=(const strPoint &p)
    {
        x += p.x; y += p.y;
        return *this;
    }
};

typedef strPoint<int> Point_t;
typedef strPoint<float> RealPoint_t;

typedef struct
{
    int x, y; /// coordinates of the xmin,ymin corner
    int width, height;
} Rectangle_t;

#endif /* IMGALT_TYPES_H_ */
