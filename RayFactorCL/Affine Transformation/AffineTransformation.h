//
//  AffineTransformation.h
//  RayFactor
//
//  Created by Trevor Walker on 26/08/11.
//  Copyright 2011 Native Dynamics. All rights reserved.
//
#include "Matrices.h"
#include "Matrix4.h"
#include <cmath>
#include <stdio.h>

#ifndef RayFactor_AffineTransformation_h
#define RayFactor_AffineTransformation_h

namespace affineTransformation {

    static const float sEps = 0.00001f;
    static const float pi = 3.141592653589793f;
        
    void rotate(Matrix4 &m, const float angle, const float &x, const float &y, const float &z);
    void invRotate(Matrix4 &m, const float angle, const float &x, const float &y, const float &z);
    void scale(Matrix4 &m, const float sx, const float sy, const float sz);
    void invScale(Matrix4 &m, const float sx, const float sy, const float sz);
    void translate(Matrix4 &m, const float x, const float y, const float z);
    void invTranslate(Matrix4 &m, const float x, const float y, const float z);

}
#endif
