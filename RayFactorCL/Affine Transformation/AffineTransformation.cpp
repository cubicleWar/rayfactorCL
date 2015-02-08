//
//  AffineTransformation.cpp
//  RayFactor
//
//  Created by Trevor Walker on 26/08/11.
//  Copyright 2011 Native Dynamics. All rights reserved.
//

#include "AffineTransformation.h"

namespace affineTransformation {

    // Assumes that the vector (x, y, z) has already been normalised
    void rotate(Matrix4 &m, const float angle, const float &x, const float &y, const float &z)
    {
        Matrix4 rm;
        
        float ang = (float)angle*pi/180.0f;					//Convert the angle to radians
        float c = cos(ang);
        float s = sin(ang);
        float mc = 1.0f - c;
        
        // This is clockwise rotation
        rm.m[0] = c + mc*x*x;
        rm.m[1] = mc*x*y + s*z;
        rm.m[2] = mc*x*z - s*y;
        rm.m[4] = mc*y*x - s*z;
        rm.m[5] = c + mc*y*y;
        rm.m[6] = mc*y*z + s*x;
        rm.m[8] = mc*z*x + s*y;
        rm.m[9] = mc*z*y - s*x;
        rm.m[10] = c + mc*z*z;
        
        postMultiply(m, rm);
    }

    // Assumes that the vector (x, y, z) has already been normalised
    void invRotate(Matrix4 &m, const float angle, const float &x, const float &y, const float &z)
    {
        Matrix4 invRm;
        float ang = (float)angle*pi/180.0f;					//Convert the angle to radians
        float c = cos(ang);
        float s = sin(ang);
        float mc = 1.0f - c;
        
        // This is counter clockwise rotation
        invRm.m[0] = c + mc*x*x;
        invRm.m[1] = mc*x*y - s*z;
        invRm.m[2] = mc*x*z + s*y;
        invRm.m[4] = mc*y*x + s*z;
        invRm.m[5] = c + mc*y*y;
        invRm.m[6] = mc*y*z - s*x;
        invRm.m[8] = mc*z*x - s*y;
        invRm.m[9] = mc*z*y + s*x;
        invRm.m[10] = c + mc*z*z;
        
        preMultiply(m, invRm);
    }

    void scale(Matrix4 &m, const float sx, const float sy, const float sz)
    {
        Matrix4 scale;
        
        if(fabs(sx) < sEps || fabs(sy) < sEps || fabs(sz) < sEps) {
            fprintf(stderr, "Degenerate scaling transformation!\n");
        }
        
        scale.m[0] = sx;
        scale.m[5] = sy;
        scale.m[10] =  sz;
        
        postMultiply(m, scale);
    }

    void invScale(Matrix4 &m, const float sx, const float sy, const float sz)
    {
        Matrix4 invScale;
        
        if(fabs(sx) < sEps || fabs(sy) < sEps || fabs(sz) < sEps)
        {
            fprintf(stderr, "Degenerate scaling transformation!\n");
        }
        
        invScale.m[0] = 1.0f/sx;
        invScale.m[5] = 1.0f/sy;
        invScale.m[10] = 1.0f/sz;
        
        preMultiply(m, invScale);
    }

    void translate(Matrix4 &m, const float x, const float y, const float z)
    {
        Matrix4 tr;
        
        tr.m[3] = x;
        tr.m[7] = y;
        tr.m[11] = z;
        
        postMultiply(m, tr);
    }

    void invTranslate(Matrix4 &m, const float x, const float y, const float z)
    {
        Matrix4 invTr;
        
        invTr.m[3] = -x;
        invTr.m[7]= -y;
        invTr.m[11] = -z;
        
        preMultiply(m, invTr);
    }
    
}