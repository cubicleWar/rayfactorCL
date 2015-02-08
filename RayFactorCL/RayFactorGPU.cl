#include <Random123/philox.h>

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable

#define R123_0x1p_24f (1.f/16777216.f)
#define EPS 1E-3f           // Cannot be smaller than 1E-3 due to use of native_sqrt in hit routines
__constant float M_2PI = 2.0f*M_PI_F;


#pragma mark -
#pragma mark Typedefs
typedef struct Vector {
    float x;
    float y;
    float z;
} Vector;


typedef struct Ray {
    Vector s;
    Vector c;
} Ray;

enum PrimitiveType {
    cylinder,
    frustum,
    sphere,
    annulus,
    rectangle,
    disc,
    none
};

// need to check this is aligned
// instead use util2 to hold boundMod and get rid of BoundMod and id to get size of 112 bytes
typedef struct Primitive {
    float m[12];
    float invm[12];
    enum PrimitiveType type;
    int rayDensity;
    float util1;                // Used to store utility data (Tapered cylinder = top radius; Annulus = small radius)
    float util2;                // Used to store utility data (Annulus = large radius, Bounding modification for 3D shapes)
} Primitive;

typedef struct BoundingVolume {
    float invm[12];
    enum PrimitiveType type;
    int objStartIndex;
    int objEndIndex;
    float util1;                // Used to store utility data (Tapered cylinder = top radius; Annulus = small radius^2)
    float util2;                // Used to store utility data (Annulus = large radius^2, Bounding modification for 3D shapes)
} BoundingVolume;

#pragma mark -
#pragma mark Utility methods

inline float4 u01_closed_open_4x32_24(uint4 i)
{
    return convert_float4(i >> 8)* R123_0x1p_24f;
}

inline void transformPoint(global const float *m, const Vector *osPoint, Vector *wsPoint)
{
    float4 osPointf4 = (float4)(osPoint->x,osPoint->y,osPoint->z,1.f);
    wsPoint->x = dot((float4)(m[0],m[1],m[2],m[3]), osPointf4);
    wsPoint->y = dot((float4)(m[4],m[5],m[6],m[7]), osPointf4);
    wsPoint->z = dot((float4)(m[8],m[9],m[10],m[11]), osPointf4);
    
}

inline void transformNormal(global const float *invm, const Vector *osNormal, Vector *wsNormal)
{
    wsNormal->x = mad(invm[0],osNormal->x, mad(invm[4],osNormal->y, invm[8]*osNormal->z));
    wsNormal->y = mad(invm[1],osNormal->x, mad(invm[5],osNormal->y, invm[9]*osNormal->z));
    wsNormal->z = mad(invm[2],osNormal->x, mad(invm[6],osNormal->y, invm[10]*osNormal->z));
    
    
    const float t = native_rsqrt(mad(wsNormal->x,wsNormal->x, mad(wsNormal->y,wsNormal->y, wsNormal->z*wsNormal->z)));
    
    wsNormal->x *= t;
    wsNormal->y *= t;
    wsNormal->z *= t;
}

inline void alignRayDirection(Vector *wsDir, const Vector *osDir, const Vector *wsNormal)
{
    float d = native_recip(1.0f + wsNormal->z);
    
    if(isinf(d)) d = 0.5f;
    
    wsDir->x = wsNormal->z*osDir->x + wsNormal->y*d*(wsNormal->y*osDir->x - wsNormal->x*osDir->y) + wsNormal->x*osDir->z;
    wsDir->y = wsNormal->x*d*(wsNormal->x*osDir->y - wsNormal->y*osDir->x) + wsNormal->z*osDir->y + wsNormal->y*osDir->z;
    wsDir->z = -wsNormal->x*osDir->x - wsNormal->y*osDir->y + wsNormal->z*osDir->z;
    
}

inline void xfrmRay(global const float *invm, Ray *osRay, const Ray *wsRay )
{
    osRay->s.x = mad(invm[0],wsRay->s.x,mad(invm[1],wsRay->s.y,mad(invm[2],wsRay->s.z,invm[3])));
    osRay->c.x = mad(invm[0],wsRay->c.x,mad(invm[1],wsRay->c.y,invm[2]*wsRay->c.z));
    
    osRay->s.y = mad(invm[4],wsRay->s.x,mad(invm[5],wsRay->s.y,mad(invm[6],wsRay->s.z,invm[7])));
    osRay->c.y = mad(invm[4],wsRay->c.x,mad(invm[5],wsRay->c.y,invm[6]*wsRay->c.z));
    
    osRay->s.z = mad(invm[8],wsRay->s.x,mad(invm[9],wsRay->s.y,mad(invm[10],wsRay->s.z,invm[11])));
    osRay->c.z = mad(invm[8],wsRay->c.x,mad(invm[9],wsRay->c.y,invm[10]*wsRay->c.z));
}

void inline osRayDirection(const float4 *rands, Vector *osVec)
{
    osVec->y = sincos((*rands).z*M_2PI, &osVec->x);
    float temp1 = (*rands).w;
    osVec->z = native_sqrt(1.0f-temp1);
    
    temp1 = native_sqrt(temp1);
    
    osVec->x = temp1*osVec->x;
    osVec->y = temp1*osVec->y;
}

#pragma mark -
#pragma mark Object intersection methods

inline float rectangleHit(global const float *invm, const Ray *wsRay)
{
    Ray osRay;
    
    xfrmRay(invm, &osRay, wsRay);
    
    const float th = -osRay.s.z*native_recip(osRay.c.z);
    
    const float hx = fabs(osRay.s.x + osRay.c.x*th);
    const float hy = fabs(osRay.s.y + osRay.c.y*th);
    
    if(hx > 1.f || hy > 1.f || th < EPS)
    {
        return INFINITY;
    }
    else
    {
        return th;
    }
}

inline float discHit(global const float *invm, const Ray *wsRay)
{
    Ray osRay;
    
    xfrmRay(invm, &osRay, wsRay);
    
    const float th = -osRay.s.z * native_recip(osRay.c.z);
    const float hx = osRay.s.x + osRay.c.x*th;
    const float hy = osRay.s.y + osRay.c.y*th;
    const float exp = hx * hx + hy * hy;
    
    if(exp <= 1.f && th > EPS)
    {
        return th;
    }
    else
    {
        return INFINITY;
    }
}

inline float annulusHit(global const float *invm, const Ray *wsRay, const float iDiameter, const float oDiameter)
{
    Ray osRay;
    
    xfrmRay(invm, &osRay, wsRay);
    
    const float th = -osRay.s.z * native_recip(osRay.c.z);
    const float hx = osRay.s.x + osRay.c.x*th;
    const float hy = osRay.s.y + osRay.c.y*th;
    const float exp = mad(hx, hx, hy * hy);
    
    if(exp >= iDiameter && exp <= oDiameter && th > EPS)
    {
        return th;
    }
    else
    {
        return INFINITY;
    }
}

inline float cylinderHit(global const float *invm, const Ray *wsRay)
{
    Ray osRay;
    
    xfrmRay(invm, &osRay, wsRay);
    
    const float A = mad(osRay.c.x, osRay.c.x, osRay.c.y * osRay.c.y);
    const float B = -mad(osRay.s.x, osRay.c.x, osRay.s.y * osRay.c.y);
    const float C = mad(osRay.s.x, osRay.s.x, mad(osRay.s.y, osRay.s.y, -1.f));
    
    const float discrim = native_sqrt(B*B - A*C);
    
    float t1, t2;
    t1 = t2 = native_recip(A);
    t1 *= (B - discrim);
    t2 *= (B + discrim);
    
    const float zHit1 = osRay.s.z + osRay.c.z*t1;
    const float zHit2 = osRay.s.z + osRay.c.z*t2;
    
    if(t1 < EPS || zHit1 < 0.f || zHit1 > 1.f)
    {
        t1 = INFINITY;
    }
    
    if(t2 < t1 && zHit2 > 0.f && zHit2 < 1.0f && t2 > EPS)
    {
        return t2;
    }
    else
    {
        return t1;
    }
}

inline float frustumHit(global const float *invm, const Ray *wsRay, const float tDiameter)
{
    Ray osRay;
    
    xfrmRay(invm, &osRay, wsRay);
    
    const float fDir = (tDiameter - 1.0f)*osRay.c.z;
    const float fStart = (tDiameter - 1.0f)*osRay.s.z + 1.0f;
    
    const float A = mad(osRay.c.x, osRay.c.x, mad(osRay.c.y, osRay.c.y, -fDir*fDir));
    const float B = -mad(osRay.s.x, osRay.c.x, mad(osRay.s.y, osRay.c.y, -fStart*fDir));
    const float C = mad(osRay.s.x, osRay.s.x, mad(osRay.s.y, osRay.s.y, -fStart*fStart));
    
    const float discrim = native_sqrt(B*B - A*C);     // This was native_sqrt but it introduces alot of error in the intersection check
    
    float t1, t2;
    t1 = t2 = native_recip(A);
    t1 *= (B - discrim);
    t2 *= (B + discrim);
    
    const float zHit1 = osRay.s.z + osRay.c.z * t1;
    const float zHit2 = osRay.s.z + osRay.c.z * t2;
    
    if(t1 < EPS || zHit1 < 0.f || zHit1 > 1.f)
    {
        t1 = INFINITY;
    }
    
    if(t2 < t1 && zHit2 > 0.f && zHit2 < 1.0f && t2 > EPS)
    {
        return t2;
    }
    else
    {
        return t1;
    }
}

inline float sphereHit(global const float *invm, const Ray *wsRay)
{
    Ray osRay;
    
    xfrmRay(invm, &osRay, wsRay);
    
    const float A = osRay.c.x*osRay.c.x + osRay.c.y*osRay.c.y + osRay.c.z*osRay.c.z;
    const float B = osRay.c.x*osRay.s.x + osRay.c.y*osRay.s.y + osRay.c.z*osRay.s.z;
    const float C = osRay.s.x*osRay.s.x + osRay.s.y*osRay.s.y + osRay.s.z*osRay.s.z - 1.0f;
    
    float discrim = B*B - A*C;
    discrim = select(discrim, INFINITY, isless(discrim, 0.0f));
    
    if(discrim > 0.f)
    {
        float t1 = (-B - discrim)*native_recip(A);
        float t2 = (-B + discrim)*native_recip(A);
        
        if(t2 < t1 && t2 > EPS)
        {
            return t2;
        }
        else if(t1 > EPS)
        {
            return t1;
        }
    }
    
    return INFINITY;
}

#pragma mark -
#pragma mark Ray formulation methods


// better to move the ray shooting loop into this function
inline void discRay(const float4 *rands, Ray *ray, global const Primitive *p)
{
    Vector osVec, osNormal, wsNormal;
    
    float temp1 = (*rands).x * M_2PI;
    
    osVec.x = sincos(temp1, &osVec.y);
    
    temp1 = native_sqrt((*rands).y);
    
    osVec.x = temp1 * osVec.x;
    osVec.y = temp1 * osVec.y;
    osVec.z = 0.f;
    
    transformPoint(&p->m[0], &osVec, &ray->s);
    
    osRayDirection(rands, &osVec);
    // Get the world space Normal
    osNormal.x = osNormal.y = 0.f;
    osNormal.z = 1.f;
    
    transformNormal(&p->invm[0], &osNormal, &wsNormal);
    alignRayDirection(&ray->c, &osVec, &wsNormal);
}

inline void annulusRay(const float4 *rands, Ray *ray, global const struct Primitive *p)
{
    Vector osVec, osNormal, wsNormal;
    
    osVec.x = sincos((*rands).x * M_2PI, &osVec.y);
    
    float temp1 = (p->util2 - p->util1) * (*rands).y + p->util1;
    temp1 = native_sqrt(temp1);
    
    osVec.x = temp1 * osVec.x;
    osVec.y = temp1 * osVec.y;
    osVec.z = 0.f;
    
    transformPoint(&p->m[0], &osVec, &ray->s);
    
    osRayDirection(rands, &osVec);
    // Get the world space Normal
    osNormal.x = osNormal.y = 0.f;
    osNormal.z = 1.f;
    
    transformNormal(&p->invm[0], &osNormal, &wsNormal);
    alignRayDirection(&ray->c, &osVec, &wsNormal);
}

inline void cylinderRay(const float4 *rands, Ray *ray, global const struct Primitive *p)
{
    Vector osVec, osNormal, wsNormal;
    
    osVec.y = sincos((*rands).x * M_2PI, &osVec.x);
    osVec.z = (*rands).y;
    transformPoint(&p->m[0], &osVec, &ray->s);
    
    // Normal is based on the os surface points
    osNormal.x = osVec.x * p->util2; // Make isBounding a float and just times by it
    osNormal.y = osVec.y * p->util2;
    osNormal.z = 0.f;
    
    //wsNormal = osNormal;
    transformNormal(&p->invm[0], &osNormal, &wsNormal);
    
    osRayDirection(rands, &osVec);
    alignRayDirection(&ray->c, &osVec, &wsNormal);
}



inline void frustumRay(const float4 *rands, Ray *ray, global const struct Primitive *p)
{
    Vector osVec, osNormal, wsNormal;
    
    osVec.z = (native_sqrt(1.0f - (*rands).x + (*rands).x*p->util1*p->util1) - 1.0f) * native_recip(p->util1 - 1.0f);
    
    float temp1 = 1.0f + osVec.z*(p->util1 - 1.0f);        // the radius at the z cordinate
    
    osVec.y = sincos(M_2PI*(*rands).y, &osVec.x);
    
    osVec.x *= temp1;
    osVec.y *= temp1;
    
    transformPoint(&p->m[0], &osVec, &ray->s);
    
    // Normal is based on the os surface points
    osNormal.x = osVec.x * p->util2; // Make isBounding a float and just times by it
    osNormal.y = osVec.y * p->util2;
    osNormal.z = (1.f - p->util1)*(1.f + (p->util1 - 1.f)*osVec.z) * p->util2;
    
    transformNormal(&p->invm[0], &osNormal, &wsNormal);
    
    osRayDirection(rands, &osVec);
    alignRayDirection(&ray->c, &osVec, &wsNormal);
}



// Here osVec is used as both osPoint and osRay to save on memory of have two os
// Maybe pass everything to this method
inline void rectangleRay(const float4 *rands, Ray *ray, global const struct Primitive *p)
{
    Vector osVec, osNormal, wsNormal;
    
    osVec.x = 2.0f*(*rands).x - 1.0f;
    osVec.y = 2.0f*(*rands).y - 1.0f;
    osVec.z = 0.f;
    
    // Transform the point
    transformPoint(&p->m[0], &osVec, &ray->s);
    osRayDirection(rands, &osVec);
    
    // Get the world space Normal
    osNormal.x = osNormal.y = 0.f;
    osNormal.z = 1.f;
    
    transformNormal(&p->invm[0], &osNormal, &wsNormal);
    alignRayDirection(&ray->c, &osVec, &wsNormal);
}

inline void sphereRay(const float4 *rands, Ray *ray, global const struct Primitive *p)
{
    Vector osVec, osNormal, wsNormal;
    
    osVec.z = 2.0f*(*rands).x - 1.0f;
    
    float temp1 = native_sqrt(1.0f - osVec.z*osVec.z);
    
    osVec.y = sincos(M_2PI*(*rands).y, &osVec.x);
    
    osVec.x *= temp1;
    osVec.y *= temp1;
    
    transformPoint(&p->m[0], &osVec, &ray->s);
    
    osNormal.x = osVec.x * p->util2; // Make isBounding a float and just times by it
    osNormal.y = osVec.y * p->util2;
    osNormal.z = osVec.z * p->util2;
    
    transformNormal(&p->invm[0], &osNormal, &wsNormal);
    
    osRayDirection(rands, &osVec);
    alignRayDirection(&ray->c, &osVec, &wsNormal);
    
}

inline int findObject(const int startIndex, const int endIndex, __global const Primitive *objects, const Ray *ray)
{
    int objectNo = -1;
    float tbest = INFINITY, t;
    
    for(int w = startIndex; w < endIndex; w++)
    {
        switch(objects[w].type)
        {
            case rectangle  : t = rectangleHit(objects[w].invm, ray); break;
            case disc       : t = discHit(objects[w].invm, ray); break;
            case annulus    : t = annulusHit(objects[w].invm, ray, objects[w].util1, objects[w].util2); break;
            case cylinder   : t = cylinderHit(objects[w].invm, ray); break;
            case sphere     : t = sphereHit(objects[w].invm, ray); break;
            case frustum    : t = frustumHit(objects[w].invm, ray, objects[w].util1); break;
            default: break;
        }
        
        if(t < tbest)
        {
            tbest = t;
            objectNo = w;
        }
        
    }
    
    return objectNo;
}

#pragma mark -
#pragma mark OpenCL Kernels

// Need to put primitives into constant memory - good for broadcasts constant memory is around 16 - 64 KB's so 136 - 546 at 120 bytes a primitive
// Could load the primitives into local memory rather than private memory

// Version which writes view factor results diretly to global memory (appears 35% faster than using local memory on smaller geometries)
__kernel void runGPU(int fromIndex, const int n_objects, global const Primitive *objects, __local int *x_res, __global int *hits)
{
    int gid = get_global_id(0);

    
    philox4x32_key_t k = {{gid, 0xdecafbad}};
    philox4x32_ctr_t c = {{0, 0xf00dcafe, 0xdeadbeef, 0xbeeff00d}};
    
    float4 rands;
    int objectNo;
    Ray ray;

    union {
        philox4x32_ctr_t c;
        uint4 i;
    } u;
    c.v[0]++;
    
    u.c = philox4x32(c, k);   // Appears random number generation is around 30% of the run time
    rands = u01_closed_open_4x32_24(u.i);
    
    for (int v = 0; v < 4; v++)
    {
        objectNo = -1;
        //tbest = INFINITY;
        rands = rands.yzwx; //rotate the random numbers
        
        switch(objects[fromIndex].type)
        {
            case rectangle  : rectangleRay(&rands, &ray, &objects[fromIndex]); break;
            case disc       : discRay(&rands, &ray, &objects[fromIndex]); break;
            case annulus    : annulusRay(&rands, &ray, &objects[fromIndex]); break;
            case cylinder   : cylinderRay(&rands, &ray, &objects[fromIndex]); break;
            case sphere     : sphereRay(&rands, &ray, &objects[fromIndex]); break;
            case frustum    : frustumRay(&rands, &ray, &objects[fromIndex]); break;
            default: break;
        }
        
        objectNo = findObject(0, n_objects, objects, &ray);
        
        
        if(objectNo != -1)
        {
            atomic_inc(&hits[fromIndex*n_objects + objectNo]);
        }
    } // end for sub rays
}

// This kernel uses a simple bounding volume method instead
__kernel void runBVGPU(const int fromElementIndex,
                       const int n_objects,
                       global const Primitive *objects,
                       __local int *x_res,
                       __global int *hits,
                       const int n_bvs,
                       global const BoundingVolume *bvs)
{
    float4 rands;
    float tbest, t;
    int bvIndex;
    Ray ray;
    
    // Set the PRNG key and seed
    int gid = get_global_id(0);
    philox4x32_key_t k = {{gid, 0xdecafbad}};
    philox4x32_ctr_t c = {{0, 0xf00dcafe, 0xdeadbeef, 0xbeeff00d}};
    
    // Mask the PRN to 4 x 32 bit floats
    union {
        philox4x32_ctr_t c;
        uint4 i;
    } u;
    c.v[0]++;
    
    // Generate the PRN's
    u.c = philox4x32(c, k);   // Appears random number generation is around 30% of the run time
    rands = u01_closed_open_4x32_24(u.i);
    
    for (int v = 0; v < 4; v++)
    {
        bvIndex = -1;
        tbest = INFINITY;
        rands = rands.yzwx; //rotate the random numbers
        
        switch(objects[fromElementIndex].type)
        {
            case rectangle  : rectangleRay(&rands, &ray, &objects[fromElementIndex]); break;
            case disc       : discRay(&rands, &ray, &objects[fromElementIndex]); break;
            case annulus    : annulusRay(&rands, &ray, &objects[fromElementIndex]); break;
            case cylinder   : cylinderRay(&rands, &ray, &objects[fromElementIndex]); break;
            case sphere     : sphereRay(&rands, &ray, &objects[fromElementIndex]); break;
            case frustum    : frustumRay(&rands, &ray, &objects[fromElementIndex]); break;
            default: break;
        }
        
        
        
        // Test for intersection against the bounding volumes
        for (int w = 0; w < n_bvs; w++)
        {
            switch(bvs[w].type)
            {
                // Check intersection with everything in the same volume,
                //
                case rectangle  : t = rectangleHit(bvs[w].invm, &ray); break;
                case disc       : t = discHit(bvs[w].invm, &ray); break;
                case annulus    : t = annulusHit(bvs[w].invm, &ray, objects[w].util1, objects[w].util2); break;
                case cylinder   : t = cylinderHit(bvs[w].invm, &ray); break;
                case sphere     : t = sphereHit(bvs[w].invm, &ray); break;
                case frustum    : t = frustumHit(bvs[w].invm, &ray, objects[w].util1); break;
                default: break;
            }
            
            
            
            if(t < tbest)
            {
                tbest = t;
                bvIndex = w;
            }
        }
        
        if (bvIndex != -1)
        {
            int objectIndex = findObject(bvs[bvIndex].objStartIndex, bvs[bvIndex].objEndIndex, objects, &ray);
          
            if (objectIndex > -1)
            {
                atomic_inc(&hits[fromElementIndex*n_objects + objectIndex]);
            }
        }
    } // end for sub rays
}
