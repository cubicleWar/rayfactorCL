/*
 *  RayfactorConstants.h
 *  RayFactor 0.5
 *
 *  Created by Trevor Walker on 20/02/10.
 *  Copyright 2010 University Of Sydney. All rights reserved.
 *
 */

#include <limits>

#define PI 3.141592f       // for float use 7 sig fig, for double 16

//when working do the same for maximum doubles etc
const double inf = std::numeric_limits<float>::infinity();

#define kDefaultRayDensity 100000

// Parsing constants used when parsing the input xml file in Scene
#define kGlobalRayDensity	"globalRayDensity"
#define kDescription		"description"
#define kValue				"value"

#define noOfPrimitveTypes   6
enum PrimitiveType {
    cylinder,
    frustum,
    sphere,
    annulus,
    rectangle,
    disc,
    none
};

//need to make sure this is aligned
// Could make this a 4 x 3 (rows) and save 4 * 4 = 16 bytes - Done check for implications
// Current size is 156 bytes (byte = size of a char = 8 bits (but could be different depending on the platform)
typedef struct Primitive {
    float m[12];
    float invm[12];
    enum PrimitiveType type;
    int rayDensity;
    float util1;                // Used to store utility data (Tapered cylinder = top radius; Annulus = small radius^2)
    float util2;                // Used to store utility data (Annulus = large radius^2, Bounding modification for 3D shapes)
} Primitive;
