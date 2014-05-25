/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	Object : Scene
//
//	Description : The Scene object stores all the data/objects for the system to be analysed
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __SCENE_H__
#define __SCENE_H__

#include <vector>
#include <cstdlib>
#include <cstring>

#include "RayfactorConstants.h"
#include "tinyxml2.h"

#include "Matrix4.h"
#include "AffineTransformation.h"


class Scene {
    private:
		
    public:
        float globalRayDensity;
        std::vector<Primitive> primitives;
        Scene();
		// Methods for reading a scene saved in an xml document
		void readScene(const char* filename);
		bool loadSettings(tinyxml2::XMLElement *settings);
		bool loadGeometry(tinyxml2::XMLElement *geometry);
		bool initPrimitive(tinyxml2::XMLElement *xmlElement, Primitive &primitive);
};

#endif
