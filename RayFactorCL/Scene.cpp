#include "Scene.h"
#include <iostream>


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::Scene
//
//	Comments : Default Constructor
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	01/02/06	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Scene::Scene()
{
    
}

//Note put a deconstructor: see free Scene

#pragma mark readScene

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::readScene
//
//	Comments : The function for reading an xml NDRay project file and creating the scene
//
//	Arguments : filename is a c string that holds the file name of the xml project file
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	23/01/10	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void Scene::readScene(const char* filename)
{
    tinyxml2::XMLDocument doc;
    
    FILE* fp = fopen(filename, "r");
    if( !fp) {
        printf("Failed to load file \"%s\"\n", filename);
        exit(1);
    }
    fclose(fp);
    doc.LoadFile(filename);
    
    tinyxml2::XMLElement *settings = doc.FirstChildElement( "project" )->FirstChildElement( "settings" );
    tinyxml2::XMLElement *geometry = doc.FirstChildElement( "project" )->FirstChildElement( "geometry" );

    this->loadSettings( settings );
    this->loadGeometry( geometry );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::loadSettings
//
//	Comments : Parses the project settings from the <settings> xml block
//
//	Arguments : settings is a pointer to a TiXmlElement representing the <settings> block.
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	23/01/10	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Scene::loadSettings( tinyxml2::XMLElement *settings )
{
	tinyxml2::XMLElement *child = settings->FirstChild()->ToElement();
	
	//For each element in the <settings> block
	for( ; child; child=child->NextSiblingElement())
	{
		const char* typeStr = child->Name();
		if(strcmp(typeStr, kGlobalRayDensity) == 0) {
			globalRayDensity = child->FloatAttribute(kValue);
		} 
	}
	
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::loadGeometry
//
//	Comments : Is the higher level handler for the primitives stored in the <geometry> xml block. It 
//			   will determine what type of primitive is specified, create the primitive object, extract
//			   any unique attributes than pass it on for parsing of the generic primitive properties
//
//	Arguments : geometry is pointer to a TiXmlElement which contains the <geometry> block.
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	23/01/10	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// This method is cut down alot and probably needs some error checking
bool Scene::loadGeometry( tinyxml2::XMLElement *geometry )
{
	tinyxml2::XMLElement *primitive = geometry->FirstChild()->ToElement();
    
	//For each <primitive> block
	for( ; primitive; primitive=primitive->NextSiblingElement() )
	{
        Primitive *object = new Primitive;
		//string prType = "";
        
		bool willAnalyse = true;
		float boundMod = 1.f;
		const char* prType;
        
		if(primitive->Attribute( "type") != NULL) {
			prType = primitive->Attribute("type");
		}
		
		if(primitive->Attribute("isBounding") != NULL) {
			const char* isBoundingStr = primitive->Attribute("isBounding");
			boundMod = (strcmp(isBoundingStr, "true") == 0) ? -1.f : 1.f;
		}
		
		if(primitive->Attribute("analyse") != NULL) {
			const char* analyseStr = primitive->Attribute("analyse");
			willAnalyse = (strcmp(analyseStr, "true") == 0);
		}
		
		if(strcmp(prType, "cylinderSurface") == 0 || strcmp(prType, "taperedCylinderSurface") == 0 || strcmp(prType, "frustum") == 0) {
			
			float smallRadius;
			primitive->FirstChildElement( "smallRadius" )->QueryFloatAttribute( "value", &smallRadius );
            
            // Handle whether to use a faster primitive based on its type i.e cone, cylinder, frustum
            
			if(smallRadius == 1.0) {
                object->type = cylinder;
            } else {
                object->type = frustum;
                object->util1 = smallRadius;
            }
            object->util2 = boundMod;
		} else if(strcmp(prType, "annulus" ) == 0) {
			float smallRadius, largeRadius;
            object->type = annulus;
            
			primitive->FirstChildElement( "smallRadius" )->QueryFloatAttribute( "value", &smallRadius );
			primitive->FirstChildElement( "largeRadius" )->QueryFloatAttribute( "value", &largeRadius );
			
			//Assumes largeRaidus and smallRadius are defined and correct
            object->util1 = smallRadius*smallRadius;
            object->util2 = largeRadius*largeRadius;
		}
		else if(strcmp(prType, "rectangle" ) == 0) {
            object->type = rectangle;
		}
		else if(strcmp(prType, "disc" ) == 0) {
			object->type = disc;
		} else if(strcmp(prType, "sphere" ) == 0) {
			object->type = sphere;
            object->util2 = boundMod;
		} else {
            object->type = none;
        }
        
        if(object->type != none) {
            this->initPrimitive( primitive, *object);
			this->primitives.push_back(*object);
        }
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Scene::initPrimitive
//
//	Comments : Parsers all the attributes common to all primitives and applies them to the primitive 
//			   object
//
//	Arguments : xmlElement is a pointer to a TiXmlElement representing the <primitive> block.
//				primative is a reference to a ScnObject to which the parsed attributes will be applied.
//
//	Date		Developer		Action
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//	23/01/10	Trevor Walker	Created
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Scene::initPrimitive( tinyxml2::XMLElement *xmlElement, Primitive &object )
{
	if(!(xmlElement->FirstChild())) {
		return true;
	}
    
	tinyxml2::XMLElement *child = xmlElement->FirstChild()->ToElement();
	bool foundRayDensity = false;
    Matrix4 *trm = new Matrix4();
    Matrix4 *invtrm = new Matrix4();
    float rayDensity = 1;
    object.rayDensity = 1.0f;
    float xScale = 1.0f;
    float yScale = 1.0f;
    float zScale = 1.0f;
    
	// For each child in a <primitive> block
	for( ; child; child=child->NextSiblingElement())
	{
		const char* typeStr = child->Name();

		//Should be done in the order of scale, rotate, translate
		if(strcmp(typeStr, "rayDensity") == 0) 
		{
			rayDensity = child->FloatAttribute( "value");
			foundRayDensity = true;	
		}
		else if(strcmp(typeStr, "translation") == 0) 
		{
			float xTrans = child->FloatAttribute("x");
			float yTrans = child->FloatAttribute( "y");
			float zTrans = child->FloatAttribute( "z");
			
            affineTransformation::translate(*trm, xTrans, yTrans, zTrans);
            affineTransformation::invTranslate(*invtrm, xTrans, yTrans, zTrans);
		} 
		else if(strcmp(typeStr, "scale") == 0) 
		{
			xScale = child->FloatAttribute( "x");
			yScale = child->FloatAttribute( "y");
			zScale = child->FloatAttribute( "z");
            
            affineTransformation::scale(*trm, xScale, yScale, zScale);
            affineTransformation::invScale(*invtrm, xScale, yScale, zScale);
		} 
		else if(strcmp(typeStr, "rotate") == 0) 
		{
			float degRot = child->FloatAttribute( "degrees");
			float xAxis = child->FloatAttribute( "x");
			float yAxis = child->FloatAttribute( "y");
			float zAxis = child->FloatAttribute( "z");
			
			affineTransformation::rotate(*trm, degRot, xAxis, yAxis, zAxis);
            affineTransformation::invRotate(*invtrm, degRot, xAxis, yAxis, zAxis);
		}
	}
	
    // Copying only the first 12 elements as the last 4 elements are not used with the current transformations
    for(int i = 0; i < 12; i++) {
        object.m[i] = trm->m[i];
        object.invm[i] = invtrm->m[i];
    }

    // Set the ray density using the scaling and 
    if(object.type == cylinder) {
        object.rayDensity = (int)rayDensity*2.0*PI*zScale*pow(0.5*(pow(xScale,2)+pow(yScale, 2)),0.5); 
    } else if(object.type == frustum) {
        float sm = object.util1 + 1E-5f;            // Fix for when the small radius is 0
        object.rayDensity = (int)rayDensity*(sqrt((xScale*xScale+yScale*yScale)*sm*sm)*(1.f+sm)*zScale*PI)/(sqrt(2.f)*sm); 
    } else if(object.type == sphere) {
        float p = 1.6075; // Knud Thomens's forumla (max error +-1.061%)
        object.rayDensity = (int)rayDensity*(4.0f*PI*pow((pow(xScale,p)*pow(yScale,p)+pow(xScale,p)*pow(zScale,p)+pow(yScale,p)*pow(zScale,p))/3.0f,(1.0f/p)));
    } else if(object.type == annulus) {
        object.rayDensity = (int)(rayDensity*PI*(object.util2 - object.util1) + 0.5f);             // Doesn't account for additonal scalling NOTE : utils are already squared radius
    } else if(object.type == rectangle) {
        object.rayDensity = (int)rayDensity*2.0f*xScale*2.0f*yScale;
    } else if(object.type == disc) {
        object.rayDensity = (int)rayDensity*PI*xScale*yScale;       // Assumes it is scales uniformly
    }

	return true;
}
