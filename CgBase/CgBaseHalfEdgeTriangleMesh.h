#ifndef CGBASE_HALFEDGE_TRIANGLEMESH_H
#define CGBASE_HALFEDGE_TRIANGLEMESH_H


#include <vector>
#include <glm/glm.hpp>
#include "CgBase/CgBaseRenderableObject.h"
#include "CgBase/CgBaseHalfdgePrimitives.h"





class CgBaseHalfEdgeTriangleMesh : public CgBaseRenderableObject
{
public:
    virtual ~CgBaseHalfEdgeTriangleMesh(){};


    virtual const std::vector<CgBaseHeFace*>& getFaces()     const =0;

};

#endif // CGBASE_HALFEDGE_TRIANGLEMESH_H
