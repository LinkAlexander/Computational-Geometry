#ifndef CGHALFEDGETRIANGLEMESH_H
#define CGHALFEDGETRIANGLEMESH_H

#include "CgBase/CgBaseHalfEdgeTriangleMesh.h"
#include "CgHalfEdgePrimitives.h"

#include <vector>
#include <glm/glm.hpp>
#include <string>

class CgHalfEdgeTriangleMesh : public CgBaseHalfEdgeTriangleMesh
{
public:
    CgHalfEdgeTriangleMesh();
    CgHalfEdgeTriangleMesh(int id);

    ~CgHalfEdgeTriangleMesh();

    //inherited from CgBaseRenderableObject
    Cg::ObjectType getType() const;
    unsigned int getID() const;

    //inherited from CgBaseHalfEdgeTriangleMesh

    const std::vector<CgBaseHeFace*>& getFaces() const;

    //own stuff

    void init( std::string filename);
    const glm::vec3 getCenter() const;


private:

    std::vector<CgBaseHeFace*> m_faces;
    std::vector<CgBaseHeVert*> m_verts;
    std::vector<CgBaseHeEdge*> m_edges;


    const Cg::ObjectType m_type;
    const unsigned int m_id;

};


inline Cg::ObjectType  CgHalfEdgeTriangleMesh::getType() const {return m_type;}
inline unsigned int CgHalfEdgeTriangleMesh::getID() const {return m_id;}



#endif // CGHALFEDGETRIANGLEMESH_H
