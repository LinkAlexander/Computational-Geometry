#ifndef CGHALFEDGETRIANGLEMESH_H
#define CGHALFEDGETRIANGLEMESH_H

#include "CgBase/CgBaseHalfEdgeTriangleMesh.h"
#include "CgHalfEdgePrimitives.h"

#include <glm/glm.hpp>
#include <string>
#include <vector>

class CgHalfEdgeTriangleMesh : public CgBaseHalfEdgeTriangleMesh {
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

    void init(std::string filename);
    const glm::vec3 getCenter() const;
    void loopSubdivision();

    void applyPickRay(glm::vec3 pickRayStart, glm::vec3 pickRayDirection);

private:
    std::vector<CgBaseHeFace*> m_faces;
    std::vector<CgBaseHeVert*> m_verts;
    std::vector<CgBaseHeEdge*> m_edges;

    const Cg::ObjectType m_type;
    const unsigned int m_id;

    void loadDemoTriangles();
    void loadFromVertexList(std::vector<glm::vec3> temp_vertices, std::vector<unsigned int> temp_indices);

    float calculateLoopSubdivisionBeta(int n);

    size_t getClosestPointToRay(glm::vec3 rayStart, glm::vec3 rayDirection);

    float distanceRayToPoint(glm::vec3 rayStart, glm::vec3 rayDirection, glm::vec3 point);

};

inline Cg::ObjectType CgHalfEdgeTriangleMesh::getType() const { return m_type; }
inline unsigned int CgHalfEdgeTriangleMesh::getID() const { return m_id; }

#endif // CGHALFEDGETRIANGLEMESH_H
