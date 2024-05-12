#ifndef CGEXAMPLETRIANGLE_H
#define CGEXAMPLETRIANGLE_H

#include <vector>
#include <glm/glm.hpp>
#include <string>
#include "CgBase/CgBaseTriangleMesh.h"


class CgTriangleMesh : public CgBaseTriangleMesh
{

public:
    CgTriangleMesh();
    CgTriangleMesh(int id);

    ~CgTriangleMesh();

    //inherited from CgBaseRenderableObject
    Cg::ObjectType getType() const;
    unsigned int getID() const;

    //inherited from CgBaseTriangleMesh
    const std::vector<glm::vec3>& getVertices() const;
    const std::vector<glm::vec3>& getVertexNormals() const;
    const std::vector<glm::vec3>& getVertexColors() const;
    const std::vector<glm::vec2>& getVertexTexCoords() const;

    const std::vector<unsigned int>& getTriangleIndices() const;

    const std::vector<glm::vec3>& getFaceNormals() const;
    const std::vector<glm::vec3>& getFaceColors() const;

    //own stuff

    void init( std::string filename);
    const glm::vec3 getCenter() const;



private:

    std::vector<glm::vec3> m_vertices;
    std::vector<glm::vec3> m_vertex_normals;
    std::vector<glm::vec3> m_vertex_colors;
    std::vector<glm::vec2> m_tex_coords;

    std::vector<unsigned int>  m_triangle_indices;

    std::vector<glm::vec3> m_face_normals;
    std::vector<glm::vec3> m_face_colors;

    const Cg::ObjectType m_type;
    const unsigned int m_id;

};


inline Cg::ObjectType  CgTriangleMesh::getType() const {return m_type;}
inline unsigned int CgTriangleMesh::getID() const {return m_id;}

#endif // CGEXAMPLETRIANGLE_H
