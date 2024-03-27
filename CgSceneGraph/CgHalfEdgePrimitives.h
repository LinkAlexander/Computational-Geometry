#ifndef CGHALFEDGEPRIMITIVES_H
#define CGHALFEDGEPRIMITIVES_H

#include <glm/glm.hpp>
#include "CgBase/CgBaseHalfdgePrimitives.h"
#include <vector>

//forward declaration
class CgHeVert;
class CgHeEdge;

class CgHeFace : public CgBaseHeFace
{
public:
    CgHeFace();
    ~CgHeFace();

    const CgBaseHeEdge* edge() const;
    const glm::vec3 normal() const;
    void calculateNormal();
    CgHeEdge* m_edge;

    glm::vec3 m_normal;


};





class CgHeVert : public CgBaseHeVert
{
public:
    CgHeVert();
    ~CgHeVert();

    const CgBaseHeEdge* edge() const;
    const glm::vec3 position() const;
    const glm::vec3 color() const;
    void calculateNormal();

    const std::vector<CgHeFace*> getFacesOfVertex();
    const std::vector<CgHeEdge*> getEdgesOfVertex();
    const std::vector<CgHeVert*> getNeighborVerts(bool *isBorderVert);

    CgHeEdge* m_edge;
    glm::vec3 m_position;
    glm::vec3 m_color;
    glm::vec3 m_normal;

};

class CgHeEdge : public CgBaseHeEdge
{

public:

    CgHeEdge();
    ~CgHeEdge();

    const CgBaseHeEdge* next() const;
    const CgBaseHeEdge* pair() const;
    const CgBaseHeVert* vert() const;
    const CgBaseHeFace* face() const;


    CgHeEdge* m_next;
    CgHeEdge* m_pair;
    CgHeVert* m_vert;
    CgHeFace* m_face;

};




#endif // CGHALFEDGEPRIMITIVES_H
