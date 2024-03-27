#include "CgHalfEdgePrimitives.h"
#include <vector>

CgHeFace::CgHeFace()
{
    m_edge=NULL;
}

CgHeVert::CgHeVert()
{
    m_edge=NULL;
}

CgHeEdge::CgHeEdge()
{
    m_next=NULL;
    m_pair=NULL;
    m_vert=NULL;
    m_face=NULL;
}

CgHeEdge::~CgHeEdge()
{

}

CgHeVert::~CgHeVert()
{

}

CgHeFace::~CgHeFace()
{

}

const CgBaseHeEdge* CgHeFace::edge() const
{
    return (CgBaseHeEdge*) m_edge;
}

const glm::vec3 CgHeFace::normal() const
{
    return m_normal;
}

const CgBaseHeVert* CgHeEdge::vert() const
{
    return (CgBaseHeVert*) m_vert;
}

const CgBaseHeEdge* CgHeEdge::next() const
{
    return (CgBaseHeEdge*) m_next;
}

const CgBaseHeEdge* CgHeEdge::pair() const
{
    return (CgBaseHeEdge*) m_pair;
}

const CgBaseHeFace* CgHeEdge::face() const
{
    return (CgBaseHeFace*) m_face;
}


const CgBaseHeEdge* CgHeVert::edge() const
{
    return (CgBaseHeEdge*) m_edge;
}

const glm::vec3 CgHeVert::position() const
{
    return m_position;
}

const glm::vec3 CgHeVert::color() const
{
    return m_color;
}

void CgHeFace::calculateNormal() {
    CgHeFace* face = this;

    const CgHeEdge* edge0 = (CgHeEdge*)face->edge();
    const CgHeEdge* edge1 = (CgHeEdge*)edge0->next();
    const CgHeEdge* edge2 = (CgHeEdge*)edge1->next();

    const CgHeVert* vert0 = (CgHeVert*)edge0->vert();
    const CgHeVert* vert1 = (CgHeVert*)edge1->vert();
    const CgHeVert* vert2 = (CgHeVert*)edge2->vert();

    glm::vec3 vecA = vert1->position() - vert0->position();
    glm::vec3 vecB = vert2->position() - vert0->position();
    glm::vec3 face_normal = glm::normalize(glm::cross(vecA, vecB));
    face->m_normal = face_normal;
}

/**
 * @brief CgHeVert::calculateNormal calculates the vertex normal
 * and saves the result into the m_normal attribute of the vertex.
 */
void CgHeVert::calculateNormal()
{
    const std::vector<CgHeFace*> faces = getFacesOfVertex();
    glm::vec3 vertex_normal = { 0, 0, 0 };

    // sum up all face normals to get the vertex normal
    for (const CgHeFace* f : faces) {
        vertex_normal += f->normal();
    }
    m_normal = glm::normalize(vertex_normal);
}


/**
 * @brief CgHeVert::getFacesOfVertex finds all faces of a given vertex
 * @return std::vector of all faces
 */
const std::vector<CgHeFace*> CgHeVert::getFacesOfVertex()
{
    // get all edges of a given vert
    std::vector<CgHeEdge*> vert_edges = getEdgesOfVertex();
    std::vector<CgHeFace*> faces;

    // find all faces of the given edges
    for (CgHeEdge* e : vert_edges) {
        CgHeFace* face = (CgHeFace*)e->face();
        faces.push_back(face);
    }
    return faces;
}

/**
 * @brief CgHeVert::getEdgesOfVertex searches all edges which are directly
 * connected to the given vertex.
 * Has a routine for border cases.
 * @return a std::vector of all found edges
 */
const std::vector<CgHeEdge*> CgHeVert::getEdgesOfVertex()
{
    CgHeVert* vert = this;

    // define starting edge
    CgHeEdge* edge = (CgHeEdge*)vert->edge();
    // for border cases
    CgHeEdge* edge_temp = (CgHeEdge*)vert->edge();
    std::vector<CgHeEdge*> edges;
    std::vector<CgHeEdge*> edges_forward;
    std::vector<CgHeEdge*> edges_backward;

    bool isBorder = false;

    // get all edges until you are at the starting point again or at a border
    do {
        if (edge->pair() != NULL) {
            edge = (CgHeEdge*)edge->pair()->next();
            edges_forward.push_back(edge);
        } else {
            isBorder = true;
            break;
        }
    } while (edge != vert->edge());

    // border case: walk into other direction from original starting point
    // to get all other edges
    // Walk into other direction when reaching pair == 0
    // instead of edge-pair-next go edge-next-next-pair
    if (isBorder) {
        do {
            if (edge_temp->next()->next()->pair() != NULL) {

                edge_temp = (CgHeEdge*)edge_temp->next()->next()->pair();
                edges_backward.push_back(edge_temp);
            } else {
                break;
            }
        } while (edge_temp != vert->edge());
    }

    for (int i = edges_backward.size() - 1; i >= 0; i--) {
        CgHeEdge* edge_backward = edges_backward[i];
        edges.push_back(edge_backward);
    }
    edges.push_back((CgHeEdge*)vert->edge());
    for (CgHeEdge* edge_forward : edges_forward) {
        edges.push_back(edge_forward);
    }
    return edges;
}
