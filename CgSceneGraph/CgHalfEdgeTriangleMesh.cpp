#include "CgHalfEdgeTriangleMesh.h"

#include "CgBase/CgEnums.h"
#include "CgUtils/ObjLoader.h"
#include <map>
#include <utility>

CgHalfEdgeTriangleMesh::CgHalfEdgeTriangleMesh()
    : m_type(Cg::HalfEdgeTriangleMesh)
    , m_id(42)
{
    this->loadDemoTriangles();
}
/**
 * @brief CgHalfEdgeTriangleMesh::loadDemoTriangles Function which loads the demo triangles from the beginning */
void CgHalfEdgeTriangleMesh::loadDemoTriangles()
{

    CgHeFace* nf1 = new CgHeFace();
    CgHeFace* nf2 = new CgHeFace();

    CgHeVert* nv1 = new CgHeVert();
    CgHeVert* nv2 = new CgHeVert();
    CgHeVert* nv3 = new CgHeVert();
    CgHeVert* nv4 = new CgHeVert();

    CgHeEdge* n1 = new CgHeEdge();
    CgHeEdge* n2 = new CgHeEdge();
    CgHeEdge* n3 = new CgHeEdge();
    CgHeEdge* n4 = new CgHeEdge();
    CgHeEdge* n5 = new CgHeEdge();
    CgHeEdge* n6 = new CgHeEdge();

    n1->m_next = n2;
    n2->m_next = n3;
    n3->m_next = n1;
    nf1->m_edge = n1;

    n4->m_next = n5;
    n5->m_next = n6;
    n6->m_next = n4;
    nf2->m_edge = n4;

    nv1->m_position = glm::vec3(0.0, 0.0, 0.0);
    nv2->m_position = glm::vec3(0.0, 1.0, 0.0);
    nv3->m_position = glm::vec3(1.0, 1.0, 0.0);
    nv4->m_position = glm::vec3(1.0, 0.0, 0.0);

    nv1->m_edge = n1;
    nv2->m_edge = n2;
    nv3->m_edge = n6;
    nv4->m_edge = n3;

    n4->m_vert = nv4;
    n5->m_vert = nv2;
    n6->m_vert = nv3;

    n1->m_vert = nv1;
    n2->m_vert = nv2;
    n3->m_vert = nv4;

    n2->m_pair = n4;
    n4->m_pair = n2;

    // store into lists
    m_faces.push_back(nf1);
    m_faces.push_back(nf2);

    m_verts.push_back(nv1);
    m_verts.push_back(nv2);
    m_verts.push_back(nv3);
    m_verts.push_back(nv4);

    m_edges.push_back(n1);
    m_edges.push_back(n2);
    m_edges.push_back(n3);
    m_edges.push_back(n4);
    m_edges.push_back(n5);
    m_edges.push_back(n6);

    //attributes
    nv1->m_color = glm::vec3(0.0, 1.0, 1.0);
    nv2->m_color = glm::vec3(0.0, 1.0, 1.0);
    nv3->m_color = glm::vec3(0.0, 1.0, 1.0);
    nv4->m_color = glm::vec3(0.0, 1.0, 1.0);

    nf1->m_normal = glm::vec3(0.0, 0.0, 1.0);
    nf2->m_normal = glm::vec3(0.0, 0.0, 1.0);
}

/**
 * @brief CgHalfEdgeTriangleMesh::loadFromVertexList
 * Loads the vertices and indices of the temporary variables into the half edge data structure.
 * Calculates the vert normals.
 * @param temp_vertices temporary verts which will be saved into the new data structure
 * @param temp_indices temporary indices which will be saved into the new data structure
 */
void CgHalfEdgeTriangleMesh::loadFromVertexList(std::vector<glm::vec3> temp_vertices,
    std::vector<unsigned int> temp_indices)
{
    // iterate through the List of given vertices.
    for (std::size_t i = 0; i < temp_vertices.size(); i++) {
        // save verts from temp list into m_verts. Get position and set color (green for now)
        CgHeVert* vert = new CgHeVert();
        vert->m_position = temp_vertices[i];
        vert->m_color = { 0, 255, 255 };
        m_verts.push_back(vert);
    }

    // Helper structure to accelerate pair lookup
    // Stores the indices of both vertices for an edge as the key to look up the edge
    std::map<std::pair<std::size_t, std::size_t>, CgHeEdge*> edges_by_vertices;

    // iterate through the list of indices to build new half edge data structure
    //
    for (std::size_t i = 0; i < temp_indices.size(); i += 3) {

        // edges of a triangle
        CgHeEdge* edge0 = new CgHeEdge();
        CgHeEdge* edge1 = new CgHeEdge();
        CgHeEdge* edge2 = new CgHeEdge();

        // indices
        unsigned vert_index_0 = temp_indices[i];
        unsigned vert_index_1 = temp_indices[i + 1];
        unsigned vert_index_2 = temp_indices[i + 2];

        // vertices of a triangle
        CgHeVert* vert0 = (CgHeVert*)m_verts[vert_index_0];
        CgHeVert* vert1 = (CgHeVert*)m_verts[vert_index_1];
        CgHeVert* vert2 = (CgHeVert*)m_verts[vert_index_2];

        // face of a triangle
        CgHeFace* face = new CgHeFace();

        // assign vertices to edges
        edge0->m_vert = vert0;
        edge1->m_vert = vert1;
        edge2->m_vert = vert2;

        // assign edges to vertices if not already happend
        if (!vert0->edge())
            vert0->m_edge = edge0;
        if (!vert1->edge())
            vert1->m_edge = edge1;
        if (!vert2->edge())
            vert2->m_edge = edge2;

        // connect the edges
        // Clockwise TODO Maybe change this
        edge0->m_next = edge1;
        edge1->m_next = edge2;
        edge2->m_next = edge0;

        // assign edges to face
        edge0->m_face = face;
        edge1->m_face = face;
        edge2->m_face = face;

        // assign starting edge to face
        face->m_edge = edge0;

        // save edges 1-3 into m_edges
        m_edges.push_back(edge0);
        m_edges.push_back(edge1);
        m_edges.push_back(edge2);

        // Register this edge and its vertices in the pair lookup table
        edges_by_vertices[std::make_pair(vert_index_0, vert_index_1)] = edge0;
        edges_by_vertices[std::make_pair(vert_index_1, vert_index_2)] = edge1;
        edges_by_vertices[std::make_pair(vert_index_2, vert_index_0)] = edge2;

        // Try to load pairs for each new edge from the pair lookup table
        CgHeEdge* pair0 = edges_by_vertices[std::make_pair(vert_index_1, vert_index_0)];
        CgHeEdge* pair1 = edges_by_vertices[std::make_pair(vert_index_2, vert_index_1)];
        CgHeEdge* pair2 = edges_by_vertices[std::make_pair(vert_index_0, vert_index_2)];

        // Connect each new edge with its pair if it exists
        if (pair0) {
            edge0->m_pair = pair0;
            pair0->m_pair = edge0;
        }
        if (pair1) {
            edge1->m_pair = pair1;
            pair1->m_pair = edge1;
        }
        if (pair2) {
            edge2->m_pair = pair2;
            pair2->m_pair = edge2;
        }

        // calculating face normal
        face->calculateNormal();

        // // save face into m_faces
        m_faces.push_back(face);
    }

    // calculate the normal for each vert
    for (CgBaseHeVert* base_vert : m_verts) {
        CgHeVert* vert = (CgHeVert*)base_vert;
        vert->calculateNormal();
    }
}

// Do loopDivision

void CgHalfEdgeTriangleMesh::loopSubdivision()
{
    // for additional verts, edges and faces
    //temp datastructures
    std::vector<CgHeVert*> a_verts;
    std::vector<CgHeEdge*> a_edges;
    std::vector<CgHeFace*> a_faces;

    //=======================Split=======================

    // Clear all temp_split_verts for when used multiple times / intialize with null
    for (CgBaseHeEdge* e : m_edges) {
        CgHeEdge* edge = (CgHeEdge*)e;
        edge->temp_split_vert = NULL;
    }

    // Create a new splitting vertex for all edges
    for (CgBaseHeEdge* base_edge : m_edges) {
        CgHeEdge* edge = (CgHeEdge*)base_edge;

        // if temp_split_vert is not empty
        if (edge->temp_split_vert) {
            continue;
        }

        /* find the four relevant vertices around the edge to split (stencil)
        *     o
        *    / \
        *   /   \
        *  /     \
        *  o--o--o
        *  \     /
        *   \   /
        *    \ /
        *     o
        */
        const CgHeVert* n1 = (CgHeVert*)edge->vert(); // vert connected to the edge
        const CgHeVert* n2 = (CgHeVert*)edge->next()->vert(); // vert connected to the edge

        CgHeEdge* pair = (CgHeEdge*)edge->pair();

        CgHeVert* split_vert = new CgHeVert();

        // if pair is not empty
        if (pair) {
            // Find the two relevant references to the vertices (Start and end of an Edge through the other edges)
            const CgHeVert* f1 = (CgHeVert*)edge->next()->next()->vert(); // vert not directly connected to the edge but relevant
            const CgHeVert* f2 = (CgHeVert*)edge->pair()->next()->next()->vert(); // vert not directly connected to the edge but relevant

            // calculate the Position of the new vert
            split_vert->m_position = (3.0f / 8.0f) * n1->position() + (3.0f / 8.0f) * n2->position() + (1.0f / 8.0f) * f1->position() + (1.0f / 8.0f) * f2->position();

            // calculate the color of the new vert
            split_vert->m_color = (3.0f / 8.0f) * n1->color() + (3.0f / 8.0f) * n2->color() + (1.0f / 8.0f) * f1->color() + (1.0f / 8.0f) * f2->color();

            // save the new vert into the edge
            edge->temp_split_vert = split_vert;

            // save the new vert as a pair
            pair->temp_split_vert = split_vert;

            // if pair is empty
        } else {
            // calculate the Position of the new vert
            split_vert->m_position = (1.0f / 2.0f) * n1->position() + (1.0f / 2.0f) * n2->position();

            // calculate the color of the new vert
            split_vert->m_color = (1.0f / 2.0f) * n1->color() + (1.0f / 2.0f) * n2->color();

            // save the new vert into the edge
            edge->temp_split_vert = split_vert;
        }

        // save new vert into temp datastructure
        a_verts.push_back(split_vert);
    }
    //=======================Average=======================

    // calculate adjusted positions for old vertices
    std::map<int, float> betas; // cache for already calculated betas
    for (CgBaseHeVert* base_vert : m_verts) {
        CgHeVert* vert = (CgHeVert*)base_vert;

        // get all neighbor vertices
        bool isBorderVert = false;
        std::vector<CgHeVert*> n_verts = vert->getNeighborVerts(&isBorderVert);
        // if vert is at the border
        if (isBorderVert) {
            // get first and last verts of all neighboring verts
            const CgHeVert* first_vert = n_verts[0];
            const CgHeVert* last_vert = n_verts[n_verts.size() - 1];

            // calculate new position
            // therefore use positions of first and last verts of neighbours
            // as well as the position of current vert
            vert->temp_position = (1.0f / 8.0f) * first_vert->position() + (1.0f / 8.0f) * last_vert->position() + (6.0f / 8.0f) * vert->position();

            // if vert is not at the border
        } else {
            // number of neighbors
            int n = n_verts.size();
            // if betas[n] is empty
            if (!betas[n]) {
                // calculate beta
                betas[n] = calculateLoopSubdivisionBeta(n);
            }
            float beta = betas[n];

            // calculate new position of old verts
            glm::vec3 new_pos = { 0, 0, 0 };
            for (const CgHeVert* n_vert : n_verts) {
                // beta = 1/n * (5/8 - (3/8 + 1/4 * cos((2*pi)/n))²)
                new_pos += n_vert->position() * beta;
            }
            new_pos += vert->position() * (1.0f - (n * beta));
            vert->temp_position = new_pos;
        }
    }

    // connect splitting verts, one face at a time
    for (CgBaseHeFace* base_face : m_faces) {
        CgHeFace* face = (CgHeFace*)base_face;

        CgHeEdge* oe1 = (CgHeEdge*)face->edge(); // original outer edge (outer vert1 to outer vert2)
        CgHeEdge* oe2 = (CgHeEdge*)oe1->next(); // original outer edge (outer vert2 to outer vert3)
        CgHeEdge* oe3 = (CgHeEdge*)oe2->next(); // original outer edge (outer vert3 to outer vert1)

        CgHeEdge* oe_ext1 = new CgHeEdge(); // outer edge from inner vert1 to vert2
        CgHeEdge* oe_ext2 = new CgHeEdge(); // outer edge from inner vert2 to vert3
        CgHeEdge* oe_ext3 = new CgHeEdge(); // outer edge from inner vert3 to vert1

        CgHeVert* ov1 = (CgHeVert*)oe1->vert(); // An original outer vertex
        CgHeVert* ov2 = (CgHeVert*)oe2->vert(); // An original outer vertex
        CgHeVert* ov3 = (CgHeVert*)oe3->vert(); // An original outer vertex

        CgHeVert* iv1 = oe1->temp_split_vert; // new vertex splitting outer edge1
        CgHeVert* iv2 = oe2->temp_split_vert; // new vertex splitting outer edge2
        CgHeVert* iv3 = oe3->temp_split_vert; // new vertex splitting outer edge3

        CgHeEdge* ie1 = new CgHeEdge(); // inner edge from inner vert1 to split vert2
        CgHeEdge* ie2 = new CgHeEdge(); // inner edge from inner vert2 to split vert3
        CgHeEdge* ie3 = new CgHeEdge(); // inner edge from inner vert3 to split vert1

        CgHeEdge* ie_pair1 = new CgHeEdge(); // inner edge from inner vert1 to split vert2
        CgHeEdge* ie_pair2 = new CgHeEdge(); // inner edge from inner vert2 to split vert3
        CgHeEdge* ie_pair3 = new CgHeEdge(); // inner edge from inner vert3 to split vert1

        CgHeFace* face1 = new CgHeFace(); // new face next to outer vert 1
        CgHeFace* face2 = new CgHeFace(); // new face next to outer vert 2
        CgHeFace* face3 = new CgHeFace(); // new face next to outer vert 3
        CgHeFace* center_face = face; // original face that now is the new center face

        // add elements to our temp datastructures
        a_edges.push_back(oe_ext1);
        a_edges.push_back(oe_ext2);
        a_edges.push_back(oe_ext3);

        a_edges.push_back(ie1);
        a_edges.push_back(ie2);
        a_edges.push_back(ie3);

        a_edges.push_back(ie_pair1);
        a_edges.push_back(ie_pair2);
        a_edges.push_back(ie_pair3);

        a_faces.push_back(face1);
        a_faces.push_back(face2);
        a_faces.push_back(face3);

        // configuration of face elements (+ reconfiguration of existing elements)
        iv1->m_edge = oe_ext1;
        iv2->m_edge = oe_ext2;
        iv3->m_edge = oe_ext3;

        oe1->m_face = face1;
        oe2->m_face = face2;
        oe3->m_face = face3;

        oe1->m_next = ie_pair3;
        oe2->m_next = ie_pair1;
        oe3->m_next = ie_pair2;

        oe_ext1->m_vert = iv1;
        oe_ext2->m_vert = iv2;
        oe_ext3->m_vert = iv3;

        oe_ext1->m_pair = NULL;
        oe_ext2->m_pair = NULL;
        oe_ext3->m_pair = NULL;

        // The old edges maintain their pairing, thus no NULL

        oe_ext1->m_face = face2;
        oe_ext2->m_face = face3;
        oe_ext3->m_face = face1;

        oe_ext1->m_next = oe2;
        oe_ext2->m_next = oe3;
        oe_ext3->m_next = oe1;

        ie1->m_vert = iv1;
        ie2->m_vert = iv2;
        ie3->m_vert = iv3;

        ie1->m_pair = ie_pair1;
        ie2->m_pair = ie_pair2;
        ie3->m_pair = ie_pair3;

        ie1->m_face = center_face;
        ie2->m_face = center_face;
        ie3->m_face = center_face;

        ie1->m_next = ie2;
        ie2->m_next = ie3;
        ie3->m_next = ie1;

        ie_pair1->m_vert = iv2;
        ie_pair2->m_vert = iv3;
        ie_pair3->m_vert = iv1;

        ie_pair1->m_pair = ie1;
        ie_pair2->m_pair = ie2;
        ie_pair3->m_pair = ie3;

        ie_pair1->m_face = face2;
        ie_pair2->m_face = face3;
        ie_pair3->m_face = face1;

        ie_pair1->m_next = oe_ext1;
        ie_pair2->m_next = oe_ext2;
        ie_pair3->m_next = oe_ext3;

        face1->m_edge = oe1;
        face2->m_edge = oe2;
        face3->m_edge = oe3;

        center_face->m_edge = ie1;
    }

    // fix pair connections for external edges
    // starting with the original edges ==> still have useful pair edge set (now pair for extension edge)
    for (CgBaseHeEdge* base_edge : m_edges) {
        CgHeEdge* edge = (CgHeEdge*)base_edge;

        CgHeEdge* ep = (CgHeEdge*)edge->pair(); // get original pair before split (extensionpair)-

        // if original pair is empty, there will not be a new pair
        if (!ep) {
            continue;
        }

        CgHeEdge* e_edge = (CgHeEdge*)edge->next()->pair()->next()->pair()->next(); // walk through internals of new triangle to reach extension edge
        CgHeEdge* pair = (CgHeEdge*)ep->next()->pair()->next()->pair()->next(); // walk through internals of new triangle to reach extension edge

        edge->m_pair = pair;
        e_edge->m_pair = ep;
    }

    // position adjustment for all original verts
    for (CgBaseHeVert* base_vert : m_verts) {
        CgHeVert* vert = (CgHeVert*)base_vert;

        vert->m_position = vert->temp_position;
    }

    m_verts.insert(m_verts.end(), a_verts.begin(), a_verts.end());
    m_edges.insert(m_edges.end(), a_edges.begin(), a_edges.end());
    m_faces.insert(m_faces.end(), a_faces.begin(), a_faces.end());

    // re-calculate all face normals
    for (CgBaseHeFace* base_face : m_faces) {
        CgHeFace* face = (CgHeFace*)base_face;
        face->calculateNormal();
    }

    // re-calculate all vertex normals
    for (CgBaseHeVert* base_vert : m_verts) {
        CgHeVert* vert = (CgHeVert*)base_vert;
        vert->calculateNormal();
    }
}

CgHalfEdgeTriangleMesh::CgHalfEdgeTriangleMesh(int id)
    : m_type(Cg::HalfEdgeTriangleMesh)
    , m_id(id)
{
}

CgHalfEdgeTriangleMesh::~CgHalfEdgeTriangleMesh()
{
    // delete pointers and clear the vector

    for (CgBaseHeFace* face : m_faces) {
        delete face;
    }
    m_faces.clear();
    for (CgBaseHeEdge* edge : m_edges) {
        delete edge;
    }
    m_edges.clear();
    for (CgBaseHeVert* vert : m_verts) {
        delete vert;
    }
    m_verts.clear();
}

const std::vector<CgBaseHeFace*>& CgHalfEdgeTriangleMesh::getFaces() const
{
    return m_faces;
}

void CgHalfEdgeTriangleMesh::init(std::string filename)
{
    std::vector<glm::vec3> temp_vertices;
    std::vector<glm::vec3> temp_vertnormals;
    std::vector<unsigned int> temp_indices;

    m_verts.clear();
    m_edges.clear();
    m_faces.clear();

    ObjLoader loader;
    loader.load(filename);

    loader.getPositionData(temp_vertices);
    loader.getNormalData(temp_vertnormals);
    loader.getFaceIndexData(temp_indices);

    this->loadFromVertexList(temp_vertices, temp_indices);
}

const glm::vec3 CgHalfEdgeTriangleMesh::getCenter() const
{
    glm::vec3 center(0.);
    for (unsigned int i = 0; i < m_verts.size(); i++) {
        center += m_verts[i]->position();
    }
    center /= (double)m_verts.size();
    return center;
}

/**
 * @brief CgHalfEdgeTriangleMesh::calculateLoopSubdevisionBeta
 * Calculates the beta value of the Loop Subdevision for a given n.
 * @param n number of old Vertices
 * @return the calculation result (beta)
 */
float CgHalfEdgeTriangleMesh::calculateLoopSubdivisionBeta(int n)
{
    return (1.0 / n) * ((5.0 / 8.0) - pow((3.0 / 8.0) + (1.0 / 4.0) * cos(2.0 * M_PIf64 / n), 2.0));
}
