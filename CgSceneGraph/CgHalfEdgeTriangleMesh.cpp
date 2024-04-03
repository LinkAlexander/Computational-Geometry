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
    nv1->m_color = glm::vec3(0.0, 1.0, 0.0);
    nv2->m_color = glm::vec3(0.0, 1.0, 0.0);
    nv3->m_color = glm::vec3(0.0, 1.0, 0.0);
    nv4->m_color = glm::vec3(0.0, 1.0, 0.0);

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

        // save face into m_faces
        m_faces.push_back(face);
    }

    // calculate the normal for each vert
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
