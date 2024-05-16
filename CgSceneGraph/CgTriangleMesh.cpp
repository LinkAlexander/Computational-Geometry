#include "CgTriangleMesh.h"
#include "CgBase/CgEnums.h"
#include "CgUtils/ObjLoader.h"

CgTriangleMesh::CgTriangleMesh():
m_type(Cg::TriangleMesh),
m_id(42)
{
}

CgTriangleMesh::CgTriangleMesh(int id):
m_type(Cg::TriangleMesh),
m_id(id)
{
    // ellipse (a* cos(t) , b*sin (t)) if axes are x and y with 0 <= t <= 2 pi

    m_vertices.push_back(glm::vec3(0.0,0.0,0.0));
    m_vertex_normals.push_back(glm::vec3(0.0,0.0,1.0));
    m_vertex_colors.push_back(glm::vec3(0.0,1.0,0.0));

    unsigned int num_steps=20;
    double stepsize=2.0*M_PI/(double)num_steps;


    for(unsigned int i=0;i<num_steps;i++)
    {
        m_vertices.push_back(glm::vec3(cos(i*stepsize),sin(i*stepsize),0.0));
        m_vertex_normals.push_back(glm::vec3(0.0,0.0,1.0));
    }

    for(unsigned int i=0;i<num_steps;i++)
    {
        m_triangle_indices.push_back(0);
        m_triangle_indices.push_back(i);
        m_triangle_indices.push_back(i+1);
      }

    m_triangle_indices.push_back(0);
    m_triangle_indices.push_back(num_steps);
    m_triangle_indices.push_back(1);


    for(unsigned int i=0;i<m_triangle_indices.size();i+=3)
    {
         glm::vec3 v1 = m_vertices[m_triangle_indices[i+1]]-m_vertices[m_triangle_indices[i]];
         glm::vec3 v2 = m_vertices[m_triangle_indices[i+2]]-m_vertices[m_triangle_indices[i]];

         glm::vec3 normal = glm::cross(v1,v2);

         m_face_normals.push_back(normal);
    }


}



CgTriangleMesh::~CgTriangleMesh()
{
    m_vertices.clear();
    m_vertex_normals.clear();
    m_vertex_colors.clear();
    m_tex_coords.clear();
    m_triangle_indices.clear();
    m_face_normals.clear();
    m_face_colors.clear();
    m_vertex_colors.clear();
}


void CgTriangleMesh::init( std::string filename)
{
    m_vertices.clear();
    m_vertex_normals.clear();
    m_triangle_indices.clear();
    m_triangle_indices.clear();
    m_face_normals.clear();


    ObjLoader loader;
    loader.load(filename);

    loader.getPositionData(m_vertices);
    loader.getNormalData(m_vertex_normals);
    loader.getFaceIndexData(m_triangle_indices);

    for(unsigned int i=0;i<m_triangle_indices.size();i+=3)
    {
         glm::vec3 v1 = m_vertices[m_triangle_indices[i+1]]-m_vertices[m_triangle_indices[i]];
         glm::vec3 v2 = m_vertices[m_triangle_indices[i+2]]-m_vertices[m_triangle_indices[i]];

         glm::vec3 normal = glm::cross(v1,v2);

         m_face_normals.push_back(normal);
    }

    for(unsigned int i=0;i<m_vertices.size();i++)
    {
      m_vertex_colors.push_back(glm::vec3(0.0,1.0,1.0));
    }
}


const glm::vec3 CgTriangleMesh::getCenter() const
{
  glm::vec3 center(0.);
  for(unsigned int i=0;i<m_vertices.size();i++)
    {
      center+=m_vertices[i];
    }
  center/=(double)m_vertices.size();
  return center;
}

const std::vector<glm::vec3>& CgTriangleMesh::getVertices() const
{
    return m_vertices;
}

const std::vector<glm::vec3>& CgTriangleMesh::getVertexNormals() const
{
    return m_vertex_normals;
}

const std::vector<glm::vec3>& CgTriangleMesh::getVertexColors() const
{
     return m_vertex_colors;
}

const std::vector<glm::vec2>& CgTriangleMesh:: getVertexTexCoords() const
{
    return m_tex_coords;
}

const std::vector<unsigned int>& CgTriangleMesh::getTriangleIndices() const
{
    return m_triangle_indices;
}

const std::vector<glm::vec3>& CgTriangleMesh::getFaceNormals() const
{
    return m_face_normals;
}

const std::vector<glm::vec3>& CgTriangleMesh::getFaceColors() const
{
    return m_face_colors;
}

void CgTriangleMesh::addQuadrangle(glm::vec3 pos1, glm::vec3 pos2, glm::vec3 pos3, glm::vec3 pos4, glm::vec3 color)
{
    // Face1: pos1, pos2, pos3
    // Face2: pos3, pos4, pos1
    glm::vec3 faceNormal1 = glm::normalize(glm::cross(pos2 - pos1, pos3 - pos1));
    glm::vec3 faceNormal2 = glm::normalize(glm::cross(pos4 - pos3, pos1 - pos3));

    size_t indexOffset = m_vertices.size() - 1;

    m_vertices.push_back(pos1);
    m_vertices.push_back(pos2);
    m_vertices.push_back(pos3);
    m_vertices.push_back(pos4);
    m_vertex_normals.push_back((faceNormal1 + faceNormal2) / 2.0f);
    m_vertex_normals.push_back(faceNormal1);
    m_vertex_normals.push_back(faceNormal1);
    m_vertex_normals.push_back((faceNormal1 + faceNormal2) / 2.0f);
    m_vertex_colors.push_back(color);
    m_vertex_colors.push_back(color);
    m_vertex_colors.push_back(color);
    m_vertex_colors.push_back(color);

    m_triangle_indices.push_back(indexOffset + 1);
    m_triangle_indices.push_back(indexOffset + 2);
    m_triangle_indices.push_back(indexOffset + 3);
    m_triangle_indices.push_back(indexOffset + 3);
    m_triangle_indices.push_back(indexOffset + 4);
    m_triangle_indices.push_back(indexOffset + 1);
    m_face_normals.push_back(faceNormal1);
    m_face_normals.push_back(faceNormal2);
    m_face_colors.push_back(color);
    m_face_colors.push_back(color);
}

/**
 * @brief CgTriangleMesh::addQuadrangleMesh Add the given vertices and add edges to connect them in a mesh, according to the given matrix from
 *
 * Structure of the grid:
 * (doubled lines represent the responsibility for creation of edges for vertex O)
 *
 * X----X----X----X----X
 * |  / | //||  / |  / |
 * | /  |// || /  | /  |
 * X----O====X----X----X
 * |  /|| // |  / |  / |
 * | / ||//  | /  | /  |
 * X----X----X----X----X
 *
 * @param vertices a rectangular set of vertices to add as a mesh
 * @param rows number of rows for given vertices
 * @param cols number of columns for given vertices
 */
void CgTriangleMesh::addQuadrangleMesh(glm::vec3* vertices, size_t rows, size_t cols, glm::vec3 normal, glm::vec3 color)
{
    size_t indexOffset = m_vertices.size();

    for (size_t row = 0; row < rows; row++) {
        for (size_t col = 0; col < cols; col++) {
            glm::vec3 vertex = vertices[row * cols + col];
            m_vertices.push_back(vertex);
            m_vertex_normals.push_back(normal);
            m_vertex_colors.push_back(color);

            size_t currentIndex = row * cols + col + indexOffset;

            // Skip triangle creation for last column,
            // as the triangles are always created to the right
            if (col != cols - 1) {
                // Create up-right triangle
                // Skip for first row
                if (row != 0) {
                    m_triangle_indices.push_back(currentIndex);
                    m_triangle_indices.push_back(currentIndex + 1);
                    m_triangle_indices.push_back(currentIndex - cols + 1);

                    m_face_normals.push_back(normal);
                    m_face_colors.push_back(color);
                }

                // Create down-triangle
                // Skip for last row
                if (row != rows - 1) {
                    m_triangle_indices.push_back(currentIndex);
                    m_triangle_indices.push_back(currentIndex + cols);
                    m_triangle_indices.push_back(currentIndex + 1);

                    m_face_normals.push_back(normal);
                    m_face_colors.push_back(color);
                }
            }
        }
    }
}