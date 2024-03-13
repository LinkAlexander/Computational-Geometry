#include "CgPolyLine.h"
#include "CgBase/CgEnums.h"

CgPolyLine::CgPolyLine():
m_type(Cg::Polyline),
m_id(23)
{
}

CgPolyLine::CgPolyLine(int id):
m_type(Cg::Polyline),
m_id(id),
m_linewidth(1.0)
{
}

CgPolyLine::CgPolyLine(int id,std::vector<glm::vec3> points):
m_type(Cg::Polyline),
m_id(id),
m_linewidth(1.0)
{
  m_vertices.clear();
  m_vertices=points;

  for(unsigned int i=0;i<m_vertices.size();i++)
    {
    m_vertex_colors.push_back(glm::vec3(0.0,0.0,1.0));
    }
}

CgPolyLine::~CgPolyLine()
{
    m_vertices.clear();
    m_vertex_colors.clear();
}

const std::vector<glm::vec3>& CgPolyLine::getVertices() const {return m_vertices;}
const std::vector<glm::vec3>& CgPolyLine::getVertexColors() const{return m_vertex_colors;}
unsigned int CgPolyLine::getLineWidth() const {return m_linewidth;}


