#ifndef CG_POLYLINE_H
#define CG_POLYLINE_H

#include <vector>
#include <glm/glm.hpp>
#include <string>
#include "CgBase/CgBasePolyline.h"

class CgPolyLine : public CgBasePolyline
{
public:
  CgPolyLine();
  CgPolyLine(int);
  CgPolyLine(int, std::vector<glm::vec3>);

  ~CgPolyLine();

  //inherited from CgBaseRenderableObject
  Cg::ObjectType getType() const;
  unsigned int   getID()   const;


  //inherited from CgBasePolyline
  const std::vector<glm::vec3>& getVertices()  const;
  const std::vector<glm::vec3>& getVertexColors()     const;
  unsigned int                  getLineWidth() const;

private:

  const Cg::ObjectType m_type;
  const unsigned int m_id;

  std::vector<glm::vec3> m_vertices;
  std::vector<glm::vec3> m_vertex_colors;
  unsigned int m_linewidth;
};


inline Cg::ObjectType  CgPolyLine::getType() const {return m_type;}
inline unsigned int CgPolyLine::getID() const {return m_id;}




#endif
