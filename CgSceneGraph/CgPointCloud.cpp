#include "CgPointCloud.h"
#include "CgBase/CgEnums.h"
#include "CgUtils/ObjLoader.h"
#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>


CgPointCloud::CgPointCloud():
CgPointCloud::CgPointCloud(51)
{

}

CgPointCloud::CgPointCloud(int id):
m_type(Cg::PointCloud),
m_id(id)
{

    m_vertices.push_back(glm::vec3(0.0,0.0,0.0));
    m_vertex_normals.push_back(glm::vec3(0.0,0.0,1.0));
    m_vertex_colors.push_back(glm::vec3(0.0,0.0,1.0));

    calculateSplatOrientations();

}


CgPointCloud::~CgPointCloud()
{
    m_vertices.clear();
    m_vertex_normals.clear();
    m_vertex_colors.clear();
    m_splat_indices.clear();
}


void CgPointCloud::calculateSplatOrientations()
{
  // calculate local coordinate system for splats (arbitrary orientation of ellipse in plane)
  // replace this if you have the real coordinate system, use up vector = y-Axis of your local coordinate system instead of getPerpendicularVector(...)

  m_splat_orientations.clear();
  m_splat_scaling.clear();
  m_splat_indices.clear();

  for(unsigned int i=0;i<m_vertices.size();i++)
    {
      glm::mat4 lookAt_matrix(glm::lookAt(glm::vec3(m_vertices[i]),glm::vec3(m_vertices[i]-m_vertex_normals[i]),getPerpendicularVector(m_vertex_normals[i])));
      m_splat_orientations.push_back(lookAt_matrix);
      m_splat_scaling.push_back(glm::vec2(0.02,0.005));

      // use all points for splatting by default
      m_splat_indices.push_back(i);
    }


}


void CgPointCloud::init( std::string filename, bool cheat_normals)
{
    m_vertices.clear();
    m_vertex_normals.clear();
    m_vertex_colors.clear();
    m_splat_orientations.clear();
    m_splat_scaling.clear();
    m_splat_indices.clear();

    // load obj File
    ObjLoader loader;
    loader.load(filename);
    loader.getPositionData(m_vertices);


    // do this for cheating with the normals
    // you need to replace this by a normal estimation algorithm
    if(cheat_normals)
      loader.getNormalData(m_vertex_normals);


    // calculate local coordinate system for splats (arbitrary orientation of ellipse in plane)
    // replace this if you have the real coordinate system, use up vector = y-Axis of your local coordinate system instead of getPerpendicularVector(...)

    calculateSplatOrientations();

    //add a standard color for each point if lighting turned off
    for(unsigned int i=0;i<m_vertices.size();i++)
      {
          m_vertex_colors.push_back(glm::vec3(0.0,1.0,0.0));
      }


    //test of getNeartestNeighbors(..) method
    // generates blue dots on the tail of the bunny

    // unsigned int k=50;
    // std::vector<int> neighbors = getNearestNeighbors(0,k);

    // for(unsigned int i=0;i<k || i < m_vertex_colors.size();i++)
    //   {
    //     m_vertex_colors[neighbors[i]]=glm::vec3(0.0,0.0,1.0);
    //   }
}


std::vector<int> CgPointCloud::getNearestNeighbors(int current_point,unsigned int k)
{

  glm::vec3 q= m_vertices[current_point];

  std::vector<std::pair<double,int>> distances;

  // very inefficient, just to show that it works for rendering colored neighborhood
  // use min heap for real purposes


  for(unsigned int i=0;i<m_vertices.size();i++)
    {
      double dist=glm::distance(m_vertices[i],q);

      distances.push_back(std::make_pair(dist,i));
    }

    std::sort(distances.begin(),distances.end());

    std::vector<int> erg;

   for(unsigned int i=0;i<k;i++)
    {
       erg.push_back(distances[i].second);
     }

    return erg;
  }


// calculates an arbitrary verctor perpendicular to the given one
glm::vec3 CgPointCloud::getPerpendicularVector(glm::vec3 arg)
{
  if((arg[0]==0.0)&&(arg[1]==0.0))
    {
    if(arg[2]==0.0)
      return glm::vec3(0.);

    return glm::vec3(0.0,1.0,0.0);
    }
  return glm::normalize(glm::vec3(-arg[1],arg[0],0.0));
}




const glm::vec3 CgPointCloud::getCenter() const
{
  glm::vec3 center(0.);
  for(unsigned int i=0;i<m_vertices.size();i++)
    {
      center+=m_vertices[i];
    }
  center/=(double)m_vertices.size();
  return center;
}

const std::vector<glm::vec2>& CgPointCloud::getSplatScalings() const
{
  return m_splat_scaling;
}
