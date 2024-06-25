#include "CgPointCloud.h"
#include "CgBase/CgEnums.h"
#include "CgUtils/ObjLoader.h"
#include "CgKdTree.h"
#include "CgMovingLeastSquares.h"
#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>
#include <iostream>


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
    m_vertex_colors.push_back(glm::vec3(0.0,1.0,1.0));

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


    // initialising kdTree
    kdTree = new CgKdTree(this, 9);

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
          m_vertex_colors.push_back(glm::vec3(0.0,1.0,1.0));
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




std::vector<size_t> CgPointCloud::getNearestNeighbors(int current_point, unsigned int k)
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

    std::vector<size_t> erg;

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

/**
 * @brief CgPointCloud::applyPickRay Apply the given pick-ray to this point cloud
 * @param pickRayStart The start of the pick-ray
 * @param pickRayDirection The direction of the pick-ray
 */
void CgPointCloud::applyPickRay(glm::vec3 pickRayStart, glm::vec3 pickRayDirection)
{
    // Get the selected point closest to the pick ray
    size_t centerIndex = kdTree->getClosestPointToRay(pickRayStart, pickRayDirection);

    unsigned int k = 50;

    if(this->m_vertices.size() < 50) {
        k = 3;
    }
    // Get the nearest neighbors of the selected point
    std::vector<size_t> neighbors = getNearestNeighbors(centerIndex, k);

    // Color all vertices cyan
    for (glm::vec3& color : m_vertex_colors) {
        color = { 0.0, 1.0, 1.0 };
    }

    // Color all neighbors red
    for (unsigned int i = 0; i < k; i++) {
        m_vertex_colors[neighbors[i]] = glm::vec3(1.0, 0.0, 0.0);
    }
}

/**
 * @brief distanceRayToPoint Calculates the distance of the picking ray to the selected point
 * @param rayStart initial point of the ray
 * @param rayDirection direction of the ray
 * @param point selected point
 * @return the distance between point and ray
 */
float CgPointCloud::distanceRayToPoint(glm::vec3 rayStart, glm::vec3 rayDirection, glm::vec3 point)
{
    return glm::length(glm::cross(rayDirection, point - rayStart));
}

/**
 * @brief getClosestPointToRay Returns the index of the point which is the nearest to the picking ray
 * @param rayStart initial point of the ray
 * @param rayDirection direction of the ray
 * @return index of the point which is the nearest to the picking ray
 */
size_t CgPointCloud::getClosestPointToRay(glm::vec3 rayStart, glm::vec3 rayDirection)
{
    // index and distance which are closest to the ray
    size_t closestIndex = 0;
    float closestDist = INFINITY;

    // go through all positions and look if there are any which are closer to the ray than the
    // previously found
    for (size_t i = 0; i < m_vertices.size(); i++) {
        glm::vec3 position = m_vertices[i];
        float dist = distanceRayToPoint(rayStart, rayDirection, position);
        if (dist < closestDist) {
            closestIndex = i;
            closestDist = dist;
        }
    }
    return closestIndex;
}

const std::vector<glm::vec2>& CgPointCloud::getSplatScalings() const
{
  return m_splat_scaling;
}

/**
 * @brief CgPointCloud::getSplitPlanes Get the split planes of the KD-Tree for this point cloud
 * @param maxDepth The max depth of split planes to get
 * @return The split planes of the KD-Tree for this point cloud
 */
std::vector<SplitPlane*> CgPointCloud::getSplitPlanes(size_t maxDepth)
{
    return kdTree->getSplitPlanes(maxDepth);
}

/**
 * @brief CgPointCloud::smoothPoint Smooth the given point based on its neighbors
 * @param centerIndex Index of the center point to smooth
 * @param neighborCount The amount of nearest neighbors to use for smoothing
 * @param bivariateFunctionDegree The degree for the bivariate function to smooth with
 * @return The smoothed position of the given point
 */
glm::vec3 CgPointCloud::smoothPoint(size_t centerIndex, size_t neighborCount, size_t bivariateFunctionDegree)
{
    // Get the neirest neighbors
    std::vector<size_t> n_indices = getNearestNeighbors(centerIndex, neighborCount);

    // Get the position of the center point to smooth
    glm::vec3 center_pos = m_vertices[centerIndex];

    // Load the neighbor positions
    std::vector<glm::vec3> neighborPositions(n_indices.size());

    for (size_t i = 0; i < n_indices.size(); i++) {
        size_t n_index = n_indices[i];
        glm::vec3 n_pos = m_vertices[n_index];
        neighborPositions[i] = n_pos;
    }

    // Calculate the smoothed position for the given center point
    CgMovingLeastSquares movingLeastSquares;
    glm::vec3 centerPositionUpdated = movingLeastSquares.smoothBivariate(
            center_pos, neighborPositions, bivariateFunctionDegree);
//    std::cout << centerPositionUpdated - center_pos << std::endl;
    return centerPositionUpdated;
}

/**
 * @brief CgPointCloud::smoothSelectedPoint Generate a triangle mesh to plot a fitted bivariate function for a selected point
 * @param pickRayStart The start of the pick-ray
 * @param pickRayDirection The direction of the pick-ray
 * @param neighborCount The amount of nearest neighbors to use for smoothing
 * @param bivariateFunctionDegree The degree for the bivariate function to smooth with
 * @return The plotted bivariate function for the given selected point
 */
CgTriangleMesh* CgPointCloud::smoothSelectedPoint(glm::vec3 pickRayStart, glm::vec3 pickRayDirection, size_t neighborCount, size_t bivariateFunctionDegree, size_t plottingSteps)
{
    // Get the selected point closest to the pick ray
    size_t centerIndex = kdTree->getClosestPointToRay(pickRayStart, pickRayDirection);

    // Get the nearest neighbor for the selected point
    std::vector<size_t> neighborIndices = getNearestNeighbors(centerIndex, neighborCount);

    // Get the position of the center point to smooth
    glm::vec3 center = m_vertices[centerIndex];

    // Load the neighbor positions
    std::vector<glm::vec3> neighbors(neighborIndices.size());

    for (size_t i = 0; i < neighborIndices.size(); i++) {
        size_t neighborIndex = neighborIndices[i];
        glm::vec3 neighborPosition = m_vertices[neighborIndex];
        neighbors[i] = neighborPosition;
    }

    Plane tangentPlane = CgMovingLeastSquares::covarianceAnalysis(center, neighbors);

    // Generate the sampling points and values for fitting based on the neighbor positions
    std::vector<glm::vec2> samplingPoints(neighbors.size());
    std::vector<float> samplingValues(neighbors.size());
    for (size_t i = 0; i < neighbors.size(); i++) {
        glm::vec3 neighbor = neighbors[i];
        glm::vec3 neighborProjected = tangentPlane.projectPoint(neighbor);

        glm::vec2 samplingPoint {
                neighborProjected[0],
                neighborProjected[1]
        };
        float samplingValue = neighborProjected[2];

        samplingPoints[i] = samplingPoint;
        samplingValues[i] = samplingValue;
    }

    // Fit a bivariate function based on the sampling points and values
    CgMovingLeastSquares movingLeastSquares;
    movingLeastSquares.fitBivariate(samplingPoints, samplingValues, bivariateFunctionDegree);

    // Plot the fitted function
    // Draw a bounding box aligned with the axes of the tangent plane around the sampling points
    glm::vec2 min = { INFINITY, INFINITY };
    glm::vec2 max = { -INFINITY, -INFINITY };

    for (glm::vec2 samplingPoint : samplingPoints) {
        if (samplingPoint[0] < min[0]) {
            min[0] = samplingPoint[0];
        }
        if (samplingPoint[1] < min[1]) {
            min[1] = samplingPoint[1];
        }
        if (samplingPoint[0] > max[0]) {
            max[0] = samplingPoint[0];
        }
        if (samplingPoint[1] > max[1]) {
            max[1] = samplingPoint[1];
        }
    }

    // Step through the bounding box in the given amount of plotting steps
    // The step size is the size in each dimension divided by the step amount minus one,
    // because step 0 should capture the minimum axis values
    glm::vec2 stepSize;
    stepSize[0] = std::abs(max[0] - min[0]) / (plottingSteps - 1);
    stepSize[1] = std::abs(max[1] - min[1]) / (plottingSteps - 1);

    // Capture the plotted points as a 2-dimensional matrix
    glm::vec3 plottedPoints[plottingSteps][plottingSteps];

    for (size_t row = 0; row < plottingSteps; row++) {
        for (size_t col = 0; col < plottingSteps; col++) {

            // Calculate the sampling point for plotting
            glm::vec2 pointProjected;
            pointProjected[0] = min[0] + stepSize[0] * row;
            pointProjected[1] = min[1] + stepSize[1] * col;

            // Calculate the fitted value for the sampling point
            float fittedValue = movingLeastSquares.evaluateBivariate(pointProjected);
            glm::vec3 updatedPointProjected = {
                    pointProjected[0],
                    pointProjected[1],
                    fittedValue
            };

            // Unproject the point from tangent plane coordinates into global coordinates
            glm::vec3 plottedPoint = tangentPlane.unprojectPoint(updatedPointProjected);

            // Add the plotted point to the plotting matrix
            plottedPoints[row][col] = plottedPoint;

            //std::cout << "Projected:   " << row << ", " << col << ": "  << std::endl;
            // std::cout << "Unprojected: " << row << ", " << col << ": " << plottedPoint << std::endl;
        }
        //project the selected point to the plane
//        m_vertices[centerIndex] = smoothPoint(centerIndex, neighborCount, bivariateFunctionDegree);
//        m_vertex_colors[centerIndex] = {1,1,1};
    }

    // Create a triangle mesh based on the plotted points and color them red
    CgTriangleMesh* plottedPointMesh = new CgTriangleMesh();
    glm::vec3 color = { 0, 255, 0 };
    plottedPointMesh->addQuadrangleMesh((glm::vec3*)plottedPoints, plottingSteps, plottingSteps, tangentPlane.getNormal(), color);

    return plottedPointMesh;
}

/**
 * @brief CgPointCloud::smoothSurface Smooth the surface of this point cloud
 * @param neighborCount The amount of neighbors to base smoothing on
 * @param bivariateFunctionDegree The degree of the bivariate function to use for smoothing
 */
void CgPointCloud::smoothSurface(size_t neighborCount, size_t bivariateFunctionDegree)
{
    std::vector<glm::vec3> smoothedVertices(m_vertices.size());
    // Calculate all smoothed vertex positions
    for (size_t i = 0; i < m_vertices.size(); i++) {
        std::cout << "Smooth Vertex " << i << " Of "  << m_vertices.size() << std::endl;
        glm::vec3 updatedPosition = smoothPoint(i, neighborCount, bivariateFunctionDegree);
        if(updatedPosition != m_vertices[i]) {
            m_vertex_colors[i] = { 1.0, 0.0, 0.0};
        }
        smoothedVertices[i] = updatedPosition;
    }
    // Update all vertex positions to the smoothed version
    for (size_t i = 0; i < m_vertices.size(); i++) {
        m_vertices[i] = smoothedVertices[i];

    }
}
