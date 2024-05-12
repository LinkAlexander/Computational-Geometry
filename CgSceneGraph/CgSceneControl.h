#ifndef CGSCENECONTROL_H
#define CGSCENECONTROL_H

#include "CgBase/CgObserver.h"
#include "CgBase/CgBaseSceneControl.h"
#include <glm/glm.hpp>

class CgBaseEvent;
class CgBaseRenderer;
class CgPointCloud;
class CgTriangleMesh;
class CgPolyLine;
class CgHalfEdgeTriangleMesh;

class CgSceneControl : public CgObserver, public CgBaseSceneControl
{
public:
    CgSceneControl();
    ~CgSceneControl();

    void handleEvent(CgBaseEvent* e);
    void setRenderer(CgBaseRenderer* r);
    void renderObjects();

private:

    // since there is only one object, the select ray is in "bunny"-coordinates
    void calculatePickRay(double x, double y);
    bool m_show_pickray;
    CgPolyLine* m_select_ray;

    // example: how to calculate an Eigen decomosition of a 3x3 matrix
    void calculateEigenDecomposition3x3();

    // example: calculate a singular value decomposition an a Moore Penrose inverse of a (n,m)-matrix
    void calculateSingularValueDecomposition();

    // one single Disc used as Object to render the spats (in own local coordinates)
    CgTriangleMesh* m_disc;
    bool m_show_splats;

    // one simple HalfEdgeMEsh
    CgHalfEdgeTriangleMesh* m_he_mesh;

    // there is no scenegraph, either a point cloud or a mesh is used
    CgPointCloud* m_pointcloud;
    CgTriangleMesh* m_triangle_mesh;
    CgHalfEdgeTriangleMesh* m_HETriangle_mesh;

    // all the things that you would normally find in a scenegraph
    glm::mat4 m_current_transformation;
    glm::mat4 m_trackball_rotation;
    glm::mat4 m_lookAt_matrix;
    glm::mat4 m_proj_matrix;
    glm::vec3 m_center;

    // interface for OpenGl Rendering
    CgBaseRenderer* m_renderer;

    glm::vec3 pickRayStart;
    glm::vec3 pickRayDirection;
};

#endif // CGSCENECONTROL_H
