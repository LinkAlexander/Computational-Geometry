#define GLM_ENABLE_EXPERIMENTAL

#include "CgSceneControl.h"
#include "CgBase/CgEnums.h"
#include "CgEvents/CgMouseEvent.h"
#include "CgEvents/CgKeyEvent.h"
#include "CgEvents/CgButtonPressedEvent.h"
#include "CgUtils/CgEventEnums.h"
#include "CgEvents/CgWindowResizeEvent.h"
#include "CgEvents/CgLoadMeshEvent.h"
#include "CgEvents/CgLoadHalfEdgeMeshEvent.h"
#include "CgEvents/CgLoadPointCloudEvent.h"
#include "CgEvents/CgTrackballEvent.h"
#include "CgEvents/CgSplatEvent.h"
#include "CgEvents/CgPickRayEvent.h"
#include "CgBase/CgBaseRenderer.h"
#include "CgPointCloud.h"
#include "CgPolyLine.h"
#include "CgTriangleMesh.h"
#include "CgHalfEdgeTriangleMesh.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>
#include "CgUtils/ObjLoader.h"
#include "CgMath/CgEigenDecomposition3x3.h"
#include <string>
#include <CgMath/Eigen/SVD>
#include <CgMath/Eigen/Core>
using namespace Eigen;

CgSceneControl::CgSceneControl()
{
    m_pointcloud=nullptr;
    m_current_transformation=glm::mat4(1.);
    m_lookAt_matrix= glm::lookAt(glm::vec3(0.0,0.0,1.0),glm::vec3(0.0,0.0,0.0),glm::vec3(0.0,1.0,0.0));
    m_proj_matrix= glm::mat4x4(glm::vec4(1.792591, 0.0, 0.0, 0.0), glm::vec4(0.0, 1.792591, 0.0, 0.0), glm::vec4(0.0, 0.0, -1.0002, -1.0), glm::vec4(0.0, 0.0, -0.020002, 0.0));

    m_trackball_rotation=glm::mat4(1.);

    m_select_ray=NULL;
    m_pointcloud = NULL;
    m_disc = NULL;
    m_he_mesh= NULL;

    m_he_mesh = new CgHalfEdgeTriangleMesh();
    m_center = m_he_mesh->getCenter();

    m_triangle_mesh= NULL;
    m_show_splats=false;
    m_show_pickray=false;
}


CgSceneControl::~CgSceneControl()
{
    if(m_pointcloud!=NULL)
        delete m_pointcloud;

    if(m_he_mesh!=NULL)
        delete m_he_mesh;

    if(m_triangle_mesh!=NULL)
        delete m_triangle_mesh;
    if(m_disc!=NULL)
        delete m_disc;
    if(m_select_ray!=NULL)
        delete m_select_ray;
}



void CgSceneControl::calculatePickRay(double x, double y)
{
    double w = (double)m_renderer->getViewportWidth();
    double h = (double)m_renderer->getViewportHeight();

    // normalize into [-1;1]
    x = 2.0 * x / w - 1.0;
    if (x < -1.0)
        x = -1.0;
    if (x > 1.0)
        x = 1.0;

    y = 2.0 * y / h - 1.0;
    if (y < -1.0)
        y = -1.0;
    if (y > 1.0)
        y = 1.0;

    // change to right handed coordinate system
    y = -y;

    glm::mat4 inverse_proj = glm::inverse(m_proj_matrix);

    //unproject point on front clipping plane
    glm::vec4 p(x, y, -0.01, 1);
    glm::vec4 q = inverse_proj * p;
    q /= q.w;

    //unproject point on back clipping plane
    glm::vec4 r(x, y, 1.0, 1);
    glm::vec4 s = inverse_proj * r;
    s /= s.w;

    //construct current modelview matrix by hand, since there is no scenegraph
    m_current_transformation = glm::translate(m_current_transformation, -m_center);
    glm::mat4 mv_matrix = m_lookAt_matrix * m_trackball_rotation * m_current_transformation;
    m_current_transformation = glm::translate(m_current_transformation, m_center);

    // convert pick ray into local "bunny coordinates"
    glm::mat4 mv_matrix_inv = glm::inverse(mv_matrix);
    glm::vec4 raystart = mv_matrix_inv * q;
    glm::vec4 rayend = mv_matrix_inv * s;

    // init new CgPolyline to draw the pick ray
    if (m_select_ray != NULL)
        delete m_select_ray;
    std::vector<glm::vec3> pointlist;
    pointlist.push_back(raystart);
    pointlist.push_back(rayend);
    m_select_ray = new CgPolyLine(33, pointlist);
    m_renderer->init(m_select_ray);

    pickRayStart = raystart;
    pickRayDirection = rayend - raystart;

    // Apply the pick ray to point cloud
    if (m_pointcloud) {
        glm::vec3 pickRayStart = raystart;
        glm::vec3 pickRayDirection = rayend - raystart;
        m_pointcloud->applyPickRay(pickRayStart, pickRayDirection);
        m_renderer->init(m_pointcloud);
    }

    m_renderer->redraw();
}


void CgSceneControl::calculateEigenDecomposition3x3()
{
    glm::mat3 testmat = glm::mat3(3.,2.,4.,2.,0.,2.,4.,2.,3.);
    std::cout << glm::to_string(testmat) << std::endl;

    CgEigenDecomposition3x3 eigen(testmat);
    glm::mat3 eigenvectors= eigen.getEigenvectors();
    glm::vec3 eigenvalues = eigen.getEigenvalues();

    std::cout << glm::to_string(eigenvalues) << std::endl;
    std::cout << glm::to_string(eigenvectors) << std::endl;
}

void CgSceneControl::calculateSingularValueDecomposition()
{
      //init some arbitrary test matrix
      MatrixXd C(2,3);
      C(0,0)=3.0;C(0,1)=2.0;C(0,2)=2.0;
      C(1,0)=2.0;C(1,1)=3.0;C(1,2)=-2.0;

      std::cout << C << std::endl;

      //compute decomposition
      JacobiSVD<Eigen::MatrixXd> svd(C,ComputeThinU|ComputeThinV);
      Eigen::MatrixXd U=svd.matrixU();
      Eigen::MatrixXd V=svd.matrixV();
      Eigen::VectorXd SV=svd.singularValues();

      //build a diagonal matrix out of the singular values
      Eigen::MatrixXd S(2,2);
      S.setZero();
      S(0,0)=SV(0);S(1,1)=SV(1);

      std::cout << U << std::endl;
      std::cout << S << std::endl;
      std::cout << V.transpose() << std::endl;


    // compute Moore-Penrose inverse now

     S(0,0)=1.0/S(0,0);S(1,1)=1.0/S(1,1);
     Eigen::MatrixXd C_plus = V*S*U.transpose();

     std::cout << std::endl;
     std::cout << C_plus << std::endl;
}

void CgSceneControl::handleEvent(CgBaseEvent* e)
{
    // any kind of MouseEvents can be checked via CgEvent::CgMouseEvent a
    // see CgEvent enums in CgEnums.h

    if(e->getType() & Cg::CgMouseEvent)
    {
        CgMouseEvent* ev = (CgMouseEvent*)e;
        if(ev->button()==Cg::RightButton)
            calculatePickRay((double)ev->x(),(double)ev->y());
    }

    if(e->getType() & Cg::CgTrackballEvent)
    {
        CgTrackballEvent* ev = (CgTrackballEvent*)e;

        m_trackball_rotation=ev->getRotationMatrix();
        m_renderer->redraw();
    }

    // any kind of KeyEvents can be checked via CgEvent::CgKeyEvent
    // see CgEvent enums in CgEnums.h

    if(e->getType() & Cg::CgKeyEvent)
    {
        CgKeyEvent* ev = (CgKeyEvent*)e;

        if(ev->text()=="e")
        {
            // example usage of eigen library to compute the eigen-decomposition of a matrix
            calculateEigenDecomposition3x3();
        }

        if(ev->text()=="s")
        {
          // example usage of eigen library to compute the eigen-decomposition of a matrix
          calculateSingularValueDecomposition();
        }

        if(ev->text()=="+")
        {
            //scaling of scene, i.e. of the one single object
            glm::mat4 scalemat = glm::mat4(1.);
            scalemat = glm::scale(scalemat,glm::vec3(1.1,1.1,1.1));
            m_current_transformation=m_current_transformation*scalemat;
            m_renderer->redraw();
        }
        if(ev->text()=="-")
        {
            //scaling of scene, i.e. of the one single object
            glm::mat4 scalemat = glm::mat4(1.);
            scalemat = glm::scale(scalemat,glm::vec3(0.9,0.9,0.9));
            m_current_transformation=m_current_transformation*scalemat;
            m_renderer->redraw();
        }
    }


    if(e->getType() & Cg::CgSplatEvent)
    {
        // checkbox in GUI changed for showing/hiding splats
        CgSplatEvent* ev = (CgSplatEvent*)e;
        m_show_splats=ev->useSplatting();
        m_renderer->redraw();
    }

    if(e->getType() & Cg::CgPickRayEvent)
    {
        // checkbox in GUI changed for showing/hiding pick ray
        CgPickRayEvent* ev = (CgPickRayEvent*)e;
        m_show_pickray=ev->showPickRay();
        renderObjects();
        m_renderer->redraw();
    }


    if(e->getType() & Cg::WindowResizeEvent)
    {
         CgWindowResizeEvent* ev = (CgWindowResizeEvent*)e;
         std::cout << *ev <<std::endl;
         m_proj_matrix=glm::perspective(45.0f, (float)(ev->w()) / ev->h(), 0.01f, 100.0f);
    }

    if(e->getType() & Cg::CgLoadPointCloudEvent)
    {

        CgLoadPointCloudEvent* ev = (CgLoadPointCloudEvent*)e;
        ObjLoader* loader= new ObjLoader();

        std::string filename =ev->FileName();
        if(filename !="")
          {
            loader->load(filename);

            if(m_triangle_mesh!=NULL)
              {
                delete m_triangle_mesh;
                m_triangle_mesh=NULL;
              }
            if(m_pointcloud!=NULL)
              {
                delete m_pointcloud;
                m_pointcloud=NULL;
              }

            if(m_he_mesh!=NULL)
              {
                delete m_he_mesh;
                m_he_mesh=NULL;
              }


            m_pointcloud= new CgPointCloud();
            m_pointcloud->init(ev->FileName(),true);

            m_center=m_pointcloud->getCenter();
            m_renderer->init(m_pointcloud);
            m_renderer->redraw();
          }
    }


    if(e->getType() & Cg::CgLoadMeshEvent)
    {

        CgLoadMeshEvent* ev = (CgLoadMeshEvent*)e;
        ObjLoader* loader= new ObjLoader();
        std::string filename =ev->FileName();

        if(filename !="")
          {
            loader->load(filename);

            if(m_triangle_mesh!=NULL)
              {
                delete m_triangle_mesh;
                m_triangle_mesh=NULL;
              }
            if(m_pointcloud!=NULL)
              {
                delete m_pointcloud;
                m_pointcloud=NULL;
              }

            if(m_he_mesh!=NULL)
              {
                delete m_he_mesh;
                m_he_mesh=NULL;
              }


            m_triangle_mesh= new CgTriangleMesh();
            m_triangle_mesh->init(ev->FileName());
            m_center=m_triangle_mesh->getCenter();
            m_renderer->init(m_triangle_mesh);
            m_renderer->redraw();
          }
    }


    if(e->getType() & Cg::CgLoadHalfEdgeMeshEvent)
    {

        CgLoadHalfEdgeMeshEvent* ev = (CgLoadHalfEdgeMeshEvent*)e;
        ObjLoader* loader= new ObjLoader();
        std::string filename =ev->FileName();

        if(filename !="")
          {
            loader->load(filename);

            if(m_triangle_mesh!=NULL)
              {
                delete m_triangle_mesh;
                m_triangle_mesh=NULL;
              }
            if(m_pointcloud!=NULL)
              {
                delete m_pointcloud;
                m_pointcloud=NULL;
              }

            if(m_he_mesh!=NULL)
              {
                delete m_he_mesh;
                m_he_mesh=NULL;
              }

            m_he_mesh= new CgHalfEdgeTriangleMesh();
            m_he_mesh->init(ev->FileName());

            m_center=m_he_mesh->getCenter();
            m_renderer->init(m_he_mesh);
            m_renderer->redraw();
          }
    }
    // Subdivision
    if (e->getType() & Cg::CgButtonPressedEvent) {
        CgButtonPressedEvent* ev = (CgButtonPressedEvent*)e;
        int code = ev->getCode();

        switch (code) {
        case SUBDIVISION:
            m_he_mesh->loopSubdivision();
            m_renderer->init(m_he_mesh);
            m_renderer->redraw();
            break;
            /*case RESET:
            delete m_he_mesh;
            m_renderer->init(m_he_mesh);
            m_renderer->redraw();
            break;*/
        }
    }




    // delete event
    delete e;


}





// rendering stuff
/*******************************************************************************************************/

void CgSceneControl::setRenderer(CgBaseRenderer* r)
{
    m_renderer=r;
    m_renderer->setSceneControl(this);

    if(m_pointcloud!=NULL)
    m_renderer->init(m_pointcloud);

    if(m_disc!=NULL)
    m_renderer->init(m_disc);

    if(m_select_ray!=NULL)
    m_renderer->init(m_select_ray);

    if(m_triangle_mesh!=NULL)
    m_renderer->init(m_triangle_mesh);

    if(m_he_mesh!=NULL)
    m_renderer->init(m_he_mesh);
}


void CgSceneControl::renderObjects()
{

    // Materialeigenschaften setzen


    m_renderer->setUniformValue("mycolor",glm::vec4(0.0,1.0,0.0,1.0));

    m_renderer->setUniformValue("matDiffuseColor",glm::vec4(0.780392f, 0.568627f, 0.113725f, 1.0f ));
    m_renderer->setUniformValue("lightDiffuseColor",glm::vec4(1.0,1.0,1.0,1.0));

    m_renderer->setUniformValue("matAmbientColor",glm::vec4(0.329412f, 0.223529f, 0.027451f,1.0f ));
    m_renderer->setUniformValue("lightAmbientColor",glm::vec4(1.0,1.0,1.0,1.0));

    m_renderer->setUniformValue("matSpecularColor",glm::vec4(0.992157f, 0.941176f, 0.807843f, 1.0f));
    m_renderer->setUniformValue("lightSpecularColor",glm::vec4(1.0,1.0,1.0,1.0));

    m_renderer->setUniformValue("shininess",27.8974f);



    if(m_he_mesh!=NULL)
      {
        m_current_transformation = glm::translate(m_current_transformation,-m_center);
        glm::mat4 mv_matrix = m_lookAt_matrix * m_trackball_rotation* m_current_transformation ;
        m_renderer->render(m_he_mesh,mv_matrix,m_proj_matrix);
        m_current_transformation = glm::translate(m_current_transformation,m_center);
      }





    if(m_triangle_mesh!=NULL)
      {
        m_current_transformation = glm::translate(m_current_transformation,-m_center);
        glm::mat4 mv_matrix = m_lookAt_matrix * m_trackball_rotation* m_current_transformation ;
        m_renderer->render(m_triangle_mesh,mv_matrix,m_proj_matrix);
        m_current_transformation = glm::translate(m_current_transformation,m_center);
      }

    if((m_pointcloud!=NULL)&&(!m_show_splats))
      {
        m_current_transformation = glm::translate(m_current_transformation,-m_center);
        glm::mat4 mv_matrix = m_lookAt_matrix * m_trackball_rotation* m_current_transformation ;
        m_renderer->render(m_pointcloud,mv_matrix,m_proj_matrix);
        m_current_transformation = glm::translate(m_current_transformation,m_center);
      }




    if((m_select_ray!=NULL)&&(m_show_pickray))
      {
        m_current_transformation = glm::translate(m_current_transformation,-m_center);
        glm::mat4 mv_matrix = m_lookAt_matrix * m_trackball_rotation* m_current_transformation ;
        m_renderer->render(m_select_ray,mv_matrix,m_proj_matrix);
        m_current_transformation = glm::translate(m_current_transformation,m_center);
      }


    // you would never ever render splats like this, a splat is never a scenegraph object
    // just to have not much effort to show the shape, inefficient rendering!
    if((m_pointcloud!=NULL)&&(m_show_splats))
      {

        m_renderer->setUniformValue("rendersplats",1);

    const std::vector<glm::vec3>& vertex_colors = m_pointcloud->getVertexColors();
    const std::vector<glm::vec3>& vertices = m_pointcloud->getVertices();
    const std::vector<unsigned int>& splat_indices= m_pointcloud->getSplatIndices();
    const std::vector<glm::mat4>& orientations = m_pointcloud->getSplatOrientations();
    const std::vector<glm::vec2>& scalings = m_pointcloud->getSplatScalings();

    glm::mat4 mv_matrix;

      for(unsigned int i=0;i<splat_indices.size();i++)
      {
        //m_renderer->setUniformValue("mycolor",glm::vec4(vertex_colors[splat_indices[i]],1.0));

        m_current_transformation = glm::translate(m_current_transformation,-m_center);
        m_current_transformation = m_current_transformation * glm::inverse(orientations[splat_indices[i]]);

        m_current_transformation=glm::scale(m_current_transformation,glm::vec3(scalings[splat_indices[i]][0],scalings[splat_indices[i]][1],1.0));

        mv_matrix = m_lookAt_matrix * m_trackball_rotation* m_current_transformation ;

        if(m_disc!=NULL)
             m_renderer->render(m_disc,mv_matrix,m_proj_matrix);

        m_current_transformation=glm::scale(m_current_transformation,glm::vec3(1.0/scalings[splat_indices[i]][0],1.0/scalings[splat_indices[i]][1],1.0));

        m_current_transformation = m_current_transformation * orientations[splat_indices[i]];
        m_current_transformation = glm::translate(m_current_transformation,m_center);

        }
    }

}

