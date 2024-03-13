#include "CgQtGlBufferObject.h"
#include "CgBase/CgBasePolyline.h"
#include "CgBase/CgBaseHalfEdgeTriangleMesh.h"
#include "CgBase/CgBaseTriangleMesh.h"
#include "CgBase/CgBasePointCloud.h"

CgQtGlBufferObject::CgQtGlBufferObject():
 indexbuffer(QOpenGLBuffer::IndexBuffer)
{
       vertexbuffer.create();
       normalbuffer.create();
       colorbuffer.create();
       indexbuffer.create();
}


CgQtGlBufferObject::CgQtGlBufferObject(QOpenGLShaderProgram* program):
 indexbuffer(QOpenGLBuffer::IndexBuffer)
{
    m_program = program;
    vertexbuffer.create();
     colorbuffer.create();
    normalbuffer.create();
    indexbuffer.create();
}






void CgQtGlBufferObject::initPolyline(CgBasePolyline* obj)
{

    vertexbuffer.bind();
    vertexbuffer.setUsagePattern(QOpenGLBuffer::StaticDraw);
    vertexbuffer.allocate(obj->getVertices().data(),  obj->getVertices().size()* sizeof(glm::vec3));

    if(obj->getVertexColors().size()>0)
      {
    colorbuffer.bind();
    colorbuffer.allocate(obj->getVertexColors().data(),  obj->getVertexColors().size()* sizeof(glm::vec3));
      }


}



void CgQtGlBufferObject::initTriangleMesh(CgBaseTriangleMesh* obj)
{


    vertexbuffer.bind();
    vertexbuffer.allocate(obj->getVertices().data(),  obj->getVertices().size()* sizeof(glm::vec3));

    normalbuffer.bind();
    normalbuffer.allocate(obj->getVertexNormals().data(),  obj->getVertexNormals().size()* sizeof(glm::vec3));

    if(obj->getVertexColors().size()>0)
      {
    colorbuffer.bind();
    colorbuffer.allocate(obj->getVertexColors().data(),  obj->getVertexColors().size()* sizeof(glm::vec3));
      }


    indexbuffer.bind();
    indexbuffer.allocate(obj->getTriangleIndices().data(),  obj->getTriangleIndices().size()* sizeof(unsigned int));

}

void CgQtGlBufferObject::initTriangleMeshFlatShading(CgBaseTriangleMesh* obj)
{
    std::vector<glm::vec3> verts = obj->getVertices();
    std::vector<unsigned int> indices=obj->getTriangleIndices();
    std::vector<glm::vec3> normals = obj->getFaceNormals();
    std::vector<glm::vec3> colors = obj->getVertexColors();

    std::vector<glm::vec3> multiple_verts;
    std::vector<glm::vec3> multiple_colors;
    std::vector<glm::vec3> multiple_normals;
    std::vector<unsigned int> multiple_indices;

    for(unsigned int i=0;i<indices.size();i++)
    {
        multiple_verts.push_back(verts[indices[i]]);
        multiple_colors.push_back(colors[indices[i]]);
        multiple_normals.push_back(normals[i/3]);
        multiple_indices.push_back(i);
    }

    vertexbuffer.bind();
    vertexbuffer.allocate(multiple_verts.data(),  multiple_verts.size()* sizeof(glm::vec3));

    normalbuffer.bind();
    normalbuffer.allocate(multiple_normals.data(),  multiple_normals.size()* sizeof(glm::vec3));

    if(obj->getVertexColors().size()>0)
      {
    colorbuffer.bind();
    colorbuffer.allocate(multiple_colors.data(),  multiple_colors.size()* sizeof(glm::vec3));
      }


    indexbuffer.bind();
    indexbuffer.allocate(multiple_indices.data(),  multiple_indices.size()* sizeof(unsigned int));

}


void CgQtGlBufferObject::initHalfEdgeTriangleMesh(CgBaseHalfEdgeTriangleMesh* obj)
{
    std::vector<CgBaseHeFace*> faces = obj->getFaces();

    std::vector<glm::vec3> multiple_verts;
    std::vector<glm::vec3> multiple_colors;
    std::vector<glm::vec3> multiple_normals;
    std::vector<unsigned int> multiple_indices;

    for(unsigned int i=0;i<faces.size();i++)
    {
        multiple_verts.push_back(faces[i]->edge()->vert()->position());
        multiple_colors.push_back(faces[i]->edge()->vert()->color());
        multiple_indices.push_back(3*i);

        multiple_verts.push_back(faces[i]->edge()->next()->vert()->position());
        multiple_colors.push_back(faces[i]->edge()->next()->vert()->color());
        multiple_indices.push_back(3*i+1);

        multiple_verts.push_back(faces[i]->edge()->next()->next()->vert()->position());
        multiple_colors.push_back(faces[i]->edge()->next()->next()->vert()->color());
        multiple_indices.push_back(3*i+2);

        multiple_normals.push_back(faces[i]->normal());
        multiple_normals.push_back(faces[i]->normal());
        multiple_normals.push_back(faces[i]->normal());
    }

    vertexbuffer.bind();
    vertexbuffer.allocate(multiple_verts.data(),  multiple_verts.size()* sizeof(glm::vec3));

    normalbuffer.bind();
    normalbuffer.allocate(multiple_normals.data(),  multiple_normals.size()* sizeof(glm::vec3));

    colorbuffer.bind();
    colorbuffer.allocate(multiple_colors.data(),  multiple_colors.size()* sizeof(glm::vec3));

    indexbuffer.bind();
    indexbuffer.allocate(multiple_indices.data(),  multiple_indices.size()* sizeof(unsigned int));

}

void CgQtGlBufferObject::initPointCloud(CgBasePointCloud* obj)
{
  vertexbuffer.bind();
  vertexbuffer.setUsagePattern(QOpenGLBuffer::StaticDraw);
  vertexbuffer.allocate(obj->getVertices().data(),  obj->getVertices().size()* sizeof(glm::vec3));

  if(obj->getVertexNormals().size()>0)
    {
  normalbuffer.bind();
  normalbuffer.allocate(obj->getVertexNormals().data(),  obj->getVertexNormals().size()* sizeof(glm::vec3));
    }

  if(obj->getVertexColors().size()>0)
    {
  colorbuffer.bind();
  colorbuffer.allocate(obj->getVertexColors().data(),  obj->getVertexColors().size()* sizeof(glm::vec3));
    }


}

void CgQtGlBufferObject::draw(CgBaseRenderableObject* obj)
{

    if(obj->getType()==Cg::TriangleMesh)
        drawTriangleMesh((CgBaseTriangleMesh*) obj);
    if(obj->getType()==Cg::HalfEdgeTriangleMesh)
        drawHalfEdgeTriangleMesh((CgBaseHalfEdgeTriangleMesh*) obj);
    if(obj->getType()==Cg::Polyline)
        drawPolyline((CgBasePolyline*)obj);
    if(obj->getType()==Cg::PointCloud)
        drawPointCloud((CgBasePointCloud*) obj);
}

void CgQtGlBufferObject::drawTriangleMesh(CgBaseTriangleMesh* obj)
{

    m_program->bind();

    if(obj->getVertices().size()!=0)
    {
        vertexbuffer.bind();
        int vertexlocation = m_program->attributeLocation("vertex");
        m_program->enableAttributeArray(vertexlocation);
        m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));
    }

    if(obj->getVertexColors().size()!=0)
    {
        colorbuffer.bind();
        int colorlocation = m_program->attributeLocation("color");
        m_program->enableAttributeArray(colorlocation);
        m_program->setAttributeBuffer(colorlocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));
    }

    if(obj->getVertexNormals().size()!=0)
    {
        normalbuffer.bind();
        int normallocation = m_program->attributeLocation("normal");
        m_program->enableAttributeArray(normallocation);
        m_program->setAttributeBuffer(normallocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));
    }

    if(obj->getTriangleIndices().size()!=0)
    {
        indexbuffer.bind();
        glDrawElements(GL_TRIANGLES,obj->getTriangleIndices().size(), GL_UNSIGNED_INT,0);
    }
}




void CgQtGlBufferObject::drawPolyline(CgBasePolyline* obj)
{
    vertexbuffer.bind();
    int vertexlocation = m_program->attributeLocation("vertex");
    m_program->enableAttributeArray(vertexlocation);
    m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));

    if(obj->getVertexColors().size()!=0)
    {
        colorbuffer.bind();
        int colorlocation = m_program->attributeLocation("color");
        m_program->enableAttributeArray(colorlocation);
        m_program->setAttributeBuffer(colorlocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));
    }

    glLineWidth((GLfloat) obj->getLineWidth());
    glDrawArrays(GL_LINE_STRIP,0,obj->getVertices().size());
}


void CgQtGlBufferObject::drawHalfEdgeTriangleMesh(CgBaseHalfEdgeTriangleMesh* obj)
{
    m_program->bind();

        vertexbuffer.bind();
        int vertexlocation = m_program->attributeLocation("vertex");
        m_program->enableAttributeArray(vertexlocation);
        m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));

        colorbuffer.bind();
        int colorlocation = m_program->attributeLocation("color");
        m_program->enableAttributeArray(colorlocation);
        m_program->setAttributeBuffer(colorlocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));



        normalbuffer.bind();
        int normallocation = m_program->attributeLocation("normal");
        m_program->enableAttributeArray(normallocation);
        m_program->setAttributeBuffer(normallocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));


    if(obj->getFaces().size()!=0)
    {
        indexbuffer.bind();
        glDrawElements(GL_TRIANGLES,obj->getFaces().size()*3, GL_UNSIGNED_INT,0);
    }
}

void CgQtGlBufferObject::drawPointCloud(CgBasePointCloud* obj)
{
    if(obj->getVertices().size()!=0)
    {
        vertexbuffer.bind();
        int vertexlocation = m_program->attributeLocation("vertex");
        m_program->enableAttributeArray(vertexlocation);
        m_program->setAttributeBuffer(vertexlocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));


        if(obj->getVertexColors().size()!=0)
        {
            colorbuffer.bind();
            int colorlocation = m_program->attributeLocation("color");
            m_program->enableAttributeArray(colorlocation);
            m_program->setAttributeBuffer(colorlocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));
        }


        if(obj->getVertexNormals().size()!=0)
        {
            normalbuffer.bind();
            int normallocation = m_program->attributeLocation("normal");
            m_program->enableAttributeArray(normallocation);
            m_program->setAttributeBuffer(normallocation, GL_FLOAT, 0, 3 , sizeof(glm::vec3));
        }

        glDrawArrays(GL_POINTS,0,obj->getVertices().size());
    }
}









