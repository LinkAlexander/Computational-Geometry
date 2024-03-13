#ifndef CGQTGLBUFFEROBJECT_H
#define CGQTGLBUFFEROBJECT_H


class CgBasePolyline;
class CgBaseTriangleMesh;
class CgBaseHalfEdgeTriangleMesh;
class CgBasePointCloud;
class CgBaseRenderableObject;

#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLShaderProgram>

class CgQtGlBufferObject
{
public:
    CgQtGlBufferObject();
    CgQtGlBufferObject(QOpenGLShaderProgram *);


    QOpenGLBuffer vertexbuffer;
    QOpenGLBuffer colorbuffer;
    QOpenGLBuffer normalbuffer;
    QOpenGLBuffer indexbuffer;

    QOpenGLShaderProgram *m_program;

    void draw(CgBaseRenderableObject*);

    void initPolyline(CgBasePolyline* obj);
    void initTriangleMesh(CgBaseTriangleMesh* obj);
    void initTriangleMeshFlatShading(CgBaseTriangleMesh* obj);
    void initHalfEdgeTriangleMesh(CgBaseHalfEdgeTriangleMesh* obj);
    void initPointCloud(CgBasePointCloud* obj);

private:

    void drawPolyline(CgBasePolyline* obj);
    void drawTriangleMesh(CgBaseTriangleMesh* obj);
    void drawHalfEdgeTriangleMesh(CgBaseHalfEdgeTriangleMesh* obj);
    void drawPointCloud(CgBasePointCloud* obj);

};

#endif // CGQTGLBUFFEROBJECT_H
