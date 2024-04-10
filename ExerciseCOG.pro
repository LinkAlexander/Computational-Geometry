QT       += core gui opengl widgets
TEMPLATE = app
TARGET = PointViewer
QMAKE_CXXFLAGS += -std=c++17
CONFIG += c++17

SOURCES += main.cpp \
    CgEvents/CgButtonPressedEvent.cpp \
    CgEvents/CgLoadHalfEdgeMeshEvent.cpp \
    CgEvents/CgLoadMeshEvent.cpp \
    CgEvents/CgLoadPointCloudEvent.cpp \
    CgEvents/CgPickRayEvent.cpp \
    CgEvents/CgSplatEvent.cpp \
    CgMath/CgEigenDecomposition3x3.cpp \
    CgQtViewer/CGQtGLRenderWidget.cpp \
    CgQtViewer/CgQtGui.cpp \
    CgBase/CgObservable.cpp \
    CgEvents/CgMouseEvent.cpp \
    CgQtViewer/CgQtMainApplication.cpp \
    CgSceneGraph/CgHalfEdgePrimitives.cpp \
    CgSceneGraph/CgHalfEdgeTriangleMesh.cpp \
    CgSceneGraph/CgPointCloud.cpp \
    CgSceneGraph/CgPolyLine.cpp \
    CgSceneGraph/CgSceneControl.cpp \
    CgEvents/CgKeyEvent.cpp \
    CgQtViewer/CgQtGlBufferObject.cpp \
    CgQtViewer/CgTrackball.cpp \
    CgEvents/CgWindowResizeEvent.cpp \
    CgSceneGraph/CgTriangleMesh.cpp \
    CgUtils/ObjLoader.cpp \
    CgEvents/CgTrackballEvent.cpp

HEADERS += \
    CgBase/CgBaseHalfEdgeTriangleMesh.h \
    CgBase/CgBaseHalfdgePrimitives.h \
    CgEvents/CgButtonPressedEvent.h \
    CgEvents/CgLoadHalfEdgeMeshEvent.h \
    CgEvents/CgLoadMeshEvent.h \
    CgEvents/CgLoadPointCloudEvent.h \
    CgEvents/CgPickRayEvent.h \
    CgEvents/CgSplatEvent.h \
    CgMath/CgEigenDecomposition3x3.h \
    CgMath/Eigen/Core \
    CgMath/Eigen/Eigen \
    CgMath/Eigen/SVD \
    CgQtViewer/CgQtGLRenderWidget.h \
    CgQtViewer/CgQtGui.h \
    CgBase/CgObserver.h \
    CgBase/CgObservable.h \
    CgBase/CgBaseEvent.h \
    CgBase/CgEnums.h \
    CgEvents/CgMouseEvent.h \
    CgQtViewer/CgQtMainApplication.h \
    CgSceneGraph/CgHalfEdgePrimitives.h \
    CgSceneGraph/CgHalfEdgeTriangleMesh.h \
    CgSceneGraph/CgPointCloud.h \
    CgSceneGraph/CgPolyLine.h \
    CgSceneGraph/CgSceneControl.h \
    CgEvents/CgKeyEvent.h \
    CgBase/CgBaseRenderer.h \
    CgBase/CgBaseRenderableObject.h \
    CgBase/CgBasePointCloud.h \
    CgBase/CgBaseTriangleMesh.h \
    CgBase/CgBasePolyline.h \
    CgBase/CgBaseSceneControl.h \
    CgQtViewer/CgQtGlBufferObject.h \
    CgQtViewer/CgTrackball.h \
    CgEvents/CgWindowResizeEvent.h \
    CgSceneGraph/CgTriangleMesh.h \
    CgUtils/CgEventEnums.h \
    CgUtils/ObjLoader.h \
    CgBase/CgBaseImage.h \
    CgEvents/CgTrackballEvent.h

