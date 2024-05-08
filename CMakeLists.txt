cmake_minimum_required(VERSION 3.27)
project(PointViewer)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)


find_package(Qt5 COMPONENTS
        Core
        Gui
        Widgets
        OpenGL
        REQUIRED)
find_package(OpenGL REQUIRED)

set(SOURCES
        CgBase/CgBaseEvent.h
        CgBase/CgBaseHalfEdgeTriangleMesh.h
        CgBase/CgBaseHalfdgePrimitives.h
        CgBase/CgBaseImage.h
        CgBase/CgBasePointCloud.h
        CgBase/CgBasePolyline.h
        CgBase/CgBaseRenderableObject.h
        CgBase/CgBaseRenderer.h
        CgBase/CgBaseSceneControl.h
        CgBase/CgBaseTriangleMesh.h
        CgBase/CgEnums.h
        CgBase/CgObservable.cpp CgBase/CgObservable.h
        CgBase/CgObserver.h
        CgEvents/CgKeyEvent.cpp CgEvents/CgKeyEvent.h
        CgEvents/CgLoadHalfEdgeMeshEvent.cpp CgEvents/CgLoadHalfEdgeMeshEvent.h
        CgEvents/CgLoadMeshEvent.cpp CgEvents/CgLoadMeshEvent.h
        CgEvents/CgLoadPointCloudEvent.cpp CgEvents/CgLoadPointCloudEvent.h
        CgEvents/CgMouseEvent.cpp CgEvents/CgMouseEvent.h
        CgEvents/CgPickRayEvent.cpp CgEvents/CgPickRayEvent.h
        CgEvents/CgSplatEvent.cpp CgEvents/CgSplatEvent.h
        CgEvents/CgTrackballEvent.cpp CgEvents/CgTrackballEvent.h
        CgEvents/CgWindowResizeEvent.cpp CgEvents/CgWindowResizeEvent.h
        CgMath/CgEigenDecomposition3x3.cpp CgMath/CgEigenDecomposition3x3.h
        CgMath/Eigen/Core
        CgMath/Eigen/Eigen
        CgMath/Eigen/SVD
        CgQtViewer/CGQtGLRenderWidget.cpp
        CgQtViewer/CgQtGLRenderWidget.h
        CgQtViewer/CgQtGlBufferObject.cpp CgQtViewer/CgQtGlBufferObject.h
        CgQtViewer/CgQtGui.cpp CgQtViewer/CgQtGui.h
        CgQtViewer/CgQtMainApplication.cpp CgQtViewer/CgQtMainApplication.h
        CgQtViewer/CgTrackball.cpp CgQtViewer/CgTrackball.h
        CgSceneGraph/CgHalfEdgePrimitives.cpp CgSceneGraph/CgHalfEdgePrimitives.h
        CgSceneGraph/CgHalfEdgeTriangleMesh.cpp CgSceneGraph/CgHalfEdgeTriangleMesh.h
        CgSceneGraph/CgPointCloud.cpp CgSceneGraph/CgPointCloud.h
        CgSceneGraph/CgPolyLine.cpp CgSceneGraph/CgPolyLine.h
        CgSceneGraph/CgSceneControl.cpp CgSceneGraph/CgSceneControl.h
        CgSceneGraph/CgTriangleMesh.cpp CgSceneGraph/CgTriangleMesh.h
        CgUtils/ObjLoader.cpp CgUtils/ObjLoader.h
        main.cpp
)



include_directories(.) #  or include_directory(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(CgEvents)
include_directories(CgMath)
include_directories(CgQtViewer)
include_directories(CgBase)
include_directories(CgSceneGraph)
include_directories(CgUtils)
include_directories(CgMath/Eigen)

add_executable(PointViewer ${SOURCES})

target_link_libraries(PointViewer
        Qt5::Core
        Qt5::Gui
        Qt5::Widgets
        Qt5::OpenGL
        ${OPENGL_LIBRARIES}
)

