#ifndef CGBASEPOINTCLOUD_H
#define CGBASEPOINTCLOUD_H

#include <vector>
#include <glm/glm.hpp>
#include "CgBase/CgBaseRenderableObject.h"

class CgBasePointCloud : public CgBaseRenderableObject
{
public:

    virtual ~CgBasePointCloud(){};

    // vertex positions in local coordinates
    virtual const std::vector<glm::vec3>& getVertices()     const =0;

    // normal directions in local coordinates (normalized)
    virtual const std::vector<glm::vec3>& getVertexNormals() const =0;

    // vertex colors
    virtual const std::vector<glm::vec3>& getVertexColors() const =0;

    // a local "look at matrix" for normal direction and the spanning directions of the tangent plane
    virtual const std::vector<glm::mat4>& getSplatOrientations() const =0;

    // for rendering ellipes use a different scaling in the two spanning directions of the tangent plane
    virtual const std::vector<glm::vec2>& getSplatScalings() const =0;



};

#endif // CGBASEPOINTCLOUD_H
