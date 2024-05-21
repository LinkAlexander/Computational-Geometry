#ifndef CGKDTREE_H
#define CGKDTREE_H



#include <set>
#include <stddef.h>
#include <vector>

#include "CgPointCloud.h"
#include "CgTriangleMesh.h"
#include "LimitedSizePriorityQueue.h"

class CgKdTree;
class CgKdTreeNode;
class CgPointCloud;
class SplitPlane;

typedef std::pair<float, size_t> DistAndIndex;

class CgKdTree
{
    public:
        typedef enum axis3d {
            x = 0,
            y = 1,
            z = 2
        } axis3d;

        CgKdTree(CgPointCloud *pointCloud, int maxDepth);

        glm::vec3 getVertPos(size_t vertIndex);
        std::vector<size_t> getNearestNeighbors(size_t centerIndex, unsigned int neighborCount);
        size_t getClosestPointToRay(glm::vec3 rayStart, glm::vec3 rayDirection);
        std::vector<SplitPlane*> getSplitPlanes(size_t maxDepth);

    private:
        CgPointCloud *pointCloud;
        const std::vector<glm::vec3> &vertPositions;
        const std::vector<glm::vec3> &vertNormals;
        const std::vector<glm::vec3> &vertColors;

        CgKdTreeNode *root;
};

class CgKdTreeNode {
    public:
        CgKdTreeNode(CgKdTree *kdTree, CgKdTreeNode *parent, std::vector<size_t> &vertIndices, int depth, int maxDepth);
        LimitedSizePriorityQueue<DistAndIndex> getNearestNeighbors(size_t centerIndex, unsigned int neighborCount);
        std::vector<SplitPlane*> getSplitPlanes(size_t currentDepth, size_t maxDepth);

    private:
        CgKdTree::axis3d splitAxis;
        float splitNode;

        std::vector<size_t> vertIndices;

        CgKdTree *kdTree;

        CgKdTreeNode *parent;
        CgKdTreeNode *lesserChild; // left child
        CgKdTreeNode *greaterChild; // right child

        float distanceToSplitPlane(glm::vec3 pos);
        std::pair<glm::vec3, glm::vec3> getBoundingBox();
};

class SplitPlane {
    public:
        SplitPlane(glm::vec3 min, glm::vec3 max, CgKdTree::axis3d splitAxis, float splitValue);

        const glm::vec3 getPos1()const;
        const glm::vec3 getPos2()const;
        const glm::vec3 getPos3()const;
        const glm::vec3 getPos4()const;

        CgKdTree::axis3d getSplitAxis() const;
        float getSplitValue() const;

    private:
        glm::vec3 pos1, pos2, pos3, pos4;
        CgKdTree::axis3d splitAxis;
        float splitValue;
};

inline const glm::vec3 SplitPlane::getPos1() const
{
    return pos1;
}

inline const glm::vec3 SplitPlane::getPos2() const
{
    return pos2;
}

inline const glm::vec3 SplitPlane::getPos3() const
{
    return pos3;
}

inline const glm::vec3 SplitPlane::getPos4() const
{
    return pos4;
}

inline CgKdTree::axis3d SplitPlane::getSplitAxis() const
{
    return splitAxis;
}

inline float SplitPlane::getSplitValue() const
{
    return splitValue;
}

#endif // CGKDTREE_H
