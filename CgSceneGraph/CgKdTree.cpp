#include "CgKdTree.h"

#include <numeric>

#include <QDebug>

template class LimitedSizePriorityQueue<DistAndIndex>;

/* KdTree */

/**
 * @brief CgKdTree::CgKdTree KdTree constructor
 * @param pointCloud the point cloud to be processed
 * @param maxDepth the maximum depth of the KdTree
 */
CgKdTree::CgKdTree(CgPointCloud* pointCloud, int maxDepth)
    : pointCloud(pointCloud)
    , vertPositions(pointCloud->getVertices())
    , vertNormals(pointCloud->getVertexNormals())
    , vertColors(pointCloud->getVertexColors())
{
    // init the indices of the root
    std::vector<size_t> vertIndices(vertPositions.size());
    std::iota(vertIndices.begin(), vertIndices.end(), 0); //fill vector with values
    root = new CgKdTreeNode(this, NULL, vertIndices, 0, maxDepth);
}

/**
 * @brief CgKdTree::getVertPos Returns the position of the wanted vert.
 * @param vertIndex Index of the wanted vert
 * @return the position of the wanted vert
 */
glm::vec3 CgKdTree::getVertPos(size_t vertIndex)
{
    return vertPositions[vertIndex];
}

/**
 * @brief CgKdTree::getNearestNeighbor Returns all indices of the nearest neighbors of one point.
 * @param centerIndex the index of the point for which we want to get the nearest neighbors
 * @param neighborCount the number of nearest neighbors for which the indices are to be returned
 * @return a vector with all indices of the nearest neighbors
 */
std::vector<size_t> CgKdTree::getNearestNeighbors(size_t centerIndex, unsigned int neighborCount)
{
    // get neighbors starting at the root node (recursive)
    LimitedSizePriorityQueue<DistAndIndex> neighbors = root->getNearestNeighbors(centerIndex, neighborCount);

    std::vector<size_t> neighborIndices(neighborCount);

    size_t i = 0;

    // write results into a vector (sorted->smallest to largest)
    while (!neighbors.empty()) {
        DistAndIndex neighbor = neighbors.popSmall();
        neighborIndices[i++] = neighbor.second;
    }
    return neighborIndices;
}

/**
 * @brief distanceRayToPoint Calculates the distance of the picking ray to the selected point
 * @param rayStart initial point of the ray
 * @param rayDirection direction of the ray
 * @param point selected point
 * @return the distance between point and ray
 */
float distanceRayToPoint(glm::vec3 rayStart, glm::vec3 rayDirection, glm::vec3 point)
{
    return glm::length(glm::cross(rayDirection, point - rayStart));
}

/**
 * @brief CgKdTree::getClosestPointToRay Returns the index of the point which is the nearest to the picking ray
 * @param rayStart initial point of the ray
 * @param rayDirection direction of the ray
 * @return index of the point which is the nearest to the picking ray
 */
size_t CgKdTree::getClosestPointToRay(glm::vec3 rayStart, glm::vec3 rayDirection)
{
    // index and distance which are closest to the ray
    size_t closestIndex = 0;
    float closestDist = INFINITY;

    // go through all positions and look if there are any which are closer to the ray than the
    // previously found
    for (size_t i = 0; i < vertPositions.size(); i++) {
        glm::vec3 position = vertPositions[i];
        float dist = distanceRayToPoint(rayStart, rayDirection, position);
        if (dist < closestDist) {
            closestIndex = i;
            closestDist = dist;
        }
    }
    return closestIndex;
}

/**
 * @brief CgKdTree::getSplitPlanes Returns the split planes of a KdTree for a given depth
 * @param maxDepth given max depth
 * @return all split planes of the KdTree
 */
std::vector<SplitPlane*> CgKdTree::getSplitPlanes(size_t maxDepth)
{
    return root->getSplitPlanes(0, maxDepth);
}

/* CgKdTreeNode*/

/**
 * @brief CgKdTreeNode::CgKdTreeNode Constructor for creating new KdTree nodes
 * @param kdTree KdTree for which the node will be initialised
 * @param parent Parent of the new node
 * @param vertexIndices
 * @param depth depth of the new node
 * @param maxDepth max depth of the tree
 */
CgKdTreeNode::CgKdTreeNode(CgKdTree* kdTree, CgKdTreeNode* parent, std::vector<size_t>& vertexIndices, int depth, int maxDepth)
    : kdTree(kdTree)
    , parent(parent)
{
    this->vertIndices = vertexIndices;

    // look up which axis to use for further splitting the vertices of this node
    // modulo 3: there are three axes, x, y and z, which we cycle through as we go deeper into the tree
    splitAxis = static_cast<CgKdTree::axis3d>(depth % 3);

    // get split values based on median position along split axis (normal of split plane)
    if (depth < maxDepth) {
        std::vector<size_t> lesserIndices; // left indices
        std::vector<size_t> greaterIndices; // right indices

        std::vector<float> vertAxisPositions;
        for (size_t vertIndex : vertexIndices) {

            // get Vertex Position via the kdTree
            glm::vec3 vertPos = kdTree->getVertPos(vertIndex);

            // position on split axis
            float axisPos = vertPos[splitAxis];
            vertAxisPositions.push_back(axisPos);
        }

        // sort all Positions on axis from smallest to largest
        std::sort(vertAxisPositions.begin(), vertAxisPositions.end(), std::less<float>());
        if(vertAxisPositions.size()>0){
            splitNode = vertAxisPositions[(int)vertAxisPositions.size() / 2]; // median
        }


        // look which positions are smaller and larger (on the axis) than the median (splitNode)
        for (size_t vertIndex : vertexIndices) {
            glm::vec3 vertPos = kdTree->getVertPos(vertIndex);
            if (vertPos[splitAxis] <= splitNode) {
                lesserIndices.push_back(vertIndex);
            } else {
                greaterIndices.push_back(vertIndex);
            }
        }

        // creat new nodes (one for smaller one for larger)
        lesserChild = new CgKdTreeNode(kdTree, this, lesserIndices, depth + 1, maxDepth);
        greaterChild = new CgKdTreeNode(kdTree, this, greaterIndices, depth + 1, maxDepth);

    } else {
        // no children
        lesserChild = NULL;
        greaterChild = NULL;
    }
}

/**
 * @brief CgKdTreeNode::getNearestNeighbors Returns the nearest neighbors which are located in the current node.
 * @param centerIndex index of the point for which the nearest neighbors are sought
 * @param neighborCount number of nearest neighbors searched for
 * @return a priority queue with the nearest neighbors
 */
LimitedSizePriorityQueue<DistAndIndex> CgKdTreeNode::getNearestNeighbors(size_t centerIndex, unsigned int neighborCount)
{
    // main vertex of neighbors
    glm::vec3 centerPos = kdTree->getVertPos(centerIndex);

    // distance between centerPos and SplitPlane
    float distCenterToSplitPlane = distanceToSplitPlane(centerPos);

    CgKdTreeNode* firstCheckedChild = NULL;
    CgKdTreeNode* secondCheckedChild = NULL;

    // check where the neighbors are
    if (centerPos[splitAxis] <= splitNode) {
        if (lesserChild) {
            firstCheckedChild = lesserChild;
        }
        if (greaterChild) {
            secondCheckedChild = greaterChild;
        }
    } else {
        if (greaterChild) {
            firstCheckedChild = greaterChild;
        }
        if (lesserChild) {
            secondCheckedChild = lesserChild;
        }
    }

    // if it exists
    if (firstCheckedChild) {
        // Check, whether the split plane for the current node is closer than then furthest away neighbor found so far
        // If this is the case, we also need to search for neighbors on the other side of the split plane
        LimitedSizePriorityQueue<DistAndIndex> firstCheckedNeighbors = lesserChild->getNearestNeighbors(centerIndex, neighborCount);
        DistAndIndex firstCheckedFurthestNeighbors = firstCheckedNeighbors.peekLarge();

        // if distance between centerPos and Splitplane is smaller than the distance of the first furthest neighbor
        if (distCenterToSplitPlane < firstCheckedFurthestNeighbors.first) {
            if (secondCheckedChild) {
                // Also get the nearest neighbors for the other side of the split plane and combine the results
                LimitedSizePriorityQueue<DistAndIndex> secondCheckedNeighbors = greaterChild->getNearestNeighbors(centerIndex, neighborCount);
                while (!secondCheckedNeighbors.empty()) {
                    firstCheckedNeighbors.push(secondCheckedNeighbors.popLarge());
                }
            }
        }
        return firstCheckedNeighbors;
    } else {
        if (secondCheckedChild) {
            return secondCheckedChild->getNearestNeighbors(centerIndex, neighborCount);
        } else {
            // Search in elements of this node
            LimitedSizePriorityQueue<DistAndIndex> neighbours(neighborCount);

            for (size_t potentialNeighborIndex : vertIndices) {
                if (potentialNeighborIndex == centerIndex) {
                    continue;
                }
                // get new neighbor position and calculate distance to main position
                glm::vec3 neighborPos = kdTree->getVertPos(potentialNeighborIndex);
                float distCenterNeighbor = glm::length(neighborPos - centerPos);
                neighbours.push(std::make_pair(distCenterNeighbor, potentialNeighborIndex));
            }

            return neighbours;
        }
    }
}

/**
 * @brief CgKdTreeNode::getSplitPlanes Return all split planes for this node and all child nodes
 * @param currentDepth The current depth into the kd tree
 * @param maxDepth The maximum depth to go into the kd tree (i.e. the maximum depth of split planes to return)
 * @return All split planes for this node and all child nodes
 */
std::vector<SplitPlane*> CgKdTreeNode::getSplitPlanes(size_t currentDepth, size_t maxDepth)
{
    // Check, whether this node actually further splits data, or is just a leaf holding un-split data
    if (!(lesserChild || greaterChild)) {
        return {};
    }

    // Collect all split planes for this node and all child nodes
    std::vector<SplitPlane*> splitPlanes;

    // Get a bounding box for the vertex positions in this node
    std::pair<glm::vec3, glm::vec3> boundingBox = getBoundingBox();

    // Read the minimum and maximum corner of the bounding box
    glm::vec3 min = boundingBox.first;
    glm::vec3 max = boundingBox.second;

    // Create a new split plane for the current nodes
    // The split plane center is based on the split value
    // The split plane orientation is based on the split axis
    // The split plane size is based on the bounding box
    SplitPlane* splitPlane = new SplitPlane(min, max, splitAxis, splitNode);
    splitPlanes.push_back(splitPlane);

    // Go deeper into the kd tree and add all split planes for all child nodes
    if (currentDepth + 1 < maxDepth) {
        if (lesserChild) {
            std::vector<SplitPlane*> lesserSplitPlanes = lesserChild->getSplitPlanes(currentDepth + 1, maxDepth);
            splitPlanes.insert(splitPlanes.end(), lesserSplitPlanes.begin(), lesserSplitPlanes.end());
        }
        if (greaterChild) {
            std::vector<SplitPlane*> greaterSplitPlanes = greaterChild->getSplitPlanes(currentDepth + 1, maxDepth);
            splitPlanes.insert(splitPlanes.end(), greaterSplitPlanes.begin(), greaterSplitPlanes.end());
        }
    }

    return splitPlanes;
}

/**
 * @brief CgKdTreeNode::distanceToSplitPlane Calculate the distance of the given position to the split plane of the current node
 * @param pos The given position
 * @return The distance of the given position to the split plane of the current node
 */
float CgKdTreeNode::distanceToSplitPlane(glm::vec3 pos)
{
    // The distance is only based on the difference in the split axis
    return abs(splitNode - pos[splitAxis]);
}

/**
 * @brief CgKdTreeNode::getBoundingBox Get a bounding box for the vertex positions in this node
 * @return The bounding box as a pair of minimum and maximum corner for the vertex positions in the current node
 */
std::pair<glm::vec3, glm::vec3> CgKdTreeNode::getBoundingBox()
{
    // Start with infinite values for minimum and maximum
    glm::vec3 min = { INFINITY, INFINITY, INFINITY };
    glm::vec3 max = { -INFINITY, -INFINITY, -INFINITY };

    // Update the minimum and maximum values for all vertices
    for (size_t i : vertIndices) {
        glm::vec3 pos = kdTree->getVertPos(i);
        // Update the minimum and maximum for each dimension individually
        for (int dim = 0; dim < 3; dim++) {
            if (pos[dim] < min[dim]) {
                min[dim] = pos[dim];
            }
            if (pos[dim] > max[dim]) {
                max[dim] = pos[dim];
            }
        }
    }
    return std::make_pair(min, max);
}

template <typename ElementType>
/**
 * @brief LimitedSizePriorityQueue<ElementType>::LimitedSizePriorityQueue Initialize the queue with empty data and the given maximum size
 * @param maxSize The maximum amount of elements to store in this queue before removing old ones upon insertion
 */
LimitedSizePriorityQueue<ElementType>::LimitedSizePriorityQueue(size_t maxSize)
    : data()
    , maxSize(maxSize)
{
}

/* Split Plane */

/**
 * @brief SplitPlane::SplitPlane Create a new split plane based on a bounding box and information about the split
 * @param min The minimum position of the split plane
 * @param max The maximum position of the split plane
 * @param splitAxis The axis on which this split plane splits the values
 * @param splitValue The value on the split axis to split based on
 */
SplitPlane::SplitPlane(glm::vec3 min, glm::vec3 max, CgKdTree::axis3d splitAxis, float splitValue)
    : splitAxis(splitAxis)
    , splitValue(splitValue)
{
    // Ensure, that the min and max are on the same axis-aligned split plane
    min[splitAxis] = splitValue;
    max[splitAxis] = splitValue;

    /* Use the minimum and maximum value as two opposing corners
     *
     *  1 ---- 2
     *  |      |
     *  |      |
     *  4 ---- 3
     */

    pos1 = min;
    pos3 = max;

    // Set the other two corners based on the split axis
    switch (splitAxis) {
    case CgKdTree::axis3d::x:
        pos2 = { splitValue, min[1], max[2] };
        pos4 = { splitValue, max[1], min[2] };
        break;
    case CgKdTree::axis3d::y:
        pos2 = { min[0], splitValue, max[2] };
        pos4 = { max[0], splitValue, min[2] };
        break;
    case CgKdTree::axis3d::z:
        pos2 = { min[0], max[1], splitValue };
        pos4 = { max[0], min[1], splitValue };
        break;
    }
}
