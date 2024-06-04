#ifndef CGMOVINGLEASTSQUARES_H
#define CGMOVINGLEASTSQUARES_H

#include <CgMath/Eigen/Core>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>
#include <vector>

class Plane;

/**
 * @brief The CgMovingLeastSquares class
 */
class CgMovingLeastSquares {
private:
    Eigen::VectorXd coefficients; // The coefficients for the fitted function
    size_t bivariateFunctionDegree; // The bivariate function degree of the  function

public:
    void fitBivariate(
            std::vector<glm::vec2> samplingPoints,
            std::vector<float> samplingValues,
            size_t bivariateFunctionDegree);
    float evaluateBivariate(glm::vec2 samplingPoint);
    glm::vec3 smoothBivariate(glm::vec3 center,
                              Plane tangentPlane,
                              std::vector<glm::vec3> neighbors,
                              size_t bivariateFunctionDegree);
    glm::vec3 smoothBivariate(
            glm::vec3 center,
            std::vector<glm::vec3> neighbors,
            size_t bivariateFunctionDegree);
    static Plane covarianceAnalysis(
            glm::vec3 center,
            std::vector<glm::vec3> neighbors);

private:
    static Eigen::MatrixXd calculateMoorePenrosePseudoInverse(Eigen::MatrixXd A);
    static Eigen::VectorXd evaluateBivariateBaseFunctions(float s, float t, size_t bivariateFunctionDegree);
};

/**
 * @brief Represent a geometric vector-based plane
 */
class Plane {
private:
    glm::vec3 origin; // Origin of the plane in global coordinates, which is at (0,0,0) in plane coordinates
    glm::vec3 normal; // Normal of the plane in global coordinates, which is (0,0,1) in plane coordinates

    glm::vec3 spanningVectorS; // First spanning vector of the plane
    glm::vec3 spanningVectorT; // Second spanning vector of the plane

    glm::mat3x3 base; // The base of the plane, based on the spanning vectors and normal
    glm::mat3x3 baseInverted; // The inverted base of the plane, based on the spanning vectors and normal

public:
    Plane(glm::vec3 origin, glm::vec3 normal);
    Plane(glm::vec3 origin, glm::vec3 normal, glm::vec3 spanningVectorS, glm::vec3 spanningVectorT);
    glm::vec3 projectPoint(glm::vec3 point);
    glm::vec3 unprojectPoint(glm::vec3 point);

    glm::vec3 getNormal();

private:
    glm::vec3 getPerpendicularVector(glm::vec3 vector);
};

#endif // CGMOVINGLEASTSQUARES_H
