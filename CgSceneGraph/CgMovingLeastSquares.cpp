#include "CgMovingLeastSquares.h"

#include "CgMath/Eigen/Eigenvalues"
#include <CgMath/Eigen/Core>
#include <CgMath/Eigen/SVD>
#include <glm/ext.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/io.hpp>
#include <iostream>
#include <map>
#include <string>

/**
 * @brief CgMovingLeastSquares::fitBivariate Fit a bivariate function to the given sampling points and values
 * @param samplingPoints The sampling point to fit based on
 * @param samplingValues The sampling values to fit based on
 * @param bivariateFunctionDegree The degree of the bivariate function to fit
 */
void CgMovingLeastSquares::fitBivariate(
        std::vector<glm::vec2> samplingPoints,
        std::vector<float> samplingValues,
        size_t bivariateFunctionDegree)
{
    // Create the fitting matrix for all function values of the bivariate base functions evaluated at the sampling points
    // Rows:    Sampling points
    // Columns: Bivariate base function component evaluated for sampling point
    Eigen::MatrixXd A(samplingPoints.size(), (int)pow(bivariateFunctionDegree + 1, 2));
    for (size_t j = 0; j < samplingPoints.size(); j++) {
        glm::vec2 samplingPoint = samplingPoints[j];
        Eigen::VectorXd functionValues = evaluateBivariateBaseFunctions(samplingPoint[0], samplingPoint[1], bivariateFunctionDegree);
        for (long i = 0; i < functionValues.rows(); i++) {
            A(j, i) = functionValues[i];
        }
    }

    // Calculate the (pseudo) inverse for the fitting matrix
    Eigen::MatrixXd A_plus = calculateMoorePenrosePseudoInverse(A);

    // Transform the sampling values into the format used by the Eigen library
    Eigen::VectorXd samplingValuesEigen(samplingValues.size());
    for (size_t i = 0; i < samplingValues.size(); i++) {
        samplingValuesEigen[i] = samplingValues[i];
    }

    // Store the bivariate function degree used
    this->bivariateFunctionDegree = bivariateFunctionDegree;

    // Store the fitted coefficients
    coefficients = A_plus * samplingValuesEigen;
}

/**
 * @brief CgMovingLeastSquares::evaluateBivariate Evaluate the fitted bivariate function and coefficients for the given sampling point
 * Fitting has to be done before calling this function
 * @param samplingPoint The sampling point to evaluate
 * @return The function value of the sampling point
 */
float CgMovingLeastSquares::evaluateBivariate(glm::vec2 samplingPoint)
{
    Eigen::VectorXd baseFunctionValues = evaluateBivariateBaseFunctions(samplingPoint[0], samplingPoint[1], bivariateFunctionDegree);
    return coefficients.dot(baseFunctionValues);
}

/**
 * @brief CgMovingLeastSquares::smoothBivariate Smooth the given point based on the neighbor positions projected onto the normal plane
 * @param center The given center point to smooth
 * @param tangentPlane The tangent plane of the given center point, which is used for the tangent plane projection
 * @param neighbors The neighbors for the given center point to use for smoothing
 * @param bivariateFunctionDegree The degree of the bivariate function to use for smoothing
 * @return The smoothed position of the given point
 */
glm::vec3 CgMovingLeastSquares::smoothBivariate(glm::vec3 center,
                                                Plane tangentPlane, std::vector<glm::vec3> neighbors, size_t bivariateFunctionDegree)
{
    // Collect the sampling points and values by projecting the neighbor points onto the tangent plane
    // Sampling points: Position on the tangent plane
    // Sampling values: Height above the tangent plane
    std::vector<glm::vec2> samplingPoints(neighbors.size());
    std::vector<float> samplingValues(neighbors.size());
    for (size_t i = 0; i < neighbors.size(); i++) {
        glm::vec3 neighbor = neighbors[i];
        glm::vec3 neighborProjected = tangentPlane.projectPoint(neighbor);

        glm::vec2 samplingPoint {
                neighborProjected[0],
                neighborProjected[1]
        };
        float samplingValue = neighborProjected[2];

        samplingPoints[i] = samplingPoint;
        samplingValues[i] = samplingValue;
    }

    // Fit a bivariate function to the sampling points and values
    fitBivariate(samplingPoints, samplingValues, bivariateFunctionDegree);

    glm::vec3 pointProjected = tangentPlane.projectPoint(center);
    float fittedValue = evaluateBivariate(pointProjected);

    glm::vec3 updatedPointProjected = {
            pointProjected[0],
            pointProjected[1],
            fittedValue
    };
    glm::vec3 updatedPoint = tangentPlane.unprojectPoint(updatedPointProjected);

    return updatedPoint;
}

/**
 * @brief CgMovingLeastSquares::Bivariate Smooth the given point based on the neighbor positions projected onto the normal plane
 * @param center The given center point to smooth
 * @param neighbors The neighbors for the given center point to use for smoothing
 * @param bivariateFunctionDegree The degree of the bivariate function to use for smoothing
 * @return The smoothed position of the given center point
 */
glm::vec3 CgMovingLeastSquares::smoothBivariate(
        glm::vec3 center,
        std::vector<glm::vec3> neighbors,
        size_t bivariateFunctionDegree)
{
    Plane plane = covarianceAnalysis(center, neighbors);

    return smoothBivariate(center, plane, neighbors, bivariateFunctionDegree);
}

/**
 * @brief CgMovingLeastSquares::covarianceAnalysis Estimate the tangent plane for a given center point based on its neighbors
 * @param center The given center point to estimate a tangent plane for
 * @param neighbors The neighbors of the given center point to base the tangent plane estimation on
 * @return An estimation of the tangent plane for a given point
 */
Plane CgMovingLeastSquares::covarianceAnalysis(glm::vec3 center, std::vector<glm::vec3> neighbors)
{
    // Calculate the centroid of the neighborhood (P Strich (Schwerpunkt Nachbarschaft))
    glm::vec3 neighborhoodCentroid = { 0, 0, 0 };
    for (glm::vec3 neighbor : neighbors) {
        neighborhoodCentroid += neighbor;
    }
    neighborhoodCentroid /= neighbors.size();

    // Set up the covariance matrix
    Eigen::MatrixXd covarianceMatrixPartial(neighbors.size(), 3);
    for (size_t row = 0; row < neighbors.size(); row++) {
        glm::vec3 neighbor = neighbors[row];
        glm::vec3 rowValues = neighbor - neighborhoodCentroid;
        // Convert the rows from GLM to Eigen format
        for (size_t col = 0; col < 3; col++) {
            covarianceMatrixPartial(row, col) = rowValues[col];
        }
    }
    Eigen::Matrix3d covarianceMatrix = covarianceMatrixPartial.transpose() * covarianceMatrixPartial;

    // Calculate the eigenvectors and eigenvalues for the covariance matrix
    Eigen::EigenSolver<Eigen::Matrix3d> solver;
    solver.compute(covarianceMatrix, true);
    Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXcd eigenvectors = solver.eigenvectors();

    // Sort the eigenvectors and eigenvectors by eigenvalue from smallest to largest
    std::vector<std::pair<double, glm::vec3>> eigenvaluesSorted(3);
    for (size_t i = 0; i < 3; i++) {
        double eigenvalue = eigenvalues(i).real();
        Eigen::Vector3cd eigenvector = eigenvectors.col(i);
        glm::vec3 eigenvectorGlm;
        // Convert the eigenvectors from Eigen to GLM format
        for (size_t j = 0; j < 3; j++) {
            eigenvectorGlm[j] = eigenvector[j].real();
        }
        eigenvaluesSorted[i] = std::make_pair(eigenvalue, eigenvectorGlm);
    }
    std::sort(eigenvaluesSorted.begin(), eigenvaluesSorted.end(), [](const std::pair<double, glm::vec3> a, const std::pair<double, glm::vec3> b) {
        return a.first < b.first;
    });

    // Create the estimated tangent plane based on the eigenvectors
    // The eigenvector for the smallest eigenvalue is assumed to be the normal
    glm::vec3 normal = eigenvaluesSorted[0].second;
    glm::vec3 spanningVectorS = eigenvaluesSorted[1].second;
    glm::vec3 spanningVectorT = eigenvaluesSorted[2].second;
    Plane plane(center, normal, spanningVectorS, spanningVectorT);

    return plane;
}

/**
 * @brief CgMovingLeastSquares::calculateMoorePenrosePseudoInverse Calculate the Moore-Penrose-Pseudo-Inverse for the given matrix
 * @param A The given matrix
 * @return The Moore-Penrose-Pseudo-Inverse for the given matrix
 */
Eigen::MatrixXd CgMovingLeastSquares::calculateMoorePenrosePseudoInverse(Eigen::MatrixXd A)
{
    // Compute SVD decomposition
    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::VectorXd SV = svd.singularValues();

    // Build a diagonal matrix out of the singular values
    Eigen::MatrixXd S(A.rows(), A.cols());
    S.setZero(A.rows(), A.cols());
    for (int i = 0; i < SV.rows(); i++) {
        S(i, i) = SV(i);
    }

    // Compute Moore-Penrose inverse
    Eigen::MatrixXd S_inv = S.transpose();
    for (int i = 0; i < SV.rows(); i++) {
        S_inv(i, i) = 1.0 / S_inv(i, i);
    }

    Eigen::MatrixXd A_plus = V * S_inv * U.transpose();

    /* DEBUG OUTPUT */
    /*
    {
        std::cout << "A:\n"
                  << A << std::endl;
        std::cout << "U:\n"
                  << U << std::endl;
        std::cout << "S:\n"
                  << S << std::endl;
        std::cout << "V^T:\n"
                  << V.transpose() << std::endl;
        std::cout << "A=USV^T:\n"
                  << U * S * V.transpose() << std::endl;
        std::cout << "S^-1:\n"
                  << S_inv << std::endl;
        std::cout << "A+:\n"
                  << A_plus << std::endl;
    }
    */

    return A_plus;
}

/**
 * @brief CgMovingLeastSquares::evaluateBivariateBaseFunctions Evaluate a bivariate base function of the given degree for the given parameters
 * @param s First parameter for the bivariate function
 * @param t Second parameter for the bivariate function
 * @param bivariateFunctionDegree Degree of the bivariate function
 * @return The bivariate function value for the given degree and parameters
 */
Eigen::VectorXd CgMovingLeastSquares::evaluateBivariateBaseFunctions(float s, float t, size_t bivariateFunctionDegree)
{
    Eigen::VectorXd functionValues((int)pow(bivariateFunctionDegree + 1, 2));

    // Iterate through all combinations of potencies for s and t smaller than or equal to the degree of the bivariate function
    for (size_t power_t = 0; power_t <= bivariateFunctionDegree; power_t++) {
        for (size_t power_s = 0; power_s <= bivariateFunctionDegree; power_s++) {
            // The unique position in the function values array is calculated in the same way
            // as you would access a matrix of size (bivariateFunctionDegree + 1)^2.
            // There are bivariateFunctionDegree values of t for each value of s,
            // resulting in a total of (bivariateFunctionDegree + 1)^2 values.
            // The function value is calculated as: s^power_s * t^power_t
            functionValues[power_t * (bivariateFunctionDegree + 1) + power_s] = pow(s, power_s) * pow(t, power_t);
        }
    }

    return functionValues;
}

/**
 * @brief Plane::Plane Create a plane based on the origin and normal
 * @param origin Origin of the plane
 * @param normal Normal of the plane
 */
Plane::Plane(glm::vec3 origin, glm::vec3 normal)
        : origin(origin)
        , normal(normal)
{
    // Create two arbitrary vectors spanning the plane
    // They need to be perpendicular to each other as well as to the normal of the plane
    spanningVectorS = getPerpendicularVector(normal);
    spanningVectorT = glm::cross(normal, spanningVectorS);

    // Create the base and inverted base for describing points on the plane
    base[0] = spanningVectorS;
    base[1] = spanningVectorT;
    base[2] = normal;
    baseInverted = glm::inverse(base);
}

/**
 * @brief Plane::Plane Create a plane based on the origin, normal and spanning vectors
 * @param origin Origin of the plane
 * @param normal Normal of the plane
 * @param spanningVectorS Spanning vector of the plane
 * @param spanningVectorS Spanning vector of the plane
 */
Plane::Plane(glm::vec3 origin, glm::vec3 normal, glm::vec3 spanningVectorS, glm::vec3 spanningVectorT)
        : origin(origin)
        , normal(normal)
        , spanningVectorS(spanningVectorS)
        , spanningVectorT(spanningVectorT)
{
    // Create the base and inverted base for describing points on the plane
    base[0] = spanningVectorS;
    base[1] = spanningVectorT;
    base[2] = normal;
    baseInverted = glm::inverse(base);
}

/**
 * @brief Plane::projectPoint Project a point onto the plane
 * @param point The point in cartesian coordinates
 * @return The point in 'plane coordinates' relative to the plane origin
 */
glm::vec3 Plane::projectPoint(glm::vec3 point)
{
    // Calculate the position of the point relative to the plane origin and change into the plane base
    // Given E: Cartesian base; B: Plane base
    // A_E->B = (A_B->E)^-1 = B^-1
    return baseInverted * (point - origin);
}

/**
 * @brief Plane::unprojectPoint Reverse-project a point from the plane
 * @param point The point in 'plane coordinates' relative to the plane origin
 * @return The point in cartesian coordinates
 */
glm::vec3 Plane::unprojectPoint(glm::vec3 point)
{
    // Change base into cartesian coordinates and adjust point for plane origin
    // Given E: Cartesian base; B: Plane base
    // A_B->E = B
    return (base * point) + origin;
}

glm::vec3 Plane::getNormal()
{
    return normal;
}

/**
 * @brief Plane::getPerpendicularVector Calculates an arbitrary verctor perpendicular to the given one
 * @param vector The given vector
 * @return An abritrary vector perpendicular to the given vector
 */
glm::vec3 Plane::getPerpendicularVector(glm::vec3 vector)
{
    if ((vector[0] == 0.0) && (vector[1] == 0.0)) {
        if (vector[2] == 0.0)
            return glm::vec3(0.);

        return glm::vec3(0.0, 1.0, 0.0);
    }
    return glm::normalize(glm::vec3(-vector[1], vector[0], 0.0));
}
