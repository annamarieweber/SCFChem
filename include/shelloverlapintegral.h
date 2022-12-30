/**
 * @class ShellOverlapIntegral
 * @brief A class for calculating electron shell overlap integrals
 *
 * The ShellOverlapIntegral class provides a way to calculate electron shell overlap
 * integrals between two shells, represented by the Shell class.
 *
 * @note The ShellOverlapIntegral class depends on the Shell class, which should be
 * included in the "shell.h" header file.
 *
 * @author Anna Weber (anna@scfchem.com)
 * @date 2022-12-17
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef SHELL_OVERLAP_INTEGRAL
#define SHELL_OVERLAP_INTEGRAL
#include "shell.h"

using arma::cube;
using std::vector;

class ShellOverlapIntegral
{
private:
    /**
     * @brief The first Shell object
     *
     * The first Shell object for which the overlap integral will be calculated.
     */
    Shell _s_a;

    /**
     * @brief The second Shell object
     *
     * The second Shell object for which the overlap integral will be calculated.
     */
    Shell _s_b;

    /**
     * @brief The overlap integral matrix
     *
     * The overlap integral matrix between the two Shell objects.
     */
    mat _overlap;

    /**
     * @brief Calculate the product of the alpha values for each exponent in the two Shell objects
     *
     * Calculate the product of the alpha values for each exponent in the two Shell objects.
     *
     * @return A vector of alpha value products
     */
    vec alphas_product();

    /**
     * @brief Calculate the sum of the alpha values for each exponent in the two Shell objects
     *
     * Calculate the sum of the alpha values for each exponent in the two Shell objects.
     *
     * @return A vector of alpha value sums
     */
    vec alphas_sum();

    /**
     * @brief Calculate the squared distance between the centers of the two Shell objects
     *
     * Calculate the squared distance between the centers of the two Shell objects.
     *
     * @return A vector of squared distances
     */
    vec dim_dist_sqr();

    /**
     * @brief Calculate the exponential prefactor for each exponent in the two Shell objects
     *
     * Calculate the exponential prefactor for each exponent in the two Shell objects.
     *
     * @return A vector of exponential prefactors
     */
    vec exponential_prefactor();

    /**
     * @brief Calculate the root term for each exponent in the two Shell objects
     *
     * Calculate the root term for each exponent in the two Shell objects.
     *
     * @return A vector of root terms
     */
    vec root_term();

    /**
     * @brief Calculate the overlap summation for each exponent in the two Shell objects
     *
     * Calculate the overlap summation for each exponent in the two Shell objects.
     *
     * @param[in] x_p A vector of alpha value products
     * @param[in] l_pair_a The angular momentum pair for the first Shell object
     * @param[in] l_pair_b The angular momentum pair for the second Shell object
     *
     * @return A vector of overlap summations
     */
    vec overlap_summation(vec x_p, int lpair_a, int l_pair_b);

    /**
     * @brief Calculate the overlap integral matrix between the two Shell objects
     *
     * Calculate the overlap integral matrix between the two Shell objects.
     *
     * @return The overlap integral matrix
     */
    mat overlap();

    /**
     * @brief Calculate the product center for the two Shell objects
     *
     * Calculate the product center for the two Shell objects.
     *
     * @return A vector representing the product center
     */
    vec product_center();

public:
    /**
     * @brief Default constructor for the ShellOverlapIntegral class
     *
     * Default constructor for the ShellOverlapIntegral class. This constructs an empty
     * ShellOverlapIntegral object with default-initialized Shell objects.
     */
    ShellOverlapIntegral();

    /**
     * @brief Constructor for the ShellOverlapIntegral class
     *
     * Constructor for the ShellOverlapIntegral class. This constructs a ShellOverlapIntegral
     * object with the specified Shell objects.
     *
     * @param[in] s_a The first Shell object
     * @param[in] s_b The second Shell object
     */
    ShellOverlapIntegral(Shell s_a, Shell s_b);

    /**
     * @brief Overloaded function call operator to retrieve the overlap integral matrix
     *
     * Overloaded function call operator to calculate the overlap integral matrix between the two
     * Shell objects.
     *
     * @return The overlap integral matrix
     */
    mat operator()();

    /**
     * @brief Calculate the derivative integral matrix with respect to the first Shell object
     *
     * Calculate the derivative integral matrix with respect to the first Shell object.
     *
     * @return The derivative integral matrix
     */
    mat saDerivativeIntegral();

    /**
     * @brief Calculate the derivative integral matrix with respect to the second Shell object
     *
     * Calculate the derivative integral matrix with respect to the second Shell object.
     *
     * @return The derivative integral matrix
     */
    mat sbDerivativeIntegral();
};
#endif
