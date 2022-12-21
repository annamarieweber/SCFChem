/**
 * @file algoutils.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief Defines utility functions used in algorithm calculations
 * @version 0.1
 * @date 2022-12-21
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef ALGO_UTILS
#define ALGO_UTILS
#include <armadillo>
#include "cluster.h"

using namespace arma;
/**
 * @brief namespace for all algorithms used in SCF chem calculations
 *
 */
namespace algo
{

    /**
     * @brief namespace for utility functions used in algorithms
     *
     */
    namespace utils
    {

        /**
         * @brief calculates the gamma matrix used to calculate the fock matrix
         *
         * @return mat
         */
        mat gammaMatrix(Cluster cluster);

        /**
         * @brief gradient of gamma matrix, rows are x y z, columns are atom pairs (A, B)
         *
         * @return mat
         */
        mat gammaMatrixRA(Cluster cluster);

        /**
         * @brief overlap matrix
         *
         * @return mat
         */
        mat overlapMatrix(Cluster cluster);

        /**
         * @brief gradient of overlap matrix, rows are x y z, columns are overlap matrix elements
         *
         * @return mat
         */
        mat overlapMatrixRA(Cluster cluster);

    }
}
#endif