/**
 * @file scf.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief definition of SCF algorithm methods
 * @version 0.1
 * @date 2022-12-21
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SCF
#define SCF
#include <armadillo>
#include "cluster.h"
namespace algo
{

    using namespace arma;
    /**
     * @brief namespace contains functions for running SCF algorithms
     *
     */
    namespace scf
    {
        /**
         * @brief calculate the x component of the gradient calculation
         *
         * @param p_a
         * @param p_b
         * @param bonding_parameters
         * @return mat
         */
        mat x_component(mat p_a, mat p_b, mat bonding_parameters);

        /**
         * @brief calculates and prints the nuclear, electron, and total energy for the molecule
         * @details calculates the energy using the SCF algorithm running until the input threshold is reached
         *
         * @param threshold (float): the threshold to run the algorithm to
         */
        void calcSCFEnergy(float threshold, Cluster cluster);
    }
}
#endif