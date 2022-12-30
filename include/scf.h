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
#include "timedfunctional.h"
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
         * @brief calls eigsym and updates values for alpha and beta
         *
         * @param e_a
         * @param c_a
         * @param f_a
         * @param e_b
         * @param c_b
         * @param f_b
         */
        int run_eig_sym(vec &e_a, mat &c_a, mat &f_a, vec &e_b, mat &c_b, mat &f_b);

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
         * @param cluster (Cluster): the cluster to run the scf algorithm on
         */
        void calcSCFEnergy(float threshold, Cluster cluster);
    }
}
#endif