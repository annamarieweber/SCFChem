/**
 * @file cndo2.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief Definition of the CNDO2 algorithm methods
 * @version 0.1
 * @date 2022-12-21
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef CNDO2
#define CNDO2

#include <armadillo>
#include "cluster.h"
namespace algo
{

    /**
     * @brief namespace for calculating fockmatrix using cndo2 algorithm
     *
     */
    namespace cndo2
    {
        using namespace arma;
        /**
         * @brief an object for storing output of the cndo2 calculation
         */
        struct CNDO2Result
        {
            mat fock;              /**< an n*n matrix representing the fock matrix calculated by the cndo algorithm */
            mat bonding_param_tot; /**< an n*n matrix with the bonding parameter totals calculated while running cndo */
        };

        /**
         * @brief returns the CNDO2 Fock Matrix for the molecule
         *
         * @param p_a (mat): the alpha density matrix
         * @param p_b (mat): the beta density matrix
         *
         * @return CNDO2Result
         */
        CNDO2Result fockMatrix(mat p_a, mat p_b, Cluster cluster);
    }
}
#endif