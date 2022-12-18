/**
 * @file molutils.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief Defines utility methods for calculations on molecules
 * @version 0.1
 * @date 2022-12-17
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef UTILS_H
#define UTILS_H
#include <armadillo>
#include "constants.h"
#include "cluster.h"

using namespace arma;
using namespace constants;

/**
 * @brief calculates the distance between two atoms represented by 1d matrices with 3 values for x,y,and z respectively
 * @param a1 (mat): matrix representing the first atom
 * @param a2 (mat): matrix representing the second atom
 * @return double: the distance between the two atoms described as having positions a1 and a2
 **/
double calcDistance(mat a1, mat a2);

/**
 * @brief resizes matrix x to be sized by electrons rather than atoms
 *
 * @param x (mat): an n * n matrix where n = a number of atoms
 * @param atomMatrix (atoms): the matrix of atoms
 * @param electronCount (int): the number of electrons in the cluster
 * @return mat
 */
mat broadcastToOrbitals(mat x, mat atomMatrix, int electronCount);

/**
 * @brief resizes matrix x with dimensional flattened matrices to be sized by electron count rather than atoms
 *
 * @param x (mat): an 3 * (n * n) matrix where n = a number of atoms
 * @param atomMatrix (atoms): the matrix of atoms
 * @param electronCount (int): the number of electrons in the cluster
 * @return mat
 */
mat broadcastToOrbitalsRA(mat x, mat atomMatrix, int electronCount);

/**
 * @brief Calculate the number of electron pairs n for the atomMatrix
 * @param atomMatrix (mat): the matrix of atoms
 * @param expectEvenElectrons (bool): if the molecule should have an even number of electrons defaults to false
 * @return int
 *
 * @throws BaseException invalid electron configuration if expectsEvenElectrons and electron count is odd
 */
int countElectronPairs(mat atomMatrix, bool expectEvenElectrons);

/**
 * @brief Calculate the number of spin up electrons
 * @param atomMatrix (mat): the matrix of atoms
 * @return int
 */
int spinUpElectrons(mat atomMatrix);

/**
 * @brief Calculate the number of spin down electrons
 * @param atomMatrix (mat): the matrix of atoms
 * @return int
 */
int spinDownElectrons(mat atomMatrix);

/**
 * @brief Calculate the number of basis functions N for the described cluster
 * @param atomMatrix (mat): the matrix of atoms
 * @return int
 */
int countBasisFunctions(mat atomMatrix);

/**
 * @brief Check if all atoms in matrix match provided atomic number
 * @param num (int): the atomic number to match
 * @param a mat the matrix to check
 * @return bool
 **/
bool allAtomsMatch(int num, mat a);

/**
 * @brief Get a vector of valence atomic numbers for each atom in atomMatrix
 * @param atomMatrix (mat): the matrix of atoms
 * @return vec
 */
vec valenceElectronCounts(mat atomMatrix);

#endif