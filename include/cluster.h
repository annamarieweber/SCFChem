/**
 * @file cluster.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief Definition of the cluster class and its methods
 * @version 0.1
 * @date 2022-12-17
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef CLUSTER
#define CLUSTER
#include <armadillo>
#include <iostream>
#include <fstream>
#include <string>

using arma::mat;
using arma::vec;
using std::string;

class Cluster
{
private:
  vec epsilons;
  vec sigmas;
  int K;
  int _p; // spin up electrons
  int _q; // spin down electrons

public:
  mat atomMatrix;

  Cluster();

  /**
   * @brief Cluster constructor
   * @detail Creates a cluster of atoms represented as a matrix with numAtoms rows and 4 columns
   * where the first column is the atomic number and the remaining there columns are the x, y, and z coordinates
   * of the given atom respectively.
   * @param numAtoms int: the number of atoms to be included in the cluster
   **/
  Cluster(int numAtoms);

  /**
   * @brief Check if all atoms in cluster match provided atomic number
   * @param num int the atomic number to match
   * @return bool
   **/
  bool allAtomsMatch(int num);

  /**
   * @brief Calculate the number of basis functions N for the described cluster
   * @return int
   */
  int countBasisFunctions();

  /**
   * @brief Calculate the number of electron pairs n for the described cluster
   * @return int
   */
  int countElectronPairs();

  /**
   * @brief Returns a matrix representing the basis functions for all of the overlaps
   * @return mat
   */
  mat basisFunctions();

  /**
   * @brief Returns a matrix representing the s basis functions for all of the overlaps
   * @return mat
   */
  mat sBasisFunctions();

  mat broadcastToOrbitals(mat x);

  mat broadcastToOrbitalsRA(mat x);

  /**
   * @brief calculate the x component of the gradient calculation
   *
   * @param p_a
   * @param p_b
   * @param bonding_parameters
   * @return mat
   */
  mat x(mat p_a, mat p_b, mat bonding_parameters);

  /**
   * @brief returns the CNDO2 Fock Matrix for the molecule
   *
   * @return mat
   */
  mat cndo2FockMatrix(mat p_a, mat p_b, mat &bonding_params);

  /**
   * @brief returns the derivative of CNDO2 Fock Matrix for the molecule with respect to Ra
   *
   * @return mat
   */
  mat cndo2FockMatrixRA(mat p_a, mat p_b);

  mat x(mat p_a, mat p_b, vec b_a, vec b_b);

  void calcSCFEnergy(float threshold);

  /**
   * @brief overlap matrix
   *
   * @return mat
   */
  mat overlapMatrix();

  /**
   * @brief gradient of overlap matrix, rows are x y z, columns are overlap matrix elements
   *
   * @return mat
   */
  mat overlapMatrixRA();

  /**
   * @brief  gamma matrix,
   *
   * @return vec
   */
  mat gammaMatrix();

  /**
   * @brief  gradient of gamma matrix, rows are x y z, columns are atom pairs (A, B)
   *
   * @return vec
   */
  mat gammaMatrixRA();

  /**
   * @brief gradient of energy, rows are x y z, columns are atoms
   *
   * @return vec
   */
  void calcSCFEnergyRA(float threshold);

  vec z_vals();

  mat molecularOrbitalCoefficients();

  vec eigenvalues();

  /*;
   * @brief used eigenvalues to calculate total molecular energy
   * @return double
   */
  double molecularEnergy();

  /**
   * @brief calculates the distance between two atoms represented by 1d matrices with 3 values for x,y,and z respectively
   * @param a1 (mat): matrix representing the first atoms
   * @param a2 (mat): matrix representing the second value
   * @return double: the distance between the two atoms described as having positions a1 and a2
   **/
  double calcDistance(mat a1, mat a2);

  /**
   * @brief inserts an atom into the cluster
   * @param index (int): the index where the atom should be inserted
   * @param atomNum (int): the atomic number representing the atom type
   * @param x (double): the x coordinate of the atom
   * @param y (double): the y coordinate of the atom
   * @param z (double): the z coordinate of the atom
   * @param e (double): the epsilon corresponding to the atom at the specified index
   * @param s (double): the sigma corresponding to the atom at the specified index
   **/
  void addAtom(int index, int atomNum, double x, double y, double z, double e, double s);

  /**
   * @brief overloading the << operator
   **/
  friend std::ostream &operator<<(std::ostream &os, const Cluster &c);
};
#endif
