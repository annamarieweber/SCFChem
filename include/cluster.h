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

public:
  int K;
  mat atomMatrix;
  int p; // spin up electrons
  int q; // spin down electrons
  int numElectronPairs;
  int numValenceElectrons;
  int numBasisFunctions;
  vec valenceElectronCountsVec;

  Cluster();

  /**
   * @brief Cluster constructor
   * @details Creates a cluster of atoms represented as a matrix with numAtoms rows and 4 columns
   * where the first column is the atomic number and the remaining there columns are the x, y, and z coordinates
   * of the given atom respectively.
   * @param numAtoms int: the number of atoms to be included in the cluster
   **/
  Cluster(int numAtoms);

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

  mat molecularOrbitalCoefficients();

  vec eigenvalues();

  /*;
   * @brief uses eigenvalues to calculate total molecular energy
   * @return double
   */
  double molecularEnergy();

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
