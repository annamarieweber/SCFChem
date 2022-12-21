#include <armadillo>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include "cluster.h"
#include "baseexception.h"
#include "constants.h"
#include "combinatorics.h"
#include "shell.h"
#include "shelloverlapintegral.h"
#include "molutils.h"
#include <cmath>

using std::string;
using namespace constants;
using namespace arma;

Cluster::Cluster(){};

Cluster::Cluster(int numAtoms)
{
  atomMatrix = mat(numAtoms, 4);
  K = 3;
  epsilons = vec(numAtoms);
  sigmas = vec(numAtoms);
  numBasisFunctions = countBasisFunctions(atomMatrix);
  numElectronPairs = countElectronPairs(atomMatrix, false);
  valenceElectronCountsVec = valenceElectronCounts(atomMatrix.col(0));
  numValenceElectrons = sum(valenceElectronCountsVec);
  p = numValenceElectrons / 2 + numValenceElectrons % 2;
  q = numValenceElectrons / 2;
}

mat Cluster::basisFunctions()
{
  vec l(numBasisFunctions);

  mat basisfns(numBasisFunctions, 5 * K + 2);

  int fn = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    vec l_vals(ATOM_ANGULAR_MOMENTUM_MAP[atomMatrix(i, 0) - 1]);
    for (int j = 0; j < l_vals.n_elem; j++)
    {
      Shell s(atomMatrix(i, 0), atomMatrix.row(i).cols(1, atomMatrix.n_cols - 1).t(), l_vals(j));
      ShellOverlapIntegral s_aa(s, s);
      mat norm_consts = 1.0 / sqrt(s_aa());
      for (int p = 0; p < s.num_quantum_arrangements(); p++)
      {
        basisfns.row(fn).cols(0, K) = atomMatrix.row(i).cols(0, K);
        basisfns.row(fn).cols(K + 1, 2 * K) = s.l_a(p).t();
        basisfns.row(fn).cols(2 * K + 1, 3 * K) = s.alpha().t();
        basisfns.row(fn).cols(3 * K + 1, 4 * K) = s.d_k().t();
        basisfns.row(fn).cols(4 * K + 1, 5 * K) = (l_vals(j) + 1) * norm_consts.col(p).t();
        basisfns.row(fn).col(5 * K + 1) = j;
        fn++;
      }
    }
  }
  return basisfns;
}

mat Cluster::sBasisFunctions()
{

  mat basisfns(atomMatrix.n_rows, 5 * K + 2);
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    rowvec l_vals("0 0 0");
    int atomNum = atomMatrix(i, 0);
    mat coeffs(ATOM_COEFFS_MAP[atomNum - 1]);
    Shell s(atomMatrix(i, 0), atomMatrix.row(i).cols(1, atomMatrix.n_cols - 1).t(), coeffs.col(1), coeffs.col(0), l_vals);
    ShellOverlapIntegral s_aa(s, s);
    mat norm_consts = 1.0 / sqrt(s_aa());
    basisfns.row(i).cols(0, K) = atomMatrix.row(i).cols(0, K);
    basisfns.row(i).cols(K + 1, 2 * K) = s.l_a(0).t();
    basisfns.row(i).cols(2 * K + 1, 3 * K) = s.alpha().t();
    basisfns.row(i).cols(3 * K + 1, 4 * K) = s.d_k().t();
    basisfns.row(i).cols(4 * K + 1, 5 * K) = (0 + 1) * norm_consts.col(0).t();
    basisfns.row(i).col(5 * K + 1) = 0;
  }
  return basisfns;
}


mat Cluster::molecularOrbitalCoefficients()
{
  // vec eigval;
  // mat eigvec;
  // mat overlap_matrix = overlapMatrix();
  // mat hamiltonian = extendedHuckelHamiltonian();

  // eig_sym(eigval, eigvec, overlap_matrix);

  // // make the orthogonalization transformation
  // mat x_diag = arma::diagmat(arma::pow(eigval, -0.5));
  // mat x = eigvec * x_diag * arma::trans(eigvec);
  // // form the hamiltonian in the orthogonalized basis
  // mat h_p = arma::trans(x) * hamiltonian * x;
  // // diagonalize
  // vec e;
  // mat C_p;
  // eig_sym(e, C_p, h_p);
  // // Form the MO coefficients:
  // mat C = x * C_p;
  // mat mo_overlap = arma::trans(C) * overlap_matrix * C;

  // return mo_overlap;
}

vec Cluster::eigenvalues()
{
  // vec eigval;
  // mat eigvec;
  // mat overlap_matrix = overlapMatrix();
  // mat hamiltonian = extendedHuckelHamiltonian();

  // eig_sym(eigval, eigvec, overlap_matrix);

  // // make the orthogonalization transformation
  // mat x_diag = arma::diagmat(arma::pow(eigval, -0.5));
  // mat x = eigvec * x_diag * arma::trans(eigvec);
  // // form the hamiltonian in the orthogonalized basis
  // mat h_p = arma::trans(x) * hamiltonian * x;
  // // diagonalize
  // vec e;
  // mat C_p;
  // eig_sym(e, C_p, h_p);

  // return e;
}

void Cluster::addAtom(int index, int atomNum, double x, double y, double z, double e, double s)
{
  atomMatrix(index, 0) = atomNum;
  atomMatrix(index, 1) = x;
  atomMatrix(index, 2) = y;
  atomMatrix(index, 3) = z;
  epsilons(index) = e;
  sigmas(index) = s;
  numBasisFunctions = countBasisFunctions(atomMatrix);
  numElectronPairs = countElectronPairs(atomMatrix, false);
  valenceElectronCountsVec = valenceElectronCounts(atomMatrix.col(0));
  numValenceElectrons = sum(valenceElectronCountsVec);
  p = numValenceElectrons / 2 + numValenceElectrons % 2;
  q = numValenceElectrons / 2;
}

std::ostream &operator<<(std::ostream &os, const Cluster &c)
{
  os << "Atom Matrix: \n"
     << c.atomMatrix << "\n";

  return os;
}
