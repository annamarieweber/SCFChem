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
  int numElectrons = sum(z_vals());
  _p = numElectrons / 2 + numElectrons % 2;
  _q = numElectrons / 2;
}

int Cluster::countBasisFunctions()
{
  int sum = 0;
  vec atoms = atomMatrix.col(0);
  for (int i = 0; i < atoms.size(); i++)
  {
    sum += ATOM_BASIS_FN_MAP[atoms[i] - 1];
  }
  return sum;
}

int Cluster::countElectronPairs()
{
  int n = countBasisFunctions();
  if (n / 2 * 2 != n)
  {
    return n / 2;
    // throw BaseException("InvalidConfig: number of electron pairs must be an integer");
  }
  else
  {
    return n / 2;
  }
}

mat Cluster::basisFunctions()
{
  int b = countBasisFunctions();
  vec l(b);

  mat basisfns(b, 5 * K + 2);

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
    mat coeffs(ATOM_COEFS_MAP[atomNum - 1]);
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

vec Cluster::z_vals()
{
  mat valenceElectrons = atomMatrix.col(0);
  valenceElectrons.for_each([](vec::elem_type &val)
                            { val = VALENCE_ATOMIC_NUM[val - 1]; });
  return valenceElectrons;
}

mat Cluster::cndo2FockMatrix(mat p_a, mat p_b, mat &bonding_params)
{
  mat gamma = gammaMatrix();

  mat p_tot = p_a + p_b;

  mat gamma_expanded = broadcastToOrbitals(gamma);
  mat z = z_vals();
  int electron_count = sum(sum(z));

  vec gamma_z_tot = (gamma - diagmat(diagvec(gamma))) * z;

  vec diag_gamma_z = (z.col(0) - (vec(z.n_rows, fill::ones) / 2.0)) % diagvec(gamma);

  vec expanded_diag_gamma_z(electron_count);
  vec expanded_gamma_z_tot(electron_count);
  vec ion_energy_electron_affinity(electron_count);
  vec atoms = atomMatrix.col(0);
  mat h_off_diag(electron_count, electron_count, fill::ones);
  h_off_diag -= diagmat(vec(electron_count, fill::ones));
  vec density_tot(electron_count);
  mat g_off_diag(electron_count, electron_count, fill::ones);
  g_off_diag -= diagmat(vec(electron_count, fill::ones));
  vec density_gamma_off_diag(electron_count);
  vec density_tot_by_atom(atoms.n_elem);

  int k = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    ion_energy_electron_affinity.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) = vec(ATOM_TO_IONIZATION_ENERGY_ELECTRON_AFFINITY_PARAMS_MAPPING[atoms(i) - 1]);
    expanded_diag_gamma_z.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(diag_gamma_z(i));
    expanded_gamma_z_tot.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(gamma_z_tot(i));
    bonding_params.rows(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) += mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], electron_count).fill(ATOMIC_BONDING_PARAMETERS[atoms(i) - 1]);
    bonding_params.cols(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) += mat(electron_count, VALENCE_ATOMIC_NUM[atoms(i) - 1]).fill(ATOMIC_BONDING_PARAMETERS[atoms(i) - 1]);
    float density_a = sum(diagvec(p_tot.submat(k, k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1)));
    density_tot.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(density_a);
    density_tot_by_atom(i) = density_a;

    k += VALENCE_ATOMIC_NUM[atoms(i) - 1];
  }

  vec gamma_density_sum = (gamma - diagmat(diagvec(gamma))) * density_tot_by_atom;
  vec g_diag = (density_tot - diagvec(p_a)) % diagvec(gamma_expanded) + diagvec(broadcastToOrbitals(diagmat(gamma_density_sum)));
  g_off_diag = g_off_diag % (-1.0 * p_a % gamma_expanded);

  h_off_diag %= (1.0 / 2.0 * bonding_params % overlapMatrix());

  vec h_diag = -1 * ion_energy_electron_affinity - expanded_diag_gamma_z - expanded_gamma_z_tot;
  return diagmat(h_diag) + diagmat(g_diag) + h_off_diag + g_off_diag;
}

mat Cluster::cndo2FockMatrixRA(mat p_a, mat p_b)
{
  mat gamma = gammaMatrixRA();

  mat p_tot = p_a + p_b;

  mat gamma_expanded = broadcastToOrbitals(gamma);
  mat z = z_vals();
  int electron_count = sum(sum(z));

  std::cout << gamma_expanded << std::endl;
  std::cout << diagmat(diagvec(gamma)) << std::endl;

  mat x_gamma = resize(mat(gamma.row(0)), atomMatrix.n_rows, atomMatrix.n_rows);
  mat y_gamma = resize(mat(gamma.row(1)), atomMatrix.n_rows, atomMatrix.n_rows);
  mat z_gamma = resize(mat(gamma.row(2)), atomMatrix.n_rows, atomMatrix.n_rows);

  mat x_gamma_exp = resize(mat(gamma_expanded.row(0)), electron_count, electron_count);
  mat y_gamma_exp = resize(mat(gamma_expanded.row(1)), electron_count, electron_count);
  mat z_gamma_exp = resize(mat(gamma_expanded.row(2)), electron_count, electron_count);

  vec gamma_z_tot_x = (x_gamma - diagmat(diagvec(x_gamma))) * z;
  vec gamma_z_tot_y = (y_gamma - diagmat(diagvec(y_gamma))) * z;
  vec gamma_z_tot_z = (z_gamma - diagmat(diagvec(z_gamma))) * z;

  vec diag_gamma_z_x = (z.col(0) - (vec(z.n_rows, fill::ones) / 2.0)) % diagvec(x_gamma);
  vec diag_gamma_z_y = (z.col(0) - (vec(z.n_rows, fill::ones) / 2.0)) % diagvec(y_gamma);
  vec diag_gamma_z_z = (z.col(0) - (vec(z.n_rows, fill::ones) / 2.0)) % diagvec(z_gamma);

  vec expanded_diag_gamma_z_x(electron_count);
  vec expanded_diag_gamma_z_y(electron_count);
  vec expanded_diag_gamma_z_z(electron_count);

  vec expanded_gamma_z_tot_x(electron_count);
  vec expanded_gamma_z_tot_y(electron_count);
  vec expanded_gamma_z_tot_z(electron_count);
  vec atoms = atomMatrix.col(0);

  mat bonding_params(electron_count, electron_count);
  mat h_off_diag(electron_count, electron_count, fill::ones);
  h_off_diag -= diagmat(vec(electron_count, fill::ones));
  vec density_tot(electron_count);
  mat g_off_diag(electron_count, electron_count, fill::ones);
  g_off_diag -= diagmat(vec(electron_count, fill::ones));
  vec density_gamma_off_diag(electron_count);
  vec density_tot_by_atom(atoms.n_elem);

  int k = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    expanded_diag_gamma_z_x.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(diag_gamma_z_x(i));
    expanded_diag_gamma_z_y.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(diag_gamma_z_y(i));
    expanded_diag_gamma_z_z.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(diag_gamma_z_z(i));
    expanded_gamma_z_tot_x.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(gamma_z_tot_x(i));
    expanded_gamma_z_tot_y.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(gamma_z_tot_y(i));
    expanded_gamma_z_tot_y.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(gamma_z_tot_z(i));
    bonding_params.rows(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) += mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], electron_count).fill(ATOMIC_BONDING_PARAMETERS[atoms(i) - 1]);
    bonding_params.cols(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1) += mat(electron_count, VALENCE_ATOMIC_NUM[atoms(i) - 1]).fill(ATOMIC_BONDING_PARAMETERS[atoms(i) - 1]);
    float density_a = sum(diagvec(p_tot.submat(k, k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1)));
    density_tot.subvec(k, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1).fill(density_a);
    density_tot_by_atom(i) = density_a;
    std::cout << "11i: " << i << std::endl;

    k += VALENCE_ATOMIC_NUM[atoms(i) - 1];
    std::cout << "kkkk" << std::endl;
  }

  mat gamma_density_sum(3, pow(electron_count, 2));

  std::cout << x_gamma - diagmat(diagvec(x_gamma)) << std::endl;
  std::cout << density_tot_by_atom << std::endl;

  gamma_density_sum.row(0) = (x_gamma - diagmat(diagvec(x_gamma))) * density_tot_by_atom;
  gamma_density_sum.row(0) = (y_gamma - diagmat(diagvec(y_gamma))) * density_tot_by_atom;
  gamma_density_sum.row(0) = (z_gamma - diagmat(diagvec(z_gamma))) * density_tot_by_atom;

  // neeed to make shapes compatible
  mat gamma_density_sum_expanded = broadcastToOrbitals(gamma_density_sum);

  vec g_diag_x = gamma_density_sum_expanded.row(0);
  vec g_diag_y = gamma_density_sum_expanded.row(1);
  vec g_diag_z = gamma_density_sum_expanded.row(2);
  // g_off_diag = g_off_diag % (-1.0 * p_a % gamma_expanded);

  // h_off_diag %= (-1.0 / 2.0 * bonding_params % overlapMatrixRA());

  vec h_diag_x = -1.0 * expanded_gamma_z_tot_x;
  vec h_diag_y = -1.0 * expanded_gamma_z_tot_y;
  vec h_diag_z = -1.0 * expanded_gamma_z_tot_z;

  // diagmat(h_diag_x) + diagmat(g_diag_z) + h_off_diag + g_off_diag;

  mat result(3, pow(electron_count, 2));

  result.row(0) = vectorise(diagmat(h_diag_x) + diagmat(g_diag_x));
  result.row(1) = vectorise(diagmat(h_diag_y) + diagmat(g_diag_y));
  result.row(2) = vectorise(diagmat(h_diag_z) + diagmat(g_diag_z));

  return result;
}

mat Cluster::overlapMatrixRA()
{
  mat basisfns = basisFunctions();
  int b = countBasisFunctions();

  mat result(3, pow(b, 2));
  for (int m = 0; m < b; m++)
  {
    for (int v = 0; v < b; v++)
    {
      vec dk_m = basisfns.row(m).cols(3 * K + 1, 4 * K).t();
      rowvec dk_v = basisfns.row(v).cols(3 * K + 1, 4 * K);
      vec nk_m = basisfns.row(m).cols(4 * K + 1, 5 * K).t();
      rowvec nk_v = basisfns.row(v).cols(4 * K + 1, 5 * K);

      for (int k = 0; k < K; k++)
      {
        for (int l = 0; l < K; l++)
        {
          Shell s_kx(basisfns(m, 0), basisfns.row(m).cols(1, 1).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 1, K + 1));
          Shell s_lx(basisfns(v, 0), basisfns.row(v).cols(1, 1).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 1, K + 1));
          ShellOverlapIntegral s_klx(s_kx, s_lx);

          Shell s_ky(basisfns(m, 0), basisfns.row(m).cols(2, 2).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 2, K + 2));
          Shell s_ly(basisfns(v, 0), basisfns.row(v).cols(2, 2).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 2, K + 2));
          ShellOverlapIntegral s_kly(s_ky, s_ly);

          Shell s_kz(basisfns(m, 0), basisfns.row(m).cols(3, 3).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 3, K + 3));
          Shell s_lz(basisfns(v, 0), basisfns.row(v).cols(3, 3).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 3, K + 3));
          ShellOverlapIntegral s_klz(s_kz, s_lz);

          result(0, m * b + v) += dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l) * sum(sum(s_klx.saDerivativeIntegral() * s_kly() * s_klz()));
          result(1, m * b + v) += dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l) * sum(sum(s_klx() * s_kly.saDerivativeIntegral() * s_klz()));
          result(2, m * b + v) += dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l) * sum(sum(s_klx() * s_kly() * s_klz.saDerivativeIntegral()));
        }
      }
    }
  }

  return result;
}

mat Cluster::overlapMatrix()
{
  mat basisfns = basisFunctions();
  int b = countBasisFunctions();

  mat result(b, b);
  for (int m = 0; m < b; m++)
  {
    for (int v = 0; v < b; v++)
    {
      vec dk_m = basisfns.row(m).cols(3 * K + 1, 4 * K).t();
      rowvec dk_v = basisfns.row(v).cols(3 * K + 1, 4 * K);
      vec nk_m = basisfns.row(m).cols(4 * K + 1, 5 * K).t();
      rowvec nk_v = basisfns.row(v).cols(4 * K + 1, 5 * K);

      for (int k = 0; k < K; k++)
      {
        for (int l = 0; l < K; l++)
        {
          Shell s_kx(basisfns(m, 0), basisfns.row(m).cols(1, 1).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 1, K + 1));
          Shell s_lx(basisfns(v, 0), basisfns.row(v).cols(1, 1).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 1, K + 1));
          ShellOverlapIntegral s_klx(s_kx, s_lx);

          Shell s_ky(basisfns(m, 0), basisfns.row(m).cols(2, 2).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 2, K + 2));
          Shell s_ly(basisfns(v, 0), basisfns.row(v).cols(2, 2).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 2, K + 2));
          ShellOverlapIntegral s_kly(s_ky, s_ly);

          Shell s_kz(basisfns(m, 0), basisfns.row(m).cols(3, 3).t(), basisfns.row(m).cols(3 * K + 1 + k, 3 * K + 1 + k).t(), basisfns.row(m).cols(2 * K + 1 + k, 2 * K + 1 + k).t(), basisfns.row(m).cols(K + 3, K + 3));
          Shell s_lz(basisfns(v, 0), basisfns.row(v).cols(3, 3).t(), basisfns.row(v).cols(3 * K + 1 + l, 3 * K + 1 + l).t(), basisfns.row(v).cols(2 * K + 1 + l, 2 * K + 1 + l).t(), basisfns.row(v).cols(K + 3, K + 3));
          ShellOverlapIntegral s_klz(s_kz, s_lz);

          result(m, v) += dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l) * sum(sum(s_klx() * s_kly() * s_klz()));
        }
      }
    }
  }

  return result;
}

mat Cluster::gammaMatrixRA()
{
  mat basisfns = sBasisFunctions();
  int atomCount = atomMatrix.n_rows;
  mat gammas(3, pow(atomCount, 2));

  for (int m = 0; m < atomCount; m++)
  {
    for (int v = 0; v < atomCount; v++)
    {
      vec dk_m = basisfns.row(m).cols(3 * K + 1, 4 * K).t();
      vec dk_v = basisfns.row(v).cols(3 * K + 1, 4 * K).t();
      vec nk_m = basisfns.row(m).cols(4 * K + 1, 5 * K).t();
      vec nk_v = basisfns.row(v).cols(4 * K + 1, 5 * K).t();

      vec dk_m_prime = dk_m % nk_m;
      vec dk_v_prime = dk_v % nk_v;

      mat dk_m_prime_combos = dk_m_prime * dk_m_prime.t();
      mat dk_v_prime_combos = dk_v_prime * dk_v_prime.t();

      vec dk_m_prime_combos_vec = arma::vectorise(dk_m_prime_combos.t());
      vec dk_v_prime_combos_vec = arma::vectorise(dk_v_prime_combos.t());

      float r_12 = calcDistance(atomMatrix.row(m).cols(1, 3), atomMatrix.row(v).cols(1, 3));

      mat gamma_summation_matrix = dk_m_prime_combos_vec * dk_v_prime_combos_vec.t();

      vec alphas = basisfns.row(m).cols(2 * K + 1, 3 * K).t();
      mat alpha_copy(alphas.n_elem, alphas.n_elem, arma::fill::ones);
      alpha_copy.each_col() %= alphas;

      mat alpha_prime = alpha_copy + alpha_copy.t();

      vec betas = basisfns.row(v).cols(2 * K + 1, 3 * K).t();
      mat beta_copy(betas.n_elem, betas.n_elem, arma::fill::ones);
      beta_copy.each_col() %= betas;

      mat beta_prime = beta_copy + beta_copy.t();

      vec alpha_prime_vec = arma::vectorise(alpha_prime, 1).t();
      vec beta_prime_vec = arma::vectorise(beta_prime, 1).t();

      vec sigma_a = arma::pow(alpha_prime_vec, -1.0);
      mat sigma_a_copy(sigma_a.n_elem, sigma_a.n_elem, arma::fill::ones);
      sigma_a_copy.each_col() %= sigma_a;
      vec u_a = arma::pow(sigma_a * M_PI, 3.0 / 2.0);
      vec sigma_b = arma::pow(beta_prime_vec, -1.0);
      mat sigma_b_copy(sigma_a.n_elem, sigma_a.n_elem, arma::fill::ones);
      sigma_b_copy.each_col() %= sigma_b;
      vec u_b = arma::pow(sigma_b * M_PI, 3.0 / 2.0);
      mat v_sq = arma::pow(sigma_a_copy + sigma_b_copy.t(), -1.0);
      mat u = u_a * u_b.t();
      float r_ab_dist = calcDistance(atomMatrix.row(m).cols(1, 3), atomMatrix.row(v).cols(1, 3));
      mat t = v_sq * pow(r_ab_dist, 2);
      mat erf_t = v_sq * pow(r_ab_dist, 2);
      erf_t.for_each([](mat::elem_type &val)
                     { val = erf(sqrt(val)); });

      if (m == v)
      {
        gammas(0, m * atomCount + v) = 0.0;
        gammas(1, m * atomCount + v) = 0.0;
        gammas(2, m * atomCount + v) = 0.0;
      }
      else
      {
        gammas(0, m * atomCount + v) = CONVERSION_FACTOR * accu(gamma_summation_matrix % ((u / pow(r_ab_dist, 2)) % ((2 * sqrt(v_sq / M_PI)) % exp(-1.0 * t) - (erf_t / abs(r_ab_dist))) * (atomMatrix.row(m).cols(1, 3)(0) - atomMatrix.row(v).cols(1, 3)(0))));
        gammas(1, m * atomCount + v) = CONVERSION_FACTOR * accu(gamma_summation_matrix % ((u / pow(r_ab_dist, 2)) % ((2 * sqrt(v_sq / M_PI)) % exp(-1.0 * t) - (erf_t / abs(r_ab_dist))) * (atomMatrix.row(m).cols(1, 3)(1) - atomMatrix.row(v).cols(1, 3)(1))));
        gammas(2, m * atomCount + v) = CONVERSION_FACTOR * accu(gamma_summation_matrix % ((u / pow(r_ab_dist, 2)) % ((2 * sqrt(v_sq / M_PI)) % exp(-1.0 * t) - (erf_t / abs(r_ab_dist))) * (atomMatrix.row(m).cols(1, 3)(2) - atomMatrix.row(v).cols(1, 3)(2))));
      }
    }
  }

  return gammas;
}

mat Cluster::x(mat p_a, mat p_b, mat bonding_params)
{
  mat p_tot = p_a + p_b;
  p_tot %= bonding_params;
  return p_tot;
}

mat Cluster::gammaMatrix()
{
  mat basisfns = sBasisFunctions();
  int atomCount = atomMatrix.n_rows;
  mat gammas(atomCount, atomCount);
  mat result(atomCount, atomCount);

  for (int m = 0; m < atomCount; m++)
  {
    for (int v = 0; v < atomCount; v++)
    {
      vec dk_m = basisfns.row(m).cols(3 * K + 1, 4 * K).t();
      vec dk_v = basisfns.row(v).cols(3 * K + 1, 4 * K).t();
      vec nk_m = basisfns.row(m).cols(4 * K + 1, 5 * K).t();
      vec nk_v = basisfns.row(v).cols(4 * K + 1, 5 * K).t();

      vec dk_m_prime = dk_m % nk_m;
      vec dk_v_prime = dk_v % nk_v;

      mat dk_m_prime_combos = dk_m_prime * dk_m_prime.t();
      mat dk_v_prime_combos = dk_v_prime * dk_v_prime.t();

      vec dk_m_prime_combos_vec = arma::vectorise(dk_m_prime_combos.t());
      vec dk_v_prime_combos_vec = arma::vectorise(dk_v_prime_combos.t());

      float r_12 = calcDistance(atomMatrix.row(m).cols(1, 3), atomMatrix.row(v).cols(1, 3));

      mat gamma_summation_matrix = dk_m_prime_combos_vec * dk_v_prime_combos_vec.t();

      vec alphas = basisfns.row(m).cols(2 * K + 1, 3 * K).t();
      mat alpha_copy(alphas.n_elem, alphas.n_elem, arma::fill::ones);
      alpha_copy.each_col() %= alphas;

      mat alpha_prime = alpha_copy + alpha_copy.t();

      vec betas = basisfns.row(v).cols(2 * K + 1, 3 * K).t();
      mat beta_copy(betas.n_elem, betas.n_elem, arma::fill::ones);
      beta_copy.each_col() %= betas;

      mat beta_prime = beta_copy + beta_copy.t();

      vec alpha_prime_vec = arma::vectorise(alpha_prime, 1).t();
      vec beta_prime_vec = arma::vectorise(beta_prime, 1).t();

      vec sigma_a = arma::pow(alpha_prime_vec, -1.0);
      mat sigma_a_copy(sigma_a.n_elem, sigma_a.n_elem, arma::fill::ones);
      sigma_a_copy.each_col() %= sigma_a;
      vec u_a = arma::pow(sigma_a * M_PI, 3.0 / 2.0);
      vec sigma_b = arma::pow(beta_prime_vec, -1.0);
      mat sigma_b_copy(sigma_a.n_elem, sigma_a.n_elem, arma::fill::ones);
      sigma_b_copy.each_col() %= sigma_b;
      vec u_b = arma::pow(sigma_b * M_PI, 3.0 / 2.0);
      mat v_sq = arma::pow(sigma_a_copy + sigma_b_copy.t(), -1.0);
      mat u = u_a * u_b.t();
      float r_ab_dist = calcDistance(atomMatrix.row(m).cols(1, 3), atomMatrix.row(v).cols(1, 3));
      mat t = v_sq * pow(r_ab_dist, 2);
      t.for_each([](mat::elem_type &val)
                 { val = erf(sqrt(val)); });

      if (m == v)
      {
        gammas(m, v) = CONVERSION_FACTOR * sum(sum(gamma_summation_matrix % (u % sqrt(2.0 * v_sq) * sqrt(2.0 / M_PI))));
      }
      else
      {
        gammas(m, v) += CONVERSION_FACTOR * sum(sum(gamma_summation_matrix % (u % t * sqrt(1.0 / pow(r_ab_dist, 2)))));
      }
    }
  }

  return gammas;
}

void Cluster::calcSCFEnergyRA(float threshold)
{
  bool mats_eq = false;

  mat z = z_vals();
  int electron_count = accu(z);

  // store the alpha and beta density matrices
  mat p_a(electron_count, electron_count, fill::zeros);
  mat p_b(electron_count, electron_count, fill::zeros);

  // store the old alpha and beta density matrices
  mat p_a_old(electron_count, electron_count, fill::zeros);
  mat p_b_old(electron_count, electron_count, fill::zeros);

  // store the output fo the cndo2 fock matrix calculation
  mat f_a(electron_count, electron_count, fill::zeros);
  mat f_b(electron_count, electron_count, fill::zeros);

  // store the eigenvectors
  mat c_a(electron_count, electron_count, fill::zeros);
  mat c_b(electron_count, electron_count, fill::zeros);

  // store the eigenvalues
  vec e_a(electron_count, fill::zeros);
  vec e_b(electron_count, fill::zeros);

  mat x(electron_count, electron_count);
  mat y(atomMatrix.n_rows, atomMatrix.n_rows);

  while (!mats_eq)
  {
    f_a = cndo2FockMatrixRA(p_a, p_b);
    f_b = cndo2FockMatrixRA(p_b, p_a);

    p_a_old = p_a;
    p_b_old = p_b;

    eig_sym(e_a, c_a, f_a);
    eig_sym(e_b, c_b, f_b);

    p_a = (c_a.cols(0, _p - 1) * c_a.cols(0, _p - 1).t());
    p_b = (c_b.cols(0, _q - 1) * c_b.cols(0, _q - 1).t());

    mats_eq = approx_equal(p_a_old, p_a, "absdiff", threshold) && approx_equal(p_b_old, p_b, "absdiff", threshold);
  }

  mat h_core = cndo2FockMatrixRA(mat(electron_count, electron_count, fill::zeros), mat(electron_count, electron_count, fill::zeros));

  x = (p_a + p_b);

  float electron_energy = (1.0 / 2.0) * sum(sum((p_a % (h_core + f_a) + p_b % (h_core + f_b))));

  // my original goal here was to eliminate the loop and do a sum of the triangular matrix but i couldnt figure out how to do that with the atom distances and my brain was tired
  float nuc_repulsion = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    for (int j = 0; j < i; j++)
    {
      nuc_repulsion += (CONVERSION_FACTOR * ((VALENCE_ATOMIC_NUM[atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[atomMatrix(j, 0) - 1]) / calcDistance(atomMatrix.row(i).cols(1, 3), atomMatrix.row(j).cols(1, 3))));
    }
  }

  std::cout << "Nuclear Repulsion Energy is " << nuc_repulsion << " eV." << std::endl;
  std::cout << "Electronic Energy is " << electron_energy << " eV." << std::endl;
  std::cout << "Total Energy is " << electron_energy + nuc_repulsion << " eV." << std::endl;
}

void Cluster::calcSCFEnergy(float threshold)
{
  bool mats_eq = false;

  vec atoms = atomMatrix.col(0);

  mat z = z_vals();
  int electron_count = sum(sum(z));

  // store the alpha and beta density matrices
  mat p_a(electron_count, electron_count, fill::zeros);
  mat p_b(electron_count, electron_count, fill::zeros);

  // store the old alpha and beta density matrices
  mat p_a_old(electron_count, electron_count, fill::zeros);
  mat p_b_old(electron_count, electron_count, fill::zeros);

  // store the output fo the cndo2 fock matrix calculation
  mat f_a(electron_count, electron_count, fill::zeros);
  mat f_b(electron_count, electron_count, fill::zeros);

  // store the eigenvectors
  mat c_a(electron_count, electron_count, fill::zeros);
  mat c_b(electron_count, electron_count, fill::zeros);

  // store the eigenvalues
  vec e_a(electron_count, fill::zeros);
  vec e_b(electron_count, fill::zeros);

  // matrix that stores the bonding param sum
  mat bonding_params(electron_count, electron_count, fill::zeros);

  while (!mats_eq)
  {

    bonding_params = mat(electron_count, electron_count, fill::zeros);

    f_a = cndo2FockMatrix(p_a, p_b, bonding_params);

    bonding_params = mat(electron_count, electron_count, fill::zeros);
    f_b = cndo2FockMatrix(p_b, p_a, bonding_params);

    p_a_old = p_a;
    p_b_old = p_b;

    eig_sym(e_a, c_a, f_a);
    eig_sym(e_b, c_b, f_b);

    p_a = (c_a.cols(0, _p - 1) * c_a.cols(0, _p - 1).t());
    p_b = (c_b.cols(0, _q - 1) * c_b.cols(0, _q - 1).t());

    mats_eq = approx_equal(p_a_old, p_a, "absdiff", threshold) && approx_equal(p_b_old, p_b, "absdiff", threshold);
  }

  bonding_params = mat(electron_count, electron_count, fill::zeros);
  mat h_core = cndo2FockMatrix(mat(electron_count, electron_count, fill::zeros), mat(electron_count, electron_count, fill::zeros), bonding_params);

  float electron_energy = (1.0 / 2.0) * sum(sum((p_a % (h_core + f_a) + p_b % (h_core + f_b))));

  vec p_a_tot(atomMatrix.n_rows);
  vec p_b_tot(atomMatrix.n_rows);

  int idx = 0;

  vec valence_atomic_num_vec(electron_count);

  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    p_a_tot(i) = sum(diagvec(p_a.submat(idx, idx, idx + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1)));
    p_b_tot(i) = sum(diagvec(p_b.submat(idx, idx, idx + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1)));
    valence_atomic_num_vec(i) = VALENCE_ATOMIC_NUM[atoms(i) - 1];

    idx += VALENCE_ATOMIC_NUM[atoms(i) - 1];
  }

  vec p_tot = p_a_tot + p_b_tot;

  mat y(atomMatrix.n_rows, atomMatrix.n_rows);

  int a_idx = 0;
  int b_idx = 0;

  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    b_idx = 0;
    for (int j = 0; j < atomMatrix.n_rows; j++)
    {
      y(i, j) = p_tot(i) * p_tot(j) - valence_atomic_num_vec(j) * p_tot(i) - valence_atomic_num_vec(i) * p_tot(j) - accu(pow(p_a.submat(a_idx, b_idx, a_idx + valence_atomic_num_vec(i) - 1, b_idx + valence_atomic_num_vec(j) - 1), 2) + pow(p_a.submat(a_idx, b_idx, a_idx + valence_atomic_num_vec(i) - 1, b_idx + valence_atomic_num_vec(j) - 1), 2));
      b_idx += VALENCE_ATOMIC_NUM[atoms(j) - 1];
    }
    a_idx += VALENCE_ATOMIC_NUM[atoms(i) - 1];
  }

  mat gradient_electron(3, atomMatrix.n_rows, fill::zeros);
  mat gradient_nuc(3, atomMatrix.n_rows, fill::zeros);

  mat ra_overlap = overlapMatrixRA();

  mat ra_gamma = gammaMatrixRA();

  mat xs_mask = mat(electron_count, electron_count, fill::ones) - broadcastToOrbitals(diagmat(vec(atomMatrix.n_rows, fill::ones)));

  mat density_term_expanded_x = xs_mask % x(p_a, p_b, bonding_params) % reshape(mat(ra_overlap.row(0)), electron_count, electron_count);
  mat density_term_x(atomMatrix.n_rows, atomMatrix.n_rows);
  mat density_term_expanded_y = xs_mask % x(p_a, p_b, bonding_params) % reshape(mat(ra_overlap.row(1)), electron_count, electron_count);
  mat density_term_y(atomMatrix.n_rows, atomMatrix.n_rows);
  mat density_term_expanded_z = xs_mask % x(p_a, p_b, bonding_params) % reshape(mat(ra_overlap.row(2)), electron_count, electron_count);
  mat density_term_z(atomMatrix.n_rows, atomMatrix.n_rows);

  int idx_i = 0;
  int idx_j = 0;

  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    int idx_j = 0;
    for (int j = 0; j < atomMatrix.n_rows; j++)
    {
      density_term_x(i, j) = accu(density_term_expanded_x.submat(idx_i, idx_j, idx_i + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx_j + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1));
      density_term_y(i, j) = accu(density_term_expanded_y.submat(idx_i, idx_j, idx_i + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx_j + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1));
      density_term_z(i, j) = accu(density_term_expanded_z.submat(idx_i, idx_j, idx_i + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, idx_j + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1));

      idx_j += VALENCE_ATOMIC_NUM[atoms(j) - 1];
    }
    idx_i += VALENCE_ATOMIC_NUM[atoms(i) - 1];
  }

  gradient_electron.row(0) = sum(density_term_x) + sum((y - diagmat(diagvec(y))) % reshape(mat(ra_gamma.row(0)), atomMatrix.n_rows, atomMatrix.n_rows));
  gradient_electron.row(1) = sum(density_term_y) + sum((y - diagmat(diagvec(y))) % reshape(mat(ra_gamma.row(1)), atomMatrix.n_rows, atomMatrix.n_rows));
  gradient_electron.row(2) = sum(density_term_z) + sum((y - diagmat(diagvec(y))) % reshape(mat(ra_gamma.row(2)), atomMatrix.n_rows, atomMatrix.n_rows));

  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    for (int j = 0; j < atomMatrix.n_rows; j++)
    {
      if (j != i)
      {
        gradient_nuc.row(0).col(i) += (CONVERSION_FACTOR * -1 * (VALENCE_ATOMIC_NUM[atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[atomMatrix(j, 0) - 1]) * (atomMatrix.row(i).col(1) - atomMatrix.row(j).col(1)) * pow(calcDistance(atomMatrix.row(i).cols(1, 3), atomMatrix.row(j).cols(1, 3)), -3));
        gradient_nuc.row(1).col(i) += (CONVERSION_FACTOR * -1 * (VALENCE_ATOMIC_NUM[atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[atomMatrix(j, 0) - 1]) * (atomMatrix.row(i).col(2) - atomMatrix.row(j).col(2)) * pow(calcDistance(atomMatrix.row(i).cols(1, 3), atomMatrix.row(j).cols(1, 3)), -3));
        gradient_nuc.row(2).col(i) += (CONVERSION_FACTOR * -1 * (VALENCE_ATOMIC_NUM[atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[atomMatrix(j, 0) - 1]) * (atomMatrix.row(i).col(3) - atomMatrix.row(j).col(3)) * pow(calcDistance(atomMatrix.row(i).cols(1, 3), atomMatrix.row(j).cols(1, 3)), -3));
      }
    }
  }

  std::cout << "gradient (Nuclear part)" << std::endl
            << gradient_nuc << std::endl;

  std::cout << "gradient (Electron part)" << std::endl
            << gradient_electron << std::endl;

  std::cout << "gradient" << std::endl
            << gradient_electron + gradient_nuc << std::endl;

  // my original goal here was to eliminate the loop and do a sum of the triangular matrix but i couldnt figure out how to do that with the atom distances and my brain was tired
  float nuc_repulsion = 0;
  for (int i = 0; i < atomMatrix.n_rows; i++)
  {
    for (int j = 0; j < i; j++)
    {
      nuc_repulsion += (CONVERSION_FACTOR * ((VALENCE_ATOMIC_NUM[atomMatrix(i, 0) - 1] * VALENCE_ATOMIC_NUM[atomMatrix(j, 0) - 1]) / calcDistance(atomMatrix.row(i).cols(1, 3), atomMatrix.row(j).cols(1, 3))));
    }
  }

  std::cout << "Nuclear Repulsion Energy is " << nuc_repulsion << " eV." << std::endl;
  std::cout << "Electronic Energy is " << electron_energy << " eV." << std::endl;
  std::cout << "Total Energy is " << electron_energy + nuc_repulsion << " eV." << std::endl;
}

mat Cluster::broadcastToOrbitals(mat x)
{
  mat z = z_vals();
  int electron_count = sum(sum(z));
  vec atoms = atomMatrix.col(0);

  mat orbitalRepresentation(electron_count, electron_count);
  int k = 0;
  int l = 0;
  for (int i = 0; i < x.n_rows; i++)
  {
    l = 0;

    for (int j = 0; j < x.n_cols; j++)
    {
      orbitalRepresentation.submat(k, l, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, l + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1) = mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], VALENCE_ATOMIC_NUM[atoms(j) - 1]).fill(x(i, j));
      l += VALENCE_ATOMIC_NUM[atoms(j) - 1];
    }
    k += VALENCE_ATOMIC_NUM[atoms(i) - 1];
  }
  return orbitalRepresentation;
}

mat Cluster::broadcastToOrbitalsRA(mat x)
{
  mat z = z_vals();
  int electron_count = accu(z);
  vec atoms = atomMatrix.col(0);

  mat xOrbitalRepresentation(electron_count, electron_count);
  mat yOrbitalRepresentation(electron_count, electron_count);
  mat zOrbitalRepresentation(electron_count, electron_count);

  mat orbitalRepresentation(3, pow(electron_count, 2));

  int k = 0;
  int l = 0;

  for (int i = 0; i < sqrt(x.n_cols); i++)
  {
    l = 0;

    for (int j = 0; j < sqrt(x.n_cols); j++)
    {
      xOrbitalRepresentation.submat(k, l, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, l + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1) = mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], VALENCE_ATOMIC_NUM[atoms(j) - 1]).fill(x(0, i * sqrt(electron_count) + j));
      yOrbitalRepresentation.submat(k, l, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, l + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1) = mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], VALENCE_ATOMIC_NUM[atoms(j) - 1]).fill(x(1, i * sqrt(electron_count) + j));
      zOrbitalRepresentation.submat(k, l, k + VALENCE_ATOMIC_NUM[atoms(i) - 1] - 1, l + VALENCE_ATOMIC_NUM[atoms(j) - 1] - 1) = mat(VALENCE_ATOMIC_NUM[atoms(i) - 1], VALENCE_ATOMIC_NUM[atoms(j) - 1]).fill(x(2, i * sqrt(electron_count) + j));

      l += VALENCE_ATOMIC_NUM[atoms(j) - 1];
    }
    k += VALENCE_ATOMIC_NUM[atoms(i) - 1];
  }

  std::cout << xOrbitalRepresentation << std::endl;
  orbitalRepresentation.row(0) = vectorise(xOrbitalRepresentation, 1);
  orbitalRepresentation.row(1) = vectorise(yOrbitalRepresentation, 1);
  orbitalRepresentation.row(2) = vectorise(zOrbitalRepresentation, 1);

  return orbitalRepresentation;
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

double Cluster::calcDistance(mat a1, mat a2)
{
  double distance = 0;
  rowvec coordDists = square(a2.row(0) - a1.row(0));
  return sqrt(sum(coordDists));
}

void Cluster::addAtom(int index, int atomNum, double x, double y, double z, double e, double s)
{
  atomMatrix(index, 0) = atomNum;
  atomMatrix(index, 1) = x;
  atomMatrix(index, 2) = y;
  atomMatrix(index, 3) = z;
  epsilons(index) = e;
  sigmas(index) = s;
  int numElectrons = sum(z_vals());
  _p = numElectrons / 2 + numElectrons % 2;
  _q = numElectrons / 2;
}

std::ostream &operator<<(std::ostream &os, const Cluster &c)
{
  os << "Atom Matrix: \n"
     << c.atomMatrix << "\n";

  return os;
}
