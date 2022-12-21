#include <armadillo>
#include "cluster.h"
#include "algoutils.h"
#include "molutils.h"
#include "shell.h"
#include "shelloverlapintegral.h"

using namespace arma;

mat algo::utils::overlapMatrixRA(Cluster cluster)
{
  int K = cluster.K;
  int numBasisFunctions = cluster.numBasisFunctions;
  mat basisfns = cluster.basisFunctions();

  mat result(3, pow(numBasisFunctions, 2));
  for (int m = 0; m < numBasisFunctions; m++)
  {
    for (int v = 0; v < numBasisFunctions; v++)
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

          result(0, m * numBasisFunctions + v) += dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l) * sum(sum(s_klx.saDerivativeIntegral() * s_kly() * s_klz()));
          result(1, m * numBasisFunctions + v) += dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l) * sum(sum(s_klx() * s_kly.saDerivativeIntegral() * s_klz()));
          result(2, m * numBasisFunctions + v) += dk_m(k) * dk_v(l) * nk_m(k) * nk_v(l) * sum(sum(s_klx() * s_kly() * s_klz.saDerivativeIntegral()));
        }
      }
    }
  }

  return result;
}

mat algo::utils::overlapMatrix(Cluster cluster)
{
  int K = cluster.K;
  int numBasisFunctions = cluster.numBasisFunctions;
  mat basisfns = cluster.basisFunctions();

  mat result(numBasisFunctions, numBasisFunctions);
  for (int m = 0; m < numBasisFunctions; m++)
  {
    for (int v = 0; v < numBasisFunctions; v++)
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

mat algo::utils::gammaMatrixRA(Cluster cluster)
{
  int K = cluster.K;
  int numBasisFunctions = cluster.numBasisFunctions;
  mat basisfns = cluster.sBasisFunctions();
  int atomCount = cluster.atomMatrix.n_rows;
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

      float r_12 = calcDistance(cluster.atomMatrix.row(m).cols(1, 3), cluster.atomMatrix.row(v).cols(1, 3));

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
      float r_ab_dist = calcDistance(cluster.atomMatrix.row(m).cols(1, 3), cluster.atomMatrix.row(v).cols(1, 3));
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
        gammas(0, m * atomCount + v) = CONVERSION_FACTOR * accu(gamma_summation_matrix % ((u / pow(r_ab_dist, 2)) % ((2 * sqrt(v_sq / M_PI)) % exp(-1.0 * t) - (erf_t / abs(r_ab_dist))) * (cluster.atomMatrix.row(m).cols(1, 3)(0) - cluster.atomMatrix.row(v).cols(1, 3)(0))));
        gammas(1, m * atomCount + v) = CONVERSION_FACTOR * accu(gamma_summation_matrix % ((u / pow(r_ab_dist, 2)) % ((2 * sqrt(v_sq / M_PI)) % exp(-1.0 * t) - (erf_t / abs(r_ab_dist))) * (cluster.atomMatrix.row(m).cols(1, 3)(1) - cluster.atomMatrix.row(v).cols(1, 3)(1))));
        gammas(2, m * atomCount + v) = CONVERSION_FACTOR * accu(gamma_summation_matrix % ((u / pow(r_ab_dist, 2)) % ((2 * sqrt(v_sq / M_PI)) % exp(-1.0 * t) - (erf_t / abs(r_ab_dist))) * (cluster.atomMatrix.row(m).cols(1, 3)(2) - cluster.atomMatrix.row(v).cols(1, 3)(2))));
      }
    }
  }

  return gammas;
}

mat algo::utils::gammaMatrix(Cluster cluster)
{

  int K = cluster.K;
  int numBasisFunctions = cluster.numBasisFunctions;
  mat basisfns = cluster.sBasisFunctions();
  int atomCount = cluster.atomMatrix.n_rows;
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

      float r_12 = calcDistance(cluster.atomMatrix.row(m).cols(1, 3), cluster.atomMatrix.row(v).cols(1, 3));

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
      float r_ab_dist = calcDistance(cluster.atomMatrix.row(m).cols(1, 3), cluster.atomMatrix.row(v).cols(1, 3));
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
