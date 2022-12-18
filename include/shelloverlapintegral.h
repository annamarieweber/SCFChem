/**
 * @file shelloverlapintegral.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief definition of electron shell overlap integral
 * @version 0.1
 * @date 2022-12-17
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef SHELL_OVERLAP_INTEGRAL
#define SHELL_OVERLAP_INTEGRAL
#include "shell.h"

using arma::cube;
using std::vector;

class ShellOverlapIntegral
{
private:
    Shell _s_a;
    Shell _s_b;

    vec alphas_product();

    vec alphas_sum();

    vec dim_dist_sqr();

    vec exponential_prefactor();

    vec root_term();

    vec overlap_summation(vec x_p, int l_pair_a, int l_pair_b);

    vec product_center();

public:
    ShellOverlapIntegral();

    ShellOverlapIntegral(Shell s_a, Shell s_b);

    mat operator()();

    mat saDerivativeIntegral();

    mat sbDerivativeIntegral();
};
#endif