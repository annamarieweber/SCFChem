/**
 * @file combinatorics.h
 * @author Anna Weber (anna@scf.com)
 * @brief This file contains an API to all combinatorics functions
 * @version 0.1
 * @date 2022-12-17
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef COMBO
#define COMBO
#include <string>
#include <armadillo>

using arma::vec;

int factorial(int n);

// calculates the a factorial where the output is the product of all integers from 1 to n where x%a == n%a
int factorial(int n, int a);

int ncr(int n, int r);

double calc_binomial(int m, int n);

vec vec_factorial(int n, int a);

vec vec_factorial(int n);

vec calc_binomial(vec m, vec n);
#endif
