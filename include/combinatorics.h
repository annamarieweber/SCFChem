/**
 * @file combinatorics_iter.h
 * @author Anna Weber (anna@scf.com)
 * @brief This file contains an API to all iterative combinatorics functions
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

/**
 * @brief combinatorics functions
 *
 */
namespace combo
{

    /**
     * @brief Iterative combinatorics implementation.
     */

    namespace iter
    {
        /**
         * @brief Calculates the factorial of a given integer. Uses an iterative implementation.
         * @param n The integer to calculate the factorial for.
         * @return The factorial of n.
         */
        int factorial(int n);

        /**
         * @brief Calculates the factorial of a given integer where the output is the product of all integers from 1 to n where x%a == n%a. Uses an iterative implementation.
         * @param n The integer to calculate the factorial for.
         * @param a The modulo value.
         * @return The factorial of n.
         */
        int factorial(int n, int a);

        /**
         * @brief Calculates the number of combinations of r items from a set of n items. Uses an iterative implementation.
         * @param n The number of items in the set.
         * @param r The number of items to choose.
         * @return The number of combinations.
         */
        int ncr(int n, int r);

        /**
         * @brief Calculates the binomial coefficient for the given values of m and n. Uses an iterative implementation.
         * @param m The first value.
         * @param n The second value.
         * @return The binomial coefficient.
         */
        double calc_binomial(int m, int n);

        /**
         * @brief Calculates the factorial of a given vector of integers. Uses an iterative implementation.
         * @param n The vector of integers to calculate the factorial for.
         * @param a The modulo value.
         * @return The factorial of n.
         */
        vec vec_factorial(int n, int a);

        /**
         * @brief Calculates the factorial of a given vector of integers. Uses an iterative implementation.
         * @param n The vector of integers to calculate the factorial for.
         * @return The factorial of n.
         */
        vec vec_factorial(int n);
        /**
         * @brief Calculates the binomial coefficient for the given vectors of m and n. Uses an iterative implementation.
         * @param m The first vector.
         * @param n The second vector.
         * @return The binomial coefficient.
         */
        vec calc_binomial(vec m, vec n);
    }

    /**
     * @brief Recursive combinatorics implementation
     *
     */
    namespace recur
    {
        /**
         * @brief Calculates the factorial of a given integer. Uses a recursive implementation.
         * @param n The integer to calculate the factorial for.
         * @return The factorial of n.
         */
        int factorial(int n);

        /**
         * @brief Calculates the factorial of a given integer where the output is the product of all integers from 1 to n where x%a == n%a. Uses a recursive implementation.
         * @param n The integer to calculate the factorial for.
         * @param a The modulo value.
         * @return The factorial of n.
         */
        int factorial(int n, int a);

        /**
         * @brief Calculates the number of combinations of r items from a set of n items. Uses a recursive implementation.
         * @param n The number of items in the set.
         * @param r The number of items to choose.
         * @return The number of combinations.
         */
        int ncr(int n, int r);

        /**
         * @brief Calculates the binomial coefficient for the given values of m and n. Uses a recursive implementation.
         * @param m The first value.
         * @param n The second value.
         * @return The binomial coefficient.
         */
        double calc_binomial(int m, int n);

        /**
         * @brief Calculates the factorial of a given vector of integers. Uses a recursive implementation.
         * @param n The vector of integers to calculate the factorial for.
         * @param a The modulo value.
         * @return The factorial of n.
         */
        vec vec_factorial(int n, int a);

        /**
         * @brief Calculates the factorial of a given vector of integers. Uses a recursive implementation.
         * @param n The vector of integers to calculate the factorial for.
         * @return The factorial of n.
         */
        vec vec_factorial(int n);
        /**
         * @brief Calculates the binomial coefficient for the given vectors of m and n. Uses a recursive implementation.
         * @param m The first vector.
         * @param n The second vector.
         * @return The binomial coefficient.
         */
        vec calc_binomial(vec m, vec n);

    }

}
#endif
