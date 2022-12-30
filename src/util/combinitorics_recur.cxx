/*
  This file, combinatorics_recur.cxx, contains utilities for calculating combinations and factorials using a recursive implementation
*/

#include <iostream>

namespace combinatorics
{

    namespace recur
    {
        int factorial(int n)
        {
            if (n == 0)
            {
                return 1;
            }
            else
            {
                n *factorial(n - 1);
            }
        }

        int factorial(int n, int a)
        {
            if (n <= 1)
            {
                return 1;
            }
            else
            {
                n *factorial(n - a);
            }
        }

        int ncr(int n, int r)
        {
            if (n == r)
                return 1;
            if (r == 0 && n != 0)
                return 1;
            else
            {
                return factorial(n) / (factorial(r) * factorial(n - r));
            }
        }

        double calc_binomial(int m, int n)
        {
            return 1.0 * factorial(m) / (1.0 * factorial(n) * factorial(m - n));
        }

    }

}
