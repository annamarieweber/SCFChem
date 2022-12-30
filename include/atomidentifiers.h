/**
 * @file atomidentifiers.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief Defines constants for mapping atom symbols to atomic number
 * @version 0.1
 * @date 2022-12-26
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef ATOM_ID_CONSTANTS_H
#define ATOM_ID_CONSTANTS_H

#include <vector>
#include <string>
#include <armadillo>

using arma::cube;
using arma::mat;
using arma::vec;
using std::map;
using std::string;

namespace constants
{
    namespace atomidentifiers
    {
        const map<string, int> ATOM_NUM = {
            {"H", 1},
            {"C", 6},
            {"N", 7},
            {"O", 8},
        };
    }
}

#endif