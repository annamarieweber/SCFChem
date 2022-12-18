/**
 * @file filereader.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief Definitions for file reading functions for SCFChem
 * @version 0.1
 * @date 2022-12-17
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef FILEREADER
#define FILEREADER
#include <string>
#include <cstdlib>
#include "cluster.h"

enum CompareMethod
{
  LTE,
  EQ
};

/**
 * @brief Validates that the number of atoms counted is <= the expected number of atoms.
 *
 * @param atom_count (int) the number of atoms read
 * @param expected_count (int) the expected number of atoms
 * @param method the method of comparison to use
 * @throws ClusterConfigException(AtomCountMismatch)
 */
void validateAtomCount(int atom_count, int expected_count, CompareMethod method);

/**
 * @brief Converts string identifier used for atoms into atomic number
 *
 * @param atom_id
 * @returns int
 */
int getAtomNum(string atom_id);

/**
 * @brief reads files with molecule data
 *
 * @details reads molecule files where the first line is the number of molecules and
 *    		the remaining lines follow the format below:
 *         	atomic_number x_coord y_coord z_coord
 *
 *         	atomic_numbers must be > 0 and < MAX_ATOMIC_NUM (Set by MAX_ATOMIC_NUM environment variable)
 *        	atomic_numbers can be restricted with ALLOWED_ATOMIC_NUM environment variable which is a list of allowed atomic_numbers
 *
 *         	the number of atoms described by the file must match the specified number of atoms in the first line
 *
 * @param filename (string) the name of the file containing molecule data for the cluster
 * @returns Cluster
 * @author Anna Weber
 * @version 0.0.1
 **/
Cluster readfile(string &filename);
#endif
