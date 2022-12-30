#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "cluster.h"
#include "constants.h"
#include "exceptionconstants.h"
#include "atomidentifiers.h"
#include "filereader.h"
#include "clusterconfigexception.h"

using std::string;
using std::vector;

using namespace constants::exceptions;
using namespace constants::atomidentifiers;

namespace util
{
	namespace io
	{

		void validateAtomCount(int atom_count, int expected_count, CompareMethod method)
		{
			if (method == LTE && atom_count > expected_count)
			{
				throw ClusterConfigException(ATOM_COUNT_MISMATCH);
			}
			if (method == EQ && atom_count != expected_count)
			{
				throw ClusterConfigException(ATOM_COUNT_MISMATCH);
			}
		}

		int getAtomNum(string atom_id)
		{

			for (int i = 0; i < atom_id.length(); i++)
			{
				if (!isdigit(atom_id[i]))
				{
					try
					{
						return ATOM_NUM.at(atom_id);
					}
					catch (const std::exception &e)
					{
						throw ClusterConfigException(INVALID_ATOM_IDENTIFIER);
					}
				}
			}
			return stoi(atom_id);
		}

		Cluster readfile(string &filename)
		{
			std::ifstream infile(filename);
			if (infile.is_open())
			{
				string atom_id;
				int atom_num;
				double x_coord;
				double y_coord;
				double z_coord;

				// read first line
				int expected_atoms;
				int miscnum;
				infile >> expected_atoms;
				infile >> miscnum;
				int atoms_counted = 0;
				Cluster cluster(expected_atoms);

				while (infile >> atom_id >> x_coord >> y_coord >> z_coord)
				{
					util::io::validateAtomCount(atoms_counted, expected_atoms, util::io::LTE);
					atom_num = getAtomNum(atom_id);
					cluster.addAtom(atoms_counted, atom_num, x_coord, y_coord, z_coord, 0, 0);
					atoms_counted++;
				}

				util::io::validateAtomCount(atoms_counted, expected_atoms, util::io::EQ);
				infile.close();
				return cluster;
			}
			else
			{
				throw std::invalid_argument("Can't open file to read.");
			}
		}
	}
}
