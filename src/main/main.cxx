#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <string>
#include "cluster.h"
#include "filereader.h"
#include "baseexception.h"
#include "clusterconfigexception.h"
#include "algoutils.h"
#include "cndo2.h"
#include "scf.h"
#include <exception> // std::set_terminate
using std::string;
using std::vector;

int main(int argc, char *argv[])
{
	std::set_terminate([]() -> void
					   {
        std::cerr << "terminate called after throwing an instance of ";
        try
        {
            std::rethrow_exception(std::current_exception());
        }
        catch (const std::exception &ex)
        {
            std::cerr << typeid(ex).name() << std::endl;
            std::cerr << ex.what() << std::endl;
        }
        catch (...)
        {
            std::cerr << typeid(std::current_exception()).name() << std::endl;
            std::cerr << "an unknown error occured" << std::endl;
        }
        std::cerr << "errno: " << errno << ": " << std::strerror(errno) << std::endl;
        std::abort(); });

	if (argc < 2)
	{
		std::cout << "Please input the file name!" << std::endl;
	}

	try
	{
		string filename = argv[1];
		Cluster cluster = util::io::readfile(filename);
		std::cout << cluster << std::endl;

		std::cout << "numBasisFunctions: " << cluster.numBasisFunctions << std::endl;
		std::cout << "numElectronPairs: " << cluster.numElectronPairs << std::endl;
		std::cout << "BasisFunctions:" << std::endl
				  << cluster.basisFunctions() << std::endl;
		// std::cout << "OV_mat_Ra: " << std::endl
		// 		  << algo::utils::overlapMatrixRA(cluster) << std::endl;
		// std::cout << "H_mat: " << std::endl
		// 		  << cluster.extendedHuckelHamiltonian() << std::endl;
		// std::cout << "MO Coeffs:" << std::endl
		// 		  << cluster.molecularOrbitalCoefficients() << std::endl;
		// std::cout << "The molecule in file " << filename << " has energy: " << cluster.molecularEnergy() << std::endl;
		// std::cout << "sBasisFunctions:" << std::endl
		// 		  << cluster.sBasisFunctions() << std::endl;
		// std::cout << "Gamma RA: " << std::endl
		// 		  << algo::utils::gammaMatrixRA(cluster) << std::endl;
		std::cout << "OV_mat: " << std::endl
				  << algo::utils::overlapMatrix(cluster) << std::endl;
		// std::cout << "Gamma RA: " << std::endl
		// 		  << algo::utils::gammaMatrixRA(cluster) << std::endl;
		// std::cout << "FockMatrix: " << std::endl
		// 		  << cluster.cndo2FockMatrix(mat(testpa), mat(testpb)) << std::endl;

		algo::scf::calcSCFEnergy(10e-4, cluster);
	}
	catch (std::invalid_argument &e)
	{
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch (BaseException &e)
	{
		e.displayError();
		return EXIT_FAILURE;
	}
	catch (ClusterConfigException &e)
	{
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	// catch (const std::exception &ex)
	// {
	// 	std::cout << ex.what() << std::endl;
	// 	return EXIT_FAILURE;
	// }
	// catch (...)
	// {
	// 	std::cout << "An unknown Error occurred" << std::endl;
	// 	return EXIT_FAILURE;
	// }

	return EXIT_SUCCESS;
}
