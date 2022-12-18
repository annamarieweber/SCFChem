/**
 * @file clusterconfigexception.h
 * @author Anna Weber (anna.weber@scfchem.com)
 * @brief Definition for cluster configuration exceptions
 * @version 0.1
 * @date 2022-12-17
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef CLUSTER_CONFIG_EXCEPTION
#define CLUSTER_CONFIG_EXCEPTION
#include <exception>
#include <string>
#include <iostream>
#include "exceptionconstants.h"
using std::string;

using namespace exception_constants;

class ClusterConfigException : public std::exception
{
private:
  string message;

public:
  ClusterConfigException(Error_Type);
  void displayError();
  string what();
};

#endif
