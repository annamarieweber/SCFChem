#include "clusterconfigexception.h"
#include "exceptionconstants.h"
#include <string>
#include <iostream>

using std::string;

using namespace constants::exceptions;

ClusterConfigException::ClusterConfigException(std::string t)
{
  message = t;
}

void ClusterConfigException::displayError()
{
  std::cout << message << std::endl;
}

string ClusterConfigException::what()
{
  return message;
}