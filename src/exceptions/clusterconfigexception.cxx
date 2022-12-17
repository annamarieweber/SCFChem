#include "clusterconfigexception.h"
#include "exceptionconstants.h"
#include <string>
#include <iostream>

using std::string;

using namespace exception_constants;

ClusterConfigException::ClusterConfigException(Error_Type t)
{
  std::cout << "making exception " << ERROR_MESSAGE.at(t) << std::endl;
  message = ERROR_MESSAGE.at(t);
}

void ClusterConfigException::displayError()
{
  std::cout << message << std::endl;
}

string ClusterConfigException::what()
{
  return message;
}