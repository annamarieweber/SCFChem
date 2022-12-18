/**
 * @file baseexception.h
 * @author Anna Weber (anna@scfchem.com)
 * @brief Definition for base Exception
 * @version 0.1
 * @date 2022-12-17
 *
 * @copyright Copyright (c) 2022
 *
 */
#ifndef BASE_EXCEPTION
#define BASE_EXCEPTION
#include <exception>
#include <string>
#include <iostream>

using std::string;

class BaseException : public std::exception
{
private:
  string message;

public:
  BaseException(string msg);
  void displayError();
};

#endif
