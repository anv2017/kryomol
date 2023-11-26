/*****************************************************************************************
                            exception.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#ifndef QRYOMOLEXCEPTION_H
#define QRYOMOLEXCEPTION_H

#include <exception>
#include <string>


namespace kryomol
{
/**
  A simple exception class*/
class Exception : public std::exception
{
  public:
  Exception(const std::string& e) throw() : _error(e) {}
  virtual ~Exception() throw() {}
  const char* what() const throw() { return _error.c_str(); }
  private:
  std::string _error;
  
  
};
}
#endif

