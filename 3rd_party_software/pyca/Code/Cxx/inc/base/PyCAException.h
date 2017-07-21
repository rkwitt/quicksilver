/* ================================================================
 *
 * PyCA Project
 *
 * Copyright (c) J. Samuel Preston, Linh K. Ha, Sarang C. Joshi. All
 * rights reserved.  See Copyright.txt or for details.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the above copyright notice for more information.
 *
 * ================================================================ */


#ifndef __PYCA_EXCEPTION_H__
#define __PYCA_EXCEPTION_H__

#ifndef SWIG

#include <string>
#include <exception>

#endif // SWIG

#include "StringUtils.h"

namespace PyCA {

/**
 * Class representing an error thrown by the PyCA libraries
 */
class PyCAException : public std::exception
{
public:
  PyCAException(const char *file="unknown file", int line=0,
		   const std::string& text = "undefined exception",
		   const std::exception *cause = NULL) :
    mTypeDescription("PyCAException")
  {
    mText = text;
    if(cause){
      mText = mText + "\n   Caused by:\n" + cause->what();
    }
    mLocation = StringUtils::strPrintf("%s:%d", file, line);
  }

  virtual ~PyCAException() throw() {}

  std::string Text() const { return mText; }
  std::string Location() const { return mLocation; }
  std::string TypeDescription() const { return mTypeDescription; }

  virtual const char* what() const throw (){
    static std::string ex; 
    ex = TypeDescription() + " : From " + Location() + " : " + Text();
    return ex.c_str();
  }
  
protected:

  std::string mText;
  std::string mLocation;
  std::string mTypeDescription;

};

} // end namespace PyCA

#endif // __ATLASWERKS_EXCEPTION_H__
