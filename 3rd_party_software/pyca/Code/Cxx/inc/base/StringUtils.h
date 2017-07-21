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

#ifndef __STRING_UTILS_H__
#define __STRING_UTILS_H__

#ifndef SWIG

#include <string>
#include <list>
#include <cstdlib>
#include <sstream>
#include <stdarg.h>

#endif // SWIG

#ifdef WIN32
#define STRINGUTILS_PATH_SEPERATOR "\\"
#else
#define STRINGUTILS_PATH_SEPERATOR "/"
#endif

#define STRPRINTF_BUFFSIZE 1024
#define STRPRINTF_MAXSIZE 1000000

namespace PyCA {

class StringUtils
{
public:
  // make all alpha chars upper case
  static std::string toUpper(const std::string& s);

  // make all alpha chars lower case
  static std::string toLower(const std::string& s);

  // remove whitespace from either end of string
  static std::string trimWhitespace(const std::string& s);

  // fill a list with the tokens in s.  delimiters are thrown out.
  static void tokenize(const std::string& s, 
                       const std::string& delimiters, 
                       std::list<std::string>& tokens);

  // convert to bool
  static bool toBool(const std::string& s)
  {
    std::string trimmed = trimWhitespace(s);
    if(toUpper(trimmed) == "TRUE"){
      return true;
    }else if(toUpper(trimmed) == "FALSE"){
      return false;
    }else if(trimmed == "1"){
      return true;
    }else if(trimmed == "0"){
      return false;
    }else if(trimmed == "T"){
      return true;
    }else if(trimmed == "F"){
      return false;
    }

    return (toDouble(trimmed) != 0.0);
  }

  // convert to integer
  static
  int
  toInt(const std::string& s)
  {
    return atoi(s.c_str());
  }

  // convert to double
  static
  float
  toFloat(const std::string& s)
  {
    return (float)atof(s.c_str());
  }

  // convert to double
  static
  double
  toDouble(const std::string& s)
  {
    return atof(s.c_str());
  }

  // convert int to string
  static
  std::string
  toString(int i)
  {
    std::ostringstream sstr;
    sstr << i;
    return sstr.str();
  }

  // convert double to string
  static
  std::string
  toString(double i)
  {
    std::ostringstream sstr;
    sstr << i;
    return sstr.str();
  }

  // convert bool to string
  static
  std::string
  toString(bool i)
  {
    std::string s;
    if(i){
      s = "True";
    }else{
      s = "False";
    }
    return s;
  }

  static 
  std::string
  strPrintf(const char*format, ...);

  // return the path-to-directory part of a path
  // e.g., /afs/radonc/pkg/test.txt yields /afs/radonc/pkg/
  static std::string getPathDirectory(const std::string& s);

  // return the filename part of a path
  // e.g.,  /afs/radonc/pkg/test.txt yields test.txt
  static std::string getPathFile(const std::string& s);

  // return the filename part of a path
  // e.g.,  /afs/radonc/pkg/test.txt yields txt
  static std::string getPathExtension(const std::string& s);

  // remove extension from path
  // e.g.,  /afs/radonc/pkg/test.txt yields /afs/radonc/pkg/test
  static std::string stripPathExtension(const std::string& s);

  // require a particular extension 
  // forceExtension("foo.png", "png") yields "foo.png".
  // forceExtension("foo.eps", "png") yields "foo.eps.png".
  static std::string forcePathExtension(const std::string& filename, 
                                        const std::string& extension);
};

} // end namespace PyCA

#endif // __STRING_UTILS_H__
