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

#include <stdio.h>
#include <cstdio>

#include "StringUtils.h"
#include "PyCAException.h"
#include "PyCAThread.h"

namespace PyCA {

std::string
StringUtils::
toUpper(const std::string& s)
{
  std::string t(s);
  for (unsigned int i = 0; i < s.size(); ++i)
  {
    if (t[i] >= 'a' && t[i] <= 'z')
    {
      t[i] -= 'a' - 'A';
    } 
  }
  return t;      
}

std::string
StringUtils::
toLower(const std::string& s)
{
  std::string t(s);
  for (unsigned int i = 0; i < s.size(); ++i)
  {
    if (t[i] >= 'A' && t[i] <= 'Z')
    {
      t[i] += 'a' - 'A';
    } 
  }
  return t;    
}

std::string
StringUtils::
trimWhitespace(const std::string& s)
{
  int length = s.size();
  int spos = 0;
  while (spos < length && isspace(s.at(spos))) spos++;
  int epos = length - 1;
  while (epos > spos && isspace(s.at(epos))) epos--;
  return s.substr(spos, epos - spos + 1);    
}

void
StringUtils::
tokenize(const std::string& s, 
         const std::string& delimiters, 
         std::list<std::string>& tokens)
{
  unsigned int currIndex = 0;
  while (currIndex < s.size())
  {
    int currSize = 0;
    while (currIndex + currSize < s.size()
           && delimiters.find(s[currIndex + currSize]) == std::string::npos)
    {
      currSize++;
    }
    if (currSize > 0)
    {
      tokens.push_back(s.substr(currIndex, currSize));
    }
    currIndex += currSize + 1;
  }  
}

std::string
StringUtils::
strPrintf(const char *format, ...){
  va_list argptr;
  va_start(argptr, format);
  static char buff[STRPRINTF_BUFFSIZE];
  static Mutex mutex = PYCA_MUTEX_INITIALIZER;
  int buffSz = STRPRINTF_BUFFSIZE;
  // attempt to format to buffer
  
  // have to synchronize use of static resource
  mutex_lock(&mutex);
  int nChars = vsnprintf(buff, buffSz, format, argptr);
  
  if(nChars >= 0 && nChars < buffSz){
    va_end(argptr);
    std::string result(buff);
    mutex_unlock(&mutex);
    return result;
  }else{
    mutex_unlock(&mutex);
  }
  
  // if that didn't work we'll try allocating memory on the heap
  char *heapbuff = NULL;
  while(nChars < 0 || nChars >= buffSz){
    buffSz *= 2;
    if(buffSz > STRPRINTF_MAXSIZE){
      if(heapbuff) free(heapbuff);
      throw PyCAException(__FILE__, __LINE__, "strPrintf string exceeded maximum size");
    }
    heapbuff = (char*)realloc(heapbuff, buffSz);
    // reset the argptr
    va_end(argptr);
    va_start(argptr, format);
    nChars = vsnprintf(heapbuff, buffSz, format, argptr);
  }

  std::string result = heapbuff;
  
  free(heapbuff);
  va_end(argptr);

  return result;
}

std::string
StringUtils::
getPathDirectory(const std::string& s)
{
  std::string::size_type pos = s.find_last_of(STRINGUTILS_PATH_SEPERATOR);
  if (pos == std::string::npos)
  {
    return "";
  }
  else
  {
    return s.substr(0, pos + 1);
  }
}

std::string
StringUtils::
getPathFile(const std::string& s)
{
  std::string::size_type pos = s.find_last_of(STRINGUTILS_PATH_SEPERATOR);
  if (pos == std::string::npos)
  {
    return s;
  }
  else
  {
    return s.substr(pos + 1);
  }
}

std::string
StringUtils::
getPathExtension(const std::string& s)
{
  std::string extension = s;
  int pos = extension.find_last_of(".");
  extension.erase(extension.begin(),extension.begin() + pos + 1 );
  return extension;
}

std::string 
StringUtils::
stripPathExtension(const std::string& s)
{
  std::string result = s;
  int pos = result.find_last_of(".");
  result.erase(result.begin() + pos, result.end());
  return result;
}

std::string
StringUtils::
forcePathExtension(const std::string& filename, const std::string& extension)
{
  std::string result = filename;
  if (getPathExtension(result) != extension)
  {
    result += std::string(".") + extension;
  }
  return result;
}

} // end namespace PyCA








