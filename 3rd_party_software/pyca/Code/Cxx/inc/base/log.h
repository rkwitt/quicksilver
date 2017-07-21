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

/**
 *  Note: most of this code was taken from Dr. Dobb's Journal October
 *  2007 issue, node and thread info added
 */

#ifndef __LOG_H__
#define __LOG_H__

#include <sstream>
#include <string>
#include <cstdio>
#include <map>
#include "StringUtils.h"
#include "PyCAException.h"
#include "PyCAThread.h"

#ifdef MPI_ENABLED
#include <mpi.h>
#endif // MPI_ENABLED

namespace PyCA {

inline std::string NowTime();

enum TLogLevel {logERROR, logWARNING, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4};

template <typename T>
class Log
{
  typedef std::map<ThreadID, std::string> ThreadNameMap;
public:
  Log(const char *file="unknown file", int line=0);
  virtual ~Log();
  std::ostringstream& Get(TLogLevel level = logINFO, bool printNode = false, bool printThread = false);
public:
  static TLogLevel& ReportingLevel();
  static std::string ToString(TLogLevel level);
  static TLogLevel FromString(const std::string& level);
  static void SetNodeID(int nodeId);
  static void SetNodeName(const std::string &name);
  static void SetThreadID(int threadId);
  static void SetThreadName(const std::string &name);
  static std::string ThreadName();
  static bool IsSingleNode();
protected:
  std::ostringstream mOs;
  TLogLevel mLevel;
  const char *mFile;
  int mLine;
protected:
  static std::string sNodeName;
  static ThreadNameMap sThreadNameMap;
  static unsigned int sUnnamedThreadCount;
  static bool sIsSingleNode;
  static bool sIsSingleNodeTested;
private:
  // purposefully unimplemented
  Log(const Log&);
  Log& operator=(const Log&);
};

template <typename T>
std::string Log<T>::sNodeName = "UnnamedNode";

template <typename T> 
typename Log<T>::ThreadNameMap Log<T>::sThreadNameMap;

template <typename T>
unsigned int Log<T>::sUnnamedThreadCount = 0;

template <typename T>
void 
Log<T>::
SetNodeID(int nodeId)
{
  std::string name  = StringUtils::strPrintf("Node%02d",nodeId);
  SetNodeName(name);
}

template <typename T>
void 
Log<T>::
SetNodeName(const std::string &name)
{
  sNodeName = name;
}

template <typename T>
void 
Log<T>::
SetThreadID(int threadId)
{
  std::string name  = StringUtils::strPrintf("Thread%02d",threadId);
  SetThreadName(name);
}

template <typename T>
void 
Log<T>::
SetThreadName(const std::string &name)
{
  ThreadID tid = current_thread_id();
  ThreadNameMap::iterator it = sThreadNameMap.begin();
  for(;it != sThreadNameMap.end(); ++it){
    if(compare_threads(it->first, tid)){
      it->second = name;
      break;
    }
  }
  // if tid isn't in map yet
  if(it == sThreadNameMap.end()){
    sThreadNameMap[tid] = name;
  }
}

template <typename T>
std::string
Log<T>::
ThreadName()
{
  ThreadID tid = current_thread_id();
  ThreadNameMap::iterator it = sThreadNameMap.begin();
  std::string name;
  for(;it != sThreadNameMap.end(); ++it){
    if(compare_threads(it->first, tid)){
      name = it->second;
      break;
    }
  }
  // if tid isn't in map yet
  if(it == sThreadNameMap.end()){
    name = StringUtils::strPrintf("UnnamedThread%02d",sUnnamedThreadCount++);
    sThreadNameMap[tid] = name;
  }
  return name;
}

template <typename T>
bool
Log<T>::
IsSingleNode(){
  if(sIsSingleNodeTested){
    return sIsSingleNode;
  }else{
#ifdef MPI_ENABLED
    int nodeId=0;
    MPI_Comm_rank( MPI_COMM_WORLD, &nodeId);
    sIsSingleNode = (nodeId == 0);
#else
    sIsSingleNode = true;
#endif
    sIsSingleNodeTested = true;
  }
  return sIsSingleNode;
}

template <typename T>
Log<T>::Log(const char *file, int line)
  :  mFile(file),
     mLine(line)
{
}

template <typename T>
std::ostringstream& Log<T>::Get(TLogLevel level, bool printNode, bool printThread)
{
  mLevel = level;
  mOs << "- " << NowTime();
  mOs << " " << ToString(mLevel);
  if(printNode){
    mOs << " " << sNodeName;
  }
  if(printThread){
    mOs << " " << ThreadName();
  }
  mOs << ": ";
  // tab indention for debug messages
  mOs << std::string(level > logDEBUG ? level - logDEBUG : 0, '\t');
  return mOs;
}

template <typename T>
Log<T>::~Log()
{
  mOs << std::endl;
  T::Output(mOs.str());
  if(mLevel == logERROR){
    throw PyCAException(mFile, mLine, mOs.str());
  }
}

template <typename T>
TLogLevel& Log<T>::ReportingLevel()
{
  static TLogLevel reportingLevel = logDEBUG4;
  return reportingLevel;
}

template <typename T>
std::string Log<T>::ToString(TLogLevel level)
{
  static const char* const buffer[] = {"ERROR", "WARNING", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"};
  return buffer[level];
}

template <typename T>
TLogLevel Log<T>::FromString(const std::string& level)
{
  if (level == "DEBUG4")
    return logDEBUG4;
  if (level == "DEBUG3")
    return logDEBUG3;
  if (level == "DEBUG2")
    return logDEBUG2;
  if (level == "DEBUG1")
    return logDEBUG1;
  if (level == "DEBUG")
    return logDEBUG;
  if (level == "INFO")
    return logINFO;
  if (level == "WARNING")
    return logWARNING;
  if (level == "ERROR")
    return logERROR;
  Log<T>().Get(logWARNING) << "Unknown logging level '" << level << "'. Using INFO level as default.";
  return logINFO;
}

class Output2STDERR
{
public:
  static void Output(const std::string& msg);
};

inline void Output2STDERR::Output(const std::string& msg)
{   
  fprintf(stderr, "%s", msg.c_str());
  fflush(stderr);
}

class ErrLog : public Log<Output2STDERR> {
public:
  ErrLog(const char *file, int line)
    : Log<Output2STDERR>(file, line)
  {
  }
};
//typedef Log<Output2FILE> FILELog;

#ifndef ERRLOG_MAX_LEVEL
#define ERRLOG_MAX_LEVEL logDEBUG4
#endif

#define LOG(level)							\
  if (level > ERRLOG_MAX_LEVEL) ;					\
  else if (level > ErrLog::ReportingLevel()) ;				\
  else ErrLog(__FILE__,__LINE__).Get(level)

#define LOGSINGLENODE(level)						\
  if (level > ERRLOG_MAX_LEVEL) ;					\
  else if (level > ErrLog::ReportingLevel()) ;				\
  else if (IsSingleNode())                                              \
    ErrLog(__FILE__,__LINE__).Get(level, true)

#define LOGNODE(level)				        		\
  if (level > ERRLOG_MAX_LEVEL) ;					\
  else if (level > ErrLog::ReportingLevel()) ;				\
  else ErrLog(__FILE__,__LINE__).Get(level, true)

#define LOGTHREAD(level)			        		\
  if (level > ERRLOG_MAX_LEVEL) ;					\
  else if (level > ErrLog::ReportingLevel()) ;				\
  else ErrLog(__FILE__,__LINE__).Get(level, false, true)

#define LOGNODETHREAD(level)			        		\
  if (level > ERRLOG_MAX_LEVEL) ;					\
  else if (level > ErrLog::ReportingLevel()) ;				\
  else ErrLog(__FILE__,__LINE__).Get(level, true, true)

class Output2FILE
{
public:
  static FILE*& Stream();
  static void Output(const std::string& msg);
};

inline FILE*& Output2FILE::Stream()
{
  static FILE* pStream = stderr;
  return pStream;
}

inline void Output2FILE::Output(const std::string& msg)
{   
  FILE* pStream = Stream();
  if (!pStream)
    return;
  fprintf(pStream, "%s", msg.c_str());
  fflush(pStream);
}

class FILELog : public Log<Output2FILE> {};
//typedef Log<Output2FILE> FILELog;

#ifndef FILELOG_MAX_LEVEL
#define FILELOG_MAX_LEVEL logDEBUG4
#endif

#define FILE_LOG(level)							\
  if (level > FILELOG_MAX_LEVEL) ;					\
  else if (level > FILELog::ReportingLevel() || !Output2FILE::Stream()) ; \
  else FILELog().Get(level)

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)

#include <windows.h>

inline std::string NowTime()
{
  const int MAX_LEN = 200;
  char buffer[MAX_LEN];
  if (GetTimeFormatA(LOCALE_USER_DEFAULT, 0, 0, 
		     "HH':'mm':'ss", buffer, MAX_LEN) == 0)
    return "Error in NowTime()";

  char result[100] = {0};
  static DWORD first = GetTickCount();
  std::sprintf(result, "%s.%03ld", buffer, (long)(GetTickCount() - first) % 1000); 
  return result;
}

#else

#include <sys/time.h>

inline std::string NowTime()
{
  char buffer[11];
  time_t t;
  time(&t);
  tm r = {0};
  strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
  struct timeval tv;
  gettimeofday(&tv, 0);
  char result[100] = {0};
  std::sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);
  std::string rtnString(result);
  return rtnString;
}

#endif //WIN32

} // end namespace PyCA

#endif //__LOG_H__
