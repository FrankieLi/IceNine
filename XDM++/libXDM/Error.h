/////////////////////////////////////////////////////////////////
//
//  File:    Error.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose: Runtime error and assert routines
//
//
/////////////////////////////////////////////////////////////////

#ifndef _ERROR_H_
#define _ERROR_H_

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstring>
////////////////////////////////////////////////////
//
//  Debug Inline Functions
//
//
////////////////////////////////////////////////////
inline void RUNTIME_ASSERT(bool pred, const string msg)
{
  if(!pred){
    std::cerr << "RUNTIME_ASSERT\n";
    std::cerr << msg;
    exit(0);
  }
}

////////////////////////////////////////////////////
inline void RUNTIME_WARNING(bool pred, const string msg)
{
  if(!pred){
    std::cerr << "RUNTIME_WARNING\n";
    std::cerr << msg;
  }
}

#endif
