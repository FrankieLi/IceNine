/////////////////////////////////////////////////////////////////
//
//  File:    Error.h
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//  Purpose:   Define all debug functions/macros here so they can be removed
//             with the appropiate compilation flags.  Each functional unit has
//             its own debugging flag so that specific debugging functions can
//             be removed at compile time.
//
/////////////////////////////////////////////////////////////////

#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstring>


#if defined(DEBUG_LEVEL_MAX)

   #define DEBUG_LEVEL_1
   #define DEBUG_LEVEL_2
   #define DEBUG_LEVEL_3
#endif

#ifndef UNUSED
#define UNUSED( x ) ( NULL )
#endif

using std::string;
using std::cerr;


////////////////////////////////////////////
//
//  Inline function prototypes
//
/////////////////////////////////////////////
//void DEBUG_ASSERT(bool pred, const string &msg);


////////////////////////////////////////////////////
//
//  Debug Inline Functions
//
//
////////////////////////////////////////////////////
#ifdef DEBUG_LEVEL_1	// DEBUG mode
//inline void DEBUG_ASSERT(bool pred, const string &msg = "")
inline void DEBUG_ASSERT(bool pred, const char* msg = NULL)
{
	if(!pred){
		cerr << "DEBUG_ASSERT\n";
		cerr << msg;
		exit(0);
	}
}
#else
//inline void DEBUG_ASSERT(bool pred, const string &msg = "")
inline void DEBUG_ASSERT(bool pred, const char * msg = NULL)
{
  UNUSED( bPred );  UNUSED( sMsg );
}

#endif

#ifdef DEBUG_LEVEL_1	// DEBUG mode
inline void DEBUG_ALERT(bool pred, const string &msg = " ")
//inline void DEBUG_ALERT(bool pred, const char * msg = NULL)
{
	if(!pred){
		cerr << msg;
	}
}
#else
//inline void DEBUG_ALERT(bool pred, const char * msg = NULL)
inline void DEBUG_ALERT(bool bPred, const string &sMsg = " ")
{
  UNUSED( bPred ); UNUSED( sMsg );
}
#endif


// ********************************************  //
//
//  Master Debug Flag.  Disable this to remove all 
//  debug statements
//
// ********************************************  // 
#if defined(DEBUG_LEVEL_1) || defined(DEBUG_LEVEL_MAX)

// Define debug functions here //
#define DEBUG_CALL( fn ) (fn)

#else

// Define disabled debug functions //
#define DEBUG_CALL( fn ) (NULL)
#endif 

// ************* end DEBUG_LEVEL_1 ************** //





// ********************************************  //
//
//  Debug functions for ConfigFile.cpp
//
// ********************************************  // 

#if defined(DEBUG_CONFIG_FILE) ||defined(DEBUG_LEVEL_MAX)

#define CONFIG_DEBUG( fn ) ( DEBUG_CALL(fn) )

#else
#define CONFIG_DEBUG( fn ) ( NULL )

#endif
// ************* end DEBUG_CONFIG_FILE ************** //




// ********************************************  //
//
//  Debug functions for MicIO.cpp
//
// ********************************************  // 

#if defined(DEBUG_MIC_IO) || defined(DEBUG_LEVEL_MAX)

#define MIC_DEBUG( fn ) ( DEBUG_CALL(fn) )

#else
#define MIC_DEBUG( fn ) ( NULL )

#endif
// ************* end DEBUG_SAMPLE ************** //


// ********************************************  //
//
//  Debug functions for Sample.cpp
//
// ********************************************  // 

#if defined(DEBUG_SAMPLE)  || defined(DEBUG_LEVEL_MAX)

#define SAMPLE_DEBUG( fn ) ( DEBUG_CALL(fn) )

#else
#define SAMPLE_DEBUG( fn ) ( NULL )

#endif
// ************* end DEBUG_SAMPLE************** //



// ********************************************  //
//
//  Debug functions for QuadTree.h
//
// ********************************************  // 

#if defined(DEBUG_QUADTREE)  || defined(DEBUG_LEVEL_MAX)

#define QUADTREE_DEBUG( fn ) ( DEBUG_CALL(fn) )

#else
#define QUADTREE_DEBUG( fn ) ( NULL )

#endif
// ************* end DEBUG_SAMPLE************** //


////////////////////////////////////////////////////
//
//  Debug Inline Functions
//
//
////////////////////////////////////////////////////




#endif 
//**************  end __DEBUG_H__ *************** //

