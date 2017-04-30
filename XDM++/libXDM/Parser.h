///////////////////////////////////////////////////////
//
//  Parser.h
//
//  Purpose:  Conatins generalized parsing tools, such as tokenizer
//            parsers, and lexer
//       
///////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "Debug.h"
#include "Types.h"
#include "Error.h"

#ifndef PARSER_H_
#define PARSER_H_

using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;


namespace GeneralLib
{
  void Tokenize( vector<vector<string> > &vsTokens, const string &sBuf,
                 const string &sDelimiters);

}

#endif
