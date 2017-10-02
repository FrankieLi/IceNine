///////////////////////////////////////////////////////
//
//  Parser.cpp
//
//  Purpose:  Implements generalized parsing tools, such as tokenizer
//            parsers, and lexer
//       
///////////////////////////////////////////////////////

#include "Parser.h"

namespace GeneralLib
{

  ///----------------------------------------------------------------------------
  ///  Tokenize - given the character string buffer sBuf and delimiters,
  ///  a set of tokens will be outputted as vectors of vectors of strings (lines of strings)
  ///
  ///----------------------------------------------------------------------------
  void Tokenize( vector<vector<string> > &vsTokens, const string &sBuf,  const string &sDelimiters)
  {
    
    Size_Type iCurrentPos = 0;
    Size_Type iDelimiterPos = 0;
    Size_Type iEndLine = 0;


    while(iCurrentPos < sBuf.size()){
      vector<string> vsCurrentLine;
		
		
      iEndLine = sBuf.find_first_of("\n", iEndLine + 1);
		
      if(iEndLine == string::npos)
        return;
		
      while(iCurrentPos < iEndLine){
        string currentToken;
	
        // find either space of tab or endline
        iDelimiterPos = sBuf.find_first_of(sDelimiters.c_str(), iDelimiterPos);
        iCurrentPos = sBuf.find_first_not_of(sDelimiters.c_str(), iCurrentPos);
			
        if( (iCurrentPos == string::npos) || (iDelimiterPos == string::npos) )
          break;

        for(Size_Type i = iCurrentPos; i < iDelimiterPos; i++)
          currentToken.append(sizeof(char), sBuf.at(i));
			

        if (currentToken.size() > 0)
          vsCurrentLine.push_back(currentToken);
			
        iCurrentPos = iDelimiterPos + 1;
        iDelimiterPos ++;
      }
      //      if(vsCurrentLine.size() > 0){
        vsTokens.push_back(vsCurrentLine);
        /*     }else{
        vsCurrentLine.push_back("#");
        vsCurrentLine.push_back("EMPTY LINE");
        vsTokens.push_back(vsCurrentLine);
        }*/
    }
  }




}
