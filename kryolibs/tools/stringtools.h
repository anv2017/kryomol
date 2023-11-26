/*****************************************************************************************
                            stringtools.h  -  description
                             -------------------
This file is part of the KryoMol project.
For more information, see <http://kryomol.sourceforge.io/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
******************************************************************************************/

#include <cstdlib>
#include <string>
#include <vector>
#include <locale>
#include <iostream>
#include <cmath>

#ifndef STRINGTOOLS_H
#define STRINGTOOLS_H

namespace kryomol
{
  /** @brief a STL based tokenizer 
    build a vector of tokens separated by the characters specified in delimiters*/
  class  StringTokenizer : public std::vector<std::string>
  {
  public:
    StringTokenizer(const std::string &rStr, const std::string &rDelimiters = " ,\n");

  };

  inline StringTokenizer::StringTokenizer(const std::string &rStr, const std::string &rDelimiters)
  {
    std::string::size_type lastPos(rStr.find_first_not_of(rDelimiters, 0));
    std::string::size_type pos(rStr.find_first_of(rDelimiters, lastPos));
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
      push_back(rStr.substr(lastPos, pos - lastPos));
      lastPos = rStr.find_first_not_of(rDelimiters, pos);
      pos = rStr.find_first_of(rDelimiters, lastPos);
    }
  }
  
  class TokByString : public std::vector <std::string>
  {
  public:
    /** build a vector of tokens divided exactly by the string specified in rDelimiters*/
    TokByString(const std::string& rStr, const std::string& rDelimiters = "\n");

  };

  inline TokByString::TokByString(const std::string& rStr, const std::string& rDelimiters)
  {
    std::string::size_type last(rStr.find(rDelimiters,0));

    std::string::size_type first=0;
    while ( last != std::string::npos  )
    {
      push_back(rStr.substr(first, last-first));
      first=last+rDelimiters.size();
      last =rStr.find(rDelimiters,first);
    }


  }
  /** transform each character into string s into its uppercase equivalent*/
  inline std::string toupper(const std::string& s)
  {
    std::string cs=s;
    std::string::iterator it;
    for(it=cs.begin();it!=cs.end();it++)
    {
      (*it)=std::toupper(*it,std::cout.getloc());
    }

    return cs;
  }

  /** transform each character into string s into its lowercase equivalent*/
  inline std::string tolower(const std::string& s)
  {
    std::string ls=s;
    std::string::iterator it;
    for(it=ls.begin();it!=ls.end();++it)
    {
      (*it)=std::tolower(*it,std::cout.getloc());
    }

    return ls;
  }
  /** an specialized atof with support for scientific notation*/
  inline double atof(std::string& s)
  {
    StringTokenizer tok(s,"DE");
    double mantisa=std::atof(tok[0].c_str());
    if(tok.size() > 1)
      mantisa*=pow((double)10.0,std::atoi(tok[1].c_str()));
    return mantisa;
  }
  
  /** remove character ch from string s*/
  inline void remove_character(char ch, std::string& s)
  {
    size_t pos = 0;

    while (true)
    {

      pos = s.find(ch, pos);

      // check if there are no more ch in str

      if (pos == std::string::npos)
        break;

      s.replace(pos, 1, "");

    }
  }

  /** remove all numbers in string s*/
  inline std::string& remove_numbers(std::string& s)
  {
    char numbers[10] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };
    for(int i=0;i<10;i++)
      remove_character(numbers[i],s);

    return s;
  }

  /** true if string snumber represents an integer number*/
 inline bool isinteger(const std::string& snumber)
  {
#ifdef __GNUC__
//#warning implement "9-" case
#endif
    return ( snumber.find_first_not_of("0123456789-") == std::string::npos
             && snumber.find_first_of("0123456789") != std::string::npos );
  }

  /** true if snumber represents a integer of real number, supports scientific notation*/
   inline bool isnum(const std::string& snumber)
  {
    return (snumber.find_first_not_of("0123456789.,-+eEdD") == std::string::npos );


  }


}

#endif
