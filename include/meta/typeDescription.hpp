/** \file typeDescription.hpp
  * \brief Type description helper classes
  * \ingroup utils_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of mxlib.
//
// mxlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mxlib is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with mxlib.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//

#ifndef typeDescription_hpp
#define typeDescription_hpp

#include <string>
#include <vector>

namespace mx
{
namespace meta 
{
   
/// Struct which contains static members describing a type.
/** Specializations are provided for the fundamental types, std::string, and std::vector of fundamental types.
  */
template<typename T>
struct typeDescription
{
   typedef T type; ///< The type iteself.
   
   /// Returns a unique numeric code for this type.
   static  constexpr int code() { return  -1;}
   
   /// Returns the name of this type.
   static constexpr const char * name(){ return "unknown";}
};

template<>
struct typeDescription<bool>
{
   typedef bool type;
   
   static  constexpr int code() { return  0;}
   
   static constexpr const char * name(){ return "bool";}
};

template<>
struct typeDescription<signed char>
{
   typedef signed char type;
   
   static  constexpr int code() { return  1;}
   
   static constexpr const char * name(){ return "signed char";}
};

template<>
struct typeDescription<unsigned char>
{
   typedef unsigned char type;
   
   static  constexpr int code() { return  2;}
   
   static constexpr const char * name(){ return "unsigned char";}
};

template<>
struct typeDescription<char>
{
   typedef char type;
   
   static  constexpr int code() { return  3;}
   
   static constexpr const char * name(){ return "char";}
};

template<>
struct typeDescription<wchar_t>
{
   typedef wchar_t type;
   
   static  constexpr int code() { return  4;}
   
   static constexpr const char * name(){ return "wchar_t";}
};

template<>
struct typeDescription<char16_t>
{
   typedef char16_t type;
   
   static  constexpr int code() { return 5;}
   
   static constexpr const char * name(){ return "char16_t";}
};

template<>
struct typeDescription<char32_t>
{
   typedef char32_t type;
   
   static  constexpr int code() { return 6;}
   
   static constexpr const char * name(){ return "char32_t";}
};

template<>
struct typeDescription<int>
{
   typedef int type;
   
   static constexpr int code() { return  7;}
   
   static constexpr const char * name(){ return "int";}
};

template<>
struct typeDescription<unsigned int>
{
   typedef unsigned int type;
   
   static  constexpr int code() { return  8;}
   
   static constexpr const char * name(){ return "unsigned int";}
};

template<>
struct typeDescription<short int>
{
   typedef short int type;
   
   static  constexpr int code() { return  9;}
   
   static constexpr const char * name(){ return "short int";}
};

template<>
struct typeDescription<short unsigned int>
{
   typedef short unsigned int type;
   
   static  constexpr int code() { return  10;}
   
   static constexpr const char * name(){ return "short unsigned int";}
};

template<>
struct typeDescription<long int>
{
   typedef long int type;
   
   static  constexpr int code() { return  11;}
   
   static constexpr const char * name(){ return "long int";}
};

template<>
struct typeDescription<long unsigned int>
{
   typedef long unsigned int type;
   
   static  constexpr int code() { return  12;}
   
   static constexpr const char * name(){ return "long unsigned int";}
};

template<>
struct typeDescription<long long int>
{
   typedef long long int type;
   
   static  constexpr int code() { return  13;}
   
   static constexpr const char * name(){ return "long long int";}
};

template<>
struct typeDescription<long long unsigned int>
{
   typedef long long unsigned int type;
   
   static  constexpr int code() { return  14;}
   
   static constexpr const char * name(){ return "long long unsigned int";}
};

template<>
struct typeDescription<float>
{
   typedef float type;
   
   static  constexpr int code() { return  15;}
   
   static constexpr const char * name(){ return "float";}
};

template<>
struct typeDescription<double>
{
   typedef double type;
   
   static  constexpr int code() { return  16;}
   
   static constexpr const char * name(){ return "double";}
};

template<>
struct typeDescription<long double>
{
   typedef long double type;
   
   static  constexpr int code() { return  17;}
   
   static constexpr const char * name(){ return "long double";}
};

template<>
struct typeDescription<std::string>
{
   typedef std::string type;
   
   static  constexpr int code() { return  100;}
   
   static constexpr const char * name(){ return "std::string";}
};


template<>
struct typeDescription<std::vector<bool>>
{
   typedef std::vector<bool> type;
   
   static  constexpr int code() { return  1000;}
   
   static constexpr const char * name(){ return "std::vector<bool>";}
};

template<>
struct typeDescription<std::vector<signed char>>
{
   typedef std::vector<signed char> type;
   
   static  constexpr int code() { return  1001;}
   
   static constexpr const char * name(){ return "std::vector<signed char>";}
};

template<>
struct typeDescription<std::vector<unsigned char>>
{
   typedef std::vector<unsigned char> type;
   
   static  constexpr int code() { return  1002;}
   
   static constexpr const char * name(){ return "std::vector<unsigned char>";}
};

template<>
struct typeDescription<std::vector<char>>
{
   typedef std::vector<char> type;
   
   static  constexpr int code() { return  1003;}
   
   static constexpr const char * name(){ return "std::vector<char>";}
};

template<>
struct typeDescription<std::vector<wchar_t>>
{
   typedef std::vector<wchar_t> type;
   
   static  constexpr int code() { return  1004;}
   
   static constexpr const char * name(){ return "std::vector<wchar_t>";}
};

template<>
struct typeDescription<std::vector<char16_t>>
{
   typedef std::vector<char16_t> type;
   
   static  constexpr int code() { return  1005;}
   
   static constexpr const char * name(){ return "std::vector<char16_t>";}
};

template<>
struct typeDescription<std::vector<char32_t>>
{
   typedef std::vector<char32_t> type;
   
   static  constexpr int code() { return  1006;}
   
   static constexpr const char * name(){ return "std::vector<char32_t>";}
};

template<>
struct typeDescription<std::vector<int>>
{
   typedef std::vector<int> type;
   
   static  constexpr int code() { return  1007;}
   
   static constexpr const char * name(){ return "std::vector<int>";}
};

template<>
struct typeDescription<std::vector<unsigned int>>
{
   typedef std::vector<unsigned int> type;
   
   static  constexpr int code() { return  1008;}
   
   static constexpr const char * name(){ return "std::vector<unsigned int>";}
};

template<>
struct typeDescription<std::vector<short int>>
{
   typedef std::vector<short int> type;
   
   static  constexpr int code() { return  1009;}
   
   static constexpr const char * name(){ return "std::vector<short int>";}
};

template<>
struct typeDescription<std::vector<short unsigned int>>
{
   typedef std::vector<short unsigned int> type;
   
   static  constexpr int code() { return  1010;}
   
   static constexpr const char * name(){ return "std::vector<short unsigned int>";}
};

template<>
struct typeDescription<std::vector<long int>>
{
   typedef std::vector<long int> type;
   
   static  constexpr int code() { return  1011;}
   
   static constexpr const char * name(){ return "std::vector<long int>";}
};

template<>
struct typeDescription<std::vector<long unsigned int>>
{
   typedef std::vector<long unsigned int> type;
   
   static  constexpr int code() { return  1012;}
   
   static constexpr const char * name(){ return "std::vector<long unsigned int>";}
};

template<>
struct typeDescription<std::vector<long long int>>
{
   typedef std::vector<long long int> type;
   
   static  constexpr int code() { return  1013;}
   
   static constexpr const char * name(){ return "std::vector<long long int>";}
};

template<>
struct typeDescription<std::vector<long long unsigned int>>
{
   typedef std::vector<long long unsigned int> type;
   
   static  constexpr int code() { return  1014;}
   
   static constexpr const char * name(){ return "std::vector<long long unsigned int>";}
};

template<>
struct typeDescription<std::vector<float>>
{
   typedef std::vector<float> type;
   
   static  constexpr int code() { return  1015;}
   
   static constexpr const char * name(){ return "std::vector<float>";}
};

template<>
struct typeDescription<std::vector<double>>
{
   typedef std::vector<double> type;
   
   static  constexpr int code() { return  1016;}
   
   static constexpr const char * name(){ return "std::vector<double>";}
};

template<>
struct typeDescription<std::vector<long double>>
{
   typedef std::vector<long double> type;
   
   static  constexpr int code() { return  1017;}
   
   static constexpr const char * name(){ return "std::vector<long double>";}
};

template<>
struct typeDescription<std::vector<std::string>>
{
   typedef std::vector<std::string> type;
   
   static  constexpr int code() { return  1100;}
   
   static constexpr const char * name(){ return "std::vector<std::string>";}
};

//For additions:
/*template<>
struct typeDescription<>
{
   typedef  type;
   
   static  constexpr int code() { return  ;
   
   static constexpr const char * name(){ return "";}
};*/

}//namespace mx
}//namespace meta

#endif //typeDescription_hpp
