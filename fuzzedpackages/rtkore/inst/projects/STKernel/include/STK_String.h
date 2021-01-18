/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  STKernel::Base
 * Purpose:  Define the fundamental type String.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_String.h
 *  @brief In this file we define the fundamental type String.
 **/

#ifndef STK_STRING_H
#define STK_STRING_H

#include <string>
#include <algorithm>
#include <ctype.h>

#include "STK_Stream.h"
#include "STK_Arithmetic.h"
#include "STK_IdTypeImpl.h"

namespace STK
{
/** @ingroup Base
 *  @brief STK fundamental type of a String.
 *
 *  The type String is defined for the internal representation
 *  of the string variables (strings).
 **/
typedef std::basic_string<Char> String;

/** @ingroup Arithmetic
 *  @brief Specialization for String.
 *
 * The STK fundamental type String use the empty String to represent
 * NA String.
 */
template<>
struct Arithmetic<String>: public std::numeric_limits<String>
{
  /** Adding a Non Available (NA) special String (the empty String) */
  static String NA() throw() { return String();}
  /** True if the type has a representation for a "Not Available." */
  static const bool hasNA = true;
  /** Test if x is a Non Available (NA) String.
   *  @param x the String to test.
   *  @return @c true if @c x is a NA String, @c false otherwise
   **/
  static bool isNA(String const& x) throw() { return (x == NA());}
  /** Test if x is  infinite.
   *  @param x the String to test.
   *  @return @c false
   **/
  static bool isInfinite(String const& x) throw() { return false; }
  /** Test if x is  finite.
   *  @param x the String to test.
   *  @return @c true if @c x not a NA string, @c false otherwise
   **/
  static bool isFinite(String const& x) throw() { return (!isNA(x));}
};

/** @ingroup RTTI
 *  @brief Specialization of the IdTypeImpl for the Type String.
 *
 *  @return the IdType of a String.
 **/
template<>
struct IdTypeImpl<String>
{
  /** @return the IdType string_ */
  static Base::IdType returnType() { return(Base::string_);}
};

#ifndef IS_RTKPP_LIB
/** @ingroup Base
  * @brief Representation of a Not Available value.
  *
  * By default we represent a Not Available value of any type as a "." (like in
  * (SAS(R))) for the end-user. This value can be overloaded at runtime.
  * @note if the value is modified at runtime, the value of @c stringNaSize
  * have to be modified accordingly. It is safer to use the @ref setStringNa
  * function.
  **/
extern String stringNa;

/** @ingroup Base
  * @brief Size (in number of Char) of a Not Available value.
  * We represent a Not Available value of any type as a "." (like in
  * (SAS(R))) for the end-user.
  **/
extern int stringNaSize;

/** @ingroup Base
  * @brief Set a new value to the na String representation and modify
  * stringNaSize accordingly.
  **/
inline void setStringNa(String const& na)
{ stringNa = na; stringNaSize = na.size();}

#else /* is rtkpp lib */

/** @ingroup Base
  * @brief Representation of a Not Available value.
  *
  * We represent a Not Available value of any type as a "NA" (like in
  * (R language) for the end-user. This value cannot be overloaded at runtime.
  **/
static const String stringNa = _T("NA");
/** @ingroup Base
  * @brief Size (in number of Char) of the Not Available value.
  **/
static const int stringNaSize = 2;

#endif

/** @ingroup Base
 *  @brief convert the characters of the String to upper case
 *
 *  @param s The String to convert
 *  @return the string converted to upper case
 **/
inline String const& toUpperString( String& s)
{
  // iterate along the String
  for (String::iterator it = s.begin(); it != s.end(); ++it)
  { *it = std::toupper(*it);}
  // return upper cased string
 return s;
}

/** @ingroup Base
 *  @brief convert the characters of the String to upper case
 *
 *  @param s The String to convert
 *  @return the string converted to upper case
 **/
inline String toUpperString( String const& s)
{
  // iterate along the String
  String str = s;
  for (String::iterator it = str.begin(); it != str.end(); ++it)
  { *it = std::toupper(*it);}
  // return upper cased string
 return str;
}

/** @ingroup Base
 *  @brief convert the characters of the String to upper case
 *
 *  @param s The String to convert
 *  @return the string converted to upper case
 **/
inline String removeWhiteSpaces( String const& s)
{
  // iterate along the String
  String str = s;
  str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
  return str;
}

}

#include "STK_Proxy.h"

namespace STK
{
/** @ingroup Base
 *  @brief convert a String to Type
 *
 *  This method return true if the String @c s could be converted into
 *  a correct Type t.
 *  http://www.codeguru.com/forum/showpost.php?p=678440&postcount=1
 *  http://c.developpez.com/faq/cpp/?page=strings#STRINGS_is_type
 *
 *  @note The operator >> have been overloaded for the @c Proxy class in order
 *  to return a NA value if the conversion fail.
 *
 *  @param t The value to get from the String
 *  @param s the String to convert
 *  @param f flags
 *  @return @c true if the conversion succeed, @c false otherwise
 **/
template <class Type>
bool stringToType( Type &t, String const& s
                 , std::ios_base& (*f)(std::ios_base&) = std::dec
                 )
{ istringstream iss(s);
  return !( iss >> f >> Proxy<Type>(t)).fail();
}

/** @ingroup Base
 *  @brief convert a String to Type without error check
 *
 *  This method return the String @c s converted into a correct Type t
 *  without formatting.
 *  http://www.codeguru.com/forum/showpost.php?p=678440&postcount=1
 *  http://c.developpez.com/faq/cpp/?page=strings#STRINGS_is_type
 *
 *  @note if the conversion fail, the method return a NA value if available
 *  for this @c Type.
 *  @param s the String to convert
 *  @return The value to get from the String
 **/
template <class Type>
Type stringToType( String const& s)
{
  Type t;
  istringstream iss(s);
  iss >> Proxy<Type>(t);
  return(t);
}

/** @ingroup Base
 *  @brief convert a Type to String
 *
 *  This method return the Type t into a String s.
 *  @see http://www.codeguru.com/forum/showpost.php?p=678440&postcount=1
 *  @see http://c.developpez.com/faq/cpp/?page=strings#STRINGS_convertform
 *
 *  @param t The value to convert to String
 *  @param f flag, by default write every number in decimal
 **/
template <class Type>
String typeToString( Type const& t, std::ios_base& (*f)(std::ios_base&) = std::dec)
{
  if (Arithmetic<Type>::isNA(t)) return stringNa;
  ostringstream oss;
  return static_cast<ostringstream&>(oss << f << Proxy<Type>(t)).str();
}

} // namespace STK

#endif /*STK_STRING_H*/
