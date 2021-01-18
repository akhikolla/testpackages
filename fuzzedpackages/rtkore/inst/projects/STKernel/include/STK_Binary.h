/*--------------------------------------------------------------------*/
/*     Copyright (C) 2007-2015  Serge Iovleff

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
 * Project:  Base
 * Purpose:  Define the fundamental type Binary.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Binary.h
 *  @brief In this file we define the fundamental type Binary.
 **/

#ifndef STK_BINARY_H
#define STK_BINARY_H

#include <map>
#include "STK_String.h"

namespace STK
{
/** @ingroup Base
 *  @brief STK fundamental type of a binary.
 *
 *  The type Binary is a representation of dichotomic variables.
 **/
enum Binary
{ zero_ =0, ///< 0 value
  one_  =1, ///< 1 value
  binaryNA_=INT_MIN ///< Not Available value
};

 /** @ingroup Arithmetic
  *  @brief Specialization for Binary.
  *  NA (not available) numbers is part of the @c enum Binary.
  */
 template<>
 struct Arithmetic<Binary> : public std::numeric_limits<Binary>
 {
   /** Adding a Non Available (NA) special number. */
   static inline Binary NA() throw() { return binaryNA_;}
   /** True if the type has a representation for a "Not Available". */
   static const bool hasNA = true;
   /** Test if x is a Non Available (NA) special number
    *  @param x the Binary number to test.
    **/
   static inline bool isNA(const Binary& x) throw() { return (x==binaryNA_);}
   /** test if x is  infinite.
    *  @param x the Binary number to test.
    **/
   static inline bool isInfinite(const Binary& x) throw() { return false; }
   /** test if x is  finite.
    *  @param x the Binary number to test.
    **/
   static inline bool isFinite(const Binary& x) throw() { return (!isNA(x) && !isInfinite(x));}
 };

/** @ingroup RTTI
 *  @brief Specialization of the IdTypeImpl for the Type Binary.
 **/
template<>
struct IdTypeImpl<Binary>
{
  /** @return the IdType of the type Binary. */
  static inline Base::IdType returnType() { return(Base::binary_);}
};
/** @ingroup stream
 *  @brief Overloading of the ostream << for the type Binary.
 *  @param os the output stream
 *  @param value the value to send to the stream
 **/
inline ostream& operator << (ostream& os, Binary const& value)
{
  return Arithmetic<Binary>::isNA(value) ? (os <<  stringNa)
                                         : (os << static_cast<int>(value));
}

/** @ingroup stream
 *  @brief Overloading of the istream >> for the type Binary.
 *  @param is the input stream
 *  @param value the value to get from the stream
 **/
inline istream& operator >> (istream& is, Binary& value)
{
  int res;
  // try to read an integer
  if (!(is >> res).fail())
  {
    switch (res)
    {
      case 0:
        value = zero_;
        break;
      case 1:
        value = one_;
        break;
      default:
        value = binaryNA_;
        break;
    }
  }
  else
  { value = binaryNA_; is.clear(); is.setstate(std::ios::failbit);}
  return is;
}

/** @ingroup Base
 *  Convert a String to a Binary.
 *  @param str the String we want to convert
 *  @return the Binary represented by the String @c str. If the string
 *  does not match any known name, the @c binaryNA_ value  is returned.
 **/
inline Binary stringToBinary( String const& str)
{
  if (toUpperString(str) == toUpperString(_T("0"))) return zero_;
  if (toUpperString(str) == toUpperString(_T("1"))) return one_;
  return binaryNA_;
}

/** @ingroup Base
 *  Convert a String to a Binary using a map.
 *  @param str the String we want to convert
 *  @param mapping the mapping between the String and the Binary
 *  @return the Binary represented by the String @c str. If the string
 *  does not match any known name, the @c binaryNA_ type is returned.
 **/
inline Binary stringToBinary( String const& str, std::map<String, Binary> const& mapping)
{
  std::map<String, Binary>::const_iterator it=mapping.find(str);
  return (it == mapping.end()) ? binaryNA_ : it->second;
}

/** @ingroup Base
 *  Convert a Binary to a String.
 *  @param value the Binary we want to convert
 *  @param f format, by default write every number in decimal
 *  @return the string associated to this type.
 **/
inline String binaryToString( Binary const& value, std::ios_base& (*f)(std::ios_base&) = std::dec)
{
  if (Arithmetic<Binary>::isNA(value)) return stringNa;
  ostringstream os;
  os << f << static_cast<int>(value);
  return os.str();
}

/** @ingroup Base
 *  Convert a Binary to a String.
 *  @param value the Binary we want to convert
 *  @param mapping the mapping between the Binary and the String
 *  @return the string associated to this value.
 **/
inline String binaryToString( Binary const& value, std::map<Binary, String> const& mapping)
{
  std::map<Binary, String>::const_iterator it=mapping.find(value);
  return (it == mapping.end()) ? Arithmetic<String>::NA() : it->second;
}

/** @ingroup Base
 *  @brief specialization for Binary
 *  @param str the String to convert
 *  @return The value to get from the String
 **/
template<>
inline Binary stringToType<Binary>( String const& str)
{ return stringToBinary(str);}

/** @ingroup Base
 *  @brief specialization for Binary
 *  @param t The Binary to convert to String
 *  @param f format, by default write every number in decimal
 **/
template<>
inline String typeToString<Binary>( Binary const& t, std::ios_base& (*f)(std::ios_base&))
{ return binaryToString(t);}

} // namespace STK

#endif /*STK_BINARY_H*/
