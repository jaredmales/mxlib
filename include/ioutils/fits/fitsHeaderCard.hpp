/** \file fitsHeaderCard.hpp
  * \brief A class to work with a FITS header card
  * \ingroup fits_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */

//***********************************************************************//
// Copyright 2015-2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef ioutils_fits_fitsHeaderCard_hpp
#define ioutils_fits_fitsHeaderCard_hpp

#include "fitsUtils.hpp"

namespace mx
{
namespace fits 
{
   

/// Simple structure to hold the three components of a FITS header card
/** Is a struct so all members are public for easy access.  Several string-to-number conversions are provided.
  *
  * \ingroup fits_processing 
  */
struct fitsHeaderCard
{
   ///The keyword
   std::string keyword;

   ///The value
   std::string value;

   ///The comment
   std::string comment;

   ///The FITS type
   int type;

   /// \name Constructors
   /**
     */
   //@{

   ///Basic c'tor
   fitsHeaderCard()
   {
      type = getFitsType<fitsUnknownType>();
   }

   ///Construct from the three components, when value's type is unknown
   /**
     */
   fitsHeaderCard( const std::string & k, ///< [in] the keyword
                   const void * v,        ///< [in] the value, must be castable to char*
                   const std::string & c  ///< [in] the comment
                 );

   ///Construct from the three components, when value's type is unknown
   /**
     */
   fitsHeaderCard( const std::string & k, ///< [in] the keyword
                   void * v,              ///< [in] the value, must be castable to char*
                   const std::string & c  ///< [in] the comment
                 );

   ///Construct from the three components, when it's really a comment card
   /** This overload is provided to facilitate handling of comments when re-writing the file.
     *
     */
   fitsHeaderCard( const std::string & k, ///< [in] the keyword
                   fitsCommentType v,     ///< [in] an object of type fitsCommentType
                   const std::string & c  ///< [in] the comment
                 );

   ///Construct from the three components, when it's really a history card
   /** This overload is provided to facilitate handling of history when re-writing the file.
     *
     */
   fitsHeaderCard( const std::string & k, ///< [in] the keyword
                   fitsHistoryType v,     ///< [in] an object of type fitsHistoryType
                   const std::string & c  ///< [in] the comment
                 );

   ///Construct from the three components for a known type.
   /**
     */
   template<typename typeT>
   fitsHeaderCard( const std::string & k, ///< [in] they keyword
                   const typeT v,         ///< [in] the value
                   const std::string & c  ///< [in] the commend
                 );

   ///Construct from just two components, when value's type is unknown
   /**
     */
   fitsHeaderCard( const std::string & k, ///< [in] the keyword
                   void* v                ///< [in] the value, must be cast-able to char*
                 );

   ///Construct from just two components, when value's type is unknown
   /**
     */
   fitsHeaderCard( const std::string & k, ///< [in] the keyword
                   const void* v          ///< [in] the value, must be cast-able to char*
                 );

   ///Construct from just two components when value's type is known.
   /**
     */
   template<typename typeT>
   fitsHeaderCard( const std::string & k, ///< [in] the keyword
                   const typeT v          ///< [in] the value
                 );

   ///Construct from just keyword, when value's type is unknown
   /**
     */
   explicit fitsHeaderCard(const std::string & k  /**< [in] the keyword*/);


   //@}

   ///\name Setting Components
   /**
     */
   //@{

   ///Set the value string, applying the appropriate conversion and setting the FITS type.
   /**
     * \param v the value to convert to string and set
     *
     * \tparam typeT is the type to convert from
     */
   template<typename typeT>
   void setValue(const typeT v);

   //@}

   ///\name Conversions
   /**
     */
   //@{

   ///Get value converted to the specified type
   /** Uses \ref convertFromString.  typeT can be any <a href=\isfundamental>fundamental type</a>,
     * or any type which can be typecast from std::string.  This is specialized for typeT=std::string so
     * that leading and trailing apostrophes are removed.
     */
   template<typename typeT>
   typeT Value() const;

   ///Convert value to an string, stripping apostrophes if necessary.
   /** \returns the result of converting the value string to a string
    */
   std::string String();

   ///Convert value to an integer
   /** \returns the result of converting the value string to an integer
    */
   int Int();

   ///Convert value to an unsigned integer
   /** \returns the result of converting the value string to an unsigned integer
    */
   unsigned int unsignedInt();

   ///Convert value to a long integer
   /** \returns the result of converting the value string to a long integer
    */
   long Long();

   ///Convert value to an unsigned long integer
   /** \returns the result of converting the value string to an unsigned long integer
    */
   unsigned long unsignedLong();

   ///Convert value to a long long integer
   /** \returns the result of converting the value string to a long long integer
    */
   long long longLong();

   ///Convert value to an unsigned long long integer
   /** \returns the result of converting the value string to anunsigned long long  integer
    */
   unsigned long long unsignedLongLong();

   ///Convert value to a float
   /** \returns the result of converting the value string to a float
    */
   float Float();

   ///Convert value to a double
   /** \returns the result of converting the value string to a double
    */
   double Double();

   ///Convert value to a long double
   /** \returns the result of converting the value string to a long double
    */
   long double longDouble();

   //@}

   ///\name Output
   /**
     */
   //@{

   /// Writes this card to a FITS file, using \ref mx::improc::fits_write_key<typename typeT>(fitsfile * fptr, char * keyword, void * value, char * comment).
   /**
     */
   int write(fitsfile * fptr);

   //@}

}; //fitsHeaderCard



template<typename typeT>
fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                const typeT v,
                                const std::string & c
                              )
{
   keyword = k;
   value = ioutils::convertToString<typeT>(v);
   type = getFitsType<typeT>();
   comment = c;
}

template<typename typeT>
fitsHeaderCard::fitsHeaderCard( const std::string & k,
                                const typeT v
                              )
{
   keyword = k;
   value = ioutils::convertToString<typeT>(v);
   type = getFitsType<typeT>();
}

template<typename typeT>
void fitsHeaderCard::setValue(const typeT v)
{
   value = ioutils::convertToString<typeT>(v);
   type = getFitsType<typeT>();
}

template<typename typeT>
typeT fitsHeaderCard::Value() const
{
   return ioutils::convertFromString<typeT>(value);
}

} //namespace fits
} //namespace mx



#endif //ioutils_fits_fitsHeaderCard_hpp
