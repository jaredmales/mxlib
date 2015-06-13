/** \file fitsHeaderCard.hpp
  * \brief Declares and defines a class to work with a FITS header card
  * \ingroup fits_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
  *
  */
  
#ifndef __fitsHeaderCard_hpp__
#define __fitsHeaderCard_hpp__

#include "stringUtils.hpp"

#include "fitsUtils.hpp"

namespace mx
{

   
/** \ingroup fits_processing
  *@{
  */
      

/// Simple structure to hold the three components of a FITS header card 
/** Is a struct so all members are public for easy access.  Several string-to-number conversions are provided.
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
      type = getFitsType<char *>();
   }
   
   ///Construct from the three components.
   /**
     * \param k is the keyword
     * \param v is the value
     * \param c is the comment
     */ 
   template<typename typeT> fitsHeaderCard(std::string k, const typeT v, std::string c);
   
   ///Construct from just two components.
   /**
     * \param k is the keyword
     * \param v is the value
     */ 
   template<typename typeT> fitsHeaderCard(std::string k, const typeT v);

   //@}
   
   ///\name Conversions
   /**
     */
   //@{

   ///Get value converted to the specified type
   /** Uses \ref convertFromString.  typeT can be any <a href=\isfundamental>fundamental type</a>, 
     * or any type which can be typecast from std::string.
     */
   template<typename typeT> typeT Value();
   
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
   
   /// Writes this card to a FITS file, using \ref fits_update_key.
   /** 
     */
   int write(fitsfile * fptr);
   
   //@}
   
}; //fitsHeaderCard


template<typename typeT> fitsHeaderCard::fitsHeaderCard(std::string k, const typeT v, std::string c)
{
   keyword = k;
   value = convertToString<typeT>(v);
   type = getFitsType<typeT>();
   comment = c;
}

template<typename typeT> fitsHeaderCard::fitsHeaderCard(std::string k, const typeT v)
{
   keyword = k;
   value = convertToString<typeT>(v);
   type = getFitsType<typeT>();
}

template<typename typeT> typeT fitsHeaderCard::Value()
{
   return convertFromString<typeT>(value);
}



///@}

} //namespace mx



#endif //__fitsHeaderCard__
