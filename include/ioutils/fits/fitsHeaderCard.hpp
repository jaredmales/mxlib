/** \file fitsHeaderCard.hpp
 * \brief A class to work with a FITS header card
 * \ingroup fits_processing_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2015-2022 Jared R. Males (jaredmales@gmail.com)
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

#include "../../mxException.hpp"
#include "fitsUtils.hpp"

namespace mx
{
   namespace fits
   {

      /// Class to manage the three components of a FITS header card
      /** Since FITS does not provide the type in keyword=value pairs in a FITS header, it is up to the user
       * to determine the type.  Futhermore, since we want to read values from files, type must
       * be set at runtime.   The result is that we must be able to accept a string, which is converted
       * to a given type on demand at determined at runtime.
       *
       * Conversion from string to native type, or vice versa, only occurs when needed.  So if you set the value to,
       * say, a double, the value is not converted to string format unless specifically requested.  If the write function is called
       * when in this state, the cfitsio routine is called directly.  This converson only on demand is most important for values read
       * from a file, then written to another file.  In this case, no conversion to its double (etc) representation.
       * occurs.
       *
       * Note that because of the checks to determine the type and appropriate return values, accessing the value in a card
       * is possibly slower than accessing a variable due to various if statements.  This means that you should typically do 
       * so once and use a local variable for repeated use.
       *
       * \ingroup fits_processing
       */
      class fitsHeaderCard
      {

      protected:
         /// The keyword
         std::string m_keyword;

         /// The FITS type of the value, and indicates which member of m_values to access.
         int m_type{fitsType<fitsUnknownType>()};

         /// The native type is held in a union.
         union values
         {
            char Char;                          ///< the char value
            unsigned char UChar;                ///< the unsigned char value
            short Short;                        ///< the short value
            unsigned short UShort;              ///< the unsigned short value
            int Int;                            ///< the int value
            unsigned int UInt;                  ///< the unsigned int value
            long Long;                          ///< the long value
            unsigned long ULong;                ///< the unsigned long value
            long long LongLong;                 ///< the long long value
            unsigned long long ULongLong;       ///< the unsigned long long value
            float Float;                        ///< the float value
            std::complex<float> complexFloat;   ///< the std::complex<float> value
            double Double;                      ///< the double value
            std::complex<double> complexDouble; ///< the std::complex<double> value

            /// c'tor.  have to specify due to inclusion of std::complex types.
            values()
            {
               return;
            }

         } m_value;

         /// The value in string form
         std::stringstream m_valueStr;

         bool m_valueGood{false};    ///< Flag indicating if the value is valid
         bool m_valueStrGood{false}; ///< Flag indicating if the value string is valid

         /// The comment
         std::string m_comment;

      public:
         /// \name Constructors
         /**
          */
         //@{

         /// Basic c'tor
         fitsHeaderCard()
         {
         }

         /// Construct from the three components for a value of string type
         /**
          */
         fitsHeaderCard(const std::string &k,     ///< [in] the keyword
                        const std::string &v,     ///< [in] the value string
                        const std::string &c = "" ///< [in] the comment
         );

         /// Construct from the three components, when already in a string format
         /** Use this when the value is not a string
          */
         fitsHeaderCard(const std::string &k,     ///< [in] the keyword
                        const std::string &v,     ///< [in] the value string
                        const int &type,          ///< [in] the type of the value
                        const std::string &c = "" ///< [in] the comment
         );

         /// Construct from the three components, when it's really a comment card
         /** This overload is provided to facilitate handling of comments when re-writing the file.
          *
          */
         fitsHeaderCard(const std::string &k, ///< [in] the keyword
                        fitsCommentType v,    ///< [in] an object of type fitsCommentType
                        const std::string &c  ///< [in] the comment
         );

         /// Construct from the three components, when it's really a history card
         /** This overload is provided to facilitate handling of history when re-writing the file.
          *
          */
         fitsHeaderCard(const std::string &k, ///< [in] the keyword
                        fitsHistoryType v,    ///< [in] an object of type fitsHistoryType
                        const std::string &c  ///< [in] the comment
         );

         /// Construct from just keyword, when value's type is unknown
         /**
          */
         explicit fitsHeaderCard(const std::string &k /**< [in] the keyword*/);

         /// Construct from just keyword, when value's type known
         /**
          */
         fitsHeaderCard(const std::string &k, ///< [in] the keyword
                        const int type        ///< [in] the type
         );

         /// Construct from the three components for a char.
         /**
          */
         template <typename typeT>
         fitsHeaderCard(const std::string &k,     ///< [in] they keyword
                        const typeT v,            ///< [in] the value
                        const std::string &c = "" ///< [in] the comment
         );

         /// Copy constructor
         fitsHeaderCard(const fitsHeaderCard &card);

         //@}

         fitsHeaderCard &operator=(const fitsHeaderCard &card);

      protected:
         ///\name Converters
         /** @{
          */

         /// Convert from the type to a string.
         /** This populates m_valueStr and sets m_valueStrGood so that this conversion
          * only occurs once.
          */
         void convertToString();

         /// Convert from string to the type
         /** This populates the appropriate union field and sets m_valueGood so that
          * this conversion only occurs once.
          */
         template <typename typeT>
         void convertFromString();

         /// Get the value from its type converted to a different type.
         template <typename typeT>
         typeT convertedValue();

         /// Convert the value from its type to a different type.
         void convertValue(int newtype /**< [in] the new type */);

         ///@}

      public:
         ///\name Accessors
         /** @{
          */

         /// Get the keyword
         /** \returns a const reference to m_keyword
          */
         const std::string &keyword();

         /// Set the keyword
         void keyword(const std::string &kw /**< [in] the new keyword */);

         /// Get the type
         /** \returns the value of m_type
          */
         int type();

         /// Set the type
         /** If this is a change in type and the native type is set in m_value (indicated by m_valueGood == true)
          * then it is converted to the new type.  Otherwise, no conversion occurs.
          */
         void type(const int &t /**< [in] the new type */);

         /// Get the value
         /** Returns the value as typeT.  Conversions occur
          * automatically if necessary.
          *
          * \returns the value converted to typeT as necessary
          *
          * \throws if the value can't be converted to typeT
          */
         template <typename typeT>
         typeT value();

         /// Get the value as a string
         /** This calls value<string>().
          *
          * \returns the value converted to string as necessary
          *
          * \throws if the value can't be converted to a string
          */
         std::string String();

         /// Get the value as a char
         /** This calls value<char>().
          *
          * \returns the value converted to char as necessary
          *
          * \throws if the value can't be converted to char
          */
         char Char();

         /// Get the value as an unsigned char
         /** This calls value<unsigned char>().
          *
          * \returns the value converted to unsigned char as necessary
          *
          * \throws if the value can't be converted to unsigned char
          */
         unsigned char UChar();

         /// Get the value as a short
         /** This calls value<short>().
          *
          * \returns the value converted to short as necessary
          *
          * \throws if the value can't be converted to short
          */
         short Short();

         /// Get the value as an unsigned short
         /** This calls value<unsigned short>().
          *
          * \returns the value converted to unsigned short as necessary
          *
          * \throws if the value can't be converted to unsigned short
          */
         unsigned short UShort();

         /// Get the value as a int
         /** This calls value<int>().
          *
          * \returns the value converted to int as necessary
          *
          * \throws if the value can't be converted to int
          */
         int Int();

         /// Get the value as an unsigned int
         /** This calls value<unsigned int>().
          *
          * \returns the value converted to unsigned int as necessary
          *
          * \throws if the value can't be converted to unsigned int
          */
         unsigned int UInt();

         /// Get the value as a long
         /** This calls value<long>().
          *
          * \returns the value converted to long as necessary
          *
          * \throws if the value can't be converted to long
          */
         long Long();

         /// Get the value as an unsigned long
         /** This calls value<unsigned long>().
          *
          * \returns the value converted to unsigned long as necessary
          *
          * \throws if the value can't be converted to unsigned long
          */
         unsigned long ULong();

         /// Get the value as a long long
         /** This calls value<long long>().
          *
          * \returns the value converted to long long as necessary
          *
          * \throws if the value can't be converted to long long
          */
         long long LongLong();

         /// Get the value as an unsigned long long
         /** This calls value<unsigned long long>().
          *
          * \returns the value converted to unsigned long long as necessary
          *
          * \throws if the value can't be converted to unsigned long long
          */
         unsigned long long ULongLong();

         /// Get the value as a float
         /** This calls value<float>().
          *
          * \returns the value converted to float as necessary
          *
          * \throws if the value can't be converted to float
          */
         float Float();

         /// Get the value as a std::complex<float>
         /** This calls value<std::complex<float>>().
          *
          * \returns the value converted to std::complex<float> as necessary
          *
          * \throws if the value can't be converted to std::complex<float>
          */
         std::complex<float> complexFloat();

         /// Get the value as a double
         /** This calls value<double>().
          *
          * \returns the value converted to double as necessary
          *
          * \throws if the value can't be converted to double
          */
         double Double();

         /// Get the value as a std::complex<double>
         /** This calls value<std::complex<double>>().
          *
          * \returns the value converted to std::complex<double> as necessary
          *
          * \throws if the value can't be converted to std::complex<double>
          */
         std::complex<double> complexDouble();

         /// Set the value to a char * string
         void value(const char *v /**< [in] a character string*/);

         /// Set the value to a std::string
         /** \overload
          */
         void value(const std::string &v /**< [in] a std::string*/);

         /// Set the value to a char
         /** \overload
          */
         void value(const char &v /**< [in] a char*/);

         /// Set the value to an unsigned char
         /** \overload
          */
         void value(const unsigned char &v /**< [in] an unsigned char */);

         /// Set the value to a short int
         /** \overload
          */
         void value(const short int &v /**< [in] a short int*/);

         /// Set the value to an unsigned short int
         /** \overload
          */
         void value(const unsigned short int &v /**< [in] an unsigned short int*/);

         /// Set the value to an int
         /** \overload
          */
         void value(const int &v /**< [in] an int*/);

         /// Set the value to an unsigned int
         /** \overload
          */
         void value(const unsigned int &v /**< [in] an unsigned int*/);

         /// Set the value to a long int
         /** \overload
          */
         void value(const long &v /**< [in] a long int*/);

         /// Set the value to an unsigned long int
         /** \overload
          */
         void value(const unsigned long int &v /**< [in] an unsigned long int*/);

         /// Set the value to a long long int
         /** \overload
          */
         void value(const long long &v /**< [in] a long long int*/);

         /// Set the value to an unsigned long long int
         /** \overload
          */
         void value(const unsigned long long int &v /**< [in] an unsigned long long int*/);

         /// Set the value to a float
         /** \overload
          */
         void value(const float &v /**< [in] a float*/);

         /// Set the value to a complex float
         /** \overload
          */
         void value(const std::complex<float> &v /**< [in] a complex float*/);

         /// Set the value to a double
         /** \overload
          */
         void value(const double &v /**< [in] a double*/);

         /// Set the value to a complex double
         /** \overload
          */
         void value(const std::complex<double> &v /**< [in] a complex double*/);

         std::string valueStr();

         bool valueGood();

         bool valueStrGood();

         /// Get the comment
         /** \returns the value of m_comment
          */
         const std::string &comment();

         /// Set the comment
         void comment(const std::string &c /**< [in] the new comment */);

         //@}

         ///\name Output
         /**
          */
         //@{

         /// Writes this card to a FITS file, using \ref mx::improc::fits_write_key<typename typeT>(fitsfile * fptr, char * keyword, void * value, char * comment).
         /**
          */
         int write(fitsfile *fptr);

         //@}

      }; // fitsHeaderCard

      template <typename typeT>
      fitsHeaderCard::fitsHeaderCard(const std::string &k,
                                     const typeT v,
                                     const std::string &c)
      {
         m_keyword = k;
         value(v);
         m_comment = c;
      }

   } // namespace fits
} // namespace mx

#endif // ioutils_fits_fitsHeaderCard_hpp
