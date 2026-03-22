/** \file fitsHeaderCard.hpp
 * \brief A class to work with a FITS header card
 * \ingroup fits_processing_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2015-2025 Jared R. Males (jaredmales@gmail.com)
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

#include <format>

#include "../../mxlib.hpp"

#include "../../meta/tagT.hpp"
#include "fitsUtils.hpp"

namespace mx
{
namespace fits
{

/// \cond
// strip leading and trailing whitespace and then opening and closing ''.  leaves spaces between ''.
inline void stripApostWS( std::string &str )
{
    if( str.size() == 0 )
    {
        return;
    }

    if( str[0] != '\'' && str[0] != ' ' && str.back() != ' ' ) // get out fast if we can
    {
        return;
    }

    // strip white space at front
    size_t ns = str.find_first_not_of( " \t\r\n" );
    if( ns != std::string::npos && ns != 0 )
    {
        str.erase( 0, ns );

        if( str.size() == 0 )
        {
            return;
        }
    }
    else if( ns == std::string::npos ) // the rare all spaces
    {
        str = "";
        return;
    }

    // strip white space at back
    ns = str.find_last_not_of( " \t\r\n" );
    if( ns != std::string::npos && ns != str.size() - 1 )
    {
        str.erase( ns + 1 );

        if( str.size() == 0 )
        {
            return;
        }
    }

    if( str[0] == '\'' && str.back() == '\'' )
    {
        if( str.size() == 1 || str.size() == 2 )
        {
            str = "";
            return;
        }
        str.erase( str.size() - 1, 1 );
        str.erase( 0, 1 );
    }
}
/// \endcond

/// Class to manage the three components of a FITS header card
/** Since FITS does not provide the type in keyword=value pairs in a FITS header, it is up to the user
 * to determine the type.  Furthermore, since we want to read values from files, type conversions must
 * be done at runtime.   The result is that we must be able to accept a string, which is converted
 * to a given type on demand as determined at runtime.
 *
 * Conversion from string to native type, or vice versa, only occurs when needed.  So if you set the value to,
 * say, a double, the value is not converted to string format unless specifically requested.  If the write function is
 * called when in this state, the cfitsio routine is called directly.  This conversion only on demand is most important
 * for values read from a file, then written to another file.  In this case, no conversion to its double (etc)
 * representation occurs.
 *
 * Note that because of the checks to determine the type and appropriate return values, accessing the value in a card
 * is possibly slower than accessing a variable due to various if statements and error checking.
 * This means that you should typically do
 * so once and use a local variable for repeated use.
 *
 * \ingroup fits_processing
 */
template <class verboseT = verbose::d>
class fitsHeaderCard
{

  protected:
    /// The keyword
    std::string m_keyword;

    /// The FITS type of the value, and indicates which member of m_values to access.
    int m_type{ fitsType<fitsUnknownType>() };

    /// The native type is held in a union.
    union values
    {
        bool Bool;                                   ///< the bool value
        char Char;                                   ///< the char value
        unsigned char UChar;                         ///< the unsigned char value
        short Short;                                 ///< the short value
        unsigned short UShort;                       ///< the unsigned short value
        int Int;                                     ///< the int value
        unsigned int UInt;                           ///< the unsigned int value
        long Long;                                   ///< the long value
        unsigned long ULong;                         ///< the unsigned long value
        long long LongLong;                          ///< the long long value
        unsigned long long ULongLong;                ///< the unsigned long long value
        float Float;                                 ///< the float value
        std::complex<float> complexFloat;            ///< the std::complex<float> value
        double Double;                               ///< the double value
        std::complex<double> complexDouble;          ///< the std::complex<double> value
        long double LongDouble;                      ///< the long double value
        std::complex<long double> complexLongDouble; ///< the std::complex<long double> value

        // clang-format off
        #ifdef MXLIB_HAS_QUAD
        __float128 Quad;                      ///< the float128 value
        std::complex<__float128> complexQuad; ///< the std::complex<__float128> value
        #endif
        // clang-format on

        /// c'tor.  have to specify due to inclusion of std::complex types.
        values()
        {
            return;
        }

        bool &member( meta::tagT<bool> )
        {
            return Bool;
        }

        char &member( meta::tagT<char> )
        {
            return Char;
        }

        unsigned char &member( meta::tagT<unsigned char> )
        {
            return UChar;
        }

        short &member( meta::tagT<short> )
        {
            return Short;
        }

        unsigned short &member( meta::tagT<unsigned short> )
        {
            return UShort;
        }

        int &member( meta::tagT<int> )
        {
            return Int;
        }

        unsigned int &member( meta::tagT<unsigned int> )
        {
            return UInt;
        }

        long &member( meta::tagT<long> )
        {
            return Long;
        }

        unsigned long &member( meta::tagT<unsigned long> )
        {
            return ULong;
        }

        long long &member( meta::tagT<long long> )
        {
            return LongLong;
        }

        unsigned long long &member( meta::tagT<unsigned long long> )
        {
            return ULongLong;
        }

        float &member( meta::tagT<float> )
        {
            return Float;
        }

        std::complex<float> &member( meta::tagT<std::complex<float>> )
        {
            return complexFloat;
        }

        double &member( meta::tagT<double> )
        {
            return Double;
        }

        std::complex<double> &member( meta::tagT<std::complex<double>> )
        {
            return complexDouble;
        }

        long double &member( meta::tagT<long double> )
        {
            return LongDouble;
        }

        std::complex<long double> &member( meta::tagT<std::complex<long double>> )
        {
            return complexLongDouble;
        }

        // clang-format off
        #ifdef MXLIB_HAS_QUAD
        __float128 &member( meta::tagT<__float128> )
        {
            return Quad;
        }

        std::complex<__float128> &member( meta::tagT<std::complex<__float128>> )
        {
            return complexQuad;
        }
        #endif
        // clang-format on

        template <typename typeT>
        typeT &member()
        {
            return member( meta::tagT<typeT>() );
        }

    } m_value;

    /// The value in string form
    std::stringstream m_valueStr;

    bool m_valueGood{ false };    ///< Flag indicating if the value is valid
    bool m_valueStrGood{ false }; ///< Flag indicating if the value string is valid

    /// The comment
    std::string m_comment;

  public:
    /// \name Constructors
    /**
     */
    //@{

    /// Basic c'tor
    fitsHeaderCard();

    /// Construct from the three components for a value of string type
    /**
     */
    fitsHeaderCard( const std::string &k,     ///< [in] the keyword
                    const std::string &v,     ///< [in] the value string
                    const std::string &c = "" ///< [in] the comment
    );

    /// Construct from the three components for a value of string type
    /** Have to provide overload for char * to avoid template version
     */
    fitsHeaderCard( const std::string &k,     ///< [in] the keyword
                    char *v,                  ///< [in] the value string
                    const std::string &c = "" ///< [in] the comment
    );

    /// Construct from the three components for a value of string type
    /** Have to provide overload for const char * to avoid template version
     */
    fitsHeaderCard( const std::string &k,     ///< [in] the keyword
                    const char *v,            ///< [in] the value string
                    const std::string &c = "" ///< [in] the comment
    );

    /// Construct from the three components, when already in a string format
    /** Use this when the value is not a string
     */
    fitsHeaderCard( const std::string &k,     ///< [in] the keyword
                    const std::string &v,     ///< [in] the value string
                    const int &type,          ///< [in] the type of the value
                    const std::string &c = "" ///< [in] the comment
    );

    /// Construct from the three components, when it's really a comment card
    /** This overload is provided to facilitate handling of comments when re-writing the file.
     *
     */
    fitsHeaderCard( const std::string &k, ///< [in] the keyword
                    fitsCommentType v,    ///< [in] an object of type fitsCommentType
                    const std::string &c  ///< [in] the comment
    );

    /// Construct from the three components, when it's really a history card
    /** This overload is provided to facilitate handling of history when re-writing the file.
     *
     */
    fitsHeaderCard( const std::string &k, ///< [in] the keyword
                    fitsHistoryType v,    ///< [in] an object of type fitsHistoryType
                    const std::string &c  ///< [in] the comment
    );

    /// Construct from just keyword, when value's type is unknown
    /**
     */
    explicit fitsHeaderCard( const std::string &k /**< [in] the keyword*/ );

    /// Construct from just keyword, when value's type known
    /**
     */
    fitsHeaderCard( const std::string &k, ///< [in] the keyword
                    const int type        ///< [in] the type
    );

    /// Construct from the three components for a char.
    /**
     */
    template <typename typeT>
    fitsHeaderCard( const std::string &k,     ///< [in] they keyword
                    const typeT v,            ///< [in] the value
                    const std::string &c = "" ///< [in] the comment
    );

    /// Copy constructor
    fitsHeaderCard( const fitsHeaderCard &card );

    //@}

    /// Assignment
    fitsHeaderCard &operator=( const fitsHeaderCard &card );

  protected:
    ///\name Converters
    /** @{
     */

    /// Convert from the type to a string.
    /** This populates m_valueStr and sets m_valueStrGood so that this conversion
     * only occurs once.
     */
    mx::error_t convertToString();

    /// Convert from string to the type
    /** This populates the appropriate union field and sets m_valueGood so that
     * this conversion only occurs once.
     */
    template <typename typeT>
    error_t convertFromString();

    /// Get the value from its type converted to a different type.
    template <typename typeT>
    error_t convertedValue( typeT &cval /**< [out] the converted value */ );

    /// Convert the value from its type to a different type.
    error_t convertValue( int newtype /**< [in] the new type */ );

    ///@}

  public:
    ///\name Accessors
    /** @{
     */

    /// Get the keyword
    /**
     * \returns a const reference to m_keyword
     */
    const std::string &keyword() const;

    /// Set the keyword
    /**
     * \returns error_t::noerror on success
     */
    error_t keyword( const std::string &kw /**< [in] the new keyword */ );

    /// Get the type
    /**
     * \returns the value of m_type
     */
    int type() const;

    /// Set the type
    /** If this is a change in type and the native type is set in m_value (indicated by m_valueGood == true)
     * then it is converted to the new type.  Otherwise, no conversion occurs.
     *
     * \returns error_t::noerror on success
     */
    error_t type( const int &t /**< [in] the new type */ );

  protected:
    /** \name tag dispatching for getting value
     * @{
     */

    /// Get value for anything not a string
    template <typename typeT>
    typeT valueNonString( mx::error_t &errc );

    /// Get value tag dispatcher for std::string
    /**
     */
    std::string value( meta::tagT<std::string>, mx::error_t &errc );

    /// Get value tag dispatcher for char
    /** Calls valueNonString
     */
    char value( meta::tagT<char>, mx::error_t &errc );

    /// Get value tag dispatcher for unsigned char
    /** Calls valueNonString
     */
    unsigned char value( meta::tagT<unsigned char>, mx::error_t &errc );

    /// Get value tag dispatcher for short
    /** Calls valueNonString
     */
    short value( meta::tagT<short>, mx::error_t &errc );

    /// Get value tag dispatcher for unsigned short
    /** Calls valueNonString
     */
    unsigned short value( meta::tagT<unsigned short>, mx::error_t &errc );

    /// Get value tag dispatcher for int
    /** Calls valueNonString
     */
    int value( meta::tagT<int>, mx::error_t &errc );

    /// Get value tag dispatcher for unsigned int
    /** Calls valueNonString
     */
    unsigned int value( meta::tagT<unsigned int>, mx::error_t &errc );

    /// Get value tag dispatcher for long
    /** Calls valueNonString
     */
    long value( meta::tagT<long>, mx::error_t &errc );

    /// Get value tag dispatcher for unsigned long
    /** Calls valueNonString
     */
    unsigned long value( meta::tagT<unsigned long>, mx::error_t &errc );

    /// Get value tag dispatcher for long long
    /** Calls valueNonString
     */
    long long value( meta::tagT<long long>, mx::error_t &errc );

    /// Get value tag dispatcher for unsigned long long
    /** Calls valueNonString
     */
    unsigned long long value( meta::tagT<unsigned long long>, mx::error_t &errc );

    /// Get value tag dispatcher for float
    /** Calls valueNonString
     */
    float value( meta::tagT<float>, mx::error_t &errc );

    /// Get value tag dispatcher for double
    /** Calls valueNonString
     */
    double value( meta::tagT<double>, mx::error_t &errc );

    ///@}

  public:
    /// Get the value
    /** Returns the value as typeT.  Conversions occur
     * automatically if necessary.
     *
     * \returns the value converted to typeT as necessary
     *
     * \b Errors
     * - mx::error_t::noerror on success
     * - other errors possible due to conversions
     *
     */
    template <typename typeT>
    typeT value( mx::error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a string
    /** This calls value<string>().
     *
     * \returns the value converted to string as necessary
     *
     */
    std::string String( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a char
    /** This calls value<char>().
     *
     * \returns the value converted to char as necessary
     *
     */
    char Char( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as an unsigned char
    /** This calls value<unsigned char>().
     *
     * \returns the value converted to unsigned char as necessary
     *
     */
    unsigned char UChar( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a short
    /** This calls value<short>().
     *
     * \returns the value converted to short as necessary
     *
     */
    short Short( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as an unsigned short
    /** This calls value<unsigned short>().
     *
     * \returns the value converted to unsigned short as necessary
     *
     */
    unsigned short UShort( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a int
    /** This calls value<int>().
     *
     * \returns the value converted to int as necessary
     *
     */
    int Int( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as an unsigned int
    /** This calls value<unsigned int>().
     *
     * \returns the value converted to unsigned int as necessary
     *
     */
    unsigned int UInt( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a long
    /** This calls value<long>().
     *
     * \returns the value converted to long as necessary
     *
     */
    long Long( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as an unsigned long
    /** This calls value<unsigned long>().
     *
     * \returns the value converted to unsigned long as necessary
     *
     */
    unsigned long ULong( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a long long
    /** This calls value<long long>().
     *
     * \returns the value converted to long long as necessary
     *
     */
    long long LongLong( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as an unsigned long long
    /** This calls value<unsigned long long>().
     *
     * \returns the value converted to unsigned long long as necessaryvalue(
     *
     */
    unsigned long long ULongLong( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a float
    /** This calls value<float>().
     *
     * \returns the value converted to float as necessary
     *
     */
    float Float( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a std::complex<float>
    /** This calls value<std::complex<float>>().
     *
     * \returns the value converted to std::complex<float> as necessary
     *
     */
    std::complex<float> complexFloat( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a double
    /** This calls value<double>().
     *
     * \returns the value converted to double as necessary
     *
     */
    double Double( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Get the value as a std::complex<double>
    /** This calls value<std::complex<double>>().
     *
     * \returns the value converted to std::complex<double> as necessary
     *
     */
    std::complex<double> complexDouble( error_t *errc = nullptr /**< [in] [optional] error code */ );

    /// Set the value to a char * string
    error_t value( const char *v /**< [in] a character string*/ );

    /// Set the value to a std::string
    /**
     */
    error_t value( const std::string &v /**< [in] a std::string*/ );

    /// Set the value for a non-string type
    template <typename typeT>
    mx::error_t value( const typeT &v /**< [in] the value to set */ );

    /// Get the current value string
    /**
     * \returns m_valueStr;
     */
    std::string valueStr();

    /// Get the current value good flag
    /**
     * \returns m_valueGood;
     */
    bool valueGood();

    /// Get the current value string good flag
    /**
     * \returns m_valueStrGood;
     */
    bool valueStrGood();

    /// Get the comment
    /** \returns the value of m_comment
     */
    const std::string &comment();

    /// Set the comment
    /**
     * \returns error_t::noerror on success
     */
    error_t comment( const std::string &c /**< [in] the new comment */ );

    //@}

    error_t appendContinue( const fitsHeaderCard &card );

    ///\name Output
    /**
     */
    //@{

    /// Writes this card to a FITS file, using \ref mx::improc::fits_write_key<typename typeT>(fitsfile * fptr, char *
    /// keyword, void * value, char * comment).
    /**
     */
    mx::error_t write( fitsfile *fptr );

    //@}

}; // fitsHeaderCard

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard()
{
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k, const std::string &v, const std::string &c )
{
    m_keyword = k;

    std::string str = v;
    stripApostWS( str );
    m_valueStr.str( str );
    m_valueGood = false;
    m_valueStrGood = true;
    m_type = fitsType<std::string>();
    m_comment = c;
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k, char *v, const std::string &c )
{
    m_keyword = k;
    std::string str = v;
    stripApostWS( str );
    m_valueStr.str( str );
    m_valueGood = false;
    m_valueStrGood = true;
    m_type = fitsType<std::string>();
    m_comment = c;
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k, const char *v, const std::string &c )
{
    m_keyword = k;
    std::string str = v;
    stripApostWS( str );
    m_valueStr.str( str );
    m_valueGood = false;
    m_valueStrGood = true;
    m_type = fitsType<std::string>();
    m_comment = c;
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k,
                                          const std::string &v,
                                          const int &type,
                                          const std::string &c )
{
    m_keyword = k;
    std::string str = v;
    stripApostWS( str );
    m_valueStr.str( str );
    m_valueGood = false;
    m_valueStrGood = true;
    m_type = type;
    m_comment = c;
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k, fitsCommentType v, const std::string &c )
{
    m_keyword = k;
    m_valueGood = false;
    m_valueStrGood = false;
    m_type = fitsType<fitsCommentType>();
    m_comment = c;
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k, fitsHistoryType v, const std::string &c )
{
    m_keyword = k;
    m_valueGood = false;
    m_valueStrGood = false;
    m_type = fitsType<fitsHistoryType>();
    m_comment = c;
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k )
{
    m_keyword = k;
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k, const int type )
{
    m_keyword = k;
    m_type = type;
}

template <class verboseT>
template <typename typeT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const std::string &k, const typeT v, const std::string &c )
{
    m_keyword = k;
    value( v );
    m_comment = c;
}

template <class verboseT>
fitsHeaderCard<verboseT>::fitsHeaderCard( const fitsHeaderCard &card )
{
    m_keyword = card.m_keyword;
    m_type = card.m_type;
    memcpy( &m_value, &card.m_value, sizeof( values ) );
    m_valueStr.str( card.m_valueStr.str() );
    m_valueGood = card.m_valueGood;
    m_valueStrGood = card.m_valueStrGood;
    m_comment = card.m_comment;
}

template <class verboseT>
fitsHeaderCard<verboseT> &fitsHeaderCard<verboseT>::operator=( const fitsHeaderCard &card )
{
    m_keyword = card.m_keyword;
    m_type = card.m_type;
    memcpy( &m_value, &card.m_value, sizeof( values ) );
    m_valueStr.str( card.m_valueStr.str() );
    m_valueGood = card.m_valueGood;
    m_valueStrGood = card.m_valueStrGood;
    m_comment = card.m_comment;

    return *this;
}

template <class verboseT>
mx::error_t fitsHeaderCard<verboseT>::convertToString()
{
    if( !m_valueGood )
    {
        return mx::error_t::paramnotset;
    }

    if( m_type == fitsType<char *>() || m_type == fitsType<std::string>() )
    {
        m_valueStrGood = true; // It should be hard to get here, but just in case.
        return mx::error_t::noerror;
    }

    m_valueStr.str( "" );
    m_valueStr.precision( 10 );

    switch( m_type )
    {
        case fitsType<char>():
            m_valueStr << static_cast<int>( m_value.Char );
            break;
        case fitsType<unsigned char>():
            m_valueStr << static_cast<int>( m_value.UChar );
            break;
        case fitsType<short>():
            m_valueStr << m_value.Short;
            break;
        case fitsType<unsigned short>():
            m_valueStr << m_value.UShort;
            break;
        case fitsType<int>():
            m_valueStr << m_value.Int;
            break;
        case fitsType<unsigned int>():
            m_valueStr << m_value.UInt;
            break;
        case fitsType<long>():
            m_valueStr << m_value.Long;
            break;
        case fitsType<unsigned long>():
            m_valueStr << m_value.ULong;
            break;
        case fitsType<long long>():
            m_valueStr << m_value.LongLong;
            break;
        case fitsType<unsigned long long>():
            m_valueStr << m_value.ULongLong;
            break;
        case fitsType<float>():
            m_valueStr << m_value.Float;
            break;
        case fitsType<std::complex<float>>():
            m_valueStr << m_value.complexFloat;
            break;
        case fitsType<double>():
            m_valueStr << m_value.Double;
            break;
        case fitsType<std::complex<double>>():
            m_valueStr << m_value.complexDouble;
            break;
        case fitsType<fitsCommentType>():
            return mx::error_t::noerror;
        case fitsType<fitsHistoryType>():
            return mx::error_t::noerror;
        case fitsType<fitsContinueType>():
            return mx::error_t::noerror;
        default:
            return internal::mxlib_error_report<verboseT>( mx::error_t::invalidarg,
                                                           "Unknown FITS type for " + m_keyword );
    }

    m_valueStrGood = true;
    return mx::error_t::noerror;
}

template <class verboseT>
template <typename typeT>
error_t fitsHeaderCard<verboseT>::convertFromString()
{
    m_type = fitsType<typeT>();

    error_t errc;

    m_value.template member<typeT>() = ioutils::stoT<typeT>( m_valueStr.str(), &errc );

    if( errc != mx::error_t::noerror )
    {
        m_valueGood = false;
        return errc;
    }

    m_valueGood = true;

    return error_t::noerror;
}

template <class verboseT>
template <typename typeT>
error_t fitsHeaderCard<verboseT>::convertedValue( typeT &cval )
{
    switch( m_type )
    {
        case fitsType<unsigned char>():
        {
            cval = m_value.UChar;
            return error_t::noerror;
        }
        case fitsType<char>():
        {
            cval = m_value.Char;
            return error_t::noerror;
        }
        case fitsType<short>():
        {
            cval = m_value.Short;
            return error_t::noerror;
        }
        case fitsType<unsigned short>():
        {
            cval = m_value.UShort;
            return error_t::noerror;
        }
        case fitsType<int>():
        {
            cval = m_value.Int;
            return error_t::noerror;
        }
        case fitsType<unsigned int>():
        {
            cval = m_value.UInt;
            return error_t::noerror;
        }
        case fitsType<long>():
        {
            cval = m_value.Long;
            return error_t::noerror;
        }
        case fitsType<unsigned long>():
        {
            cval = m_value.ULong;
            return error_t::noerror;
        }
        case fitsType<long long>():
        {
            cval = m_value.LongLong;
            return error_t::noerror;
        }
        case fitsType<unsigned long long>():
        {
            cval = m_value.ULongLong;
            return error_t::noerror;
        }
        case fitsType<float>():
        {
            cval = m_value.Float;
            return error_t::noerror;
        }
        case fitsType<std::complex<float>>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::notimpl,
                                                           "can't convert complex type for " + m_keyword );
        }
        case fitsType<double>():
        {
            cval = m_value.Double;
            return error_t::noerror;
        }
        case fitsType<std::complex<double>>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::notimpl,
                                                           "can't convert complex type for " + m_keyword );
        }
        case fitsType<fitsCommentType>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "cannot convert comment to numeric type for " + m_keyword );
        }
        case fitsType<fitsHistoryType>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "cannot convert history to numeric type for " + m_keyword );
        }
        case fitsType<fitsContinueType>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "cannot convert continue to numeric type for " + m_keyword );
        }
        case TSTRING:
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "cannot convert string to numeric type for " + m_keyword );
        }
        default:
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "invalid FITS type conversion for " + m_keyword );
        }
    }
}

template <class verboseT>
error_t fitsHeaderCard<verboseT>::convertValue( int newtype )
{
    if( !m_valueGood )
    {
        m_type = newtype;
        return error_t::noerror;
    }

    mx::error_t errc;
    switch( newtype )
    {
        case fitsType<unsigned char>():
        {
            errc = convertedValue<unsigned char>( m_value.UChar );
            break;
        }
        case fitsType<char>():
        {
            errc = convertedValue<char>( m_value.Char );
            break;
        }
        case fitsType<short>():
        {
            errc = convertedValue<short>( m_value.Short );
            break;
        }
        case fitsType<unsigned short>():
        {
            errc = convertedValue<unsigned short>( m_value.UShort );
            break;
        }
        case fitsType<int>():
        {
            errc = convertedValue<int>( m_value.Int );
            break;
        }
        case fitsType<unsigned int>():
        {
            errc = convertedValue<unsigned int>( m_value.UInt );
            break;
        }
        case fitsType<long>():
        {
            errc = convertedValue<long>( m_value.Long );
            break;
        }
        case fitsType<unsigned long>():
        {
            errc = convertedValue<unsigned long>( m_value.ULong );
            break;
        }
        case fitsType<long long>():
        {
            errc = convertedValue<long long>( m_value.LongLong );
            break;
        }
        case fitsType<unsigned long long>():
        {
            errc = convertedValue<unsigned long long>( m_value.ULongLong );
            break;
        }
        case fitsType<float>():
        {
            errc = convertedValue<float>( m_value.Float );
            break;
        }
        case fitsType<std::complex<float>>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::notimpl,
                                                           "can't convert complex type for " + m_keyword );
        }
        case fitsType<double>():
        {
            errc = convertedValue<double>( m_value.Double );
            break;
        }
        case fitsType<std::complex<double>>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::notimpl,
                                                           "can't convert complex type for " + m_keyword );
        }
        case fitsType<fitsCommentType>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "cannot convert comment to numeric type for " + m_keyword );
        }
        case fitsType<fitsHistoryType>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "cannot convert history to numeric type for " + m_keyword );
        }
        case fitsType<fitsContinueType>():
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "cannot convert continue to numeric type for " + m_keyword );
        }
        case TSTRING:
        {
            convertToString();
            m_type = newtype;
            m_valueGood = false;
            return error_t::noerror;
        }
        default:
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg,
                                                           "invalid FITS type conversion for " + m_keyword );
        }
    }

    if( !!errc )
    {
        return errc;
    }

    m_type = newtype;
    m_valueGood = true;

    return error_t::noerror;
}

template <class verboseT>
const std::string &fitsHeaderCard<verboseT>::keyword() const
{
    return m_keyword;
}

template <class verboseT>
error_t fitsHeaderCard<verboseT>::keyword( const std::string &kw )
{
    m_keyword = kw;
    return error_t::noerror;
}

template <class verboseT>
int fitsHeaderCard<verboseT>::type() const
{
    return m_type;
}

template <class verboseT>
error_t fitsHeaderCard<verboseT>::type( const int &t )
{
    if( t == m_type )
    {
        return error_t::noerror;
    }

    if( m_valueGood )
    {
        error_t errc = convertValue( t );
        if( errc != error_t::noerror )
        {
            return errc;
        }
    }
    else
    {
        m_type = t;
    }

    // Need to reconvert, always favor the actual value.
    if( m_valueGood && m_valueStrGood )
    {
        m_valueStrGood = false;
    }

    return error_t::noerror;
}

template <class verboseT>
template <typename typeT>
typeT fitsHeaderCard<verboseT>::valueNonString( mx::error_t &errc )
{
    errc = mx::error_t::noerror; // to be changed if needed

    if( m_valueGood == false )
    {
        errc = convertFromString<typeT>();
    }

    if( m_type != fitsType<typeT>() )
    {
        typeT val;
        errc = convertedValue<typeT>( val );
        return val;
    }

    return m_value.template member<typeT>();
}

template <class verboseT>
std::string fitsHeaderCard<verboseT>::value( meta::tagT<std::string>, mx::error_t &errc )
{
    errc = mx::error_t::noerror;

    if( m_valueStrGood == false )
    {
        errc = convertToString();
        if( !!errc )
        {
            return "";
        }
    }

    // Strip ' from beginning and end if present
    std::string str = m_valueStr.str();

    // Reload it so it's there for next time:
    m_valueStr.str( str );

    return str;
}

template <class verboseT>
char fitsHeaderCard<verboseT>::value( meta::tagT<char>, mx::error_t &errc )
{
    return valueNonString<char>( errc );
}

template <class verboseT>
unsigned char fitsHeaderCard<verboseT>::value( meta::tagT<unsigned char>, mx::error_t &errc )
{
    return valueNonString<unsigned char>( errc );
}

template <class verboseT>
short fitsHeaderCard<verboseT>::value( meta::tagT<short>, mx::error_t &errc )
{
    return valueNonString<short>( errc );
}

template <class verboseT>
unsigned short fitsHeaderCard<verboseT>::value( meta::tagT<unsigned short>, mx::error_t &errc )
{
    return valueNonString<unsigned short>( errc );
}

template <class verboseT>
int fitsHeaderCard<verboseT>::value( meta::tagT<int>, mx::error_t &errc )
{
    return valueNonString<int>( errc );
}

template <class verboseT>
unsigned int fitsHeaderCard<verboseT>::value( meta::tagT<unsigned int>, mx::error_t &errc )
{
    return valueNonString<unsigned int>( errc );
}

template <class verboseT>
long fitsHeaderCard<verboseT>::value( meta::tagT<long>, mx::error_t &errc )
{
    return valueNonString<long>( errc );
}

template <class verboseT>
unsigned long fitsHeaderCard<verboseT>::value( meta::tagT<unsigned long>, mx::error_t &errc )
{
    return valueNonString<unsigned long>( errc );
}

template <class verboseT>
long long fitsHeaderCard<verboseT>::value( meta::tagT<long long>, mx::error_t &errc )
{
    return valueNonString<long long>( errc );
}

template <class verboseT>
unsigned long long fitsHeaderCard<verboseT>::value( meta::tagT<unsigned long long>, mx::error_t &errc )
{
    return valueNonString<unsigned long long>( errc );
}

template <class verboseT>
float fitsHeaderCard<verboseT>::value( meta::tagT<float>, mx::error_t &errc )
{
    return valueNonString<float>( errc );
}

template <class verboseT>
double fitsHeaderCard<verboseT>::value( meta::tagT<double>, mx::error_t &errc )
{
    return valueNonString<double>( errc );
}

template <class verboseT>
template <typename typeT>
typeT fitsHeaderCard<verboseT>::value( mx::error_t *errc )
{
    mx::error_t _errc;

    typeT val = value( meta::tagT<typeT>(), _errc );

    if( errc )
    {
        *errc = _errc;
    }

    if( _errc != mx::error_t::noerror )
    {
        internal::mxlib_error_report<verboseT>( _errc, "getting value for " + m_keyword );
    }

    return val;
}

template <class verboseT>
std::string fitsHeaderCard<verboseT>::String( error_t *errc )
{
    return value<std::string>( errc );
}

template <class verboseT>
char fitsHeaderCard<verboseT>::Char( error_t *errc )
{
    return value<char>( errc );
}

template <class verboseT>
unsigned char fitsHeaderCard<verboseT>::UChar( error_t *errc )
{
    return value<unsigned char>( errc );
}

template <class verboseT>
short fitsHeaderCard<verboseT>::Short( error_t *errc )
{
    return value<short>( errc );
}

template <class verboseT>
unsigned short fitsHeaderCard<verboseT>::UShort( error_t *errc )
{
    return value<unsigned short>( errc );
}

template <class verboseT>
int fitsHeaderCard<verboseT>::Int( error_t *errc )
{
    return value<int>( errc );
}

template <class verboseT>
unsigned int fitsHeaderCard<verboseT>::UInt( error_t *errc )
{
    return value<unsigned int>( errc );
}

template <class verboseT>
long fitsHeaderCard<verboseT>::Long( error_t *errc )
{
    return value<long>( errc );
}

template <class verboseT>
unsigned long fitsHeaderCard<verboseT>::ULong( error_t *errc )
{
    return value<unsigned long>( errc );
}

template <class verboseT>
long long fitsHeaderCard<verboseT>::LongLong( error_t *errc )
{
    return value<long long>( errc );
}

template <class verboseT>
unsigned long long fitsHeaderCard<verboseT>::ULongLong( error_t *errc )
{
    return value<unsigned long long>( errc );
}

template <class verboseT>
float fitsHeaderCard<verboseT>::Float( error_t *errc )
{
    return value<float>( errc );
}

/*template <class verboseT>
std::complex<float> fitsHeaderCard<verboseT>::complexFloat()
{
    return value<std::complex<float>>();
}*/

template <class verboseT>
double fitsHeaderCard<verboseT>::Double( error_t *errc )
{
    return value<double>( errc );
}

/*
template <class verboseT>
std::complex<double> fitsHeaderCard<verboseT>::complexDouble()
{
    return value<std::complex<double>>();
}*/

template <class verboseT>
mx::error_t fitsHeaderCard<verboseT>::value( const char *v )
{
    std::string str = v;

    return value( str );
}

template <class verboseT>
mx::error_t fitsHeaderCard<verboseT>::value( const std::string &v )
{
    // Strip ' from beginning and end if present
    std::string str = v;

    stripApostWS( str );

    // Reload it so it's there for next time:
    m_valueStr.str( str );
    m_valueGood = false;
    m_valueStrGood = true;
    m_type = fitsType<std::string>();

    return mx::error_t::noerror;
}

template <class verboseT>
template <typename typeT>
mx::error_t fitsHeaderCard<verboseT>::value( const typeT &v )
{
    m_type = fitsType<typeT>();
    m_value.template member<typeT>() = v;
    m_valueGood = true;
    m_valueStrGood = false;

    return mx::error_t::noerror;
}

template <class verboseT>
std::string fitsHeaderCard<verboseT>::valueStr()
{
    if( !m_valueStrGood )
    {
        convertToString();
    }

    std::string s = m_valueStr.str();
    return s;
}

template <class verboseT>
bool fitsHeaderCard<verboseT>::valueGood()
{
    return m_valueGood;
}

template <class verboseT>
bool fitsHeaderCard<verboseT>::valueStrGood()
{
    return m_valueStrGood;
}

template <class verboseT>
const std::string &fitsHeaderCard<verboseT>::comment()
{
    return m_comment;
}

template <class verboseT>
error_t fitsHeaderCard<verboseT>::comment( const std::string &c )
{
    m_comment = c;
    return error_t::noerror;
}

template <class verboseT>
error_t fitsHeaderCard<verboseT>::appendContinue( const fitsHeaderCard<verboseT> &card )
{
    // Check if m_type is string
    if( m_type != fitsType<char *>() && m_type != fitsType<std::string>() )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "attempt to continue a non-string card" );
    }

    // Check if m_valueStrGood is true
    if( !m_valueStrGood )
    {
        return internal::mxlib_error_report<verboseT>( error_t::invalidconfig,
                                                       "attempt to continue a card with no value" );
    }

    std::string newstr = card.m_comment;

    // Check if this is the last one
    size_t slash = newstr.find_last_of( '\'' );
    if( slash != std::string::npos )
    {
        slash = newstr.find_first_of( '/', slash );

        if( slash != std::string::npos )
        {
            // extract comment
            size_t comm = newstr.find_first_not_of( ' ', slash + 1 );
            if( comm != std::string::npos )
            {
                m_comment = newstr.substr( comm );
            }

            // and then erase it
            newstr.erase( slash );
        }
    }
    stripApostWS( newstr );

    // have to remove & if needed
    std::string vstr = m_valueStr.str();
    if( vstr.back() == '&' )
    {
        vstr.erase( vstr.size() - 1 );
    }

    m_valueStr.str( vstr + newstr );

    return error_t::noerror;
}

template <class verboseT>
mx::error_t fitsHeaderCard<verboseT>::write( fitsfile *fptr )
{
    if( m_type == fitsType<char *>() || m_type == fitsType<std::string>() )
    {
        return fits_write_key<char *>( fptr,
                                       (char *)m_keyword.c_str(),
                                       (void *)m_valueStr.str().c_str(),
                                       (char *)m_comment.c_str() );
    }

    // If the string is good, meaning already converted.
    if( m_valueStrGood == true )
    {
        // This populates the card directly.
        return fits_write_key<fitsUnknownType>( fptr,
                                                (char *)m_keyword.c_str(),
                                                (void *)m_valueStr.str().c_str(),
                                                (char *)m_comment.c_str() );
    }

    // Ok, now we write the type directly using fitsio routines because it hasn't been converted.
    switch( m_type )
    {
        case fitsType<bool>():
        {
            return fits_write_key<bool>( fptr, (char *)m_keyword.c_str(), &m_value.Bool, (char *)m_comment.c_str() );
        }
        case fitsType<char>():
        {
            return fits_write_key<char>( fptr, (char *)m_keyword.c_str(), &m_value.Char, (char *)m_comment.c_str() );
        }
        case fitsType<unsigned char>():
        {
            return fits_write_key<unsigned char>( fptr,
                                                  (char *)m_keyword.c_str(),
                                                  &m_value.UChar,
                                                  (char *)m_comment.c_str() );
        }
        case fitsType<short>():
        {
            return fits_write_key<short>( fptr, (char *)m_keyword.c_str(), &m_value.Short, (char *)m_comment.c_str() );
        }
        case fitsType<unsigned short>():
        {
            return fits_write_key<unsigned short>( fptr,
                                                   (char *)m_keyword.c_str(),
                                                   &m_value.UShort,
                                                   (char *)m_comment.c_str() );
        }
        case fitsType<int>():
        {
            return fits_write_key<int>( fptr, (char *)m_keyword.c_str(), &m_value.Int, (char *)m_comment.c_str() );
        }
        case fitsType<unsigned int>():
        {
            return fits_write_key<unsigned int>( fptr,
                                                 (char *)m_keyword.c_str(),
                                                 &m_value.UInt,
                                                 (char *)m_comment.c_str() );
        }
        case fitsType<long>():
        {
            return fits_write_key<long>( fptr, (char *)m_keyword.c_str(), &m_value.Long, (char *)m_comment.c_str() );
        }
        case fitsType<unsigned long>():
        {
            return fits_write_key<unsigned long>( fptr,
                                                  (char *)m_keyword.c_str(),
                                                  &m_value.ULong,
                                                  (char *)m_comment.c_str() );
        }
        case fitsType<long long>():
        {
            return fits_write_key<long long>( fptr,
                                              (char *)m_keyword.c_str(),
                                              &m_value.LongLong,
                                              (char *)m_comment.c_str() );
        }
        case fitsType<unsigned long long>():
        {
            return fits_write_key<unsigned long long>( fptr,
                                                       (char *)m_keyword.c_str(),
                                                       &m_value.ULongLong,
                                                       (char *)m_comment.c_str() );
        }
        case fitsType<float>():
        {
            return fits_write_key<float>( fptr, (char *)m_keyword.c_str(), &m_value.Float, (char *)m_comment.c_str() );
        }
        case fitsType<double>():
        {
            return fits_write_key<double>( fptr,
                                           (char *)m_keyword.c_str(),
                                           &m_value.Double,
                                           (char *)m_comment.c_str() );
        }
        case fitsType<fitsCommentType>():
        {
            return fits_write_comment( fptr, (char *)m_comment.c_str() );
        }
        case fitsType<fitsHistoryType>():
        {
            return fits_write_history( fptr, (char *)m_comment.c_str() );
        }
        default:
        {
            return internal::mxlib_error_report<verboseT>( error_t::invalidarg, "invalid FITS type for " + m_keyword );
        }
    }
}

extern template class fitsHeaderCard<verbose::d>;

} // namespace fits
} // namespace mx

#endif // ioutils_fits_fitsHeaderCard_hpp
