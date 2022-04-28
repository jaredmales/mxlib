/** \file mxError.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines the mxlib error reporting system.
  * \ingroup error_handling_files
  * 
*/

//***********************************************************************//
// Copyright 2015, 2016, 2017 Jared R. Males (jaredmales@gmail.com)
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

#ifndef __mxError__
#define __mxError__

#include <cerrno>
#include <cstring>
#include <sstream>

/** \addtogroup error_macros 
  * 
  * mxlib components use these macros to report errors to the user, which by default are wrappers for mx::error_report and mx::errno_report.
  * This behavior can be changed with preprocessor defines.
  * 
  * To completely suppress all mxlib error reporting, include the following before includng any mxlib headers
  * \code
  * #define MX_NO_ERROR_REPORTS
  * \endcode
  */

  /** \ingroup error_handling
  * @{
  */

/** \def MX_NO_ERROR_REPORTS
  * \brief If this is defined, then no errors are reported.
  * \ingroup error_macros
  */
#ifdef MX_NO_ERROR_REPORTS

//Just defining these as empty macros, though this is currently untested.
/** \todo test MX_NO_ERROR_REPORTS */

#define mxError(esrc,ecode,expl) 
#define mxPError(esrc,ecode,expl)

#else


#ifndef mxError

#include <iostream>

/** \def mxError
  * \brief This reports an mxlib specific error
  * 
  * Can be predefined to change the error reporting of mxlib.
  *
  * \param[in] esrc is intended to identify the component (i.e. the class name)
  * \param[in] ecode is an errno erro code
  * \param[in] expl [optional] if more information can be provided, use this to inform the user.
  * 
  * \ingroup error_macros
 */
#define mxError(esrc,ecode,expl) std::cerr << "\n" << mx::error_report(esrc, ecode, __FILE__, __LINE__, expl) << "\n";

#endif //mxError

#ifndef mxPError

/** \def mxPError
  * \brief This reports a standard library error, taking after perror.
  * 
  * Can be pre-defined to change the error reporting of mxlib.
  * 
  * \param[in] esrc is intended to identify the component (i.e. the class name)
  * \param[in] ecode is an errno erro code
  * \param[in] expl [optional] if more information can be provided, use this to inform the user.
  * 
  * \addtogroup error_macros
  */
#define mxPError(esrc,ecode,expl) std::cerr << "\n" << mx::errno_report(esrc, ecode, __FILE__, __LINE__, expl) << "\n";


#endif //mxPError

#endif //MX_NO_ERROR_REPORTS

/// @}


namespace mx
{
/** \defgroup mxe_errors mxlib Error Codes
  * \ingroup error_handling
  */  

/** \def MXE_INVALIDARG
  * \brief An argument was invalid.
  * \ingroup mxe_errors
  */   
#define MXE_INVALIDARG 25
#define MXE_INVALIDARG_NAME "MXE_INVALIDARG"
#define MXE_INVALIDARG_MSG "An argument was invalid."

/** \def MXE_INVALIDCONFIG
  * \brief A config setting was invalid.
  * \ingroup mxe_errors
  */   
#define MXE_INVALIDCONFIG 27
#define MXE_INVALIDCONFIG_NAME "MXE_INVALIDCONFIG"
#define MXE_INVALIDCONFIG_MSG "A config setting was invalid."

/** \def MXE_NOTIMPL
  * \brief A component or technique is not implemented.
  * \ingroup mxe_errors
  */   
#define MXE_NOTIMPL 30
#define MXE_NOTIMPL_NAME "MXE_NOTIMPL"
#define MXE_NOTIMPL_MSG "A component or technique is not implemented."

/** \def MXE_PARAMNOTSET
  * \brief A parameter was not set
  * \ingroup mxe_errors
  */   
#define MXE_PARAMNOTSET 35
#define MXE_PARAMNOTSET_NAME "MXE_PARAMNOTSET"
#define MXE_PARAMNOTSET_MSG "A parameter was not set."

/** \def MXE_ENVNOTSET
  * \brief An environment variable is not set
  * \ingroup mxe_errors
  */   
#define MXE_ENVNOTSET 36
#define MXE_ENVNOTSET_NAME "MXE_ENVNOTSET"
#define MXE_ENVNOTSET_MSG "An environment variable is not set."

/** \def MXE_NOTFOUND
  * \brief An item was not found
  * \ingroup mxe_errors
  */   
#define MXE_NOTFOUND 40
#define MXE_NOTFOUND_NAME "MXE_NOTFOUND"
#define MXE_NOTFOUND_MSG "An item was not found."

/** \def MXE_SIZEERR
  * \brief A size was invalid or calculated incorrectly
  * \ingroup mxe_errors
  */   
#define MXE_SIZEERR 55
#define MXE_SIZEERR_NAME "MXE_SIZEERR"
#define MXE_SIZEERR_MSG "A size was invalid or calculated incorrectly."

/** \def MXE_ALLOCERR
  * \brief An error occurred during memory allocation.
  * \ingroup mxe_errors
  */   
#define MXE_ALLOCERR 60
#define MXE_ALLOCERR_NAME "MXE_ALLOCERR"
#define MXE_ALLOCERR_MSG "An error occurred during memory allocation."

/** \def MXE_FREEERR
  * \brief An error occurred during memory de-allocation.
  * \ingroup mxe_errors
  */   
#define MXE_FREEERR 65
#define MXE_FREEERR_NAME "MXE_FREEERR"
#define MXE_FREEERR_MSG "An error occurred during memory de-allocation."


/** \def MXE_PARSEERR
  * \brief A parsing error occurred.
  * \ingroup mxe_errors
  */   
#define MXE_PARSEERR 75
#define MXE_PARSEERR_NAME "MXE_PARSEERR"
#define MXE_PARSEERR_MSG "A parsing error occurred."


/** \def MXE_FILEOERR
  * \brief An error occurred while opening a file.
  * \ingroup mxe_errors
  */   
#define MXE_FILEOERR 1034
#define MXE_FILEOERR_NAME "MXE_FILEOERR"
#define MXE_FILEOERR_MSG "An error occurred while opening a file."

/** \def MXE_FILEWERR
  * \brief An error occurred while writing to a file.
  * \ingroup mxe_errors
  */   
#define MXE_FILEWERR 1044
#define MXE_FILEWERR_NAME "MXE_FILEWERR"
#define MXE_FILEWERR_MSG "An error occurred while writing to a file."

/** \def MXE_FILERERR
  * \brief An error occurred while reading from a file.
  * \ingroup mxe_errors
  */   
#define MXE_FILERERR 1049
#define MXE_FILERERR_NAME "MXE_FILERERR"
#define MXE_FILERERR_MSG "An error occurred while reading from a file."

/** \def MXE_FILECERR
  * \brief An error occurred while closing a file.
  * \ingroup mxe_errors
  */   
#define MXE_FILECERR 1054
#define MXE_FILECERR_NAME "MXE_FILECERR"
#define MXE_FILECERR_MSG "An error occurred while closing a file."

/** \def MXE_FILENOTFOUND
  * \brief The file was not found.
  * \ingroup mxe_errors
  */   
#define MXE_FILENOTFOUND 1059
#define MXE_FILENOTFOUND_NAME "MXE_FILENOTFOUND"
#define MXE_FILENOTFOUND_MSG "The file was not found."

/** \def MXE_PROCERR
  * \brief An error occrred while starting a process.
  * \ingroup mxe_errors
  */   
#define MXE_PROCERR 2001
#define MXE_PROCERR_NAME "MXE_PROCERR"
#define MXE_PROCERR_MSG "An error occured while starting a process."

/** \def MXE_TIMEOUT
  * \brief A timeout occurred.
  * \ingroup mxe_errors
  */   
#define MXE_TIMEOUT 2322
#define MXE_TIMEOUT_NAME "MXE_TIMEOUT"
#define MXE_TIMEOUT_MSG "A timeout occurred."

/** \def MXE_LIBERR
  * \brief An error was returned by a library.
  * \ingroup mxe_errors
  */   
#define MXE_LIBERR 4000
#define MXE_LIBERR_NAME "MXE_LIBERR"
#define MXE_LIBERR_MSG "An error was returned by a library."


/** \def MXE_GNUPLOTERR
  * \brief An error was returned by gnuplot.
  * \ingroup mxe_errors
  */   
#define MXE_GNUPLOTERR 4567
#define MXE_GNUPLOTERR_NAME "MXE_GNUPLOTERR"
#define MXE_GNUPLOTERR_MSG "An error was returned by gnuplot."

/** \def MXE_LAPACKERR
  * \brief An error was returned by Lapack.
  * \ingroup mxe_errors
  */   
#define MXE_LAPACKERR 6890
#define MXE_LAPACKERR_NAME "MXE_LAPACKERR"
#define MXE_LAPACKERR_MSG "An error was returned by Lapack."


///Return the name for an mxlib error code
/**
  * 
  * \returns the name of the macro corresponding to the error code.
  * 
  * \ingroup error_handling
  */
std::string MXE_CodeToName( int ec /**< [in] the error code */ );

///Return the description for an mxlib error code
/**
  * 
  * \returns the description for and error code.
  * 
  * \ingroup error_handling
  */
std::string MXE_CodeToDescription( int ec /**< [in] the error code */ );

///Return the macro name and a message for a standard errno code
/**
  * 
  * \returns the name of the macro corresponding to the code.
  * 
  * \ingroup error_handling
  */
std::string errno_CodeToName( int ec /**< [in] the error code */ );

///Construct a rich error report given an mxlib error code
/**
  * 
  * \return the formatted error report.
  * 
  * \ingroup error_handling
  */ 
std::string error_report( const std::string & source,   ///< [in] is intended to identify the mxlib component (i.e. the class name)
                          const int & code,             ///< [in] is an MXE_* error code
                          const std::string & file,     ///< [in] should be passed the __FILE__ macro
                          const int & line,             ///< [in] should be passed the __LINE__ macro
                          const std::string & expl = "" ///< [in] [optional] if more information can be provided, use this to inform the user.
                        );


///Construct a rich error report given a standard errno error code
/**
  * \return the formatted error report.
  * 
  * \ingroup error_handling
  */ 
std::string errno_report( const std::string & source,   ///< [in] intended to identify the component (i.e. the class name)
                          int ec,                       ///< [in] an errno erro code
                          const std::string & file,     ///< [in] file should be passed the __FILE__ macro
                          const int & line,             ///< [in] line should be passed the __LINE__ macro
                          const std::string & expl = "" ///< [in] [optional] if more information can be provided, use this to inform the user.
                        );






} //namespace mx

#endif //__mxError__
