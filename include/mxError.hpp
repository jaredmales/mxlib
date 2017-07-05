/** \file mxError.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Declares and defines the mxlib error reporting system.
  * \ingroup error_handling
  * 
*/

#ifndef __mxError__
#define __mxError__

#include <cerrno>
#include <cstring>
#include <sstream>

/** \defgroup error_macros Error Handling Macros
  * \brief Macros controlling how mxlib reports errors
  * 
  * mxlib components use these macros to report errors to the user, which by default are wrappers for mx::error_report and mx::errno_report.
  * This behavior can be changed with preprocessor defines.
  * 
  * To completely suppress all mxlib error reporting, include the following before includng any mxlib headers
  * \code
  * #define MX_NO_ERROR_REPORTS
  * \endcode
  * 
  * \ingroup error_handling
  * @{
  */

/** \def MX_NO_ERROR_REPORTS
  * \brief If this is defined, then no errors are reported.
  * \addtogroup error_macros
  */
#ifdef MX_NO_ERROR_REPORTS


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
  * \addtogroup error_macros
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
#define MXE_INVALIDARG_MSG "An argument was invalid."

/** \def MXE_NOTIMPL
  * \brief A component or technique is not implemented.
  * \ingroup mxe_errors
  */   
#define MXE_NOTIMPL 30
#define MXE_NOTIMPL_MSG "A component or technique is not implemented."

/** \def MXE_PARAMNOTSET
  * \brief A parameter was not set
  * \ingroup mxe_errors
  */   
#define MXE_PARAMNOTSET 35
#define MXE_PARAMNOTSET_MSG "A parameter was not set."

/** \def MXE_NOTFOUND
  * \brief An item was not found
  * \ingroup mxe_errors
  */   
#define MXE_NOTFOUND 40
#define MXE_NOTFOUND_MSG "An item was not found."


/** \def MXE_SIZEERR
  * \brief A size was calculated incorrectly
  * \ingroup mxe_errors
  */   
#define MXE_SIZEERR 55
#define MXE_SIZEERR_MSG "A size was calculated incorrectly."


/** \def MXE_FILEOERR
  * \brief An error occurred while opening a file.
  * \ingroup mxe_errors
  */   
#define MXE_FILEOERR 1034
#define MXE_FILEOERR_MSG "An error occurred while opening a file."

/** \def MXE_FILEWERR
  * \brief An error occurred while writing to a file.
  * \ingroup mxe_errors
  */   
#define MXE_FILEWERR 1044
#define MXE_FILEWERR_MSG "An error occurred while writing to a file."

/** \def MXE_FILERERR
  * \brief An error occurred while reading from a file.
  * \ingroup mxe_errors
  */   
#define MXE_FILERERR 1049
#define MXE_FILERERR_MSG "An error occurred while reading from a file."

/** \def MXE_FILECERR
  * \brief An error occurred while closing a file.
  * \ingroup mxe_errors
  */   
#define MXE_FILECERR 1054
#define MXE_FILECERR_MSG "An error occurred while closing a file."

/** \def MXE_FILENOTFOUND
  * \brief The file was not found.
  * \ingroup mxe_errors
  */   
#define MXE_FILENOTFOUND 1059
#define MXE_FILENOTFOUND_MSG "The file was not found."

/** \def MXE_PROCERR
  * \brief An error occrred while starting a process.
  * \ingroup mxe_errors
  */   
#define MXE_PROCERR 2001
#define MXE_PROCERR_MSG "An error occured while starting a process."

/** \def MXE_TIMEOUT
  * \brief A timeout occurred.
  * \ingroup mxe_errors
  */   
#define MXE_TIMEOUT 2322
#define MXE_TIMEOUT_MSG "A timeout occurred."


/** \def MXE_GNUPLOTERR
  * \brief An error was returned by gnuplot.
  * \ingroup mxe_errors
  */   
#define MXE_GNUPLOTERR 4567
#define MXE_GNUPLOTERR_MSG "An error was returned by gnuplot."
   
///Return the macro name and a message for an mxlib error code
/**
  * \param [in] ec is the error code
  * \param [out] message is the human friendly message
  * 
  * \returns the name of the macro corresponding to the error code.
  * 
  * \ingroup error_handling
  */
std::string MXE_CodeToName( int ec, std::string & message)
{
   
   switch(ec)
   {
      case MXE_INVALIDARG:
         message = MXE_INVALIDARG_MSG;
         return "MXE_INVALIDARG";
      case MXE_NOTIMPL:
         message = MXE_NOTIMPL;
         return "MXE_NOTIMPL";
      case MXE_PARAMNOTSET:
         message = MXE_PARAMNOTSET_MSG;
         return "MXE_PARAMNOTSET";
      case MXE_NOTFOUND:
         message = MXE_NOTFOUND_MSG;
         return "MXE_NOTFOUND";
      case MXE_SIZEERR:
         message = MXE_SIZEERR_MSG;
         return "MXE_SIZEERR";
      case MXE_FILEOERR:
         message = MXE_FILEOERR_MSG;
         return "MXE_FILEOERR";
      case MXE_FILEWERR:
         message = MXE_FILEWERR_MSG;
         return "MXE_FILEWERR";
      case MXE_FILERERR:
         message = MXE_FILERERR_MSG;
         return "MXE_FILERERR";
      case MXE_FILECERR:
         message = MXE_FILECERR_MSG;
         return "MXE_FILECERR";
      case MXE_FILENOTFOUND:
         message = MXE_FILENOTFOUND_MSG;
         return "MXE_FILENOTFOUND";
      case MXE_PROCERR:
         message = MXE_PROCERR_MSG;
         return "MXE_PROCERR";
      case MXE_TIMEOUT:
         message = MXE_TIMEOUT_MSG;
         return "MXE_TIMEOUT";
      case MXE_GNUPLOTERR:
         message = MXE_GNUPLOTERR_MSG;
         return "MXE_GNUPLOTERR";
      default:
         message = "Unknown mxlib error code.";
         return "?";
   }
}

///Return the macro name and a message for a standard errno code
/**
  * \param[in] ec is the error code
  * \param[out] message is the human friendly message retrieved using std::strerror
  * 
  * \returns the name of the macro corresponding to the code.
  * 
  * \ingroup error_handling
  */
std::string errno_CodeToName( int ec, std::string & message)
{
   message = std::strerror(ec);
   
   switch(ec)
   {
      #ifdef E2BIG
      case E2BIG:
          return "E2BIG";
      #endif 
      #ifdef EACCES
      case EACCES:
          return "EACCES";
      #endif 
      #ifdef EADDRINUSE
      case EADDRINUSE:
          return "EADDRINUSE";
      #endif 
      #ifdef EADDRNOTAVAIL
      case EADDRNOTAVAIL:
          return "EADDRNOTAVAIL";
      #endif 
      #ifdef EAFNOSUPPORT
      case EAFNOSUPPORT:
          return "EAFNOSUPPORT";
      #endif 
      #ifdef EAGAIN
      case EAGAIN:
          #if (EWOULDBLOCK == EAGAIN)
             return "EAGIAN / EWOULDBLOCK";
          #else
             return "EAGAIN";
          #endif
      #endif 
      #ifdef EALREADY
      case EALREADY:
          return "EALREADY";
      #endif 
      #ifdef EBADF
      case EBADF:
          return "EBADF";
      #endif 
      #ifdef EBADMSG
      case EBADMSG:
          return "EBADMSG";
      #endif 
      #ifdef EBUSY
      case EBUSY:
          return "EBUSY";
      #endif 
      #ifdef ECANCELED
      case ECANCELED:
          return "ECANCELED";
      #endif 
      #ifdef ECHILD
      case ECHILD:
          return "ECHILD";
      #endif 
      #ifdef ECONNABORTED
      case ECONNABORTED:
          return "ECONNABORTED";
      #endif 
      #ifdef ECONNREFUSED
      case ECONNREFUSED:
          return "ECONNREFUSED";
      #endif 
      #ifdef ECONNRESET
      case ECONNRESET:
          return "ECONNRESET";
      #endif 
      #ifdef EDESTADDRREQ
      case EDESTADDRREQ:
          return "EDESTADDRREQ";
      #endif 
      #ifdef EDOM
      case EDOM:
          return "EDOM";
      #endif 
      #ifdef EEXIST
      case EEXIST:
          return "EEXIST";
      #endif 
      #ifdef EFAULT
      case EFAULT:
          return "EFAULT";
      #endif 
      #ifdef EFBIG
      case EFBIG:
          return "EFBIG";
      #endif 
      #ifdef EHOSTUNREACH
      case EHOSTUNREACH:
          return "EHOSTUNREACH";
      #endif 
      #ifdef EIDRM
      case EIDRM:
          return "EIDRM";
      #endif 
      #ifdef EILSEQ
      case EILSEQ:
          return "EILSEQ";
      #endif 
      #ifdef EINPROGRESS
      case EINPROGRESS:
          return "EINPROGRESS";
      #endif 
      #ifdef EINTR
      case EINTR:
          return "EINTR";
      #endif 
      #ifdef EINVAL
      case EINVAL:
          return "EINVAL";
      #endif 
      #ifdef EIO
      case EIO:
          return "EIO";
      #endif 
      #ifdef EISCONN
      case EISCONN:
          return "EISCONN";
      #endif 
      #ifdef EISDIR
      case EISDIR:
          return "EISDIR";
      #endif 
      #ifdef ELOOP
      case ELOOP:
          return "ELOOP";
      #endif 
      #ifdef EMFILE
      case EMFILE:
          return "EMFILE";
      #endif 
      #ifdef EMLINK
      case EMLINK:
          return "EMLINK";
      #endif 
      #ifdef EMSGSIZE
      case EMSGSIZE:
          return "EMSGSIZE";
      #endif 
      #ifdef ENAMETOOLONG
      case ENAMETOOLONG:
          return "ENAMETOOLONG";
      #endif 
      #ifdef ENETDOWN
      case ENETDOWN:
          return "ENETDOWN";
      #endif 
      #ifdef ENETRESET
      case ENETRESET:
          return "ENETRESET";
      #endif 
      #ifdef ENETUNREACH
      case ENETUNREACH:
          return "ENETUNREACH";
      #endif 
      #ifdef ENFILE
      case ENFILE:
          return "ENFILE";
      #endif 
      #ifdef ENOBUFS
      case ENOBUFS:
          return "ENOBUFS";
      #endif 
      #ifdef ENODATA
      case ENODATA:
          return "ENODATA";
      #endif 
      #ifdef ENODEV
      case ENODEV:
          return "ENODEV";
      #endif 
      #ifdef ENOENT
      case ENOENT:
          return "ENOENT";
      #endif 
      #ifdef ENOEXEC
      case ENOEXEC:
          return "ENOEXEC";
      #endif 
      #ifdef ENOLCK
      case ENOLCK:
          return "ENOLCK";
      #endif 
      #ifdef ENOLINK
      case ENOLINK:
          return "ENOLINK";
      #endif 
      #ifdef ENOMEM
      case ENOMEM:
          return "ENOMEM";
      #endif 
      #ifdef ENOMSG
      case ENOMSG:
          return "ENOMSG";
      #endif 
      #ifdef ENOPROTOOPT
      case ENOPROTOOPT:
          return "ENOPROTOOPT";
      #endif 
      #ifdef ENOSPC
      case ENOSPC:
          return "ENOSPC";
      #endif 
      #ifdef ENOSR
      case ENOSR:
          return "ENOSR";
      #endif 
      #ifdef ENOSTR
      case ENOSTR:
          return "ENOSTR";
      #endif 
      #ifdef ENOSYS
      case ENOSYS:
          return "ENOSYS";
      #endif 
      #ifdef ENOTCONN
      case ENOTCONN:
          return "ENOTCONN";
      #endif 
      #ifdef ENOTDIR
      case ENOTDIR:
          return "ENOTDIR";
      #endif 
      #ifdef ENOTEMPTY
      case ENOTEMPTY:
          return "ENOTEMPTY";
      #endif 
      #ifdef ENOTRECOVERABLE
      case ENOTRECOVERABLE:
          return "ENOTRECOVERABLE";
      #endif 
      #ifdef ENOTSOCK
      case ENOTSOCK:
          return "ENOTSOCK";
      #endif 
      #ifdef ENOTSUP
      #if (ENOTSUP != EOPNOTSUPP)
      case ENOTSUP:
          return "ENOTSUP";
      #endif
      #endif 
      #ifdef ENOTTY
      case ENOTTY:
          return "ENOTTY";
      #endif 
      #ifdef ENXIO
      case ENXIO:
          return "ENXIO";
      #endif 
      #ifdef EOPNOTSUPP
      case EOPNOTSUPP:
          return "EOPNOTSUPP";
      #endif 
      #ifdef EOVERFLOW
      case EOVERFLOW:
          return "EOVERFLOW";
      #endif 
      #ifdef EOWNERDEAD
      case EOWNERDEAD:
          return "EOWNERDEAD";
      #endif 
      #ifdef EPERM
      case EPERM:
          return "EPERM";
      #endif 
      #ifdef EPIPE
      case EPIPE:
          return "EPIPE";
      #endif 
      #ifdef EPROTO
      case EPROTO:
          return "EPROTO";
      #endif 
      #ifdef EPROTONOSUPPORT
      case EPROTONOSUPPORT:
          return "EPROTONOSUPPORT";
      #endif 
      #ifdef EPROTOTYPE
      case EPROTOTYPE:
          return "EPROTOTYPE";
      #endif 
      #ifdef ERANGE
      case ERANGE:
          return "ERANGE";
      #endif 
      #ifdef EROFS
      case EROFS:
          return "EROFS";
      #endif 
      #ifdef ESPIPE
      case ESPIPE:
          return "ESPIPE";
      #endif 
      #ifdef ESRCH
      case ESRCH:
          return "ESRCH";
      #endif 
      #ifdef ETIME
      case ETIME:
          return "ETIME";
      #endif 
      #ifdef ETIMEDOUT
      case ETIMEDOUT:
          return "ETIMEDOUT";
      #endif 
      #ifdef ETXTBSY
      case ETXTBSY:
          return "ETXTBSY";
      #endif 
      #ifdef EWOULDBLOCK
      #if (EWOULDBLOCK != EAGAIN)
      case EWOULDBLOCK:
          return "EWOULDBLOCK";
      #endif
      #endif 
      #ifdef EXDEV
      case EXDEV:
          return "EXDEV";
      #endif 


      default:
         message = "Unknown errno code.";
         return "?";
   }
}


///Construct a rich error report given an mxlib error code
/**
  * \param[in] source is intended to identify the mxlib component (i.e. the class name)
  * \param[in] code is an MXE_* error code
  * \param[in] file should be passed the __FILE__ macro
  * \param[in] line should be passed the __LINE__ macro
  * \param[in] expl [optional] if more information can be provided, use this to inform the user.
  * 
  * \return the formatted error report.
  * 
  * \ingroup error_handling
  */ 
std::string error_report(const std::string & source, 
                           const int & code, 
                           const std::string & file, 
                           const int & line,
                           const std::string & expl = "")
{
   
   std::string codeName, codeMessage;
   
   codeName = MXE_CodeToName(code, codeMessage);
   
   std::ostringstream s;
   s.str("");
   
   s << "An error has occured in an mxlib component.\n";
   s << "      source: " << source << "\n";
   s << "        code: " << codeName << "(" << code << ")\n";
   s << "    code msg: " << codeMessage << "\n";
   s << "     in file: " << file << "\n";
   s << "     at line: " << line << "\n";
   if(expl != "") 
   s << " explanation: " << expl << "\n";
   
   return s.str();
}

///Construct a rich error report given a standard errno error code
/**
  * \param[in] source is intended to identify the component (i.e. the class name)
  * \param[in] code is an errno erro code
  * \param[in] file should be passed the __FILE__ macro
  * \param[in] line should be passed the __LINE__ macro
  * \param[in] expl [optional] if more information can be provided, use this to inform the user.
  * 
  * \return the formatted error report.
  * 
  * \ingroup error_handling
  */ 
std::string errno_report( const std::string & source,
                          int ec,
                          const std::string & file, 
                          const int & line,
                          const std::string & expl = "" )
{
   std::string codeName, codeMessage;
   
   codeName = errno_CodeToName(ec, codeMessage);
   
   std::ostringstream s;
   s.str("");
   
   s << "An  error has occured in an mxlib component.\n";
   s << "      source: " << source << "\n";
   s << "  errno code: " << codeName << "(" << ec << ")\n";
   s << "    code msg: " << codeMessage << "\n";
   s << "     in file: " << file << "\n";
   s << "     at line: " << line << "\n";
   if(expl != "") 
   s << " explanation: " << expl << "\n";
   
   return s.str();
}






} //namespace mx

#endif //__mxError__
