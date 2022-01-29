/** \file mxError.cpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Function definitions for the mxlib error reporting system.
  * \ingroup error_handling_files
  * 
*/

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#include "mxError.hpp"

namespace mx
{

std::string MXE_CodeToName( int ec )
{
   
   switch(ec)
   {
      case MXE_INVALIDARG:
         return MXE_INVALIDARG_NAME;
      case MXE_INVALIDCONFIG:
         return MXE_INVALIDCONFIG_NAME;
      case MXE_NOTIMPL:
         return MXE_NOTIMPL_NAME;
      case MXE_PARAMNOTSET:
         return MXE_PARAMNOTSET_NAME;
      case MXE_ENVNOTSET:
         return MXE_ENVNOTSET_NAME;
      case MXE_NOTFOUND:
         return MXE_NOTFOUND_NAME;
      case MXE_SIZEERR:
         return MXE_SIZEERR_NAME;
      case MXE_ALLOCERR:
         return MXE_ALLOCERR_NAME;
      case MXE_FREEERR:
         return MXE_FREEERR_NAME;
      case MXE_PARSEERR:
         return MXE_PARSEERR_NAME;
      case MXE_FILEOERR:
         return MXE_FILEOERR_NAME;
      case MXE_FILEWERR:
         return MXE_FILEWERR_NAME;
      case MXE_FILERERR:
         return MXE_FILERERR_NAME;
      case MXE_FILECERR:
         return MXE_FILECERR_NAME;
      case MXE_FILENOTFOUND:
         return MXE_FILENOTFOUND_NAME;
      case MXE_PROCERR:
         return MXE_PROCERR_NAME;
      case MXE_TIMEOUT:
         return MXE_TIMEOUT_NAME;
      case MXE_LIBERR:
         return MXE_LIBERR_NAME;
      case MXE_GNUPLOTERR:
         return MXE_GNUPLOTERR_NAME;
      case MXE_LAPACKERR:
         return MXE_LAPACKERR_NAME;   
      default:
         return "?";
   }
}

std::string MXE_CodeToDescription( int ec )
{
   
   switch(ec)
   {
      case MXE_INVALIDARG:
         return MXE_INVALIDARG_MSG;
      case MXE_INVALIDCONFIG:
         return MXE_INVALIDCONFIG_MSG;
      case MXE_NOTIMPL:
         return MXE_NOTIMPL_MSG;
      case MXE_PARAMNOTSET:
         return MXE_PARAMNOTSET_MSG;
      case MXE_ENVNOTSET:
         return MXE_ENVNOTSET_MSG;
      case MXE_NOTFOUND:
         return MXE_NOTFOUND_MSG;
      case MXE_SIZEERR:
         return MXE_SIZEERR_MSG;
      case MXE_ALLOCERR:
         return MXE_ALLOCERR_MSG;
      case MXE_FREEERR:
         return MXE_FREEERR_MSG;
      case MXE_PARSEERR:
         return MXE_PARSEERR_MSG;
      case MXE_FILEOERR:
         return MXE_FILEOERR_MSG;
      case MXE_FILEWERR:
         return MXE_FILEWERR_MSG;
      case MXE_FILERERR:
         return MXE_FILERERR_MSG;
      case MXE_FILECERR:
         return MXE_FILECERR_MSG;
      case MXE_FILENOTFOUND:
         return MXE_FILENOTFOUND_MSG;
      case MXE_PROCERR:
         return MXE_PROCERR_MSG;
      case MXE_TIMEOUT:
         return MXE_TIMEOUT_MSG;
      case MXE_LIBERR:
         return MXE_LIBERR_MSG;
      case MXE_GNUPLOTERR:
         return MXE_GNUPLOTERR_MSG;
      case MXE_LAPACKERR:
         return MXE_LAPACKERR_MSG;   
      default:
         return "Unknown mxlib error code.";
   }
}

std::string errno_CodeToName( int ec )
{
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
         return "?";
   }
}


std::string error_report( const std::string & source,   
                          const int & code,             
                          const std::string & file,     
                          const int & line,             
                          const std::string & expl
                        )
{
   
   std::string codeName, codeMessage;
   
   codeName = MXE_CodeToName(code);
   codeMessage = MXE_CodeToDescription(code);
   
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

std::string errno_report( const std::string & source,
                          int ec,
                          const std::string & file, 
                          const int & line,
                          const std::string & expl
                        )
{
   std::string codeName, codeMessage;
   
   codeName = errno_CodeToName(ec);
   codeMessage = strerror(ec);

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

