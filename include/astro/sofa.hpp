/** \file sofa.hpp
  * \brief Wrapper for the sofa library headers, adding a namespace
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup timeutils
  *
  */
  
#ifndef __sofa_hpp__
#define __sofa_hpp__

namespace sofa
{
#ifndef MX_GLOBAL_SOFA
   #include "sofa/sofa.h"
#else
   #include <sofa.h>
#endif
}

//Have to undefine DC to avoid symbol conflict with Eigen.
#undef DC



#endif //__sofa_hpp__
