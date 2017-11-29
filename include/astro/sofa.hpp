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
   #include <sofa.h>
}

//Have to undefine DC to avoid symbol conflict with Eigen.
#undef DC



#endif //__sofa_hpp__
