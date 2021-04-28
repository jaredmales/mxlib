/** \file jincFuncs.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief The +/- Jinc functions for analyzing contrast
  * \ingroup mxAO_files
  *
  */

#ifndef jincFuncs_hpp
#define jincFuncs_hpp

#include "../../math/constants.hpp"
#include "../../math/func/jinc.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{

/// The +/- Jinc functions for analyzing contrast
/** See Males and Guyon
  *
  */ 
template<typename realT>
int jincFuncs( realT & Jp, ///< [out]
               realT & Jm, ///< [out]
               int mpos,   ///< [in]
               int npos,   ///< [in]
               int mmode,  ///< [in]
               int nmode   ///< [in]
             )
{
   realT ku = mpos;
   realT kv = npos;
   
   realT kp = sqrt( pow(ku + mmode,2) + pow(kv + nmode,2) );
   realT km = sqrt( pow(ku - mmode,2) + pow(kv - nmode,2) );

   Jp = math::func::jinc(math::pi<realT>()*kp);

   Jm = math::func::jinc(math::pi<realT>()*km);
   
   return 0;
}




} //namespace analysis
} //namespace AO
} //namespace mx

#endif //jincFuncs_hpp
