/** \file fourierBasis.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Generating a fourier basis.
  * \ingroup mxAO_files
  * 
  */

#ifndef __fourierBasis_hpp__
#define __fourierBasis_hpp__

#include <mx/fourierModes.hpp>
#include <mx/fitsFile.hpp>


namespace mx
{
   
namespace AO
{
   
///Make the modified Fourier basis
/** 
  * \param [in] basisName the name of the basis (not including the mx::AO path)
  * \param [in] dim the linear size of the maps, that is they will be dimxdim in size.
  * \param [in] N is the number of degrees of freedom.  Number of modes will be (N+1)(N+1) - 1.
  *
  * \tparam realT the real numeric type for calculations
  */  
template<typename realT>
void makeModfBasis( const std::string & basisName,
                    int dim,
                    int N )
{
   mx::eigenCube<realT> modes;

   makeFourierBasis_Rect( modes, dim, N, MX_FOURIER_MODIFIED);

//    realT p2v;
//    for(int i=0; i<N; ++i)
//    {
//       p2v = modes.image(i).maxCoeff() - modes.image(i).minCoeff();
//       modes.image(i) /= p2v;
//    }

   mx::fitsFile<realT> ff;
   
   std::string fName = mx::AO::path::basis::modes(basisName, true);
      
   ff.write(fName, modes);
}

} //namespace AO

   
} //namespace mx

   
#endif //__fourierBasis_hpp__
