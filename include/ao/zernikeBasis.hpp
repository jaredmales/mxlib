/** \file zernikeBasis.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Generating a Zernike basis.
  * \ingroup mxAO_files
  * 
  */

#ifndef __zernikeBasis_hpp__
#define __zernikeBasis_hpp__

#include "../sigproc/zernike.hpp"
#include "../improc/fitsFile.hpp"


namespace mx
{
   
namespace AO
{
   
///Make the Zernike basis
/** 
  * \param [in] basisName the name of the basis (not including the mx::AO path)
  * \param [in] dim the linear size of the maps, that is they will be dimxdim in size.
  * \param [in] N is the number of degrees of freedom.  Number of modes will be (N+1)(N+1) - 1.
  *
  * \tparam realT the real numeric type for calculations
  */  
template<typename realT>
void makeZernikeBasis( const std::string & basisName,
                       const std::string & pupilName,
                       int dim,
                       int N )
{
   improc::eigenCube<realT> rawModes;

   rawModes.resize(dim, dim, N);
   zernikeBasis( rawModes);

   std::string pupilFName = mx::AO::path::pupil::pupilFile(pupilName);
   Eigen::Array<realT, -1, -1> pupil;
   
   mx::improc::fitsFile<realT> ff;
   ff.read(pupil, pupilFName);
   
   realT psum = pupil.sum();
   realT norm;
   
   for(int i=0; i< rawModes.planes(); ++i)
   {
      rawModes.image(i) *= pupil;
      
      norm = rawModes.image(i).square().sum()/psum;
      
      rawModes.image(i)/= sqrt(norm);
   }
 
 
 
//    realT p2v;
//    for(int i=0; i<N; ++i)
//    {
//       p2v = rawModes.image(i).maxCoeff() - rawModes.image(i).minCoeff();
//       rawModes.image(i) /= p2v;
//    }
   
   
   
   
   
   
   
   std::string fName = mx::AO::path::basis::modes(basisName, true);
      
   ff.write(fName, rawModes);
}

} //namespace AO

   
} //namespace mx

   
#endif //__fourierBasis_hpp__
