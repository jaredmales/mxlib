/** \file varmapToImage.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief A utility to convert a wavefront variance map to an intensity image
  * \ingroup mxAOAnalytic_files
  * 
  */

#ifndef __varmapToImage_hpp__
#define __varmapToImage_hpp__

namespace mx
{
namespace AO
{
   
///Convert a wavefront variance map to an intensity image by convolving with the PSF.
/**  
  * The PSF should have odd dimensions, and the peak pixel should be in the center pixel
  * defined by [0.5*psf.rows(), 0.5*psf.cols()].  The platescale (lambda/D per pixel) of the PSF should match that of
  * the variance map. The PSF should be normalized such that the center/peak pixel has value 1.
  * 
  * \param [out] im is the intensity image, resized to match varmap
  * \param [in] varmap is the wavefront variance map
  * \param [in] psf is the point spread function
  * 
  * \tparam imageT is an Eigen-array-like type.
  */
template< typename imageT >
void varmapToImage( imageT & im,
                    imageT & varmap,
                    imageT & psf )
{
   typedef typename imageT::Scalar floatT;

   im.resize(varmap.rows(), varmap.cols());
   im.setZero();


   floatT psfNorm = 1.0;//psf.sum();
   
   int psf_i, psf_j;
   floatT psfVal;
   for(int i=0; i< im.rows(); ++i)
   {
      for(int j=0; j< im.cols(); ++j)
      {
         for( int ii =0; ii < im.rows(); ++ii)
         {
            psf_i = 0.5*(psf.rows()) - (ii - i);
            
            if( psf_i < 0 || psf_i >= psf.rows()) continue;
            
            for( int jj =0; jj < im.cols(); ++jj)
            {
               psf_j = 0.5*(psf.cols()) - (jj - j);

               if( psf_j >= 0 and psf_j < psf.cols() )
               {
                  psfVal = psf(psf_i, psf_j)/psfNorm;
                  im(i, j) += varmap(ii, jj) * psfVal;
               }
            }
         }
      }
   }
}

} //namespace AO
} //namespace mx

#endif //__varmapToImage_hpp__

