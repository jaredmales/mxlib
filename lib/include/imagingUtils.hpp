#ifndef __imaging_hpp__
#define __imaging_hpp__

namespace mx
{

#include <cmath>


/// Fill in an Eigen Array with a circular pupil mask.
/** 
  * \ingroup imaging
  *
  * \param m is the allocated Array
  * \param eps [optional] is the central obscuration.  0-1, default is 0.
  * \param rad [optional] is the desired radius. if 0 the maximum radius is used.
  */  
template<class arrayT> 
void circular_pupil( arrayT & m, 
                     typename arrayT::Scalar eps=0, 
                     typename arrayT::Scalar rad=0 
                   )
{
   size_t l0 = m.rows();
   size_t l1 = m.cols();
   
   typename arrayT::Scalar r;
   typename arrayT::Scalar xc = 0.5*(l0-1);
   typename arrayT::Scalar yc = 0.5*(l1-1);
   
   if(rad == 0) rad = 0.5*std::min(l0-1, l1-1);
   
   for(size_t i=0; i < l0; i++)
   {
      for(size_t j=0; j < l1; j++)
      {
         r = std::sqrt( std::pow(i-xc, 2) + std::pow(j-yc, 2) );
         
         if(r <= rad+0.5 && r >= eps*rad) m(i,j) = 1;
         else m(i,j) = 0;
      }
   }
}

// template<class carrayT>
// void fraunpropFocal(carrayT & cmplxAmpFocal, carrayT & cmplxAmpPupil, int Nfft=-1)
// {
//    typedef typename carrayT::Scalar complexT;
//    
//    if(Nfft == -1) Nfft = cmplxAmpPupil.rows();
//    
//    carrayT pad(Nfft, Nfft);
//    pad.setZero();
//    
//    pad.topLeftCorner(cmplxAmpPupil.rows(), cmplxAmpPupil.cols()) = cmplxAmpPupil;
//    
//    //Center the result by multiplying the input by (-1)^(ii+jj) and 1/2 pixel phase slope
//    int one  = 1.0;
//    float xc = .5*(Nfft-1);
//    float yc = .5*(Nfft-1);
//    for(int ii=0; ii< cmplxAmpIn.rows(); ++ii)
//    {
//       for(int jj=0; jj<cmplxAmpIn.cols(); ++jj)
//       {     
//          pad(ii,jj) = pad(ii,jj)*complexT(one,one)* exp(complexT(0.,-6.28*((ii-xc)+(jj-yc))*.5/Nfft));
//          one *= -1;
//       }
//       one*=-1.;
//    }
//    
//    
//    cmplxAmpOut.resize(Nfft, Nfft);
//    
//    fft(cmplxAmpOut, pad);
// }

template<typename arithT>
void makeComplexPupil(Eigen::Array<std::complex<arithT>, Eigen::Dynamic, Eigen::Dynamic> & complexPupil, 
                       Eigen::Array<arithT, Eigen::Dynamic, Eigen::Dynamic> & realPupil, int wavefrontSizePixels)
{
   
   complexPupil.resize(wavefrontSizePixels, wavefrontSizePixels);
   complexPupil.setZero();
   
   //complexPupil.bottomRightCorner(realPupil.rows(), realPupil.cols()) = realPupil*std::complex<arithT>(1,0);
   
   int bl = 0.5*(complexPupil.rows()-1) - 0.5*(realPupil.rows()-1.);
   //int ur = 0.5*(complexPupil.rows()-1) + 0.5*(realPupil.rows()-1.);
   
   complexPupil.block(bl, bl, realPupil.rows(), realPupil.rows()) = realPupil*std::complex<arithT>(1,0);

}

} //namespace mx

#endif //__imaging_hpp__

