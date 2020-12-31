/** \file gramSchmidt.hpp
  * \brief Procedures to orthogonalize vector basis sets
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing_files
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

#ifndef gramSchmidt_hpp
#define gramSchmidt_hpp

namespace mx
{
   
namespace sigproc 
{
   
///Perform Gram-Schmidt ortogonalization of a basis set, and normalize the result.
/** Performs the stabilized Gram-Schmidt procedure on the input basis set, followed
  * by normalization of the result.
  *
  * \param out [out] is the orthonormal basis set constructed from the input
  * \param int [in] is a basis set, where each column represents one vector.
  * 
  * \tparam progress if true, then the loop index is printed for progress reporting
  * \tparam eigenTout is the Eigen array type of the desired output
  * \tparam eigenTin is the Eigen array type of the input
  * 
  * \ingroup signal_processing 
  */ 
template<int progress=0, typename eigenTout, typename eigenTin>
void gramSchmidt(eigenTout & out, const eigenTin & in)
{
   out.resize(in.rows(), in.cols());

   out.col(0) = in.col(0);
   
   for(int i=1;i< in.cols(); ++i)
   {
      if(progress)
      {
         std::cout << i+1 << "/" << in.cols() << "\n";
      }
      
      //out.col(i) = in.col(i);
      
      out.col(i) = in.col(i) - ((in.col(i).matrix().dot(out.col(0).matrix())) / (out.col(0).matrix().dot(out.col(0).matrix())) )* out.col(0);

      for(int j=1;j<i; ++j)
      {       
         out.col(i) = out.col(i) - ((out.col(i).matrix().dot(out.col(j).matrix()))/(out.col(j).matrix().dot(out.col(j).matrix())))* out.col(j);
      }
   }
   
   for(int i=0; i<out.cols(); ++i)
   {
      out.col(i) = out.col(i)/ out.col(i).matrix().norm();
   }
}


///Perform Gram-Schmidt ortogonalization of a basis set on a window, and normalize the result.
/** Performs the stabilized Gram-Schmidt procedure on the input basis set over a window (or
  * weight function), followed by normalization of the result.
  *
  * \param out [out] is the orthonormal basis set constructed from the input
  * \param in [in] is a basis set, where each column represents one vector.
  * \param window [in] is the window, or weighting function
  * 
  * \tparam progress if true, then the loop index is printed for progress reporting
  * \tparam eigenTout is the Eigen array type of the desired output
  * \tparam eigenTin is the Eigen array type of the input
  * \tparam eigenTWin is the Eigen array type of the window
  * 
  * \ingroup signal_processing 
  */ 
template<int progress=0, typename eigenTout, typename eigenTin, typename eigenTWin>
void gramSchmidt(eigenTout & out, const eigenTin & in, const eigenTWin & window)
{
   //out.resize(in.rows(), in.cols());

   out.col(0) = in.col(0);
   
   for(int i=1;i< in.cols(); ++i)
   {
      if(progress)
      {
         std::cout << i+1 << "/" << in.cols() << "\n";
      }
      
      //out.col(i) = in.col(i);
      
      out.col(i) = in.col(i) - (((in.col(i)*window).matrix().dot(out.col(0).matrix())) / ( (out.col(0)*window).matrix().dot(out.col(0).matrix())) )* out.col(0);

      for(int j=1;j<i; ++j)
      {       
         out.col(i) = out.col(i) - ((  (out.col(i)*window).matrix().dot(out.col(j).matrix()))/((out.col(j)*window).matrix().dot(out.col(j).matrix())))* out.col(j);
      }
   }
   
   for(int i=0; i<out.cols(); ++i)
   {
      out.col(i) = out.col(i)/ (out.col(i)*window.sqrt()).matrix().norm();
   }
}



//Unwraps the gram schmidt coefficients to give a spectrum in terms of the original basis set
//Helper function for gramSchmidtSpectrum (below)
template<int progress=0, typename eigenT>
void baseSpectrum(eigenT & bspect, eigenT & gsspect)
{
   bspect.resize(gsspect.rows(), gsspect.cols());
   bspect.setZero();
   
//#pragma omp parallel for
   for(int i=0; i< gsspect.rows(); ++i)
   {      
      bspect(i,i) = gsspect(i,i);
      
      for(int j=i-1; j >=0; --j)
      {
         bspect.row(i) -= gsspect(i,j)*bspect.row(j);
      }
   }
}

///Perform Gram-Schmidt ortogonalization of a basis set, and normalize the result, while recording the spectrum.
/** Performs the stabilized Gram-Schmidt procedure on the input basis set, followed
  * by normalization of the result.  Also records the spectrum, that is the coefficients of the linear expansion
  * in the orginal basis set for the resultant basis set.
  *
  * 
  * \tparam progress if true, then the loop index is printed for progress reporting
  * \tparam eigenTout is the Eigen array type of the output orthogonalized array
  * \tparam eigenTout2 is the Eigen array type of the spectrum
  * \tparam eigenTin is the Eigen array type of the input
  * 
  * \ingroup signal_processing 
  */ 
template<int progress=0, typename eigenTout, typename eigenTout2, typename eigenTin>
void gramSchmidtSpectrum( eigenTout & out,                        ///< [out] the orthonormal basis set constructed from the input
                          eigenTout2 & spect,                     ///< [out] the spectrum 
                          const eigenTin & in,                    ///< [in] a basis set, where each column represents one vector
                          typename eigenTin::Scalar normPix = 1.0 ///< [in] [optional] area of (usually number of pixels in) the orthogonal region for normalization.
                        )
{
   typedef typename eigenTout::Scalar Scalar;
   
   out.resize(in.rows(), in.cols());
   
   eigenTout2 gsspect;
   
   gsspect.resize(in.cols(), in.cols());
   gsspect.setZero();
   
   out.col(0) = in.col(0);
   gsspect(0,0) = 1;
   
   for(int i=1;i< in.cols(); ++i)
   {
      if(progress)
      {
         std::cout << i+1 << "/" << in.cols() << "\n";
      }
      
      gsspect(i,i)=1;
      
      gsspect(i,0) = ((in.col(i).matrix().dot(out.col(0).matrix())) / (out.col(0).matrix().dot(out.col(0).matrix())) );
      out.col(i) = in.col(i) - gsspect(i,0) * out.col(0);

      for(int j=1;j<i; ++j)
      {       
         gsspect(i,j) = ((out.col(i).matrix().dot(out.col(j).matrix()))/(out.col(j).matrix().dot(out.col(j).matrix())));
         out.col(i) = out.col(i) - gsspect(i,j)* out.col(j);
      }
   }
   
   //Here we unwrap the gram schmidt coefficients, giving us the coefficients in the original basis
   baseSpectrum(spect, gsspect);
   
   Scalar norm;
      
   for(int i=0; i<out.cols(); ++i)
   {
      norm = sqrt(out.col(i).square().sum() / normPix);
      
      out.col(i) /= norm;
      //spect.row(i) /= norm;
   }
}



} //namespace sigproc 
} //namespace mx

#endif //gramSchmidt_hpp


