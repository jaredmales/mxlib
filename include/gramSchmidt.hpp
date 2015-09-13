/** \file gramSchmidt.hpp
  * \brief Procedures to orthogonalize vector basis sets
  * 
  * \author Jared R. Males (jaredmales@gmail.com)
  * 
  * \ingroup signal_processing
  *
  */

#ifndef __gramSchmidt_hpp__
#define __gramSchmidt_hpp__

namespace mx
{
   
/** \addtogroup signal_processing
  * @{
  */
     
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
  */ 
template<int progress=0, typename eigenTout, typename eigenTin, typename eigenTWin>
void gramSchmidt(eigenTout & out, const eigenTin & in, const eigenTWin & window)
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


///@}
} //namespace mx

#endif //__gramSchmidt_hpp__


