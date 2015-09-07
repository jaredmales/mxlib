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

///Perform Gram-Schmidt ortogonalization of a basis set, and normalize the result.
/** Performs the stabilized Gram-Schmidt procedure on the input basis set, followed
  * by normalization of the result.
  *
  * \param out [out] is the orthonormal basis set constructed from the input
  * \param int [in] is a basis set, where each column represents one vector.
  * 
  * \tparam progress if true, then the loop index is printed for progress reporting
  * \tparam eigenTin is the Eigen array type of the desired output
  * \tparam eigenTout is the Eigen array type of the input
  */ 
template<int progress=0, typename eigenTin, typename eigenTout>
void gramSchmidt(eigenTin & out, const eigenTout & in)
{
   out.resize(in.rows(), in.cols());

   out.col(0) = in.col(0);
   
   for(int i=1;i< in.cols(); ++i)
   {
      if(progress)
      {
         std::cout << i+1 << "/" << in.cols() << "\n";
      }
      
      out.col(i) = in.col(i);
      
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


#endif //__gramSchmidt_hpp__


