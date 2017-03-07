/** \file basis.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Utilities for working with a modal basis.
  * \ingroup mxAO_files
  * 
  */

#ifndef __basis_hpp__
#define __basis_hpp__

#include <mx/fourierModes.hpp>
#include <mx/fitsFile.hpp>
#include <mx/gramSchmidt.hpp>
#include <mx/eigenUtils.hpp>
#include <mx/imageFilters.hpp>
#include <mx/imagePads.hpp>

#include "aoPaths.hpp"

namespace mx
{
   
namespace AO
{

///Multiply a raw modal basis by a pupil mask.   
template<typename realT>
void applyPupil2Basis( mx::eigenCube<realT> & modes,
                       const std::string & basisName,
                       const std::string & pupilName,
                       realT fwhm = 0  )
{
   mx::fitsFile<realT> ff;
   //mx::eigenCube<realT> modes;
   
   
   std::string rawFName =  mx::AO::path::basis::modes(basisName);
   ff.read(rawFName, modes);
   
   std::string pupilFName = mx::AO::path::pupil::pupilFile(pupilName);
   Eigen::Array<realT, -1, -1> pupil;
   ff.read(pupilFName, pupil);
   
   Eigen::Array<realT, -1, -1> im, pim, fim;

   for(int i=0; i< modes.planes(); ++i)
   {
      
      modes.image(i) = modes.image(i) * pupil;
      
      if(fwhm > 0)
      {
         im = modes.image(i);
         
         padImage(pim, im, pupil,3);
         
         
         filterImage(fim, pim, gaussKernel<Eigen::Array<realT, -1, -1>>(fwhm));
         
         modes.image(i) += fim* (pupil - 1).abs();
      }
   }

   //std::string procFName = mx::AO::path::basis::pupilModes(basisName, pupilName, true);

   //mx::fitsHeader head;
   //head.append( "FWHM", fwhm, "FWHM of smoothing kernel used for slaving");
   //ff.write(procFName, rawModes, head);
   
}

/** \def MXAO_ORTHO_METHOD_SGS
  * \brief Constant to specify using the stabilized Gramm Schmidt (SGS) orthogonalization procedure.
  */ 
#define MXAO_ORTHO_METHOD_SGS 0

/** \def MXAO_ORTHO_METHOD_SVD
  * \brief Constant to specify using the singular value decomposition (SVD) for orthogonalizing a modal basis.
  */ 
#define MXAO_ORTHO_METHOD_SVD 1


template<typename realT>
int orthogonalizeBasis( mx::eigenCube<realT> & ortho,
                        Eigen::Array<realT, -1, -1>  & spect,
                         mx::eigenCube<realT> & modes,
                         int method )
{
      
   
   if(method == MXAO_ORTHO_METHOD_SGS)
   {
      ortho.resize(modes.rows(), modes.cols(), modes.planes());
      Eigen::Map< Eigen::Array<realT, -1, -1>> gsout(ortho.data(), ortho.rows()*ortho.cols(), ortho.planes());
      
      gramSchmidtSpectrum<1>(gsout, spect, modes.asVectors());
    
      std::cout << "out : " << "\n";
      
      //ortho.asVectors() = gsout;
      
   }
   else if (method == MXAO_ORTHO_METHOD_SVD)
   {
      Eigen::Array<realT, -1, -1> U, S, VT, A;
      int info;
      A = modes.asVectors();
      info = eigenGESDD(U,S,VT,A);
   
      if(info != 0)
      {
         std::cerr  << "Some error occurred in SVD\n";
         return -1;
      }
      
      //This maps to a cube, but does not copy or take ownership.
      
      mx::eigenCube<realT> omodes( U.data(), modes.rows(), modes.cols(), modes.planes());
      
      ortho.resize(modes.rows(), modes.cols(), modes.planes());
      ortho = omodes;
      
   }
   else
   {
      std::cerr << "Invalid orthogonalization method.\n";
      return -1;
   }
   
/*   
   realT p2v;
   for(int i=0; i<ortho.planes(); ++i)
   {
      p2v = ortho.image(i).maxCoeff() - ortho.image(i).minCoeff();
      ortho.image(i) /= p2v;
   }*/
   
   
   return 0;
   
}

///Calculate an orthogonal projection of a basis on a pupil.
/** Can use either the stabilized Gramm Schmidt (SGS) procedure, or Singular Value Decomposition (SVD).
  * This calls \ref applyPupil2Basis as a first step..
  * 
  * \param [in] basisName the name of the basis to orthogonalize
  * \param [in] pupilName the name of the pupil on which to orthogonalize.  
  * \param [in] method either \ref MXAO_ORTHO_METHOD_SGS (for SGS) or \ref MXAO_ORTHO_METHOD_SVD (for SVD)
  *
  * \tparam realT
  */ 
template<typename realT>
void orthogonalizeBasis( const std::string & orthoName,
                         const std::string & basisName,
                         const std::string & pupilName,
                         int method )
{
   mx::fitsFile<realT> ff;
   mx::fitsHeader head;
   
   
   mx::eigenCube<realT> modes;
   applyPupil2Basis(modes, basisName, pupilName);
   
   mx::eigenCube<realT> ortho;
   Eigen::Array<realT, -1, -1> spect;
   
   orthogonalizeBasis(ortho, spect, modes, method);
   
   modes.resize(1,1,1);
   
   std::string pupilFName = mx::AO::path::pupil::pupilFile(pupilName);
   Eigen::Array<realT, -1, -1> pupil;
   ff.read(pupilFName, pupil);
   
   realT psum = pupil.sum();
   realT norm;
   
   for(int i=0; i< ortho.planes(); ++i)
   {
      norm = ortho.image(i).square().sum()/psum;
      
      ortho.image(i)/= sqrt(norm);
      
      
      if(method == MXAO_ORTHO_METHOD_SGS)
      {
         spect.row(i) /= sqrt(norm);
      }
   }
 
   std::string orthoFName =  mx::AO::path::basis::modes(orthoName, true);
   
   std::string orthm = "UNK";
   if(method == MXAO_ORTHO_METHOD_SGS) orthm = "SGS";
   if(method == MXAO_ORTHO_METHOD_SVD) orthm = "SVD";
   
   head.append("ORTHMETH", orthm, "Orthogonalization method.");
   head.append("ORTHPUP", pupilName, "Pupil used for orthogonalization");
   
   
   ff.write(orthoFName, ortho, head);
   
   if(method == MXAO_ORTHO_METHOD_SGS)
   {
      std::string orthoFName =  mx::AO::path::basis::spectrum(orthoName, true);
      ff.write(orthoFName, spect, head);
   }
   
}

template<typename realT>
int subtractBasis( const std::string basisName,
                   const std::string basisName1,
                   const std::string basisName2,
                   const std::string pupilName,
                   int method  )
{
   mx::fitsFile<realT> ff;
   mx::fitsHeader head;
   
   std::string basisFName1 = mx::AO::path::basis::modes(basisName1);//, pupilName); 
   mx::eigenCube<realT> modes1;
   ff.read(basisFName1, modes1);
   
   std::string basisFName2 = mx::AO::path::basis::modes(basisName2);//, pupilName); 
   mx::eigenCube<realT> modes2;
   ff.read(basisFName2, modes2, head);
   realT fwhm = head["FWHM"].Value<realT>();
   
   mx::eigenCube<realT> modes;
   modes.resize( modes1.rows(), modes1.cols(), modes1.planes() + modes2.planes());
   
   std::string pupilFName = mx::AO::path::pupil::pupilFile(pupilName);
   Eigen::Array<realT, -1, -1> pupil;
   ff.read(pupilFName, pupil);
   
   for(int i=0; i < modes1.planes(); ++i)
   {
      modes.image(i) = modes1.image(i)*pupil;
   }
   
   for(int i=0; i < modes2.planes(); ++i)
   {
      modes.image( modes1.planes() + i) = modes2.image(i)*pupil;
   }
   
   mx::eigenCube<realT> ortho;
      
   orthogonalizeBasis(ortho, modes, method);

   //Now copy out just the modes corresponding to the 2nd basis set.
   //Apply the gaussian smooth FWHM again.
   Eigen::Array<realT, -1, -1> im, ppim, pim, fim, ppupil;
   
   modes.resize(ortho.rows(), ortho.cols(), modes2.planes());
   
   realT fsmooth = 6.0;
   
   padImage(ppupil, pupil, 10);
   
   for(int  i=0; i< modes2.planes(); ++i)
   {
      modes.image(i) = ortho.image( modes1.planes() + i)*pupil;
      
      if(fwhm > 0)
      {
         im = modes.image(i);
         padImage(ppim, im, 10);
         padImage(pim, im, 10);//ppupil,4);
         
         filterImage(fim, pim, gaussKernel<Eigen::Array<realT, -1, -1>>(fwhm));
         fim = ppim + fim*(ppupil - 1).abs();
         
         if(fsmooth > 0)
         {
            pim = fim;
            filterImage(fim, pim, gaussKernel<Eigen::Array<realT, -1, -1>>(fsmooth));
         }
         
         cutPaddedImage(im, fim, 10);
         modes.image(i) = im; 
      }
   }
   
   
   
   std::string orthoFName =  mx::AO::path::basis::modes(basisName, true);
   ff.write(orthoFName, modes);
   
   return 0;
   
}

} //namespace mx
} //namespace AO

   
#endif //__basis_hpp__
