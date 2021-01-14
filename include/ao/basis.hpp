/** \file basis.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Utilities for working with a modal basis.
  * \ingroup mxAO_files
  * 
  */

#ifndef __basis_hpp__
#define __basis_hpp__

#include "../ioutils/fits/fitsFile.hpp"
using namespace mx::fits;

#include "../sigproc/fourierModes.hpp"
#include "../sigproc/gramSchmidt.hpp"
//#include "../improc/eigenUtils.hpp"
#include "../improc/imageFilters.hpp"
#include "../improc/imagePads.hpp"
#include "../improc/eigenCube.hpp"

#include "../sigproc/signalWindows.hpp"
using namespace mx::improc;
using namespace mx::sigproc;

#include "../math/eigenLapack.hpp"
using namespace mx::math;

#include "aoPaths.hpp"

namespace mx
{
   
namespace AO
{

///Multiply a raw modal basis by a pupil mask.   
template<typename realT>
void applyPupil2Basis( eigenCube<realT> & modes,
                       const std::string & basisName,
                       const std::string & pupilName,
                       realT fwhm = 0  )
{
   fitsFile<realT> ff;
   
   std::string rawFName =  mx::AO::path::basis::modes(basisName);
   ff.read(modes, rawFName);
   
   std::string pupilFName = mx::AO::path::pupil::pupilFile(pupilName);
   Eigen::Array<realT, -1, -1> pupil;
   ff.read(pupil, pupilFName);
   
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


template<typename realT, typename spectRealT>
int orthogonalizeBasis( eigenCube<realT> & ortho,
                        Eigen::Array<spectRealT, -1, -1>  & spect,
                        eigenCube<realT> & modes,
                        int method,
                        realT psum = 1.0
                      )
{
      
   
   if(method == MXAO_ORTHO_METHOD_SGS)
   {
      ortho.resize(modes.rows(), modes.cols(), modes.planes());
      Eigen::Map< Eigen::Array<realT, -1, -1>> gsout(ortho.data(), ortho.rows()*ortho.cols(), ortho.planes());
      
      gramSchmidtSpectrum<1>(gsout, spect, modes.asVectors(), psum);
    
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
      
      eigenCube<realT> omodes( U.data(), modes.rows(), modes.cols(), modes.planes());
      
      ortho.resize(modes.rows(), modes.cols(), modes.planes());
      ortho = omodes;
      
   }
   else
   {
      std::cerr << "Invalid orthogonalization method.\n";
      return -1;
   }
   
   
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
   fitsFile<realT> ff;
   fitsHeader head;
   
   
   eigenCube<realT> modes;
   applyPupil2Basis(modes, basisName, pupilName);
   
   eigenCube<realT> ortho;
   Eigen::Array<realT, -1, -1> spect;
   
   orthogonalizeBasis(ortho, spect, modes, method);
   
   modes.resize(1,1,1);
   
   std::string pupilFName = mx::AO::path::pupil::pupilFile(pupilName);
   Eigen::Array<realT, -1, -1> pupil;
   ff.read(pupil, pupilFName);
   
   realT psum = pupil.sum();
   realT norm;
   
  for(int i=0; i< ortho.planes(); ++i)
  {
     norm = ortho.image(i).square().sum()/psum;
     ortho.image(i)/= sqrt(norm);
      
      if(method == MXAO_ORTHO_METHOD_SGS)
      {
         //if(i == 0) std::cout << sqrt(norm) << "\n";
         
         //spect.row(i) /= sqrt(norm);
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

template<typename spectRealT, typename realT>
void normalizeSpectrum( mx::improc::eigenImage<spectRealT> & spectrum,
                        mx::improc::eigenCube<realT> & modes,
                        mx::improc::eigenCube<realT> & ortho,
                        mx::improc::eigenImage<realT> & pupil
                      )
{
   int nModes = modes.planes();
   
   std::vector<realT> w(nModes);
   
   size_t psum = pupil.sum();
   
   #pragma omp parallel
   {
      std::vector<realT> amps( nModes, 0.0 );
      std::vector<realT> uwAmps(nModes);
      std::vector<realT> rat(psum);
      
      #pragma omp for
      for(int k=0; k< nModes; ++k)
      {
         amps[k] = 1.0;
      
         for(int j=0; j< nModes; ++j)
         {
            uwAmps[j] = 0;// amps(j, k)*spectrum(j,j);
         
            for(int i= j; i < nModes;  ++i)
            {
               uwAmps[j] += amps[i]*spectrum(i,j);
            }
         }
      
         mx::improc::eigenImage<realT> uwp(modes.rows(), modes.cols());
         uwp.setZero();
   
         for(int i=0;i<nModes; ++i) uwp+= uwAmps[i]*modes.image(i);
    
         int idx = 0;
         for(int r=0; r< pupil.rows(); ++r)
         {
            for(int c=0;c < pupil.cols(); ++c)
            {
               if(pupil(r,c) == 1) rat[idx++] = ( uwp(r,c) / ortho.image(k)(r,c));
            }
         }  
      
         w[k] = mx::math::vectorMedianInPlace(rat);
      
         if(k == 1200) std::cout << w[k] << "\n";
         
         amps[k] = 0.0;
      }
   }
    
   for(int i=0; i< nModes; ++i) spectrum.row(i)/= w[i];
   
   return;
}

template<typename realT>
int slaveBasis( const std::string & outputBasisN,
                const std::string & inputBasisN,
                const std::string & pupilN,
                const std::string & dmN,
                realT fwhm,
                realT fsmooth,
                int firstMode = 0
              )
{
   fitsFile<realT> ff;
   fitsHeader head;

   //get DM size and clear memory first
   std::string dmFName = mx::AO::path::dm::influenceFunctions(dmN);
   eigenCube<realT> inf;
   ff.read(dmFName, inf);

   int dmSz = inf.rows();
   
   inf.resize(0,0);

   //Now load basis
   std::string basisFName = mx::AO::path::basis::modes(inputBasisN);//, pupilName); 
   eigenCube<realT> modes;
   ff.read(basisFName, modes);

   //Now load pupil
   std::string pupilFName = mx::AO::path::pupil::pupilFile(pupilN);
   Eigen::Array<realT, -1, -1> pupil;
   ff.read(pupilFName, pupil);

   Eigen::Array<realT, -1, -1> im, ppim, pim, fim, ppupil;
   
   int padw = 2*fwhm;
   
   if(fwhm ==0 ) padw = 2*fsmooth;

   int cutw = 0.5*(dmSz - modes.rows());
   
   if( padw < cutw )
   {
      padw = cutw;
   }

   eigenCube<realT> oModes;
   
   oModes.resize( modes.rows()+2*cutw, modes.cols()+2*cutw, modes.planes());
   
   
   if(modes.rows() > pupil.rows())
   {
      padImage(ppupil, pupil, 0.5*(modes.rows()-pupil.rows()));
      pupil = ppupil;
   }
   
   
   padImage(ppupil, pupil, padw);

   
   for(int i=0;i<firstMode; ++i)
   {
      im = modes.image(i);
      padImage(ppim, im, padw);
      cutPaddedImage(im, ppim, 0.5*ppim.rows() - (0.5*modes.rows()+cutw));
                  
      oModes.image(i) = im;
   }
   
   for(int  i=firstMode; i< modes.planes(); ++i)
   {

      if(fwhm > 0 || fsmooth > 0)
      {
         
         im = modes.image(i)*pupil;
         padImage(ppim, im, padw);
         padImage(pim, im, ppupil, padw);//ppupil,4);
         
         
         if( fwhm > 0 )
         {
            filterImage(fim, pim, gaussKernel<Eigen::Array<realT, -1, -1>>(fwhm));
            fim = ppim + fim*(ppupil - 1).abs();
         }
         else
         {
            fim = ppim;
         }
         
         if(fsmooth > 0)
         {
            pim = fim;
            
            filterImage(fim, pim, gaussKernel<Eigen::Array<realT, -1, -1>>(fsmooth));
         }
         
         
         cutPaddedImage(im, fim, 0.5*fim.rows() - (0.5*modes.rows()+cutw));
                  
         oModes.image(i) = im; 
      }
   }
   
   std::string oFName =  mx::AO::path::basis::modes(outputBasisN, true);
   ff.write(oFName, oModes);
   
   return 0;
}


template<typename realT>
int apodizeBasis( const std::string & outputBasisN,
                  const std::string & inputBasisN,
                  realT tukeyAlpha,
                  realT centralObs,
                  realT overScan,
                  int firstMode )
{
   fitsFile<realT> ff;
   fitsHeader head;
   
   std::string basisFName = mx::AO::path::basis::modes(inputBasisN);//, pupilName); 
   eigenCube<realT> modes;
   ff.read(basisFName, modes);

   realT cen = 0.5*(modes.rows() - 1.0);
   
   Eigen::Array<realT, -1, -1> pupil;
   pupil.resize( modes.rows(), modes.cols());
   
   if(centralObs == 0)
   {
      window::tukey2d<realT>(pupil.data(), modes.rows(), modes.rows() + overScan, tukeyAlpha, cen,cen);
   }
   else
   {
      window::tukey2dAnnulus<realT>(pupil.data(), modes.rows(), modes.rows() + overScan, centralObs, tukeyAlpha, cen,cen);
   }
   
   
   
   for(int  i=firstMode; i< modes.planes(); ++i)
   {
      modes.image(i) = modes.image(i)*pupil;
   }
   
   std::string oFName =  mx::AO::path::basis::modes(outputBasisN, true);
   ff.write(oFName, modes);
   
   return 0;
}


template<typename realT>
int subtractBasis( const std::string & basisName,
                   const std::string & basisName1,
                   const std::string & basisName2,
                   const std::string & pupilName,
                   int method  )
{
   fitsFile<realT> ff;
   fitsHeader head;
   
   std::string basisFName1 = mx::AO::path::basis::modes(basisName1);//, pupilName); 
   eigenCube<realT> modes1;
   ff.read(basisFName1, modes1);
   
   std::string basisFName2 = mx::AO::path::basis::modes(basisName2);//, pupilName); 
   eigenCube<realT> modes2;
   ff.read(basisFName2, modes2, head);
   realT fwhm = head["FWHM"].Value<realT>();
   
   eigenCube<realT> modes;
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
   
   eigenCube<realT> ortho;
      
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
