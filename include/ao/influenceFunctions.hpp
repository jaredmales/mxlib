/** \file influenceFunctions.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Utilities for generating and analyzing deformable mirror influence functions.
  * \ingroup mxAO_files
  * 
  */

#ifndef __influenceFunctions_hpp__
#define __influenceFunctions_hpp__


#include "../improc/eigenCube.hpp"
#include "../math/func/gaussian.hpp"
#include "../improc/fitsFile.hpp"

#include <mx/eigenUtils.hpp>
#include <mx/gnuPlot.hpp>
#include <mx/pout.hpp>
#include <mx/randomT.hpp>

#include <vector>

#include "aoPaths.hpp"

namespace mx
{
   
namespace AO
{
 
template<typename realT>
struct influenceFunctionGaussianSpec
{
   std::string dmName;
   realT dmSz;
   int linNAct;
   realT diameter;
   realT coupling;
   realT pupilSz;
   bool offsetOdd;
   
   std::vector<int> badActi;
   std::vector<int> badActj;
   std::vector<int> badActType;
   
   influenceFunctionGaussianSpec()
   {
      pupilSz = 0;
      offsetOdd =false;
   }
    
   void influenceFunctionsGaussian()
   {
      if(pupilSz <= 0) pupilSz = dmSz;
   
      realT pupilScale = pupilSz/dmSz;
   
      realT act_space = (dmSz - 1.0) / linNAct;  
   
      std::vector<realT> xcoords, ycoords;
   
      realT xc, yc, xcen, ycen;
      realT _rad;
   
   
      //The center of the pupil
      xcen = 0.5*(pupilSz - 1); //0.5*(linNAct*pupilScale - 1);
      ycen = 0.5*(pupilSz - 1); //0.5*(linNAct*pupilScale - 1);   
   
      realT xoff = 0.0;
   
      int q = 0;
      for(int i=0; i< linNAct; ++i)
      {
         xc = (i+0.5)*act_space-0.5*(dmSz-1); //coordinate in pixels w.r.t. dm, not pupil
      
         for(int j=0; j< linNAct; ++j)
         {
            int skip = 0;
               
            if(badActi.size() > 0 && badActi.size() == badActj.size())
            {
               for(int k=0; k< badActi.size(); ++k)
               {
                  if( badActi[k] == i && badActj[k] == j)
                  {
                     skip = 1;
                  }
               }
            }
            
            if(skip) continue;
            
            yc = (j+0.5)*act_space-0.5*(dmSz-1); //coordinate in pixels w.r.t dm, not pupil
       
            if(offsetOdd)
            {
               if( j % 2) xoff = 0.0;
               else xoff = 0.5*act_space;
            }
         
            _rad = sqrt( pow(xc + xoff,2) + pow(yc,2));
   
            if(_rad/act_space <= 0.5*diameter +.01  || diameter == 0)
            {
               xcoords.push_back(xcen + (xc+xoff) + 0.0*(dmSz-pupilSz));
               ycoords.push_back(ycen + yc + 0.0*(dmSz-pupilSz));
               ++q;
            }
         }
      }

      std::cout << "Number of actuators: " << q << " " << xcoords.size() << "\n";
   
      gnuPlot gp;
      gp.command("set size square");
      gp.plot(xcoords, ycoords, "w p ps 0.75");
      gp.point(xcen, ycen);
      gp.circle(xcen, ycen, 0.5*(pupilSz-1), "w l", "Pupil");
   
      realT fwhm = act_space * 2.*sqrt(log2(2))/sqrt(-log2(coupling));

      mx::eigenCube<realT>  act_inf(pupilSz, pupilSz, xcoords.size());

      for(int i=0; i < xcoords.size(); ++i)
      {
         xc = xcoords[i];
         yc = ycoords[i];

         mx::gaussian2D( act_inf.image(i).data(),pupilSz, pupilSz, 0.0, 1.0, xc, yc, fwhm2sigma(fwhm)); 
      }

      if(badActi.size() == badActj.size() && badActi.size() == badActType.size())
      {
         Eigen::Array<realT,-1,-1> im;
         
         std::cerr << "ok trying bad acts\n";
         std::cerr << badActi.size() << "\n";
         for(int k=0; k< badActi.size(); ++k)
         {
            if(badActType[k] == 0) continue;
            
            im.resize(pupilSz, pupilSz);
            
            xc = (badActi[k]+0.5)*act_space-0.5*(dmSz-1); //coordinate in pixels w.r.t. dm, not pupil
            yc = (badActj[k]+0.5)*act_space-0.5*(dmSz-1); //coordinate in pixels w.r.t dm, not pupil
         
            if(offsetOdd)
            {
               if( badActj[k] % 2) xoff = 0.0;
               else xoff = 0.5*act_space;
            }
         
            _rad = sqrt( pow(xc + xoff,2) + pow(yc,2));
   
            std::cerr << xc << " " << yc << "\n";
            
            if(_rad/act_space <= 0.5*diameter +.01  || diameter == 0)
            {
               xc = (xcen + (xc+xoff) + 0.0*(dmSz-pupilSz));
               yc = (ycen + yc + 0.0*(dmSz-pupilSz));
            
               std::cerr << xc << " " << yc << "\n";
            
               mx::gaussian2D( im.data() ,pupilSz, pupilSz, 0.0, -150.0, xc, yc, fwhm2sigma(fwhm));
            
               for(int i=0; i< xcoords.size(); ++i)
               {
                  act_inf.image(i) += im;
               }
            }
            
         }
      }
      else
      {
         std::cerr << "Didn't do bad atuators!\n";
      }
      
      
      
      
      
      
      
      
      
      mx::fitsFile<realT> ff;
      mx::fitsHeader head;

      head.append("LINNACT", linNAct, "Linear number of actuators");
      head.append("DIAMACT", diameter, "Diameter of DM in actuators");
      head.append("DIAMPIX", diameter*act_space, "Diameter of DM in pixels");
      head.append("ACTSPACE", act_space, "Actuator spacing, in pixels");
      head.append("IFTYPE", "GAUSSIAN", "Type of influence function");
      head.append("COUPLING", coupling, "Coupling at 1 actuator spacing");
      head.append("FWHM", fwhm/act_space, "Resultant FWHM in actuators");
      head.append("OFFSTODD", offsetOdd, "Whether odd rows are offset by 1/2");
   
      std::string fName =  mx::AO::path::dm::influenceFunctions(dmName, true);
      ff.write( fName, act_inf, head);
   
      //--------------------------------------------
      //-- Write coords file --
      //--------------------------------------------
      Eigen::Array<realT, -1, -1> coords(2, xcoords.size());
      for(int i=0; i< xcoords.size(); ++i)
      {
         coords(0,i) = xcoords[i];
         coords(1,i) = ycoords[i];
      }
   
      fName =  mx::AO::path::dm::actuatorPositions(dmName, true);
      ff.write( fName, coords);
   
   }
};


///Create a set of Gaussian influence functions (IFs) on a square grid.
/** The width of the IFs is specified by the coupling, which is the height of the IF at 1 actuator separation.  For a FWHM = 1 actuator, coupling = 0.0625.
  * The IFs are normalized to a height of 1.
  * 
  * \param [in] dmName the name of the deformable mirror (DM) (the mx::AO path will be constructed from this).
  * \param [in] dmSz the linear size of the DM in pixels.
  * \param [in] linNAct the linear number of actuators across the DM.
  * \param [in] diameter is the diameter of the DM, in actuators.  If 0 then the DM is a square.
  * \param [in] coupling is the height of the IF at 1 actuator separation.  E.G., for a FWHM = 1 actuator, set this to 0.0625. 
  * \param [in] pupilSz is the size of the pupil, and therefore the size of the maps (pupilSz x pupilSz),  Is set to dmSz if 0, can be pupilSz \< dmSx for an oversized DM.  The pupil is assumed to be centered.
  * 
  * \tparam realT is the real floating point type in which to do calculations.
  */
template<typename realT>
void influenceFunctionsGaussian( const std::string & dmName,
                                 realT dmSz,
                                 int linNAct,
                                 realT diameter,
                                 realT coupling,
                                 realT couplingRange, ///< Uniformly distributed uncertainy in coupling, in fraction of the coupling.
                                 realT pupilSz = 0,
                                 bool offsetOdd = false )
{
   if(pupilSz <= 0) pupilSz = dmSz;
   
   realT pupilScale = pupilSz/dmSz;
   
   realT act_space = (dmSz - 1.0) / linNAct;  
   
   std::vector<realT> xcoords, ycoords;
   
   realT xc, yc, xcen, ycen;
   realT _rad;
   
   
   //The center of the pupil
   xcen = 0.5*(pupilSz - 1); //0.5*(linNAct*pupilScale - 1);
   ycen = 0.5*(pupilSz - 1); //0.5*(linNAct*pupilScale - 1);   
   
   realT xoff = 0.0;
   
   mx::uniDistT<realT> uniVar; ///< Uniform deviate, used in shiftRandom.
   
   if(couplingRange > 0)
   {
      uniVar.seed();
   }
   
   for(int i=0; i< linNAct; ++i)
   {
      xc = (i+0.5)*act_space-0.5*(dmSz-1); //coordinate in pixels w.r.t. dm, not pupil
      
      for(int j=0; j< linNAct; ++j)
      {
         yc = (j+0.5)*act_space-0.5*(dmSz-1); //coordinate in pixels w.r.t dm, not pupil
         
         if(offsetOdd)
         {
            if( j % 2) xoff = 0.0;
            else xoff = 0.5*act_space;
         }
         
         _rad = sqrt( pow(xc + xoff,2) + pow(yc,2));
   
         if(_rad/act_space <= 0.5*diameter +.01  || diameter == 0)
         {
            xcoords.push_back(xcen + (xc+xoff) + 0.0*(dmSz-pupilSz));
            ycoords.push_back(ycen + yc + 0.0*(dmSz-pupilSz));
         }
      }
   }

   std::cout << "Number of actuators: " << xcoords.size() << "\n";
   
   gnuPlot gp;
   gp.command("set size square");
   gp.plot(xcoords, ycoords, "w p ps 0.75");
   gp.point(xcen, ycen);
   gp.circle(xcen, ycen, 0.5*(pupilSz-1), "w l", "Pupil");
   
   realT fwhm = act_space * 2.*sqrt(log2(2))/sqrt(-log2(coupling));

   mx::eigenCube<realT>  act_inf(pupilSz, pupilSz, xcoords.size());

   realT fw;
   
   for(int i=0; i < xcoords.size(); ++i)
   {
      xc = xcoords[i];
      yc = ycoords[i];

      if(couplingRange > 0)
      {
         fw = fwhm + fwhm*(1.0 - 2.0*uniVar)*couplingRange;
      }
      else
      {
         fw = fwhm;
      }
      
      mx::gaussian2D( act_inf.image(i).data(),pupilSz, pupilSz, 0.0, 1.0, xc, yc, fwhm2sigma(fw)); 

   }

   mx::fitsFile<realT> ff;
   mx::fitsHeader head;

   head.append("LINNACT", linNAct, "Linear number of actuators");
   head.append("DIAMACT", diameter, "Diameter of DM in actuators");
   head.append("DIAMPIX", diameter*act_space, "Diameter of DM in pixels");
   head.append("ACTSPACE", act_space, "Actuator spacing, in pixels");
   head.append("IFTYPE", "GAUSSIAN", "Type of influence function");
   head.append("COUPLING", coupling, "Coupling at 1 actuator spacing");
   head.append("FWHM", fwhm/act_space, "Resultant FWHM in actuators");
   head.append("OFFSTODD", offsetOdd, "Whether odd rows are offset by 1/2");
   
   std::string fName =  mx::AO::path::dm::influenceFunctions(dmName, true);
   ff.write( fName, act_inf, head);
   
   //--------------------------------------------
   //-- Write coords file --
   //--------------------------------------------
   Eigen::Array<realT, -1, -1> coords(2, xcoords.size());
   for(int i=0; i< xcoords.size(); ++i)
   {
      coords(0,i) = xcoords[i];
      coords(1,i) = ycoords[i];
   }
   
   fName =  mx::AO::path::dm::actuatorPositions(dmName, true);
   ff.write( fName, coords);
   
}

///Calculate the pseudo-inverse and mirror modes for a set of influence functions.
/** The influence functions are a cube of 2D maps of mirror deformations, and the matrix \f$ A \f$ is formed by vectorizing the influence functions. 
  * The pseudo inverse is then calculated from the SVD as follows:
  * \f$ A = U S V^T \f$
  * which gives the Moore-Penrose pseudo-inverse:
  * \f$ A^+ = V S^+ U^T\f$.
  * The columns of U contain the orthogonal mirror modes, and S contains the singular values.  These are both written to disk along with the pseudo-inverse.
  * 
  * \param[in] dmName is the name of the deformable mirror.
  * \param[in] maxCondition [optional] the maximum condition number to accept for the pseudo-inverse.  If < 0 [default] then the user is interactively asked what to use.
  * 
  * \tparam realT is the real floating point type in which to do calculations.
  */ 
template<typename realT>
void ifPInv( const std::string & dmName,
             realT maxCondition = -1)
{
   std::string pinvName = mx::AO::path::dm::pseudoInverse(dmName);
   std::string mmodesName = mx::AO::path::dm::mirrorModes(dmName, true);
   std::string svsName =  mx::AO::path::dm::singularValues(dmName);
   std::string ifName = mx::AO::path::dm::influenceFunctions(dmName);
   
   std::cout << "mmodesName: " <<  mmodesName << "\n";
   
   
   
   mx::eigenCube<realT>  actInf;
         
   mx::fitsFile<realT> ff;
   ff.read(ifName, actInf);
   
   int nrows = actInf.rows();
   int ncols = actInf.cols();
   int nplanes = actInf.planes();
   
   Eigen::Array<realT, -1, -1> rowInf=actInf.asVectors();
   actInf.resize(0,0,0);
   
   Eigen::Array<realT, -1,-1> PInv;
   
   realT condition;
   int nRejected;
   
   Eigen::Array<realT, -1, -1> U, S, VT;
   
   
   int interact = MX_PINV_PLOT | MX_PINV_ASK;
   if(maxCondition >= 0) interact = MX_PINV_NO_INTERACT;
   
   mx::eigenPseudoInverse(PInv, condition, nRejected, U, S, VT, rowInf, maxCondition, interact );
      
   mx::fitsHeader head;
   head.append("MAXCONDN", maxCondition, "Max. condition no. in pseudo inverse");
   head.append("CONDN", condition, "Actual condition number of pseudo inverse");
   head.append("NREJECT", nRejected, "Number of s.v.s rejected");
   
   //std::string pinvFile = pinvName;//ifName + "_pinv.fits";
   
   ff.write(pinvName, PInv, head);
   
   
   //This maps to a cube, but does not copy or take ownership.
   mx::eigenCube<realT> mmodes( U.data(), nrows, ncols, nplanes);
   
   
   ff.write(mmodesName, mmodes, head);
   
   
   std::ofstream fout;
   fout.open(svsName);
   
   for(int i=0;i<S.rows(); ++i)
   {
      fout << S(i) << "\n";
   }
   
   fout.close();
   
}

///Calculate the modes-to-commands matrix for a set of modes
/**
  * Given the pseudo-inverse \f$ A^+ \f$ of the influence functions, the commands to the actuators to take a mirror shape \f$ \vec{s} \f$ are calculated as
  * \f[
   \vec{c} = \mathbf{A^+} \vec{s}.
  * \f]
  * Now given a basis set \f$ \mathbf{M} \f$, the mirror shape is specified as
   \f[
   \vec{s}= \sum_i m_i M_i
   \f]
   * where \f$ \vec{m} = m_0, m_1, \cdot m_i \cdot \f$ are the modal coefficients. If the mirror-to-commands matrix, \f$ \mathbf{M2c} \f$ gives commands as
  * \f[
   \vec{c} = \mathbf{M2c} \vec{m}
  * \f]
  * then we can calculate \f$ \mathbf{M2c} \f$ as:
  * \f[
    \mathbf{M2c} = \mathbf{A}^{+T} \mathbf{M}
   \f]
  * 
  * \param[out] M2c is the M2c matrix, allocated during the calculationg
  * \param[in] Ap is the pseudo-inverse of the influence function basis set
  * \param[in] M is the modal basis set, a cube of shapes
  *
  * \tparam realT is the real floating point type used for all calculations.
  */ 
template<typename realT>
void m2cMatrix( Eigen::Array<realT, -1, -1> & M2c,
                Eigen::Array<realT, -1, -1> & Ap,
                mx::eigenCube<realT> & M )
{
   M2c = Ap.matrix().transpose() * M.asVectors().matrix();
}

///Calculate the modes-to-commands matrix for a set of modes
/**
  * A wrapper for \ref m2cMatrix, where the various matrices are here specified with names, which
  * in turn generate mx::AO paths
  * 
  * \param[in] dmName is the name of the deformable mirror
  * \param[in] basisName is the name of the modal basis
  * \param[in] pupilName is the name of the pupil
  * 
  * \tparam realT is the real floating point type in which to do calculations.
  */ 
template<typename realT>
void m2cMatrix( const std::string & dmName,
                const std::string & basisName )
{
   std::string M2cFName;
   mx::fitsFile<realT> ff;
   
   Eigen::Array<realT, -1,-1> Ap; //The pseudo-inverse
   std::string pinvFName = mx::AO::path::dm::pseudoInverse(dmName);   
   ff.read(pinvFName, Ap);//read from the file
   
   mx::eigenCube<realT>  M; //The modal basis, an image cube of shapes
   std::string basisFName;

   basisFName = mx::AO::path::basis::modes(basisName);
   ff.read(basisFName, M); //read from file
   
   
   Eigen::Array<realT, -1,-1> M2c; //The M2c matrix
   m2cMatrix(M2c, Ap, M);
   
   M2cFName = mx::AO::path::dm::M2c(dmName, basisName, true);
      
   ff.write(M2cFName, M2c);
}

///Generate the modes-to-commands matrix for a set of modes on the modal DM
/**
  * Generates an Identity matrix of the appropriate size.
  * 
  * \param[in] basisName is the name of the modal basis
  * 
  * \tparam realT is the real floating point type in which to do calculations.
  */ 
template<typename realT>
void modalDMM2cMatrix( const std::string & basisName )
{
   mx::fitsFile<realT> ff;
   
   
   mx::eigenCube<realT>  M; //The modal basis, an image cube of shapes

   std::string basisFName;
   
   basisFName = mx::AO::path::basis::modes(basisName);
   
   ff.read(basisFName, M); //read from file
   
   
   Eigen::Array<realT, -1,-1> M2c; //The M2c matrix

   M2c.resize(M.planes(), M.planes());
   M2c.setZero();
   
   for(int i=0; i< M2c.rows(); ++i)
   {
      M2c(i,i) = 1;
   }
         
   
   std::string M2cFName;

   M2cFName = mx::AO::path::dm::M2c("modalDM", basisName, true);
      
   ff.write(M2cFName, M2c);
}

///Calculate the basis set as it will be reproduced by the DM.
/** The projected modes are written to the path specified by mx::AO::path::dm::projectedModes().
 * 
  * \param [in] dmName is the name of the DM.
  * \param [in] basisName is the name of the basis.
  * 
  * \tparam realT is the real floating point type in which to do calculations.
  */ 
template<typename realT>
void modesOnDM( const std::string & dmName,
                const std::string & basisName )
{
   Eigen::Array<realT, -1, -1> m2c ;//, m, c;
   
   
   std::string m2cName;

   m2cName = mx::AO::path::dm::M2c(dmName, basisName );
   
   mx::fitsFile<realT> ff;
   ff.read(m2cName, m2c);
   
   std::cout << m2c.rows() << "\n";
   std::cout << m2c.cols() << "\n";
   
   int nmodes = m2c.cols();
   int nact = m2c.rows();
   
//   m.resize(nmodes,1);
   
   
   mx::eigenCube<realT> act_inf, projModes;
   
   std::string fName =  mx::AO::path::dm::influenceFunctions(dmName, true);
   ff.read(fName, act_inf);
   
   projModes.resize(act_inf.rows(), act_inf.cols(), nmodes);
   projModes.setZero();
   
#pragma omp parallel for
   for(int i=0; i< nmodes; ++i)
   {
      Eigen::Array<realT, -1, -1> m, c;
         
      std::cout << i+1 << "/" << nmodes << "\n";
      m.resize(nmodes,1);

      m.setZero();
      m(i,0) = 1.0;
      c = m2c.matrix() * m.matrix();
    
      
      //#pragma omp parallel
      //{
    //     Eigen::Array<realT, -1, -1> tmp(act_inf.rows(), act_inf.cols());
     //    tmp.setZero();
         
        // #pragma omp for
         for(int j=0; j< nact; ++j)
         {
            projModes.image(i) += c(j,0) * act_inf.image(j);
  //          tmp += c(j,0) * act_inf.image(j);
         }
         //#pragma omp critical
//         projModes.image(i) += tmp;
      //}
   }
   
   std::string oName;
   oName = mx::AO::path::dm::projectedModes(dmName, basisName, true );
   
   ff.write(oName, projModes);
}

} //namespace AO
} //namespace mx

#endif //__influenceFunctions_hpp__


