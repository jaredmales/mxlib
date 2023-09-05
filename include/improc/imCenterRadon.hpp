/** \file imCenterRadon.hpp
  * \brief A class to find PSF centers using circular symmetry.
  * \ingroup image_processing_files
  * \author Jared R. Males (jaredmales@gmail.com)
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

#ifndef __imCenterRadon_hpp__
#define __imCenterRadon_hpp__

#include "../math/vectorUtils.hpp"
#include "../math/fit/fitGaussian.hpp"
#include "imageFilters.hpp"
#include "eigenImage.hpp"
#include "imagePeakInterp.hpp"

//#define ICCS_OMP
namespace mx
{
namespace improc
{
   

///Find the center of the image of a point source using circular symmetry of the PSF.
/** Conducts a grid search through possible center positions.  Finds the point with the minimum
  * total variance around concentric circles about that point.  The cost grid is fit with a Gaussian.
  * 
  * \tparam realT is a real floating point type.
  * 
  * \ingroup image_reg
  */ 
template<typename transformT>
struct imCenterRadon
{

   typedef typename transformT::arithT realT;

   eigenImage<realT> m_grid; ///< The cost grid.
   
   realT m_maxPix {2}; ///<The maximum range of the grid search, the grid will be +/- maxPix.  Default: 5.0.
   realT m_dPix {0.1}; ///< The spacing of the points in the grid search.  Default: 0.1.
   realT m_peakTol {0.1}; ///< The tolerance for the peak finding.  This is multiplicative with m_dPix.

   realT m_minRad {1.0}; ///< The minimum radius to include in the transform.   Default: 1.0.
   realT m_maxRad {10}; ///< The maximum radius to include in the transform.  Default: 10.0.
   realT m_dRad {0.1}; ///< The spacing of the radii.   Default: 0.1.

   realT m_angTolPix {1.0}; ///< The angle tolerance in pixels at the maximum radius, sets the angular step size. Default is 1.0.

   realT m_guessWidth {0.5}; ///< Guess in original pixels for the FWHM of the cost grid.  Default is 0.5, which seems to work well in most cases.
   
   ///The fitter for finding the centroid of the grid.  You can use this to inspect the details of the fit if desired.
   //mx::math::fit::fitGaussian2D<mx::gaussian2D_gen_fitter<realT>> fit; 
   mx::math::fit::fitGaussian2Dgen<realT> m_fit;
   
   realT m_x0; ///< Internal storage of the guess center position
   realT m_y0; ///< Internal storage of the guess center position.
   
   imCenterRadon()
   {
   }
   

   ///Peform the grid search and fit.
   int center( realT & x, ///< [out] the x-coordinate of the estimated center of rotational symmetry
               realT & y, ///< [out] the y-coordinate estimated center of rotational symmetry
               eigenImage<realT> & im,  ///< [in] the image to centroid
               realT x0,  ///< [in] an initial guess at the x-coordinate, the grid is centered about this point.
               realT y0  ///< [in] an initial guess at the y-coordinate, the grid is centered about this point.
             )
   {
      m_x0 = x0;
      m_y0 = y0;
      
      const int lbuff = transformT::lbuff;
      const int width = transformT::width;

      transformT trans;

      m_grid.resize( 2*(m_maxPix/m_dPix) + 1, 2*(m_maxPix/m_dPix) + 1);
 
      int N = 2*(m_maxPix/m_dPix)+1;

      realT startX, startY;

      cubicConvolTransform<realT> cct;

      for(int i =0; i<N; ++i)
      {
         startX = x0 - m_maxPix + i*m_dPix;
         for(int j=0; j<N;++j)
         {
            startY = y0 - m_maxPix + j*m_dPix;
            
            m_grid(i,j) = 0;

            realT dang = math::rtod(atan(1.0/m_maxRad));
            int Nang = 360./dang + 1;
            dang = 360./Nang;

            #ifdef MX_IMPROC_IMCENTERRADON_OMP
            #pragma omp parallel for
            #endif
            for(int iang =0; iang < Nang; ++iang)
            {
               realT ang = iang*dang;

               eigenImage<realT> kern;
               kern.resize(width,width);

               realT cosang = cos(ang*3.14159/180.);
               realT sinang = sin(ang*3.14159/180.);
               realT lineint = 0;
               for(realT rad=m_minRad; rad < m_maxRad; rad += m_dRad)
               {
                  realT xx = startX + rad * cosang;
                  realT yy = startY + rad * sinang;

                  int i0 = xx;
                  int j0 = yy;
                  realT dx = xx - i0;
                  realT dy = yy - j0;

                  trans(kern, dx, dy);
                  lineint += (im.block(i0-lbuff,j0-lbuff, width, width) * kern).sum();
               }
               #pragma omp critical
               m_grid(i,j) += lineint*lineint*pow(m_dRad,2);
            }
            m_grid(i,j) /= Nang;
         } //j
      } //i

      return fitGrid(x,y, m_grid);
      
   }

   /// Fit the transform grid
   /** A separate function for testing, does not normally need to be called independently
     */
   int fitGrid( realT & x, ///< [out] the x-coordinate of the estimated center of rotational symmetry
                realT & y, ///< [out] the y-coordinate estimated center of rotational symmetry
                eigenImage<realT> grid
              )
   {
      imagePeakInterp<cubicConvolTransform<realT>> ipi(m_peakTol);

      ipi(x,y,grid);

      x = m_x0 - m_maxPix + x*m_dPix;
      y = m_y0 - m_maxPix + y*m_dPix;

      return 0;

      //Get initial guess for the fit.
      realT mx, mn, gx, gy;

      mx = grid.maxCoeff(&gx, &gy);
      mn = grid.minCoeff();


      //------------- Now fit the grid to a 2D elliptical Gaussian ------------
      /** \todo This needs to be able to handle highly non-symmetric peaks much better.
        */
      m_fit.setArray(grid.data(),  grid.rows(), grid.cols());
      m_fit.setGuess(mn, mx-mn, gx, gy, m_guessWidth/m_dPix, m_guessWidth/m_dPix, 0);

      m_fit.fit();

      //The estimated enter of circular symmetry from the fit.
      x = m_x0 - m_maxPix + m_fit.x0()*m_dPix;
      y = m_y0 - m_maxPix + m_fit.y0()*m_dPix;

      return 0;
   }
   
   ///Output the results to a stream
   /** Prints a result summary to the input stream.
     *
     * \tparam iosT is a std::ostream-like type.
     * \tparam comment is a comment character to start each line.  Can be '\0'.
     */ 
   template<typename iosT, char comment='#'>
   iosT & dumpResults( iosT & ios )
   {
      char c[] = {comment, '\0'};
      
      ios << c << "--------------------------------------\n";
      ios << c << "mx::improc::imCenterRadon Results \n";
      ios << c << "--------------------------------------\n";
      ios << c << "Estimated x-center: " << m_x0 - m_maxPix + m_fit.x0()*m_dPix << "\n";
      ios << c << "Estimated y-center: " << m_y0 - m_maxPix + m_fit.y0()*m_dPix << "\n";
      ios << c << "Cost half-width: " << 0.5*math::func::sigma2fwhm<realT>( sqrt( pow(m_fit.sigma_x(),2) + pow(m_fit.sigma_y(),2)) * m_dPix) << " pixels\n";
      ios << c << "--------------------------------------\n";
      ios << c << "Setup:\n";
      ios << c << "--------------------------------------\n";
      ios << c << "maxPix: " << m_maxPix << "\n";
      ios << c << "dPix: " << m_dPix << "\n";
      ios << c << "minRad: " << m_minRad << "\n";
      ios << c << "maxRad: " << m_maxRad << "\n";
      ios << c << "dRad: " << m_dRad << "\n";
      ios << c << "angTolPix: " << m_angTolPix << "\n";
      ios << c << "guessWidth: " << m_guessWidth << "\n";
      ios << c << "--------------------------------------\n";
      ios << c << "Fit results:\n"   ;
      m_fit.dumpReport(ios);
      ios << c << "--------------------------------------\n";
      
      return ios;
   }
   
   ///Output the results to std::cout
   /** Prints a result summary to std::cout
     *
     * \tparam comment is a comment character to start each line.  Can be '\0'.
     */
   template<char comment ='#'>
   std::ostream & dumpResults()
   {
      return dumpResults<std::ostream, comment>(std::cout);
   }
   
};
   
} //improc
} //mx 

#endif //__imCenterRadon_hpp__
