/** \file imCenterCircleSym.hpp
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

#ifndef __imCenterCircleSym_hpp__
#define __imCenterCircleSym_hpp__

#include "../math/vectorUtils.hpp"
#include "../math/fit/fitGaussian.hpp"
#include "imageFilters.hpp"
#include "eigenImage.hpp"

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
template<typename realT>
struct imCenterCircleSym
{
   
   eigenImage<realT> * mask; ///< An optional pointer to a mask.  If not NULL and the same size as the input image, any 0 pixels in this image are ignored.
   
   eigenImage<realT> grid; ///< The cost grid.
   eigenImage<realT> smGrid; ///< The smoothed cost grid.
   
   
   realT maxPix; ///<The maximum range of the grid search, the grid will be +/- maxPix.  Default: 5.0.
   realT dPix; ///< The spacing of the points in the grid search.  Default: 0.1.
   
   realT minRad; ///< The minimum radius to include in the cost function.   Default: 1.0.
   realT maxRad; ///< The maximum radius to include in the cost function  Default: 20.0.
   realT dRad; ///< The spacing of the radii.   Default: 1.0.
   
   realT smWidth; ///< FWHM of the Guassian smoothing kernel.  Default is 10 (seems to work well regardless of dPix). Set to 0 for no smoothing.
   realT guessWidth; ///< Guess in original pixels for the FWHM of the cost grid.  Default is 0.5, which seems to work well in most cases.
   
   ///The fitter for finding the centroid of the grid.  You can use this to inspect the details of the fit if desired.
   //mx::math::fit::fitGaussian2D<mx::gaussian2D_gen_fitter<realT>> fit; 
   mx::math::fit::fitGaussian2Dgen<realT> fit;
   
   realT _x0; ///< Internal storage of the guess center position
   realT _y0; ///< Internal storage of the guess center position.
   
   imCenterCircleSym()
   {

      mask = 0;
      
      init( 5, 0.1, 1.0, 20.0, 1.0);
            
      smWidth = 10;
      guessWidth = 0.5;
            
   }
   
   ///Initialize the grid.
   int init( realT maxP, ///< [in] the new value of maxPix
             realT dP, ///< [in] the new value of dPix
             realT minR, ///< [in] the new value of minRad
             realT maxR, ///< [in] the new value of maxRad
             realT dR ///< [in] the new value of dRad
           )
   {
      maxPix = maxP;
      dPix = dP;
      minRad = minR;
      maxRad = maxR;
      dRad = dR;
      return 0;
   }

   ///Peform the grid search and fit.
   int center( realT & x, ///< [out] the x-coordinate of the estimated center of rotational symmetry
               realT & y, ///< [out] the y-coordinate estimated center of rotational symmetry
               eigenImage<realT> & im,  ///< [in] the image to centroid
               realT x0,  ///< [in] an initial guess at the x-coordinate, the grid is centered about this point.
               realT y0  ///< [in] an initial guess at the y-coordinate, the grid is centered about this point.
             )
   {
      _x0 = x0;
      _y0 = y0;
      
      bool doMask = false;
      if( mask )
      {
         if( mask->rows() == im.rows() && mask->cols() == im.cols()) doMask = true;
      }
      
      grid.resize( 2*(maxPix/dPix) + 1, 2*(maxPix/dPix) + 1);
 
#ifdef ICCS_OMP      
      #pragma omp parallel
      {
#endif
         int nRad = (maxRad - minRad)/dRad+1;
         std::vector< std::vector<realT> > rads;
         rads.resize( nRad );
      
         int N = 2*(maxPix/dPix)+1;
      
         realT startX, startY;
#ifdef ICCS_OMP
         #pragma omp for
#endif
         for(int i =0; i<N; ++i)
         {
            startX = x0 - maxPix + i*dPix;
            for(int j=0; j<N;++j)
            {
               startY = y0 - maxPix + j*dPix;
            
               for(int k =0; k< rads.size(); ++k)
               {
                  rads[k].clear();
               }
            
               for( int ii = startX - maxRad; ii <= startX+maxRad; ++ii)
               {      
                  for( int jj = startY - maxRad; jj <= startY+maxRad; ++jj)   
                  {
                  
                     if(doMask)
                     {
                        if( (*mask)(ii,jj) == 0) continue;
                     }
                  
                     realT rad;
                     rad = sqrt(pow( ii-startX, 2) + pow( jj-startY,2));
                  
                     if(rad < minRad || rad >= maxRad) continue;
               
                     rads[ (rad-minRad)/dRad ].push_back( im(ii,jj) );
                  }//jj
               }//ii

               grid(i,j) = 0;
            
               realT mean, var;
               for(int k =0; k< rads.size(); ++k)
               {
                  if(rads[k].size() <= 1) continue;
                  mean = math::vectorMean(rads[k]);
                  var = math::vectorVariance(rads[k], mean);
               
                  grid(i,j) += var/pow(mean,2);
               } 
            } //j
         } //i
#ifdef ICCS_OMP
      } //omp parallel
#endif
      //-------- Now Smooth the grid, and multiply by -1 for fitting -------------
      
      if(smWidth > 0)
      {
         filterImage(smGrid, grid, gaussKernel<eigenImage<realT>,2>(10), 0);
      }
      else
      {
         smGrid = grid;
      }
      
      smGrid *= -1;
      
      
      //Get initial guess for the fit.
      realT mx, mn, gx, gy;
      
      mx = smGrid.maxCoeff(&gx, &gy);
      mn = smGrid.minCoeff();
      

      //------------- Now fit the smoothed cost grid to a 2D elliptical Gaussian ------------
      /** \todo This needs to be able to handle highly non-symmetric peaks much better.
        */
      fit.setArray(smGrid.data(), smGrid.rows(), smGrid.cols());
      fit.setGuess(mn, mx-mn, gx, gy, guessWidth/dPix, guessWidth/dPix, 0);
         
      fit.fit();
      
      //The estimated enter of circular symmetry from the fit.
      x = x0 - maxPix + fit.x0()*dPix;
      y = y0 - maxPix + fit.y0()*dPix;
      
     
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
      ios << c << "mx::improc::imCenterCircleSym Results \n";
      ios << c << "--------------------------------------\n";
      ios << c << "Estimated x-center: " << _x0 - maxPix + fit.x0()*dPix << "\n";
      ios << c << "Estimated y-center: " << _y0 - maxPix + fit.y0()*dPix << "\n";
      ios << c << "Cost half-width: " << 0.5*sigma2fwhm( sqrt( pow(fit.sigma_x(),2) + pow(fit.sigma_y(),2)) * dPix) << " pixels\n";
      ios << c << "--------------------------------------\n";
      ios << c << "Setup:\n";
      ios << c << "--------------------------------------\n";
      ios << c << "maxPix: " << maxPix << "\n";
      ios << c << "dPix: " << dPix << "\n";
      ios << c << "minRad: " << minRad << "\n";
      ios << c << "maxRad: " << maxRad << "\n";
      ios << c << "dRad: " << dRad << "\n"; 
      ios << c << "smWidth: " << smWidth << "\n";
      ios << c << "guessWidth: " << guessWidth << "\n";
      ios << c << "--------------------------------------\n";
      ios << c << "Fit results:\n"   ;
      fit.dumpReport(ios);
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

#endif //__imCenterCircleSym_hpp__
