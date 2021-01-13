/** \file fourierModes.hpp
  * \author Jared R. Males
  * \brief Functions for generating 2D Fourier modes
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

#ifndef fourierModes_hpp
#define fourierModes_hpp

#include <vector>
#include <algorithm>


#include "../mxError.hpp"

#include "../math/constants.hpp"
#include "../math/geo.hpp"

namespace mx
{

namespace sigproc
{

/**
  * \ingroup fourier_basis
  *
  * \todo need to handle references and r-values so this works directly with Eigen arrays.
  *
  * @{
  */

/** \def MX_FOURIER_BASIC
  * \brief Signifies the basic sine/cosine Fourier basis.
  */
#define MX_FOURIER_BASIC 0

/** \def MX_FOURIER_MODIFIED
  * \brief Signifies the modified Fourier basis.
  */
#define MX_FOURIER_MODIFIED 1

///Populate a single cosine mode
/** Makes a 2D image of the cosine mode:
    \f[
    M_c(\vec{q}) = \cos\left( 2\pi\frac{m}{D}u + 2\pi\frac{n}{D}v \right)
    \f]
  * where
    \f[
    \vec{q} = u \hat{u} + v \hat{v}
    \f]
  * and  \f$ D \f$ is taken to be the maximum of the number of columns and rows in the image.
  *
  * \retval 0 on success
  *
  * \tparam typeN is an Eigen-like reference type.
  */
template<class typeN>
int makeFourierModeC( typeN im,                   ///<  [out] is an Eigen-like image
                      typename typeN::Scalar m,   ///<  [in] specifies the spatial frequency in the u direction
                      typename typeN::Scalar n    ///<  [in] n specifies the spatial frequency in the v direction
                    )
{
   typedef typename typeN::Scalar realT;

   int dim1 = im.cols();
   int dim2 = im.rows();

   realT uc, vc, u, v, D;

   uc = 0.5*(dim1-1.0);
   vc = 0.5*(dim2-1.0);

   D = std::max( dim1, dim2);

   for(int i=0;i<dim1; ++i)
   {
      u = i - uc;
      for(int j=0;j<dim2; ++j)
      {
         v = j-vc;

         im(i,j) = cos(math::two_pi<realT>()/D*(m*u + n*v));
      }
   }

   return 0;
}

///Populate a single sine mode
/** Makes a 2D image of the cosine mode:
    \f[
    M_s(\vec{q}) = \sin\left( 2\pi\frac{m}{D}u + 2\pi\frac{n}{D}v \right)
    \f]
  * where
    \f[
    \vec{q} = u \hat{u} + v \hat{v}
    \f]
  * and  \f$ D \f$ is taken to be the maximum of the number of columns and rows in the image.
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  *
  * \retval 0 on success
  *
  * \tparam typeN is an Eigen-like reference type.
  */
template<class typeN>
int makeFourierModeS(typeN  im, typename typeN::Scalar m, typename typeN::Scalar n)
{
   typedef typename typeN::Scalar realT;

   int dim1 = im.cols();
   int dim2 = im.rows();

   realT uc, vc, u, v, D;

   uc = 0.5*(dim1-1.0);
   vc = 0.5*(dim2-1.0);

   D = std::max(dim1, dim2);

   for(int i=0;i<dim1; ++i)
   {
      u = i - uc;
      for(int j=0;j<dim2; ++j)
      {
         v = j-vc;

         im(i,j) = sin(math::two_pi<realT>()/D*(m*u + n*v));
      }
   }

   return 0;
}

///Populate a single basic Fourier mode
/** Makes a 2D image of the basic sine or cosine mode.  Calls \ref makeFourierModeC and \ref makeFourierModeS.
  *
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  * \param [in] p species sine (-1) or cosine (+1)
  *
  * \retval 0 on success
  *
  * \tparam typeN is an Eigen-like reference type.
  */
template<class typeN>
int makeFourierMode(typeN  im, typename typeN::Scalar m, typename typeN::Scalar n, int p)
{
   if( p == -1) return makeFourierModeS(im, m, n);

   if( p == 1) return makeFourierModeC(im, m, n);

   mxError("makeFourierMode", MXE_INVALIDARG, "p must be +1 (cosine) or -1 (sine).");

   return -1;
}

///Populate a single modified Fourier mode
/** Makes a 2D image of the modified fourier mode:
    \f[
    M_p(\vec{q}) = \cos\left( 2\pi \frac{m}{D}u + 2\pi\frac{n}{D}v \right) + p  \sin\left( 2\pi\frac{m}{D}u + 2\pi\frac{n}{D}v \right)
    \f]
  * where \f$ p = \pm 1 \f$, \f$ \vec{q} = u \hat{u} + v \hat{v} \f$,
  * and  \f$ D \f$ is taken to be the maximum of the number of columns in the image.
  *
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  * \param [in] p is +/- 1 specifying which modified Fourier mode
  *
  * \retval 0 on success
  *
  * \tparam typeN is an Eigen-like reference type.
  */
template<class typeN>
int makeModifiedFourierMode( typeN im,
                             typename typeN::Scalar m,
                             typename typeN::Scalar n,
                             int p,
                             typename typeN::Scalar ang = 0
                           )
{
   typedef typename typeN::Scalar realT;

   if(p != 1 && p != -1)
   {
      mxError("makeModifiedFourierMode", MXE_INVALIDARG, "p must be +1 or -1.");
   }

   int dim1 = im.cols();
   int dim2 = im.rows();

   realT xc, yc, x, y, arg, D;

   xc = 0.5*(dim1-1.0);
   yc = 0.5*(dim2-1.0);

   D = std::max(dim1, dim2);

   realT c_ang, s_ang;
   if(ang != 0)
   {
      c_ang = cos( math::dtor(ang));
      s_ang = sin( math::dtor(ang));
   }

   for(int i=0;i<dim1; ++i)
   {

      for(int j=0;j<dim2; ++j)
      {
         x = i - xc;
         y = j - yc;

         if(ang != 0)
         {
            realT x0 = x*c_ang - y*s_ang;
            realT y0 = x*s_ang + y*c_ang;

            x=x0;
            y=y0;
         }

         arg = math::two_pi<realT>()/D*(m*x + n*y);

         im(i,j) = cos(arg) + p*sin(arg);
      }
   }

   return 0;

}

///Populate a single modified Fourier +1 mode
/** See \ref makeModifiedFourierMode.
  *
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  *
  * \tparam typeN is an Eigen-like reference type.
  */
template<class typeN>
int makeModifiedFourierModeP( typeN im,
                              typename typeN::Scalar m,
                              typename typeN::Scalar n,
                              typename typeN::Scalar ang = 0
                            )
{
   return makeModifiedFourierMode(im, m, n, 1, ang);
}

///Populate a single modified Fourier -1 mode
/** See \ref makeModifiedFourierMode.
  *
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  *
  * \tparam typeN is an Eigen-like reference type.
  */
template<class typeN>
int makeModifiedFourierModeM( typeN im,
                              typename typeN::Scalar m,
                              typename typeN::Scalar n,
                              typename typeN::Scalar ang = 0
                            )
{
   return makeModifiedFourierMode(im, m, n, -1, ang);
}


///Structure to contain the parameters of a Fourier mode.
struct fourierModeDef
{
   int m; ///< Spatial frequency m index
   int n; ///< Spatial frequency n index
   int p; ///< +/- 1 for modified Fourier modes, -1==>sine, +1==>cosine for basic Fourier modes.

   ///Constructor
   fourierModeDef()
   {
      m = 0;
      n = 0;
      p = 0;
   }
};


/// Sorting functor for modes according to the Poyner and Veran (2005) standard.
/** Follows the conventions of  Poyneer & Veran (2005) \cite poyneer_and_veran_2005
  *
  * \param[in] spf1 is an object of type \ref fourierModeDef to compare with spf2
  * \param[in] spf2 is an object of type \ref fourierModeDef to compare with spf1
  *
  * \returns the result of (spf1 < spf2) according to the above rules.
  */
bool comp_fourierModeDef_PandV (const fourierModeDef & spf1, const fourierModeDef & spf2)
{
   double fr1 = spf1.m*spf1.m + spf1.n*spf1.n;
   double fr2 = spf2.m*spf2.m + spf2.n*spf2.n;

   if(fr1 == fr2)
   {
      if(spf1.n == spf2.n)
      {
         if(spf1.m == spf2.m) return (spf1.p < spf2.p);
         else   return (spf1.m < spf2.m);
      }
      return (spf1.n < spf2.n);
   }

   //Otherwise sort by lowest absolute value
   return (fr1 < fr2);
}


///Generate a Poyneer and Veran spatial frequency grid.
/** Follows the convention of Poyneer and Veran (2005) \cite poyneer_and_veran_2005
  * \todo use push_back instead of resize
  * \param spf
  * \param N is the linear number of degrees of freedom.
  */
int makeFourierModeFreqs_PandV(std::vector<fourierModeDef> & spf, int N)
{
   int Nmodes = N*N;

   spf.resize(Nmodes);

   int _k, _l;
   int modeCount = 0;
   for(int ll=0; ll<=N; ++ll)
   {
      for(int l=ll, k=0; l>=0; --l, ++k)
      {
         if(modeCount >= Nmodes)
         {
            mxError("makeFourierModeFreqs_PandV", MXE_SIZEERR, "mode count exceeded expected number of modes");
            return -1;
         }


         if(k==0 && l > .5*N) continue;

         if(k > .5*N && (l == 0 || l>= N-k)) continue;

         if(k > 0.5*N) _k = k-N;
         else _k = k;

         if(l > 0.5*N) _l = l-N;
         else _l = l;

         //There is always a cosine mode
         spf[modeCount].m = _k;
         spf[modeCount].n = _l;
         spf[modeCount].p = 1;

         ++modeCount;

         //Make sure it isn't one of the corners, then populate the sine mode
         if(! ((k==0 && l==0) || (k ==0.5*N && l==0) || (k==0 && l==0.5*N) || (k==0.5*N && l==0.5*N)))
         {
            spf[modeCount].m = _k;
            spf[modeCount].n = _l;
            spf[modeCount].p = -1;
            ++modeCount;

         }
      }
   }

   std::sort(spf.begin(), spf.end(), comp_fourierModeDef_PandV);

   return 0;
}

/// Sorting functor for modes according to the mx::AO standard.
/** Given a spatial frequency vector \f$ \vec{k} = k_u \hat{u} + k_v \hat{v} \f$, sorts first by
  * \f$ | \vec{k} | \f$, then by the angle from the u axis.  Results in mode numbers increasing radially
  * and counter-clockwise from the \f$ u \f$ axis.
  *
  * \param[in] spf1 is an object of type \ref fourierModeDef to compare with spf2
  * \param[in] spf2 is an object of type \ref fourierModeDef to compare with spf1
  *
  * \returns the result of (spf1 < spf2) according to the above rules.
  */
bool comp_fourierModeDef (const fourierModeDef & spf1, const fourierModeDef & spf2)
{
   double k1 = spf1.m*spf1.m + spf1.n*spf1.n;
   double k2 = spf2.m*spf2.m + spf2.n*spf2.n;

   if(k1 == k2)
   {
      //Now compare by angle
      double a1 = atan2(spf1.n, spf1.m);
      if(a1 < 0) a1 += math::two_pi<double>();

      double a2 = atan2(spf2.n, spf2.m);
      if(a2 < 0) a2 += math::two_pi<double>();

      if( a1 == a2 ) return (spf1.p < spf2.p);

      return ( a1 < a2);

   }

   //Otherwise sort by lowest absolute value
   return (k1 < k2);
}


///Generate a circular spatial frequency grid.
/** \todo use push_back instead of resize
  *
  * \param spf
  * \param N is the linear number of degrees of freedom.  The number of modes is ((N+1)(N+1) - 1).
  */
int makeFourierModeFreqs_Circ( std::vector<fourierModeDef> & spf,
                               int N )
{
   int krad;
   int Nmodes;

   krad = floor(0.5*N);

   Nmodes = 0.25*math::pi<double>()*N*N*1.1;//The 1.1 is a buffer

   spf.resize(Nmodes);

   int modeCount = 0;

   for(int m=-krad; m <= krad; ++m)
   {
      int nmax = sqrt(krad*krad-m*m);

      for(int n=0; n <= nmax; ++n)
      {
         if( n==0 && m <=0 ) continue;

         for(int p=-1; p<=1; p+=2)
         {
            if(modeCount >= Nmodes)
            {
               mxError("makeFourierModeFreqs_Circ", MXE_SIZEERR, "mode count exceeded expected number of modes");
               return -1;
            }

            spf[modeCount].m = m;
            spf[modeCount].n = n;
            spf[modeCount].p = p;
            ++modeCount;
         }
      }
   }

   //Erase any extra modes (there will bedue to the buffer).
   spf.erase(spf.begin()+modeCount, spf.end());

   //And now sort it
   std::sort(spf.begin(), spf.end(), comp_fourierModeDef);

   return 0;
}

///Calculate the index of a Fourier mode given its (m,n,p) coordinates.
/** The index increases counter-clockwise from the m axis in squares, with p==-1 being even and p==1 being odd.
  *
  * \param [in] m is the m spatial frequency coordinate
  * \param [in] n is the n spatial frequency coordinage
  * \param [in] p is the parity, +/-1, for sine or cosine or specifiying the modified Fourier mode.
  *
  * \returns the index of the Fourier mode on success
  * \returns -1 on an error.
  */
int fourierModeNumber(int m, int n, int p)
{

   if( p != -1 && p != 1)
   {
      mxError("fourierModeCoordinates", MXE_INVALIDARG, "p must +/-1");
      return -1;
   }

   if(m ==0 && n == 0)
   {
      mxError("fourierModeCoordinates", MXE_INVALIDARG, "piston is ignored");
      return -1;
   }


   int m0 = std::max(abs(m),abs(n));

   int i0 = 2*(m0*m0 - m0);

   //int di0= 4*m0;

   int di;

   if( m == m0 )
   {
      di = n;
   }
   else if( m == - m0 )
   {
      di = 4*m0 - n;
   }
   else
   {
      di = 2*m0 - m;

   }

   return 2*(i0 + di) + (1 + p)/2;
}


///Calculate the (m,n,p) coordinates of a Fourier mode given its index.
/** The index increases counter-clockwise from the m axis in squares, with p==-1 being even and p==1 being odd.
  *
  * \returns 0 on success, a negative number otherwise
  */
int fourierModeCoordinates( int & m, ///< [out] is the m spatial frequency coordinate
                            int & n, ///< [out] is the n spatial frequency coordinage
                            int & p, ///< [out] is the parity, +/-1, for sine or cosine or specifiying the modified Fourier mode.
                            int i )  ///< [in] is the mode index.
{
   if(i < 0)
   {
      mxError("fourierModeCoordinates", MXE_INVALIDARG, "i must be >= 0");
      m = 0;
      n = 0;
      p = 0;
      return -1; /// \retval -1 if i < 0
   }


   int ii = (i)/2;

   int m0 = 1;

   while( 2*(m0*m0-m0) <= ii) ++m0;
   --m0;

   int di = ii - 2*(m0*m0-m0);

   if( di <= m0 )
   {
      m = m0;
      n = di;
   }
   else if( di <= 3*m0 )
   {
      m = 2*m0 - di;
      n = m0;
   }
   else
   {
      m = -m0;
      n = 4*m0 - di;
   }

   p = -1 + 2*(i%2);

   return 0;
}








///Generate a rectangular spatial frequency grid.
/**
  * \param spf is a vector of fourierModeDef structures, which will be resized and populated.
  * \param N is the linear number of degrees of freedom.  The number of modes is ((N+1)(N+1) - 1).
  */
int makeFourierModeFreqs_Rect( std::vector<fourierModeDef> & spf,
                                int N )
{
   int Ndx = floor(0.5*(N));

   int Nmodes = 2*((2*Ndx + 1)*(Ndx + 1) - (Ndx+1));

   spf.resize(Nmodes);

   int m, n, p;

   for(int i=0; i< Nmodes; ++i)
   {
      fourierModeCoordinates(m, n, p, i);

      spf[i].m = m;
      spf[i].n = n;
      spf[i].p = p;
   }

   return 0;


}


/// Fill in a cube with a Fourier basis.
/** Fills the cube with either the basic or modified Fourier basis for the modes specified.
  *
  * \param[out] cube will be allocated to hold and will be filled with the modes
  * \param[in] dim is the linear size of the maps, each is
  * \param[in] spf is a vector of mode definitions to use for each mode
  * \param[in] basisType is either MX_FOURIER_MODIFIED or MX_FOURIER_BASIC
  *
  * \tparam cubeT is an eigen-like cube, e.g. \ref mx::eigenCube.
  */
template<typename cubeT>
int fillFourierBasis( cubeT & cube,
                      int dim,
                      std::vector<fourierModeDef> spf,
                      int basisType,
                      typename cubeT::Scalar ang = 0
                    )
{
   int Nmodes = spf.size();

   cube.resize(dim,dim,Nmodes);

   for(int j=0; j < Nmodes; ++j)
   {
      if(basisType == MX_FOURIER_MODIFIED)
      {
         if( makeModifiedFourierMode( cube.image(j), spf[j].m, spf[j].n, spf[j].p, ang) < 0) return -1;
      }
      else
      {
         if( makeFourierMode(cube.image(j), spf[j].m, spf[j].n, spf[j].p) < 0) return -1;
      }
   }

   return 0;
}

/// Generate a rectangular modified Fourier basis with the Poyneer and Veran ordering.
/** Fills the cube with modes generated by \ref makeFourierMode or \ref makeModifiedFourierMode using the mode grid generated by \ref makeFourierModeFreqs_PandV.
  * This follows the convention of Poyneer and Veran (2005) \cite poyneer_and_veran_2005 for mode numbering.
  *
  * \param[out] cube will be allocated to hold and will be filled with the modes
  * \param[in] dim is the linear size of the maps, each is
  * \param[in] N is the linear number of degrees of freedom.  The number of modes is ((N+1)(N+1) - 1).
  * \param[in] basisType is either MX_FOURIER_MODIFIED or MX_FOURIER_BASIC
  *
  * \tparam cubeT is an eigen-like cube, e.g. \ref mx::eigenCube.
  */
template<typename cubeT>
int makeFourierBasis_PandV( cubeT & cube,
                            int dim,
                            int N,
                            int basisType )
{
   std::vector<fourierModeDef> spf;

   if( makeFourierModeFreqs_PandV(spf, N) < 0)
   {
      return -1;
   }

   return fillFourierBasis(cube, dim, spf, basisType);
}

/// Generate a circular Fourier basis.
/** Fills the cube with modes generated by \ref makeFourierMode or \ref makeModifiedFourierMode using the mode grid generated by \ref makeFourierModeFreqs_Circ.
  *
  * \param[out] cube will be allocated to hold and will be filled with the modes
  * \param[in] dim is the linear size of the maps, each is
  * \param[in] N is the linear number of degrees of freedom.  The number of modes is ((N+1)(N+1) - 1).
  * \param[in] basisType is either MX_FOURIER_MODIFIED or MX_FOURIER_BASIC
  *
  * \tparam cubeT is an eigen-like cube, e.g. \ref mx::eigenCube.
  */
template<typename cubeT>
int makeFourierBasis_Circ(cubeT & cube,
                          int dim,
                          int N,
                          int basisType )
{
   std::vector<fourierModeDef> spf;

   if( makeFourierModeFreqs_Circ(spf, N) < 0 )
   {
      return -1;
   }

   return fillFourierBasis(cube, dim, spf, basisType);
}


/// Generate a rectangular Fourier basis.
/** Fills the cube with modes generated by \ref makeFourierMode or \ref makeModifiedFourierMode using the mode grid generated by \ref makeFourierModeFreqs_Rect.
  *
  * \param[out] cube will be allocated to hold and will be filled with the modes
  * \param[in] dim is the linear size of the maps, each is
  * \param[in] N is the linear number of degrees of freedom.  The number of modes is ((N+1)(N+1) - 1).
  * \param[in] basisType is either MX_FOURIER_MODIFIED or MX_FOURIER_BASIC
  *
  * \tparam cubeT is an eigen-like cube, e.g. \ref mx::eigenCube.
  */
template<typename cubeT>
int makeFourierBasis_Rect( cubeT & cube,
                           int dim,
                           int N,
                           int basisType,
                           typename cubeT::Scalar ang = 0
                         )
{
   std::vector<fourierModeDef> spf;

   if( makeFourierModeFreqs_Rect(spf, N) < 0 )
   {
      return -1;
   }

   return fillFourierBasis(cube, dim, spf, basisType, ang);
}

///@}

} //namespace sigproc
} //namespace mx

#endif //fourierModes_hpp
