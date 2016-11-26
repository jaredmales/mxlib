/** \file fourierModes.hpp
  * \author Jared R. Males
  * \brief Functions for generating 2D Fourier modes
  * \ingroup signal_processing
  *
  */

#ifndef __fourierModes_hpp__
#define __fourierModes_hpp__

#include <vector>

#include <mx/eigenImage.hpp>

namespace mx
{
 
/** \addtogroup signal_processing
  * @{
  */

///Populate a single cosine mode
/** Makes a 2D image of the cosine mode:
    \f[ 
    M_c(\vec{q}) = \cos\left( 2\pi\frac{m}{D}u + 2\pi\frac{n}{D}v \right) 
    \f]
  * where
    \f[
    \vec{q} = u \hat{u} + v \hat{v}
    \f]
  * and  \f$ D \f$ is taken to be the number of columns in the image.
  *
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  *
  * \tparam typeN is an Eigen-like reference type.
  */  
template<class typeN>
void make_fourier_mode_c(typeN im, typename typeN::Scalar m, typename typeN::Scalar n)
{
   int dim1 = im.cols();
   int dim2 = im.rows();

   typename typeN::Scalar uc, vc, u, v;

   uc = 0.5*(dim1-1.0);
   vc = 0.5*(dim2-1.0);
   
   for(int i=0;i<dim1; ++i)
   {
      u = i - uc;
      for(int j=0;j<dim2; ++j)
      {
         v = j-vc;
         
         im(i,j) = cos(2.*3.14159/(dim1)*(m*u + n*v));
      }
   }
   
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
  * and  \f$ D \f$ is taken to be the number of columns in the image.
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  *
  * \tparam typeN is an Eigen-like reference type.
  */  
template<class typeN>
void make_fourier_mode_s(typeN  im, typename typeN::Scalar m, typename typeN::Scalar n)
{
   int dim1 = im.cols();
   int dim2 = im.rows();

   typename typeN::Scalar uc, vc, u, v;

   uc = 0.5*(dim1-1.0);
   vc = 0.5*(dim2-1.0);
   
   for(int i=0;i<dim1; ++i)
   {
      u = i - uc;
      for(int j=0;j<dim2; ++j)
      {
         v = j-vc;
         
         im(i,j) = sin(2.*3.14159/(dim1)*(m*u + n*v));
      }
   }
   
}

///Populate a single modified Fourier mode
/** Makes a 2D image of the modified fourier mode:
    \f[ 
    M_p(\vec{q}) = \cos\left( 2\pi \frac{m}{D}u + 2\pi\frac{n}{D}v \right) + p  \sin\left( 2\pi\frac{m}{D}u + 2\pi\frac{n}{D}v \right)
    \f]
  * where \f$ p = \pm 1 \f$, \f$ \vec{q} = u \hat{u} + v \hat{v} \f$,
  * and  \f$ D \f$ is taken to be the number of columns in the image.
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  * \param [in] p is +/- 1 specifying which modified Fourier mode
  * 
  * \tparam typeN is an Eigen-like reference type.
  */  
template<class typeN>
void make_modified_fourier_mode(typeN im, typename typeN::Scalar m, typename typeN::Scalar n, int p)
{
   int dim1 = im.cols();
   int dim2 = im.rows();

   typename typeN::Scalar xc, yc, x, y, arg;

   xc = 0.5*(dim1-1.0);
   yc = 0.5*(dim2-1.0);
   
   
   for(int i=0;i<dim1; ++i)
   {
      x = i - xc;
      for(int j=0;j<dim2; ++j)
      {
         y = j-yc;
         
         arg = 2.*3.14159/(dim1)*(m*x + n*y);
         
         im(i,j) = cos(arg) + p*sin(arg);
      }
   }
   
}

///Populate a single modified Fourier +1 mode
/** See \ref make_modified_fourier_mode.
  * 
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  * 
  * \tparam typeN is an Eigen-like reference type.
  */
template<class typeN>
void make_modified_fourier_mode_p(typeN im, typename typeN::Scalar m, typename typeN::Scalar n)
{
   make_modified_fourier_mode(im, m, n, 1);
   
}
 
///Populate a single modified Fourier -1 mode
/** See \ref make_modified_fourier_mode.
  * 
  * \param [out] im is an Eigen-like image
  * \param [in] m specifies the spatial frequency in the u direction
  * \param [in] n specifies the spatial frequency in the v direction
  * 
  * \tparam typeN is an Eigen-like reference type.
  */
template<class typeN>
void make_modified_fourier_mode_m(typeN im, typename typeN::Scalar m, typename typeN::Scalar n)
{
   make_modified_fourier_mode(im, m, n, -1);
}


template<typename floatT>
struct spfreq
{
   floatT k;
   floatT l;
   int cs;
   
   spfreq()
   {
      k = 0;
      l = 0;
      cs = 0;
   }
};

template<typename floatT>
bool comp_spfreq_PandV (const spfreq<floatT> & spf1, const spfreq<floatT> & spf2) 
{ 
   floatT fr1 = spf1.k*spf1.k + spf1.l*spf1.l;
   floatT fr2 = spf2.k*spf2.k + spf2.l*spf2.l;
   
   if(fr1 == fr2)
   {
      if(spf1.l == spf2.l) 
      {
         if(spf1.k == spf2.k) return (spf1.cs < spf2.cs);
         else   return (spf1.k < spf2.k);
      }
      return (spf1.l < spf1.l);
   }
   
   //Otherwise sort by lowest absolute value
   return (fr1 < fr2); 
}


///Fill in a vector of x and y spatial frequencies for a Fourier modal basis
/** Follows the conventions of Poyneer & Veran, JOSAA, 22:8, 1515, (2005)
  * 
  * \param ks [output] is a vector to fill in with the x spatial frequencies
  * \param ls [output] is a vector to fill in with the y spatial frequencies
  * \param cs [output] is a vector to fill in with 0 if only a cosine mode, 1 if both cosine and sine modes for this pair
  * \param N [input] the 1-D sampling, e.g. these modes sample NxN
  */
template<typename floatT>
void make_fourier_mode_freqs_PandV(std::vector<floatT> & ks, std::vector<floatT> & ls, std::vector<bool> & cs, int N)
{
   int Nmodes = N*N;
   
   ks.resize(Nmodes);
   ls.resize(Nmodes);
   cs.resize(Nmodes);
   
   floatT _k, _l;
   int modeCount = 0;
   for(int ll=0; ll<=N; ++ll)
   {
      for(int l=ll, k=0; l>=0; --l, ++k)
      {
         if(k==0 && l > .5*N) continue;
         
         if(k > .5*N && (l == 0 || l>= N-k)) continue;

         if(k > 0.5*N) _k = k-N;
         else _k = k;
         
         if(l > 0.5*N) _l = l-N;
         else _l = l;
         
         //There is always a cosine mode
         ks[modeCount] = _k;
         ls[modeCount] = _l;        
         cs[modeCount] = 0;
         
         ++modeCount;
        
         //Make sure it isn't one of the corners, then populate the sine mode
         if(! ((k==0 && l==0) || (k ==0.5*N && l==0) || (k==0 && l==0.5*N) || (k==0.5*N && l==0.5*N)))
         {
            ks[modeCount] = _k;
            ls[modeCount] = _l;
            cs[modeCount] = 1;
            ++modeCount;

         }
         
      }
   }
   
   std::cout << "modeCount = " << modeCount << "\n";
   
   std::vector<spfreq<floatT>> freqs;
   freqs.resize(ks.size());
   
   for(int i=0;i<ks.size();++i)
   {
      freqs[i].k = ks[i];
      freqs[i].l = ls[i];
      freqs[i].cs = cs[i];
   }
   
   std::sort(freqs.begin(), freqs.end(), comp_spfreq_PandV<floatT>);
   
   for(int i=0;i<freqs.size(); ++i)
   {
      ks[i] = freqs[i].k;
      ls[i] = freqs[i].l;
      cs[i] = freqs[i].cs;
   }
   
}
   

template<typename floatT>
bool comp_spfreq (const spfreq<floatT> & spf1, const spfreq<floatT> & spf2) 
{ 
   floatT fr1 = spf1.k*spf1.k + spf1.l*spf1.l;
   floatT fr2 = spf2.k*spf2.k + spf2.l*spf2.l;
   
   if(fr1 == fr2)
   {
      if(spf1.k == spf2.k) 
      {
         if(spf1.l == spf2.l) return (spf1.cs < spf2.cs);
         else   return (spf1.l < spf2.l);
      }
      return (spf1.l < spf1.l);
   }
   
   //Otherwise sort by lowest absolute value
   return (fr1 < fr2); 
}

template<typename floatT>
void make_fourier_mode_freqs_Circular(std::vector<floatT> & ks, std::vector<floatT> & ls, int N)
{
   floatT frad;
   int Nmodes;

   frad = floor( (floatT) 0.5*N);
   
   Nmodes = 0.25*3.14159*N*N*1.1;//The 1.1 is a buffer

   ks.resize(Nmodes);
   ls.resize(Nmodes);
   
   floatT _k, _l;
   floatT kmax;
   
   int modeCount = 0;
   
   for(int l=0; l<=frad; ++l)
   {
      kmax = sqrt(frad*frad-l*l);
      
      for(int k=( (int) -kmax); k <= kmax; ++k)
      {
         if(l==0 && k <=0) continue;
                  
         ks[modeCount] = k;
         ls[modeCount] = l;        
         ++modeCount;
      }
   }
   
   ks.erase(ks.begin()+modeCount, ks.end());
   ls.erase(ls.begin()+modeCount, ls.end());
   
   std::vector<spfreq<floatT>> freqs;
   freqs.resize(ks.size());
   
   for(int i=0;i<ks.size();++i)
   {
      freqs[i].k = ks[i];
      freqs[i].l = ls[i];
      freqs[i].cs = 1;
   }
   
   std::sort(freqs.begin(), freqs.end(), comp_spfreq<floatT>);
   
   for(int i=0;i<freqs.size(); ++i)
   {
      ks[i] = freqs[i].k;
      ls[i] = freqs[i].l;
   }
}
 
template<typename floatT>
void make_fourier_mode_freqs_Rect(std::vector<floatT> & ks, std::vector<floatT> & ls, int N)
{
   int Nmodes;

   
   Nmodes = 0.5*(N+1)*(N+1);

   ks.resize(Nmodes);
   ls.resize(Nmodes);
   
   floatT _k, _l;
   floatT kmax;
   
   int Ndx = floor(0.5*(N));
   
   int modeCount = 0;
   
   for(int l=-Ndx; l <= Ndx; ++l)
   {      
      for(int k= 0; k <= Ndx; ++k)
      {
         if( k==0 && l <=0 ) continue;
       
         //if(abs(l) == Ndx && abs(k) >= 0.5*Ndx) continue;
         //if(abs(k) == Ndx && abs(l) >= 0.5*Ndx) continue;
            
         ks[modeCount] = k;
         ls[modeCount] = l;        
         ++modeCount;
      }
   }
   
   ks.erase(ks.begin()+modeCount, ks.end());
   ls.erase(ls.begin()+modeCount, ls.end());
   
   std::vector<spfreq<floatT>> freqs;
   freqs.resize(ks.size());
   
   for(int i=0;i<ks.size();++i)
   {
      freqs[i].k = ks[i];
      freqs[i].l = ls[i];
      freqs[i].cs = 1;
   }
   
   std::sort(freqs.begin(), freqs.end(), comp_spfreq<floatT>);
   
   for(int i=0;i<freqs.size(); ++i)
   {
      ks[i] = freqs[i].k;
      ls[i] = freqs[i].l;
   }
} 
 
template<typename eigenT, typename floatT>
void make_freq_map(eigenT & map, std::vector<floatT> & ks, std::vector<floatT> & ls, int N)
{
   map.setZero(N+1, N+1);
   floatT mid = 0.5*(N);
   
   for(int i=0; i<ks.size(); ++i)
   {
      map( mid+ks[i], mid+ls[i]) = sqrt(pow(ks[i],2) + pow(ls[i],2));
   }
}

   
template<typename cubeT>
void make_fourier_basis_PandV(cubeT & cube, int dim, int N)
{
   std::vector<float> ks, ls;
   std::vector<bool> cs;
   
   int Nmodes = N*N;
   
   make_fourier_mode_freqs_PandV(ks, ls, cs, N);
   
   cube.resize(dim,dim,Nmodes-1);
   
   
   for(int n=1; n<Nmodes; ++n)
   {
      if(!cs[n])
         make_fourier_mode_c(cube.image(n-1), ks[n], ls[n]);
      else
         make_fourier_mode_s(cube.image(n-1), ks[n], ls[n]);
   }
}

template<typename cubeT>
void make_fourier_basis_Circular(cubeT & cube, int dim, int N)
{
   typedef typename cubeT::Scalar floatT;
   
   std::vector<floatT> ks, ls; 
   
   make_fourier_mode_freqs_Circular<floatT>(ks, ls, N);

   int Nmodes = 2*ks.size();
   
   cube.resize(dim,dim,Nmodes);
   
   
   for(int n=0; n < ks.size(); ++n)
   {
      make_fourier_mode_c(cube.image(2*n), ks[n], ls[n]);
      
      make_fourier_mode_s(cube.image(2*n+1), ks[n], ls[n]);
   }
}


template<typename cubeT>
void make_fourier_basis_Rect(cubeT & cube, int dim, int N)
{
   typedef typename cubeT::Scalar floatT;
   
   std::vector<floatT> ks, ls; 
   
   
   make_fourier_mode_freqs_Rect<floatT>(ks, ls, N);

   
   int Nmodes = 2*ks.size();
   
   cube.resize(dim,dim,Nmodes);
   
   
   for(int n=0; n < ks.size(); ++n)
   {
      make_fourier_mode_c(cube.image(2*n), ks[n], ls[n]);
      
      make_fourier_mode_s(cube.image(2*n+1), ks[n], ls[n]);
   }
}

template<typename cubeT>
void make_modified_fourier_basis_Rect(cubeT & cube, int dim, int N)
{
   typedef typename cubeT::Scalar floatT;
   
   std::vector<floatT> ks, ls; 
      
   make_fourier_mode_freqs_Rect<floatT>(ks, ls, N);

   int Nmodes = 2*ks.size();
   
   cube.resize(dim,dim,Nmodes);
   
   
   for(int n=0; n < ks.size(); ++n)
   {
      make_modified_fourier_mode_p(cube.image(2*n), ks[n], ls[n]);
      
      make_modified_fourier_mode_m(cube.image(2*n+1), ks[n], ls[n]);
   }
}

///@}

} //namespace mx

#endif //__fourierModes_hpp__
