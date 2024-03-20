/** \file fourierBasis.hpp
  * \author Jared R. Males (jaredmales@gmail.com)
  * \brief Generating a fourier basis.
  * \ingroup mxAO_files
  * 
  */

#ifndef __fourierBasis_hpp__
#define __fourierBasis_hpp__

#include "../sigproc/fourierModes.hpp"
#include "../sigproc/zernike.hpp"
#include "../improc/eigenCube.hpp"
#include "../ioutils/fits/fitsFile.hpp"

#include "aoPaths.hpp"

namespace mx
{
namespace AO
{
   
///Make the modified Fourier basis
/** 
  * \param [in] basisName the name of the basis (not including the mx::AO path)
  * \param [in] dim the linear size of the maps, that is they will be dimxdim in size.
  * \param [in] N is the number of degrees of freedom.  Number of modes will be (N+1)(N+1) - 1.
  *
  * \tparam realT the real numeric type for calculations
  */  
template<typename realT>
void makeModfBasis( const std::string & basisName,
                    int dim,
                    int N,
                    realT ang,
                    int nZern = 0
                  )
{
    improc::eigenCube<realT> modes;

    sigproc::makeFourierBasis_Rect( modes, dim, N, MX_FOURIER_MODIFIED, ang);

    improc::eigenCube<realT> * fmodes;
    bool fmodes_allocated;

    if(nZern > 0)
    {
        improc::eigenCube<realT> zModes;
        
        zModes.resize(dim, dim, nZern);
        sigproc::zernikeBasis<improc::eigenCube<realT>, double>( zModes);

        fmodes = new improc::eigenCube<realT>;
        fmodes_allocated = true;
        fmodes->resize(dim,dim,modes.planes()+nZern);

        for(int p=0; p < nZern; ++p)
        {
            fmodes->image(p) = zModes.image(p);
        }

        for(int p=0; p < modes.planes(); ++p)
        {
            fmodes->image(p+nZern) = modes.image(p);
        }

        modes.resize(0,0,0);
        zModes.resize(0,0,0);
    }
    else
    {
        fmodes = &modes;
        fmodes_allocated = false;
    }

    fits::fitsFile<realT> ff;
   
    std::string fName = mx::AO::path::basis::modes(basisName, true);
      
    ff.write(fName, *fmodes);

    if(fmodes_allocated)
    {
        delete fmodes;
    }
}

} //namespace AO
} //namespace mx

   
#endif //__fourierBasis_hpp__
