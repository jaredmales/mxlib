/** \file imageXCorrFFT.hpp
 * \brief A class to register images using the Fourier cross correlation with a peak fit.
 * \ingroup image_processing_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2020 Jared R. Males (jaredmales@gmail.com)
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

#ifndef imageXCorrFFT_hpp
#define imageXCorrFFT_hpp

#include "../mxError.hpp"
#include "../mxException.hpp"
#include "../math/fft/fft.hpp"
#include "../math/fit/fitGaussian.hpp"

#include "imageUtils.hpp"
#include "imageTransforms.hpp"
#include "imageXCorr.hpp"

namespace mx
{
namespace improc
{

/// Find the optimum shift to align two images using the FFT cross correlation.
/** The reference image must be the same size as the target image.  Both the reference image and the section of
 * target image being analyzed are mean subtracted and variance normalized.  An optional mask can be supplied,
 * which limits the pixels used for the mean and variance calculation.
 *
 * Typical usage will be to set the mask, then the reference, then repeatedly call operator() to
 * determine the shifts for a sequence of imaages.  No new heap allocations take place on these calls
 * to operator(), and the reference image is not re-normalized on each call.
 *
 * The shift is reported in pixels such that if the mxlib imageShift function is used
 * to shift the input image by the negative of the shifts, it will align with the
 * reference at the center of the array.
 *
 * Three peak finding methods are provided.  xcorrPeakMethod::centerOfLight uses center of light,
 * xcorrPeakMethod::gaussFit uses Gaussian centroiding, and xcorrPeakMethod::interpPeak uses interpolation
 * to find the peak to a given tolerance.  For well behaved images gaussFit is most accurate.  interpPeak
 * should be robust against weird shapes and other problems, but accuracy is somewhat worse.
 *
 * \tparam _ccImT is the Eigen-like array type used for image processing.  See typedefs.
 *
 * \ingroup image_reg
 */
template <class _ccImT>
class imageXCorrFFT
{
  public:
    typedef _ccImT ccImT; ///< the Eigen-like array type used for image processing

    typedef typename _ccImT::Scalar Scalar; ///< the scalar type of the image type

    typedef typename _ccImT::Scalar realT; ///< the scalar type of the image type

    typedef std::complex<realT> complexT; ///< Complex floating point type.

    typedef Eigen::Array<realT, Eigen::Dynamic, Eigen::Dynamic> realArrayT; ///< Real eigen array type

    typedef Eigen::Array<complexT, Eigen::Dynamic, Eigen::Dynamic>
        complexArrayT; ///< Complex eigen array type with Scalar==complexT

  protected:
    int m_rows{ 0 };

    int m_cols{ 0 };

    /** \name Working Memory
     * @{
     */

    ccImT m_refIm; ///< The normalized reference image.

    realT m_refX0{ 0 }; ///< The shift of the reference image to itself, used as coordinate origint

    realT m_refY0{ 0 }; ///< The shift of the reference image to itself, used as coordinate origint

    ccImT m_maskIm; ///< Mask image to use, may be needed for proper normalization even if refIm has 0 mask applied.

    bool m_haveMask{ false }; ///< Flag indicating that a mask has been provided.

    bool m_normalize{ true };

    ccImT m_normIm; ///< The normalized image.

    ccImT m_ccIm; ///< The cross-correlation image

    ccImT m_magIm; ///< The magnified image, used if m_peakMethod == xcorrPeakMethod::interp

    complexArrayT m_ftIm0; ///< Working memory for the FT of the reference image.

    complexArrayT m_ftWork; ///< Working memory for the FFT.

    math::fft::fftT<realT, complexT, 2, 0> m_fft_fwd; ///< FFT object for the forward transform.

    math::fft::fftT<complexT, realT, 2, 0> m_fft_back; ///< FFT object for the backward transfsorm.

    ///@}

    math::fit::fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<Scalar>> m_fitter;

    int m_maxLag{ 5 }; ///< The maximum lag to consider in the initial cross-correlation.  Default is 5.

    Scalar m_tol{ 0.1 }; ///< The tolerance of the interpolated-magnified image, in pixels.

    Scalar m_magSize{ 0 }; ///< Magnified size of the ccIm when using interp.  Set as function of m_tol and m_maxLag.

    xcorrPeakMethod m_peakMethod{ xcorrPeakMethod::centerOfLight };

  public:
    /// Default c'tor
    imageXCorrFFT();

    /// Constructor setting maxLag.
    explicit imageXCorrFFT( int maxLag );

    /// Set the size of the cross-correlation images.
    /** This resizes all working memory and conducts fftw planning.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int resize( int nrows, ///< [in] the number of rows in the images to register
                int ncols  ///< [in] the number of columns in the images to register
    );

    /// Get the current value of the normalize flag
    /**
     * \returns the current value of m_normalize
     */
    bool normalize();

    /// Set the normalize flag
    /**
     * \todo check if this has side effects that should be triggered
     */
    void normalize( bool no /**< [in] the new normalize flag value */ );

    /// Get the current maximum lag
    /**
     * \returns the current value of m_maxLag
     */
    int maxLag();

    /// Set the maximum lag
    void maxLag( int ml /**< [in] the new maximum lag */ );

    /// Get the tolerance of the interpolated-magnified image, in pixels.
    /**
     * \returns the current value of m_tol.
     */
    Scalar tol();

    /// Set the tolerance of the interpolated-magnified image, in pixels.
    void tol( Scalar nt /**< [in] The new value of the interpolation tolerance. */ );

    xcorrPeakMethod peakMethod();

    void peakMethod( xcorrPeakMethod xpm );

    /// Set the mask image
    /**
     * \returns 0 on success
     * \returns -1 on error
     */
    int maskIm( const ccImT &mask /**< [in] the new mask image */ );

    /// Get a reference to the mask image.
    /**
     * \returns a const referance to the mask image.
     */
    const ccImT &maskIm();

    /// Set the reference image
    /** Normalizes the reference image by mean subtraction and variance division.  Applies
     * the mask first if supplied. Then Fourier transform and conjugates the result.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int refIm( const ccImT &im0 );

    /// Get a reference to the reference image.
    /**
     * \returns a const referent to m_refIm.
     */
    const ccImT &refIm();

    /// Get a reference to the normalized image.
    /**
     * \returns a const referent to m_normIm.
     */
    const ccImT &normIm();

    /// Get a reference to the cross correlation image.
    /**
     * \returns a const referent to m_ccIm.
     */
    const ccImT &ccIm();

    /// Get a reference to the magnified image.
    /**
     * \returns a const referent to m_magIm.
     */
    const ccImT &magIm();

  protected:
    void findPeak( Scalar &xShift, ///< [out] the x shift of im w.r.t. im0, in pixels
                   Scalar &yShift  ///< [out] the y shift of im w.r.t. im0, in pixels
    );

  public:
    /// Conduct the cross correlation to a specified tolerance
    /**
     * \returns 0 on success
     * \returns -1 on error
     */
    template <class imT>
    int operator()( Scalar &xShift, ///< [out] the x shift of im w.r.t. im0, in pixels
                    Scalar &yShift, ///< [out] the y shift of im w.r.t. im0, in pixels
                    const imT &im   ///< [in] the image to cross-correlate with the reference
    );

    /// Conduct the cross correlation to a specified tolerance
    /**
     * \overload
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    template <class im0T, class imT>
    int operator()( Scalar &xShift, ///< [out] the x shift of im w.r.t. im0, in pixels
                    Scalar &yShift, ///< [out] the y shift of im w.r.t. im0, in pixels
                    im0T &im0,      ///< [in] a new reference image
                    imT &im         ///< [in] the image to cross-correlate with the reference
    );
};

template <class ccImT>
imageXCorrFFT<ccImT>::imageXCorrFFT()
{
    maxLag( m_maxLag ); // this sets up m_magSize.
}

template <class ccImT>
imageXCorrFFT<ccImT>::imageXCorrFFT( int mL )
{
    maxLag( mL );
}

template <class ccImT>
int imageXCorrFFT<ccImT>::resize( int nrows, int ncols )
{
    if( m_rows == nrows && m_cols == ncols )
    {
        return 0;
    }

    m_rows = nrows;

    m_cols = ncols;

    m_ccIm.resize( m_rows, m_cols );

    m_ftIm0.resize( (int)( 0.5 * m_rows ) + 1, m_cols );

    m_ftWork.resize( (int)( 0.5 * m_rows ) + 1, m_cols );

    // fftw is row-major, eigen defaults to column-major
    m_fft_fwd.plan( m_cols, m_rows, MXFFT_FORWARD, false );

    m_fft_back.plan( m_cols, m_rows, MXFFT_BACKWARD, false );

    return 0;
}

template <class ccImT>
bool imageXCorrFFT<ccImT>::normalize()
{
    return m_normalize;
}

template <class ccImT>
void imageXCorrFFT<ccImT>::normalize( bool no )
{
    m_normalize = no;
}

template <class ccImT>
int imageXCorrFFT<ccImT>::maxLag()
{
    return m_maxLag;
}

template <class ccImT>
void imageXCorrFFT<ccImT>::maxLag( int ml )
{
    m_maxLag = ml;
    tol( m_tol );
}

template <class ccImT>
typename ccImT::Scalar imageXCorrFFT<ccImT>::tol()
{
    return m_tol;
}

template <class ccImT>
void imageXCorrFFT<ccImT>::tol( Scalar nt )
{
    m_magSize = ceil( ( ( 2. * m_maxLag + 1 ) - 1.0 ) / nt ) + 1;

    Scalar mag = ( m_magSize - 1.0 ) / ( ( 2. * m_maxLag + 1 ) - 1.0 );

    m_tol = 1.0 / mag;

    if( m_refIm.rows() != m_rows || m_refIm.cols() != m_cols || m_rows == 0 || m_cols == 0 )
        return;

    // Find the reference shift
    m_refX0 = 0;
    m_refY0 = 0;
    operator()( m_refX0, m_refY0, m_refIm );
}

template <class ccImT>
xcorrPeakMethod imageXCorrFFT<ccImT>::peakMethod()
{
    return m_peakMethod;
}

template <class ccImT>
void imageXCorrFFT<ccImT>::peakMethod( xcorrPeakMethod xpm )
{
    m_peakMethod = xpm;

    if( m_refIm.rows() != m_rows || m_refIm.cols() != m_cols || m_rows == 0 || m_cols == 0 )
        return;

    // Find the reference shift
    m_refX0 = 0;
    m_refY0 = 0;
    operator()( m_refX0, m_refY0, m_refIm );
}

template <class ccImT>
int imageXCorrFFT<ccImT>::maskIm( const ccImT &mask )
{
    m_maskIm = mask;
    m_haveMask = true;

    return 0;
}

template <class ccImT>
const ccImT &imageXCorrFFT<ccImT>::maskIm()
{
    return m_maskIm;
}

template <class ccImT>
int imageXCorrFFT<ccImT>::refIm( const ccImT &im )
{
    ccImT im0;

    // Mask if needed
    if( m_haveMask )
    {
        if( im.rows() != m_maskIm.rows() && im.cols() != m_maskIm.cols() )
        {
            mxError( "imageXCorFit::setReference", MXE_SIZEERR, "reference and mask are not the same size" );
            return -1;
        }
        im0 = im * m_maskIm;
    }
    else
    {
        im0 = im;
    }

    // Setup the FFTW space
    resize( im0.rows(), im0.cols() );

    // Now normalize
    if( m_normalize )
    {
        realT m = imageMean( im0 );
        realT v = imageVariance( im0, m );
        m_refIm = ( im0 - m ) / sqrt( v );
    }
    else
    {
        m_refIm = im0;
    }
    // We save refIm as the un-shifted version

    // Now shift so center pixel is 0,0
    imageShiftWP( im0, m_refIm, 0.5 * m_rows + 1, 0.5 * m_cols + 1 );

    // Then FT
    m_fft_fwd( m_ftIm0.data(), im0.data() );

    // Conjugate and normalize for FFTW scaling.
    for( int c = 0; c < m_ftIm0.cols(); ++c )
    {
        for( int r = 0; r < m_ftIm0.rows(); ++r )
        {
            complexT val = std::conj( m_ftIm0( r, c ) );
            val /= ( m_rows * m_cols );
            m_ftIm0( r, c ) = val;
        }
    }

    // And finally find the reference shift
    m_refX0 = 0;
    m_refY0 = 0;

    operator()( m_refX0, m_refY0, m_refIm );

    return 0;
}

template <class ccImT>
const ccImT &imageXCorrFFT<ccImT>::refIm()
{
    return m_refIm;
}

template <class ccImT>
const ccImT &imageXCorrFFT<ccImT>::normIm()
{
    return m_normIm;
}

template <class ccImT>
const ccImT &imageXCorrFFT<ccImT>::ccIm()
{
    return m_ccIm;
}

template <class ccImT>
const ccImT &imageXCorrFFT<ccImT>::magIm()
{
    return m_magIm;
}

template <class ccImT>
void imageXCorrFFT<ccImT>::findPeak( Scalar &xShift, Scalar &yShift )
{
    int xLag0, yLag0;
    Scalar pk = m_ccIm.maxCoeff( &xLag0, &yLag0 );
    Scalar mn = m_ccIm.minCoeff();

    if( xLag0 - m_maxLag < 0 )
        m_maxLag = xLag0;
    if( xLag0 + 2 * m_maxLag + 1 >= m_ccIm.rows() )
        m_maxLag = ( m_ccIm.rows() - 1 - xLag0 ) / 2;
    if( yLag0 - m_maxLag < 0 )
        m_maxLag = yLag0;
    if( yLag0 + 2 * m_maxLag + 1 >= m_ccIm.cols() )
        m_maxLag = ( m_ccIm.cols() - 1 - yLag0 ) / 2;

    realT x0 = xLag0 - m_maxLag;
    realT y0 = yLag0 - m_maxLag;

    if( m_peakMethod == xcorrPeakMethod::gaussFit )
    {
        m_magIm = m_ccIm.block( x0, y0, 2 * m_maxLag + 1, 2 * m_maxLag + 1 );

        m_fitter.setArray( m_magIm.data(), m_magIm.rows(), m_magIm.cols() );
        m_magIm.maxCoeff( &xLag0, &yLag0 );

        m_fitter.setGuess( mn, pk - mn, m_maxLag, m_maxLag, 3, 3, 0 );
        m_fitter.fit();

        xShift = ( m_fitter.x0() + x0 ) - m_refX0;
        yShift = ( m_fitter.y0() + y0 ) - m_refY0;
        ;
    }
    else if( m_peakMethod == xcorrPeakMethod::interpPeak )
    {
        int x, y;
        m_magIm.resize( m_magSize, m_magSize );

        imageMagnify(
            m_magIm, m_ccIm.block( x0, y0, 2 * m_maxLag + 1, 2 * m_maxLag + 1 ), cubicConvolTransform<Scalar>() );

        m_magIm.maxCoeff( &x, &y );
        xShift = ( x * m_tol + x0 ) - m_refX0;
        yShift = ( y * m_tol + y0 ) - m_refY0;
    }
    else if( m_peakMethod == xcorrPeakMethod::centerOfLight )
    {
        Scalar x, y;

        m_magIm = m_ccIm.block( x0, y0, 2 * m_maxLag + 1, 2 * m_maxLag + 1 );
        m_magIm -= m_magIm.minCoeff(); // Must sum to > 0.

        imageCenterOfLight( x, y, m_magIm );

        xShift = ( x + x0 ) - m_refX0;
        yShift = ( y + y0 ) - m_refY0;
    }
    else if( m_peakMethod == xcorrPeakMethod::none )
    {
        int x, y;
        m_ccIm.maxCoeff( &x, &y );
        xShift = x;
        yShift = y;
    }
    else
    {
        mxThrowException( mx::err::invalidconfig, "imageXCorrFFT<ccImT>::operator()", "unknown peak finding method" );
    }
}

template <class ccImT>
template <class imT>
int imageXCorrFFT<ccImT>::operator()( Scalar &xShift, Scalar &yShift, const imT &im )
{
    if( im.rows() != m_rows )
    {
        mxThrowException( mx::err::sizeerr, "imageXCorrFFT", "image must be same size as reference (rows)" );
    }

    if( im.cols() != m_cols )
    {
        mxThrowException( mx::err::sizeerr, "imageXCorrFFT", "image must be same size as reference (rows)" );
    }

    int maxLag = m_maxLag;
    if( maxLag == 0 )
    {
        maxLag = 0.25 * m_rows - 1;
    }

    float maxLag_r = 0.5 * ( 1.0 * m_rows - 1.0 );
    float maxLag_c = 0.5 * ( 1.0 * m_cols - 1.0 );

    // Once ccIm is resized for maxLag appropriately:
    // if(maxLag_r > m_maxLag && m_maxLag != 0) maxLag_r = m_maxLag;
    // if(maxLag_c > m_maxLag && m_maxLag != 0) maxLag_c = m_maxLag;

    // Mask if needed
    if( m_haveMask )
    {
        m_normIm = im * m_maskIm;
    }
    else
    {
        m_normIm = im;
    }

    if( m_normalize )
    {
        realT m = imageMean( m_normIm );
        realT v = imageVariance( m_normIm, m );
        m_normIm = ( m_normIm - m ) / sqrt( v );
    }

    m_fft_fwd( m_ftWork.data(), m_normIm.data() );

    // m_ftIm0 is the conjugated, fftw-normalized, reference image in the Fourier domain
    // So this is the FT of the cross-correlation:
    m_ftWork *= m_ftIm0;

    m_fft_back( m_ccIm.data(), m_ftWork.data() );

    findPeak( xShift, yShift );

    /*
    if(m_peakMethod == xcorrPeakMethod::gaussFit)
    {
       int xLag0,yLag0;
       Scalar pk = m_ccIm.maxCoeff(&xLag0, &yLag0);
       Scalar mn = m_ccIm.minCoeff();

       if(xLag0 - m_maxLag < 0) m_maxLag = xLag0;
       if(xLag0 + 2*m_maxLag + 1 >= m_ccIm.rows()) m_maxLag = (m_ccIm.rows() - 1 - xLag0)/2;
       if(yLag0 - m_maxLag < 0) m_maxLag = yLag0;
       if(yLag0 + 2*m_maxLag + 1 >= m_ccIm.cols()) m_maxLag = (m_ccIm.cols() - 1 - yLag0)/2;

       realT x0 = xLag0-m_maxLag;
       realT y0 = yLag0-m_maxLag;

       m_magIm = m_ccIm.block( x0, y0, 2*m_maxLag+1, 2*m_maxLag+1);

       m_fitter.setArray(m_magIm.data(), m_magIm.rows(), m_magIm.cols());
       m_magIm.maxCoeff(&xLag0, &yLag0);

       m_fitter.setGuess(mn, pk-mn, m_maxLag, m_maxLag, 3, 3, 0);
       m_fitter.fit();

       xShift = (m_fitter.x0() + x0) - m_refX0;
       yShift = (m_fitter.y0() + y0) - m_refY0;;
    }
    else if(m_peakMethod == xcorrPeakMethod::interpPeak)
    {
       int x, y;
       m_magIm.resize(m_magSize, m_magSize);
       m_ccIm.maxCoeff(&x, &y);

       imageMagnify(m_magIm, m_ccIm.block(x-m_maxLag,y-m_maxLag, 2*m_maxLag+1, 2*m_maxLag+1),
    cubicConvolTransform<Scalar>());

       m_magIm.maxCoeff(&x,&y);
       xShift = x*m_tol - m_maxLag;
       yShift = y*m_tol - m_maxLag;
    }
    else if(m_peakMethod == xcorrPeakMethod::centerOfLight)
    {
       Scalar x,y;
       m_ccIm -= m_ccIm.minCoeff(); //Must sum to > 0.

       imageCenterOfLight(x,y,m_ccIm);

       xShift = x - maxLag_r;
       yShift = y - maxLag_c;
    }
    else if(m_peakMethod == xcorrPeakMethod::none)
    {
       xShift = nan("");
       yShift = nan("");
    }
    else
    {
       mxThrowException(mx::err::invalidconfig, "imageXCorrFFT<ccImT>::operator()", "unknown peak finding method");
    }*/

    return 0;
}

template <class ccImT>
template <class im0T, class imT>
int imageXCorrFFT<ccImT>::operator()( Scalar &xShift, Scalar &yShift, im0T &im0, imT &im )
{
    setReference( im0 );
    return operator()( xShift, yShift, im );
}

} // namespace improc
} // namespace mx

#endif // imageXCorrFFT_hpp
