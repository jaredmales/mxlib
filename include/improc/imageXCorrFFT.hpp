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
#include "../math/fft/ft.hpp"
#include "../math/fit/fitGaussian.hpp"

#include "eigenImage.hpp"
#include "imageUtils.hpp"
#include "imageTransforms.hpp"
#include "imageXCorr.hpp"

namespace mx
{
namespace improc
{

/// Find the optimum shift to align two images using the FFT cross correlation.
/** The reference image must be the same size as the target image.  Both the reference image and the
 * target image being analyzed can be masked, normalized by mean subtraction and variance division,
 * and windowed.  Seperate masks and windows for the reference and target image are used.
 *
 * The reference mask, reference window, and normalization flag must be set before setting the reference
 * (unless using the defaults).  If the reference window, reference window, or normalization flag are changed
 * then the reference must be set again. After configuring the reference, then operator() can be called
 * repeatedly to determine the shifts for a sequence of imaages.  The mask and window can be changed between
 * calls without requiring resetting the reference. No new heap allocations take place on these calls
 * to operator(), and the reference image is not re-normalized on each call.
 *
 * The shift is reported in pixels such that if the mxlib imageShift function is used
 * to shift the input image by the negative of the shifts, it will align with the
 * reference at the center of the array.
 *
 * Several peak finding methods are provided:
 *  - xcorrPeakMethod::centerOfLight uses center of light
 *  - xcorrPeakMethod::gaussFit uses Gaussian centroiding
 *  - xcorrPeakMethod::interpPeak uses interpolation to find the peak to a given tolerance.
 *  - xcorrPeakMethod::none will just return the coordinate of the max pixel
 *
 *
 * \tparam _realImageT is the Eigen-like array type used for image processing.  See typedefs.
 *
 * \ingroup image_reg
 */
template <class _realImageT>
class imageXCorrFFT
{
  public:
    typedef _realImageT realImageT; ///< the Eigen-like array type used for image processing

    typedef typename _realImageT::Scalar realT; ///< the scalar type of the image type

    typedef std::complex<realT> complexT; ///< Complex floating point type.

    typedef eigenImage<complexT> complexArrayT; ///< Complex eigen array type with Scalar==complexT

  protected:
    /** \name Configuration Member Data
     * @{
     */
    int m_rows{ 0 }; ///< The number of rows in the images

    int m_cols{ 0 }; ///< The number of columns in the images

    realT m_oversamp{ 1 }; ///< The oversampling factor for the CC

    int m_rowsPadded{ 0 }; ///< The number of rows in the padded CC

    int m_colsPadded{ 0 }; ///< The number of columns in the padded CC

    realImageT m_refMaskIm; /**< Mask image to use for the reference, may be needed for proper normalization
                            even if refIm has 0 mask applied.*/

    bool m_haveRefMask{ false }; ///< Flag indicating that a referece mask has been provided.

    realImageT m_maskIm; /**< Mask image to use, may be needed for proper normalization even if
                         refIm has 0 mask applied.*/

    bool m_haveMask{ false }; ///< Flag indicating that a mask has been provided.

    realImageT m_refWinIm; ///< Window image to use for the reference

    bool m_haveRefWindow{ false }; ///< Flag indicating that a window has been provided for the reference

    realImageT m_winIm; ///< Window image for the target image.

    bool m_haveWindow{ false }; ///< Flag indicating that a window has been provided.

    bool m_normalize{ true }; /**< Flag specifying whether images should be normalized prior to performing
                                   the cross-correlation*/

    realImageT m_refIm; ///< The normalized reference image.

    bool m_refValid{ false }; /**< Flag indicating whether the reference is valid.  It can be invalidated
                                   if a mask, window, or normalize is set after the reference is set.*/

    int m_maxLag{ 5 }; ///< The maximum lag to consider in the initial cross-correlation.  Default is 5.

    realT m_tol{ 0.1 }; ///< The tolerance of the interpolated-magnified image, in pixels.

    realT m_magSize{ 0 }; ///< Magnified size of the ccIm when using interp.  Set as function of m_tol and m_maxLag.

    xcorrPeakMethod m_peakMethod{ xcorrPeakMethod::interpPeak };

    ///@}

    /** \name Working Memory Member Data
     * @{
     */

    realImageT m_refACIm; ///< The auto-correlation image of the reference.

    realT m_refX0{ 0 }; /**< The x-shift of the reference image to itself using the selected
                             algorithm, used as coordinate origin*/

    realT m_refY0{ 0 }; /**< The x-shift of the reference image to itself using the selected
                             algorithm, used as coordinate origin*/

    realImageT m_normIm; ///< The normalized image.

    realImageT m_ccIm; ///< The cross-correlation image

    realImageT m_magIm; ///< The magnified image, used if m_peakMethod == xcorrPeakMethod::interpPeak

    complexArrayT m_ftIm0; ///< Working memory for the FT of the reference image.

public:
    complexArrayT m_ftWork; ///< Working memory for the FFT.

    #ifdef XCFFT_C2C
    complexArrayT m_ftWorkIn; ///< Working memory for the FFT input.
    #endif

    complexArrayT m_ftWorkPadded; ///< Working memory for the padded inverse FFT input.

    #ifdef XCFFT_C2C
    complexArrayT m_ftWorkPaddedOut; ///< Working memory for the FFT output.
    #endif
protected:

    #ifdef XCFFT_C2C

    math::ft::fftT<complexT, complexT, 2, 0> m_fft_fwd; ///< FFT object for the forward transform.

    math::ft::fftT<complexT, complexT, 2, 0> m_fft_bwd; ///< FFT object for the backward transfsorm.

    #else

    math::ft::fftT<realT, complexT, 2, 0> m_fft_fwd; ///< FFT object for the forward transform.

    math::ft::fftT<complexT, realT, 2, 0> m_fft_bwd; ///< FFT object for the backward transfsorm.

    #endif

    ///@}

    /** \name Peak Finding Tools
     * @{
     */
    math::fit::fitGaussian2D<mx::math::fit::gaussian2D_gen_fitter<realT>> m_fitter;

    ///@}

  public:
    /** \name Construction and Destruction
     * @{
     */

    /// Default c'tor
    imageXCorrFFT();

    ///@}

    /** \name Configuration Member Access
     * @{
     */

    /// Get the number of rows in the input images
    /** This can only be set by a call to resize() or refIm().
     *
     * \returns the current value of m_rows
     */
    int rows();

    /// Get the number of columns in the input images
    /** This can only be set by a call to resize() or refIm().
     *
     * \returns the current value of m_cols
     */
    int cols();

    /// Set the oversampling factor
    /** This can also be set by a call to resize() or refIm().
     *  The oversampling factor must be greater than or equal to 1.
     *
     *  If m_rows and m_cols are both greater than 0 this calls resize(int, int, realT).
     *
     */
    int oversamp( realT os /**< [in] the new oversampling factor*/ );

    /// Get the oversampling factor
    /**
     * \returns the current value of m_oversamp
     */
    realT oversamp();

    /// Get the number of rows in the padded cross-correlation
    /** This can only be set by a call to resize() or refIm().
     *
     * \returns the current value of m_rowsPadded
     */
    int rowsPadded();

    /// Get the number of rows in the padded cross-correlation
    /** This can only be set by a call to resize() or refIm().
     *
     * \returns the current value of m_colsPadded
     */
    int colsPadded();

    /// Set the size of the cross-correlation problem based on input image size.
    /** This resizes all working memory and conducts fftw planning.
     *  If there are no changes to size or oversampling, this returns immediately
     *  making no changes.
     *
     *  If there are changes to size or oversampling the reference is invalidated.
     *
     *  If nrows or ncols is 0, this resizes all working memory to 0.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int resize( int nrows, /**< [in] the number of rows in the images to register */
                int ncols, /**< [in] the number of columns in the images to register */
                realT oversamp /**< [in] the oversampling factor.  Must be `>=` 1. */ );

    /// Set the size of the cross-correlation problem based on input image size.
    /** Calls resize(int, int, realT) with the current value of m_oversamp.
     *
     * \overload
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int resize( int nrows, /**< [in] the number of rows in the images to register */
                int ncols /**< [in] the number of columns in the images to register */ );

    /// Set the reference mask image
    /** This invalidates the reference
     * \returns 0 on success
     * \returns -1 on error
     */
    int refMaskIm( const realImageT &mask /**< [in] the new reference mask image */ );

    /// Get a reference to the reference mask image.
    /**
     * \returns a const reference to the reference mask image.
     */
    const realImageT &refMaskIm();

    /// Get the value of the haveRefMask flag
    /**
     * \returns the current value of m_haveRefMask
     */
    bool haveRefMask();

    /// Set the mask image
    /** This does not invalidate the reference
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int maskIm( const realImageT &mask /**< [in] the new mask image */ );

    /// Get a reference to the mask image.
    /**
     * \returns a const referance to the mask image.
     */
    const realImageT &maskIm();

    /// Get the value of the haveMask flag
    /**
     * \returns the current value of m_haveMask
     */
    bool haveMask();

    /// Set the reference window image
    /** This invalidates the reference.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int refWinIm( const realImageT &win /**< [in] the new reference window image */ );

    /// Get a reference to the reference window image.
    /**
     * \returns a const referance to the reference window image.
     */
    const realImageT &refWinIm();

    /// Get the value of the haveRefWindow flag
    /**
     * \returns the current value of m_haveRefWindow
     */
    bool haveRefWindow();

    /// Set the window image
    /** This does not invalidate the reference.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int winIm( const realImageT &win /**< [in] the new window image */ );

    /// Get a reference to the window image.
    /**
     * \returns a const referance to the window image.
     */
    const realImageT &winIm();

    /// Get the value of the haveWindow flag
    /**
     * \returns the current value of m_haveWindow
     */
    bool haveWindow();

    /// Set the normalize flag
    /** Changing this invalidates the reference.
     */
    void normalize( bool no /**< [in] the new normalize flag value */ );

    /// Get the current value of the normalize flag
    /**
     * \returns the current value of m_normalize
     */
    bool normalize();

    /// Set the reference image
    /** Performs the following steps
     *    - Applies the mask if supplied.
     *    - If desired, normalizes the reference image by mean subtraction and variance division.
     *    - Applies the window if supplied.
     *    - Stores the result as m_refIm
     *    - Resizes working memory with a call to resize()
     *    - Shifts to FFT order and Fourier transforms and conjugates the result.
     *    - Calls operator() to measure the reference shifts
     *    - Stores the resultant auto-correlation in m_refACIm
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int refIm( const realImageT &im0, /**< [in] The new reference image */
               realT oversamp /**< [in] [optional] the oversampling factor to use */ );

    /// Set the reference image
    /** Calls refIm(const realImageT &, realT) with the current value of m_oversamp.
     *
     *  \overload
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int refIm( const realImageT &im0 /**< [in] The new reference image */ );

    /// Get a reference to the reference image.
    /**
     * \returns a const referent to m_refIm.
     */
    const realImageT &refIm();

    /// Get the value of the reference valid flag
    /**
     * \returns the value of m_refValid
     */
    bool refValid();

    /// Set the maximum lag
    void maxLag( int ml /**< [in] the new maximum lag */ );

    /// Get the current maximum lag
    /**
     * \returns the current value of m_maxLag
     */
    int maxLag();

    /// Set the tolerance of the interpolated-magnified image, in pixels.
    void tol( realT nt /**< [in] The new value of the interpolation tolerance. */ );

    /// Get the tolerance of the interpolated-magnified image, in pixels.
    /**
     * \returns the current value of m_tol.
     */
    realT tol();

    void peakMethod( xcorrPeakMethod xpm );

    xcorrPeakMethod peakMethod();

    ///@}

    /** \name Working Memory Member Access
     * @{
     */

    /// Get a reference to the reference auto-correlation image.
    /**
     * \returns a const referent to m_refACIm.
     */
    const realImageT &refACIm();

    /// Get the x-shift of the reference image to itself.
    /**
     * \returns m_refX0
     */
    realT refX0();

    /// Get the y-shift of the reference image to itself.
    /**
     * \returns m_refY0
     */
    realT refY0();

    /// Get a reference to the normalized image.
    /**
     * \returns a const referent to m_normIm.
     */
    const realImageT &normIm();

    /// Get a reference to the cross correlation image.
    /**
     * \returns a const referent to m_ccIm.
     */
    const realImageT &ccIm();

    /// Get a reference to the magnified image.
    /**
     * \returns a const referent to m_magIm.
     */
    const realImageT &magIm();

    ///@}

    /** \name Peak Finding
     * @{
     */
  protected:
    void findPeak( realT &xShift, /**< [out] the x shift of im w.r.t. im0, in pixels */
                   realT &yShift /**< [out] the y shift of im w.r.t. im0, in pixels */ );

    ///@}

    /** \name Cross Correlation Operator
     * @{
     */

  public:
    /// Conduct the cross correlation to a specified tolerance
    /**
     * \returns 0 on success
     * \returns -1 on error
     */
    template <class imT>
    int operator()( realT &xShift, ///< [out] the x shift of im w.r.t. im0, in pixels
                    realT &yShift, ///< [out] the y shift of im w.r.t. im0, in pixels
                    const imT &im  ///< [in] the image to cross-correlate with the reference
    );

    /// Conduct the cross correlation to a specified tolerance
    /**
     * \overload
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    template <class im0T, class imT>
    int operator()( realT &xShift, ///< [out] the x shift of im w.r.t. im0, in pixels
                    realT &yShift, ///< [out] the y shift of im w.r.t. im0, in pixels
                    im0T &im0,     ///< [in] a new reference image
                    imT &im        ///< [in] the image to cross-correlate with the reference
    );

    ///@}
};

template <class realImageT>
imageXCorrFFT<realImageT>::imageXCorrFFT()
{
    maxLag( m_maxLag ); // this sets up m_magSize.
}

template <class realImageT>
int imageXCorrFFT<realImageT>::rows()
{
    return m_rows;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::cols()
{
    return m_cols;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::oversamp( realT os )
{
    if( os < 1 )
    {
        mxError( "imageXCorrFFT::oversamp", MXE_INVALIDARG, "oversampling factor can't be less than 1" );
        return MXE_INVALIDARG;
    }

    if( m_rows > 0 && m_cols > 0 )
    {
        return resize( m_rows, m_cols, os );
    }
    else
    {
        m_oversamp = os;
        return 0;
    }
}

template <class realImageT>
typename imageXCorrFFT<realImageT>::realT imageXCorrFFT<realImageT>::oversamp()
{
    return m_oversamp;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::rowsPadded()
{
    return m_rowsPadded;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::colsPadded()
{
    return m_colsPadded;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::resize( int nrows, int ncols, realT oversamp )
{
    // No-op if no change
    if( m_rows == nrows && m_cols == ncols && m_oversamp == oversamp )
    {
        return 0;
    }

    if( oversamp < 1 )
    {
        mxError( "imageXCorrFFT::resize", MXE_INVALIDARG, "oversampling factor can't be less than 1" );
        return MXE_INVALIDARG;
    }

    m_rows = nrows;

    m_cols = ncols;

    m_oversamp = oversamp;

    /* Any change here invalidates the reference.
     * This is trivially true for rows or cols.  Also must be
     * handled for oversampling which will change the reference shifts.
     */
    m_refValid = false;

    if( m_rows == 0 || m_cols == 0 )
    {
        m_ftIm0.resize( 0, 0 );
        m_ftWork.resize( 0, 0 );
        #ifdef XCFFT_C2C
        m_ftWorkIn.resize( 0, 0 );
        #endif
        m_rowsPadded = 0;
        m_colsPadded = 0;
        m_ccIm.resize( 0, 0 );
        m_ftWorkPadded.resize( 0, 0 );

        #ifdef XCFFT_C2C
        m_ftWorkPaddedOut.resize( 0, 0 );
        #endif

        return 0;
    }

    #ifdef XCFFT_C2C

    m_ftIm0.resize( m_rows, m_cols );

    m_ftWork.resize( m_rows, m_cols );

    m_ftWorkIn.resize( m_rows, m_cols );

    #else

    m_ftIm0.resize( (int)( 0.5 * m_rows ) + 1, m_cols );

    m_ftWork.resize( (int)( 0.5 * m_rows ) + 1, m_cols );

    #endif

    // fftw is row-major, eigen defaults to column-major
    m_fft_fwd.plan( m_cols, m_rows, math::ft::dir::forward, false );

    m_rowsPadded = ceil( m_rows * m_oversamp );
    m_colsPadded = ceil( m_cols * m_oversamp );

    m_ccIm.resize( m_rowsPadded, m_colsPadded );

    #ifdef XCFFT_C2C

    m_ftWorkPadded.resize( m_rowsPadded, m_colsPadded );
    m_ftWorkPaddedOut.resize( m_rowsPadded, m_colsPadded );

    #else

    m_ftWorkPadded.resize( (int)( 0.5 * m_rowsPadded ) + 1, m_colsPadded );

    #endif

    m_fft_bwd.plan( m_colsPadded, m_rowsPadded, math::ft::dir::backward, false );

    return 0;
} // resize(int, int, realT)

template <class realImageT>
int imageXCorrFFT<realImageT>::resize( int nrows, int ncols )
{
    return resize( nrows, ncols, m_oversamp );
}

template <class realImageT>
int imageXCorrFFT<realImageT>::refMaskIm( const realImageT &mask )
{
    m_refMaskIm = mask;
    m_haveRefMask = true;

    m_refValid = false;

    return 0;
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::refMaskIm()
{
    return m_refMaskIm;
}

template <class realImageT>
bool imageXCorrFFT<realImageT>::haveRefMask()
{
    return m_haveRefMask;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::maskIm( const realImageT &mask )
{
    m_maskIm = mask;
    m_haveMask = true;

    return 0;
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::maskIm()
{
    return m_maskIm;
}

template <class realImageT>
bool imageXCorrFFT<realImageT>::haveMask()
{
    return m_haveMask;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::refWinIm( const realImageT &win )
{
    m_refWinIm = win;
    m_haveRefWindow = true;

    m_refValid = false;

    return 0;
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::refWinIm()
{
    return m_refWinIm;
}

template <class realImageT>
bool imageXCorrFFT<realImageT>::haveRefWindow()
{
    return m_haveRefWindow;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::winIm( const realImageT &win )
{
    m_winIm = win;
    m_haveWindow = true;

    return 0;
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::winIm()
{
    return m_winIm;
}

template <class realImageT>
bool imageXCorrFFT<realImageT>::haveWindow()
{
    return m_haveWindow;
}

template <class realImageT>
void imageXCorrFFT<realImageT>::normalize( bool no )
{
    if( no != m_normalize )
    {
        m_refValid = false;
        m_normalize = no;
    }
}

template <class realImageT>
bool imageXCorrFFT<realImageT>::normalize()
{
    return m_normalize;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::refIm( const realImageT &im, realT oversamp )
{
    realImageT im0;

    // Mask if needed
    if( m_haveRefMask )
    {
        if( im.rows() != m_refMaskIm.rows() && im.cols() != m_refMaskIm.cols() )
        {
            mxError( "imageXCorrFFT::refIm", MXE_SIZEERR, "reference and reference mask are not the same size" );
            return MXE_SIZEERR;
        }

        // Now normalize
        if( m_normalize )
        {
            realT m = imageMean<realT>( im, m_refMaskIm );
            realT v = imageVariance<realT>( im, m, m_refMaskIm );
            im0 = ( ( im - m ) * m_refMaskIm ) / sqrt( v );
        }
        else
        {
            im0 = im * m_refMaskIm;
        }
    }
    else
    {
        // Normalize
        if( m_normalize )
        {
            realT m = imageMean<realT>( im );
            realT v = imageVariance<realT>( im, m );
            im0 = ( im - m ) / sqrt( v );
        }
        else
        {
            im0 = im;
        }
    }

    // Window if needed
    if( m_haveRefWindow )
    {
        if( im0.rows() != m_refWinIm.rows() && im.cols() != m_refWinIm.cols() )
        {
            mxError( "imageXCorrFFT::refIm", MXE_SIZEERR, "reference and reference window are not the same size" );
            return MXE_SIZEERR;
        }

        im0 *= m_refWinIm;
    }

    // We save refIm as the un-shifted version
    m_refIm = im0;

    /* Setup the FFTW space
     * Note: we don't do this earlier in case we ever implement something (e.g. a filter)
     * that changes size before this point
     */
    resize( im0.rows(), im0.cols(), oversamp );

    // Now shift so center pixel is 0,0
    imageShiftWP( im0, m_refIm, 0.5 * m_rows + 1, 0.5 * m_cols + 1 );

    // Then FT
    #ifdef XCFFT_C2C
    for(int cc=0; cc < im0.cols(); ++cc)
    {
        for(int rr=0; rr < im0.rows(); ++rr)
        {
            m_ftWork(rr,cc) = im0(rr,cc);
        }
    }
    m_fft_fwd( m_ftIm0.data(), m_ftWork.data() );
    #else
    m_fft_fwd( m_ftIm0.data(), im0.data() );
    #endif

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

    m_refValid = true;

    operator()( m_refX0, m_refY0, m_refIm );

    m_refACIm = m_ccIm;

    return 0;
}

template <class realImageT>
int imageXCorrFFT<realImageT>::refIm( const realImageT &im )
{
    return refIm( im, m_oversamp );
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::refIm()
{
    return m_refIm;
}

template <class realImageT>
bool imageXCorrFFT<realImageT>::refValid()
{
    return m_refValid;
}

template <class realImageT>
void imageXCorrFFT<realImageT>::maxLag( int ml )
{
    m_maxLag = ml;
    tol( m_tol );
}

template <class realImageT>
int imageXCorrFFT<realImageT>::maxLag()
{
    return m_maxLag;
}

template <class realImageT>
typename imageXCorrFFT<realImageT>::realT imageXCorrFFT<realImageT>::tol()
{
    return m_tol;
}

template <class realImageT>
void imageXCorrFFT<realImageT>::tol( realT nt )
{
    m_magSize = ceil( ( ( 2. * m_maxLag + 1 ) - 1.0 ) / nt ) + 1;

    realT mag = ( m_magSize - 1.0 ) / ( ( 2. * m_maxLag + 1 ) - 1.0 );

    m_tol = 1.0 / mag;

    if( m_refIm.rows() != m_rows || m_refIm.cols() != m_cols || m_rows == 0 || m_cols == 0 )
        return;

    // Find the reference shift
    m_refX0 = 0;
    m_refY0 = 0;
    operator()( m_refX0, m_refY0, m_refIm );
}

template <class realImageT>
void imageXCorrFFT<realImageT>::peakMethod( xcorrPeakMethod xpm )
{
    m_peakMethod = xpm;

    if( m_refIm.rows() != m_rows || m_refIm.cols() != m_cols || m_rows == 0 || m_cols == 0 )
        return;

    // Find the reference shift
    m_refX0 = 0;
    m_refY0 = 0;
    operator()( m_refX0, m_refY0, m_refIm );
}

template <class realImageT>
xcorrPeakMethod imageXCorrFFT<realImageT>::peakMethod()
{
    return m_peakMethod;
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::refACIm()
{
    return m_refACIm;
}

template <class realImageT>
typename imageXCorrFFT<realImageT>::realT imageXCorrFFT<realImageT>::refX0()
{
    return m_refX0;
}

template <class realImageT>
typename imageXCorrFFT<realImageT>::realT imageXCorrFFT<realImageT>::refY0()
{
    return m_refY0;
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::normIm()
{
    return m_normIm;
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::ccIm()
{
    return m_ccIm;
}

template <class realImageT>
const realImageT &imageXCorrFFT<realImageT>::magIm()
{
    return m_magIm;
}

template <class realImageT>
void imageXCorrFFT<realImageT>::findPeak( realT &xShift, realT &yShift )
{
    int xLag0, yLag0;

    /*realImageT sm = m_ccIm;
    filterImage(sm, m_ccIm, gaussKernel<eigenImage<realT>, 2>(3));
    m_ccIm = sm;*/

    // m_ccIm = m_ccIm.pow(3);
    realT pk = m_ccIm.maxCoeff( &xLag0, &yLag0 );
    realT mn = m_ccIm.minCoeff();

    if( xLag0 - m_maxLag < 0 )
    {
        m_maxLag = xLag0;
    }
    if( xLag0 + 2 * m_maxLag + 1 >= m_ccIm.rows() )
    {
        m_maxLag = ( m_ccIm.rows() - 1 - xLag0 ) / 2;
    }
    if( yLag0 - m_maxLag < 0 )
    {
        m_maxLag = yLag0;
    }
    if( yLag0 + 2 * m_maxLag + 1 >= m_ccIm.cols() )
    {
        m_maxLag = ( m_ccIm.cols() - 1 - yLag0 ) / 2;
    }

    realT x0 = xLag0 - m_maxLag;
    realT y0 = yLag0 - m_maxLag;

    if( m_peakMethod == xcorrPeakMethod::gaussFit )
    {
        m_magIm = m_ccIm.block( x0, y0, 2 * m_maxLag + 1, 2 * m_maxLag + 1 );

        m_fitter.setArray( m_magIm.data(), m_magIm.rows(), m_magIm.cols() );
        m_magIm.maxCoeff( &xLag0, &yLag0 );

        m_fitter.setGuess( mn, pk - mn, m_maxLag, m_maxLag, 3, 3, 0 );
        m_fitter.fit();

        x0 = ( m_fitter.x0() + x0 );
        y0 = ( m_fitter.y0() + y0 );
    }
    else if( m_peakMethod == xcorrPeakMethod::interpPeak )
    {
        int x, y;
        m_magIm.resize( m_magSize, m_magSize );

        imageMagnify(
            m_magIm, m_ccIm.block( x0, y0, 2 * m_maxLag + 1, 2 * m_maxLag + 1 ), cubicConvolTransform<realT>() );

        m_magIm.maxCoeff( &x, &y );
        x0 = ( x * m_tol + x0 );
        y0 = ( y * m_tol + y0 );
    }
    else if( m_peakMethod == xcorrPeakMethod::centerOfLight )
    {
        realT x, y;

        m_magIm = m_ccIm.block( x0, y0, 2 * m_maxLag + 1, 2 * m_maxLag + 1 );
        m_magIm -= m_magIm.minCoeff(); // Must sum to > 0.

        imageCenterOfLight( x, y, m_magIm );

        x0 = ( x + x0 );
        y0 = ( y + y0 );
    }
    else if( m_peakMethod == xcorrPeakMethod::none )
    {
        int x, y;
        m_ccIm.maxCoeff( &x, &y );
        x0 = x;
        y0 = y;
    }
    else
    {
        mxThrowException(
            mx::err::invalidconfig, "imageXCorrFFT<realImageT>::operator()", "unknown peak finding method" );
    }

    //--> unpad here, scaling the shifts
    xShift = x0 * static_cast<realT>( m_rows ) / static_cast<realT>( m_rowsPadded ) - m_refX0;
    yShift = y0 * static_cast<realT>( m_cols ) / static_cast<realT>( m_colsPadded ) - m_refY0;
}

template <class realImageT>
template <class imT>
int imageXCorrFFT<realImageT>::operator()( realT &xShift, realT &yShift, const imT &im )
{
    if( !m_refValid )
    {
        mxThrowException( mx::err::sizeerr, "imageXCorrFFT", "reference image is not valid" );
    }

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

    // Mask and normalize as needed
    if( m_haveMask )
    {
        if( m_normalize )
        {
            realT m = imageMean<realT>( im, m_maskIm );
            realT v = imageVariance<realT>( im, m, m_maskIm );

            m_normIm = ( im - m ) * m_maskIm / sqrt( v );
        }
        else
        {
            m_normIm = im * m_maskIm;
        }
    }
    else
    {

        if( m_normalize )
        {
            realT m = imageMean<realT>( m_normIm );
            realT v = imageVariance<realT>( m_normIm, m );
            m_normIm = ( im - m ) / sqrt( v );
        }
        else
        {
            m_normIm = im;
        }
    }

    if( m_haveWindow )
    {
        m_normIm *= m_winIm;
    }

    #ifdef XCFFT_C2C
    for(int cc=0; cc < m_normIm.cols(); ++cc)
    {
        for(int rr=0; rr < m_normIm.rows(); ++rr)
        {
            m_ftWorkIn(rr,cc) = m_normIm(rr,cc);
        }
    }
    m_fft_fwd( m_ftWork.data(), m_ftWorkIn.data() );
    #else
    m_fft_fwd( m_ftWork.data(), m_normIm.data() );
    #endif

    // m_ftIm0 is the conjugated, fftw-normalized, reference image in the Fourier domain
    // So this is the FT of the cross-correlation:
    m_ftWork *= m_ftIm0;

    //--> PAD here
    int ci = 0.5 * m_ftWork.cols() + 1;
    int cf = m_ftWork.cols() - ci;

    m_ftWorkPadded.setZero(); //fftw c2r overwrites the input
    m_ftWorkPadded.block( 0, 0, m_ftWork.rows(), ci ) = m_ftWork.block( 0, 0, m_ftWork.rows(), ci );
    m_ftWorkPadded.block( 0, m_ftWorkPadded.cols() - cf, m_ftWork.rows(), cf ) =
        m_ftWork.block( 0, m_ftWork.cols() - cf, m_ftWork.rows(), cf );

    #ifdef XCFFT_C2C
    m_fft_bwd( m_ftWorkPaddedOut.data(), m_ftWorkPadded.data() );

    for(int cc = 0; cc < m_ccIm.cols(); ++cc)
    {
        for(int rr = 0; rr < m_ccIm.rows(); ++rr)
        {
            m_ccIm(rr,cc) = std::real(m_ftWorkPaddedOut(rr,cc));
        }
    }

    #else
    m_fft_bwd( m_ccIm.data(), m_ftWorkPadded.data() );
    #endif

    findPeak( xShift, yShift );

    return 0;
}

template <class realImageT>
template <class im0T, class imT>
int imageXCorrFFT<realImageT>::operator()( realT &xShift, realT &yShift, im0T &im0, imT &im )
{
    setReference( im0 );
    return operator()( xShift, yShift, im );
}

} // namespace improc
} // namespace mx

#endif // imageXCorrFFT_hpp
