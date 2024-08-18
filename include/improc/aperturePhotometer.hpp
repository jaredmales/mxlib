/** \file aperturePhotometer.hpp
 * \brief Class for conducting aperture photometry on an image.
 * \ingroup image_processing_files
 * \author Jared R. Males (jaredmales@gmail.com)
 *
 */

//***********************************************************************//
// Copyright 2018-2022 Jared R. Males (jaredmales@gmail.com)
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

#ifndef improc_aperturePhotometer_hpp
#define improc_aperturePhotometer_hpp

#include <map>
#include <mx/improc/eigenImage.hpp>
#include <mx/improc/imageMasks.hpp>

namespace mx
{
namespace improc
{

/// Class for performing aperture photometry on images.
/** Designed to very efficiently calculate the cumulative flux as a function of radius,
 * with minimal cost for performing the measurement on multiple images.
 *
 * This initial version assumes the source is at 0.5*(rows-1) and 0.5*(cols-1).
 *
 * First call resize with the desired size.  This organizes the radius vector and
 * the index image which maps pixels to position in the radius vector.
 *
 * Then call cumPhot with the image, which will sum the flux contained within each radius.
 *
 * \ingroup improc
 *
 * \tparam realT the real floating point type of the data
 */
template <typename realT>
class aperturePhotometer
{

  protected:
    int m_sizeX{ 0 }; ///< The size of the image in X (rows).
    int m_sizeY{ 0 }; ///< The size of the image in Y (columns).

    realT m_xcen{ 0 }; ///< The x coordinate of the center of the image.  Radii are caculated relative to this point.
                       ///< Default is 0.5*(sizeX-1).
    realT m_ycen{ 0 }; ///< The y coordinate of the center of the image.  Radii are caculated relative to this point.
                       ///< Default is 0.5*(sizeY-1).

    std::vector<realT> m_radius; ///< Holds the ordered unique radii of the pixels in the image.

    eigenImage<size_t> m_indexIm; ///< Maps a pixel in the image to the position in the radius vector.

    bool m_useNanMask{ false };

    eigenImage<realT> *m_nanMask{ nullptr };

    realT m_bgMin{ 0 };
    realT m_bgMax{ 0 };

    std::vector<realT> m_bgx;
    std::vector<realT> m_bgy;

  public:
    /// Resize the photometer, recalculating the radius vector and index image.
    /** If size and center are unchanged, nothing is done.
     *
     * \note Source must be at geometric center of image.
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int resize( int sizeX,  ///< [in] The new size in rows
                int sizeY,  ///< [in] The new size in columns
                realT xcen, ///< The x coordinate of the center of the image.  Radii are caculated relative to this
                            ///< point. Default is 0.5*(sizeX-1).
                realT ycen  ///< The y coordinate of the center of the image.  Radii are caculated relative to this
                            ///< point. Default is 0.5*(sizeY-1).
    );

    /// Resize the photometer for a centered image, recalculating the radius vector and index image.
    /** If size and center are unchanged, nothing is done.
     *
     * This version uses the geometric center of the image.
     *
     * \overload
     *
     * \returns 0 on success
     * \returns -1 on error f
     */
    int resize( int sizeX, ///< [in] The new size in rows
                int sizeY  ///< [in] The new size in columns
    );

    /// Get the cumulative photometry of an image as a function of radius.
    /**
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int cumPhot(
        std::vector<realT> &cumPhot, ///< [out] the cumulative photometry at each point in the radius vector.  Resized.
        eigenImage<realT> &im,       ///< [in] the image on which to perform the calculation
        realT maxr =
            0 ///< [in] [optional] the maximum radius to which to calculate.  If <= 0 then max possible is used.
    );

    /// Get the cumulative photometry of an image as a function of radius.
    /**
     * \overload
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int cumPhot(
        std::vector<realT> &cumPhot, ///< [out] the cumulative photometry at each point in the radius vector.  Resized.
        std::vector<realT> &
            deltaPhot, ///< [out] [optional] the photometry at each point in the radius vector, not-cumulative. Resized.
        eigenImage<realT> &im, ///< [in] the image on which to perform the calculation
        realT maxr =
            0 ///< [in] [optional] the maximum radius to which to calculate.  If <= 0 then max possible is used.
    );

  protected:
    /// Get the cumulative photometry of an image as a function of radius.
    /** This is called by the public cumPhot functions.
     *
     *
     * \returns 0 on success
     * \returns -1 on error
     */
    int cumPhotWork(
        std::vector<realT> &cumPhot, ///< [out] the cumulative photometry at each point in the radius vector.  Resized.
        std::vector<realT> *deltaPhot, ///< [out] if not `nullptr`, this vector is filled with the sum at each radius.
        eigenImage<realT> &im,         ///< [in] the image on which to perform the calculation
        realT maxr =
            0 ///< [in] [optional] the maximum radius to which to calculate.  If <= 0 then max possible is used.
    );

  public:
    /// Get the radius value at an index of the vector.
    realT radius( size_t i /**< [in] the index of the radius vector */ )
    {
        return m_radius[i];
    }

    void bgMin( realT bgm )
    {
        m_bgMin = bgm;
    }

    void bgMax( realT bgx )
    {
        m_bgMax = bgx;
    }

    void bg( realT bgm, realT bgx )
    {
        m_bgMin = bgm;
        m_bgMax = bgx;
    }
};

template <typename realT>
int aperturePhotometer<realT>::resize( int sizeX, int sizeY, realT xcen, realT ycen )
{
    // Don't bother if we don't have to.
    if( sizeX == m_sizeX && sizeY == m_sizeY && xcen == m_xcen && ycen == m_ycen )
        return 0;

    m_sizeX = sizeX;
    m_sizeY = sizeY;
    m_xcen = xcen;
    m_ycen = ycen;

    m_indexIm.resize( m_sizeX, m_sizeY );

    // Get radius of each pixel up front, since we need it twice.
    eigenImage<realT> radIm;
    radIm.resize( m_sizeX, m_sizeY );
    radiusImage( radIm, m_xcen, m_ycen ); // mx function to fill in the radius array

    std::map<realT, size_t> unrads; // Map of unique radii to their indices in the radius vector.

    // First populate the map
    for( int i = 0; i < radIm.cols(); ++i )
    {
        for( int j = 0; j < radIm.rows(); ++j )
        {
            // unrads[radIm(j,i)] = 0;
            unrads.insert( std::pair<realT, size_t>( radIm( j, i ), 0 ) );
        }
    }

    // Now fill in the radius vector, and set map values to the indices
    m_radius.resize( unrads.size() );
    typename std::map<realT, size_t>::iterator it;
    int n = 0;
    for( it = unrads.begin(); it != unrads.end(); ++it )
    {
        m_radius[n] = it->first;
        it->second = n;
        ++n;
    }
    // Finally, fill in the index image with index of the radius vector for that pixel.
    for( int i = 0; i < radIm.cols(); ++i )
    {
        for( int j = 0; j < radIm.rows(); ++j )
        {
            it = unrads.find( radIm( j, i ) );
            m_indexIm( j, i ) = it->second;
        }
    }

    if( m_bgMax > m_bgMin && m_bgMin > 0 )
    {
        m_bgx.clear();
        m_bgy.clear();

        for( int cc = 0; cc < radIm.cols(); ++cc )
        {
            for( int rr = 0; rr < radIm.rows(); ++rr )
            {
                if( radIm( rr, cc ) >= m_bgMin && radIm( rr, cc ) <= m_bgMax )
                {
                    m_bgx.push_back( rr );
                    m_bgy.push_back( cc );
                }
            }
        }
    }

    return 0;
}

template <typename realT>
int aperturePhotometer<realT>::resize( int sizeX, int sizeY )
{
    realT xcen = 0.5 * ( 1.0 * sizeX - 1.0 );
    realT ycen = 0.5 * ( 1.0 * sizeY - 1.0 );

    return resize( sizeX, sizeY, xcen, ycen );
}

template <typename realT>
int aperturePhotometer<realT>::cumPhot( std::vector<realT> &cumPhot, eigenImage<realT> &im, realT maxr )
{
    return cumPhotWork( cumPhot, nullptr, im, maxr );
}

template <typename realT>
int aperturePhotometer<realT>::cumPhot( std::vector<realT> &cumPhot,
                                        std::vector<realT> &deltaPhot,
                                        eigenImage<realT> &im,
                                        realT maxr )
{
    return cumPhotWork( cumPhot, &deltaPhot, im, maxr );
}

template <typename realT>
int aperturePhotometer<realT>::cumPhotWork( std::vector<realT> &cumPhot,
                                            std::vector<realT> *deltaPhot,
                                            eigenImage<realT> &im,
                                            realT maxr )
{
    realT xcen = m_xcen; // 0.5*(m_indexIm.rows()-1);
    realT ycen = m_ycen; // 0.5*(m_indexIm.cols()-1);

    if( maxr <= 0 )
        maxr = m_radius.back();

    cumPhot.resize( m_radius.size(), 0 );

    // Make sure we stay in bounds
    int x0 = xcen - maxr;
    if( x0 < 0 )
        x0 = 0;
    int x1 = xcen + maxr;
    if( x1 > m_indexIm.rows() - 1 )
        x1 = m_indexIm.rows() - 1;

    int y0 = ycen - maxr;
    if( y0 < 0 )
        y0 = 0;
    int y1 = ycen + maxr;
    if( y1 > m_indexIm.cols() - 1 )
        y1 = m_indexIm.rows() - 1;

    realT bg = 0;
    // If using background annulus, first calc background;
    if( m_bgx.size() > 0 )
    {
        std::vector<realT> bgann( m_bgx.size() );
        for( size_t n = 0; n < m_bgx.size(); ++n )
        {
            bgann[n] = im( (int)m_bgx[n], (int)m_bgy[n] );
        }

        bg = math::vectorMedian( bgann );
        std::cerr << "background: " << bg << "\n";
    }

    // First get sum at each unique radius
    // Note: switched order of cols/rows indices for speeed
    if( maxr == m_radius.back() )
    {
        for( int i = y0; i <= y1; ++i )
        {
            for( int j = x0; j <= x1; ++j )
            {
                cumPhot[m_indexIm( j, i )] += im( j, i ) - bg;
            }
        }
    }
    else
    {
        for( int i = y0; i <= y1; ++i )
        {
            for( int j = x0; j <= x1; ++j )
            {
                if( m_radius[m_indexIm( j, i )] <= maxr )
                    cumPhot[m_indexIm( j, i )] += im( j, i );
            }
        }
    }

    // If deltaPhot is desired, output it here.
    if( deltaPhot != nullptr )
    {
        deltaPhot->assign( cumPhot.begin(), cumPhot.end() );
    }

    // Now perform cumulative sum
    for( int i = 1; i < cumPhot.size(); ++i )
    {
        cumPhot[i] += cumPhot[i - 1];
    }

    return 0;
}

} // namespace improc
} // namespace mx

#endif // improc_aperturePhotometer_hpp
