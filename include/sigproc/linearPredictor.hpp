/** \file linearPredictor.hpp
 * \brief Working with linear prediction.
 *
 * \author Jared R. Males (jaredmales@gmail.com)
 *
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

#ifndef linearPredictor_hpp
#define linearPredictor_hpp

#include <vector>
#include <complex>

#include "../math/constants.hpp"
#include "../math/eigenLapack.hpp"

#include "levinsonRecursion.hpp"

namespace mx
{
namespace sigproc
{

/// A class to support linear prediction.
/** \ingroup signal_processing
 *
 * \todo document linearPredictor
 */
template <typename _realT>
struct linearPredictor
{
    typedef _realT realT;

    std::vector<realT> m_c;

    realT _setCondition{ 0 };
    realT _actCondition{ 0 };
    int _nRejected{ 0 };

    /// Calculate the LP coefficients given an autocorrelation.
    /** If condition==0 then the levinson recursion is used.
     * Otherwise, SVD pseudo-inversion is used with the given condition number.
     */
    int calcCoefficients( const std::vector<realT> &ac,
                         size_t Nc,
                         size_t Npred = 1,
                         realT condition = 0
                         )
    {
        return calcCoefficients( ac.data(), ac.size(), Nc, Npred, condition);
    }

    int calcCoefficients( const realT *ac,
                          size_t acSz,
                          size_t Nc,
                          size_t Npred = 1,
                          realT condition = 0
                          )
    {

        if( condition == 0 )
        {
            return calcCoefficientsLevinson( ac, acSz, Nc, Npred );
        }

        Eigen::Array<realT, -1, -1> Rmat, Rvec, PInv, LPcoeff;

        Rmat.resize( Nc, Nc );
        Rvec.resize( 1, Nc );

        for( int i = 0; i < Nc && i < acSz; ++i )
        {
            for( int j = 0; j < Nc && i < acSz; ++j )
            {
                Rmat( i, j ) = ac[abs( i - j )];
            }

            Rvec( 0, i ) = ac[i + Npred];
        }

        realT tmpCond = condition;

        _setCondition = condition;
        math::eigenPseudoInverse( PInv, tmpCond, _nRejected, Rmat, condition );

        _actCondition = tmpCond;

        m_c.resize(Nc);
        Eigen::Map<Eigen::Array<realT,-1,-1>> cmap(m_c.data(), 1, m_c.size());
        cmap = Rvec.matrix() * PInv.matrix();

        return 0;
    }

    int calcCoefficientsLevinson( const std::vector<realT> &ac, /**< [in] The autocorrelation, at least
                                                                          Nc+Npred in length */
                                  size_t Nc,                    /**< [in] The number of LP coefficients desired */
                                  size_t Npred = 1              /**< [in] [optional] The prediction length,
                                                                                    default is 1 */
    )
    {
        return calcCoefficientsLevinson( ac.data(), ac.size(), Nc, Npred );
    }

    int calcCoefficientsLevinson( const realT *ac,   /**< [in] The autocorrelation, at least Nc+Npred in length */
                                  size_t acSz,       /**< [in] The length of the autocorrelation */
                                  size_t Nc,         /**< [in] The number of LP coefficients desired */
                                  unsigned Npred = 1 /**< [in] [optional] The prediction length, default is 1 */
    )
    {
        if( acSz < Nc + Npred )
        {
            std::string msg = "too many coefficients for size and prediction length\n";
            msg += "    acSz  = " + std::to_string( acSz ) + "\n";
            msg += "    Nc    = " + std::to_string( Nc ) + "\n";
            msg += "    Npred = " + std::to_string( Npred ) + "\n";
            mxThrowException( err::invalidarg, "linearPredictor::calcCoefficientsLevinson", msg );
        }

        std::vector<realT> r, x, y;

        r.resize( 2. * Nc - 1 );
        m_c.resize( Nc );
        y.resize( Nc );

        for( size_t i = 0; i < Nc; ++i )
        {
            r[i] = ac[Nc - i - 1]; // this runs from Nc-1 to 0
        }

        for( size_t i = Nc; i < 2 * Nc - 1; ++i )
        {
            r[i] = ac[i - Nc + 1]; // this runs from 1 to Nc-1
        }

        for( size_t i = 0; i < Nc; ++i )
        {
            y[i] = ac[i + Npred]; // this runs from Npred to Nc-1 + Npred
        }

        levinsonRecursion( r.data(), m_c.data(), y.data(), Nc );

        return 0;
    }

    realT c( size_t i )
    {
        return m_c[i];
    }

    size_t Nc()
    {
        return m_c.size();
    }

    realT predict( std::vector<realT> &hist, int idx )
    {
        realT x = 0;

        if( idx < m_c.size() )
        {
            return x;
        }

        for( int i = 0; i < m_c.size(); ++i )
        {
            x += m_c[i] * hist[idx - i];
        }

        return x;
    }

    realT spectralResponse( realT f, realT fs )
    {
        int n = m_c.size();

        std::complex<realT> He = 0;
        for( int j = 0; j < n; ++j )
        {
            realT s = ( j + 1.0 ) * math::two_pi<realT>();
            He += m_c[j] * exp( s * std::complex<realT>( 0, -1.0 ) * f / fs );
        }

        realT one = 1.0;
        return std::norm( one / ( one - He ) );
    }
};

} // namespace sigproc
} // namespace mx
#endif // linearPredictor_hpp
