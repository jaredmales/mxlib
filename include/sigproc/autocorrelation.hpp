/** \file autocorrelation.hpp
 * \brief Tools for working with autocorrelations
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

#ifndef autocorrelation_hpp
#define autocorrelation_hpp

#include "../math/fft/fft.hpp"

namespace mx
{

namespace sigproc
{

/// Calculate the autocorrelation of a time-series
/**
 * \tparam T is the real type of the data and autocorrelation.
 *
 * \ingroup signal_processing
 */
template <typename T>
void autocorrelation(
    T *ac,      ///< [out] is the pre-allocated array of length Nac which will contain the autocorrelation on return
    size_t Nac, ///< [in] is the length of ac
    T *sig,     ///< [in] is the input time-series (signal)
    size_t Nsig ///< [in] is the length of the input time-series
)
{

    #pragma omp parallel for
    for( int i = 0; i < Nac; ++i )
    {
        ac[i] = 0;
        for( int j = 0; j < Nsig - i; ++j )
        {
            ac[i] += sig[j] * sig[j + i];
        }
    }

    T norm = ac[0];
    for( int i = 0; i < Nac; ++i )
    {
        ac[i] /= norm;
    }
}

/// Calculate the autocorrelation of a time-series
/**
 * \tparam T is the real type of the data and autocorrelation.
 *
 * \todo this should probably re-allocate if ac.size() != sig.size(), or at least <
 *
 * \ingroup signal_processing
 */
template <typename T>
void autocorrelation( std::vector<T> &ac, ///< [out] will contain the autocorrelation on return.  If ac.size()==0 then
                                          ///< it is resized to sig.size().
                      std::vector<T> &sig ///< [in] is the input time-series (signal)
)
{
    if( ac.size() == 0 )
    {
        ac.resize( sig.size() );
    }

    autocorrelation( ac.data(), ac.size(), sig.data(), sig.size() );
}

/// Functor for calculating the autocorrelation given a PSD
/** Stores the fftT object and related working memory so that
 * repeated calls do not re-allocate or re-plan the FFT.
 *
 * \tparam T is the real type of the PSD and resultig A.C.
 *
 * \ingroup signal_processing
 */
template <typename T>
struct autocorrelationFromPSD
{
    std::vector<std::complex<T>> fftOut;
    std::vector<std::complex<T>> fftIn;

    math::fft::fftT<std::complex<T>, std::complex<T>, 1, 0> fft;

    /// Calculate the A.C. as the inverse FFT of the PSD
    /** This calculates the circular autocorrelation from the PSD.
     *
     */
    void
    operator()( T *ac, ///<  [out] pre-allocated array, on output contains the first Nac points of the autocorrelation
                size_t Nac, ///< [in] the allocated size of ac.
                T *psd,     ///< [in] the 2-sided FFT storage order PSD
                size_t Npsd ///< [in] the number of points in the PSD
    )
    {
        fft.plan( Npsd, MXFFT_FORWARD );

        fftOut.resize( Npsd );
        fftIn.resize( Npsd );

        for( size_t i = 0; i < Npsd; ++i )
            fftIn[i] = psd[i];

        fft( fftOut.data(), fftIn.data() );

        T norm = fftOut[0].real();
        for( size_t i = 0; i < Npsd && i < Nac; ++i )
            ac[i] = fftOut[i].real() / norm;
    }

    /// Calculate the A.C. as the inverse FFT of the PSD
    /** This calculates the circular autocorrelation from the PSD.
     *
     */
    void operator()( std::vector<T> &ac, ///< [out] On output contains the autocorrelation.  If ac.size()==0, it is
                                         ///< allocated to psd.size().
                     std::vector<T> &psd ///<  [in] the 2-sided FFT storage order PSD
    )
    {
        if( ac.size() == 0 )
            ac.resize( psd.size() );
        operator()( ac.data(), ac.size(), psd.data(), psd.size() );
    }
};

} // namespace sigproc
} // namespace mx

#endif // autocorrelation_hpp
