/** \file speckleAmpPSD.hpp
 * \author Jared R. Males (jaredmales@gmail.com)
 * \brief Calculates the PSD of speckle intensity given a modified Fourier mode amplitude PSD.
 * \ingroup mxAO_files
 *
 */

#ifndef speckleAmpPSD_hpp
#define speckleAmpPSD_hpp

#include <vector>

#include <Eigen/Dense>

#include "../../math/constants.hpp"
#include "../../math/randomT.hpp"
#include "../../math/vectorUtils.hpp"

#include "../../sigproc/signalWindows.hpp"
#include "../../sigproc/psdUtils.hpp"
#include "../../sigproc/psdFilter.hpp"
#include "../../sigproc/psdVarMean.hpp"
#include "../../sigproc/averagePeriodogram.hpp"

#include "jincFuncs.hpp"

namespace mx
{
namespace AO
{
namespace analysis
{

/// Calculate the PSD of the speckle intensity given the PSD of Fourier mode amplitude
/**
 * Does so by generating a time-series from the PSD using Fourier-domain convolution, and calculating
 * the average periodogram.
 *
 * A Hann window is used.
 *
 * Statistics are not currently normalized.  You need to normalize to match expected contrast if desired.
 *
 * Will also calculate the binned-variances in the generated time-series, if the arguments vars and bins are not null
 * pointers.
 *
 * \todo Figure out what's going on with psdFilter normalization!
 * \todo probably need to not do overlapped averaging in periodogram.  Could drop mean-sub then.
 */
template <typename realT>
int speckleAmpPSD(
    std::vector<realT> &spFreq, ///< [out] The frequency grid of the output PSD
    std::vector<realT>
        &spPSD, ///< [out] The speckle amplitude PSD corresponding to the freq coordinates.  Will be resized.
    const std::vector<realT> &freq,  ///< [in] The Frequency grid of the input PSD
    const std::vector<realT> &fmPSD, ///< [in] The Fourier mode PSD.
    const std::vector<std::complex<realT>>
        &fmXferFxn,                 ///< [in] The complex error transfer function, as a function of freq
    const std::vector<realT> &nPSD, ///< [in] The open-loop noise PSD
    const std::vector<std::complex<realT>> &nXferFxn, ///< [in] The noise transfer function, as a function of freq
    int N,                              ///< [in] The number of trials to use in calculating the amplitude PSD.
    std::vector<realT> *vars = nullptr, ///< [out] [optional] The binned variances of the time series generated.
    std::vector<realT> *bins = nullptr, ///< [in]  [optional] The bin sizes to use in calculating vars.
    bool noPSD = false ///< [in] [optional] if true then the PSD is not actually calculated, only the binned variances
)
{

    std::vector<realT> psd2, npsd2;
    std::vector<std::complex<realT>> xfer2, nxfer2;

    bool hasZero = false;
    if( freq[0] == 0 )
        hasZero = true;

    // First augment to two-sided DFT form
    mx::sigproc::augment1SidedPSD( psd2, fmPSD, !hasZero, 0.5 );
    mx::sigproc::augment1SidedPSD( npsd2, nPSD, !hasZero, 0.5 );

    // And augment the xfer fxns to two sided form
    sigproc::augment1SidedPSD( xfer2, fmXferFxn, !hasZero, 1.0 );
    sigproc::augment1SidedPSD( nxfer2, nXferFxn, !hasZero, 1.0 );

    // The time sampling
    realT dt = 1. / ( 2 * freq.back() );

    // Indices for getting the middle half
    int Nwd = 1.0 * psd2.size();
    int NwdStart = 0.5 * psd2.size() - 0.5 * Nwd;

    int Nsamp = 1.0 * psd2.size();
    int NsampStart = 0.5 * Nwd - 0.5 * Nsamp;

    sigproc::averagePeriodogram<realT> globalAvgPgram( Nsamp * 0.1, dt );
    spPSD.resize( globalAvgPgram.size() );
    for( size_t n = 0; n < spPSD.size(); ++n )
        spPSD[n] = 0;

    std::vector<std::vector<realT>> means;
    if( bins != nullptr && vars != nullptr )
        means.resize( bins->size() );

    // Calculate the Fourier Mode variance for normalization
    realT fmVar = sigproc::psdVar( freq, fmPSD );
    // and the noise variance
    realT nVar = sigproc::psdVar( freq, nPSD );

#pragma omp parallel
    {
        // Filters for imposing the PSDs
        mx::sigproc::psdFilter<realT, 1> filt;
        filt.psd( psd2, freq[1] - freq[0] );

        mx::sigproc::psdFilter<realT, 1> nfilt;
        nfilt.psd( npsd2, freq[1] - freq[0] );

        // FFTs for going to Fourier domain and back to time domain.
        math::fft::fftT<realT, std::complex<realT>, 1, 0> fft( psd2.size() );
        math::fft::fftT<std::complex<realT>, realT, 1, 0> fftB( psd2.size(), MXFFT_BACKWARD );

        // Fourier transform working memmory
        std::vector<std::complex<realT>> tform1( psd2.size() );
        std::vector<std::complex<realT>> tform2( psd2.size() );

        std::vector<std::complex<realT>> Ntform1( psd2.size() );
        std::vector<std::complex<realT>> Ntform2( psd2.size() );

        // Normally distributed random numbers
        math::normDistT<realT> normVar;
        normVar.seed();

        // The two modal amplitude series
        std::vector<realT> fm_n( psd2.size() );
        std::vector<realT> fm_nm( psd2.size() );

        // The two noise amplitude series
        std::vector<realT> N_n( psd2.size() );
        std::vector<realT> N_nm( psd2.size() );

        // The speckle amplitude
        std::vector<realT> vnl( Nwd );
        std::vector<realT> vn( Nsamp );

        // Periodogram averager
        sigproc::averagePeriodogram<realT> avgPgram( Nsamp * 0.1, dt ); //, 0, 1);
        avgPgram.win( sigproc::window::hann );

        // The temporary periodogram
        std::vector<realT> tpgram( avgPgram.size() ); // spPSD.size());

#pragma omp for
        for( int k = 0; k < N; ++k )
        {
            // Generate the time-series
            for( int i = 0; i < psd2.size(); ++i )
            {
                fm_n[i] = normVar;
                N_n[i] = normVar;
                N_nm[i] = normVar;
            }

            // Filter and normalize the fourier mode time series
            filt( fm_n );
            math::vectorMeanSub( fm_n );
            realT actvar = math::vectorVariance( fm_n );
            realT norm = sqrt( fmVar / actvar );
            for( size_t q = 0; q < fm_n.size(); ++q )
                fm_n[q] *= norm;

            // And move it to the Fourier domain
            fft( tform1.data(), fm_n.data() );

            // Filter and normalize the measurement mode time series
            nfilt.filter( N_n );
            nfilt.filter( N_nm );

            realT Nactvar = 0.5 * ( math::vectorVariance( N_n ) + math::vectorVariance( N_nm ) );
            norm = sqrt( nVar / Nactvar );
            for( size_t q = 0; q < fm_n.size(); ++q )
                N_n[q] *= norm;
            for( size_t q = 0; q < fm_n.size(); ++q )
                N_nm[q] *= norm;

            // And move them to the Fourier domain
            fft( Ntform1.data(), N_n.data() );
            fft( Ntform2.data(), N_nm.data() );

            std::complex<realT> scale = std::complex<realT>( ( tform1.size() ), 0 );

            // Apply the modal phase shift, and apply the measurement noise.
            for( size_t m = 0; m < tform1.size(); ++m )
            {
                // Apply the phase shift to form the 2nd time series
                tform2[m] = tform1[m] * exp( std::complex<realT>( 0, math::half_pi<realT>() ) );

                // Apply the augmented ETF to two time-series
                tform1[m] *= xfer2[m] / scale;
                tform2[m] *= xfer2[m] / scale;

                // Ntform2[m] = Ntform1[m]*exp( std::complex<realT>(0, math::half_pi<realT>() ));

                Ntform1[m] *= nxfer2[m] / scale;
                Ntform2[m] *= nxfer2[m] / scale;
            }

            //<<<<<<<<****** Transform back to the time domain.
            fftB( fm_n.data(), tform1.data() );
            fftB( fm_nm.data(), tform2.data() );
            fftB( N_n.data(), Ntform1.data() );
            fftB( N_nm.data(), Ntform2.data() );

            // Calculate the speckle amplitude and mean-subtract
            realT mn = 0;
            for( int i = 0; i < Nwd; ++i )
            {
                realT h1 = fm_n[i + NwdStart] + N_n[i + NwdStart];
                realT h2 = fm_nm[i + NwdStart] + N_nm[i + NwdStart];

                vnl[i] = ( pow( h1, 2 ) + pow( h2, 2 ) );
                mn += vnl[i];
            }
            mn /= vnl.size();

            // Extract middle sample and mean subtract
            for( int i = 0; i < Nsamp; ++i )
                vn[i] = vnl[i + NsampStart] - mn;

            // Calculate PSD of the speckle amplitude
            if( !noPSD )
                avgPgram( tpgram, vn );

// Accumulate
#pragma omp critical
            {
                if( bins != nullptr && vars != nullptr )
                    math::vectorBinMeans( means, *bins, vn );
                for( size_t i = 0; i < spPSD.size(); ++i )
                    spPSD[i] += tpgram[i];
            }
        }
    } // pragma omp parallel

    if( !noPSD )
    {
        //-- Calculate frequency grid (which is maybe not identical to input freq)
        realT df = ( 1.0 ) / ( 2.0 * spPSD.size() * dt );
        spFreq.resize( spPSD.size() );
        for( size_t i = 0; i < spFreq.size(); ++i )
        {
            spFreq[i] = globalAvgPgram[i];
            spPSD[i] /= N;
        }
        // realT spVar = sigproc::psdVar1sided(df, spPSD.data(), spPSD.size());
        // for(size_t i=0; i< spPSD.size(); ++i) spPSD[i] *= (fmVar*fmVar/spVar);
    }

    // Calculate binned variances.
    if( bins != nullptr && vars != nullptr )
    {
        vars->resize( bins->size() );
        for( size_t i = 0; i < bins->size(); ++i )
        {
            ( *vars )[i] = math::vectorVariance( means[i] );
        }
    }

    return 0;
}

/// Calculate the variance of the mean vs bin size in a speckle time-series
/**
 * Does so by generating a time-series from the PSD using Fourier-domain convolution, then
 * calculates the binned-variances in the generated time-series.
 *
 * This is just a wrapper for speckleAmpPSD, with no actual PSD calculation.
 *
 */
template <typename realT>
int speckleAmpVarMean(
    std::vector<realT> &vars,                    ///< [out] The binned variances of the time series generated.
    std::vector<realT> &bins,                    ///< [in] The bin sizes to use to calculate the variances.
    std::vector<realT> &freq,                    ///< [in] The Frequency grid of the input PSD
    std::vector<realT> &fmPSD,                   ///< [in] The open-loop Fourier mode PSD.
    std::vector<std::complex<realT>> &fmXferFxn, ///< [in] The complex error transfer function, as a function of freq
    std::vector<realT> &nPSD,                    ///< [in] The open-loop noise PSD
    std::vector<std::complex<realT>> &nXferFxn,  ///< [in] The noise transfer function, as a function of freq
    int N ///< [in] The number of trials to use in calculating the amplitude time-series.
)
{
    std::vector<realT> spFreq, spPSD;

    return speckleAmpPSD( spFreq, spPSD, freq, fmPSD, fmXferFxn, nPSD, nXferFxn, N, &vars, &bins, true );
}

} // namespace analysis
} // namespace AO
} // namespace mx

#endif // speckleAmpPSD_hpp
