/** \file fourierTemporalPSD_test.cpp
 * \brief Tests of Fourier-mode temporal power spectral densities.
 */

#include "../../../catch2/catch.hpp"

#define MX_NO_ERROR_REPORTS

#include "../../../../include/ao/analysis/fourierTemporalPSD.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

using aoSystemBaseT = mx::AO::analysis::aoSystem<double, mx::AO::analysis::vonKarmanSpectrum<double>>;

/// AO-system test type that forces local template instantiation under sanitizers.
struct aoSystemT : public aoSystemBaseT
{
};

using temporalPsdT = mx::AO::analysis::fourierTemporalPSD<double, aoSystemT>;

/// Verify temporal-PSD tail initialization at and below its nominal averaging width.
/** Exercises zero, one, 49, and 50 exactly integrated frequency bins. */
TEST_CASE( "Fourier temporal PSD handles short exact tails", "[ao::analysis::fourierTemporalPSD]" )
{
    aoSystemT aoSystem;
    aoSystem.D( 6.5 );
    aoSystem.atm.setSingleLayer( 0.16, 500e-9, 25.0, 0.0, 0.0, 10.0, 0.0 );

    temporalPsdT temporalPsd;
    temporalPsd.m_aosys = &aoSystem;
    temporalPsd._useBasis = mx::AO::analysis::basis::basic;

    constexpr std::size_t frequencyCount = 55;
    std::vector<double> frequency( frequencyCount );
    for( std::size_t index = 0; index < frequency.size(); ++index )
    {
        frequency[index] = static_cast<double>( index + 1 );
    }

    for( const std::size_t exactCount : { 0U, 1U, 49U, 50U } )
    {
        DYNAMIC_SECTION( exactCount << " exact bins" )
        {
            std::vector<double> psd( frequency.size(), 0.0 );
            const double maximumExactFrequency = exactCount == 0 ? 0.5 : frequency[exactCount - 1];

            const int result = temporalPsd.singleLayerPSD( psd, frequency, 1.0, 0.0, 0, 1, maximumExactFrequency );

            if( exactCount == 0 )
            {
                REQUIRE( result < 0 );
                continue;
            }

            REQUIRE( result == 0 );

            const std::size_t averageCount = std::min<std::size_t>( exactCount, 50 );
            const double exponent = aoSystem.atm.alpha( 0 ) + 2.0;
            double expected = 0.0;
            for( std::size_t offset = averageCount; offset > 0; --offset )
            {
                const std::size_t index = exactCount - offset;
                expected += psd[index] * std::pow( frequency[index] / frequency[exactCount], exponent );
            }
            expected /= static_cast<double>( averageCount );

            REQUIRE( psd[exactCount] == Approx( expected ).epsilon( 1e-12 ) );
            REQUIRE( std::isfinite( psd[exactCount] ) );
            REQUIRE( psd[exactCount] >= 0.0 );

            const double expectedLast =
                psd[exactCount] * std::pow( frequency[exactCount] / frequency.back(), exponent );
            REQUIRE( psd.back() == Approx( expectedLast ).epsilon( 1e-12 ) );
        }
    }
}
