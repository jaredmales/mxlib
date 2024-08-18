

#include <iostream>

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

int main()
{
    typedef cpp_dec_float_100 floatT;
    // typedef double floatT;
    std::cout.precision( 100 );

    floatT tan_arcsec;

    tan_arcsec =
        tan( static_cast<floatT>( 1.0 ) * pi<floatT>() / static_cast<floatT>( 180 ) / static_cast<floatT>( 3600 ) );

    floatT twosqrt2log2 =
        static_cast<floatT>( 2.0 ) * sqrt( static_cast<floatT>( 2.0 ) * log( static_cast<floatT>( 2.0 ) ) );

    std::cout << "              "; // label buffer
    for( int i = 0; i < 10; ++i )  // print digit scale
    {
        std::cout << i;
        for( int j = 0; j < 9; ++j )
            std::cout << ' ';
    }
    std::cout << "1\n";

    std::cout << "tan_arcsec:   " << tan_arcsec << "\n";
    std::cout << "twosqrt2log2: " << twosqrt2log2 << "\n";
    std::cout << "ln_two:       " << ln_two<floatT>() << "\n";
}
