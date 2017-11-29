


#include <iostream>

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;


int main()
{
   typedef cpp_dec_float_100 floatT;
   //typedef double floatT;
   std::cout.precision(50);
   
   floatT tan_arcsec;
         
   tan_arcsec = tan( static_cast<floatT>(1.0) * pi<floatT>()/static_cast<floatT>(180)/static_cast<floatT>(3600));
   
   for(int i=0;i<5;++i)
   {
      std::cout << i;
      for(int j=0;j<9;++j) std::cout << ' ';
   }
   std::cout << "5\n";
   
   std::cout << tan_arcsec << "\n";
}
