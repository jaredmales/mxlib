
#include <iostream>
#include "randomt"

#include <random>
#include "../vmop/MMatrix"

int main()
{
   
   //mx::randomT<unsigned long, unsigned long, mx::rangen_ISO_portable_ulong> r;
   
   //(mx::randomT<long, long, mx::rangen_ISO_portable>
   //r.set_seed(126);
   
   //for(int i=0; i< 100; i++)   std::cout << r << std::endl;
   
   
   //std::mt19937_64 r(21947293);
   //std::normal_distribution<double> normd;
   
   mx::norm_distd normr;
   mx::exp_distd expr1, expr2, expr3;
   mx::lap_distd lap1;
   
   //expr2.distribution.param(static_cast<mx::exp_distd::randistT::param_type>(0.5));
   
   normr.engine().seed(859384);
   //expr1.generator.seed(21947293);
   //expr2.generator.seed(67829832);
   //expr3.generator.seed(84029);
   
   //mx::Vectord expr_sum(10000);
      
   double e1, e2, e3, l1;
   //double N = 1000.;
   for(int i=0; i< 10000; i++)
   {  
      
      e1 = expr1;
      l1 = lap1;
      //e2 = expr2;
      //e3 = expr3;
      
      std::cout << e1 << " " <<  l1 << "\n";
   }   
   
   return 0;
}
