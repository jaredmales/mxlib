
#include "geo.h"






// template<typename imageTout, typename imageTin, typename coeffT>
// void cutImageRegion(imageTout imout, const imageTin & imin,  coeffT coeffs, bool resize = true)
// {
//    if(resize)
//    {
//       imout.resize(coeffs.size());
//    }
//    
//    #pragma omp parallel for schedule(static, 1)
//    for(int i=0;i<coeffs.size();++i)
//    {
//       imout(i) = imin(coeffs[i]);
//    }
//    
// }
//  
// template<typename imageTout, typename imageTin, typename coeffT>
// void insertImageRegion(imageTout imout, const imageTin & imin,  coeffT coeffs)
// {
//    
//    #pragma omp parallel for schedule(static, 1)
//    for(int i=0;i<coeffs.size();++i)
//    {
//       imout(coeffs[i]) = imin(i);
//    }
//    
// } 
// 
// 



