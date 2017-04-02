

#ifndef __systemSetup_hpp__
#define __systemSetup_hpp__


template<typename aosysT, typename floatT>
void setupMagAOX(aosysT & aosys, floatT F0, floatT F1, floatT lam_wfs, floatT lam_sci)
{
   aosys.load_MagAO_model();
   aosys.atm.L_0(25.);
   
   aosys.d_min(aosys.D()/48.);
   aosys.npix_wfs(120*120);
   aosys.ron_wfs(0.3);//0.3;
   
   aosys.minTauWFS(1./3630.);
   aosys.deltaTau(1.5*aosys.minTauWFS());
   
   aosys.lam_wfs(lam_wfs);
   aosys.lam_sci(lam_sci);

   aosys.ncp_wfe(30e-9);
   
   
   aosys.F0(F0);

   
}

#endif //__systemSetup_hpp__

