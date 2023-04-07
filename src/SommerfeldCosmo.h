/*
 main.cpp
 
 Author : Gaetan Faclambnetti
 mail   : gaetanfacc@gmail.com
 
 */

#include "headers/mymath.h"
#include "headers/MyUnits.h"
#include "headers/Particle.h"
#include "headers/CrossSections/CrossSection_ccgaga.h"
#include "headers/CrossSections/CrossSection_ccff.h"
#include "headers/CrossSections/CrossSection_ccss.h"
#include "headers/CrossSections/CrossSection_ccpp.h"
#include "headers/CrossSections/CrossSection_ccsp.h"
#include "headers/CrossSections/CrossSectionMajorana_cccc.h"
#include "headers/CrossSections/CrossSectionDirac_cccc.h"
#include "headers/CrossSections/CrossSectionDirac_ccbccb.h"
#include "headers/CrossSections/CrossSection_cfcf.h"
#include "headers/SpecialIntegrations.h"
#include "headers/StandardModel.h"
#include "headers/DegreeFreedom.h"
#include "headers/ChemicalDecoupling.h"
#include "headers/KineticDecoupling.h"
#include "headers/MinimalMass.h"
#include "headers/CPOddHiggsLowMass.h"
#include "headers/InputReader.h"
#include "headers/DecayRates.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "headers/ExceptionHandler.h"
#include "headers/ConstraintCouplingUniversal.h"
#include "headers/PowerSpectrum.h"
#include "headers/LocalBoost.h"
#include <ctime>
#include "headers/DirectDetection.h"
#include "headers/PrimordialBlackHoles.h"
#include "headers/MassConcentrationModel.h"
#include "headers/CosmologicalBoost.h"
#include "headers/SommerfeldEnhancement.h"

namespace sommerfeldcosmo {

void test_SommerfeldCosmologicalBoost();

int _main()
{   
    SommerfeldEnhancement somm;
    //somm.plot(1e-6, 1e-2);
    test_SommerfeldCosmologicalBoost();
    return 0;
}

void test_SommerfeldCosmologicalBoost()
{
    Cosmology cosmo(CosmoModel::PLANCKONLY);
    std::shared_ptr<TransferFunction_WIMP_GHS05> TF_WIMP = TransferFunction_WIMP_GHS05::TransferFunction_WIMP_GHS05_from_Tempmchi(cosmo, 28e-3, 1e+2);
    PowerSpectrum power_spectrum(cosmo, PSWindowType::fourier_space_top_hat, false, TF_WIMP); 
    std::shared_ptr<MassFunction_ST> mass_function = std::make_shared<MassFunction_ST>(power_spectrum, 0.73, 0.353, 0.175);
    CosmologicalBoost cosmological_boost(mass_function);
    SommerfeldEnhancement sommerfeld_enhancement;
    DarkHalo dark_halo(cosmo, 0, 1, 3, 1);

    dark_halo.Initialise_from_virial_parameters(1e-5);
    std::cout << dark_halo.circular_velocity_DM_only(dark_halo.get_rs()) << std::endl;

   
    std::vector<double> ephi_res_arr = sommerfeld_enhancement.get_ephi_res_arr();
    std::vector<double> ephi_sat_arr = sommerfeld_enhancement.get_ephi_sat_arr();
    cosmological_boost.plot_TotalBoost_with_Sommerfeld_approx(1, ephi_sat_arr[ephi_sat_arr.size()-1], 0.01, 1e-20, "test_pwave_3");  
    
    //sommerfeld_enhancement.plot_vs_v(ephi_res_arr[ephi_res_arr.size()-20], 0.01);
    //std::cout << ephi_res_arr[ephi_res_arr.size()-20] << std::endl;
    
    //std::cout << cosmological_boost.TotalBoost_with_Sommerfeld_approx(10, ephi_res_arr[50], 1, 0.01, 1e-15) << std::endl;
    //std::cout << cosmological_boost.TotalBoost(20, 1e-15, ProfileConfig::ALL_NFW, 0) << std::endl;

}


} // End of namespace
