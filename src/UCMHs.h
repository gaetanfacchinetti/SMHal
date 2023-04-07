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

namespace ucmhs {

void test_PBH(); // For the module on PBHs
void test_CosmologicalBoost(); // For the module on the Cosmological boost
void plot_fractionInUCMHs();
void plot_formation_redshift();
void plot_AverageRedshiftOfCollapse_OnSpike();

int _main()
{   

    //plot_AverageRedshiftOfCollapse_OnSpike();
    test_PBH();
    exit(0);

    Cosmology cosmo(CosmoModel::PLANCKONLY);
    //PowerSpectrum power_spectrum(cosmo, PSWindowType::real_space_top_hat);

    std::vector<double> As = {1e-6};
    std::vector<double> ks = {1e+6};
    std::vector<double> eps = {0}; 

    std::cout << cosmo.keq_rad_mat() << std::endl;
    std::cout << cosmo.get_Omega_m_h2() << std::endl;
    std::shared_ptr<TransferFunction_EH98> TF =  std::make_shared<TransferFunction_EH98>(cosmo, false);
    TF->plot();

    //plot_formation_redshift();
    double D0 = cosmo.growth_factor_D1_Carroll(0);
    double Om_m_H2 = 1e+4*cosmo.get_Omega_m_h2() / pow(C_LIGHT * 1e-3, 2) ; // in Mpc^{-2}
    double kT = 6.2e-2; // Mpc^{-1}
    double Sigma_star = 4./25. * pow(0.23*D0*kT*kT/Om_m_H2, 2);
    std::cout << cosmo.get_h() << std::endl;
    std::cout << Om_m_H2 << " " << D0 << std::endl;
    std::cout << Sigma_star << std::endl;
    std::cout << pow(1.686 * 3401 * D0, 2)/Sigma_star << std::endl;
    
    PowerSpectrum power_spectrum(cosmo, PSWindowType::fourier_space_top_hat, false); 
    power_spectrum.plot_Sigma_vs_RM(0, false);

    //power_spectrum.set_spikes(As, ks, eps); 
    
    //test_PBH();

    //MassFunction PrSc(power_spectrum);
    //std::cout << power_spectrum.Mass_vs_LagrangianRadius(1/1e+3) << std::endl;

    //power_spectrum.plot_Sigma_vs_RM(0, false, false, "ks6_As6");
    //std::cout << cosmo.critical_density(0) << std::endl;
    //plot_fractionInUCMHs();
    //std::cout << my_erfc(21) << std::endl;
    //test_CosmologicalBoost();

    //cosmo._Inverse_growth_factor_Caroll(0.01);
    //exit(0);

    //std::cout << cosmo.Inverse_growth_factor_Caroll(0.5) << " " << cosmo._Inverse_growth_factor_Caroll(0.5) << std::endl;

    return 0;
}

void test_PBH()
{
    Cosmology cosmo(CosmoModel::PLANCKONLY);
    PowerSpectrum power_spectrum(cosmo, PSWindowType::real_space_top_hat, false); 
    DegreeFreedom *degf = new DegreeFreedom(T_QCD_GeV);
    PrimordialBlackHoles PBH(power_spectrum, degf);

    std::cout << PBH.Tformation_vs_MassHorizon(1e+7) << std::endl;
    std::cout << PBH.Radius_vs_MassHorizon(10) << " " << 1/cosmo.keq_rad_mat()*1e+3 << std::endl;
    PBH.plot_MassHorizon_vs_Radius();
    PBH.plot_fPBH(1e-1, 1e+7);
    std::cout << PBH.Max_fPBH(0.2, 1e-1, 1e+7) << std::endl;
    std::cout << PBH.Max_As(1e-11, 0.2, 1e+4) << std::endl;
    PBH.plot_Max_As();
    double Mpcm1_to_GeV = 1e-3*m_to_kpc/GeV_to_mm1;
    std::cout << 4/3*PI*PLANCK_MASS*pow(T0_CMB_GeV, 2.)*pow(3.93, 2./3.)*pow(106.75, -1./6.)*pow(0.05, -2.)*pow(PI/80, 1./2.)*GeV_to_MSOL*pow(Mpcm1_to_GeV, -2.) << std::endl;
    std::cout << PBH.fPBH(100, 0.2, 1e-2, 1e+5) << std::endl;
    
}

void plot_fractionInUCMHs()
{
    Cosmology cosmo(CosmoModel::PLANCKONLY);
    PowerSpectrum power_spectrum(cosmo, PSWindowType::fourier_space_top_hat, false); 
    //DegreeFreedom *degf = new DegreeFreedom(T_QCD_GeV);

    
    std::vector<double> As = {1e-6};
    std::vector<double> ks = {1e+3};
    std::vector<double> eps = {0}; 

    power_spectrum.set_spikes(As, ks, eps); 
    //power_spectrum.plot_Sigma_vs_RM(0, false, false, "no_spike"); 
    //cosmo.plot_Der_growth_factor();
    //std::cout << power_spectrum.Mass_vs_LagrangianRadius(1e-3) << std::endl;

    //MassConcentrationModel MC(power_spectrum);

    std::shared_ptr<MassFunction_PS> mass_function = std::make_shared<MassFunction_PS>(power_spectrum);
    CosmologicalBoost cosmo_boost(mass_function);

    // Plotting the mass fraction in the spike
    std::ofstream outfile;
    outfile.open("../output/UCMHS/MassFractionInSpike.out");
    double ks_min = 1, ks_max = 1e+6; // Mpc^{-1}
    double T1_2 = pow(power_spectrum.get_TransferFunction().evaluate(ks[0]), 2);
    double T1_2_ks_min = pow(power_spectrum.get_TransferFunction().evaluate(ks_min), 2);
    double T1_2_ks_max = pow(power_spectrum.get_TransferFunction().evaluate(ks_max), 2);
    double As_min = 1e-10, As_max = 1e-2;
    //double As_eff_min = As_min*pow(ks_min, 4)*T1_2_ks_min/pow(ks[0], 4)/T1_2;
    double As_eff_min = As_min;
    //double As_eff_max = As_max*pow(ks_max, 4)*T1_2_ks_max/pow(ks[0], 4)/T1_2;

    int N = 100;
    //double dAs_eff = log10(As_eff_max/As_eff_min)/(1.0*N), As_eff = 0;
    double dAs_eff = log10(As_max/As_min)/(1.0*N), As_eff = 0;
    double z_min = 1e-2, z_max = 1e+4, dz = log10(z_max/z_min)/(1.0*N), z = 0;


    for(int i = 0; i < N+1; i++)
    {
        As_eff = As_eff_min*pow(10, i*dAs_eff);
        //std::cout << As_eff << std::endl;
        PowerSpectrum PS_bis(cosmo, PSWindowType::fourier_space_top_hat, false);
        As[0] = As_eff;
        PS_bis.set_spikes(As, ks, eps);
        std::shared_ptr<MassFunction_PS> mass_function_bis = std::make_shared<MassFunction_PS>(PS_bis);
        double ms = PS_bis.Mass_vs_LagrangianRadius(1./ks[0]); 
        std::cout << i << " " << As[0] << " " << ks[0] << " " << ms << " " << mass_function_bis->get_PowerSpectrum().get_As_spike()[0] << std::endl;
        CosmologicalBoost CB_bis(mass_function_bis);

        //double h = cosmo.get_h();
        //double D1_0 = cosmo.growth_factor_D1_Carroll(0);
        //double Omega_m_0 = cosmo.get_Omega_m_h2() / (h * h);
        //double c_over_H0_Mpc = 1e-3 * LIGHT_SPEED / (100 * h);

        for(int j = 0; j < N+1; j++)
        {
            z = z_min*pow(10, j*dz);
            //outfile << (4./25.)*As_eff*T1_2*pow(ks[0], 4)*pow(D1_0 * c_over_H0_Mpc * c_over_H0_Mpc / Omega_m_0, 2)  << "\t" << z << "\t" << CB_bis.MassFractionInSpike(z, 1e-6) << std::endl;
            outfile << As_eff << "\t" << z << "\t" << CB_bis.MassFractionInSpike(z, 1e-6) << std::endl;
        }
    }

    outfile.close();
    //std::cout << As_eff_min << " " << As_eff_max << std::endl;
}


void test_CosmologicalBoost()
{
    Cosmology cosmo(CosmoModel::PLANCKONLY);
    std::shared_ptr<TransferFunction_WIMP_GHS05> TF_WIMP = TransferFunction_WIMP_GHS05::TransferFunction_WIMP_GHS05_from_Tempmchi(cosmo, 1e-5, 1e-2);
    //std::shared_ptr<TransferFunction_WIMP_GHS05> TF_WIMP = TransferFunction_WIMP_GHS05::TransferFunction_WIMP_GHS05_from_Tempmchi(cosmo, 28e-3, 1e+2);
    //PowerSpectrum power_spectrum(cosmo, PSWindowType::fourier_space_top_hat, false, TF_WIMP); 
    PowerSpectrum power_spectrum(cosmo, PSWindowType::fourier_space_top_hat, false); 

    std::vector<double> As = {8e-4};
    std::vector<double> ks = {286};
    std::vector<double> eps = {0}; 

    //power_spectrum.set_spikes(As, ks, eps); 


    std::shared_ptr<MassFunction_ST> mass_function = std::make_shared<MassFunction_ST>(power_spectrum, 0.73, 0.353, 0.175); // Values used in 21cmFAST
    //std::shared_ptr<MassFunction_PS> PrSc = std::make_shared<MassFunction_PS>(power_spectrum);
    
    //MassConcentrationModel MC = MassConcentrationModel::MassConcentrationModel_from_PowerSpectrum(power_spectrum);
    //double ms = power_spectrum.Mass_vs_LagrangianRadius(1./ks[0]); 
    //std::cout << ms << std::endl;
    //std::cout << cosmo.redshift_horizon_entry(1e-2) << std::endl;
    //std::cout << power_spectrum.SigmaM(0.99*ms, 0) << std::endl;
    //std::cout << cosmo.cosmic_abundance(0, Species::MATTER) << std::endl;
    //std::cout << cosmo.cosmological_density(0, Species::MATTER)<< std::endl;
    //std::cout << cosmo.critical_density(0)/pow(cosmo.get_h(), 2)*1e+9 << std::endl;
    //std::cout << cosmo.get_Omega_m_h2()/pow(cosmo.get_h(), 2) << std::endl;
    //std::cout << MC.RedshiftOfCollapse_vs_Mass(0.99*ms) << std::endl;
    //power_spectrum.plot_Sigma_vs_RM(0, true, false, "another_test");
    //std::cout << mass_function->NumberDensity(1e-10, 0) << std::endl;
    //exit(0);
    mass_function->plot(0, false, "comp");
    exit(0);

    //PrSc.plot_HaloMassFunction(0,  false, "Tkd28em3_m100");
    //PrSc.plot_HaloMassFunction(10, false, "Tkd28em3_m100");
    //PrSc.plot_HaloMassFunction(20, false, "Tkd28em3_m100");
    //PrSc.plot_HaloMassFunction(30, false, "Tkd28em3_m100");

    CosmologicalBoost CB(mass_function);
  
    //double ms = power_spectrum.Mass_vs_LagrangianRadius(1. / ks[0]);
    //std::cout << CB.HaloDistribution(ms, 1000, 0, 1) << std::endl;
    //CB.plot_FormationRedshiftDistribution("ks3_As8_EST", 1);
    //CB.plot_FormationRedshiftDistribution("ks3_As8_BBKS", 2);
    
    CB.plot_TotalBoost(1e-6, 0, "ST_fsth_Maccio_Tkd1em5_m1em2");

    exit(0);
}


void plot_formation_redshift()
{
    Cosmology cosmo(CosmoModel::PLANCKONLY);
    PowerSpectrum power_spectrum(cosmo, PSWindowType::fourier_space_top_hat, false); 


    double ks_min = 1e+1, ks_max = 1e+10, nks =100, ks;
    double dks = log10(ks_max/ks_min)/(1.0*nks);
    double As_min = 1e-9, As_max = 1e-1, nAs =100, As;
    double dAs = log10(As_max/As_min)/(1.0*nAs);


    std::ofstream outfile;
    outfile.open("../output/UCMHS/Formation_redshift_vs_spike.out");

    double dc = 1.686;
    double D0 = cosmo.growth_factor_D1_Carroll(0);
    
    double S, zf; 

    outfile << 0 << "\t";


    for(int j = 0; j < nAs+1; j++)
    {
        As = As_min*pow(10, j*dAs);
        outfile << As << "\t";
    }

    outfile << std::endl;

    for(int i = 0; i < nks+1; i++)
    {
        ks = ks_min*pow(10, i*dks);
        outfile << ks << "\t";

        for(int j = 0; j < nAs+1; j++)
        {
            As = As_min*pow(10, j*dAs);
            S = power_spectrum.SigmaR2_narrow_spike(0.99/ks, As, ks, 0);
            zf = cosmo.Inverse_growth_factor_Caroll(dc*D0/sqrt(S));

            outfile << zf << "\t";
        }

        outfile << std::endl;
    }



    outfile.close();
}


void plot_AverageRedshiftOfCollapse_OnSpike()
{
    Cosmology cosmo(CosmoModel::PLANCKONLY);
    PowerSpectrum power_spectrum(cosmo, PSWindowType::fourier_space_top_hat, false); 

    std::vector<double> As = {1e-8};
    std::vector<double> ks = {1e+3};
    std::vector<double> eps = {0}; 

    power_spectrum.set_spikes(As, ks, eps); 
    
    std::shared_ptr<MassFunction_PS> mass_function = std::make_shared<MassFunction_PS>(power_spectrum);
    CosmologicalBoost cosmological_boost(mass_function);  

    double As_min = 1e-9, As_max = 1e-4, nAs =20;
    double dAs = log10(As_max/As_min)/(1.0*nAs);
    double z_min = 1e-2, z_max = 1e+4, nz =101, z;
    double dz = log10(z_max/z_min)/(1.0*nz);
    


    std::ofstream outfile, outfile_zfmp, outfile_Nz;
    outfile.open("../output/UCMHs/AverageRedshiftOfCollapse_OnSpike_ks1e+3.out");
    outfile_zfmp.open("../output/UCMHs/Zfmp_ks1e+3.out");
    outfile_Nz.open("../output/UCMHs/Nz_ks1e+3.out");

    double dc = 1.686;
    double D0 = cosmo.growth_factor_D1_Carroll(0);
    
    outfile_zfmp << "# As | Zf(M_star-) | Zf(M_star+) " << std::endl;

    std::cout << power_spectrum.Mass_vs_LagrangianRadius(1/ks[0]) << " " << cosmo.cosmological_density(0, Species::MATTER) << " " << cosmo.critical_density(0) << std::endl;

    double S_m, S_p, zf_av, zf_m, zf_p, Nz; 

    outfile << 0 << "\t";
    outfile_Nz << 0 << "\t";


    for(int j = 0; j < nz+1; j++)
    {
        z = z_min*pow(10, j*dz);
        outfile << z << "\t";
        outfile_Nz << z << "\t";
    }

    outfile << std::endl;
    outfile_Nz << std::endl;

    for(int i = 0; i < nAs+1; i++)
    {
        As[0] = As_min*pow(10, i*dAs);
        outfile << As[0] << "\t";
        outfile_Nz << As[0] << "\t";

        power_spectrum.set_spikes(As, ks, eps);
        mass_function = std::make_shared<MassFunction_PS>(power_spectrum);
        cosmological_boost = CosmologicalBoost(mass_function);  
         

        S_m = power_spectrum.SigmaR2(0.99/ks[0], 0);
        S_p = power_spectrum.SigmaR2(1.01/ks[0], 0);

        zf_m = cosmo.Inverse_growth_factor_Caroll(dc*D0/sqrt(S_m));
        zf_p = cosmo.Inverse_growth_factor_Caroll(dc*D0/sqrt(S_p));

        outfile_zfmp << As[0] << "\t" << zf_m << "\t" << zf_p << std::endl;

        for(int j = 0; j < nz+1; j++)
        {
            z = z_min*pow(10, j*dz);
            zf_av = cosmological_boost.AverageRedshiftOfCollapse_OnSpike(z);
            Nz = mass_function->Der_MassFraction_spike(As[0], ks[0], z);
            outfile << zf_av << "\t";
            outfile_Nz << Nz << "\t";
        }

        outfile << std::endl;
        outfile_Nz << std::endl;
    }



    outfile.close();
    outfile_zfmp.close();
    outfile_Nz.close();
}




} // End of the namespace