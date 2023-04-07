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

namespace simplifiedmodels {

void write_input_MinimalMassHalos_vs_massDM_correctAbundance(std::string const &name_input, Proptype type, int Npts, double massMin, double massMax, double ratio);
void write_input_MinimalMassHalos_vs_massDM_vs_massProp_correctAbundance(std::string const &name_input, Proptype type, int Npts, double massMin, double massMax);

/// Plot the annihilation cross-sections DM into scalar/pseudo-scalars 
void plot_CrossSections_Annihilation_PhiPhi_universalCouplings();
/// Plot the self inteaction cross-section
void plot_CrossSections_SelfInteraction();
/// Plot the annihilation cross-section DM into SM fermions with index "part"
void plot_CrossSections_Annihilation_SMFermions(int part);
/// Plot the annihilation cross-sections DM into all SM fermions for universal couplings
void plot_CrossSections_Annihilation_FermionFermion_universalCouplings();
/// Plot the decay rates of photons or gluons into SM fermions
void plot_decay_photons_gluons_fermions();
/// Function to evaluate the value of the (s- or p-wave) cross-section giving the right abundance
void plot_Sigma_Vs_MassDM(Cosmology Cosmo, DegreeFreedom *degFree, bool is_swave);
/// Bijection method for the function to evaluate the value of the (s- or p-wave) cross section giving the right abundance
double Dichotomie_Omega(DegreeFreedom *degFree, double Omega_h2_obs, double massDM, bool is_swave);


/// Run the chain of analysis for a given mediator and a given coupling type
void run_analysis(std::string mediator, CouplingType ctype);

/// Test the resolution of chemical decoupling
void test_calculation_cross_section();

//
//
//
//
//
//


int _main()
{   

    run_analysis("Scalar", CouplingType::UNIVERSAL_YUKAWA);
    //test_calculation_cross_section();

    exit(0);

    DarkSectorModel model_s, model_ps;
    model_s.add_darkmatter(1, 0, Fermiontype::majorana);
    model_ps.add_darkmatter(1, 0, Fermiontype::majorana);
    model_s.add_propagator(Proptype::scalar, 1, "phi_s", 10, 0);
    model_ps.add_propagator(Proptype::pseudoscalar, 1, "phi_p", 0.1, 0);

    model_s.Initialise();
    model_ps.Initialise();


    double lambda_SM = 0.1; // This value does not matter

    model_s.ForceInitCouplings();  // Function that initiate the size of the coupling vectors in model
    model_ps.ForceInitCouplings(); // Function that initiate the size of the coupling vectors in model

    double mass = 0;
    double mt = M_QUARK_TOP;

    Cosmology Cosmo(CosmoModel::PLANCKONLY);
    ConstraintCouplingUniversal ccu(CouplingType::UNIVERSAL_YUKAWA, Cosmo);

    ccu.SetNewUniversalYukawaCoupling(model_s, 0.1, 0.1);
    DecayRates dc;
    dc.set_width_phis(model_s);

    std::cout << model_s.get_phis(0).get_width()/model_s.get_phis(0).get_mass() << std::endl;

    return 0;
}


void run_analysis(std::string mediator, CouplingType ctype)
{
    //ccu.Write_input_MinimalMassHalos_vs_massDM_vs_massProp("PseudoScalar_test", Proptype::pseudoscalar, 20, 0.01, 1e+4);
    Cosmology Cosmo(CosmoModel::PLANCKONLY);
    ConstraintCouplingUniversal ccu(ctype, Cosmo);

    std::string str_ctype = "";
    if (ctype == CouplingType::UNIVERSAL_YUKAWA)
        str_ctype = "Yukawa";
    
    //ccu.Coupling_correctAbundance(mediator, "Coupling_" + mediator + "_" +  str_ctype);
    ccu.TkdAndMinMass(mediator, "Coupling_" + mediator + "_" + str_ctype, "Mass_" + mediator + "_" + str_ctype);
    //ccu.UnevolvedNsubhalos("Mass_" + mediator + "_" + str_ctype, "NumberSubhalos_" + mediator + "_" + str_ctype);
    //ccu.SelfInteractingCS(mediator, "Coupling_" + mediator + "_" + str_ctype, "CS_SI_" + mediator + "_" + str_ctype);
    //ccu.DecayTime(mediator, "Coupling_" + mediator + "_" + str_ctype, "DecayTime_" + mediator + "_" + str_ctype);
}


void plot_Sigma_Vs_MassDM(Cosmology Cosmo, DegreeFreedom *degFree, bool is_swave)
{

    double massMin = 0.1;
    double massMax = 10000;
    // double SigmvMin = 1.5e-26 / conv;
    // double SigmvMax = 5.2e-26 / conv;
    double SigmvMin = 1e-26 / GeVm2_to_cm3persec;
    double SigmvMax = 1e-25 / GeVm2_to_cm3persec;
    double nPointsMass = 101;
    double nPointsSigmv = 201;

    double deltaLogMass = abs(log10(massMax) - log10(massMin)) / (1.0 * nPointsMass);
    double deltaSigmv = (SigmvMax - SigmvMin) / (1.0 * nPointsSigmv);

    double logMass, Sigmv;
    double Y0;

    std::ofstream outfile;
    if (is_swave == true)
        outfile.open("../output/Y0_Mass_SigmV_swave.out");
    else
        outfile.open("../output/Y0_Mass_SigmV_pwave.out");
    outfile.precision(8);
    outfile << "# mass DM [GeV] | sigmav [cm^3 s^{-1}] | 1 sigma constraint - | 1 sigma constraint + | 3 sigma constraint | 3 sigma constraint + | 5 sigma constraint - | 5 sigma constraint + " << std::endl;

    double sigv = 0, sigv_1sigma_p = 0, sigv_1sigma_m = 0, sigv_3sigma_p = 0, sigv_3sigma_m = 0, sigv_5sigma_p = 0, sigv_5sigma_m = 0;

    // VWARNING : Value for the Planck-Only 18 cosmology column 6 table 2 (everything without BAO) in https://arxiv.org/abs/1807.06209
    double eps = 1.2e-3;

    for (int i = 0; i < nPointsMass + 1; i++)
    {
        logMass = log10(massMin) + deltaLogMass * i;
        std::cout << i << std::endl;

        sigv = Dichotomie_Omega(degFree, Cosmo.get_Omega_c_h2(), pow(10, logMass), is_swave);
        sigv_1sigma_p = Dichotomie_Omega(degFree, Cosmo.get_Omega_c_h2() + eps, pow(10, logMass), is_swave);
        sigv_1sigma_m = Dichotomie_Omega(degFree, Cosmo.get_Omega_c_h2() - eps, pow(10, logMass), is_swave);
        sigv_3sigma_p = Dichotomie_Omega(degFree, Cosmo.get_Omega_c_h2() + 3 * eps, pow(10, logMass), is_swave);
        sigv_3sigma_m = Dichotomie_Omega(degFree, Cosmo.get_Omega_c_h2() - 3 * eps, pow(10, logMass), is_swave);
        sigv_5sigma_p = Dichotomie_Omega(degFree, Cosmo.get_Omega_c_h2() + 5 * eps, pow(10, logMass), is_swave);
        sigv_5sigma_m = Dichotomie_Omega(degFree, Cosmo.get_Omega_c_h2() - 5 * eps, pow(10, logMass), is_swave);
        outfile << pow(10, logMass) << "\t" << sigv << "\t" << sigv_1sigma_m << "\t" << sigv_1sigma_p << "\t" << sigv_3sigma_m << "\t" << sigv_3sigma_p << "\t" << sigv_5sigma_m << "\t" << sigv_5sigma_p << std::endl;
    }

    outfile.close();
}

double Dichotomie_Omega(DegreeFreedom *degFree, double Omega_h2_obs, double massDM, bool is_swave)
//massDM in GeV
{
    DarkSectorModel DS;
    DS.add_darkmatter(massDM, 0, Fermiontype::majorana);

    CS_TYPE cs;

    if (is_swave == true)
        cs = CS_TYPE::S_WAVE_ONLY;
    else
        cs = CS_TYPE::P_WAVE_ONLY;

    double Sigmv_0 = 0, Sigmv_1 = 0, Sigmv_2 = 0;
    double Y0_0 = 0, Y0_1 = 0, Y0_2 = 0;
    double SigvMax = 1e-23 / GeVm2_to_cm3persec, SigvMin = 1e-28 / GeVm2_to_cm3persec;
    double Omega_0 = 0, Omega_1 = 0, Omega_2 = 0;
    std::vector<double> xY;

    int compt = 0;

    double coeffOmega_h2 = 16 * pow(PI, 3) / 135 * 3.93773 * pow(T0_CMB_GeV, 3) * pow(100 / kpc_to_m / GeV_to_sm1, -2) * pow(M_PLANCK, -2);

    ChemicalDecoupling *chemDec0, *chemDec1, *chemDec2;

    while (abs(Omega_1 - Omega_h2_obs) / Omega_h2_obs > 1e-4 && compt < 1e+3)
    {
        if (compt == 0)
        {
            Sigmv_0 = SigvMin;
            Sigmv_2 = SigvMax;

            chemDec0 = new ChemicalDecoupling(DS, degFree, cs, Sigmv_0);
            xY = chemDec0->SolveFreezeOut();
            Y0_0 = chemDec0->ComputeAbundanceToday(xY);
            Omega_0 = Y0_0 * coeffOmega_h2 * massDM;
            //std::cout << "xY " << xY[1] << " " << Y0_0 << std::endl;

            chemDec2 = new ChemicalDecoupling(DS, degFree, cs, Sigmv_2);
            xY = chemDec2->SolveFreezeOut();
            Y0_2 = chemDec2->ComputeAbundanceToday(xY);
            Omega_2 = Y0_2 * coeffOmega_h2 * massDM;
        }

        Sigmv_1 = sqrt(Sigmv_0 * Sigmv_2);

        chemDec1 = new ChemicalDecoupling(DS, degFree, cs, Sigmv_1);
        xY = chemDec1->SolveFreezeOut();
        Y0_1 = chemDec1->ComputeAbundanceToday(xY);
        Omega_1 = Y0_1 * coeffOmega_h2 * massDM;

        if ((Omega_0 - Omega_h2_obs) * (Omega_1 - Omega_h2_obs) <= 0 && (Omega_2 - Omega_h2_obs) * (Omega_1 - Omega_h2_obs) > 0)
        {
            Sigmv_2 = Sigmv_1;
            Omega_2 = Omega_1;
        }
        else if ((Omega_0 - Omega_h2_obs) * (Omega_1 - Omega_h2_obs) > 0 && (Omega_2 - Omega_h2_obs) * (Omega_1 - Omega_h2_obs) <= 0)
        {
            Sigmv_0 = Sigmv_1;
            Omega_0 = Omega_1;
        }
        else
        {
            std::cout << "FATAL ERROR : " << __PRETTY_FUNCTION__ << std::endl;
            exit(0);
        }

        compt++;
    }

    return Sigmv_1 * GeVm2_to_cm3persec;
}

void test_chemical_decoupling()
{
    DegreeFreedom *degFree = new DegreeFreedom(T_QCD_GeV);
    std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/Simple_CrossSections_phiphi/Scalar_electron.in");
    // Kill the reader singleton to free memory and because we do not need it anymore for the moment
    InputReader::kill();
    DecayRates dc;

    for (int i = 0; i < models.size(); i++)
    {
        dc.set_width_phip(models[i]);
        dc.set_width_phis(models[i]);
        //dc.set_width_X(models[i]);
    }

    std::vector<double> xY = {0, 0};
    double Y0, Omega_h2_comp_0, Omega_h2_comp_1;

    //std::cout << "Coucou" << std::endl;
    ChemicalDecoupling *_ChemDec = new ChemicalDecoupling(models[0], degFree, CS_TYPE::FULL_ONLY);
    //std::cout << "Ici" << std::endl;

    double mass_chi = models[0].get_chi(0).get_mass();
    double coeffOmega_h2 = 16 * pow(PI, 3) / 135 * 3.93773 * pow(T0_CMB_GeV, 3) * pow(100 / kpc_to_m / GeV_to_sm1, -2) * pow(M_PLANCK, -2);

    _ChemDec->set_saveYToPrint(true);

    try
    {
        xY = _ChemDec->SolveFreezeOut();
        Y0 = _ChemDec->ComputeAbundanceToday(xY);

        Omega_h2_comp_0 = Y0 * coeffOmega_h2 * mass_chi;
        Omega_h2_comp_1 = xY[1] * coeffOmega_h2 * mass_chi;

        std::cout << Omega_h2_comp_0 << "\t" << Omega_h2_comp_1 << std::endl;
        //_ChemDec->SigmaV(100./5.);
        //_ChemDec->plotSigmaV();
        //std::cout << "Res freeze out " << _ChemDec->SolveFreezeOut()[0] << std::endl;
        //_ChemDec->PlotYFO();
    }
    catch (ExceptionHandler &e)
    {
        std::cout << "here" << std::endl;
        std::cout << e.what() << " : line " << e.get_line() << std::endl;
    }

    //_ChemDec->plotFuncToIntegrateSigmaV(5e-3);
    delete _ChemDec;
}

void plot_CrossSections_Annihilation_FermionFermion_universalCouplings()
{
    std::ofstream outfile;
    outfile.open("../output/CrossSectionAnn/CrossSectionAnn_Universal_couplings.out");
    outfile.precision(8);

    outfile << "# mchi = 100 GeV, mprop = 500 GeV and couplings = 0.1" << std::endl;
    outfile << "# s [GeV^2] | sigma [GeV^{-2}]" << std::endl;

    // Read the input models
    std::vector<DarkSectorModel> model_scalar = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/Scalar_universal.in");
    std::vector<DarkSectorModel> model_pseudoscalar = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/PseudoScalar_universal.in");
    std::vector<DarkSectorModel> model_vector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/Vector_universal.in");
    std::vector<DarkSectorModel> model_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/AxialVector_universal.in");

    std::vector<DarkSectorModel> model_scalar_pseudoscalar = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/ScalarPseudoScalar_universal.in");
    std::vector<DarkSectorModel> model_scalar_vector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/ScalarVector_universal.in");
    std::vector<DarkSectorModel> model_scalar_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/ScalarAxialVector_universal.in");
    std::vector<DarkSectorModel> model_pseudoscalar_vector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/PseudoScalarVector_universal.in");
    std::vector<DarkSectorModel> model_pseudoscalar_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/PseudoScalarAxialVector_universal.in");
    std::vector<DarkSectorModel> model_vector_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/VectorAxialVector_universal.in");
    InputReader::kill();

    //std::cout << "Dans le modele : " << model_scalar_axialvector[0].get_b_VEC_SM(0, 3) << std::endl;

    DecayRates dc;
    //dc.set_width_phis();
    dc.set_width_phis(model_scalar[0]);
    dc.set_width_phip(model_pseudoscalar[0]);
    dc.set_width_X(model_vector[0]);
    dc.set_width_X(model_axialvector[0]);

    dc.set_width_phis(model_scalar_pseudoscalar[0]);
    dc.set_width_phip(model_scalar_pseudoscalar[0]);
    dc.set_width_phis(model_scalar_vector[0]);
    dc.set_width_X(model_scalar_vector[0]);
    dc.set_width_phis(model_scalar_axialvector[0]);
    dc.set_width_X(model_scalar_axialvector[0]);

    dc.set_width_phip(model_pseudoscalar_vector[0]);
    dc.set_width_X(model_pseudoscalar_vector[0]);
    dc.set_width_phip(model_pseudoscalar_axialvector[0]);
    dc.set_width_X(model_pseudoscalar_axialvector[0]);
    dc.set_width_X(model_vector_axialvector[0]);

    std::cout << model_scalar[0].get_phis_ptr(0)->get_width() / model_scalar[0].get_phis_ptr(0)->get_mass() << std::endl;
    std::cout << model_pseudoscalar[0].get_phip_ptr(0)->get_width() / model_pseudoscalar[0].get_phip_ptr(0)->get_mass() << std::endl;
    std::cout << model_vector[0].get_X_ptr(0)->get_width() / model_vector[0].get_X_ptr(0)->get_mass() << std::endl;
    std::cout << model_axialvector[0].get_X_ptr(0)->get_width() / model_axialvector[0].get_X_ptr(0)->get_mass() << std::endl;

    //std::cout << model_vector[0].get_a_VEC_DM(0)[0][0] << std::endl;

    int n_crossSections_ccff = StandardModel::getInstance()->get_n_couples_SM_ferm();
    std::vector<CrossSection *> crossSection_SM_scalar, crossSection_SM_pseudoscalar, crossSection_SM_axialvector, crossSection_SM_vector;
    std::vector<CrossSection *> crossSection_SM_scalar_pseudoscalar, crossSection_SM_scalar_vector, crossSection_SM_scalar_axialvector;
    std::vector<CrossSection *> crossSection_SM_pseudoscalar_vector, crossSection_SM_pseudoscalar_axialvector, crossSection_SM_vector_axialvector;

    crossSection_SM_scalar.resize(0);
    crossSection_SM_pseudoscalar.resize(0);
    crossSection_SM_axialvector.resize(0);
    crossSection_SM_vector.resize(0);

    crossSection_SM_scalar_pseudoscalar.resize(0);
    crossSection_SM_scalar_vector.resize(0);
    crossSection_SM_scalar_axialvector.resize(0);
    crossSection_SM_pseudoscalar_vector.resize(0);
    crossSection_SM_pseudoscalar_axialvector.resize(0);
    crossSection_SM_vector_axialvector.resize(0);

    for (int i = 0; i < n_crossSections_ccff; i++)
    {
        crossSection_SM_scalar.push_back(new CrossSection_ccff(model_scalar[0], model_scalar[0].get_chi(0), model_scalar[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));
        crossSection_SM_pseudoscalar.push_back(new CrossSection_ccff(model_pseudoscalar[0], model_pseudoscalar[0].get_chi(0), model_pseudoscalar[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));
        crossSection_SM_vector.push_back(new CrossSection_ccff(model_vector[0], model_vector[0].get_chi(0), model_vector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));
        crossSection_SM_axialvector.push_back(new CrossSection_ccff(model_axialvector[0], model_axialvector[0].get_chi(0), model_axialvector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));

        crossSection_SM_scalar_pseudoscalar.push_back(new CrossSection_ccff(model_scalar_pseudoscalar[0], model_scalar_pseudoscalar[0].get_chi(0), model_scalar_pseudoscalar[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));
        crossSection_SM_scalar_vector.push_back(new CrossSection_ccff(model_scalar_vector[0], model_scalar_vector[0].get_chi(0), model_scalar_vector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));
        crossSection_SM_scalar_axialvector.push_back(new CrossSection_ccff(model_scalar_axialvector[0], model_scalar_axialvector[0].get_chi(0), model_scalar_axialvector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));

        crossSection_SM_pseudoscalar_vector.push_back(new CrossSection_ccff(model_pseudoscalar_vector[0], model_pseudoscalar_vector[0].get_chi(0), model_pseudoscalar_vector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));
        crossSection_SM_pseudoscalar_axialvector.push_back(new CrossSection_ccff(model_pseudoscalar_axialvector[0], model_pseudoscalar_axialvector[0].get_chi(0), model_pseudoscalar_axialvector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));
        crossSection_SM_vector_axialvector.push_back(new CrossSection_ccff(model_vector_axialvector[0], model_vector_axialvector[0].get_chi(0), model_vector_axialvector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(i)));
    }

    //std::cout << "Dans le modele : " << model_scalar_vector[0].get_b_VEC_SM(0, 0) << std::endl;
    //crossSection_SM_scalar_axialvector[6]->EvalSigmaInterfChannelS(50000);

    //exit(0);

    int n_crossSections = crossSection_SM_scalar.size();

    /*
    double mchi = 100; // GeV
    double smin = 4 * mchi * mchi;
    double smax = 4 * mchi * mchi * 1e+3;
    double Npts = 5000;
    double ds = log10(smax / smin) / Npts;
    double s;*/

    double vrel_min = 0;
    double vrel_max = 1;
    double Npts = 5000;
    double dvrel = (vrel_max - vrel_min) / Npts;
    double vrel;

    double res_scalar = 0;
    double res_pseudoscalar = 0;
    double res_vector = 0;
    double res_axialvector = 0;

    double res_scalar_pseudoscalar = 0;
    double res_scalar_vector = 0;
    double res_scalar_axialvector = 0;

    double res_pseudoscalar_vector = 0;
    double res_pseudoscalar_axialvector = 0;
    double res_vector_axialvector = 0;

    for (int i = 0; i < Npts + 1; i++)
    {
        //s = smin * pow(10, i * ds);
        vrel = vrel_min + i * dvrel;

        res_scalar = 0;
        res_pseudoscalar = 0;
        res_axialvector = 0;
        res_vector = 0;

        res_scalar_pseudoscalar = 0;
        res_scalar_vector = 0;
        res_scalar_axialvector = 0;

        res_pseudoscalar_vector = 0;
        res_pseudoscalar_axialvector = 0;
        res_vector_axialvector = 0;

        for (int j = 0; j < n_crossSections; j++)
        {
            res_scalar += crossSection_SM_scalar[j]->EvalSigma_fVrel(vrel) * vrel;
            res_pseudoscalar += crossSection_SM_pseudoscalar[j]->EvalSigma_fVrel(vrel) * vrel;
            res_vector += crossSection_SM_vector[j]->EvalSigma_fVrel(vrel) * vrel;
            res_axialvector += crossSection_SM_axialvector[j]->EvalSigma_fVrel(vrel) * vrel;

            res_scalar_pseudoscalar += crossSection_SM_scalar_pseudoscalar[j]->EvalSigma_fVrel(vrel) * vrel;
            res_scalar_vector += crossSection_SM_scalar_vector[j]->EvalSigma_fVrel(vrel) * vrel;
            res_scalar_axialvector += crossSection_SM_scalar_axialvector[j]->EvalSigma_fVrel(vrel) * vrel;

            res_pseudoscalar_vector += crossSection_SM_pseudoscalar_vector[j]->EvalSigma_fVrel(vrel) * vrel;
            res_pseudoscalar_axialvector += crossSection_SM_pseudoscalar_axialvector[j]->EvalSigma_fVrel(vrel) * vrel;
            res_vector_axialvector += crossSection_SM_vector_axialvector[j]->EvalSigma_fVrel(vrel) * vrel;
        }

        outfile << vrel << "\t" << res_scalar << "\t" << res_pseudoscalar << "\t" << res_vector << "\t" << res_axialvector << "\t"
                << res_scalar_pseudoscalar << "\t" << res_scalar_vector << "\t" << res_scalar_axialvector << "\t"
                << res_pseudoscalar_vector << "\t" << res_pseudoscalar_axialvector << "\t" << res_vector_axialvector << std::endl;
    }

    outfile.close();
}

// part = 3 for electron positrons
void plot_CrossSections_Annihilation_SMFermions(int part)
{
    std::ofstream outfile;
    std::string filename = "../output/CrossSectionAnn/CrossSectionAnn";

    if (part == 3)
        filename += "_epem";

    filename += ".out";

    outfile.open(filename);
    outfile.precision(8);

    // Read the input models
    std::vector<DarkSectorModel> model_scalar = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/Scalar_universal.in");
    std::vector<DarkSectorModel> model_pseudoscalar = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/PseudoScalar_universal.in");
    std::vector<DarkSectorModel> model_vector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/Vector_universal.in");
    std::vector<DarkSectorModel> model_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/AxialVector_universal.in");

    std::vector<DarkSectorModel> model_scalar_pseudoscalar = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/ScalarPseudoScalar_universal.in");
    std::vector<DarkSectorModel> model_scalar_vector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/ScalarVector_universal.in");
    std::vector<DarkSectorModel> model_scalar_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/ScalarAxialVector_universal.in");
    std::vector<DarkSectorModel> model_pseudoscalar_vector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/PseudoScalarVector_universal.in");
    std::vector<DarkSectorModel> model_pseudoscalar_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/PseudoScalarAxialVector_universal.in");
    std::vector<DarkSectorModel> model_vector_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/VectorAxialVector_universal.in");
    InputReader::kill();

    //std::cout << "Dans le modele : " << model_scalar_axialvector[0].get_b_VEC_SM(0, 3) << std::endl;

    DecayRates dc;
    //dc.set_width_phis();
    dc.set_width_phis(model_scalar[0]);
    dc.set_width_phip(model_pseudoscalar[0]);
    dc.set_width_X(model_vector[0]);
    dc.set_width_X(model_axialvector[0]);

    dc.set_width_phis(model_scalar_pseudoscalar[0]);
    dc.set_width_phip(model_scalar_pseudoscalar[0]);
    dc.set_width_phis(model_scalar_vector[0]);
    dc.set_width_X(model_scalar_vector[0]);
    dc.set_width_phis(model_scalar_axialvector[0]);
    dc.set_width_X(model_scalar_axialvector[0]);

    dc.set_width_phip(model_pseudoscalar_vector[0]);
    dc.set_width_X(model_pseudoscalar_vector[0]);
    dc.set_width_phip(model_pseudoscalar_axialvector[0]);
    dc.set_width_X(model_pseudoscalar_axialvector[0]);
    dc.set_width_X(model_vector_axialvector[0]);

    std::cout << model_scalar[0].get_phis_ptr(0)->get_width() / model_scalar[0].get_phis_ptr(0)->get_mass() << std::endl;
    std::cout << model_pseudoscalar[0].get_phip_ptr(0)->get_width() / model_pseudoscalar[0].get_phip_ptr(0)->get_mass() << std::endl;
    std::cout << model_vector[0].get_X_ptr(0)->get_width() / model_vector[0].get_X_ptr(0)->get_mass() << std::endl;
    std::cout << model_axialvector[0].get_X_ptr(0)->get_width() / model_axialvector[0].get_X_ptr(0)->get_mass() << std::endl;

    //std::cout << model_vector[0].get_a_VEC_DM(0)[0][0] << std::endl;

    int n_crossSections_ccff = StandardModel::getInstance()->get_n_couples_SM_ferm();
    CrossSection *crossSection_SM_scalar, *crossSection_SM_pseudoscalar, *crossSection_SM_axialvector, *crossSection_SM_vector;
    CrossSection *crossSection_SM_scalar_pseudoscalar, *crossSection_SM_scalar_vector, *crossSection_SM_scalar_axialvector;
    CrossSection *crossSection_SM_pseudoscalar_vector, *crossSection_SM_pseudoscalar_axialvector, *crossSection_SM_vector_axialvector;

    crossSection_SM_scalar = new CrossSection_ccff(model_scalar[0], model_scalar[0].get_chi(0), model_scalar[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));
    crossSection_SM_pseudoscalar = new CrossSection_ccff(model_pseudoscalar[0], model_pseudoscalar[0].get_chi(0), model_pseudoscalar[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));
    crossSection_SM_vector = new CrossSection_ccff(model_vector[0], model_vector[0].get_chi(0), model_vector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));
    crossSection_SM_axialvector = new CrossSection_ccff(model_axialvector[0], model_axialvector[0].get_chi(0), model_axialvector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));

    crossSection_SM_scalar_pseudoscalar = new CrossSection_ccff(model_scalar_pseudoscalar[0], model_scalar_pseudoscalar[0].get_chi(0), model_scalar_pseudoscalar[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));
    crossSection_SM_scalar_vector = new CrossSection_ccff(model_scalar_vector[0], model_scalar_vector[0].get_chi(0), model_scalar_vector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));
    crossSection_SM_scalar_axialvector = new CrossSection_ccff(model_scalar_axialvector[0], model_scalar_axialvector[0].get_chi(0), model_scalar_axialvector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));

    crossSection_SM_pseudoscalar_vector = new CrossSection_ccff(model_pseudoscalar_vector[0], model_pseudoscalar_vector[0].get_chi(0), model_pseudoscalar_vector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));
    crossSection_SM_pseudoscalar_axialvector = new CrossSection_ccff(model_pseudoscalar_axialvector[0], model_pseudoscalar_axialvector[0].get_chi(0), model_pseudoscalar_axialvector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));
    crossSection_SM_vector_axialvector = new CrossSection_ccff(model_vector_axialvector[0], model_vector_axialvector[0].get_chi(0), model_vector_axialvector[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(part));

    //std::cout << "Dans le modele : " << model_scalar_vector[0].get_b_VEC_SM(0, 0) << std::endl;
    //crossSection_SM_scalar_axialvector[6]->EvalSigmaInterfChannelS(50000);

    /*
    double mchi = 100; // GeV
    double smin = 4 * mchi * mchi;
    double smax = 4 * mchi * mchi * 1e+3;
    double Npts = 5000;
    double ds = log10(smax / smin) / Npts;
    double s;*/

    double vrel_min = 0;
    double vrel_max = 1;
    double Npts = 5000;
    double dvrel = (vrel_max - vrel_min) / Npts;
    double vrel;

    double res_scalar = 0;
    double res_pseudoscalar = 0;
    double res_vector = 0;
    double res_axialvector = 0;

    double res_scalar_pseudoscalar = 0;
    double res_scalar_vector = 0;
    double res_scalar_axialvector = 0;

    double res_pseudoscalar_vector = 0;
    double res_pseudoscalar_axialvector = 0;
    double res_vector_axialvector = 0;

    outfile << "# mchi = 100 GeV, mprop = 224.14 GeV and couplings = 0.1" << std::endl;
    outfile << "# vrel | vrel * sigma [cm^3 s^{-1}]" << std::endl;

    double conv = GeVm2_to_cm3persec;

    for (int i = 0; i < Npts + 1; i++)
    {
        //s = smin * pow(10, i * ds);
        vrel = vrel_min + i * dvrel;

        res_scalar = crossSection_SM_scalar->EvalSigma_fVrel(vrel) * vrel;
        res_pseudoscalar = crossSection_SM_pseudoscalar->EvalSigma_fVrel(vrel) * vrel;
        res_vector = crossSection_SM_vector->EvalSigma_fVrel(vrel) * vrel;
        res_axialvector = crossSection_SM_axialvector->EvalSigma_fVrel(vrel) * vrel;

        res_scalar_pseudoscalar = crossSection_SM_scalar_pseudoscalar->EvalSigma_fVrel(vrel) * vrel;
        res_scalar_vector = crossSection_SM_scalar_vector->EvalSigma_fVrel(vrel) * vrel;
        res_scalar_axialvector = crossSection_SM_scalar_axialvector->EvalSigma_fVrel(vrel) * vrel;

        res_pseudoscalar_vector = crossSection_SM_pseudoscalar_vector->EvalSigma_fVrel(vrel) * vrel;
        res_pseudoscalar_axialvector = crossSection_SM_pseudoscalar_axialvector->EvalSigma_fVrel(vrel) * vrel;
        res_vector_axialvector = crossSection_SM_vector_axialvector->EvalSigma_fVrel(vrel) * vrel;

        outfile << vrel << "\t" << res_scalar * conv << "\t" << res_pseudoscalar * conv << "\t" << res_vector * conv << "\t" << res_axialvector * conv << "\t"
                << res_scalar_pseudoscalar * conv << "\t" << res_scalar_vector * conv << "\t" << res_scalar_axialvector * conv << "\t"
                << res_pseudoscalar_vector * conv << "\t" << res_pseudoscalar_axialvector * conv << "\t" << res_vector_axialvector * conv << std::endl;
    }

    outfile.close();
}

void plot_CrossSections_Annihilation_PhiPhi_universalCouplings()
{
    std::ofstream outfile;
    outfile.open("../output/CrossSectionAnn/CrossSectionAnn_phiphi_Universal_couplings.out");
    outfile.precision(8);

    outfile << "# mchi = 100 GeV, mprop = 500 GeV and couplings = 0.1" << std::endl;
    outfile << "# vrel | vrel * sigma [cm^3 s^{-1}]" << std::endl;

    // Read the input models
    std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/Simple_CrossSections_phiphi/ScalarPseudoScalar_universal.in");

    InputReader::kill();

    // NOTE THAT HERE WE CONSIDER 0 WIDTHS

    // Define the particle before (as otherwise they will not be the same in the code)
    Particle phis_0 = models[0].get_phis(0), phis_1 = models[1].get_phis(0);
    Particle phip_0 = models[0].get_phip(0), phip_1 = models[1].get_phip(0);
    Particle chi_0 = models[0].get_chi(0), chi_1 = models[1].get_chi(0);

    CrossSection_ccss *c_ccss_m = new CrossSection_ccss(models[0], chi_0, chi_0, phis_0, phis_0);
    CrossSection_ccss *c_ccss_d = new CrossSection_ccss(models[1], chi_1, chi_1, phis_1, phis_1);
    CrossSection_ccpp *c_ccpp_m = new CrossSection_ccpp(models[0], chi_0, chi_0, phip_0, phip_0);
    CrossSection_ccpp *c_ccpp_d = new CrossSection_ccpp(models[1], chi_1, chi_1, phip_1, phip_1);
    CrossSection_ccsp *c_ccsp_m = new CrossSection_ccsp(models[0], chi_0, chi_0, phis_0, phip_0);
    CrossSection_ccsp *c_ccsp_d = new CrossSection_ccsp(models[1], chi_1, chi_1, phis_1, phip_1);

    double conv = GeVm2_to_cm3persec;

    double mchi = 500; // GeV
    double vrel_min = 0;
    double vrel_max = 1;
    double Npts = 5000;
    double dvrel = (vrel_max - vrel_min) / Npts;
    double vrel;

    for (int i = 0; i < Npts + 1; i++)
    {
        vrel = vrel_min + i * dvrel; // min * pow(10, i * ds);

        outfile << vrel << "\t" << c_ccss_m->EvalSigma_fVrel(vrel) * vrel * conv
                << "\t" << c_ccss_d->EvalSigma_fVrel(vrel) * vrel * conv
                << "\t" << c_ccpp_m->EvalSigma_fVrel(vrel) * vrel * conv
                << "\t" << c_ccpp_d->EvalSigma_fVrel(vrel) * vrel * conv
                << "\t" << c_ccsp_m->EvalSigma_fVrel(vrel) * vrel * conv
                << "\t" << c_ccsp_d->EvalSigma_fVrel(vrel) * vrel * conv << std::endl;
    }

    outfile.close();
}

void plot_CrossSections_SelfInteraction()
{
    std::ofstream outfile;
    outfile.open("../output/CrossSectionScatt/CrossSection_SelfInteraction.out");
    outfile.precision(8);

    // Read the input models
    std::vector<DarkSectorModel> model_scalar = InputReader::getInstance()->Read("../input/Simple_CrossSectionsSI/Scalar_SI.in");
    std::vector<DarkSectorModel> model_pseudoscalar = InputReader::getInstance()->Read("../input/Simple_CrossSectionsSI/PseudoScalar_SI.in");
    std::vector<DarkSectorModel> model_vector = InputReader::getInstance()->Read("../input/Simple_CrossSectionsSI/Vector_SI.in");
    std::vector<DarkSectorModel> model_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSectionsSI/AxialVector_SI.in");

    InputReader::kill();

    //std::cout << "Dans le modele : " << model_scalar_axialvector[0].get_b_VEC_SM(0, 3) << std::endl;

    DecayRates dc;
    //dc.set_width_phis();
    dc.set_width_phis(model_scalar[0]);
    dc.set_width_phip(model_pseudoscalar[0]);
    dc.set_width_X(model_vector[0]);
    dc.set_width_X(model_axialvector[0]);

    //std::cout << model_scalar[0].get_phis_ptr(0)->get_width() / model_scalar[0].get_phis_ptr(0)->get_mass() << std::endl;
    //std::cout << model_pseudoscalar[0].get_phip_ptr(0)->get_width() / model_pseudoscalar[0].get_phip_ptr(0)->get_mass() << std::endl;
    //std::cout << model_vector[0].get_X_ptr(0)->get_width() / model_vector[0].get_X_ptr(0)->get_mass() << std::endl;
    //std::cout << model_axialvector[0].get_X_ptr(0)->get_width() / model_axialvector[0].get_X_ptr(0)->get_mass() << std::endl;

    //std::cout << model_vector[0].get_a_VEC_DM(0)[0][0] << std::endl;

    CrossSectionMajorana_cccc *CS_M_scalar = new CrossSectionMajorana_cccc(model_scalar[0], model_scalar[0].get_chi(0));
    CrossSectionMajorana_cccc *CS_M_pseudoscalar = new CrossSectionMajorana_cccc(model_pseudoscalar[0], model_pseudoscalar[0].get_chi(0));
    CrossSectionMajorana_cccc *CS_M_axialvector = new CrossSectionMajorana_cccc(model_axialvector[0], model_axialvector[0].get_chi(0));

    // By default DM particles are Majorana, here I set it to Direc by hand and recompute the decay rates
    model_scalar[0].get_chi_ptr(0)->set_ftype(Fermiontype::dirac);
    model_pseudoscalar[0].get_chi_ptr(0)->set_ftype(Fermiontype::dirac);
    model_vector[0].get_chi_ptr(0)->set_ftype(Fermiontype::dirac);
    model_axialvector[0].get_chi_ptr(0)->set_ftype(Fermiontype::dirac);
    dc.set_width_phis(model_scalar[0]);
    dc.set_width_phip(model_pseudoscalar[0]);
    dc.set_width_X(model_vector[0]);
    dc.set_width_X(model_axialvector[0]);

    CrossSectionDirac_cccc *CS_D_scalar = new CrossSectionDirac_cccc(model_scalar[0], model_scalar[0].get_chi(0));
    CrossSectionDirac_cccc *CS_D_pseudoscalar = new CrossSectionDirac_cccc(model_pseudoscalar[0], model_pseudoscalar[0].get_chi(0));
    CrossSectionDirac_cccc *CS_D_vector = new CrossSectionDirac_cccc(model_vector[0], model_vector[0].get_chi(0));
    CrossSectionDirac_cccc *CS_D_axialvector = new CrossSectionDirac_cccc(model_axialvector[0], model_axialvector[0].get_chi(0));

    CrossSectionDirac_ccbccb *CS_D_CB_scalar = new CrossSectionDirac_ccbccb(model_scalar[0], model_scalar[0].get_chi(0), model_scalar[0].get_chi(0));
    CrossSectionDirac_ccbccb *CS_D_CB_pseudoscalar = new CrossSectionDirac_ccbccb(model_pseudoscalar[0], model_pseudoscalar[0].get_chi(0), model_pseudoscalar[0].get_chi(0));
    CrossSectionDirac_ccbccb *CS_D_CB_vector = new CrossSectionDirac_ccbccb(model_vector[0], model_vector[0].get_chi(0), model_vector[0].get_chi(0));
    CrossSectionDirac_ccbccb *CS_D_CB_axialvector = new CrossSectionDirac_ccbccb(model_axialvector[0], model_axialvector[0].get_chi(0), model_axialvector[0].get_chi(0));

    outfile << "#" << std::endl;
    outfile << "# vrel [km/s] | sigma [GeV^{-2}]" << std::endl;

    double c_kms = C_LIGHT * 1e-3; // Light speed in km/s

    double vrel_min = 10;   // km/s
    double vrel_max = 5000; // kms/s
    double Npts = 100;
    double dvrel = log10(vrel_max / vrel_min) / Npts;
    double vrel;

    double m0 = model_scalar[0].get_chi(0).get_mass();
    double m1 = model_pseudoscalar[0].get_chi(0).get_mass();
    double m2 = model_axialvector[0].get_chi(0).get_mass();
    double m3 = model_vector[0].get_chi(0).get_mass();
    //std::cout << m0 << " " << m1 << " " << m2 << " " << m3 << std::endl;
    double conv = pow(GeV_to_cmm1, -2) * pow(GeV_to_kg * 1e+3, -1);

    //std::vector<DarkSectorModel> test = InputReader::getInstance()->Read("../input/Simple_CrossSections_fermferm/Scalar_universal.in");
    //dc.set_width_phis(test[0]);
    //std::cout << StandardModel::getInstance()->get_couples_SM_ferm(4).get_particle(0).get_name() << std::endl;
    //std::cout << test[0].get_phis_ptr(0)->get_width() / test[0].get_phis_ptr(0)->get_mass() << std::endl;
    //CrossSection_cfcf *CS_test = new CrossSection_cfcf(test[0], test[0].get_chi(0), StandardModel::getInstance()->get_couples_SM_ferm(4));
    //model_scalar[0].get_phis_ptr(0)->set_width(50);
    //std::cout << model_scalar[0].get_phis_ptr(0)->get_width() << std::endl;
    //CrossSectionDirac_cccc *CS_D_scalar_2 = new CrossSectionDirac_cccc(model_scalar[0], model_scalar[0].get_chi(0));
    //model_scalar[0].get_phis_ptr(0)->set_mass(100);
    //std::cout << "res = " <<  CS_D_scalar_2->EvalSigmaTransfer_fVrel(0.9) << " " << CS_D_scalar_2->EvalSigma_fVrel(0.9) << std::endl;

    //double smin = 4e+4, smax = 4e+6, ds = log10(smax/smin)/100.;
    //for(int i = 0; i < 101; i++)
    //    std::cout << smin*pow(10, i*ds)<< "\t" << CS_D_scalar->EvalSigmaTransfer(smin*pow(10, i*ds)) << std::endl;

    //exit(0);

    for (int i = 0; i < Npts + 1; i++)
    {
        //s = smin * pow(10, i * ds);
        vrel = vrel_min * pow(10, i * dvrel);

        //std::cout << vrel/c_kms << std::endl;

        outfile << vrel << "\t"
                << vrel * CS_M_scalar->EvalSigmaTransfer_fVrel(vrel / c_kms) / m0 * conv << "\t"
                << vrel * CS_M_pseudoscalar->EvalSigmaTransfer_fVrel(vrel / c_kms) / m1 * conv << "\t"
                << vrel * CS_M_axialvector->EvalSigmaTransfer_fVrel(vrel / c_kms) / m2 * conv << "\t"
                << vrel * CS_D_scalar->EvalSigmaTransfer_fVrel(vrel / c_kms) / m0 * conv << "\t"
                << vrel * CS_D_pseudoscalar->EvalSigmaTransfer_fVrel(vrel / c_kms) / m1 * conv << "\t"
                << vrel * CS_D_axialvector->EvalSigmaTransfer_fVrel(vrel / c_kms) / m2 * conv << "\t"
                << vrel * CS_D_vector->EvalSigmaTransfer_fVrel(vrel / c_kms) / m3 * conv << "\t"
                << vrel * CS_D_CB_scalar->EvalSigmaTransfer_fVrel(vrel / c_kms) / m0 * conv << "\t"
                << vrel * CS_D_CB_pseudoscalar->EvalSigmaTransfer_fVrel(vrel / c_kms) / m1 * conv << "\t"
                << vrel * CS_D_CB_axialvector->EvalSigmaTransfer_fVrel(vrel / c_kms) / m2 * conv << "\t"
                << vrel * CS_D_CB_vector->EvalSigmaTransfer_fVrel(vrel / c_kms) / m3 * conv << "\t"
                << std::endl;
    }

    outfile.close();
}

void write_input_MinimalMassHalos_vs_massDM_correctAbundance(std::string const &name_input, Proptype type, int Npts, double massMin, double massMax, double ratio)
{
    std::ofstream out_file;
    out_file.open("../input/MassHalos/" + name_input + ".in");

    double dl_m = log10(massMax / massMin) / (Npts - 1);

    time_t now = time(0);

    // convert now to string form
    char *dt = ctime(&now);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);

    out_file << "####################################################" << std::endl;
    out_file << "#" << std::endl;
    out_file << "# DMEvolution input file" << std::endl;
    out_file << "# Automatically generated from main" << std::endl;
    out_file << "#" << std::endl;
    out_file << "# Gaetan Facchinetti" << std::endl;
    out_file << "# gaetan.facchinetti@umontpellier.fr" << std::endl;
    out_file << "# Laboratoire Univers et Particules Montpellier" << std::endl;
    out_file << "# Created on : (UTC) " << dt;
    out_file << "#" << std::endl;
    out_file << "####################################################" << std::endl;
    out_file << std::endl;

    out_file << "## Number of points to evaluate ##" << std::endl;
    out_file << "Npoints :: " << Npts << std::endl;
    out_file << std::endl;
    out_file << std::endl;

    out_file << "## Dark sector particles ##" << std::endl;
    out_file << "DSCONT :: NDM :: 1" << std::endl;
    if (type == Proptype::scalar)
        out_file << "DSCONT :: DSPart :: (scalar, 1)" << std::endl;
    else if (type == Proptype::pseudoscalar)
        out_file << "DSCONT :: DSPart :: (pseudoscalar, 1)" << std::endl;
    else if (type == Proptype::vector)
        out_file << "DSCONT :: DSPart :: (vector, 3)" << std::endl;
    out_file << std::endl;
    out_file << std::endl;

    out_file << "## Masses dark sector [GeV] ##" << std::endl;

    double massDM;

    for (int i = 0; i < Npts; i++)
    {
        massDM = massMin * pow(10, i * dl_m);
        out_file << "MDMDS :: " << i << " :: " << massDM << ", " << ratio * massDM << std::endl;
    }

    out_file << std::endl;
    out_file << std::endl;

    out_file << "## Interactions DS-SM [GeV] ##" << std::endl;

    for (int i = 0; i < Npts; i++)
        out_file << "IntDSSM :: " << i << " :: 0 :: "
                 << "l_e=" << 0.3 << std::endl;

    out_file << std::endl;
    out_file << std::endl;

    out_file << "## Interactions DS-DM [GeV] ##" << std::endl;

    for (int i = 0; i < Npts; i++)
        out_file << "IntDSDM :: " << i << " :: 0 :: "
                 << "l_0_0=" << 0.3 << std::endl;
}

//

void plot_decay_photons_gluons_fermions()
// This function compares the decay rates into massless spin 1 vectors and fermions
{
    DarkSectorModel model_s, model_ps;
    model_s.add_propagator(Proptype::scalar, 1, "phi_s", 100, 0);
    model_ps.add_propagator(Proptype::pseudoscalar, 1, "phi_p", 100, 0);

    double lambda_SM = 0.1; // This value does not matter

    model_s.ForceInitCouplings();  // Function that initiate the size of the coupling vectors in model
    model_ps.ForceInitCouplings(); // Function that initiate the size of the coupling vectors in model

    double mass = 0;
    double mt = M_QUARK_TOP;

    for (int j = 0; j < StandardModel::getInstance()->get_n_couples_SM_ferm(); j++)
    {
        mass = StandardModel::getInstance()->get_couples_SM_ferm(j).get_particle(0).get_mass();
        model_s.set_coupling_DS_S_SM(0, j, lambda_SM * mass / mt);
        model_ps.set_coupling_DS_PS_SM(0, j, lambda_SM * mass / mt);
    }

    DecayRates dc;

    std::ofstream outfile;
    outfile.open("../output/decay_photons_gluons_fermions.out");

    int Npts = 500;
    double mmin = 1e-2, mmax = 1e+4, dm = log10(mmax / mmin) / (1. * Npts), m = 0;

    outfile << "# m [GeV] | decay scalar two photons [GeV^{-1}] | decay scalar two gluons [GeV^{-1}] | sum of decay scalar two fermions [GeV^{-1}] | ... " << std::endl;
    outfile << "# ... | decay pseudoscalar two photons [GeV^{-1}] | decay pseudoscalar two gluons [GeV^{-1}] | sum of decay pseudoscalar two fermions [GeV^{-1}] " << std::endl;

    //std::cout << fOneLoopPseudoScalarVectors(10) << std::endl;
    //std::cout << regulatorTQCD(log10(10), 1.0 / 0.3, log10(T_QCD_GeV), 1, 0) << std::endl;
    //model_s.get_phis_ptr(0)->set_mass(10.);
    //std::cout << dc.DecayRate_sgluongluon(model_s, 0)/dc.DecayRate_sff(model_s, 0)  << std::endl;
    //std::cout << abs(fOneLoopScalarVectors(4*pow(M_QUARK_TOP, 2)/100.) + fOneLoopScalarVectors(4*pow(M_QUARK_BOTTOM, 2)/100.)  +  fOneLoopScalarVectors(4*pow(M_QUARK_CHARM, 2)/100.)) << std::endl;

    //exit(0);

    for (int i = 0; i < Npts + 1; i++)
    {
        m = mmin * pow(10, i * dm);
        model_s.get_phis_ptr(0)->set_mass(m);
        model_ps.get_phip_ptr(0)->set_mass(m);
        //std::cout << "here : " << i << std::endl;
        outfile << m << "\t" << dc.DecayRate_sphotonphoton(model_s, 0) << "\t" << dc.DecayRate_sgluongluon(model_s, 0) << "\t" << dc.DecayRate_sff(model_s, 0)
                << "\t" << dc.DecayRate_pphotonphoton(model_ps, 0) << "\t" << dc.DecayRate_pgluongluon(model_ps, 0) << "\t" << dc.DecayRate_pff(model_ps, 0) << std::endl;
    }

    outfile.close();
}



void test_CrossSection()
{

    DarkSectorModel model_s, model_ps;
    model_s.add_darkmatter(1, 0, Fermiontype::majorana);
    model_ps.add_darkmatter(1, 0, Fermiontype::majorana);
    model_s.add_propagator(Proptype::scalar, 1, "phi_s", 0.1, 0);
    model_ps.add_propagator(Proptype::pseudoscalar, 1, "phi_p", 0.1, 0);

    model_s.Initialise();
    model_ps.Initialise();
    model_s.ForceInitCouplings();  // Function that initiate the size of the coupling vectors in model
    model_ps.ForceInitCouplings(); // Function that initiate the size of the coupling vectors in model

    double mass = 0;
    double mt = M_QUARK_TOP;
    double lambda = sqrt(4*PI);

    for (int j = 0; j < StandardModel::getInstance()->get_n_couples_SM_ferm(); j++)
    {
        mass = StandardModel::getInstance()->get_couples_SM_ferm(j).get_first().get_mass();
        model_s.set_coupling_DS_S_SM(0, j, lambda * mass / mt);
        model_ps.set_coupling_DS_PS_SM(0, j, lambda * mass / mt);
    }             

    model_s.set_coupling_DS_S_DM(0, 0, 0, lambda);
    model_ps.set_coupling_DS_PS_DM(0, 0, 0, lambda);


    DecayRates dc;
    dc.set_width_phis(model_s);
    dc.set_width_phip(model_ps);

    //std::cout << "here" << std::endl;


    Particle chi_s = model_s.get_chi(0);
    Particle chi_ps = model_ps.get_chi(0);

    CrossSection_ccff *cr_ccff = new CrossSection_ccff(model_s, chi_s, chi_s, StandardModel::getInstance()->get_couples_SM_ferm(3)); // electrons
    CrossSection_cfcf *cr_cfcf = new CrossSection_cfcf(model_s, chi_s,  StandardModel::getInstance()->get_couples_SM_ferm(3)); // electrons
    CrossSection_cgacga *cr_cgacga = new CrossSection_cgacga(model_ps, chi_ps); // photons

    // Check that these two quantities are the correct ones
    cr_ccff->plot(4e+7, "Test_comp_tt");
    cr_ccff->plot_fromIntegral(4e+7, "Test_comp_tt");
    cr_cfcf->plot_sigmaTransfer(4e+7, "Test_comp_tt");
    cr_cfcf->plot_sigmaTransfer_fromIntegral(4e+7, "Test_comp_tt");
    //cr_cgacga->plot_sigmaTransfer(1.0001e+4, "Test_photons", "s");
    //cr_cgacga->plot_sigmaTransfer_fromIntegral(1.0001e+4, "Test_photons", "s");

    
    DegreeFreedom *degFree = new DegreeFreedom(T_QCD_GeV);
    KineticDecoupling KinDec(model_ps, degFree, 0);

    //KinDec.plotGammaTot();
    //KinDec.plotGamma_1Loop_photons();
}




void test_calculation_cross_section()
{
    DarkSectorModel model;
    model.add_darkmatter(1e-2, 0, Fermiontype::majorana);
    model.add_propagator(Proptype::scalar, 1, "phi_s", 1e-2, 0);

    model.Initialise();
    model.ForceInitCouplings();  // Function that initiate the size of the coupling vectors in model

    ConstraintCouplingUniversal ccu(CouplingType::UNIVERSAL_YUKAWA, Cosmology(CosmoModel::PLANCKONLY));
    bool widthParticleTooLarge = ccu.SetCoupling(model, 0.00405903, 0.00405903, 3);

    DegreeFreedom *degFree = new DegreeFreedom(T_QCD_GeV, &model);

    // The chemical decoupling object is a new object for a new model - i.e. new DM mass in our case
    ChemicalDecoupling ChemDec(model, degFree, CS_TYPE::FULL_ONLY);
    double coeffOmega_h2 = 16 * pow(PI, 3) / 135 * 3.93773 * pow(T0_CMB_GeV, 3) * pow(100 / kpc_to_m / GeV_to_sm1, -2) * pow(M_PLANCK, -2);

    // Computation of the comoving number density today
    std::vector<double> xY = ChemDec.SolveFreezeOut();
    //ChemDec.plotSigmaV();
    //ChemDec.plotFuncToIntegrateSigmaV(0.78476*1e-1);
    //Particle chi = model.get_chi(0);
    //CrossSection_ccff *cs_ccff = new CrossSection_ccff(model, chi, chi, StandardModel::getInstance()->get_couples_SM_ferm(3));
    //cs_ccff->plot(30*0.78476, "Test");
    double xcd = xY[2]; // x decoupling
    //Y0 = xY[1];  // Assume no later evolution
    double Y0 = ChemDec.ComputeAbundanceToday(xY);
    double Omega_h2_comp_min = Y0 * coeffOmega_h2 * model.get_chi(0).get_mass();
    //std::vector<double> zx = {18., model.get_phis(0).get_mass()/(0.78476*1e-1)};
    //std::cout << ChemDec.fToIntegrateForSigmaV_noquarks(zx) << std::endl;

    std::cout << model.get_phis(0).get_width()/model.get_phis(0).get_mass() << std::endl;
    std::cout << Omega_h2_comp_min << " " << xcd << " " << xY[1] * coeffOmega_h2 * model.get_chi(0).get_mass() << std::endl;
}


}