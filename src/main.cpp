/*
 main.cpp
 
 Author : Gaetan Faclambnetti
 mail   : gaetanfacc@gmail.com
 
 */

#include "DDCalc.hpp"
#include "headers/DirectDetection.h"
#include "SimplifiedModels.h"
#include "UCMHs.h"
#include "SommerfeldCosmo.h"

// Initialisation of singleton classes
StandardModel *StandardModel::_SM = NULL;
InputReader *InputReader::_IR = NULL;

std::vector<double> DirectDetectionConstraintsOneModel(DarkSectorModel &model, double lambda, ConstraintCouplingUniversal ccu);
void DirectDetectionConstraints(std::string const &name_model_file_input, std::string const &name_coupling_file_input, std::string const &name_output, CouplingType cptype, Cosmology cosmo);

bool Verbose = false;
int NUM_THREADS = 1;


int main(int argc, char **argv)
{
    std::cout << std::endl;
    std::cout << " --------------------- SMHhal --------------------" << std::endl;
    std::cout << "| # Gaétan Facchinetti : gaetanfacc@gmail.com     |" << std::endl;
    std::cout << " --------------------------------------------------" << std::endl;
    std::cout << std::endl;

    // Reroute the log file
    std::ofstream outlogfile("../log/log.out");
    std::streambuf *prevclogbuf = std::clog.rdbuf(outlogfile.rdbuf());
    std::clog.rdbuf(outlogfile.rdbuf());

    // Intialisation of the Verbose parameter and the number of thread used
    for (int i = 0; i < argc; i++)
    {

        if (strncmp(argv[i], "-omp", 10) == 0)
        {
            //std::cout << argv[i+1] << std::endl;
            NUM_THREADS = atoi(argv[i + 1]);
        }

        if (strncmp(argv[i], "-v", 10) == 0)
        {
            Verbose = true;
            std::cout << "Verbose set to true !" << std::endl;
        }
            
    }

    // Indicate the total possible number of threads
    std::cout << "Maximal number of threads openMP is : " << omp_get_max_threads() << std::endl;
    std::cout << "Number of thread openMP used        : " << NUM_THREADS << std::endl;

    // gsl_set_error_handler_off();
    gsl_set_error_handler(&my_gsl_handler);


    //sommerfeldcosmo::_main();
    ucmhs::_main();

    exit(0);
   //Cosmology cosmo(CosmoModel::PLANCKONLY);
    //Galaxy gal(cosmo, 0, 0);
    //std::cout << gal.dm_mass(gal.virial_radius(200, true)) << " " << gal.virial_radius(200, true) << std::endl;


    //exit(0);
    //simplifiedmodels::_main();

  
    Cosmology cosmo(CosmoModel::PLANCKONLY);

    //std::cout << cosmo.cosmological_density(0, Species::CDM)*MSOL_to_GeV*pow(kpc_to_cm, -3) << std::endl;
    
    DarkHalo dh(cosmo, 0, 1, 3, 1);
    dh.Initialise_from_virial_parameters(1);
    std::cout << dh.virial_concentration(200, true) << std::endl;
    std::cout << dh.virial_radius(200, true) << std::endl;
    std::cout << dh.get_rhos() << std::endl;
    std::cout << dh.get_rs() << std::endl;
    

    //std::cout << 4*G_NEWTON*m_to_kpc*MSOL_to_kg << std::endl;

    exit(0);
    //FSLModel fsl(1e-10, 1e-2, dh);

    //std::cout << fsl.c_bar(1e-12) << std::endl;
    DegreeFreedom *degFree = new DegreeFreedom(T_QCD_GeV);
    
    MinimalMass MM(degFree, cosmo);
    std::cout << MM.MassesComputation(1e+4, 1e+6)[1]*GeV_to_MSOL << std::endl;

    //Galaxy gal(cosmo, 0, 0);
    //std::cout << gal.dm_density(1)*MSOL_to_GeV*pow(kpc_to_cm, -3) << std::endl;
    //std::cout << dh.dm_density(1)*MSOL_to_GeV*pow(kpc_to_cm, -3) << std::endl;

    exit(0);
    
    
    
    //DirectDetectionConstraints("Scalar", "Coupling_Scalar_Yukawa", "DirectDetectionConstraint_Scalar_Yukawa", CouplingType::UNIVERSAL_YUKAWA, cosmo);


    /* 
    exit(0);
    std::cout << pow(GeV_to_cmm1,-2) << std::endl;
    CPOddHiggsLowMass *cp = new CPOddHiggsLowMass("../input/my_input.dat");
    cp->Write("../input/CPOdd_input.dat");
    //cp->plot_MixingParameters();
    std::vector<DarkSectorModel> models_cp = InputReader::getInstance()->Read("../input/CPOdd_input.dat");
    //std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/Simple_CrossSections_phiphi/ScalarPseudoScalar_universal.in");
    // Kill the reader singleton to free memory and because we do not need it anymore for the moment
    InputReader::kill();

    Particle chi = models_cp[0].get_chi(0);

    // Suppose here that we are at T = 10*T_QCD
    models_cp[0].get_phip_ptr(0)->set_width(cp->GammaTAS(0, 10)); // Attention here the decay width is set a bit differently for A1 
    std::cout << models_cp[0].get_phip_ptr(0)->get_mass() << "\t" << models_cp[0].get_phip_ptr(0)->get_width() << std::endl; 
    std::cout << models_cp[0].get_chi_ptr(0)->get_mass() << std::endl; 
    CrossSection_ccgaga *sigma_ccgaga = new CrossSection_ccgaga(models_cp[0], chi, chi, cp, 0, 10);
    CrossSection_ccff *sigma_ccee = new CrossSection_ccff(models_cp[0], chi, chi, StandardModel::getInstance()->get_couples_SM_ferm(3)); // Coupling to fermions
    //CrossSection_ccff *sigma_ccuu = new CrossSection_ccff(models_cp[0], chi, chi, StandardModel::getInstance()->get_couples_SM_ferm(7));
    //CrossSection_ccff *sigma_ccdd = new CrossSection_ccff(models_cp[0], chi, chi, StandardModel::getInstance()->get_couples_SM_ferm(6));
    double mA1 = models_cp[0].get_phip_ptr(0)->get_mass();
    //sigma_ccgaga->plot(100*mA1*mA1, "gaga_CPOddHiggs_triangle_only");  
    //sigma_ccee->plot(100*mA1*mA1, "ee_CPOddHiggs");
    //sigma_ccuu->plot(4, "uu_CPOddHiggs");
    //sigma_ccdd->plot(4, "dd_CPOddHiggs");

    exit(0);*/

    //ps.plot_ModifiedPowerSpectrumSmallScale(0, 1e-3, 10);
    //ps.plot_ModifiedPowerSpectrumSmallScale(0, 1e-3, 100);
    //ps.plot_ModifiedPowerSpectrumSmallScale(0, 1e-3, 1000);
    ///MinMass.PlotFreeStreamingLength(0.001, 100);
    //MinMass.PlotMinimalMass(1);
    //exit(0);
    //ConstraintCouplingUniversal ccu(CouplingType::UNIVERSAL_YUKAWA);
    //ccu.Coupling_correctAbundance("Scalar", "Coupling_Scalar_Yukawa");
    //ccu.plot_approximate_Yf();
    //ccu.Write_input_MinimalMassHalos_vs_massDM_vs_massProp("PseudoScalar_test", Proptype::pseudoscalar, 20, 0.01, 1e+4);
    //ccu.WriteMassesFromModel("PseudoScalar_test");
    //ccu.Coupling_correctAbundance("Scalar", "Coupling_Scalar_single");
    //ccu.Coupling_correctAbundance("Scalar", "Coupling_Scalar_v2");
    //ccu.TkdAndMinMass("PseudoScalar", "Coupling_PseudoScalar_single_ee", "Mass_PseudoScalar_single_ee");
    //ccu.TkdAndMinMass("Scalar", "Coupling_Scalar_single_ee", "Mass_Scalar_single_ee");

    /*ccu.UnevolvedNsubhalos("Mass_PseudoScalar", "NumberSubhalos_PseudoScalar");
    ccu.SelfInteractingCS("PseudoScalar", "Coupling_PseudoScalar", "CS_SI_PseudoScalar");
    ccu.DecayTime("PseudoScalar", "Coupling_PseudoScalar", "DecayTime_PseudoScalar");
    ccu.UnevolvedNsubhalos("Mass_Scalar", "NumberSubhalos_Scalar");
    ccu.SelfInteractingCS("Scalar", "Coupling_Scalar", "CS_SI_Scalar");
    ccu.DecayTime("Scalar", "Coupling_Scalar", "DecayTime_Scalar");*/

    //-----------------------
    /* ccu.UnevolvedNsubhalos("Mass_PseudoScalar_single_ee", "NumberSubhalos_PseudoScalar_single_ee");
    ccu.SelfInteractingCS("PseudoScalar", "Coupling_PseudoScalar_single_ee", "CS_SI_PseudoScalar_single_ee");
    ccu.DecayTime("PseudoScalar", "Coupling_PseudoScalar_single_ee", "DecayTime_PseudoScalar_single_ee");
    ccu.UnevolvedNsubhalos("Mass_Scalar_single_ee", "NumberSubhalos_Scalar_single_ee");
    ccu.SelfInteractingCS("Scalar", "Coupling_Scalar_single_ee", "CS_SI_Scalar_single_ee");
    ccu.DecayTime("Scalar", "Coupling_Scalar_single_ee", "DecayTime_Scalar_single_ee");*/
    //-----------------------

    //plot_CrossSectionsScatt_universalCouplings();
    //std::vector<DarkSectorModel> model_scalar_axialvector = InputReader::getInstance()->Read("../input/Simple_CrossSections/ScalarAxialVector_universal.in");
    //std::cout << "In the model : " <<  model_scalar_axialvector[0].get_b_VEC_SM(0, 4) << std::endl;
    //std::cout << pow(GeV_to_cmm1, -2) << std::endl;
    //plot_CrossSections_Annihilation_PhiPhi_universalCouplings();
    //plot_CrossSections_Annihilation_FermionFermion_universalCouplings();
    //plot_CrossSections_SelfInteraction();

    exit(0);

    std::vector<DarkSectorModel> test = InputReader::getInstance()->Read("../input/Test_content.in");
    //test[0].get_phip_ptr(0)->set_width(50);
    //std::cout << StandardModel::getInstance()->get_couples_SM_ferm(10).get_particle(0).get_name() << std::endl;
    //std::cout << StandardModel::getInstance()->get_couples_SM_ferm(10).get_particle(0).get_mass() << std::endl;
    //std::cout << test[0].get_phis_ptr(0)->get_width() / test[0].get_phis_ptr(0)->get_mass() << std::endl;
    //std::cout << test[0].get_chi(0).get_mass() << std::endl;

    DecayRates dcc;
    dcc.set_width_phip(test[0]);

    
    KineticDecoupling KinDec(test[0], degFree, 0);
    //KinDec.plotGammaTot();

    std::cout << KinDec.SolveTemperatureEI() << std::endl;
    //KinDec.plotTemperatureEvolutionEI("test");

    exit(0);


    std::clog.rdbuf(prevclogbuf);
    outlogfile.close();

    return 0;
}


std::vector<double> DirectDetectionConstraintsOneModel(DarkSectorModel &model, double lambda, ConstraintCouplingUniversal ccu)
{

    // We set the model with the constrained value of lambda
    bool widthParticleTooLarge = ccu.SetCoupling(model, lambda, lambda, 3);

    DirectDetection DDetection(model);

    int WIMP, Halo;
    int Detector1, Detector2, Detector3;
    int Detector10, Detector20, Detector30;

    Halo = DDCalc::InitHalo();
    WIMP = DDCalc::InitWIMP();

    //DDCalc::SetSHM(Halo, 0.3, 220.0, 220.0);

    Detector1 = DDCalc::XENON1T_2018_Init();
    Detector2 = DDCalc::CDMSLite_Init();
    Detector3 = DDCalc::CRESST_II_Init();

    Detector10 = DDCalc::XENON1T_2018_Init();
    Detector20 = DDCalc::CDMSLite_Init();
    Detector30 = DDCalc::CRESST_II_Init();

    double mDM = model.get_chi(0).get_mass();
    double DM_spin = model.get_chi(0).get_spin();

    DDCalc::SetWIMP_NREffectiveTheory(WIMP, mDM, DM_spin);

    double coeff = 0;
    for (int i = 1; i < 15; i++)
    {
        if (i != 2)
        {
            coeff = DDetection.NROperatorCoefficient_ISO(i, 0);
            DDCalc::SetNRCoefficient(WIMP, i, 0, coeff);

            coeff = DDetection.NROperatorCoefficient_ISO(i, 1);
            DDCalc::SetNRCoefficient(WIMP, i, 1, coeff);
        }
    }

    //std::cout << Detector1 << " " << Detector2 << " " << Detector3 << " " << Detector10 << std::endl;

    // Computing the rates with DDCalc
    DDCalc::CalcRates(Detector1, WIMP, Halo);
    DDCalc::CalcRates(Detector2, WIMP, Halo);
    DDCalc::CalcRates(Detector3, WIMP, Halo);


    // Reference point with no DM signal to evaluate the likelihood with background only
    int WIMP0 = DDCalc::InitWIMP();
    DDCalc::SetWIMP_msigma(WIMP0, mDM, 0., 0., 0., 0.);
    DDCalc::CalcRates(Detector10, WIMP0, Halo);
    DDCalc::CalcRates(Detector20, WIMP0, Halo);
    DDCalc::CalcRates(Detector30, WIMP0, Halo);

    double limit = 1.64;

    double logL_XENON1T = 2 * (DDCalc::LogLikelihood(Detector10) - DDCalc::LogLikelihood(Detector1)) - limit;   // XENON1T
    double logL_CDMSLite = 2 * (DDCalc::LogLikelihood(Detector20) - DDCalc::LogLikelihood(Detector2)) - limit;  // CDMSLite
    double logL_CRESST_II = 2 * (DDCalc::LogLikelihood(Detector30) - DDCalc::LogLikelihood(Detector3)) - limit; // CRESST_II

    std::vector<double> res = {logL_XENON1T, logL_CDMSLite, logL_CRESST_II};

    // Clean up all the objects
    DDCalc::FreeAll();

    return res;
}

void DirectDetectionConstraints(std::string const &name_model_file_input, std::string const &name_coupling_file_input, std::string const &name_output, CouplingType cptype, Cosmology cosmo)
{
    std::ofstream outfile;
    outfile.open("../output/SimplifiedModels/DirectDetectionConstraints/" + name_output + ".out");
    outfile.precision(8);

    ConstraintCouplingUniversal ccu(cptype, cosmo);

    // Construction of the dark sector models - one model for one mass
    std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/MassHalos/" + name_model_file_input + ".in"); // We collect the models
    std::vector<std::vector<double> > lambda = ccu.ReadCouplings("../output/SimplifiedModels/ConstrainedCouplings/" + name_coupling_file_input + ".out");      // we collect the results of the coupling constraints

    int n_phis = models[0].get_n_DS_phis();
    int n_phip = models[0].get_n_DS_phip();
    int n_X = models[0].get_n_DS_X();
    int n_DM = models[0].get_n_DS_DM();

    outfile << "# Model input file is : "
            << "../input/MassHalos/" + name_model_file_input + ".in" << std::endl;
    outfile << "# Coupling input file is : "
            << "../input/MassHalos/" + name_coupling_file_input + ".in" << std::endl;
    outfile << "# Number of DM particles : " << n_DM << std::endl;
    outfile << "# Number of scalar mediators : " << n_phis << std::endl;
    outfile << "# Number of pseudoscalar mediators : " << n_phip << std::endl;
    outfile << "# Number of vector mediators : " << n_X << std::endl;
    outfile << "# Number of treated points : " << models.size() << std::endl;
    outfile << "##############" << std::endl;
    outfile << "# For every DM particles their type is" << std::endl;
    for (int i = 0; i < n_DM; i++)
        outfile << "## DM Particle " << i << " is of type : " << models[0].get_chi(i).get_fermiontype_str() << std::endl;
    outfile << "##############" << std::endl;

    outfile << "# If the LogL - limit is above 0 then the coupling is constrained by the corresponding experiment" << std::endl;
    outfile << "# Model number | LogL - limit (XENON1T | CDMSLite | CRESST_II)" << std::endl;

    std::vector<double> LogL;

    for (int i = 0; i < lambda.size(); i++)
    {
        LogL = DirectDetectionConstraintsOneModel(models[(int)lambda[i][0]], lambda[i][1], ccu);
        outfile << lambda[i][0] << "\t" << LogL[0] << "\t" << LogL[1] << "\t" << LogL[2] << std::endl;
    }

    outfile.close();
}


