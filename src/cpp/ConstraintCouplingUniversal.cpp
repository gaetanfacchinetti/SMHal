#include "../headers/ConstraintCouplingUniversal.h"

double ConstraintCouplingUniversal::Coupling_correctAbundance(std::string const &name_input, std::string const &name_output)
{

    std::ofstream outfile;
    outfile.open("../output/SimplifiedModels/ConstrainedCouplings/" + name_output + ".out");

    // Construction of the dark sector models - one model for one mass
    std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/MassHalos/" + name_input + ".in");
    InputReader::kill();

    int n_phis = models[0].get_n_DS_phis();
    int n_phip = models[0].get_n_DS_phip();
    int n_X = models[0].get_n_DS_X();
    int n_DM = models[0].get_n_DS_DM();

    outfile << "# Model input file is : "
            << "../input/MassHalos/" + name_input + ".in" << std::endl;
    outfile << "# Number of DM particles : " << n_DM << std::endl;
    outfile << "# Number of scalar mediators : " << n_phis << std::endl;
    outfile << "# Number of pseudoscalar mediators : " << n_phip << std::endl;
    outfile << "# Number of vector mediators : " << n_X << std::endl;
    outfile << "# Number of treated points : " << models.size() << std::endl;
    outfile << "##############" << std::endl;
    for (int i = 0; i < n_DM; i++)
        outfile << "## DM Particle " << i << " is of type : " << models[0].get_chi(i).get_fermiontype_str() << std::endl;
    outfile << "##############" << std::endl;

    outfile << "# Model number | coupling | xcd | widthParticleTooLarge | GSLErrorOccured | NoFreezeOut | NeedTooLargeCouplings" << std::endl;

    // Degrees of freedom initialisation (without putting the DS content into it)

    

    // Loop on all the models (use of parallelisation through open mp to make the computation faster)
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < models.size(); i++)
        {
            std::vector<double> lambda;
            std::string result = "";
            printf("Model : %d treated by thread %d \n", i, omp_get_thread_num());
            DegreeFreedom *degFree = new DegreeFreedom(T_QCD_GeV, &models[i]);
            lambda = ConstraintOneModel(models[i], degFree, Cosmo.get_Omega_c_h2());
            result = std::to_string(i) + "\t" + std::to_string(lambda[0]) + "\t" + std::to_string(lambda[1]) + "\t" + std::to_string(lambda[2]) + "\t" + std::to_string(lambda[3]) + "\t" + std::to_string(lambda[4]) + "\t" + std::to_string(lambda[5]);
            outfile << result << std::endl;
            delete degFree;
        }
    }

    return 1.;
}

double ConstraintCouplingUniversal::TkdAndMinMass(std::string const &name_model_file_input, std::string const &name_coupling_file_input, std::string const &name_output)
{

    std::ofstream outfile;
    outfile.open("../output/SimplifiedModels/MinimalMassHalos/" + name_output + ".out");
    outfile.precision(8);

    // Construction of the dark sector models - one model for one mass
    std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/MassHalos/" + name_model_file_input + ".in"); // We collect the models
    std::vector<std::vector<double> > lambda = ReadCouplings("../output/SimplifiedModels/ConstrainedCouplings/" + name_coupling_file_input + ".out");          // we collect the results of the coupling constraints
    InputReader::kill();

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

    outfile << "# Model number | Tkd [GeV] | massao [Msol] | massfs [Msol] | 1 if nothermaleq -> mass_max | 1 ifalwaysthermaleq -> mass_min | xkd (only if < xcd)" << std::endl;

    // Degrees of freedom initialisation (without putting the DS content into it)

    //std::cout << degFree->get_hSmoothInterp()(-2) << std::endl;
    //std::vector<double> lambda;

    //std::cout << "Introduction over" << std::endl;
    std::vector<double> TkdMM;

    //std::cout << lambda.size() << " " << lambda[0][0] << " " << lambda[0][1] << std::endl;

    // Loop on all the models
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < lambda.size(); i++)
        {
            
            printf("Model : %f treated by thread %d \n", lambda[i][1], omp_get_thread_num());
            DegreeFreedom *degFree = new DegreeFreedom(T_QCD_GeV, &models[(int)lambda[i][0]]);

            double kinDecBeforeChemDec = 0;
            TkdMM = TemperatureKdOneModel(models[(int)lambda[i][0]], lambda[i][1], degFree);

            // Compare xkd to xcd (stored in lambda[i][2]).
            if (models[i].get_chi(0).get_mass() / TkdMM[0] < lambda[i][2])
                kinDecBeforeChemDec = models[i].get_chi(0).get_mass() / TkdMM[0];

            outfile << lambda[i][0] << "\t" << TkdMM[0] << "\t" << TkdMM[1] << "\t" << TkdMM[2] << "\t" << TkdMM[3] << "\t" << TkdMM[4] << "\t" << kinDecBeforeChemDec << std::endl;

            delete degFree;
        }
    }
    
    outfile.close();
    return 1;
}

double ConstraintCouplingUniversal::SelfInteractingCS(std::string const &name_model_file_input, std::string const &name_coupling_file_input, std::string const &name_output)
{

    std::ofstream outfile;
    outfile.open("../output/SimplifiedModels/SelfInteractions/" + name_output + ".out");
    outfile.precision(8);

    // Construction of the dark sector models - one model for one mass
    std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/MassHalos/" + name_model_file_input + ".in"); // We collect the models
    std::vector<std::vector<double> > lambda = ReadCouplings("../output/SimplifiedModels/ConstrainedCouplings/" + name_coupling_file_input + ".out");          // we collect the results of the coupling constraints
    InputReader::kill();

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

    outfile << "# Model number |  vrel*sigma_T/mchi [cm^2 g^{-1} km s^{-1}] @ vrel = 1e+2 km/s  | vrel*sigma_T/mchi [cm^2 g^{-1} km s^{-1}] @ vrel = 1e+3 km/s | vrel*sigma_T/mchi [cm^2 g^{-1} km s^{-1}] @ vrel = 1e+4 km/s " << std::endl;

    double conv = pow(GeV_to_cmm1, -2) * pow(GeV_to_kg * 1e+3, -1);
    double c_kms = C_LIGHT * 1e-3;                      // Light speed in km/s
    double vrel_1 = 1e+2, vrel_2 = 1e+3, vrel_3 = 1e+4; // km/s
    double res_1 = 0, res_2 = 0, res_3 = 0, m0 = 0;

    for (int i = 0; i < lambda.size(); i++)
    {
        Particle chi = models[(int)lambda[i][0]].get_chi(0);

        SetCoupling(models[(int)lambda[i][0]], lambda[i][1], lambda[i][1]);

        CrossSectionMajorana_cccc *CS_M = new CrossSectionMajorana_cccc(models[(int)lambda[i][0]], chi);
        m0 = chi.get_mass();
        res_1 = vrel_1 * CS_M->EvalSigmaTransfer_fVrel(vrel_1 / c_kms) / m0 * conv;
        res_2 = vrel_2 * CS_M->EvalSigmaTransfer_fVrel(vrel_2 / c_kms) / m0 * conv;
        res_3 = vrel_3 * CS_M->EvalSigmaTransfer_fVrel(vrel_3 / c_kms) / m0 * conv;

        outfile << (int)lambda[i][0] << "\t" << res_1 << "\t" << res_2 << "\t" << res_3 << std::endl;
    }

    outfile.close();
    return 1;
}

double ConstraintCouplingUniversal::DecayTime(std::string const &name_model_file_input, std::string const &name_coupling_file_input, std::string const &name_output)
{

    std::ofstream outfile;
    outfile.open("../output/SimplfiedModels/DecayTime/" + name_output + ".out");
    outfile.precision(8);

    // Construction of the dark sector models - one model for one mass
    std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/MassHalos/" + name_model_file_input + ".in"); // We collect the models
    std::vector<std::vector<double> > lambda = ReadCouplings("../output/SimplifiedModels/ConstrainedCouplings/" + name_coupling_file_input + ".out");          // we collect the results of the coupling constraints
    InputReader::kill();

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

    outfile << "# Model number | Decay time [s] " << std::endl;

    double max_time = 0;
    double new_time = 0;

    for (int i = 0; i < lambda.size(); i++)
    {
        new_time = 0;
        max_time = 0;

        //Particle chi = models[(int)lambda[i][0]].get_chi(0);
        // We set the model with the constrained value of lambda
        SetCoupling(models[(int)lambda[i][0]], lambda[i][1], lambda[i][1]);

        //std::cout << "------------" << std::endl;

        for (int j = 0; j < n_phis; j++)
        {
            if ((models[(int)lambda[i][0]].get_phis(j).get_width()) > 0)
            {
                new_time = 1. / (models[(int)lambda[i][0]].get_phis(j).get_width());

                if (new_time > max_time)
                    max_time = new_time;
            }

            //std::cout << i << " " << (int)lambda[i][0] << " " << lambda[i][1] << " " << new_time << " " << (models[(int)lambda[i][0]].get_phis(j).get_width()) << std::endl;
            //std::cout << "########" << std::endl;
        }
        for (int j = 0; j < n_phip; j++)
        {
            if ((models[(int)lambda[i][0]].get_phip(j).get_width()) > 0)
            {

                new_time = 1. / (models[(int)lambda[i][0]].get_phip(j).get_width());

                if (new_time > max_time)
                    max_time = new_time;
            }
        }
        for (int j = 0; j < n_X; j++)
        {
            if ((models[(int)lambda[i][0]].get_X(j).get_width()) > 0)
            {
                new_time = 1. / (models[(int)lambda[i][0]].get_X(j).get_width());

                if (new_time > max_time)
                    max_time = new_time;
            }
        }

        outfile << (int)lambda[i][0] << "\t" << max_time / s_to_GeVm1 << std::endl;
    }

    outfile.close();
    return 1;
}

double ConstraintCouplingUniversal::UnevolvedNsubhalos(std::string const &name_minmass_file_input, std::string const &name_output)
{

    std::ofstream outfile;
    outfile.open("../output/SimplifiedModels/UnevolvedNsubhalos/" + name_output + ".out");
    outfile.precision(8);

    // Construction of the dark sector models - one model for one mass
    //std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/MassHalos/" + name_model_file_input + ".in"); // We collect the models
    std::vector<std::vector<double> > MinMass = ReadMinimalMass("../output/SimplifiedModels/ConstrainedCouplings/" + name_minmass_file_input + ".out"); // we collect the results of the coupling constraints
    //InputReader::kill();

    outfile << "# MinMass input file is : "
            << "../input/MassHalos/" + name_minmass_file_input + ".in" << std::endl;
    outfile << "##############" << std::endl;

    outfile << "# Model number | Nsub " << std::endl;

    double mmass = 0;
    DarkHalo dh(Cosmo, 0, 1, 3, 1);
    dh.Initialise_from_virial_parameters(1e+12);

    // Loop on all the models
    for (int i = 0; i < MinMass.size(); i++)
    {
        mmass = std::max(MinMass[i][2], MinMass[i][3]);

        if (mmass > 0)
        {
            FSLModel fsl(mmass, 0, dh);
            outfile << MinMass[i][0] << "\t" << fsl.NsubUnevolved() << std::endl;
        }
        else
            outfile << MinMass[i][0] << "\t" << 0 << std::endl;
    }

    outfile.close();
    return 1;
}

std::vector<double> ConstraintCouplingUniversal::ConstraintOneModel(DarkSectorModel &model, DegreeFreedom *degFree, double _Omega_C_h2)
{

    // ------- Initialisation of all parameters ------- //

    // Dark Matter mass
    double mass_chi = model.get_chi(0).get_mass();

    // Variables/parameters needed in the computation
    double Omega_h2_comp, Omega_h2_comp_min, Omega_h2_comp_max;
    double Y0, xcd;
    std::vector<double> xY{0, 0};
    int comptDichotomieMax = 100, comptDichotomie = 0, comptFirstEval = 0;

    // Coefficient to convert the comoving number density into abundance
    double coeffOmega_h2 = 16 * pow(PI, 3) / 135 * 3.93773 * pow(T0_CMB_GeV, 3) * pow(100 / kpc_to_m / GeV_to_sm1, -2) * pow(M_PLANCK, -2);

    // Resetting the boolean for this model
    // This go to true if an error is detected (ususally couplings too low)
    bool freezeInNotFreezeOut = false;
    bool GSL_error_occured = false;
    bool widthParticleTooLarge = false;

    // Result vector
    std::vector<double> res;
    res.resize(6, 0);

    // Setting the minimal and maxmial values of the couplings
    double lambda_min_0 = 1e-5, lambda_max_0 = sqrt(4 * PI);
    double lambda_min_1 = 1e-5, lambda_max_1 = sqrt(4 * PI);
    double lambda_new = 0;

    // ------- Trying to find the minimal coupling that does a good job ------- //
    
    do
    {

        GSL_error_occured = false;
        freezeInNotFreezeOut = false;

        try
        {
            // Check the minimal value of lambda first
            // We set the model with the constrained value of lambda
            widthParticleTooLarge = SetCoupling(model, lambda_min_1, lambda_min_1, 3);

            // The chemical decoupling object is a new object for a new model - i.e. new DM mass in our case
            ChemicalDecoupling ChemDec(model, degFree, CS_TYPE::FULL_ONLY);

            // Computation of the comoving number density today
            xY = ChemDec.SolveFreezeOut();
            xcd = xY[2]; // x decoupling
            //Y0 = xY[1];  // Assume no later evolution
            Y0 = ChemDec.ComputeAbundanceToday(xY);
            Omega_h2_comp_max = Y0 * coeffOmega_h2 * mass_chi;
            //std::cout << Omega_h2_comp_max << std::endl;
        }
        catch (const ExceptionHandler &e)
        {
            // If we have suspitions of freeze in or if we have GSL_error the it probably means that
            // the initial guess for the coupling was simply to low so we can start with a new one
            // If after several trials it does not improve, then it means that freeze in must occur anyway
        

            //std::cout << e.what() << std::endl;

            if (e.get_exception_type() == ExceptionType::GSL_error || e.get_exception_type() == ExceptionType::FreezeIn)
            {

                lambda_min_1 = lambda_min_0 * pow(10, (comptFirstEval + 1) * log10(lambda_max_0 / lambda_min_0) / 10.);
                //std::cout << "Ici : " << comptFirstEval << " " << lambda_min_1 << std::endl;

                Omega_h2_comp_max = 2 * _Omega_C_h2;

                if (e.get_exception_type() == ExceptionType::GSL_error)
                    GSL_error_occured = true;

                if (e.get_exception_type() == ExceptionType::FreezeIn)
                    freezeInNotFreezeOut = true;
            }
        }
        // If the first evlauation gives something wrong we also retry at larger lambda
        if (Omega_h2_comp_max < _Omega_C_h2)
        {
            if (NUM_THREADS == 1)
                std::cout << "First evaluation of Omega_h2 not satysfying (< Omega_h2_true)" << std::endl;

            lambda_min_1 = lambda_min_0 * pow(10, (comptFirstEval + 1) * log10(lambda_max_0 / lambda_min_0) / 10);
        }

        comptFirstEval++;

    } while ((GSL_error_occured || freezeInNotFreezeOut || Omega_h2_comp_max < _Omega_C_h2) && comptFirstEval < 11);


    std::cout << "We are here " << lambda_min_1 << " " << lambda_max_1 << std::endl;

    // If problem we do not go further
    if (GSL_error_occured)
    {

        res[3] = 1;
        return res;
    }
    if (freezeInNotFreezeOut)
    {
        res[4] = 1;
        return res;
    }

    if (Omega_h2_comp_max < _Omega_C_h2)
    {
        res[3] = 1;
        res[4] = 1;
        return res;
    }
    


    // ------- Trying to find if the maximal coupling is enough to get lower than Omega_c_h2 measured ------- //

    try
    {
        // Check the minimal value of lambda first
        // We set the model with the constrained value of lambda
        widthParticleTooLarge = SetCoupling(model, lambda_max_1, lambda_max_1, 3);

        // The chemical decoupling object is a new object for a new model - i.e. new DM mass in our case
        ChemicalDecoupling ChemDec(model, degFree, CS_TYPE::FULL_ONLY);

        // Computation of the comoving number density today
        xY = ChemDec.SolveFreezeOut();
        //ChemDec.plotSigmaV();
        xcd = xY[2]; // x decoupling
        //Y0 = xY[1];  // Assume no later evolution
        Y0 = ChemDec.ComputeAbundanceToday(xY);
        Omega_h2_comp_min = Y0 * coeffOmega_h2 * mass_chi;
        //std::cout << model.get_phis(0).get_width()/model.get_phis(0).get_mass() << std::endl;
        //std::cout << Omega_h2_comp_min << " " << xcd << " " << xY[1] * coeffOmega_h2 * mass_chi << std::endl;
    }
    catch (const ExceptionHandler &e)
    {
        if (NUM_THREADS == 1)
            std::cout << "ERROR : Problem when putting the maximum value for the couplings " << __PRETTY_FUNCTION__ << std::endl;
    }

    // If the large coupling is not enough to get the correct abundance we return a problem
    if (Omega_h2_comp_min > _Omega_C_h2)
    {
        res[5] = 1;
        return res;
    }

    comptFirstEval++;

    // We reset the new lower and upper bound of the search
    lambda_min_0 = lambda_min_1;
    lambda_max_0 = lambda_max_1;

    //std::cout << "We are second point " << lambda_min_1 << " " << lambda_max_1 << std::endl;

    // ==============================================================================================
    // If we successfully managed to compute Omega_h2 for the lowest value of the coupling we keep on
    // ==============================================================================================

    do // Loop on the couplings
    {
        //std::cout << "Ici au depart" << std::endl;
        // Initialisation of the dichotomie on the couplings
        lambda_new = sqrt(lambda_max_1 * lambda_min_1);

        // We set the model with the constrained value of lambda
        widthParticleTooLarge = SetCoupling(model, lambda_new, lambda_new, 3);

        //comptGSLerror = 0;

        try
        {
            // The chemical decoupling object is a new object for a new model - i.e. new DM mass in our case
            ChemicalDecoupling ChemDec(model, degFree, CS_TYPE::FULL_ONLY);

            //std::cout << "must be called here " << &ChemDec << std::endl;
            // Computation of the comoving number density today
            xY = ChemDec.SolveFreezeOut();
            xcd = xY[2]; // x decoupling
            //xY.erase(xY.begin() + 2, xY.end());
            //Y0 = xY[1]; // Assume no later evolution
            Y0 = ChemDec.ComputeAbundanceToday(xY);
            Omega_h2_comp = Y0 * coeffOmega_h2 * mass_chi;

            if (NUM_THREADS == 1)
                std::cout << " Omega_h2 : " << Omega_h2_comp << " and coupling is : " << lambda_new << std::endl;
            //std::cout << fabs(Omega_h2_comp - _Omega_C_h2) / _Omega_C_h2 << " " << comptDichotomie << " " << comptDichotomieMax  << std::endl;

            if (Omega_h2_comp > _Omega_C_h2)
                lambda_min_1 = lambda_new;

            if (Omega_h2_comp < _Omega_C_h2)
                lambda_max_1 = lambda_new;

            comptDichotomie++;
            //std::cout << "Ici 3" << std::endl;
        }
        catch (const ExceptionHandler &e)
        {
            std::cout << e.what() << std::endl;
        }

        // Loop for the dichotomie
    } while (fabs(Omega_h2_comp - _Omega_C_h2) / _Omega_C_h2 > 1e-2 && comptDichotomie < comptDichotomieMax && lambda_new < lambda_max_0 && lambda_new > 1.1 * lambda_min_0);

    //std::cout << "here : " << lambda_new << std::endl;

    res[0] = lambda_new;
    res[1] = xcd;

    if (widthParticleTooLarge)
        res[2] = 1;

    return res;
}

std::vector<double> ConstraintCouplingUniversal::TemperatureKdOneModel(DarkSectorModel &model, double lambda, DegreeFreedom *degFree)
{

    // We set the model with the constrained value of lambda
    bool widthParticleTooLarge = SetCoupling(model, lambda, lambda, 3);

    double mass_chi = model.get_chi(0).get_mass();

    double Tkd = 0;

    bool noThermalEq = false;
    bool AlwaysThermalEq = false;
    //bool kineticDecouplingBeforeChemicalDecoupling = false;

    std::vector<double> res;
    res.resize(5, 0);

    if (lambda == 0)
        return res;

    double massAO = 0;
    double massKD = 0;
    std::vector<double> minMass = {0, 0};

    try
    {
        // Finding the kinetic decoupling temperature or equivalentely xkd = mchi(0)/Tkd between [x=1, x=10^{5}]

        if (NUM_THREADS == 1)
            std::cout << "Treating : coupling = " << lambda << std::endl;

        KineticDecoupling KinDec(model, degFree, 0);
        Tkd = KinDec.SolveTemperatureEI(); // temperature in GeV
        //KinDec.plotGammaTot();
        //exit(0);

        MinimalMass MM(degFree, Cosmo);
        minMass = MM.MassesComputation(Tkd, mass_chi);

        //out_massHaloVsMassDM << mass_chi << "\t" << mass_phip << "\t" << Tkd << "\t" << mass_chi / xcd << "\t" << mass_chi / Tkd << "\t" << xcd << "\t" << minimalMass * GeV_to_MSOL << "\t" << (MM.MassesComputation())[0] * GeV_to_MSOL << "\t" << (MM.MassesComputation())[1] * GeV_to_MSOL << "\t" << width_phip / mass_phip << "\t" << lambdaFDM_new << "\t" << Omega_h2_comp << std::endl;
    }
    catch (ExceptionHandler &e)
    {
        std::cout << e.what() << std::endl;

        // First possible exception caught : if gamma always below H, no kinetic equilibrium on the range of temperature probed
        if (e.get_exception_type() == ExceptionType::KineticDec_gammaLessThanH)
        {
            double Tmax = mass_chi / 1.;

            noThermalEq = true;

            MinimalMass MM(degFree, Cosmo);
            minMass = MM.MassesComputation(Tmax, mass_chi);
        }
        else if (e.get_exception_type() == ExceptionType::KineticDec_gammaGreatThanH)
        {
            double Tmin = std::max(mass_chi / 1e+9, 1e-5);

            AlwaysThermalEq = true;

            //std::cout << "Coucou " << Tmin << std::endl;

            MinimalMass MM(degFree, Cosmo);
            minMass = MM.MassesComputation(Tmin, mass_chi);
        }
    }

    res[0] = Tkd;
    res[1] = minMass[0] * GeV_to_MSOL;
    res[2] = minMass[1] * GeV_to_MSOL;

    if (noThermalEq)
        res[3] = 1;
    if (AlwaysThermalEq)
        res[4] = 1;

    return res;
}

// Set all the couplings with the universal value lambda
bool ConstraintCouplingUniversal::SetNewUniversalCoupling(DarkSectorModel &model, double lambda_SM, double lambda_DM)
{
    // We evaluate what are the number of different propagators
    int n_phis = model.get_n_DS_phis();
    int n_phip = model.get_n_DS_phip();
    int n_X = model.get_n_DS_X();
    int n_DM = model.get_n_DS_DM();

    model.ForceInitCouplings(); // Function that initiate the size of the coupling vectors in model

    // ==================================
    // New coupling to the standard model
    // ==================================
    for (int i = 0; i < n_phis; i++)
        for (int j = 0; j < n_SM; j++)
            model.set_coupling_DS_S_SM(i, j, lambda_SM);

    for (int i = 0; i < n_phip; i++)
        for (int j = 0; j < n_SM; j++)
            model.set_coupling_DS_PS_SM(i, j, lambda_SM);

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_SM; j++)
            model.set_coupling_DS_a_VEC_SM(i, j, lambda_SM);

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_SM; j++)
            model.set_coupling_DS_b_VEC_SM(i, j, lambda_SM);

    // ==================================
    // New coupling to Dark Matter
    // ==================================
    for (int i = 0; i < n_phis; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_S_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_phip; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_PS_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                if (model.get_chi(i).get_fermiontype() == Fermiontype::dirac && model.get_chi(j).get_fermiontype() == Fermiontype::dirac)
                    model.set_coupling_DS_a_VEC_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_b_VEC_DM(i, j, k, lambda_DM);

    DecayRates dc;

    bool widthTooLarge = false;
    // Resetting the new width of the particles in the model
    try
    {
        dc.set_width_phis(model);
        dc.set_width_phip(model);
        dc.set_width_X(model);
        dc.set_width_chi(model);
    }
    catch (ExceptionHandler &e)
    {
        if (e.get_exception_type() == ExceptionType::WidthTooLarge)
            widthTooLarge = true;
    }
    // Here we do nothing even if we catch an exception
    // The exception is treated after in the code

    return widthTooLarge;
}

// Set all the couplings with the universal value lambda
bool ConstraintCouplingUniversal::SetNewUniversalYukawaCoupling(DarkSectorModel &model, double lambda_SM, double lambda_DM)
{
    // We evaluate what are the number of different propagators
    int n_phis = model.get_n_DS_phis();
    int n_phip = model.get_n_DS_phip();
    int n_X = model.get_n_DS_X();
    int n_DM = model.get_n_DS_DM();

    model.ForceInitCouplings(); // Function that initiate the size of the coupling vectors in model

    double mass = 0;
    double mt = M_QUARK_TOP;
    double v = 246;

    // ==================================
    // New coupling to the standard model
    // ==================================
    for (int i = 0; i < n_phis; i++)
        for (int j = 0; j < n_SM; j++)
        {
            mass = StandardModel::getInstance()->get_couples_SM_ferm(j).get_particle(0).get_mass();
            model.set_coupling_DS_S_SM(i, j, lambda_SM * mass / v);
        }

    for (int i = 0; i < n_phip; i++)
        for (int j = 0; j < n_SM; j++)
        {
            mass = StandardModel::getInstance()->get_couples_SM_ferm(j).get_particle(0).get_mass();
            model.set_coupling_DS_PS_SM(i, j, lambda_SM * mass / v);
        }

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_SM; j++)
        {
            mass = StandardModel::getInstance()->get_couples_SM_ferm(j).get_particle(0).get_mass();
            model.set_coupling_DS_a_VEC_SM(i, j, lambda_SM * mass / v);
        }

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_SM; j++)
        {
            mass = StandardModel::getInstance()->get_couples_SM_ferm(j).get_particle(0).get_mass();
            model.set_coupling_DS_b_VEC_SM(i, j, lambda_SM * mass / v);
        }

    // ==================================
    // New coupling to Dark Matter
    // ==================================
    for (int i = 0; i < n_phis; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_S_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_phip; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_PS_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                if (model.get_chi(i).get_fermiontype() == Fermiontype::dirac && model.get_chi(j).get_fermiontype() == Fermiontype::dirac)
                    model.set_coupling_DS_a_VEC_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_b_VEC_DM(i, j, k, lambda_DM);

    DecayRates dc;

    bool widthTooLarge = false;
    // Resetting the new width of the particles in the model
    try
    {
        dc.set_width_phis(model);
        dc.set_width_phip(model);
        dc.set_width_X(model);
        dc.set_width_chi(model);
    }
    catch (ExceptionHandler &e)
    {
        if (e.get_exception_type() == ExceptionType::WidthTooLarge)
            widthTooLarge = true;
    }
    // Here we do nothing even if we catch an exception
    // The exception is treated after in the code

    return widthTooLarge;
}

// Set all the couplings to species coupl with value lambda
bool ConstraintCouplingUniversal::SetNewSingleCoupling(DarkSectorModel &model, double lambda_SM, double lambda_DM, int coupl)
{
    // We evaluate what are the number of different propagators
    int n_phis = model.get_n_DS_phis();
    int n_phip = model.get_n_DS_phip();
    int n_X = model.get_n_DS_X();
    int n_DM = model.get_n_DS_DM();

    model.ForceInitCouplings(); // Function that initiate the size of the coupling vectors in model

    // ==================================
    // New coupling to the standard model
    // ==================================
    for (int i = 0; i < n_phis; i++)
        model.set_coupling_DS_S_SM(i, coupl, lambda_SM);

    for (int i = 0; i < n_phip; i++)
        model.set_coupling_DS_PS_SM(i, coupl, lambda_SM);

    for (int i = 0; i < n_X; i++)
        model.set_coupling_DS_a_VEC_SM(i, coupl, lambda_SM);

    for (int i = 0; i < n_X; i++)
        model.set_coupling_DS_b_VEC_SM(i, coupl, lambda_SM);

    // ==================================
    // New coupling to Dark Matter
    // ==================================
    for (int i = 0; i < n_phis; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_S_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_phip; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_PS_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                if (model.get_chi(i).get_fermiontype() == Fermiontype::dirac && model.get_chi(j).get_fermiontype() == Fermiontype::dirac)
                    model.set_coupling_DS_a_VEC_DM(i, j, k, lambda_DM);

    for (int i = 0; i < n_X; i++)
        for (int j = 0; j < n_DM; j++)
            for (int k = 0; k < n_DM; k++)
                model.set_coupling_DS_b_VEC_DM(i, j, k, lambda_DM);

    DecayRates dc;

    bool widthTooLarge = false;
    // Resetting the new width of the particles in the model
    try
    {
        dc.set_width_phis(model);
        dc.set_width_phip(model);
        dc.set_width_X(model);
        dc.set_width_chi(model);
    }
    catch (ExceptionHandler &e)
    {
        if (e.get_exception_type() == ExceptionType::WidthTooLarge)
            widthTooLarge = true;
    }
    // Here we do nothing even if we catch an exception
    // The exception is treated after in the code

    return widthTooLarge;
}

bool ConstraintCouplingUniversal::SetCoupling(DarkSectorModel &model, double lambda_SM, double lambda_DM, int coupl)
{
    bool widthParticleTooLarge = true;

    if (cptype == CouplingType::UNIVERSAL_CONSTANT)
        widthParticleTooLarge = SetNewUniversalCoupling(model, lambda_SM, lambda_DM);
    else if (cptype == CouplingType::SINGLE_CONSTANT)
        widthParticleTooLarge = SetNewSingleCoupling(model, lambda_SM, lambda_DM, coupl);
    else if (cptype == CouplingType::UNIVERSAL_YUKAWA)
        widthParticleTooLarge = SetNewUniversalYukawaCoupling(model, lambda_SM, lambda_DM);
    else
        widthParticleTooLarge = SetNewSingleCoupling(model, lambda_SM, lambda_DM, coupl);

    return widthParticleTooLarge;
}

double ConstraintCouplingUniversal::MaximalValueCouplingDM(DarkSectorModel &model, int coupl)
{
    int n_phis = model.get_n_DS_phis();
    int n_phip = model.get_n_DS_phip();
    int n_X = model.get_n_DS_X();
    int n_DM = model.get_n_DS_DM();

    double res;

    double lambda_max = sqrt(4 * PI);

    // If we are ok with a coupling that is lambda_max, no problems
    if (SetCoupling(model, lambda_max, lambda_max, coupl) == false)
        return lambda_max;

    double lambda_DM_max = sqrt(4 * PI);

    for (int i = 0; i < n_phis; i++)
    {
        if (sqrt(4 * PI * model.get_phis(i).get_mass() / model.get_chi(0).get_mass()) < lambda_DM_max)
            lambda_DM_max = sqrt(4 * PI * model.get_phis(i).get_mass() / model.get_chi(0).get_mass());
    }
    for (int i = 0; i < n_phip; i++)
    {
        if (sqrt(4 * PI * model.get_phip(i).get_mass() / model.get_chi(0).get_mass()) < lambda_DM_max)
            lambda_DM_max = sqrt(4 * PI * model.get_phip(i).get_mass() / model.get_chi(0).get_mass());
    }
    for (int i = 0; i < n_X; i++)
    {
        if (sqrt(4 * PI * model.get_X(i).get_mass() / model.get_chi(0).get_mass()) < lambda_DM_max)
            lambda_DM_max = sqrt(4 * PI * model.get_X(i).get_mass() / model.get_chi(0).get_mass());
    }

    double lambda_min = 1e-5;
    int Npts = 40;
    double dlambda = log10(lambda_max / lambda_min) / (1.0 * Npts);
    double lambda_DM = 0, lambda_SM = 0;

    for (int i = 0; i < Npts; i++)
    {
        lambda_SM = lambda_min * pow(10, i * dlambda);

        for (int j = 0; j < Npts; j++)
        {
            lambda_DM = lambda_min * pow(10, j * dlambda);

            if (SetCoupling(model, lambda_SM, lambda_DM, coupl) == true || lambda_DM > lambda_DM_max)
                std::cout << 0 << " ";
            else
                std::cout << 1 << " ";
        }
        std::cout << std::endl;
    }

    return res;
}

void ConstraintCouplingUniversal::Write_input_MinimalMassHalos_vs_massDM_vs_massProp(std::string const &name_input, Proptype type, int Npts, double massMin, double massMax)
{
    std::ofstream outfile;
    outfile.open("../input/MassHalos/" + name_input + ".in");

    double dl_m = log10(massMax / massMin) / (Npts - 1);

    time_t now = time(0);

    // convert now to string form
    char *dt = ctime(&now);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);

    outfile << "####################################################" << std::endl;
    outfile << "#" << std::endl;
    outfile << "# DMEvolution input file" << std::endl;
    outfile << "# Automatically generated" << std::endl;
    outfile << "#" << std::endl;
    outfile << "# Gaetan Facchinetti" << std::endl;
    outfile << "# gaetan.facchinetti@umontpellier.fr" << std::endl;
    outfile << "# Laboratoire Univers et Particules Montpellier" << std::endl;
    outfile << "# Created on : (UTC) " << dt;
    outfile << "#" << std::endl;
    outfile << "####################################################" << std::endl;
    outfile << std::endl;

    outfile << "## Number of points to evaluate ##" << std::endl;
    outfile << "Npoints :: " << Npts * Npts << std::endl;
    outfile << std::endl;
    outfile << std::endl;

    outfile << "## Dark sector particles ##" << std::endl;
    outfile << "DSCONT :: NDM :: 1" << std::endl;
    if (type == Proptype::scalar)
        outfile << "DSCONT :: DSPart :: (scalar, 1)" << std::endl;
    else if (type == Proptype::pseudoscalar)
        outfile << "DSCONT :: DSPart :: (pseudoscalar, 1)" << std::endl;
    else if (type == Proptype::vector)
        outfile << "DSCONT :: DSPart :: (vector, 3)" << std::endl;
    outfile << std::endl;
    outfile << std::endl;

    outfile << "## Masses dark sector [GeV] ##" << std::endl;

    double massDM;
    double massProp;

    int compt = 0;

    for (int i = 0; i < Npts; i++)
    {
        massDM = massMin * pow(10, i * dl_m);

        for (int j = 0; j < Npts; j++)
        {
            massProp = massMin * pow(10, j * dl_m);
            outfile << "MDMDS :: " << compt << " :: " << massDM << ", " << massProp << std::endl;
            compt++;
        }
    }

    outfile << std::endl;
}

std::vector<std::vector<double> > ConstraintCouplingUniversal::ReadCouplings(std::string const &name)
{
    std::ifstream infile;
    infile.open(name); // opens the file

    std::vector<std::vector<double> > lambda;
    std::vector<double> lambda_mod;
    lambda_mod.resize(6);

    double mod, lamb, xcd, w, f, g;

    std::string line;

    if (!infile.is_open())
    {
        std::cout << "FATAL ERROR : file not found" << std::endl;
        exit(0);
    }

    //std::cout << "here " << name << std::endl;
    while (std::getline(infile, line))
    {
        //std::cout << "here " << line << std::endl;
        if (line.size() > 0)
            if (line[0] != '#')
            {
                std::istringstream iss(line);

                if (!(iss >> mod >> lamb >> xcd >> w >> f >> g))
                {
                    break;
                } // error

                lambda_mod[0] = mod;  // Model number
                lambda_mod[1] = lamb; // lambda
                lambda_mod[2] = xcd;
                lambda_mod[3] = w; // If width was too large
                lambda_mod[4] = f; // If freeze-out was impossible
                lambda_mod[5] = g; // If gsl integration function failed

                lambda.push_back(lambda_mod);
            }
    }

    infile.close();

    return lambda;
}

std::vector<std::vector<double> > ConstraintCouplingUniversal::ReadMinimalMass(std::string const &name)
{
    std::ifstream infile;
    infile.open(name); // opens the file

    std::vector<std::vector<double> > MinMass;
    std::vector<double> MinMass_mod;
    MinMass_mod.resize(6);

    double mod, Tkd, mminao, mminfs;

    std::string line;

    if (!infile.is_open())
    {
        std::cout << "FATAL ERROR : file not found" << std::endl;
        exit(0);
    }

    //std::cout << "here " << name << std::endl;
    while (std::getline(infile, line))
    {
        //std::cout << "here " << line << std::endl;
        if (line.size() > 0)
            if (line[0] != '#')
            {
                std::istringstream iss(line);

                if (!(iss >> mod >> Tkd >> mminao >> mminfs))
                {
                    break;
                } // error

                MinMass_mod[0] = mod; // Model number
                MinMass_mod[1] = Tkd; // lambda
                MinMass_mod[2] = mminao;
                MinMass_mod[3] = mminfs; // If width was too large
                //lambda_mod[4] = f; // If freeze-out was impossible
                //lambda_mod[5] = g; // If gsl integration function failed

                MinMass.push_back(MinMass_mod);
            }
    }

    infile.close();

    return MinMass;
}

double ConstraintCouplingUniversal::WriteMassesFromModel(std::string const &name_input)
{

    std::ofstream outfile;
    outfile.open("../output/MassHalos/" + name_input + "_param.out");

    // Construction of the dark sector models - one model for one mass
    std::vector<DarkSectorModel> models = InputReader::getInstance()->Read("../input/MassHalos/" + name_input + ".in");
    InputReader::kill();

    int n_phis = models[0].get_n_DS_phis();
    int n_phip = models[0].get_n_DS_phip();
    int n_X = models[0].get_n_DS_X();
    int n_DM = models[0].get_n_DS_DM();

    outfile << "# Model input file is : "
            << "../input/MassHalos/" + name_input + ".in" << std::endl;
    outfile << "# Number of DM particles : " << n_DM << std::endl;
    outfile << "# Number of scalar mediators : " << n_phis << std::endl;
    outfile << "# Number of pseudoscalar mediators : " << n_phip << std::endl;
    outfile << "# Number of vector mediators : " << n_X << std::endl;
    outfile << "# Number of treated points : " << models.size() << std::endl;
    outfile << "##############" << std::endl;
    outfile << "# For every DM particles their type is" << std::endl;
    for (int i = 0; i < n_DM; i++)
        outfile << "## DM Particle " << i << " is of type : " << models[0].get_chi(0).get_fermiontype_str() << std::endl;
    outfile << "##############" << std::endl;

    outfile << "# Model number | masses DM (in order) | masses mediators scalar-pseudoscalar-vectors(in order)" << std::endl;

    // Loop on all the models
    for (int i = 0; i < models.size(); i++)
    {
        outfile << i;
        for (int j = 0; j < models[i].get_n_DS_DM(); j++)
            outfile << "\t" << models[i].get_chi(j).get_mass();
        for (int j = 0; j < models[i].get_n_DS_phis(); j++)
            outfile << "\t" << models[i].get_phis(j).get_mass();
        for (int j = 0; j < models[i].get_n_DS_phip(); j++)
            outfile << "\t" << models[i].get_phip(j).get_mass();
        for (int j = 0; j < models[i].get_n_DS_X(); j++)
            outfile << "\t" << models[i].get_X(j).get_mass();

        outfile << std::endl;
    }

    outfile.close();

    return 1;
}

double ConstraintCouplingUniversal::plot_approximate_Yf()
{

    std::ofstream outfile;
    outfile.open("../output/Coupling_scaling_with_sigv_and_mass.out");

    // Construction of the dark sector models - one model for one mass
    DarkSectorModel model;
    model.add_darkmatter(100, 0, Fermiontype::majorana);
    DegreeFreedom *degFree = new DegreeFreedom(T_QCD_GeV);

    int Npts = 100;
    double smin = 1e-30, smax = 1e-20, ds = log10(smax / smin) / (1. * Npts);
    double s, xf, Yeq;
    double coeffOmega_h2 = 16 * pow(PI, 3) / 135 * 3.93773 * pow(T0_CMB_GeV, 3) * pow(100 / kpc_to_m / GeV_to_sm1, -2) * pow(M_PLANCK, -2);

    for (int i = 0; i < Npts + 1; i++)
    {
        s = smin * pow(10, i * ds);
        outfile << s;
        for (int j = 0; j < 5; j++)
        {
            model.get_chi_ptr(0)->set_mass(pow(10, j));
            ChemicalDecoupling ChemDec(model, degFree, CS_TYPE::S_WAVE_ONLY, s / GeVm2_to_cm3persec);

            xf = ChemDec.DeterminationXfMethod2();
            Yeq = 45 * 2 * xf * xf * BesselK2(xf) / (4 * pow(PI, 4) * degFree->get_hSmoothInterp(j));
            outfile << "\t" << Yeq * coeffOmega_h2;
        }
        outfile << std::endl;
    }

    outfile.close();
    return 1;
}


