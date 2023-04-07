#include "../headers/CosmologicalBoost.h"

CosmologicalBoost::CosmologicalBoost(std::shared_ptr<MassFunction> mass_function) : _mass_function(mass_function)
{
    // These are defined in case there is a spike in the matter power spectrum
    // It allows to know quantities on the spike is ignored
    _power_spectrum = _mass_function->get_PowerSpectrum();

    _power_spectrum_without_spikes = _power_spectrum;
    _power_spectrum_without_spikes.reset_spikes();

    _rhom0 = _power_spectrum.get_Cosmology().cosmological_density(0, Species::MATTER) * 1e+9; // [Msol * Mpc^{-3}]

    std::shared_ptr<MassFunction> mass_function_without_spikes;

    if (instanceof<MassFunction_PS>(_mass_function.get()))
        mass_function_without_spikes = std::make_shared<MassFunction_PS>(_power_spectrum_without_spikes);
    if (instanceof<MassFunction_ST>(_mass_function.get()))
        mass_function_without_spikes = std::make_shared<MassFunction_ST>(_power_spectrum_without_spikes, _mass_function->get_extra_params());

    _mass_concentration_model = MassConcentrationModel::MassConcentrationModel_from_MassFunction(_mass_function);
    _mass_concentration_model_without_spikes = MassConcentrationModel::MassConcentrationModel_from_MassFunction(mass_function_without_spikes);

    _sommerfeld_enhancement = SommerfeldEnhancement();
}

double CosmologicalBoost::OneHaloBoost(double m, double zf, double z, DensityProfile profile)
// Dimensionless quantity
{
    double c = 0;

    if (profile == DensityProfile::NFW && zf >= 0)
        c = _mass_concentration_model.ConcentrationMaccio(zf, z);
    else if (profile == DensityProfile::MOORE && zf >= 0) // Still not implemented the correct formula in that case
        c = _mass_concentration_model.ConcentrationMaccio(zf, z);
    else // zc < 0 the halo does not exist yet so it cannot contribute
        return 0;

    DarkHalo dh(_power_spectrum.get_Cosmology(), profile);
    dh.Initialise_from_virial_parameters(m, c, 200., true);

    double Omegamz = _power_spectrum.get_Cosmology().cosmic_abundance(z, Species::MATTER);
    double luminosity = dh.luminosity_profile(c);
    double mass = dh.mass_profile(c);
    double eta = dh.rm2_rs();

    return 200. / Omegamz * pow(c, 3) / 3. * luminosity * pow(mass, -2);
}

double CosmologicalBoost::AverageRedshiftOfCollapse_OnSpike(double z)
// Here we assume that the spike is only relevant for z \gg 1
{
    double ms = 0, res = 0;

    if (_power_spectrum.get_window() == PSWindowType::fourier_space_top_hat)
    {
        std::vector<double> As_spike = _power_spectrum.get_As_spike();
        std::vector<double> ks_spike = _power_spectrum.get_ks_spike();

        if (As_spike.size() >= 1)
        {
            ms = _power_spectrum.Mass_vs_LagrangianRadius(1. / ks_spike[0]);
            // double Dz = _power_spectrum.get_Cosmology().growth_factor_D1_Carroll(z);
            double Nz = _mass_function->Der_MassFraction_spike(As_spike[0], ks_spike[0], z);
            double nup = _mass_function->nuM(1.001 * ms, z);
            double num = _mass_function->nuM(0.999 * ms, z);

            // std::cout << " The average is : " << Nz << std::endl;

            std::vector<double> xx = {0, z};
            /*
            try{
                numerator = gsl_sf_gamma_inc(0, num*num/2.) -  gsl_sf_gamma_inc(0, nup*nup/2.);
                res = std::max(-1. + numerator/Dz/Nz/sqrt(2*PI), z);
            } catch (const ExceptionHandler &e)
            {

                ExceptionHandler ee((char *)"WARNING :: GSL exception", __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::GSL_error);
                //throw ee;

                // By default if the evaluation of the numerator fails
                res = z;
            }*/

            res = GaussLegendre_Integral_Static(0, 2, gl200, num, nup, xx, this, CallBack_fToIntOnmu_AverageRedshiftOfCollapse_OnSpike) / Nz;
        }
    }

    if (res != res)
        return 0;

    return res;
}

double CosmologicalBoost::fToIntOnmu_AverageRedshiftOfCollapse_OnSpike(double nu, double z)
{
    double nuD = nu * _power_spectrum.get_Cosmology_ptr()->growth_factor_D1_Carroll(z);
    return sqrt(2 / PI) * exp(-nu * nu / 2.) * _power_spectrum.get_Cosmology_ptr()->Inverse_growth_factor_Caroll(nuD);
}






double CosmologicalBoost::fToIntOnm_for_TotalBoost(double m, double z, double n_prof, double distrib_model)
{
    if (distrib_model == 0)
    {
        // In that case we don't need to perform an integral
        // The distribution of formation redshift is a delta function
        double zf = SimpleModel_FormationRedshift_Delta(m, z);
        return fToIntOnzf_for_TotalBoost(m, zf, z, n_prof, distrib_model);
    }

    if (distrib_model >= 1)
    {
        //std::cout << "Should go " << distrib_model << std::endl;
        std::vector<double> xx = {m, 0, z, n_prof, distrib_model};
        return GaussLegendre_IntegralLn_Static(1, 5, gl200, z, 1e+5, xx, this, CallBack_fToIntOnzf_for_TotalBoost);
    }

    return 0;
}

double CosmologicalBoost::fToIntOnzf_for_TotalBoost(double m, double zf, double z, double n_prof, double distrib_model)
{
    DensityProfile density_profile;

    if (n_prof == 0)
        density_profile = DensityProfile::NFW;
    if (n_prof == 1)
        density_profile = DensityProfile::MOORE;


    // HaloDistribution dimension of Msol^{-1}*Mpc^{-3}, OneHaloBoost no dimension, 
    return HaloDistribution(m, zf, z, distrib_model) * OneHaloBoost(m, zf, z, density_profile) * m / _rhom0;
}

double CosmologicalBoost::HaloDistribution(double m, double zf, double z, double distrib_model) 
{
    double res = 0;

    // Model distributed according to the mass
    if (distrib_model == 0)
        res = _mass_function->NumberDensity(m, z);

    // Models distributed according to the formation redshift (mass is fixed only by the spike)
    if (distrib_model >= 1)
    {
         
        if (_power_spectrum.get_window() == PSWindowType::fourier_space_top_hat)
        {
            std::vector<double> As_spike = _power_spectrum.get_As_spike();
            std::vector<double> ks_spike = _power_spectrum.get_ks_spike();

            if (As_spike.size() >= 1)
            {
                //double ms = _power_spectrum.Mass_vs_LagrangianRadius(1. / ks_spike[0]);
                double nu =  _mass_function->get_delta_c() / sqrt(_power_spectrum.SigmaR2_narrow_spike(0.5 / ks_spike[0], As_spike[0], ks_spike[0], zf));
                double functional_form = 0;

                if (distrib_model == 1)
                    functional_form = nu * sqrt(2. / PI) * exp(-nu * nu / 2.) / (6 * PI * PI);
                if (distrib_model == 2)
                    functional_form = nu * exp(-nu * nu / 2.) * BBKSFunction(nu) * pow( 2*PI, -2) * pow(3., -3./2.);

                //std::cout << "here" << std::endl;
                //res = pow(ks_spike[0], 3) * functional_form / (1+zf);
                res = pow(ks_spike[0], 3) * functional_form * fabs(_power_spectrum.get_Cosmology_ptr()->derLn_growth_factor_D1_Carroll(zf));
            }
        }
    }

    return res;
}

double CosmologicalBoost::BBKSFunction(double x)
{
    return (pow(x, 3) - 3*x)*( erf( sqrt(5./2.)*x )  + erf( sqrt(5./2.)*x/2. ) )/2. 
        + sqrt(2./5./PI)*( (31.*pow(x, 2)/4. + 8./5.)*exp(-5.*x*x/8.) + (pow(x, 2)/2. - 8./5.)*exp(-5.*x*x/2.) );
}

double CosmologicalBoost::SimpleModel_FormationRedshift_Delta(double m, double z)
{
    double zf = 0, ms = 0;

    if (_power_spectrum.get_window() == PSWindowType::fourier_space_top_hat)
    {
        std::vector<double> ks_spike = _power_spectrum.get_ks_spike();

        // We only treat the case of a single spike here
        // When we have a spike we need to be carefull on how we define zc for ms < m < 100*ms in the Maccio model
        if (ks_spike.size() >= 1)
        {
            ms = _power_spectrum.Mass_vs_LagrangianRadius(1. / ks_spike[0]);

            if (m > ms && m < 100 * ms)                                      // in the case we are at the limit of the transition
                zf = _mass_concentration_model_without_spikes.RedshiftOfCollapse_vs_Mass(m / 100); // Formation redshift
            else if (m != ms)
                zf = _mass_concentration_model.RedshiftOfCollapse_vs_Mass(m / 100); // Formation redshift
            else if (m == ms)
            {

                //zf = AverageRedshiftOfCollapse_OnSpike(z); // Averaged formation redshift
                zf = _mass_concentration_model.RedshiftOfCollapse_vs_Mass(0.999*ms);
                //zf = _mass_concentration_model.RedshiftOfCollapse_vs_Mass(1.001*ms);
                // std::cout << zf << " " << _mass_concentration_model.RedshiftOfCollapse_vs_Mass(1.001*ms) << " " << _mass_concentration_model.RedshiftOfCollapse_vs_Mass(0.999*ms) << std::endl;
            }
        }
        else
            zf = _mass_concentration_model.RedshiftOfCollapse_vs_Mass(m / 100); // Formation redshift
    }
    else
        zf = _mass_concentration_model.RedshiftOfCollapse_vs_Mass(m / 100); // Formation redshift

    return zf; //std::max(zf, z);
}

double CosmologicalBoost::TotalBoost(double z, double mmin, ProfileConfig prof_config, double distrib_model)
// Dimensionless quantity
// The UCMHs are only accounted for with fourier_space_top_hat AND if prof_config = ProfileConfig::MIXED or ProfileConfig::MOORE
// In this function we perform the integral on the mass (Note that is can be an integral over a delta function)
// If distrib_model > 0 and we have no spike it will return 0. Indeed distrib_model > 0 only consider the population in the spike
{
    std::vector<double> xx_nfw = {0, z, 0, distrib_model};
    std::vector<double> xx_moore = {0, z, 1, distrib_model};

    double res = 0;
    double ms = 0;

    if (_power_spectrum.get_window() == PSWindowType::fourier_space_top_hat)
    {
        std::vector<double> As_spike = _power_spectrum.get_As_spike();
        std::vector<double> ks_spike = _power_spectrum.get_ks_spike();

        // Here we only consider at most one spike in this computation
        if (As_spike.size() >= 1)
        {
            ms = _power_spectrum.Mass_vs_LagrangianRadius(1. / ks_spike[0]);

            // Here we perform the integral on the mass
            if (distrib_model == 0)
            {
                // Start by adding the contribution of the spike
                if (ms > mmin)
                {
                    double zf = SimpleModel_FormationRedshift_Delta(ms, z);
                    res += _mass_function->NumberDensity_spike(As_spike[0], ks_spike[0], z) * OneHaloBoost(ms, zf, z, DensityProfile::NFW) * ms / _rhom0;
                }

                // Adapt the configuration to the scenario
                if (prof_config == ProfileConfig::MIXED)
                {
                    res += GaussLegendre_IntegralLn_Static(0, 4, gl100, ms, 1e+16, xx_nfw, this, CallBack_fToIntOnm_for_TotalBoost);
                    res += GaussLegendre_IntegralLn_Static(0, 4, gl100, mmin, ms, xx_moore, this, CallBack_fToIntOnm_for_TotalBoost);
                }
                else if (prof_config == ProfileConfig::ALL_NFW)
                    res += GaussLegendre_IntegralLn_Static(0, 4, gl200, mmin, 1e+16, xx_nfw, this, CallBack_fToIntOnm_for_TotalBoost);
                else if (prof_config == ProfileConfig::ALL_MOORE)
                    res += GaussLegendre_IntegralLn_Static(0, 4, gl200, mmin, 1e+16, xx_moore, this, CallBack_fToIntOnm_for_TotalBoost);
            }
            // Here we perform the integral on the formation redshift
            else if (distrib_model >= 1)
            {
                //Here we have a spike on the mass, so only consider the integrand
                if (prof_config == ProfileConfig::ALL_NFW)
                    res = fToIntOnm_for_TotalBoost(ms, z, 0, distrib_model);
                if (prof_config == ProfileConfig::ALL_MOORE)
                    res = fToIntOnm_for_TotalBoost(ms, z, 1, distrib_model);
            }
        }
        else if (distrib_model == 0)
            res += GaussLegendre_IntegralLn_Static(0, 4, gl200, mmin, 1e+16, xx_nfw, this, CallBack_fToIntOnm_for_TotalBoost);
    }
    else if (distrib_model == 0)
        res += GaussLegendre_IntegralLn_Static(0, 4, gl200, mmin, 1e+16, xx_nfw, this, CallBack_fToIntOnm_for_TotalBoost);

    return 1 + res - pow(_mass_function->MassFractionInHalos(z, mmin), 2);
}

double CosmologicalBoost::MinimalBound_TotalBoost(double z, double mmin)
{
    double Omegamz = _power_spectrum.get_Cosmology().cosmic_abundance(z, Species::MATTER);
    double nu_sqrt2 = _mass_function->get_delta_c() / sqrt(2) / _power_spectrum.SigmaM(mmin, z);
    return 1 + 200. / Omegamz * (1 - erf(nu_sqrt2));
}

void CosmologicalBoost::plot_TotalBoost(double mmin, double distrib_model, std::string add_file_name)
{
    std::ofstream outfile;

    if (add_file_name != "")
        add_file_name = "_" + add_file_name; 

    outfile.open("../output/UCMHS/TotalBoost" + add_file_name + ".out");
    outfile.precision(8);

    //double z, zmin = 1e-2, zmax = 1e+5;
    double rs, rsmin = 1, rsmax = 200;

    int Npts = 100;

    //double d_logz = log10(zmax / zmin) / Npts;
    double drs = log10(rsmax/rsmin)/Npts;

    //outfile << "# mmin = " << mmin << " Msol" << std::endl;
    outfile << "# 1+z | 1 + Boost" << std::endl;

    double z0;

    for (int i = 0; i < Npts + 1; i++)
    {
        //z = zmin * pow(10, i * d_logz);
        rs = rsmin * pow(10, i*drs);
        outfile << rs << "\t" << TotalBoost(rs-1, mmin, ProfileConfig::ALL_NFW, distrib_model) << std::endl; //<< "\t" << MinimalBound_TotalBoost(z, mmin) << std::endl;
    }

    outfile.close();
}



void CosmologicalBoost::plot_FormationRedshiftDistribution(std::string name, double distrib_model)
// For the models where the formation redshift distribution is not a simple dirac
{

    std::ofstream outfile;
    outfile.open("../output/UCMHS/FormationRedshiftDistribution/FormationRedshiftDistribution_" + name + ".out");
    outfile.precision(8);

    double zf, zmin = 1e-2, zmax = 1e+5;

    int Npts = 200;

    double d_logz = log10(zmax / zmin) / Npts;

    outfile << "# zf | ks^{-3} dn/daf" << std::endl;
    
    double ms = 0;
    
    if (_power_spectrum.get_window() == PSWindowType::fourier_space_top_hat)
    {
        std::vector<double> As_spike = _power_spectrum.get_As_spike();
        std::vector<double> ks_spike = _power_spectrum.get_ks_spike();

        // Here we only consider at most one spike in this computation
        if (As_spike.size() >= 1)
        {
            ms = _power_spectrum.Mass_vs_LagrangianRadius(1. / ks_spike[0]);

            for (int i = 0; i < Npts + 1; i++)
            {
                zf = zmin * pow(10, i * d_logz);
                outfile << zf << "\t" << pow(ks_spike[0], -3) * HaloDistribution(ms, zf, 0, distrib_model)*pow(1 + zf, 2) << std::endl;
            }
        }
    }

    outfile.close();
}



double CosmologicalBoost::MassFractionInSpike(double z, double mmin)
{

    if (_power_spectrum.get_window() == PSWindowType::fourier_space_top_hat)
    {
        std::vector<double> As_spike = _power_spectrum.get_As_spike();
        std::vector<double> ks_spike = _power_spectrum.get_ks_spike();

        //std::cout << As_spike.size() << std::endl;

        // Here we only consider at most one spike in this computation
        if (As_spike.size() >= 1)
        {
            
            return _mass_function->Der_MassFraction_spike(As_spike[0], ks_spike[0], z);
        }
    }

    return 0;
}


/// Test function of Boost with Sommerfeld enhancement

double CosmologicalBoost::fToIntOnm_for_TotalBoost_with_Sommerfeld_approx(double m, double z, double ephi, double l, double alpha)
{

    DarkHalo dh(_mass_function->get_cosmo(), 0, 1, 3, 1);
    double velocity = 0;

    if (_power_spectrum.get_As_spike().size() >= 1)
        return 0;

    // In that case we don't need to perform an integral
    // The distribution of formation redshift is a delta function
    double zf = _mass_concentration_model.RedshiftOfCollapse_vs_Mass(m / 100);
    double c = _mass_concentration_model.ConcentrationMaccio(zf, z);

    dh.Initialise_from_virial_parameters(m, c);
    velocity = dh.circular_velocity_DM_only(dh.get_rs())/(1e-3*C_LIGHT); // in units of the speed of light

    //std::cout << velocity << " " << c << std::endl;

    return _mass_function->NumberDensity(m, z) * OneHaloBoost(m, zf, z, DensityProfile::NFW) * m / _rhom0 * pow(2*velocity, 2*l) * _sommerfeld_enhancement.Sommerfeld_Hulthen_v(velocity, ephi, l, alpha);

    return 0;
}



double CosmologicalBoost::TotalBoost_with_Sommerfeld_approx(double z, double ephi, double l, double alpha, double mmin)
// Dimensionless quantity
// Only implemented in the case of a NFW halos here
{
    std::vector<double> xx_nfw = {0, z, ephi, l, alpha};
    double res = GaussLegendre_IntegralLn_Static(0, 5, gl200, mmin, 1e+16, xx_nfw, this, CallBack_fToIntOnm_for_TotalBoost_with_Sommerfeld_approx);
    
    if(l==0)
        return 1 + res - pow(_mass_function->MassFractionInHalos(z, mmin), 2);
    if(l==1)
        return res;
    
    return 0;
}


void CosmologicalBoost::plot_TotalBoost_with_Sommerfeld_approx(double l, double ephi, double alpha, double mmin, std::string add_file_name)
{

    std::ofstream outfile;
    std::vector<double> ephi_arr = _sommerfeld_enhancement.get_ephi_arr();


    if (add_file_name != "")
        add_file_name = "_" + add_file_name; 

    outfile.open("../output/CosmologicalBoost/TotalBoost_with_Sommerfeld_approx_" + add_file_name + ".out");
    outfile.precision(8);

    //double z, zmin = 1e-2, zmax = 1e+5;
    double rs, rsmin = 1, rsmax = 100;

    int Npts = 70;

    //double d_logz = log10(zmax / zmin) / Npts;
    double drs = log10(rsmax/rsmin)/Npts;

    //outfile << "# mmin = " << mmin << " Msol" << std::endl;
    outfile << "# 1+z | 1 + Boost" << std::endl;

    //double ephi;
    
    //ephi = ephi_arr[j];

    for (int i = 0; i < Npts + 1; i++)
    {
        rs = rsmin * pow(10, i*drs);
        outfile << rs << "\t" << TotalBoost_with_Sommerfeld_approx(rs-1, ephi, l, alpha, mmin) << "\t" << TotalBoost_with_Sommerfeld_approx(rs-1, 1e+5, l, alpha, mmin) << std::endl; //<< "\t" << MinimalBound_TotalBoost(z, mmin) << std::endl;
    }



    outfile.close();

}