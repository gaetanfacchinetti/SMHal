#include "../headers/MassFunction.h"


double MassFunction::NumberHalos_unitsOfVolume(double Mmin, double Mmax, double z)
// [in Mpc^{-3}]
{
    std::vector<double> Mz = {0, z};
    return GaussLegendre_IntegralLn_Static(0, 2, gl1000, Mmin, Mmax, Mz, this, CallBack_NumberDensity);
}

double MassFunction::MassInHalos(double z, double Mmin)
{
    std::vector<double> Mz = {0, z};
    return GaussLegendre_IntegralLn_Static(0, 2, gl200, Mmin, 1e+18, Mz, this, CallBack_NumberDensity_Times_M);
}

// Derivative of the mass fraction if Press Schechter
double MassFunction_PS::Der_MassFraction(double M, double z)
// [in Msol^{-1}]
{

    double nu = _delta_c / _power_spectrum.SigmaM(M, z);
    return abs(_power_spectrum.Der_LnSigmaM_vs_M(M, 0)) * nu * sqrt(2. / PI) * exp(-nu * nu / 2);
}

double MassFunction_PS::NumberDensity_spike(double As_spike, double ks_spike, double z)
// [in Msol^{-1} Mpc^{-3}]
{  
    double Ms = _power_spectrum.Mass_vs_LagrangianRadius(1./ks_spike);
    double _rho_bar = 1e+9 * _cosmo.cosmological_density(0, Species::MATTER); // [Msol*Mpc^-3]
    return Der_MassFraction_spike(As_spike, ks_spike, z) * _rho_bar / Ms;
}


double MassFunction_PS::Der_MassFraction_spike(double As_spike, double ks_spike, double z)
{
    double Ms = _power_spectrum.Mass_vs_LagrangianRadius(1./ks_spike);
    double sigmaMs_without_spike = _power_spectrum.SigmaM(1.00001*Ms, z); // Here we take the limit in M->Ms M > Ms of sigma
    double sigma2_spike = _power_spectrum.SigmaR2_narrow_spike_over_window(As_spike, ks_spike, z);

    // This is only valid for the fourier space top hat window function
    if(_window_type == PSWindowType::fourier_space_top_hat)
        return (my_erf_diff(_delta_c/sqrt(2)/sigmaMs_without_spike, _delta_c/sqrt(2)/sqrt(pow(sigmaMs_without_spike, 2) + sigma2_spike)));
    
    return 0;
}

double MassFunction_ST::Der_MassFraction(double M, double z)
{
    double nu = _delta_c / _power_spectrum.SigmaM(M, z);
    double small_a = 0.707;
    double nup = sqrt(small_a) * nu;
    double q = 0.3;
    double A = 0.3222;

    return abs(_power_spectrum.Der_LnSigmaM_vs_M(M, 0)) * nup * A * (1 + pow(nup, -2. * q)) * sqrt(2. * small_a / PI) * exp(-nup * nup / 2.);
}





//
//
// Fonction plot
//
void MassFunction::plot_NormalisedHaloMassFunction(double z)
{
    std::ofstream outfile;
    std::string name_file;
    std::string name_add = "";

    if (instanceof<MassFunction_PS>(this))
        name_add = "PS";
    if (instanceof<MassFunction_ST>(this))
        name_add = "ST";

    if (_window_type == PSWindowType::real_space_top_hat)
        name_file = "../output/CosmologicalHaloMassFunctions/Normalised_HaloMassFunction_" + name_add + "_rsth";
    if (_window_type == PSWindowType::gaussian)
        name_file = "../output/CosmologicalHaloMassFunctions/Normalised_HaloMassFunction_" + name_add + "_gaussian";
    if (_window_type == PSWindowType::fourier_space_top_hat)
        name_file = "../output/CosmologicalHaloMassFunctions/Normalised_HaloMassFunction_" + name_add + "_fsth";

    int int_z = int(z);
    name_file += ("_z" + std::to_string(int_z));
    name_file += ".out";

    outfile.open(name_file);

    outfile << "# Redshift z = " << z << std::endl;
    outfile << "# M [Msol] | HaloMassFunction [Msol^{-1} Mpc^{-3}] | Approx with mass index 2.0 [Msol^{-1} Mpc^{-3}]  | Approx with mass index 1.9 [[Msol^{-1} Mpc^{-3}]  | Approx with mass index 1.95 [[Msol^{-1} Mpc^{-3}]" << std::endl;
    int Npts = 200;
    double Mmin = 1e-12, Mmax = 1e+16, dM = log10(Mmax / Mmin) / Npts, M;
    double norm1, norm2, norm3;
    norm1 = pow(Mmin, -1);
    norm2 = (1. / 0.9) * pow(Mmin, -0.9);
    norm3 = (1. / 0.95) * pow(Mmin, -0.95);

    double norm = NumberHalos_unitsOfVolume(Mmin, Mmax, z);

    for (int i = 0; i < Npts + 1; i++)
    {
        M = Mmin * pow(10, i * dM);
        outfile << M << "\t" << NumberDensity(M, z) / norm << "\t" << pow(M, -2) / norm1 << "\t" << pow(M, -1.9) / norm2 << "\t" << pow(M, -1.95) / norm3 << std::endl;
    }
    outfile.close();
}


/*
void MassFunction::plot_NormalisedHaloMassFunction_ST(double z)
{
    std::ofstream outfile;
    std::string name_file;

    if (_window_type == PSWindowType::real_space_top_hat)
        name_file = "../output/CosmologicalHaloMassFunctions/Normalised_HaloMassFunction_ST_rsth";
    if (_window_type == PSWindowType::gaussian)
        name_file = "../output/CosmologicalHaloMassFunctions/Normalised_HaloMassFunction_ST_gaussian";
    if (_window_type == PSWindowType::fourier_space_top_hat)
        name_file = "../output/CosmologicalHaloMassFunctions/Normalised_HaloMassFunction_ST_fsth";

    int int_z = int(z);
    name_file += ("_z" + std::to_string(int_z));
    name_file += ".out";

    outfile.open(name_file);

    outfile << "# Redshift z = " << z << std::endl;
    outfile << "# M [Msol] | HaloMassFunction_ST [Msol^{-1} Mpc^{-3}] | Approx with mass index 2.0 [Msol^{-1} Mpc^{-3}]  | Approx with mass index 1.9 [[Msol^{-1} Mpc^{-3}]" << std::endl;
    int Npts = 200;
    double Mmin = 1e-12, Mmax = 1e+16, dM = log10(Mmax / Mmin) / Npts, M;
    double norm1, norm2;
    norm1 = pow(Mmin, -1);
    norm2 = (1. / 0.9) * pow(Mmin, -0.9);

    double norm = NumberHalos_unitsOfVolume_ST(Mmin, Mmax, z);

    for (int i = 0; i < Npts + 1; i++)
    {
        M = Mmin * pow(10, i * dM);
        outfile << M << "\t" << NumberDensity_ST(M, z) / norm << "\t" << pow(M, -2) / norm1 << "\t" << pow(M, -1.9) / norm2 << std::endl;
    }
    outfile.close();
}*/

void MassFunction::plot(double z, bool in_units_of_h, std::string add)
{
    std::ofstream outfile;
    std::string name_file = "";

    if (instanceof<MassFunction_PS>(this))
        name_file = "../output/HaloMassFunctions/HaloMassFunction_PS_" +  add;
    if (instanceof<MassFunction_ST>(this))
        name_file = "../output/HaloMassFunctions/HaloMassFunction_ST_" +  add;

    if(add != "")    
        name_file += "_";

    if (_window_type == PSWindowType::real_space_top_hat)
        name_file += "rsth";
    if (_window_type == PSWindowType::gaussian)
        name_file += "gaussian";
    if (_window_type == PSWindowType::fourier_space_top_hat)
        name_file += "fsth";

    int int_z = int(z);
    name_file += ("_z" + std::to_string(int_z));
    name_file += ".out";

    outfile.open(name_file);

    double h_param = 1.;
    if (in_units_of_h == true)
        h_param = _cosmo.Hubble_parameter(0) / 100.;

    outfile << "# Redshift z = " << z << std::endl;

    if (in_units_of_h == false)
        outfile << "# M [Msol] | HaloMassFunction [Msol^{-1} Mpc^{-3}]" << std::endl;
    else
        outfile << "# M [Msol/h] | HaloMassFunction [(Msol/h)^{-1} (Mpc/h)^{-3}]" << std::endl;

    int Npts = 1000;
    double xmin = 1e-12 * h_param, xmax = 1e+16 * h_param, dx = log10(xmax / xmin) / Npts, x;

    for (int i = 0; i < Npts + 1; i++)
    {
        x = xmin * pow(10, i * dx);
        outfile << x << "\t" << NumberDensity(x / h_param, z) * pow(h_param, -4) << std::endl;
    }
    outfile.close();
}



void MassFunction::plot_HaloMassFunction_vs_redshift(double M, std::string add)
{
    std::ofstream outfile;
    std::string name_file = "";

    if (instanceof<MassFunction_PS>(this))
        name_file = "../output/UCMHs/HaloMassFunction_PS_vs_redshift_" +  add;
    if (instanceof<MassFunction_ST>(this))
        name_file = "../output/UCMHs/HaloMassFunction_ST_vs_redshift_" +  add;

    if (_window_type == PSWindowType::real_space_top_hat)
        name_file += "_rsth";
    if (_window_type == PSWindowType::gaussian)
        name_file += "_gaussian";
    if (_window_type == PSWindowType::fourier_space_top_hat)
        name_file += "_fsth"; 

    outfile.open(name_file);

    int Npts = 1000;
    double zmin = 1e-2, zmax = 10, dz = log10(zmax / zmin) / Npts, z;

    for (int i = 0; i < Npts + 1; i++)
    {
        z = zmin * pow(10, i * dz);
        outfile << z << "\t" << NumberDensity(M, z) << std::endl;
    }
    outfile.close();
}


/*
void MassFunction::plot_HaloMassFunction_ST(double z, bool in_units_of_h)
{
    std::ofstream outfile;
    std::string name_file;

    if (_window_type == PSWindowType::real_space_top_hat)
        name_file = "../output/CosmologicalHaloMassFunctions/HaloMassFunction_ST_rsth";
    if (_window_type == PSWindowType::gaussian)
        name_file = "../output/CosmologicalHaloMassFunctions/HaloMassFunction_ST_gaussian";
    if (_window_type == PSWindowType::fourier_space_top_hat)
        name_file = "../output/CosmologicalHaloMassFunctions/HaloMassFunction_ST_fsth";

    int int_z = int(z);

    name_file += ("_z" + std::to_string(int_z));
    name_file += ".out";

    outfile.open(name_file);

    double h_param = 1.;
    if (in_units_of_h == true)
        h_param = Cosmo.Hubble_parameter(0) / 100.;

    outfile << "# Redshift z = " << z << std::endl;

    if (in_units_of_h == false)
        outfile << "# M [Msol] | HaloMassFunction [Msol^{-1} Mpc^{-3}]" << std::endl;
    else
        outfile << "# M [Msol/h] | HaloMassFunction [(Msol/h)^{-1} (Mpc/h)^{-3}]" << std::endl;

    int Npts = 200;
    double xmin = 1e-12 * h_param, xmax = 1e+16 * h_param, dx = log10(xmax / xmin) / Npts, x;

    for (int i = 0; i < Npts + 1; i++)
    {
        x = xmin * pow(10, i * dx);
        outfile << x << "\t" << NumberDensity_ST(x / h_param, z) * pow(h_param, -4) << std::endl;
    }
    outfile.close();
}*/


double MassFunction::NumberProgenitors_EST(double M1, double M2, double z1, double z2)
// This function relies on the EST formalism a k-space top hat is necessary
{
    if (_window_type != PSWindowType::fourier_space_top_hat)
    {
        std::cout << "FATAL ERROR : This function relies on the EST formalism a k-space top hat is necessary -> " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    double S2 = _power_spectrum.SigmaM2(M2, 0);
    double S1 = _power_spectrum.SigmaM2(M1, 0);

    double omega_1 = _delta_c * _cosmo.growth_factor_D1_Carroll(0) / _cosmo.growth_factor_D1_Carroll(z1);
    double omega_2 = _delta_c * _cosmo.growth_factor_D1_Carroll(0) / _cosmo.growth_factor_D1_Carroll(z2);

    double der_S2_M2 = fabs(_power_spectrum.Der_sigmaM2_vs_M(M2, 0));

    double dom = omega_2 - omega_1;
    double ds = S2 - S1;

    return M1 / M2 * dom / sqrt(2 * PI * ds) * exp(-dom * dom / (2 * ds)) * der_S2_M2 / ds;
}

void MassFunction::plot_NumberProgenitors_EST(double M1, double z1, double z2)
// Evaluation of the number of progenitors in the Excursion Set Theory at z2 for a given z1
{
    std::ofstream outfile;
    outfile.open("../output/Masses_Merger_Tree/Masses_at_given_step/NumberProgenitors_EST_z2.out");

    outfile << "# Mass of host : " << M1 << " Msol" << std::endl;
    outfile << "# Redshift of host : " << z1 << std::endl;
    outfile << "# Redshift of evaluation " << z2 << std::endl;
    outfile << "#" << std::endl;
    outfile << "# M2/M1 | Number progenitors" << std::endl;

    int Npts = 500;
    double M2_min = M1 / 500, M2_max = M1, dM2 = log10(M2_max / M2_min) / Npts, M2;

    for (int i = 0; i < Npts + 1; i++)
    {
        M2 = M2_min * pow(10, i * dM2);
        outfile << M2 / M1 << "\t" << NumberProgenitors_EST(M1, M2, z1, z2) << std::endl;
    }

    outfile.close();
}

double MassFunction::fToIntOns_CMF_Density_EST(double S2, double S1, double z)
{
    double M2 = _power_spectrum.Mass_vs_Sigma2(S2, z);
    return 1. / (M2 * sqrt(S2 - S1));
}

double MassFunction::CMF_Density_EST(double M2, double M1, double z)
{
    double S1 = _power_spectrum.SigmaM2(M1, z);
    double S2 = _power_spectrum.SigmaM2(M2, z);
    std::vector<double> xx = {0, S1, z};
    double integral = GaussLegendre_IntegralLn_Static(0, 3, gl500, S1, S2, xx, this, CallBack_fToIntOns_CMF_Density_EST);
    return M1 / sqrt(2 * PI) * integral;
}

double MassFunction::CMF_Approx_EST(double M2, double M1, double z)
{
    double volume_form = _power_spectrum.VolumeForm();
    double Radius = _power_spectrum.LagrangianRadius_vs_Mass(M1);

    return volume_form * pow(Radius, 3) * CMF_Density_EST(M2, M1, z);
}

void MassFunction::plot_CMF_EST(double M1, double z1)
{
    std::ofstream outfile;
    outfile.open("../output/CMF_EST.out");
    outfile << "# Mass of the host halo : " << M1 << " Msol" << std::endl;

    int Npts = 100;
    double M2_min = 1e-12, M2_max = M1, dM2 = log10(M2_max / M2_min) / Npts, M2;

    for (int i = 0; i < Npts + 1; i++)
    {
        M2 = M2_min * pow(10, i * dM2);
        outfile << M2 << "\t" << CMF_Approx_EST(M2, M1, z1) << std::endl;
    }

    outfile.close();
}

double MassFunction::MF_EST(double M2, double M1, double z1)
{
    // Be carefull that the interpolation only goes under 1e-12
    double S1 = _power_spectrum.Interp_SigmaM2(M1, 0);
    double S2 = _power_spectrum.Interp_SigmaM2(M2, 0);

    return M1 / M2 / sqrt(2. * PI) * fabs(_power_spectrum.Der_sigmaM2_vs_M(M2, 0)) / sqrt(S2 - S1);
}

void MassFunction::plot_MF_EST(double M1, double z1)
{
    std::ofstream outfile;
    outfile.open("../output/MF_EST.out");
    outfile << "# Mass of the host halo : " << M1 << " Msol" << std::endl;
    outfile << "# M2 [Msol] | MF (not normalised) [Msol^(-1)]" << std::endl;

    int Npts = 1000;
    double M2_min = 1e-12 * M1, M2_max = M1, dM2 = log10(M2_max / M2_min) / Npts, M2;

    for (int i = 0; i < Npts; i++)
    {
        M2 = M2_min * pow(10, i * dM2);
        outfile << M2 / M1 << "\t" << MF_EST(M2, M1, z1) << std::endl;
    }

    outfile.close();
}

void MassFunction::plot_DeltaS_EST(double M1, int version)
{
    std::ofstream outfile;
    outfile.open("../output/DeltaS_EST_v" + std::to_string(version) + ".out");
    outfile << "# Mass of the host halo : " << M1 << " Msol" << std::endl;
    outfile << "# m/M [Msol] | MF (not normalised) [Msol^(-1)]" << std::endl;

    int Npts = 1000;
    double M2_min = 1e-5 * M1, M2_max = M1 / 2., dM2 = log10(M2_max / M2_min) / Npts, M2;

    double S1 = _power_spectrum.SigmaM2(M1, 0);

    for (int i = 0; i < Npts; i++)
    {
        M2 = M2_min * pow(10, i * dM2);
        outfile << M2 / M1 << "\t" << M1 * M1 / M2 * pow(_power_spectrum.SigmaM2(M2, 0) - S1, -1.5) * fabs(_power_spectrum.Der_sigmaM2_vs_M(M2, 0)) / sqrt(2 * PI) << "\t" << sqrt(2 / PI) * pow(_power_spectrum.SigmaM2(M2, 0) - S1, -0.5) << std::endl;
    }

    outfile.close();
}
//
// Interpolation functions

void MassFunction::Interpolate_M2_MF_EST_vs_M(double M1, double z1)
{
    int Npts = 1000;
    double dM = log10(1e+12) / Npts;
    std::vector<double> log10_M, M2_MF_EST;
    M2_MF_EST.resize(Npts);
    log10_M.resize(Npts);
    double M2 = 0;

    for (int i = 0; i < Npts; i++)
    {
        M2 = M1 * 1e-12 * pow(10, i * dM);
        log10_M[i] = log10(M2);
        M2_MF_EST[i] = M2 * M2 * MF_EST(M2, M1, 0);
    }

    _spline_M2_MF_EST_log10_M.set_points(log10_M, M2_MF_EST);
}

double MassFunction::Interp_M2_MF_EST_vs_M(double M2, double M1, double z1)
{
    return _spline_M2_MF_EST_log10_M(log10(M2)) / M2 / M2;
}

double MassFunction::Overdensity_BryanNorman_1998(double z)
{
    double Om_m_h2 = _cosmo.get_Omega_m_h2()*pow(1+z,3);
    double E2_h2_bis = Om_m_h2 + _cosmo.get_Omega_l_h2();
    
    double x = Om_m_h2/E2_h2_bis-1;

    //std::cout << x << " " << _cosmo.get_Omega_m_h2()*pow(1+z,3) << " " << _cosmo.get_Omega_l_h2() << std::endl;

    return 18*PI*PI + 82*x - 39*x*x;
}