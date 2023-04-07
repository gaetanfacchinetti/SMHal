#include "../headers/PowerSpectrum.h"

PowerSpectrum::PowerSpectrum(Cosmology n_cosmo, PSWindowType w, bool n_normalised_s8, std::shared_ptr<TransferFunction> n_add_TF)
{
    cosmo = n_cosmo;
    window_type = w;

    TF = std::make_shared<TransferFunction_EH98>(cosmo);
    add_TF = n_add_TF;

    if(add_TF && (TF->get_cosmo() != add_TF->get_cosmo()))
    {
        std::cout << "ERROR: the additional TF in the power spectrum must have the same cosmology than n_cosmo" << std::endl;
        exit(0);
    }

    // Additional spike not considered first
    As_spike.resize(0);
    ks_spike.resize(0);
    eps_spike.resize(0);
    eps_spike_0.resize(0);

    // Parameters for the window function
    gamma_rs_th = 4 * PI / 3.;
    gamma_gaussian = pow(2 * PI, 3. / 2.);
    gamma_fs_th = 6 * PI * PI;

    // Test
    //gamma_fs_th = 4 * PI / 3.;

    ns = cosmo.get_ns();
    As = cosmo.get_As();

    // Planck normalisation
    k0 = 0.05; // Mpc^{-1}

    is_interpolated_SigmaM2 = false;
    _z_interp = 0;

    // We renormalise to the value of sigma_8 if that is asked and if possible
    norm_PS = 1.;
    double s8 = 0;
    if (n_normalised_s8 == true)
    {
        window_type = PSWindowType::real_space_top_hat;
        s8 = SigmaR(8. / cosmo.get_h(), 0);
        norm_PS = pow(cosmo.get_sigma8() / s8, 2.);
        window_type = w;
        //std::cout << norm_PS << " " << cosmo.get_sigma8() << " " << s8 << std::endl;
        //std::cout << SigmaR(8. / cosmo.get_h(), 0) << std::endl;
    }
}

void PowerSpectrum::set_spikes(std::vector<double> n_as_spike, std::vector<double> n_ks_spike, std::vector<double> n_eps_spike)
{
    As_spike = n_as_spike;
    ks_spike = n_ks_spike;
    eps_spike = n_eps_spike;

    if ((As_spike.size() != ks_spike.size()) || (As_spike.size() != eps_spike.size()) || (ks_spike.size() != eps_spike.size()))
    {
        std::cout << "ERROR : all vectors must have the same size in " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    eps_spike_0.resize(0);

    for (int i = 0; i < eps_spike.size(); i++)
        if (eps_spike[i] == 0)
        {
            eps_spike_0.push_back(i);
            //std::cout << "coucou" << std::endl;
        }
}

double PowerSpectrum::regulator_spikes(double k, double ks, double eps)
// Regulated by log-normal distrubution with dimension of [Mpc] for k in Mpc^{-1}
{
    return 1. / (k * eps * sqrt(2 * PI)) * exp(-pow(log(k) - log(ks), 2) / 2 / pow(eps, 2));
}

// ATTENTION :: k must be in Mpc^-1
double PowerSpectrum::dimensionless_curvature_power_spectrum(double k)
{
    double res = 0;
    res = As * pow(k / k0, ns - 1);

    // We can add spikes to the power spectrum is we want
    for (int i = 0; i < As_spike.size(); i++)
        if (eps_spike[i] > 0) // we do not consider the extreme case eps = 0
            res += As_spike[i] * ks_spike[i] * regulator_spikes(k, ks_spike[i], eps_spike[i]);

    return res;
}

// ATTENTION :: k must be in Mpc^-1
double PowerSpectrum::curvature_power_spectrum(double k)
{
    return 2 * PI * PI * dimensionless_curvature_power_spectrum(k) * pow(k, -3);
}

// ATTENTION :: k must be in Mpc^-1
double PowerSpectrum::matter_power_spectrum(double k, double z)
{
 
    double h = cosmo.get_h();

    //std::cout << TF << std::endl;

    //TF->test();
    //TransferFunction_EH98 *TF_test = new TransferFunction_EH98(cosmo);
    double tf = TF->evaluate(k);

    if (add_TF) 
        tf = tf * add_TF->evaluate(k);
 
    double D1_z = cosmo.growth_factor_D1_Carroll(z);
    double Omega_m_0 = cosmo.get_Omega_m_h2() / (h * h);
    double c_over_H0_Mpc = 1e-3 * LIGHT_SPEED / (100 * h);
    //return 0;

    // Formula from Nakama+18
    return norm_PS * (4. / 25.) * pow(D1_z * k * k * tf * c_over_H0_Mpc * c_over_H0_Mpc / Omega_m_0, 2) * curvature_power_spectrum(k);
    // Norm_PS is here to reajust to the value of sigma 8 if asked
}



// ATTENTION :: k must be in Mpc^-1
double PowerSpectrum::dimensionless_matter_power_spectrum(double k, double z)
{
    return matter_power_spectrum(k, z) * pow(k, 3) / (2 * PI * PI);
}

void PowerSpectrum::plot_PowerSpectrum(double z, std::string add)
{
    std::ofstream outfile;
    std::string name_file = "../output/PowerSpectrum_";
    
    if (add != "")
        name_file += "_";

    outfile.open(name_file + "_.out");
    outfile << "# z = " << z << std::endl;
    outfile << "# k [Mpc^{-1}] | Matter_Power_Spectrum | Dimensionless_Matter_Power_Spectrum | Curvature power spectrum | Dimensionless curvature power spectrum | Matter_Power_Spectrum with TF = 1 | Transfer function " << std::endl;

    int N = 1000;
    double kmax = 1e+1, kmin = 1e-4;
    //double kmax = 1e+6, kmin = 1;
    double dk = log10(kmax / kmin) / N, k;

    for (int i = 0; i < N + 1; i++)
    {
        k = kmin * pow(10, i * dk);
        double tf = TF->evaluate(k);
        outfile << k << "\t" << matter_power_spectrum(k, z) << "\t" << dimensionless_matter_power_spectrum(k, z)
                << "\t" << curvature_power_spectrum(k) << "\t" << dimensionless_curvature_power_spectrum(k)
                << "\t" << matter_power_spectrum(k, z) / (tf * tf) << "\t" << tf << std::endl;
    }

    outfile.close();
}

double PowerSpectrum::Window(double kR)
{
    if (window_type == PSWindowType::real_space_top_hat)
    {
        if(kR < 0.001)
            return 1. - pow(kR,2)/10.;
        else
            return 3. * pow(kR, -3) * (sin(kR) - kR * cos(kR));
    }
    if (window_type == PSWindowType::gaussian)
        return exp(-kR * kR / 2.);
    if (window_type == PSWindowType::fourier_space_top_hat)
    {
        if (1 - kR > 0)
            return 1.;
        else
            return 0.;
    }

    std::cout << "WARNING :: Wrong type of window in -> " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}

double PowerSpectrum::Der_Window_vs_kR(double kR)
{
    if (window_type == PSWindowType::real_space_top_hat)
        return kR < 0.1 ?  -kR/5. + pow(kR, 3)/70. : 3 * sin(kR) * pow(kR, -2) - 9 * (-kR * cos(kR) + sin(kR)) * pow(kR, -4);
    if (window_type == PSWindowType::gaussian)
        return -kR * exp(-kR * kR / 2.);
    if (window_type == PSWindowType::fourier_space_top_hat)
        return 0;

    std::cout << "WARNING :: Wrong type of window in -> " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}

double PowerSpectrum::LagrangianRadius_vs_Mass(double Mass)
// Mass in Msol and result in Mpc
{
    double rho_bar = 1e+9 * (cosmo.cosmological_density(0, Species::MATTER)); // [Msol*Mpc^-3]

    if (window_type == PSWindowType::real_space_top_hat)
        return pow(Mass / (rho_bar * gamma_rs_th), 1. / 3.);
    if (window_type == PSWindowType::gaussian)
        return pow(Mass / (rho_bar * gamma_gaussian), 1. / 3.);
    if (window_type == PSWindowType::fourier_space_top_hat)
        return pow(Mass / (rho_bar * gamma_fs_th), 1. / 3.);

    std::cout << "WARNING :: Wrong type of window in -> " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}

double PowerSpectrum::Mass_vs_LagrangianRadius(double Radius)
// Radius in Mpc^{-3} and result in solar masses
{
    double rho_bar = 1e+9 * cosmo.cosmological_density(0, Species::MATTER); // [Msol*Mpc^-3]

    if (window_type == PSWindowType::real_space_top_hat)
        return rho_bar * gamma_rs_th * pow(Radius, 3);
    if (window_type == PSWindowType::gaussian)
        return rho_bar * gamma_gaussian * pow(Radius, 3);
    if (window_type == PSWindowType::fourier_space_top_hat)
        return rho_bar * gamma_fs_th * pow(Radius, 3);

    std::cout << "WARNING :: Wrong type of window in -> " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}

double PowerSpectrum::Der_Mass_vs_LagrangianRadius(double Mass) // Result in [Msol.Mpc^{-1}]
{
    double rho_bar = 1e+9 * cosmo.cosmological_density(0, Species::MATTER); // [Msol*Mpc^-3]

    if (window_type == PSWindowType::real_space_top_hat)
        return 3 * pow(rho_bar * gamma_rs_th * Mass * Mass, 1. / 3.);
    if (window_type == PSWindowType::gaussian)
        return 3 * pow(rho_bar * gamma_gaussian * Mass * Mass, 1. / 3.);
    if (window_type == PSWindowType::fourier_space_top_hat)
        return 3 * pow(rho_bar * gamma_fs_th * Mass * Mass, 1. / 3.);

    std::cout << "WARNING :: Wrong type of window in -> " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}

double PowerSpectrum::fForDichotomie_Mass_vs_Sigma2(double M, double sM2, double z)
{
    double res = Interp_SigmaM2(M, z) - sM2;
    //std::cout << M << " " << res << std::endl;
    return res;
}

double PowerSpectrum::Mass_vs_Sigma2(double sM2, double z)
{
    std::vector<double> xx = {0, sM2, z};

    if (sM2 > SigmaM2(1e-12, 0)) // The mass of the halo is so small that it sould not even exist
        return 0;

    return DichotomieLn_Static(0, 3, xx, this, CallBack_fForDichotomie_Mass_vs_SigmaM2, 1e-12, 1e+15, 1e-10, 0);
}

double PowerSpectrum::Der_sigmaR2_vs_R(double R, double z) // Result in Mpc^{-1}
// Attention here we do not incorporate the presence of other peaks
{
    std::vector<double> kRz;
    kRz.resize(3);
    kRz[1] = R;
    kRz[2] = z;

    double k_min = exp(-7.), k_max;

    if (window_type == PSWindowType::real_space_top_hat)
        k_max = 20.0 / R;

    if (window_type == PSWindowType::gaussian)
        k_max = 4.0 / R;

    //double err;
    double res = -1;

    if (window_type == PSWindowType::real_space_top_hat)
        res = GaussLegendre_IntegralLn_Static(0, 3, gl500, k_min, k_max, kRz, this, CallBack_fToIntFor_Der_sigmaR2_vs_R);
    if (window_type == PSWindowType::gaussian)
        res =  GaussLegendre_IntegralLn_Static(0, 3, gl500, k_min, k_max, kRz, this, CallBack_fToIntFor_Der_sigmaR2_vs_R);
    if (window_type == PSWindowType::fourier_space_top_hat)
        res= -dimensionless_matter_power_spectrum(1. / R, z) / R;

    //std::cout << k_min << " " << k_max << " Here we have: " << res << std::endl; 

    return res;
}

double PowerSpectrum::fToIntFor_Der_sigmaR2_vs_R(std::vector<double> kRz)
{
    double k = kRz[0], R = kRz[1], z = kRz[2];
    double res = 2 * dimensionless_matter_power_spectrum(k, z) * Der_Window_vs_kR(k * R) * Window(k * R);
    //std::cout << k*R << " " << res  << " " << Window(k * R) << " " << Der_Window_vs_kR(k * R) << " " <<  dimensionless_matter_power_spectrum(k, z) << std::endl;
    return res;
}

double PowerSpectrum::SigmaR2(double R, double z)
// !! R in Mpc
{


    double lnk_min = -8, lnk_max;
    double res = 0;

    if (window_type == PSWindowType::real_space_top_hat)
        lnk_max = log(20) - log(R); // lnk_max = 4.0 - log(R);
    if (window_type == PSWindowType::gaussian)
        lnk_max = log(4) - log(R); // lnk_max = 1.2 - log(R);
    if (window_type == PSWindowType::fourier_space_top_hat)
        lnk_max = -log(R);

    if (lnk_max <= lnk_min)
        return 0;

    std::vector<double> lnkRz;
    lnkRz.resize(3);
    lnkRz[1] = R;
    lnkRz[2] = z;
    /*
    double err;

       {
        res = Simpson_Integral1_Static(0, 3, 0, 500, lnk_min, lnk_max, lnkRz, this, CallBack_fToIntFor_Sigma2R, err);
    }*/

    if (window_type == PSWindowType::real_space_top_hat)
        res = GaussLegendre_Integral_Static(0, 3, gl2000, lnk_min, lnk_max, lnkRz, this, CallBack_fToIntFor_Sigma2R);
    if (window_type == PSWindowType::gaussian || window_type == PSWindowType::fourier_space_top_hat)
        res = GaussLegendre_Integral_Static(0, 3, gl500, lnk_min, lnk_max, lnkRz, this, CallBack_fToIntFor_Sigma2R);

    //std::cout << "size = " << eps_spike_0.size() << std::endl;

    if (eps_spike_0.size() > 0)
        for (int j = 0; j < eps_spike_0.size(); j++)
            res += SigmaR2_narrow_spike(R, As_spike[eps_spike_0[j]], ks_spike[eps_spike_0[j]], z);

    return res;
}



double PowerSpectrum::fToIntFor_Sigma2R(std::vector<double> lnkRz)
{

    
    double lnk = lnkRz[0];
    double R = lnkRz[1];
    double z = lnkRz[2];


    double k = exp(lnk);
    return dimensionless_matter_power_spectrum(k, z) * pow(Window(k * R), 2);
}

void PowerSpectrum::Interpolate_SigmaM2_vs_M(double z)
{
    int Npts = 5000;
    double dM = log10(1e+16 / 1e-12) / Npts;
    std::vector<double> log10_M, SM2;
    log10_M.resize(Npts);
    SM2.resize(Npts);

    double M;

    for (int i = 0; i < Npts; i++)
    {
        M = 1e-12 * pow(10, i * dM);
        log10_M[i] = log10(M);
        SM2[i] = SigmaM2(M, z);
        //std::cout << pow(10, log10_M[i]) << " " << SM2[i] << std::endl;
    }

    spline_sigmaM2_vs_log10_M.set_points(log10_M, SM2);
}

/**
 * Return the result of the interpolation
 * but if the redshift is the same and the interpolation has been done
*/
double PowerSpectrum::Interp_SigmaM2(double M, double z)
{

    if (M < 1e-12)
    {
        std::cout << "FATAL ERROR : Impossible to compute Interp_SigmaM2 in " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    if (z == _z_interp && is_interpolated_SigmaM2 == true)
        return spline_sigmaM2_vs_log10_M(log10(M));

    //std::cout << "here " << z << " " << _z_interp << std::endl;
    Interpolate_SigmaM2_vs_M(z);
    is_interpolated_SigmaM2 = true;
    _z_interp = z;

    //std::cout << spline_sigmaM2_vs_log10_M(log10(M)) << std::endl;

    return spline_sigmaM2_vs_log10_M(log10(M));
}

double PowerSpectrum::VolumeForm()
{
    double volume_form = -1;

    if (window_type == PSWindowType::real_space_top_hat)
        volume_form = gamma_rs_th;
    if (window_type == PSWindowType::gaussian)
        volume_form = gamma_gaussian;
    if (window_type == PSWindowType::fourier_space_top_hat)
        volume_form = gamma_fs_th;

    return volume_form;
}

//
//
// Fonction plot
//

void PowerSpectrum::plot_Sigma_vs_RM(double z, bool with_respect_to_M, bool in_units_of_h, std::string add)
{
    std::ofstream outfile;
    std::string name_file;

    if (with_respect_to_M == true)
        name_file = "../output/UCMHs/SigmaM_" + add;
    else
        name_file = "../output/UCMHs/SigmaR_" + add;

    if (add != "")
        name_file += "_";

    if (window_type == PSWindowType::real_space_top_hat)
        name_file += "rsth.out";
    if (window_type == PSWindowType::gaussian)
        name_file += "gaussian.out";
    if (window_type == PSWindowType::fourier_space_top_hat)
        name_file += "fsth.out";

    outfile.open(name_file);

    double Npoints = 1000;

    double h_param = 1.;
    if (in_units_of_h == true)
        h_param = cosmo.get_h();

    double xmin = 0, xmax = 0;

    if (with_respect_to_M == true)
    {
        xmin = 1e-12 * h_param;
        xmax = 1e+18 * h_param;
    }
    else
    {
        // We convert the mass range into a range in R [in Mpc]
        xmin = 3 * pow(3 * 1e-19 / (4 * PI), 1. / 3.) * pow(cosmo.critical_density(0), -1. / 3.) * 1e-3 * h_param;
        xmax = 3 * pow(3 * 1e+18 / (4 * PI), 1. / 3.) * pow(cosmo.critical_density(0), -1. / 3.) * 1e-3 * h_param;
    }

    double x, d_logx = log10(xmax / xmin) / Npoints;

    if (in_units_of_h == false)
    {
        if (with_respect_to_M == true)
            outfile << "# M [Msol] | sigmaM" << std::endl;
        else
            outfile << "# R [Mpc]  | sigmaM" << std::endl;
    }
    else
    {
        if (with_respect_to_M == true)
            outfile << "# M [Msol/h] | sigmaR" << std::endl;
        else
            outfile << "# R [Mpc/h]  | sigmaR" << std::endl;
    }

    for (int i = 0; i < Npoints + 1; i++)
    {
        x = xmin * pow(10, i * d_logx);

        if (with_respect_to_M == true)
            outfile << x << "\t" << SigmaM(x / h_param, z) << "\t" << sqrt(Interp_SigmaM2(x / h_param, z)) << std::endl;
        else
            outfile << x << "\t" << SigmaR(x / h_param, z) << "\t" << sqrt(Interp_SigmaM2(x / h_param, z)) << std::endl;
        //outfile << M << " " <<  << std::endl;
    }

    outfile.close();
}

void PowerSpectrum::plot_SigmaM_vs_z(double M)
{
    std::ofstream outfile;
    std::string name_file;

    if (window_type == PSWindowType::real_space_top_hat)
        name_file = "../output/SigmaM_rsth_vs_z.out";
    if (window_type == PSWindowType::gaussian)
        name_file = "../output/SigmaM_gaussian_vs_z.out";
    if (window_type == PSWindowType::fourier_space_top_hat)
        name_file = "../output/SigmaM_fsth_vs_z.out";

    outfile.open(name_file);

    double zmin = 1e-1;
    double zmax = 1e+3;

    double Npoints = 500;

    double z, d_logz = log10(zmax / zmin) / Npoints;

    outfile << "# z | sigmaM" << std::endl;

    for (int i = 0; i < Npoints + 1; i++)
    {
        z = zmin * pow(10, i * d_logz);
        outfile << z << " " << SigmaM(M, z) << std::endl;
    }

    outfile.close();
}

void PowerSpectrum::plot_Der_LnSigma_vs_LnRM(double z, bool with_respect_to_M, bool in_units_of_h)
{
    std::ofstream outfile;
    std::string name_file;

    if (with_respect_to_M == true)
        name_file = "../output/Der_LnSigmaM_lnM";
    else
        name_file = "../output/Der_LnSigmaR_lnR";

    if (window_type == PSWindowType::real_space_top_hat)
        name_file += "_rsth.out";
    if (window_type == PSWindowType::gaussian)
        name_file += "_gaussian.out";
    if (window_type == PSWindowType::fourier_space_top_hat)
        name_file += "_fsth.out";

    outfile.open(name_file);

    double h_param = 1.;
    if (in_units_of_h == true)
        h_param = cosmo.get_h();

    double xmin = 0, xmax = 0;

    if (with_respect_to_M == true)
    {
        xmin = 1e-10 * h_param;
        xmax = 1e+18 * h_param;
    }
    else
    {
        // We convert the mass range into a range in R [in Mpc]
        xmin = 3 * pow(3 * 1e-10 / (4 * PI), 1. / 3.) * pow(cosmo.critical_density(0), -1. / 3.) * 1e-3 * h_param;
        xmax = 3 * pow(3 * 1e+18 / (4 * PI), 1. / 3.) * pow(cosmo.critical_density(0), -1. / 3.) * 1e-3 * h_param;
    }

    double Npoints = 1000;

    double x, d_logx = log10(xmax / xmin) / Npoints;

    if (in_units_of_h == false)
    {
        if (with_respect_to_M == true)
            outfile << "# M [Msol] | abs(Der_LnSigmaM_vs_lnM)" << std::endl;
        else
            outfile << "# R [Mpc]  | Der_LnSigmaM_vs_lnR" << std::endl;
    }
    else
    {
        if (with_respect_to_M == true)
            outfile << "# M [Msol/h] | abs(Der_LnSigmaM_vs_lnM)" << std::endl;
        else
            outfile << "# R [Mpc/h]  | Der_LnSigmaM_vs_lnR" << std::endl;
    }

    for (int i = 0; i < Npoints + 1; i++)
    {
        x = xmin * pow(10, i * d_logx);

        if (with_respect_to_M == true)
            outfile << x << "\t" << x / h_param * abs(Der_LnSigmaM_vs_M(x / h_param, z)) << std::endl;
        else
            outfile << x << "\t" << x / h_param * abs(Der_LnSigmaR_vs_R(x / h_param, z)) << std::endl;
    }

    outfile.close();
}

void PowerSpectrum::plot_fToIntFor_Sigma2R(double z, double R)
{
    std::ofstream outfile;
    std::string name_file;

    if (window_type == PSWindowType::real_space_top_hat)
        name_file = "../output/fToIntFor_Sigm2R_rsth.out";
    if (window_type == PSWindowType::gaussian)
        name_file = "../output/fToIntFor_Sigma2R_gaussian.out";
    if (window_type == PSWindowType::fourier_space_top_hat)
        name_file = "../output/fToIntFor_Sigma2R_fsth.out";

    outfile.open(name_file);

    outfile << "# R [Mpc] | fToIntFor_SigmaR2 [...]" << std::endl;
    int Npts = 200;
    double kmin = exp(-6), kmax = exp(4) / R, d_logk = log10(kmax / kmin) / Npts, k;

    std::vector<double> lnkRz{0, R, z};

    for (int i = 0; i < Npts + 1; i++)
    {
        k = kmin * pow(10, i * d_logk);
        lnkRz[0] = log(k);
        outfile << k << "\t" << fToIntFor_Sigma2R(lnkRz) << std::endl;
    }

    outfile.close();
}

double PowerSpectrum::SigmaR2_narrow_spike_over_window(double As_s, double ks_s, double z)
{
    double h = cosmo.get_h();
    double D1_z = cosmo.growth_factor_D1_Belloso(z);
    double Omega_m_0 = cosmo.get_Omega_m_h2() / (h * h);
    double c_over_H0_Mpc = 1e-3 * LIGHT_SPEED / (100 * h);
    double tf = TF->evaluate(ks_s);

    double C = norm_PS * (4. / 25.) * pow(D1_z * c_over_H0_Mpc * c_over_H0_Mpc / Omega_m_0, 2); // Mpc^{-4}

    return C * tf * tf * pow(ks_s, 4) * As_s;
}

double PowerSpectrum::SigmaR2_narrow_spike(double R, double As_s, double ks_s, double z)
{
    return SigmaR2_narrow_spike_over_window(As_s, ks_s, z) * pow(Window(ks_s*R), 2.);
}


