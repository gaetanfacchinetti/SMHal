#include "../headers/PrimordialBlackHoles.h"

double PrimordialBlackHoles::MassHorizon_vs_Tformation(double Tform)
// Tform in GeV result in Msol
{
    double gamma = ps.VolumeForm();
    double pref = gamma * pow(3. / (8 * PI), 2) * sqrt(80. / PI);
    double g = degFree->get_gSmoothInterp(log10(Tform));

    return pref * pow(PLANCK_MASS, 3) / sqrt(g) * pow(Tform, -2) * GeV_to_MSOL;
}

double PrimordialBlackHoles::Tformation_vs_MassHorizon(double Mass)
{
    std::vector<double> TM = {0, Mass};
    double zeq = ps.get_Cosmology().Zeq_rad_mat();
    double a0OutOfaeq = zeq + 1;
    double Teq_GeV = T0_CMB_GeV * a0OutOfaeq;
    return DichotomieLn_Static(0, 2, TM, this, CallBack_MassHorizon_vs_Tformation_toSolve, Teq_GeV, 1e+10, 1e-6, 0);
}

double PrimordialBlackHoles::Radius_vs_MassHorizon(double Mass)
// Result in Mpc -> Radius = "comoving radius" here
{
    double Tform = Tformation_vs_MassHorizon(Mass); // Temperature in GeV
    double zeq = ps.get_Cosmology().Zeq_rad_mat();
    double a0OutOfaeq = zeq + 1;
    double Teq_GeV = T0_CMB_GeV * a0OutOfaeq;
    double aformOutOfaeq = (Teq_GeV / Tform) * pow(degFree->get_hSmoothInterp(log10(Teq_GeV)), 1. / 3.) / pow(degFree->get_hSmoothInterp(log10(Tform)), 1. / 3.);
    double aform = aformOutOfaeq / a0OutOfaeq;
    double Hform = pow(PLANCK_MASS, 2) / 2. / (Mass * MSOL_to_GeV);
    //double Hform_2 = 2 * sqrt((pow(PI, 3) / 45) * degFree->get_gSmoothInterp(log10(Tform))) * pow(Tform, 2) / M_PLANCK;

    //std::cout << Hform << " " << Hform_2 << std::endl;
    //std::cout << degFree->get_hSmoothInterp(log10(Teq_GeV)) << std::endl;
    double res = 1 / (aform * Hform) / GeV_to_mm1 * m_to_kpc * 1e-3;

    //std::cout << Mass << " " << res << std::endl;

    return res;
}

double PrimordialBlackHoles::MassHorizon_vs_Radius(double Radius)
{
    std::vector<double> MR = {0, Radius};
    return DichotomieLn_Static(0, 2, MR, this, CallBack_Radius_vs_MassHorizon_toSolve, 1e-18, 1e+14, 1e-6, 0);
}

double PrimordialBlackHoles::SigmaR(double R, double As, double ks)
// For a spiked power spectrum
// R and ks in the same units / inverse unit
{

    double sigma2 = 0;
    sigma2 = 16. / 81. * As * pow(R * ks, 4) * pow(TransferFunction_RDera(ks * R / sqrt(3)), 2) * pow(ps.Window(R * ks), 2);
    return sqrt(sigma2);
}

double PrimordialBlackHoles::TransferFunction_RDera(double y)
{
    if (y < 0.01)
        return 1.;
    else
        return 3. * (sin(y) - y * cos(y)) * pow(y, -3);
}

double PrimordialBlackHoles::Beta(double R, double alpha, double As, double ks)
{
    double nu = 0;
    double sig = SigmaR(R, As, ks);

    //sig = 0.08;

    if (sig > 0)
        nu = 0.3 / sig;
    else
        return 0;

    //std::cout << "nu = " << nu << std::endl;

    if (nu > 1000)
    {
        //std::cout << alpha << " " << (sqrt(2.*PI)*nu) << " " << exp(-nu*nu/2.) << std::endl;
        return 2. * alpha / (sqrt(2. * PI) * nu) * exp(-nu * nu / 2.);
    }
    else
        return alpha * erfc(nu / sqrt(2));
}

double PrimordialBlackHoles::fPBH(double Mass, double alpha, double As, double ks)
// Here the mass input mass is not the horizon mass but the true mass
// ks in Mpc, Mass in Msol
{
    double MH = Mass / alpha;             // We first recover the horizon mass
    double R = Radius_vs_MassHorizon(MH); // in Mpc
    double b = Beta(R, alpha, As, ks);
    //std::cout << "Beta  = " << b << " " << R << std::endl;

    double Tform = Tformation_vs_MassHorizon(MH); // Temperature in GeV
    double zeq = ps.get_Cosmology().Zeq_rad_mat();
    double a0OutOfaeq = zeq + 1;
    double Teq_GeV = T0_CMB_GeV * a0OutOfaeq;
    double aformOutOfaeq = (Teq_GeV / Tform) * pow(degFree->get_hSmoothInterp(log10(Teq_GeV)), 1. / 3.) / pow(degFree->get_hSmoothInterp(log10(Tform)), 1. / 3.);
    double aform = aformOutOfaeq / a0OutOfaeq;
    double Hform = pow(PLANCK_MASS, 2) / 2. / (MH * MSOL_to_GeV);
    double H0 = 100 * ps.get_Cosmology().get_h() / kpc_to_m / s_to_GeVm1; // GeV
    double OmegaDM = ps.get_Cosmology().cosmic_abundance(0, Species::CDM);

    //std::cout << 1./(1.1e-8)*pow(MH, -1./2.) << " " << pow(Hform/H0, 2)*pow(aform,3) << std::endl;
    return pow(Hform / H0, 2) * pow(aform, 3) * b / OmegaDM;
}

double PrimordialBlackHoles::Max_fPBH(double alpha, double As, double ks)
{
    double Ms_PBH = 0.2*MassHorizon_vs_Radius(1./ks);

    int Npts = 500;
    double Mmin = std::max(Ms_PBH*1e-1, 1e-18), Mmax = std::min(1e+14, Ms_PBH*1e+2), dM = log10(Mmax / Mmin) / (1.0 * Npts), M;

    double fPBH_max = 0;
    double val = 0;

    for (int i = 0; i < Npts + 1; i++)
    {
        M = Mmin * pow(10, i * dM);
        val = fPBH(M, alpha, As, ks);
        
        if( val > fPBH_max)
            fPBH_max = val;
    }

    return fPBH_max;
}

double PrimordialBlackHoles::Max_As(double max_f, double alpha, double ks)
{
    int Npts = 100;
    double As_min = 1e-4, As_max = 1e+2, dAs = log10(As_max / As_min) / (1.0 * Npts), As;

    double val = 0;
    double As_min_2 = 0;
    int i = 0;

    while(As_min_2 == 0 && i < Npts+1)
    {
        As = As_min * pow(10, i * dAs);
        val = Max_fPBH(alpha, As, ks);
        
        if(val > 0)
            As_min_2 = As;

        //std::cout << As << std::endl;

        i++;
    }

    std::vector<double> xx={0, max_f, alpha, ks};
    return DichotomieLn_Static(0, 4, xx, this, CallBack_Max_fPBH_toSolve, As_min_2, As_max, 1e-6,  0);

}


void PrimordialBlackHoles::plot_MassHorizon_vs_Radius()
{
    std::ofstream outfile;
    outfile.open("../output/PBH/Mass_Vs_Radius.out");

    int Npts = 100;
    double Rmin = 1e-18, Rmax = 1e-4, dR = log10(Rmax / Rmin) / (1.0 * Npts), R;

    double M_approx = 0;

    for (int i = 0; i < Npts + 1; i++)
    {
        R = Rmin * pow(10, i * dR);
        M_approx = 5.6e+15 * pow(R * 0.05e-3, 2);
        outfile << R << "\t" << MassHorizon_vs_Radius(R) << "\t" << M_approx << "\t" << SigmaR(R * 1e-3, 1e-5, 1e+3) << std::endl;
    }

    outfile.close();
}

void PrimordialBlackHoles::plot_fPBH(double As, double ks)
{
    std::ofstream outfile;
    outfile.open("../output/PBH/fPBH.out");

    double Ms_PBH = 0.2*MassHorizon_vs_Radius(1./ks);
    std::cout << Ms_PBH << std::endl;
    //exit(0);
    int Npts = 500;
    double Mmin = std::max(Ms_PBH*1e-1, 1e-18), Mmax = std::min(1e+14, Ms_PBH*1e+3), dM = log10(Mmax / Mmin) / (1.0 * Npts), M;
    

    outfile << "# As = " << As << std::endl;
    outfile << "# ks = " << ks << " Mpc^{-1}" << std::endl;
    outfile << "# M_PBH [Msol] | f_PBH" << std::endl;

    for (int i = 0; i < Npts + 1; i++)
    {
        M = Mmin * pow(10, i * dM);
        outfile << M << "\t" << fPBH(M, 0.2, As, ks) << "\t" << SigmaR(Radius_vs_MassHorizon(M/0.2),As, ks) << std::endl;
    }

    outfile.close();
}

void PrimordialBlackHoles::plot_Max_As()
{
    std::ofstream outfile;
    outfile.open("../output/PBH/Max_As.out");

    int Npts = 200;
    double ks_min = 1., ks_max = 1e+15, dk = log10(ks_max / ks_min) / (1.0 * Npts), ks;
    

    outfile << "# ks [Mpc^{-1}] | f_PBH" << std::endl;


    for (int i = 0; i < Npts + 1; i++)
    {
        ks = ks_min * pow(10, i * dk);
        outfile << ks << "\t" << Max_As(1e-1, 0.2, ks) << "\t" << Max_As(1e-10, 0.2, ks) << std::endl;
    }

    outfile.close();
}