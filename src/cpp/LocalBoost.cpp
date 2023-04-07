#include "../headers/LocalBoost.h"

double LocalBoost::Luminosity_OneSubhalo(double m200, double c200, double R)
// Result in Msol^2/kpc
{
    double eta = subModel.rm2_rs();
    double xt = fslmod.rt_over_rs_DMO(R, c200);
    return 200 * fslmod.get_cosmo().critical_density(0) / 3. * m200 * pow(c200, 3)* DimensionlessLuminosity(xt) * pow(subModel.mass_profile(eta * c200), -2);
}

double LocalBoost::DimensionlessLuminosity(double xt)
{
    double alpha = subModel.get_alpha();
    double beta = subModel.get_beta();
    double gamma = subModel.get_gamma();

    if (alpha == 1 && beta == 3 && gamma == 1)
    {
        //if(xt < 1)
        return (-1. / 3.) * (pow(1 + xt, -3) - 1);
        //else
        //  return (1./3.) * (7./8.);
    }
    else
    {
        std::cout << "FATAL ERROR : trying to compute luminosity for undefined parameters in " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }
}


double LocalBoost::fToIntOnc200_LuminositySubhalos(double c200, double m200, double R)
{
    double lum = Luminosity_OneSubhalo(m200, c200, R);
    return lum * fslmod.ProbDistributionConcentration(c200, m200) * fslmod.ProbDistributionMass(m200) * fslmod.ProbDistributionSpace(R);
}

double LocalBoost::fToIntOnm200_LuminositySubhalos(double m200, double R)
{
    std::vector<double> xx = {0, m200, R};
    double epst = fslmod.get_epst();
    double c200_min = fslmod.c200_min_DMO(R, epst);
    //std::cout << "here : " << m200 << std::endl;
    return GaussLegendre_IntegralLn_Static(0, 3, gl200, c200_min, fslmod.get_c200_max(), xx, this, CallBack_fToIntOnc200_LuminositySubhalos);
}

double LocalBoost::LuminositySubhalos(double R)
{
    std::vector<double> xx = {0, R};
    double NsubOverKt = fslmod.NsubUnevolved();
    //std::cout << fslmod.get_m200_min() << " " << fslmod.get_m200_max() << std::endl;
    return NsubOverKt*GaussLegendre_IntegralLn_Static(0, 2, gl200, fslmod.get_m200_min(), fslmod.get_m200_max(), xx, this, CallBack_fToIntOnm200_LuminositySubhalos);
}


double LocalBoost::LuminositySmooth(double R)
// Result in Msol^2/kpc^6
{
    DarkHalo host = fslmod.get_hostHalo();
    return pow(host.dm_density(R)- fslmod.averageRhoSub(R), 2); 
}

double LocalBoost::LuminositySmoothSubhalos(double R)
// Result in Msol^2 /kpc^6
{
    DarkHalo host = fslmod.get_hostHalo();
    return 2*(host.dm_density(R)- fslmod.averageRhoSub(R))*fslmod.averageRhoSub(R);
}

double LocalBoost::Boost(double R)
{
    return TotalLuminosityWithSubhalos(R)/TotalLuminosityWithoutSubhalos(R);
}


double LocalBoost::TotalLuminosityWithSubhalos(double R)
// Result in Msol^2 /kpc^6
{
    std::cout << "Here = " << R << std::endl;
    return LuminositySubhalos(R) + LuminositySmooth(R) + LuminositySmoothSubhalos(R);
}


double LocalBoost::TotalLuminosityWithoutSubhalos(double R)
// Result in Msol^2 /kpc^6
{
    DarkHalo host = fslmod.get_hostHalo();
    return pow(host.dm_density(R), 2);
}



double LocalBoost::TotalBoost()
{   
    std::vector<double> xx = {0};
    double Rmax = fslmod.get_R_max();
    double TotLumWith = GaussLegendre_IntegralLn_Static(0, 1, gl100, 1e-1, Rmax, xx, this, CallBack_ftoIntOnR_TotalLuminosityWithSubhalos);
    double TotLumWithout = GaussLegendre_IntegralLn_Static(0, 1, gl100, 1e-1, Rmax, xx, this, CallBack_ftoIntOnR_TotalLuminosityWithoutSubhalos);

    return TotLumWith/TotLumWithout;
}


