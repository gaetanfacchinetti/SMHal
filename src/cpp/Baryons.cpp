#include "../headers/Baryons.h"

using namespace std;


Baryons::Baryons(){}

double Baryons::bulge_density_mcmillan(double R, double z, double q, double r0,
				       double r_cut, double rho_b, double alpha_b)
{
  double r = sqrt(R*R+z*z/(q*q));
  return rho_b*pow(1+r/r0,-alpha_b)*exp(-r*r/(r_cut*r_cut));
}

double Baryons::exponential_disc_density(double R, double z, double Rd, double zd, double Sigma_d)
{return 0.5*Sigma_d/zd*exp(-R/Rd-fabs(z)/zd);}

double Baryons::gas_disc_density_mcmillan(double R, double z, double Rd, double Rm, double zd, double Sigma_d)
{return 0.25*Sigma_d/zd*exp(-Rm/R-R/Rd)*pow(cosh(0.5*z/zd),-2);}

double Baryons::exponential_disc_surface_density(double R, double Rd, double Sigma_d)
{return Sigma_d*exp(-R/Rd);}

double Baryons::gas_disc_surface_density_mcmillan(double R, double Rd, double Rm, double Sigma_d)
{return Sigma_d*exp(-Rm/R-R/Rd);}
