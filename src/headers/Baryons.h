#ifndef BARYONS_H
#define BARYONS_H

#include "CppLib.h"
#include "mymath.h"
#include "myspline.h"
#include "MyUnits.h"

class Baryons
{
 public:

  Baryons();
  double bulge_density_mcmillan(double R, double z, double q, double r0,
				double r_cut, double rho_b, double alpha_b);
  double exponential_disc_density(double R, double z, double Rd, double zd, double Sigma_d);
  double gas_disc_density_mcmillan(double R, double z, double Rd, double Rm, double zd, double Sigma_d);
  double exponential_disc_surface_density(double R, double Rd, double Sigma_d);
  double gas_disc_surface_density_mcmillan(double R, double Rd, double Rm, double Sigma_d);

};

#endif
