#ifndef MASS_MODELS_
#define MASS_MODELS_

#include "CppLib.h"
#include "mymath.h"
#include "myspline.h"
#include "MyUnits.h"
#include "DarkHalo.h"
#include "Baryons.h"


class MassModel
{
 private:

  int m_model_dm;
  int m_model_baryons;
  
 public:

  MassModel();
  MassModel(int model_dm, int model_baryons);
  // dm
  double dm_alpha() const;
  double dm_beta() const;
  double dm_gamma() const;
  double dm_scale_density() const; // Msun/kpc^3
  double dm_scale_radius() const; // kpc
  // stellar discs
  double thick_disc_Sigma_d() const; // Msun/kpc^2
  double thick_disc_Rd() const; // kpc
  double thick_disc_zd() const; // kpc
  double thin_disc_Sigma_d() const; // Msun/kpc^2
  double thin_disc_Rd() const; // kpc
  double thin_disc_zd() const; // kpc
  // gas discs
  double h1_disc_Sigma_d() const; // Msun/kpc^2
  double h1_disc_Rm() const; // kpc
  double h1_disc_Rd() const; // kpc
  double h1_disc_zd() const; // kpc
  double h2_disc_Sigma_d() const; // Msun/kpc^2
  double h2_disc_Rm() const; // kpc
  double h2_disc_Rd() const; // kpc
  double h2_disc_zd() const; // kpc
  double gas_disc_Sigma_d() const; // Msun/kpc^2
  double gas_disc_Rd() const; // kpc
  double gas_disc_zd() const; // kpc
  // bulge
  double bulge_rho_b() const; // Msun/kpc^3
  double bulge_q() const;
  double bulge_r0() const; // kpc
  double bulge_r_cut() const; // kpc
  double bulge_alpha_b() const;
  double R_max() const; // kpc

  double bulge_density_mcmillan(double R, double z, double q, double r0,
				double r_cut, double rho_b, double alpha_b);
  double exponential_disc_density(double R, double z, double Rd, double zd, double Sigma_d);
  double gas_disc_density_mcmillan(double R, double z, double Rd, double Rm, double zd, double Sigma_d);
  double exponential_disc_surface_density(double R, double Rd, double Sigma_d);
  double gas_disc_surface_density_mcmillan(double R, double Rd, double Rm, double Sigma_d);
  double miyamoto_nagai_disc_density(double R, double z, double a, double b, double M);
  double miyamoto_nagai_disc_potential(double R, double z, double a, double b, double M);  // [(km/s)^]
};

#endif
