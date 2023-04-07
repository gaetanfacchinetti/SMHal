#ifndef GALAXY_
#define GALAXY_

#include "DarkHalo.h"
#include "Baryons.h"
#include "MassModels.h"
#include "CppLib.h"
#include "mymath.h"
#include "myspline.h"
#include "MyUnits.h"

class Galaxy : public DarkHalo
{
protected:
  MassModel model;
  int m_model_dm;
  int m_model_baryons;

public:
  Galaxy(Cosmology cosmo);
  Galaxy(Cosmology cosmo, int model_dm, int model_baryons);
  double get_rmax() const;                      // [kpc]
  double total_ISM_density(double R, double z); // [Msun/kpc^3]
  double HI_density(double R, double z);        // [Msun/kpc^3]
  double H2_density(double R, double z);        // [Msun/kpc^3]
  double der_R_HI_density(double R, double z);
  double der_z_HI_density(double R, double z);
  double der_R_H2_density(double R, double z);
  double der_z_H2_density(double R, double z);
  double total_baryon_density(double R, double z); // [Msun/kpc^3]
  double disc_surface_density(double R);           // [Msun/kpc^2]
  double integrand_spherical_density(std::vector<double> z_r);
  void compute_spherical_density();
  double interpolate_spherical_density(double r);
  double density(double r); // [Msun/kpc^3]
  double integrand_baryon_mass(std::vector<double> rr);
  void compute_baryon_mass();
  double interpolate_baryon_mass(double r);
  double mass(double r); // [Msun]
  double integrand_baryon_potential(std::vector<double> rr);
  void compute_baryon_potential();
  double interpolate_baryon_potential(double r);
  double gravitational_potential(double r); // [(km/s)^2]
  double circular_velocity(double r);       // [km/s]
  double circular_period(double r);         // [Myr]
  double escape_speed(double r);            // [km/s]
  double circular_Ncross(double r);         // Number of disk crossing
  double circular_velocity_baryons_only(double r);
  void plot_circular_velocity();

  // Getter
  double const get_model_dm() const { return m_model_dm; };
  double const get_model_baryons() const { return m_model_baryons; };
};

static double CallBack_integrand_spherical_density(void *pt2object, std::vector<double> z_r);
static double CallBack_integrand_baryon_mass(void *pt2object, std::vector<double> rr);
static double CallBack_integrand_baryon_potential(void *pt2object, std::vector<double> rr);

/* double miyamoto_nagai_potential(double R, double z, double a, double b, double M); */
/* static double integrand1(std::vector<double> z_r_a_b_M_r); */
/* static double integrand2(std::vector<double> z_r_a_b_M_r); */
/* double MN_spherical_potential(double r, double a, double b, double M); */

#endif
