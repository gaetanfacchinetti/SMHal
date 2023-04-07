#ifndef DEF_DARK_HALO
#define DEF_DARK_HALO

#include "CppLib.h"
#include "mymath.h"
#include "myspline.h"
#include "MyUnits.h"
#include "Cosmology.h"

enum class DensityProfile
{
    NFW,
    MOORE
};

class DarkHalo
{
 protected:
  
  int m_profile;
  double m_alpha;
  double m_beta;
  double m_gamma;
  double m_rhos;
  double m_rs;
  Cosmology m_cosmo;

 public:

  DarkHalo(Cosmology cosmo);
  DarkHalo(Cosmology cosmo, DensityProfile profile, double rhos = 0, double rs = 0);
  DarkHalo(Cosmology cosmo, int profile, double alpha, double beta, double gamma, double rhos=0, double rs=0);
  DarkHalo(Cosmology cosmo, double rhos, double rs);

  void Initialise_from_virial_parameters(int profile, double alpha, double beta, double gamma, double mvir, double cvir=-1, double delta=200, bool is_critical=true);
  void Initialise_from_virial_parameters(double mvir, double cvir=-1, double delta = 200, bool is_critical= true);
  double c_bar(double m200);

  int get_profile() const;
  double get_alpha() const;
  double get_beta() const;
  double get_gamma() const;
  double get_rhos() const;
  double get_rs() const;
  Cosmology get_cosmo() const {return m_cosmo;};

  void update_rhos(double rhos);
  void update_rs(double rs);
  void update_profile(int profile);
  void update_alpha(double alpha);
  void update_beta(double beta);
  void update_gamma(double gamma);

  double virial_radius(double delta, bool is_critical); // [kpc]
  double virial_mass(double delta, bool is_critical); // [Msun]
  double virial_concentration(double delta, bool is_critical);
  double scale_density(double c_vir, double delta, bool is_critical); // [Msun/kpc^3]
  double scale_radius(double c_vir, double m_vir, double delta, bool is_critical); // [kpc]
  double scale_radius_bis(double rhos, double m_vir, double delta, bool is_critical);
  
  double density_profile(double x);
  double dm_density(double r); // [Msun/kpc^3]
  double integrand_mass_profile(std::vector<double> xx);
  void compute_mass_profile();
  double interpolate_mass_profile(double x);
  double mass_profile(double x);
  double dm_mass(double r); // [Msun]
  double integrand_grav_potential(std::vector<double> xx);
  void compute_grav_potential_profile();
  double interpolate_grav_potential_profile(double x);
  double grav_potential_profile(double x);
  double dm_grav_potential(double r); // [(km/s)^2]
  double integrand_potential_energy(std::vector<double> xx);
  void compute_potential_energy_profile();
  double interpolate_potential_energy_profile(double x);
  double potential_energy_profile(double x); // negative
  double dm_potential_energy(double r); // [Msun*(km/s)^2] negative
  double integrand_binding_energy(std::vector<double> xx);
  void compute_binding_energy_profile();
  double interpolate_binding_energy_profile(double x);
  double binding_energy_profile(double x); // negative
  double dm_binding_energy(double r); // [Msun*(km/s)^2] positive
  double integrand_luminosity_profile(std::vector<double> xx);
  void compute_luminosity_profile();
  double interpolate_luminosity_profile(double x);
  double luminosity_profile(double x);
  double dm_luminosity(double r); // [Msun^2/kpc^3]

  double rm2_rs();

  double velocity_dispersion_profile(double x, double x0);
  double velocity_dispersion(double r, double r0);
  double orbital_frequency_profile(double x, double x0);
  double orbital_frequency(double r, double r0);
  double circular_velocity_DM_only(double r); // r [kpc] -> [km/s] 

  void plot_velocity_dispersion_and_grav_pot(double x0);

  // Same function than above but this time the dispersion in computed
  // using a direct integration instead of polyLog functions
  double fToInt_velocity_dispersion_integral(double r);
  double velocity_dispersion_integral(double r, double r0);
  double orbital_frequency_integral(double r, double r0);

  private:
    double f_for_velocity_disperison(double y);

};

static double CallBack_integrand_mass_profile(void *pt2object, std::vector<double> xx);
static double CallBack_integrand_grav_potential(void *pt2object, std::vector<double> xx);
static double CallBack_integrand_potential_energy(void *pt2object, std::vector<double> xx);
static double CallBack_integrand_binding_energy(void *pt2object, std::vector<double> xx);
static double CallBack_integrand_luminosity_profile(void *pt2object, std::vector<double> xx);
static double CallBack_fToInt_velocity_dispersion_integral(void *pt2object, std::vector<double> xx){ return ((DarkHalo *)pt2object)->fToInt_velocity_dispersion_integral(xx[0]); };

#endif
