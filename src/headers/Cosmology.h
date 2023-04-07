#ifndef DEF_COSMOLOGY
#define DEF_COSMOLOGY

#include "CppLib.h"
#include "MyUnits.h"
#include "mymath.h"
#include "myspline.h"

enum class Species
{
  MATTER,
  BARYONS,
  CDM,
  LAMBDA,
  RADIATION
};

enum class CosmoModel
{
  PLANCK,
  PLANCKONLY,
  WMAP5,
  WMAP5ONLY,
  EdSPLANCK,
  UNKNOWN
};

class Cosmology
{

public:
  Cosmology(CosmoModel n_model = CosmoModel::PLANCKONLY);
  //Cosmology(){_Omega_c_h2 = OMEGA_C_h2, _Omega_b_h2 = OMEGA_B_h2, _ns=SPECTRAL_INDEX, _As = SPECTRAL_AMPLITUDE,  _Omega_m_h2 = _Omega_c_h2 + _Omega_b_h2, _Omega_l_h2 = pow(HUBBLE_h,2) - _Omega_m_h2;};
  Cosmology(double n_Omega_c_h2, double n_Omega_b_h2, double n_ns, double n_As, double n_h);
  ~Cosmology(){};

  double get_Omega_m_h2() const { return _Omega_m_h2; };
  double get_Omega_c_h2() const { return _Omega_c_h2; };
  double get_Omega_b_h2() const { return _Omega_b_h2; };
  double get_Omega_l_h2() const { return _Omega_l_h2; };
  double get_sigma8() const
  {
    if (_sigma8 > 0)
      return _sigma8;
    else
    {
      std::cout << "ERROR : No sigma8 defined for this cosmology" << std::endl;
      exit(0);
    }
  };
  double get_h() const { return _h; };
  double get_ns() const { return _ns; };
  double get_As() const { return _As; };

  double cosmic_abundance(double z, Species species);
  double Ez(double z);
  double Hubble_parameter(double z); // [km/s/Mpc]
  double critical_density(double z); // [Msun/kpc^3]
  double cosmological_density(double z, Species species);
  double growth_factor_D1_Carroll(double z);
  double growth_factor_D1_Belloso(double z);
  double ftoIntOnz_growth_factor_exact(double z);
  double growth_factor_exact(double z);
  double derLn_growth_factor_exact(double z);
  double derLn_growth_factor_D1_Carroll(double z);
  double Der_growth_factor_exact(double z);
  double Der_growth_factor_D1_Carroll(double z);
  double Der_growth_factor_D1_Carroll_numerical(double z);

  // Inverse of the growth function
  double Inverse_growth_factor_Caroll(double y);
  double _Inverse_growth_factor_Caroll(double y); // direct computation (avoid using it)
  double f_forBissection_Inverse_growth_factor_Caroll(double z, double y);
  void Interpolate_inverse_growth_factor_Caroll();

  double Zeq_rad_mat() const { return _zeq_rm; };
  double Zeq_mat_lam() const { return pow(_Omega_l_h2 / _Omega_m_h2, 1. / 3.) - 1; };
  
  // k at equivalence in Mpc^{-1}
  double keq_rad_mat() const { return _keq_Mpcm1; };

  // k in Mpc^{-1}
  double redshift_horizon_entry(double k);
  double fToSolve_redshift_horizon_entry(double z, double k);

  double CosmicTime_in_yr(double z);
  double CosmicTime_times_Hubble_parameter(double z);
  double fToIntForCosmicTime(double z);

  void plot_cosmic_abundance(std::string cosmo);
  void plot_Der_growth_factor();
  void plot_Hubble_Parameter(std::string cosmo);
  void plot_growth_factor();

  static double CallBack_fToSolve_redshift_horizon_entry(void *pt2Object, std::vector<double> zz) { return ((Cosmology *)pt2Object)->fToSolve_redshift_horizon_entry(zz[0], zz[1]); };
  static double CallBack_growth_factor_D1_Carroll(void *pt2Object, std::vector<double> zz) { return ((Cosmology *)pt2Object)->growth_factor_D1_Carroll(zz[0]); };
  static double CallBack_fToIntForCosmicTime(void *pt2Object, std::vector<double> zz) { return ((Cosmology *)pt2Object)->fToIntForCosmicTime(zz[0]); };
  static double CallBack_ftoIntOnz_growth_factor_exact(void *pt2Object, std::vector<double> zz) { return ((Cosmology *)pt2Object)->ftoIntOnz_growth_factor_exact(zz[0]); };
  static double CallBack_f_forBissection_Inverse_growth_factor_Caroll(void *pt2Object, std::vector<double> xx) { return ((Cosmology *)pt2Object)->f_forBissection_Inverse_growth_factor_Caroll(xx[0], xx[1]); };

  // Overcharged comparison operators
  friend bool operator== (const Cosmology& c1, const Cosmology& c2);
  friend bool operator!= (const Cosmology& c1, const Cosmology& c2);

private:
  double _Omega_c_h2, _Omega_b_h2, _ns, _As, _h;
  double _Omega_m_h2, _Omega_l_h2, _Omega_r_h2;
  double _sigma8; ///< Value of sigma8 At R = 8/h Mpc in the real-space top hat
  CosmoModel model;

  double _keq_Mpcm1, _zeq_rm, _aeq_rm; ///< Constants at equivalence Matter/Radiation

  bool is_inverse_growth_Caroll_interpolated;
  my_spline::spline spline_inverse_log10z_Caroll_log10D1;

};


#endif