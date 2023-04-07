#ifndef DEF_POWERSPECTRUM
#define DEF_POWERSPECTRUM

#include "CppLib.h"
#include "mymath.h"
#include "Cosmology.h"
#include "TransferFunction.h"
#include "myspline.h"
#include <memory>

enum class PSWindowType
{
  real_space_top_hat,
  gaussian,
  fourier_space_top_hat
};

class PowerSpectrum
{
public:
  PowerSpectrum(){};
  PowerSpectrum(Cosmology n_Cosmo, PSWindowType w = PSWindowType::real_space_top_hat, bool n_normalised_s8 = false, std::shared_ptr<TransferFunction> n_add_TF = nullptr);
  ~PowerSpectrum(){};

  void set_spikes(std::vector<double> n_as_spike, std::vector<double> n_ks_spike, std::vector<double> n_eps_spike);
  void reset_spikes() {As_spike.resize(0); ks_spike.resize(0); eps_spike.resize(0); eps_spike_0.resize(0);};
  double regulator_spikes(double k, double ks_spike, double eps);

  double dimensionless_curvature_power_spectrum(double k);
  double curvature_power_spectrum(double k);
  double matter_power_spectrum(double k, double z);
  double dimensionless_matter_power_spectrum(double k, double z);

  void plot_PowerSpectrum(double z, std::string add = "");
  

  double Window(double kR);
  double Der_Window_vs_kR(double kR);
  double LagrangianRadius_vs_Mass(double Mass);
  double Mass_vs_LagrangianRadius(double Radius);
  double Der_Mass_vs_LagrangianRadius(double Mass);

  /** Variance of the matter power spectrum at redshift z for \f$\Lambda\f$CDM power spectrum
   *  Input  : \f$\mathcal{A}_s\f$ (dimensionless)
   *          \f$k_s~({\rm Mpc}^{-1})\f$
   *          \f$z\f$ (dimensionless)  
   *  Output : \f$\sigma_M(M)\f$ (dimensionless) */
  double SigmaM(double M, double z) { return SigmaR(LagrangianRadius_vs_Mass(M), z); }; // M in Msol
  double SigmaR(double R, double z) { return sqrt(SigmaR2(R, z)); };             // R in Mpc
  double SigmaR2(double R, double z);
  double SigmaM2(double M, double z) { return SigmaR2(LagrangianRadius_vs_Mass(M), z); };
  double fForDichotomie_Mass_vs_Sigma2(double M, double sM2, double z);
  double Mass_vs_Sigma2(double sM2, double z);
  double fToIntFor_Sigma2R(std::vector<double> lnkRz);

  double Der_LnSigmaR_vs_R(double R, double z) { return 0.5 * pow(SigmaR(R, z), -2) * Der_sigmaR2_vs_R(R, z); };
  double Der_LnSigmaM_vs_M(double M, double z) { return Der_LnSigmaR_vs_R(LagrangianRadius_vs_Mass(M), z) / Der_Mass_vs_LagrangianRadius(M); };
  double Der_sigmaR2_vs_R(double R, double z);
  double Der_sigmaM2_vs_M(double M, double z) { return Der_sigmaR2_vs_R(LagrangianRadius_vs_Mass(M), z) / Der_Mass_vs_LagrangianRadius(M); };
  double fToIntFor_Der_sigmaR2_vs_R(std::vector<double> kRz);

  /** Variance of the matter power spectrum at redshift z for spiked power spectrum :  
   *  \f$ \mathcal{P}_{\mathcal{R}} = \mathcal{A}_s * k_s * delta(k- k_s) \f$ 
   *  Input : \f$\mathcal{A}_s\f$ (dimensionless)
   *          \f$k_s~({\rm Mpc}^{-1})\f$
   *          \f$z\f$ (dimensionless)   
   * Output : \f$\sigma_M(M)\f$ (dimensionless) */
  //double SigmaM_spike(double M, double As, double ks, double z) {return sqrt(SigmaM2_spike(M, As, ks, z)); };
  double SigmaR2_narrow_spike_over_window(double As_spike, double ks_spike, double z);
  double SigmaR2_narrow_spike(double R, double As_spike, double ks_spike, double z);


  void Interpolate_SigmaM2_vs_M(double z);
  double Interp_SigmaM2(double M, double z);

  void plot_SigmaM_vs_z(double M);
  void plot_Sigma_vs_RM(double z, bool with_respect_to_M, bool in_units_of_h = false, std::string add = "");
  void plot_Der_LnSigma_vs_LnRM(double z, bool with_respect_to_M = false, bool in_units_of_h = false);
  void plot_fToIntFor_Sigma2R(double z, double R);

  Cosmology get_Cosmology() const { return cosmo; };
  Cosmology *get_Cosmology_ptr()  {return &cosmo;};
  PSWindowType get_window() const { return window_type; };
  TransferFunction_EH98 get_TransferFunction() const {return *TF;};
  
  std::vector<double> get_ks_spike() const {return ks_spike;};
  std::vector<double> get_As_spike() const {return As_spike;};

  double VolumeForm();

  static double CallBack_fToIntFor_Sigma2R(void *pt2Object, std::vector<double> xx) { return ((PowerSpectrum *)pt2Object)->fToIntFor_Sigma2R(xx); };
  static double CallBack_fToIntFor_Der_sigmaR2_vs_R(void *pt2Object, std::vector<double> xx) { return ((PowerSpectrum *)pt2Object)->fToIntFor_Der_sigmaR2_vs_R(xx); };
  static double CallBack_fForDichotomie_Mass_vs_SigmaM2(void *pt2Object, std::vector<double> xx) { return ((PowerSpectrum *)pt2Object)->fForDichotomie_Mass_vs_Sigma2(xx[0], xx[1], xx[2]); };


private:

  Cosmology cosmo;
  std::shared_ptr<TransferFunction_EH98> TF; // Usual transfer function for the matter power spectrum
  std::shared_ptr<TransferFunction> add_TF;  // Additional transfer function to account for specific properties of the dark matter

  PSWindowType window_type;
  double gamma_rs_th, gamma_gaussian, gamma_fs_th;

  double k0;     // pivot scale
  double ns, As; // spectral index and spectral amplitude
  
  std::vector<double> As_spike, ks_spike, eps_spike;
  std::vector<int> eps_spike_0; // spikes for which eps_spike = 0 -> true delta function

  my_spline::spline spline_sigmaM2_vs_log10_M;
  bool is_interpolated_SigmaM2;
  double _z_interp;

  bool normalised_s8;
  double norm_PS;
};

#endif