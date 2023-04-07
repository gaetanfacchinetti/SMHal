#ifndef DEF_FSLModel
#define DEF_FSLModel

#include "CppLib.h"
#include "mymath.h"
#include "MyUnits.h"
#include "Galaxy.h"
#include "DarkHalo.h"
#include "myspline.h"

// BE CAREFULL : It is impossible to define a c200_max_point with the new definition
// of a point subhalo

class FSLModel
{
public:
  FSLModel(double n_m200_min, double n_epst, DarkHalo n_hostHalo);
  ~FSLModel(){};

  double tidal_radius_over_scale_radius(double R, double c200);
  void plot_tidal_radius_over_scale_radius();
  void plot_tidal_radius_over_scale_radius_linear_concentration();

  void ReadFormTable_rt_over_rs_DMO();
  double Interpolation_rt_over_rs_DMO(double R, double c200);
  void plot_rt_over_rs();

  double c200_min_DMO(double R, double eps_t);
  double funcToSolvec200min_DMO(std::vector<double> Rc200eps);


  // Function to compute c_bar in terms of the physical mass mt
  double c_bar_vs_mphys(double mt, double R);
  double fForDichotomie_c_bar_vs_mphys(double c200, double mt, double R);
  void plot_c200_vs_mphys();

  //double c200_max_point(double s_i, double cosbcosl, double eta_i, double theta_d, double m200);
  //double funcToSolvec200max_point(std::vector<double> c200sRetathetam200);
  //void plot_funcToSolvec200max_point(double s_i, double cosbcosl, double eta_i, double theta_d, double m200);
  //void plot_c200_max_point(double s_i, double psi, double eta_i, double theta_d, double m200_min, double m200_max);

  double rt_over_rs_DMO(double R, double c200);
  double funcToSolve_tidal_radius_DMO(std::vector<double> xtRc200);
  void plot_dimensionless_tidal_radius_DMO();

  double fToIntOnm200AverageRhoSub(std::vector<double> yy);
  double fToIntOnc200AverageRhoSub(std::vector<double> yy);
  double UnNormalizedAverageRhoSub(double R);
  int InitialiseAverageRhoSub();

  double averageRhoSub(double R)
  {
    //return NsubUnevolved()*UnNormalizedAverageRhoSub(R);
    if (is_rhoSub_initialized == true)
      return spline_AverageRhoSub(R);
    else
    {
      is_rhoSub_initialized = true;
      InitialiseAverageRhoSub();
      return spline_AverageRhoSub(R);
    }
  };

  double get_Kt_Tot()
  {
    if (is_Kt_Tot_initialized == true)
      return _Kt_Tot;
    else
    {
      _Kt_Tot = Normalisation_Kt_Tot();
      is_Kt_Tot_initialized = true;
      return _Kt_Tot;
    }
  }

  double get_Nsub()
  {
    if (is_N_sub_initialized == true)
      return _N_sub;
    else
    {
      //std::cout << NsubUnevolved() << " " << get_Kt_Tot() << std::endl;
      _N_sub = NsubUnevolved()*get_Kt_Tot();
      is_N_sub_initialized = true;
      return _N_sub;
    }
  }

  // Getter
  double const get_alpha_m() const { return alpha_m; };
  double const get_epst() const { return epst; };
  double const get_m200_min() const { return m200_min; };
  double const get_m200_max() const { return m200_max; };
  double const get_c200_max() const { return c200_max; };
  double const get_Km() const { return Km; };
  double const get_R_max() const { return R_max; };
  DarkHalo const get_subModel() const { return subModel; };
  DarkHalo const get_hostHalo() const {return hostHalo;};
  //Galaxy const get_galaxy() const { return gal; };
  Cosmology get_cosmo() const {return cosmo;};

  double const get_vec_c200(int i) const { return vec_c200[i]; };
  double const get_vec_R(int i) const { return vec_R[i]; };
  double const get_vec_rt_over_rs(int i, int j) const { return vec_rt_over_rs[i][j]; };

  //
  //
  double MassFunction(double m200_M200, double M200);
  double TotalNumberOfSubhalosUnevolved();

  double ProbDistributionSpace(double R);
  double ProbDistributionMass(double m200);
  double ProbDistributionConcentration(double c200, double m200);
  double c_bar(double m200);

  double fToIntOnc200_MassFraction(std::vector<double> yy);
  double fToIntOnm200_MassFraction(std::vector<double> yy);
  double fToIntOnR_MassFraction(std::vector<double> yy);
  double MassFraction();
  double NsubUnevolved()
  {
    if (is_NsubUnevolved_initialised == true)
    {
      //std::cout << "There" << std::endl;
      return _N_sub_un;
    }
    else
    {
      is_NsubUnevolved_initialised = true;
      _N_sub_un = TotalNumberOfSubhalosUnevolved();
      //std::cout << "here we have : " << _N_sub_un << std::endl;
      return _N_sub_un;
    }
  };

  double fToIntOnc200_EvolvedMassFunction(double c200, double mt, double R, double c200_m);
  double EvolvedMassFunction_Rfixed(double mt, double R, double c200_m);
  double fToIntOnR_EvolvedMassFunction(double mt, double R);
  double EvolvedMassFunction(double mt);
  void plot_EvolveMassFunction(int version);
  double Nsim_Ktsim() { return 1; };

  double fToIntOnm200_Norm_Tot(std::vector<double> varNorm);
  double fToIntOnc200_Norm_Tot(std::vector<double> varNorm);
  double fToIntOnR_Norm_Tot(std::vector<double> varNorm);
  double Normalisation_Kt_Tot();

  std::vector<double> Nsub_calibrated();

  double NumberDensity_subhalos(double R);
  void plot_NumberDensity_subhalos(int id);

  //
  //

  
  static double CallBack_funcToSolvec200min_DMO(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->funcToSolvec200min_DMO(xx); };
  static double CallBack_funcToSolve_tidal_radius_DMO(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->funcToSolve_tidal_radius_DMO(xx); };

  static double CallBack_fToIntOnc200_EvolvedMassFunction(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnc200_EvolvedMassFunction(xx[0], xx[1], xx[2], xx[3]); };
  static double CallBack_fToIntOnR_EvolvedMassFunction(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnR_EvolvedMassFunction(xx[0], xx[1]); };

  static double CallBack_fToIntOnm200_Norm_Tot(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnm200_Norm_Tot(xx); };
  static double CallBack_fToIntOnc200_Norm_Tot(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnc200_Norm_Tot(xx); };
  static double CallBack_fToIntOnR_Norm_Tot(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnR_Norm_Tot(xx); };

  static double CallBack_MassFunction(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->MassFunction(xx[0], xx[1]); };

  static double CallBack_fToIntOnm200_MassFraction(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnm200_MassFraction(xx); };
  static double CallBack_fToIntOnc200_MassFraction(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnc200_MassFraction(xx); };
  static double CallBack_fToIntOnR_MassFraction(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnR_MassFraction(xx); };

  static double CallBack_fToIntOnc200AverageRhoSub(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnc200AverageRhoSub(xx); };
  static double CallBack_fToIntOnm200AverageRhoSub(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fToIntOnm200AverageRhoSub(xx); };
  static double CallBack_fForDichotomie_c_bar_vs_mphys(void *pt2Object, std::vector<double> xx) { return ((FSLModel *)pt2Object)->fForDichotomie_c_bar_vs_mphys(xx[0], xx[1], xx[2]); };

  //
  //
  //

private:


  std::vector<std::vector<double>> vec_rt_over_rs, vec_rt_over_rs_DMO, vec_c200_min;
  std::vector<double> vec_R, vec_R_DMO;
  std::vector<double> vec_c200, vec_c200_DMO;
  std::vector<double> vec_R_bis, vec_epst;

  double distance_SGC;
  double Km, alpha_m, m200_min, m200_max, c200_max, R_max;
  double epst;

  Cosmology cosmo;
  //Galaxy gal;
  DarkHalo subModel, hostHalo;
  double host_virialmass;

  std::string input_filename;

  // Average value of rho_sub
  my_spline::spline spline_AverageRhoSub;
  double _N_sub, _N_sub_un, _Kt_Tot;

  bool is_rhoSub_initialized;
  bool is_Kt_Tot_initialized;
  bool is_N_sub_initialized;

  bool is_NsubUnevolved_initialised;

  bool are_disk_effects; // variable to say if disk effects are taken into account or not
};

#endif