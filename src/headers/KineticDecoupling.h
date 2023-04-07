#ifndef DEF_KINETICDECOUPLING
#define DEF_KINETICDECOUPLING

#include "CrossSections/CrossSection_cfcf.h"
#include "CrossSections/CrossSection_cgacga.h"
#include "DegreeFreedom.h"
#include "DarkSectorModel.h"
#include "GaussLegendre.h"
#include "ExceptionHandler.h"
#include <exception>

// Exception thrown when a problem occurs in Kinetic decoupling solver
class KineticDecouplingException : public std::exception
{
public:
  KineticDecouplingException() { _msg = (char *)"WARNING << KineticDecoupling exception : "; };
  KineticDecouplingException(char *msg) { _msg = msg; };
  KineticDecouplingException(char *msg, int line, char *func, int type = +1)
  {
    _msg = msg;
    _line = line;
    _func = func;
    _type = type;
  };
  virtual const char *what() const throw() { return _msg; }
  const char *get_func() const { return _func; };
  int get_line() const { return _line; };
  int get_type() const { return _type; };

  // type +1 means that gamma < H and that the kinetic decoupling temperature is lower than the minimal value tested in the dichotomie for its estimation
  // type -1 means that gamma > H and that the kinetic decoupling temperature is higher than the maximal value tested in the dichotomie for its estimation

private:
  char *_msg, *_func;
  int _line, _type;
};

class KineticDecoupling
{

public:
  KineticDecoupling(DarkSectorModel const& _DS, DegreeFreedom *_degFree, int index_chi); // Constructor
  KineticDecoupling(){};
  ~KineticDecoupling()
  {
    for(std::vector<CrossSection *>::iterator it = crossSectionScatt.begin(); it != crossSectionScatt.end(); ++it)
      delete (*it);

    crossSectionScatt.clear(); 
    delete crossSectionScatt_photons;
  }; // Destructor


  /* Total momentum relaxation rate \f$\gamma\f$ 
  *  This function accounts for any possible tree-level scattering with standard model fermions
  *  Input  : \f$T\f$ [GeV] 
  *  Output : \f$\gammaf$ [GeV] */
  double gammaTot(double T);

  /* Momentum relaxation rate \f$\gamma\f$ 
  *  This function accounts for any possible tree-level scattering with standard model fermions
  *  Input  : index_cs (index of the particle the DM scatters onto)
  *           \f$T\f$ [GeV] 
  *  Output : \f$\gammaf$ [GeV] */
  double gamma(int index_cs, double T);

  /* One loop correction of the momentum relaxation rate \f$\gamma\f$ 
  *  This function accounts for the triangle coupling to photons with the (pseudo-)scalar mediators.
  *  Input  : \f$T\f$ [GeV] 
  *  Output : \f$\gammaf$ [GeV] */
  double gamma_1Loop_photons(double T); 

 

  void DerivativeOfTemperature(std::vector<double> var, double w, double &deriv);
  double FToEstimateDecTem(double T);
  double EstimateDecouplingTemperature();
  double SolveTemperatureEI(); // Solver using the euler implicit scheme
  double SolveTemperature();   // Solver using other numerical schemes
  double FToCheckDecTem(double T);


  void plotTemperatureEvolution();
  void plotTemperatureEvolutionEI(std::string name);
  void plotGammaTot();
  void plotGamma_1Loop_photons();
  void plotfToIntForGamma(int index_cs, double T);
  void plotfToEstimateDecTem();


  // CallBack functions
  static double CallBackFToEstimateDecTem(void *pt2Object, std::vector<double> TT) { return ((KineticDecoupling *)pt2Object)->FToEstimateDecTem(TT[0]); };
  static double CallBackFToCheckDecTem(void *pt2Object, double T) { return ((KineticDecoupling *)pt2Object)->FToCheckDecTem(T); };
  static void CallBackDerivativeOfTemperature(void *pt2Object, std::vector<double> var, double w, double &deriv) { return ((KineticDecoupling *)pt2Object)->DerivativeOfTemperature(var, w, deriv); };

  static double CallBackFunctionForImplicitEulerLinearA(void *pt2Object, std::vector<double> variables) { return ((KineticDecoupling *)pt2Object)->FunctionForImplicitEulerLinearA(variables); };
  static double CallBackFunctionForImplicitEulerLinearZ(void *pt2Object, std::vector<double> variables) { return ((KineticDecoupling *)pt2Object)->FunctionForImplicitEulerLinearZ(variables); };

  static double CallBack_fToIntForGamma(void *pt2Object, std::vector<double> xx) { return ((KineticDecoupling *)pt2Object)->fToIntForGamma(xx[0], xx[1], xx[2]); };
  static double CallBack_fToIntForGamma_1Loop_photons(void *pt2Object, std::vector<double> xx) { return ((KineticDecoupling *)pt2Object)->fToIntForGamma_1Loop_photons(xx[0], xx[1]); };

private:
  // Private functions
  void InterpolateGammaTot();
  double fToIntForGamma(double om, double index_cs, double T);
  double fToIntForGamma_1Loop_photons(double om,  double T);

  double FunctionForImplicitEulerLinearA(std::vector<double> var);
  double FunctionForImplicitEulerLinearZ(std::vector<double> var);

  // table af all the scattering cross sections we need
  std::vector<CrossSection *> crossSectionScatt;
  CrossSection *crossSectionScatt_photons;

  Particle chi;
  int n_scatter;
  std::vector<double> mass_scattered_part, spin_scattered_part;

  std::vector<double> ySavedRK6, ySavedEIL, ySavedAML, ySavedBDL, xSaved, yeqSaved;

  double yinf;
  double T_transQCD;

  my_spline::spline spline_derLogHSmooth;
  my_spline::spline spline_gSmooth;
  my_spline::spline spline_hSmooth;
  my_spline::spline spline_gStarHalf;

  my_spline::spline spline_gammaTot; // gives gammaTot vs 
};

#endif
