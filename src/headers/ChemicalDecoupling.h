#ifndef DEF_CHEMICALDECOUPLING
#define DEF_CHEMICALDECOUPLING

#include "CrossSections/CrossSection_ccff.h"
#include "CrossSections/CrossSection_ccss.h"
#include "CrossSections/CrossSection_ccsp.h"
#include "CrossSections/CrossSection_ccpp.h"
#include "DegreeFreedom.h"
#include "DarkSectorModel.h"
#include "gsl/gsl_integration.h"
#include "ExceptionHandler.h"
#include <exception>

enum class CS_TYPE
{
  FULL_ONLY,
  S_WAVE_ONLY,
  P_WAVE_ONLY,
  SIMPLIFIED_ONLY,
  FULL_SIMPLIFIED
};
// Full simplfied means that we use the expension on

// Exception thrown when freeze-in might occur
class FreezeInException : public std::exception
{
public:
  FreezeInException() { _msg = (char *)"WARNING << FreezeIn exception : "; };
  FreezeInException(char *msg) { _msg = msg; };
  FreezeInException(char *msg, int line, char *func)
  {
    _msg = msg;
    _line = line;
    _func = func;
  };
  virtual const char *what() const throw() { return _msg; }
  const char *get_func() const { return _func; };
  int get_line() const { return _line; };

private:
  char *_msg, *_func;
  int _line;
};

//
//
struct gsl_f_params;

class ChemicalDecoupling
{
public:
  ChemicalDecoupling(DarkSectorModel const &_DS, DegreeFreedom *_degreeFreedom, CS_TYPE cs_type = CS_TYPE::FULL_ONLY, double n_simpleSigmaV = -1);
  ~ChemicalDecoupling()
  {
    for (int i = 0; i < crossSection.size(); i++)
      for (int j = 0; j < crossSection[i].size(); j++)
        for (int k = 0; k < crossSection[i][j].size(); k++)
        {
          delete crossSection[i][j][k];
          crossSection[i][j][k] = NULL;
        }

    // Do not forget to delete the workspace ! Otherwise terrible memory leaks problems
    gsl_integration_workspace_free(workspace);
  };

  void set_saveYToPrint(bool const &NewSaveYToPrint)
  {
    if (saveYToPrint != NewSaveYToPrint)
    {

      saveYToPrint = NewSaveYToPrint;
      std::cout << "|| INFO    : (AbundanceComp.) SaveYToPrint has been changed to        : " << saveYToPrint << std::endl;
    }
  };

  double DeterminationXfMethod1();
  double DeterminationXfMethod2();
  double fToSolveXfModel1Method1(double const &x);
  double fToSolveXfModel1Method2(double const &x);
  double Dichotomie(double (ChemicalDecoupling::*func)(double const &), double const &xMin, double const &xMax, double const &prec);
  double static FunctionTest(std::vector<double> x);

  void OutputGammaAndH();
  void OutputCrossSection();

  /** 
   * Thermally averaed annihilation cross-section   
   * not taking into account the quarks once they   
   * are not only present as hadrons (when \f$T < T_QCD\f$)
   */
  double SigmaV(double T, bool with_quarks_after_TQCD = false);
  /// Function to integrate in order to evaluate the total thermally averaged annihilation cross-section including the quarks
  double fToIntegrateForSigmaV_all(std::vector<double> variables);
  /// Function to integrate in order to evaluate the thermally averaged annihilation cross-section without quarks after \f$T<T_QCD\f$
  double fToIntegrateForSigmaV_noquarks(std::vector<double> variables);
  // Evaluate the s-wave term of the cross-section
  void Initialise_s_wave_term();


  static double DerivativeOfComovingDensity(std::vector<double> variables, double Y);
  void DerivativeOfComovingDensityNonStiff(std::vector<double> variables, double W, double &deriv);
  std::vector<double> SolveFreezeOut();

  double ComputeAbundanceToday(std::vector<double> const &xYInitForAbundanceToday);
  double FunctionToIntForAbundanceToday(std::vector<double> variables);

  // Zone of tested functions
  void PlotYFO();
  void PlotYFI();
  void plotFuncToIntegrateSigmaV(double T);
  void plotSigmaV();
  void SolveFreezeIn();
  double ComputeYFreezeInApprox(double x);
  double FunctionToIntFreezeInApprox(std::vector<double> var);
  double FunctionForImplicitEulerA(std::vector<double> variables);
  double FunctionForImplicitEulerB(std::vector<double> variables);
  double FunctionForImplicitEulerLnA(std::vector<double> variables);
  double FunctionForImplicitEulerZ(std::vector<double> variables);

  //
  //
  // CallBack functions

  static double CallBackfToSolveXfModel1Method1(void *pt2Object, double x) { return ((ChemicalDecoupling *)pt2Object)->fToSolveXfModel1Method1(x); };
  static double CallBackfToSolveXfModel1Method2(void *pt2Object, std::vector<double> xx) { return ((ChemicalDecoupling *)pt2Object)->fToSolveXfModel1Method2(xx[0]); };
  static double CallBackfToIntegrateForSigmaV_all(void *pt2Object, std::vector<double> variables) { return ((ChemicalDecoupling *)pt2Object)->fToIntegrateForSigmaV_all(variables); };
  static double CallBackfToIntegrateForSigmaV_noquarks(void *pt2Object, std::vector<double> variables) { return ((ChemicalDecoupling *)pt2Object)->fToIntegrateForSigmaV_noquarks(variables); };
  static double CallBackFunctionToIntForAbundanceToday(void *pt2Object, std::vector<double> variables) { return ((ChemicalDecoupling *)pt2Object)->FunctionToIntForAbundanceToday(variables); };
  static void CallBackDerivativeOfAbundanceNonStiff(void *pt2Object, std::vector<double> variables, double W, double &deriv) { return ((ChemicalDecoupling *)pt2Object)->DerivativeOfComovingDensityNonStiff(variables, W, deriv); };
  static double CallBackFunctionForImplicitEulerA(void *pt2Object, std::vector<double> variables) { return ((ChemicalDecoupling *)pt2Object)->FunctionForImplicitEulerA(variables); };
  static double CallBackFunctionForImplicitEulerB(void *pt2Object, std::vector<double> variables) { return ((ChemicalDecoupling *)pt2Object)->FunctionForImplicitEulerB(variables); };
  static double CallBackFunctionForImplicitEulerLnA(void *pt2Object, std::vector<double> variables) { return ((ChemicalDecoupling *)pt2Object)->FunctionForImplicitEulerLnA(variables); };
  static double CallBackFunctionForImplicitEulerZ(void *pt2Object, std::vector<double> variables) { return ((ChemicalDecoupling *)pt2Object)->FunctionForImplicitEulerZ(variables); };
  /*




  */

private:
  void InterpolateCrossSectionP1cm2();
  void InterpolateSigmaV();

  gsl_integration_workspace *workspace;
  int n_points_int;

  // WIMPs particles
  std::vector<Particle> chi;
  std::vector<double> m_chi, m_prop, w_prop;
  std::vector<int> g_chi;
  bool is_CoAnnihilation;

  /** 
   * Vector of all cross sections to use
   * For each couple of incoming particles [k][l] one assign a cross-section  
   * that correspond to [m] the possible output
   */
  std::vector<std::vector<std::vector<CrossSection *> > > crossSection;
  double n_crossSections, n_crossSections_ccff;
  double simpleSigmaV;

  // Variable that set the method of the resolution
  CS_TYPE cs_type;

  int degLibDM;
  double couplingScalarDM, couplingPScalarDM, couplingScalarFermions, couplingPScalarFermions;
  double simpleCrossSection;

  double T_transQCD;

  int modelApproximateCrossSection;
  int useRealCrossSectionWP; // Use real cross section when possible

  std::vector<double> coeffSigmaV;
  std::vector<double> xPrint, YIEPrint, YRK6Print, YeqPrint, YFIPrint;

  double GeV2toCm3Sm1;
  bool saveYToPrint;

  my_spline::spline spline_derLogHSmooth;
  my_spline::spline spline_gSmooth;
  my_spline::spline spline_hSmooth;
  my_spline::spline spline_gStarHalf;

  std::vector<my_spline::spline> spline_crossSection;
  my_spline::spline spline_sigmaV_all;
  my_spline::spline spline_sigmaV_noquarks;
};

// ===================================================
// Parameters for integration with gsl QAGS function
// ===================================================

struct gsl_f_params
{
  double _x;
  ChemicalDecoupling *pt_MyClass;
};

static double gslClassWrapperLn_fToIntegrateForSigmaV_all(double lnz, void *pp)
{
  gsl_f_params *p = (gsl_f_params *)pp;
  std::vector<double> zx;
  zx.resize(2);
  zx[0] = exp(lnz);
  zx[1] = p->_x;
  //for(int i = 1; i<p->isize; i++) xx[i] = p->params[i-1];
  return zx[0] * (p->pt_MyClass->fToIntegrateForSigmaV_all(zx));
}

static double gslClassWrapperLn_fToIntegrateForSigmaV_noquarks(double lnz, void *pp)
{
  gsl_f_params *p = (gsl_f_params *)pp;
  std::vector<double> zx;
  zx.resize(2);
  zx[0] = exp(lnz);
  zx[1] = p->_x;
  //for(int i = 1; i<p->isize; i++) xx[i] = p->params[i-1];
  return zx[0] * (p->pt_MyClass->fToIntegrateForSigmaV_noquarks(zx));
}

/*
static double gslClassWrapperLnLn_FunctionToIntegrateForSigmaV(double lnlnz, void *pp)
{
  gsl_f_params *p = (gsl_f_params *)pp;
  std::vector<double> zx;
  zx.resize(2);
  zx[0] = exp(exp(lnlnz));
  zx[1] = p->_x;
  //for(int i = 1; i<p->isize; i++) xx[i] = p->params[i-1];
  return exp(lnlnz) * zx[0] * (p->pt_MyClass->FunctionToIntegrateForSigmaV(zx));
}*/

/*
static double gslClassWrapper_FunctionToIntegrateForSigmaV(double z, void *pp)
{
  gsl_f_params *p = (gsl_f_params *)pp;
  std::vector<double> zx;
  zx.resize(2);
  zx[0] = z;
  zx[1] = p->_x;
  //for(int i = 1; i<p->isize; i++) xx[i] = p->params[i-1];
  return (p->pt_MyClass->FunctionToIntegrateForSigmaV(zx));
}*/

#endif

//double FunctionToIntegrateForSigmaV(double z, void* parameters){ gsl_p->pt_MyClass = this; = return 0;};