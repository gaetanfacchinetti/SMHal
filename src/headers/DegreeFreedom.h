#ifndef DEF_DEGREEFREEDOM
#define DEF_DEGREEFREEDOM


#include "StandardModel.h"
#include "DarkSectorModel.h"
#include "myspline.h"
#include "Derivateur.h"

enum class DegFreeModel
{
  LATTICE2016,
  SIMPLE
};

class DegreeFreedom
{

public:
  DegreeFreedom(double T_transQCD, DarkSectorModel *_DS = NULL);
  ~DegreeFreedom(){};

  void ComputeGAndH();
  void ComputeGAndHSmooth();
  void ComputeGStarHalf();

  void InitialiseLattice();

  int Interpolate();
  std::vector<double> reverseVector(std::vector<double> const& vec);


  // Getter
  std::vector<double> get_g() const {return g_vec;};
  std::vector<double> get_h() const {return h_vec;};
  std::vector<double> get_gSmooth() const {return gSmooth_vec;};
  std::vector<double> get_hSmooth() const {return hSmooth_vec;};
  std::vector<double> get_gStarHald() const {return gStarHalf_vec;};
  std::vector<double> get_logT()      const {return logT_vec;};

  double get_T_transQCD() const {return T_transQCD;};

  // BE CAREFULL these functions take log10(T) in input and NOT T
  double get_derLogHSmoothInterp(double const& log10T) const {return spline_derLogHSmooth(log10T);};
  double get_gSmoothInterp(double const& log10T) const {return spline_gSmooth(log10T);};
  double get_hSmoothInterp(double const& log10T) const {return spline_hSmooth(log10T);};
  double get_gStarHalfInterp(double const& log10T) const {return spline_gStarHalf(log10T);};

   // BE CAREFULL the returned functions take log10(T) in input and NOT T
  my_spline::spline get_gSmoothInterp() const {return spline_gSmooth;};
  my_spline::spline get_hSmoothInterp() const {return spline_hSmooth;};
  my_spline::spline get_gStarHalfInterp() const {return spline_gStarHalf;};
  my_spline::spline get_derLogHSmoothInterp() const {return spline_derLogHSmooth;};

  // Hubble constant function of T
  //double Hubble_T(double T) {return 2*pow(PI, 1.5)*sqrt(spline_gSmooth(log10(T))) *pow(T,2) / (sqrt(45)*M_PLANCK);};

  /// Function that plot the evolution of G and H
  void plot_GAndH();
  //void QCDPhaseTransitionForGAndHv1();

private:
  Derivateur *Der;
  std::vector<Particle> all_Part;

  double T_transQCD;
  double numSpecies;

  std::vector<double> mass, stat;
  std::vector<int> degFree;
  std::vector<bool> isPresentAfterQCDPT, isPresentBeforeQCDPT;
  std::vector<Particle> SM_neutrinos;
  std::vector<int> neutrinos_index; 
  int n_neutrinos;

  std::vector<double> g_vec, h_vec, logG_vec, logH_vec, logT_vec;
  std::vector<double> gSmooth_vec, hSmooth_vec, logHSmooth_vec, gStarHalf_vec;
  std::vector<double> derLogHSmooth_vec;
  std::vector<std::string> name;

  std::ofstream out_GAndH, out_GAndHSmooth, out_GStarHalf;

  int nBessel;         // Number of terms in sums
  int numPoints;       // Number of discretisation points
  double Tmax;         // Maximal temperature (en GeV)
  double Tmin;         // Minimal temperature (en GeV)

  double logTmax,logTmin, deltaLogT;

  bool isComputeGStarHalfCalled, isInterpolateCalled;

  // Spline interpolation for the simple QCD phase transition
  my_spline::spline spline_derLogHSmooth;
  my_spline::spline spline_gSmooth;
  my_spline::spline spline_hSmooth;
  my_spline::spline spline_gStarHalf;

  // Spline interpolation from LATTICE arXiv:1606.07494
  my_spline::spline spline_gLattice;
  my_spline::spline spline_hLattice;
  my_spline::spline spline_gStarHalfLattice; 



};


#endif
