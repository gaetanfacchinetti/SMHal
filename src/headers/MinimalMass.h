#ifndef DEF_MINIMALMASS
#define DEF_MINIMALMASS

#include "CppLib.h"
#include "DegreeFreedom.h"
#include "Cosmology.h"

class MinimalMass
{

public:
  MinimalMass(DegreeFreedom *_degFree, Cosmology Cosmo);
  MinimalMass(){};
  ~MinimalMass(){};

  double FreeStreamingLength(double Tkd, double mDM);
  double FreeStreamingLength(double Tkd, double mDM, double a){return -1;};
  double FToIntToComputeKfs(std::vector<double> variables);
  std::vector<double> MassesComputation(double Tkd, double mDM);
  //void set_Tkd(double const &TkdNew);
  double VelocityAtKineticDecoupling(double Tkd, double mDM); 
  double ComovingFreeStreamingScale(double Tkd, double mDM, double z); // GeV
  double ComovingCollisionalDampingScale(double Tkd, double mDM); // GeV
  double ComovingFreeStreamingLength(double Tkd, double mDM, double z); // GeV^{-1}
  
  void PlotMinimalMass(double mDM);
  void PlotFreeStreamingLength_vs_z(double Tkd, double mDM);
  void PlotFreeStreamingLength_vs_Tkd(double z, double mDM);
  // CallBack functions
  static double CallBackFToIntToComputeKfs(void *pt2Object, std::vector<double> variables) { return ((MinimalMass *)pt2Object)->FToIntToComputeKfs(variables); }

private:


  double Tstar, Teq_GeV;
  double a0OutOfaeq, zeq, rhoEq, rhoM0, Heq, aeq;
  //double akdOutOfaeq;

  my_spline::spline spline_derLogHSmooth;
  my_spline::spline spline_gSmooth;
  my_spline::spline spline_hSmooth;
  my_spline::spline spline_gStarHalf;
};

#endif