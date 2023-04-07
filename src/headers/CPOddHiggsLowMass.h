#ifndef DEF_CPODDHIGGSLOWMASS
#define DEF_CPODDHIGGSLOWMASS

#include "CppLib.h"
#include "StandardModel.h"
#include "mymath.h"

class CPOddHiggsLowMass
{
public:
  CPOddHiggsLowMass(std::string const &n_name_file);
  ~CPOddHiggsLowMass(){};

  void Read();
  void Write(std::string const &name_file_out);
  void Diagonalise(int i);

  dcomp FF(double x);
  dcomp Cg(int i);
  dcomp sCgam(int i, double s);
  dcomp sdCgam(int i, double s);
  double ach0A(int i, double T_over_TQCD);

  dcomp sFg(int i, double s, double T_over_TQCD);

  // For numerical tests
  dcomp sFg_test(int i, double s, double T_over_TQCD);

  // Initialise the vector version of : 4 mesons and A1-3 mesons amplitudes
  void AmplitudeInitialise();

  // 4 mesons and A1-3 mesons amplitudes
  double AmpA1f(int i, int mes1, int mes2, int mes3);
  double AmpPi3f(int i, int mes1, int mes2, int mes3);
  double AmpEtaf(int i, int mes1, int mes2, int mes3);
  double AmpEtapf(int i, int mes1, int mes2, int mes3);

  // Amplitude for the decay in 3 mesons of the physical A1 tilde
  double AmpAf(int i, int mes1, int mes2, int mes3);

  // Decay Rates
  double GamaA(int i, int mes1, int mes2, int mes3); // A1->mes1+mes2+mes3
  double fToIntForGamaA(std::vector<double> ismes1mes2mes3);

  double GamaAch(int i, double T_over_TQCD);  // Total decay of A1 into three mesons with the chiral lagrangian
  double GamaAsp(int i);                      // Total decay of A1 into three mesons with the spectator lagrangian
  double GamaLe(int i, double T_over_TQCD);   // Decay of A1 into e+ and e-
  double GamaLmu(int i, double T_over_TQCD);  // Decay of A1 into mu+ and mu-
  double GamaL(int i, double T_over_TQCD);    // Total decay of A1 into e+,e- and mu+,mu-
  double Gamach0(int i, double T_over_TQCD);  // Decay of A1 into dark matter
  double GamaAt1S(int i, double T_over_TQCD); // Decay of A1 into two photons
  double GammaTAS(int i, double T_over_TQCD); // Total decay of A1

  double Sfactor(int mes1, int mes2, int mes3);
  void Order(int &i, int &j, int &k);

  void plot_MixingParameters();
  void plot_GammaA();
  void plot_GamaA(int mes1, int mes2, int mes3);
  void plot_Amp();
  void plot_sFg(int i);
  void plot_sCgam(int i);

  // CallBack functions
  static double CallBack_fToIntForGamaA(void *pt2Object, std::vector<double> xx) { return ((CPOddHiggsLowMass *)pt2Object)->fToIntForGamaA(xx); };

private:
  std::string name_file;

  // Number of points in the input file
  int Npoints;

  std::vector<double> lamb, kappa, Tbeta, mH1, mA1, mch0, mch1, mch2, Cotb, Sbeta, Cbeta;
  std::vector<double> GA1N1N1, U11, U12, U21, U22, V11, V12, V21, V22, P11, P12;
  double V, del, del1, del2, del3, del4, del5, del6, del7, del8, del9, Delta, fpi, mpi, mpip, mpim, mpi8, mpi9, meta, metp, te, Nc, Qu, Qd, Qe, Qch, g;
  double vv, mK, mKz;
  std::vector<double> OAA, OA3, OA8, OA9, OAeta, OAetp;
  std::vector<double> MM; // Meson mass vector
  std::vector<double> YllA1e, YllA1mu;

  // Wess-Zumino-Witten couplings between CP Odd mesons and photons
  double Cpi3, Ceta, Cetp;
  double Yepi3, Ymueta;

  // Mesons index
  int _pi3, _pip, _pim, _eta, _etp, _Kp, _Km, _Kz;

  // Amplitudes for decay into three mesons
  std::vector<std::vector<std::vector<std::vector<double>>>> AmpA1, AmpPi3, AmpEta, AmpEtap;

};

#endif
