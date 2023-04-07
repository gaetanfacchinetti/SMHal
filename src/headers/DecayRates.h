#ifndef DEF_DECAYRATES
#define DEF_DECAYRATES

#include "DarkSectorModel.h"
#include "ExceptionHandler.h"

class DecayRates
{
public:
  DecayRates(){  v = 246.; T_transQCD = T_QCD_GeV;};
  //DecayRates(){};
  ~DecayRates(){}; // we do not delete _DS since we

  // Decay rates for the scalar mediators
  double DecayRate_sff(DarkSectorModel const& _DS, int prop_s);
  double DecayRate_scc(DarkSectorModel const& _DS, int prop_s);
  double DecayRate_sss(DarkSectorModel const& _DS, int prop_s);
  double DecayRate_spp(DarkSectorModel const& _DS, int prop_s);
  double DecayRate_sXX(DarkSectorModel const& _DS, int prop_s);
  double DecayRate_sgluongluon(DarkSectorModel const &_DS, int prop_s); // This is a one-loop effect
  double DecayRate_sphotonphoton(DarkSectorModel const &_DS, int prop_s); // This is a one-loop effect

  // Decay rates for the pseudo scalar mediators
  double DecayRate_pff(DarkSectorModel const& _DS, int prop_p);
  double DecayRate_pcc(DarkSectorModel const& _DS, int prop_p);
  double DecayRate_psp(DarkSectorModel const& _DS, int prop_p);
  double DecayRate_pgluongluon(DarkSectorModel const &_DS, int prop_p); // This is a one-loop effect
  double DecayRate_pphotonphoton(DarkSectorModel const &_DS, int prop_p); // This is a one-loop effect

  // Decay rates for the vector mediators
  double DecayRate_Xff(DarkSectorModel const& _DS, int prop_X);
  double DecayRate_Xcc(DarkSectorModel const& _DS, int prop_X);

  // Decay rates for the heavy WIMPs
  double DecayRate_csc(DarkSectorModel const& _DS, int prop_chi);
  double DecayRate_cpc(DarkSectorModel const& _DS, int prop_chi);
  double DecayRate_cXc(DarkSectorModel const& _DS, int prop_chi);

  // Total decay rates
  double DecayRate_s(DarkSectorModel const& _DS, int prop_s);
  double DecayRate_p(DarkSectorModel const& _DS, int prop_p);
  double DecayRate_X(DarkSectorModel const& _DS, int prop_X);
  double DecayRate_chi(DarkSectorModel const& _DS, int prop_chi);

  void set_width_phis(DarkSectorModel & _DS);
  void set_width_phip(DarkSectorModel & _DS);
  void set_width_X(DarkSectorModel & _DS);
  void set_width_chi(DarkSectorModel & _DS);

  double pf(double m0, double m1, double m2);

private:

  double v, T_transQCD;
};

#endif