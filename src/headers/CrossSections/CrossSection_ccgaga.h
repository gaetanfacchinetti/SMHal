#ifndef DEF_CROSSSECTION_CCGAGA
#define DEF_CROSSSECTION_CCGAGA

#include "CrossSection.h"
#include "../DarkSectorModel.h"
#include "../CPOddHiggsLowMass.h"

class CrossSection_ccgaga : public CrossSection
{

public:
  CrossSection_ccgaga(DarkSectorModel & DSMod, Particle const& n_chii, Particle const& n_chij, CPOddHiggsLowMass *n_CPOHLM, int i_CPOddHiggsModel = 0, double T_over_TQCD=0.1);
  virtual ~CrossSection_ccgaga(){};

  dcomp Q_ssphipphip(int n, double s);

  double sFg(double s);

  dcomp Q(ParticlePair const &partPair, int n, double s, double t=0);
  int Q_order(ParticlePair const &partPair);

private:
  CPOddHiggsLowMass *CPOHLM;

  double mi, mj;
  double mc, Mc, dc;
  int i, j; // index of the chi particles

  int i_model;
  double T_over_TQCD;

  double lPScc, alpha;

};

#endif