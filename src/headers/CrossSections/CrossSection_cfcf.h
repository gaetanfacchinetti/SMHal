#ifndef CROSSSECTIONSCATT_CFCF
#define CROSSSECTIONSCATT_CFCF

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSection_cfcf : public CrossSection
{

public:
  CrossSection_cfcf(DarkSectorModel const& DSMod, Particle const& n_chii, ParticlePair const& couple_SM_ferm);
  virtual ~CrossSection_cfcf(){};

  dcomp Q_ttphisphis(int n, Particle const& s1, Particle const& s2, double s);
  dcomp Q_ttphipphip(int n, Particle const& ps1, Particle const& ps2,  double s);
  dcomp Q_ttphisphip(int n, Particle const& s1, Particle const& ps2,  double s);
  dcomp Q_ttXX(int n,  Particle const& X1, Particle const& X2, double s);
  dcomp Q_ttphisX(int n,  Particle const& s1, Particle const& X2, double s);
  dcomp Q_ttphipX(int n,  Particle const& ps1, Particle const& X2, double s);


  dcomp Q(ParticlePair const &partPair, int n, double s, double t=0);
  int Q_order(ParticlePair const &partPair);



private:
  double mi, mf;

  std::vector<double> lSff, lPSff, aXff, bXff;
  std::vector<double> lScc, lPScc, aXcc, bXcc;

  std::vector<double> mX;

  // Additional attribute for specific cases
  //double mS, mPS, wS, wPS;
};

#endif