#ifndef CROSSSECTION_CCFF
#define CROSSSECTION_CCFF

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSection_ccff : public CrossSection
{

public:
  /**
  * Cross section for the annihilation of DM (Dirac or Majorana) into SM fermions
  * WARNING : Note that this cross-section is summed over the spins but not over the colours (for the quarks).
  */
  CrossSection_ccff(DarkSectorModel const &DSMod, Particle const &n_chii, Particle const &n_chij, ParticlePair const &couple_SM_ferm);
  virtual ~CrossSection_ccff(){};

  /*
  dcomp Q_ssphisphis(int n, Particle const&s1, Particle const&s2, double s);
  dcomp Q_ssphipphip(int n, Particle const&ps1, Particle const&ps2,  double s);
  dcomp Q_ssXXgg(int n,  Particle const&X1, Particle const&X2, double s);
  dcomp Q_ssXXgppg(int n,  Particle const&X1, Particle const&X2, double s);
  dcomp Q_ssXXpp(int n,  Particle const&X1, Particle const&X2, double s);
  dcomp Q_ssXXpg(int n,  Particle const&X1, Particle const&X2, double s);
  dcomp Q_ssphisXg(int n,  Particle const&s1, Particle const&X2, double s);
  dcomp Q_ssphipXg(int n,  Particle const&ps1, Particle const&X2, double s);
  dcomp Q_ssphipXp(int n,  Particle const&ps1, Particle const&X2, double s);*/

  dcomp Q_ssphisphis(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_ssphipphip(int n, Particle const &ps1, Particle const &ps2, double s);
  dcomp Q_ssXX(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ssphisphip(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_ssphisX(int n, Particle const &ps1, Particle const &X2, double s);
  dcomp Q_ssphipX(int n, Particle const &ps1, Particle const &X2, double s);

  dcomp Q(ParticlePair const &partPair, int n, double s, double t=0);
  int Q_order(ParticlePair const &partPair);

  void printErrorMessage(std::string func_name)
  {
    std::cout << "FATAL ERROR : " << func_name << " -> trying to access non defined coeff." << std::endl;
    exit(0);
  };

  //std::vector<double> AppSigmaVChanSPropSPS();

private:
  double mi, mj, mf;
  double mc, Mc, dc;

  std::vector<double> lSff, lPSff, aXff, bXff;
  std::vector<double> lScc, lPScc, aXcc, bXcc;

  std::vector<double> mX;

  // Additional attribute for specific cases
  //double mS, mPS, wS, wPS;
};

#endif