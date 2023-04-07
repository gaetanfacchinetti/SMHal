#ifndef DEF_CROSSSECTION_CCSP
#define DEF_CROSSSECTION_CCSP

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSection_ccsp : public CrossSection
{

public:
  CrossSection_ccsp(DarkSectorModel const &DSMod, Particle const &n_chii, Particle const &n_chij, Particle const &phisa, Particle const &phipb);
  virtual ~CrossSection_ccsp(){};

  dcomp Q_ssphipphip(int n, Particle const &ps1, Particle const &p2, double s);
  dcomp Q_ttchichi(int n, Particle const &pl, Particle const &pr, double s);
  dcomp Q_uuchichi(int n, Particle const &pl, Particle const &pr, double s);
  dcomp Q_stphipchi(int n, Particle const &ps1, Particle const &pl, double s);
  dcomp Q_suphipchi(int n, Particle const &ps1, Particle const &pl, double s);
  dcomp Q_tuchichi(int n, Particle const &pl, Particle const &pr, double s);

  dcomp Q(ParticlePair const &partPair, int n, double s, double t=0);
  int Q_order(ParticlePair const &partPair);

  void printErrorMessage(std::string func_name)
  {
    std::cout << "FATAL ERROR : " << func_name << " -> trying to access non defined coeff." << std::endl;
    exit(0);
  };

private:
  double mi, mj, msa, mpb;
  double mc, Mc, dc;
  int i, j; // index of the chi particles
  int a, b;

  std::vector<std::vector<std::vector<double> > > lPScc, lScc, dpsp;
  double v;
};

#endif