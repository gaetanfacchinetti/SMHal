#ifndef DEF_CROSSSECTION_CCXX
#define DEF_CROSSSECTION_CCXX

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSection_ccXX : public CrossSection
{

public:
  CrossSection_ccXX(DarkSectorModel const &DSMod, Particle const &n_chii, Particle const &n_chij, Particle const &Xa, Particle const &Xb);
  virtual ~CrossSection_ccXX(){};

  dcomp Q_ttchichi(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_tuchichi(int n, Particle const &pl, Particle const &pr, double s);
  dcomp Q_uuchichi(int n, Particle const &pl, Particle const &pr, double s);

  dcomp Q(ParticlePair const &partPair, int n, double s, double t=0);
  int Q_order(ParticlePair const &partPair);

  void printErrorMessage(std::string func_name)
  {
    std::cout << "FATAL ERROR : " << func_name << " -> trying to access non defined coeff." << std::endl;
    exit(0);
  };

private:
  double mi, mj, mXa, mXb;
  double mc, Mc, dc;
  int i, j, a, b; // index of the chi particles

  std::vector<std::vector<std::vector<double> > > aXcc, bXcc;
};

#endif