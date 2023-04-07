#ifndef DEF_CROSSSECTION_CCSS
#define DEF_CROSSSECTION_CCSS

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSection_ccss : public CrossSection
{

public:
  CrossSection_ccss(DarkSectorModel const& DSMod, Particle const& n_chii, Particle const& n_chij, Particle const& phisa, Particle const& phisb);
  virtual ~CrossSection_ccss(){};

  dcomp Q_ssphisphis(int n, Particle const& s1, Particle const& s2, double s);
  dcomp Q_ttchichi(int n, Particle const& pl, Particle const& pr, double s);
  dcomp Q_uuchichi(int n, Particle const& pl, Particle const& pr, double s);
  dcomp Q_stphischi(int n, Particle const& s1, Particle const& pl, double s);
  dcomp Q_suphischi(int n, Particle const& s1, Particle const& pl, double s);
  dcomp Q_tuchichi(int n, Particle const& pl, Particle const& pr, double s);
 
  dcomp Q(ParticlePair const& partPair, int n, double s, double t=0);
  int Q_order(ParticlePair const& partPair);

  void printErrorMessage(std::string func_name)
  {
    std::cout << "FATAL ERROR : " << func_name << " -> trying to access non defined coeff." << std::endl;
    exit(0);
  };

  

private:
  double mi, mj, msa, msb;
  double mc, Mc, dc;
  int i,j;  // index of the chi particles
  int a,b;  // index of the phis particles
  
  std::vector<std::vector<std::vector<double>>> lScc, csss;

  double v;
};

#endif