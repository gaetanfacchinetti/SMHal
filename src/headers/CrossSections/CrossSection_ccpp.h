#ifndef DEF_CROSSSECTION_CCPP
#define DEF_CROSSSECTION_CCPP

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSection_ccpp : public CrossSection
{

public:
  CrossSection_ccpp(DarkSectorModel const& DSMod, Particle const& n_chii, Particle const& n_chij, Particle const& phipa, Particle const& phipb);
  virtual ~CrossSection_ccpp(){};

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
  double mi, mj, mpa, mpb;
  double mc, Mc, dc;
  int i,j;  // index of the chi particles
  int a,b;
  
  std::vector<std::vector<std::vector<double> > > lPScc, dspp, lScc;
  double v;


};

#endif