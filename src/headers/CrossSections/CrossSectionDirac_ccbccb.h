#ifndef DEF_CROSSSECTIONDIRAC_CCBCCB
#define DEF_CROSSSECTIONDIRAC_CCBCCB

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSectionDirac_ccbccb : public CrossSection
{

public:
  CrossSectionDirac_ccbccb(DarkSectorModel const &DSMod, Particle const &chi, Particle const & chib);
  virtual ~CrossSectionDirac_ccbccb(){};

  dcomp Q_ssphisphis(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_ssphipphip(int n, Particle const &ps1, Particle const &ps2, double s);
  dcomp Q_ssXXgg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ssXXgppg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ssXXpp(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ssXXpg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ssphisphip(int n, Particle const &s1, Particle const &ps2, double s);
  dcomp Q_ssphisXg(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_ssphisXp(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_ssphipXg(int n, Particle const &ps1, Particle const &X2, double s);
  dcomp Q_ssphipXp(int n, Particle const &ps1, Particle const &X2, double s);

  dcomp Q_ttphisphis(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_ttphipphip(int n, Particle const &ps1, Particle const &ps2, double s);
  dcomp Q_ttXXgg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ttXXgppg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ttXXpp(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ttXXpg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_ttphisphip(int n, Particle const &s1, Particle const &ps2, double s);
  dcomp Q_ttphisXg(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_ttphisXp(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_ttphipXg(int n, Particle const &ps1, Particle const &X2, double s);
  dcomp Q_ttphipXp(int n, Particle const &ps1, Particle const &X2, double s);


  dcomp Q_stphisphis(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_stphipphip(int n, Particle const &ps1, Particle const &ps2, double s);
  dcomp Q_stXXgg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_stXXgppg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_stXXpp(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_stXXpg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_stphisphip(int n, Particle const &s1, Particle const &ps2, double s);
  dcomp Q_stphisXg(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_stphisXp(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_stphipXg(int n, Particle const &ps1, Particle const &X2, double s);
  dcomp Q_stphipXp(int n, Particle const &ps1, Particle const &X2, double s);


  dcomp Q(ParticlePair const &partPair, int n, double s, double t=0);
  int Q_order(ParticlePair const &partPair);

  void printErrorMessage(std::string func_name)
  {
    std::cout << "FATAL ERROR : " << func_name << " -> trying to access non defined coeff." << std::endl;
    exit(0);
  };

private:
  double mc;
  int chi_index;

  std::vector<double> lScc, lPScc, aXcc, bXcc;
};

#endif