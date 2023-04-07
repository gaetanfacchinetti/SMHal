#ifndef DEF_CROSSSECTIONDIRAC_CCCC
#define DEF_CROSSSECTIONDIRAC_CCCC

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSectionDirac_cccc : public CrossSection
{

public:
  CrossSectionDirac_cccc(DarkSectorModel const &DSMod, Particle const &chi);
  virtual ~CrossSectionDirac_cccc(){};


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

  dcomp Q_uuphisphis(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_uuphipphip(int n, Particle const &ps1, Particle const &ps2, double s);
  dcomp Q_uuXXgg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_uuXXgppg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_uuXXpp(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_uuXXpg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_uuphisphip(int n, Particle const &s1, Particle const &ps2, double s);
  dcomp Q_uuphisXg(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_uuphisXp(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_uuphipXg(int n, Particle const &ps1, Particle const &X2, double s);
  dcomp Q_uuphipXp(int n, Particle const &ps1, Particle const &X2, double s);

  dcomp Q_tuphisphis(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_tuphipphip(int n, Particle const &ps1, Particle const &ps2, double s);
  dcomp Q_tuXXgg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_tuXXgppg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_tuXXpp(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_tuXXpg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_tuphisphip(int n, Particle const &s1, Particle const &ps2, double s);
  dcomp Q_tuphisXg(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_tuphisXp(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_tuphipXg(int n, Particle const &ps1, Particle const &X2, double s);
  dcomp Q_tuphipXp(int n, Particle const &ps1, Particle const &X2, double s);

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