#ifndef DEF_CROSSSECTIONMAJORANA_CCCC
#define DEF_CROSSSECTIONMAJORANA_CCCC

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSectionMajorana_cccc : public CrossSection
{

public:
  CrossSectionMajorana_cccc(DarkSectorModel const &DSMod, Particle const &chi);
  virtual ~CrossSectionMajorana_cccc(){};

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

  dcomp Q_suphisphis(int n, Particle const &s1, Particle const &s2, double s);
  dcomp Q_suphipphip(int n, Particle const &ps1, Particle const &ps2, double s);
  dcomp Q_suXXgg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_suXXgppg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_suXXpp(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_suXXpg(int n, Particle const &X1, Particle const &X2, double s);
  dcomp Q_suphisphip(int n, Particle const &s1, Particle const &ps2, double s);
  dcomp Q_suphisXg(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_suphisXp(int n, Particle const &s1, Particle const &X2, double s);
  dcomp Q_suphipXg(int n, Particle const &ps1, Particle const &X2, double s);
  dcomp Q_suphipXp(int n, Particle const &ps1, Particle const &X2, double s);

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

  dcomp Q(ParticlePair const &partPair, int n, double s, double t = 0);
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