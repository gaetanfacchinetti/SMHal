#ifndef CROSSSECTIONSCATT_CGACGA
#define CROSSSECTIONSCATT_CGACGA

#include "CrossSection.h"
#include "../CppLib.h"
#include "../mymath.h"
#include "../DarkSectorModel.h"

class CrossSection_cgacga : public CrossSection
{

public:
    CrossSection_cgacga(DarkSectorModel const &DSMod, Particle const &n_chii);
    virtual ~CrossSection_cgacga(){};

    dcomp Q_ttphisphis(int n, Particle const &s1, Particle const &s2, double s, double t);
    dcomp Q_ttphipphip(int n, Particle const &ps1, Particle const &ps2, double s, double t);
    dcomp Q_ttphisphip(int n, Particle const &s1, Particle const &ps2, double s, double t);

    dcomp Q(ParticlePair const &partPair, int n, double s, double t=0);
    int Q_order(ParticlePair const &partPair);

    dcomp Fphis(int i, double q2); ///< Loop factor for a coupling to scalar mediators
    dcomp Fphip(int i, double q2); ///< Loop factor for a coupling to pseudo-scalar mediators
    double aem(double q2) {return ALPHA_EM;}; ///< \f$\alpha_{\rm em}\f$ (we assume no energy dependance here for the moment)

private:
    double mi, mf;

    std::vector<std::vector<double>> lSff, lPSff;
    std::vector<double> lScc, lPScc;
};

#endif