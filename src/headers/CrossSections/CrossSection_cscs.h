#ifndef CROSSSECTIONSCATT_CSCS
#define CROSSSECTIONSCATT_CSCS

#include "CrossSection.h"
#include "../DarkSectorModel.h"

class CrossSection_cscs : public CrossSection
{

public:
    CrossSection_cscs(DarkSectorModel const &DSMod, Particle const &chii, Particle const &phisa, Particle const &chij, Particle const &phisb);
    virtual ~CrossSection_cscs(){};

    dcomp Q_sschichi(int n, Particle const &c1, Particle const &c2, double s);
    dcomp Q_ttphisphis(int n, Particle const &s1, Particle const &s2, double s);
    dcomp Q_uuchichi(int n, Particle const &c1, Particle const &c2, double s);
    dcomp Q_stchiphis(int n, Particle const &c1, Particle const &s2, double s);
    dcomp Q_suchichi(int n, Particle const &c1, Particle const &c2, double s);
    dcomp Q_tuphischi(int n, Particle const &s1, Particle const &c2, double s);

    dcomp Q(ParticlePair const &partPair, int n, double s, double t=0);
    int Q_order(ParticlePair const &partPair);

private:
    double mi, mj, msa, msb;
    double mc, Mc, dc;
    int i, j; // index of the chi particles
    int a, b; // index of the phis particles

    std::vector<std::vector<std::vector<double> > > lScc, csss;

    double w;
    // Additional attribute for specific cases
    //double mS, mPS, wS, wPS;
};

#endif