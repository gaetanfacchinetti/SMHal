#include "../../headers/CrossSections/CrossSection_ccgaga.h"

CrossSection_ccgaga::CrossSection_ccgaga(DarkSectorModel & DSMod, Particle const& chii, Particle const& chij, CPOddHiggsLowMass *n_CPOHLM, int n_i_model, double n_T_over_TQCD)
    : CrossSection(chii, chij, *(StandardModel::getInstance()->get_photon_ptr()), *(StandardModel::getInstance()->get_photon_ptr()))
{

    // std::cout << chii->get_index() << std::endl;
        
    
    int n_PS = 1;

    CPOHLM = n_CPOHLM;
    i_model = n_i_model;
    T_over_TQCD = n_T_over_TQCD;
    

    lPScc = CPOHLM->ach0A(i_model, T_over_TQCD);


    // First cross section
    std::vector<Particle> as, at, au;

    int const np = n_PS;
    int const nt = 0;
    int const nu = 0;

    as.resize(np);
    at.resize(nt);
    au.resize(nu);

    // redefine the width of the particle according to the tmeperature
    DSMod.get_phip_ptr(0)->set_width(CPOHLM->GammaTAS(i_model, T_over_TQCD));

    // COmplete the table of mediator particles
    as[0] = DSMod.get_phip(0);

    mi = _m1;
    mj = _m2;

    mc = (mi + mj) / 2;
    Mc = sqrt(mi * mj);
    dc = (mj - mi) / 2;

    Initialise(as, at, au, "ccgaga");

    alpha = ALPHA_EM;


}

// Function completed from mathematica results
dcomp CrossSection_ccgaga::Q_ssphipphip(int n, double s)
{
    if(n == 0)
        return (std::pow(alpha,2)*std::pow(lPScc,2)*std::pow(s,2)*(-4*std::pow(dc,2) + s)*std::pow(abs(sFg(s)),2))/std::pow(PI,2);

    std::cout << "FATAL ERROR : " << __PRETTY_FUNCTION__ << " -> trying to access non defined coeff." << std::endl;
    exit(0);
}

dcomp CrossSection_ccgaga::Q(ParticlePair const &pp, int n, double s, double t)
{
    if (pp.get_type() == INT_TYPE::SS)
        return Q_ssphipphip(n, s);

    std::cout << "FATAL ERROR : " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}

int CrossSection_ccgaga::Q_order(ParticlePair const &pp)
{
    if (pp.get_type() == INT_TYPE::SS)
        return 0;

    std::cout << "FATAL ERROR : " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}

double CrossSection_ccgaga::sFg(double s)
{
    return abs(CPOHLM->sFg(i_model, s, T_over_TQCD));
}