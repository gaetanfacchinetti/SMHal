#include "../../headers/CrossSections/CrossSection_ccsp.h"

CrossSection_ccsp::CrossSection_ccsp(DarkSectorModel const& DSMod, Particle const& chii, Particle const& chij, Particle const& phisa, Particle const& phipb)
    : CrossSection(chii, chij, phisa, phipb)
{

    // std::cout << chii.get_index() << std::endl;

    int n_chi = DSMod.get_n_DS_DM();
    int n_PS = DSMod.get_n_DS_phip();

    lPScc = DSMod.get_lambda_PS_DM();
    lScc = DSMod.get_lambda_S_DM();
    dpsp = DSMod.get_d_psp();

    // First cross section
    std::vector<Particle> as, at, au;

    int const np = n_PS;
    int const nt = n_chi;
    int nu = 0;

    // Do not set the u channel if the two incoming particles are not the same
    if ((&chii == &chij && chii.get_fermiontype() == Fermiontype::majorana))
        nu = n_chi;
    else
        nu = 0;

    as.resize(np);
    at.resize(nt);
    au.resize(nu);

    for (int i = 0; i < np; i++)
        as[i] = DSMod.get_phip(i);

    for (int i = 0; i < nt; i++)
        at[i] = DSMod.get_chi(i);

    if ((&chii == &chij && chii.get_fermiontype() == Fermiontype::majorana))
        for (int i = 0; i < nu; i++)
            au[i] = DSMod.get_chi(i);

    mi = _m1;
    mj = _m2;
    msa = _m3;
    mpb = _m4;

    mc = (mi + mj) / 2;
    Mc = sqrt(mi * mj);
    dc = (mj - mi) / 2;

    i = chii.get_specific_index();
    j = chij.get_specific_index();
    a = phisa.get_specific_index();
    b = phipb.get_specific_index();

    Initialise(as, at, au, "ccsp");

    v = 246.;
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccsp::Q_ssphipphip(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int c2 = p2.get_specific_index();

  	if (n == 0)
		return (-8*std::pow(mc,2) + 4*std::pow(Mc,2) + 4*mi*mj + 2*s)*std::pow(v,2)*dpsp[c1][a][b]*dpsp[c2][a][b]*lPScc[c1][i][j]*lPScc[c2][i][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccsp::Q_ttchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (-4*std::pow(Mc,4) + 8*std::pow(mc,2)*std::pow(mi,2) - 4*std::pow(Mc,2)*std::pow(mi,2) - 2*std::pow(mi,4) + 8*std::pow(mc,2)*mi*ml - 4*std::pow(Mc,2)*mi*ml - 2*std::pow(mi,3)*ml + 2*std::pow(mi,2)*mj*ml - 4*mi*std::pow(mj,2)*ml + 2*std::pow(mi,2)*std::pow(mpb,2) + 2*mi*ml*std::pow(mpb,2) + 8*std::pow(mc,2)*mi*mr - 4*std::pow(Mc,2)*mi*mr - 2*std::pow(mi,3)*mr + 2*std::pow(mi,2)*mj*mr - 4*mi*std::pow(mj,2)*mr + 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr + 2*mi*std::pow(mpb,2)*mr - 8*std::pow(mc,2)*std::pow(msa,2) + 4*std::pow(Mc,2)*std::pow(msa,2) + 2*std::pow(mi,2)*std::pow(msa,2) + 4*std::pow(mj,2)*std::pow(msa,2) - 2*mj*ml*std::pow(msa,2) - 2*std::pow(mpb,2)*std::pow(msa,2) - 2*mj*mr*std::pow(msa,2) + 2*ml*mr*s)*lPScc[b][l][j]*lPScc[b][r][j]*lScc[a][l][i]*lScc[a][r][i];

 	if (n == 1)
		return (8*std::pow(mc,2) - 2*std::pow(mi,2) - 2*std::pow(mj,2) - 2*mi*ml + 2*mj*ml + 2*std::pow(mpb,2) - 2*mi*mr + 2*mj*mr + 2*std::pow(msa,2) - 2*s)*lPScc[b][l][j]*lPScc[b][r][j]*lScc[a][l][i]*lScc[a][r][i];

 	if (n == 2)
		return -2*lPScc[b][l][j]*lPScc[b][r][j]*lScc[a][l][i]*lScc[a][r][i];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccsp::Q_uuchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (-4*std::pow(Mc,4) + 4*std::pow(Mc,2)*std::pow(mi,2) - 2*std::pow(mi,4) - 8*std::pow(mc,2)*std::pow(mj,2) + 8*std::pow(Mc,2)*std::pow(mj,2) + 2*std::pow(mi,3)*ml - 8*std::pow(mc,2)*mj*ml + 4*std::pow(Mc,2)*mj*ml - 2*std::pow(mi,2)*mj*ml + 4*mi*std::pow(mj,2)*ml - 8*std::pow(mc,2)*std::pow(mpb,2) + 8*std::pow(Mc,2)*std::pow(mpb,2) + 2*std::pow(mj,2)*std::pow(mpb,2) + 2*mi*ml*std::pow(mpb,2) + 2*std::pow(mi,3)*mr - 8*std::pow(mc,2)*mj*mr + 4*std::pow(Mc,2)*mj*mr - 2*std::pow(mi,2)*mj*mr + 4*mi*std::pow(mj,2)*mr + 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr + 2*mi*std::pow(mpb,2)*mr + 4*std::pow(Mc,2)*std::pow(msa,2) - 2*std::pow(mj,2)*std::pow(msa,2) - 2*mj*ml*std::pow(msa,2) - 2*std::pow(mpb,2)*std::pow(msa,2) - 2*mj*mr*std::pow(msa,2) - 4*std::pow(Mc,2)*s + 2*std::pow(mi,2)*s + 2*std::pow(mj,2)*s - 2*mi*ml*s + 2*mj*ml*s - 2*mi*mr*s + 2*mj*mr*s + 2*ml*mr*s)*lPScc[b][l][i]*lPScc[b][r][i]*lScc[a][l][j]*lScc[a][r][j];

 	if (n == 1)
		return (8*std::pow(mc,2) - 8*std::pow(Mc,2) + 2*std::pow(mi,2) + 2*std::pow(mj,2) - 2*mi*ml + 2*mj*ml + 2*std::pow(mpb,2) - 2*mi*mr + 2*mj*mr + 2*std::pow(msa,2) - 2*s)*lPScc[b][l][i]*lPScc[b][r][i]*lScc[a][l][j]*lScc[a][r][j];

 	if (n == 2)
		return -2*lPScc[b][l][i]*lPScc[b][r][i]*lScc[a][l][j]*lScc[a][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccsp::Q_stphipchi(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double mr = p2.get_mass();

  	if (n == 0)
		return (8*std::pow(mc,2)*mi - 4*std::pow(Mc,2)*mi - 2*std::pow(mi,3) + 2*std::pow(mi,2)*mj - 4*mi*std::pow(mj,2) + 2*mi*std::pow(mpb,2) + 4*std::pow(Mc,2)*mr - 2*std::pow(mi,2)*mr - 2*std::pow(mj,2)*mr - 2*mj*std::pow(msa,2) + 2*mr*s)*v*dpsp[c1][a][b]*lPScc[b][r][j]*lPScc[c1][i][j]*lScc[a][r][i];

 	if (n == 1)
		return (-2*mi + 2*mj)*v*dpsp[c1][a][b]*lPScc[b][r][j]*lPScc[c1][i][j]*lScc[a][r][i];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccsp::Q_suphipchi(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double mr = p2.get_mass();

  	if (n == 0)
		return (2*std::pow(mi,3) - 8*std::pow(mc,2)*mj + 4*std::pow(Mc,2)*mj - 2*std::pow(mi,2)*mj + 4*mi*std::pow(mj,2) + 2*mi*std::pow(mpb,2) + 4*std::pow(Mc,2)*mr - 2*std::pow(mi,2)*mr - 2*std::pow(mj,2)*mr - 2*mj*std::pow(msa,2) - 2*mi*s + 2*mj*s + 2*mr*s)*v*dpsp[c1][a][b]*lPScc[b][r][i]*lPScc[c1][i][j]*lScc[a][r][j];

 	if (n == 1)
		return (-2*mi + 2*mj)*v*dpsp[c1][a][b]*lPScc[b][r][i]*lPScc[c1][i][j]*lScc[a][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccsp::Q_tuchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (-16*std::pow(mc,4) + 16*std::pow(mc,2)*std::pow(Mc,2) - 4*std::pow(mc,2)*std::pow(mi,2) + 2*std::pow(mi,4) + 4*std::pow(mc,2)*std::pow(mj,2) - 4*std::pow(Mc,2)*std::pow(mj,2) + 2*std::pow(mi,3)*ml - 8*std::pow(mc,2)*mj*ml + 4*std::pow(Mc,2)*mj*ml - 2*std::pow(mi,2)*mj*ml + 4*mi*std::pow(mj,2)*ml - 4*std::pow(mc,2)*std::pow(mpb,2) + 4*std::pow(Mc,2)*std::pow(mpb,2) + std::pow(mi,2)*std::pow(mpb,2) + std::pow(mj,2)*std::pow(mpb,2) + 2*mi*ml*std::pow(mpb,2) + 8*std::pow(mc,2)*mi*mr - 4*std::pow(Mc,2)*mi*mr - 2*std::pow(mi,3)*mr + 2*std::pow(mi,2)*mj*mr - 4*mi*std::pow(mj,2)*mr + 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr + 2*mi*std::pow(mpb,2)*mr - 4*std::pow(mc,2)*std::pow(msa,2) + 4*std::pow(Mc,2)*std::pow(msa,2) + std::pow(mi,2)*std::pow(msa,2) + std::pow(mj,2)*std::pow(msa,2) - 2*mj*ml*std::pow(msa,2) - 2*std::pow(mpb,2)*std::pow(msa,2) - 2*mj*mr*std::pow(msa,2) + 8*std::pow(mc,2)*s - 4*std::pow(Mc,2)*s - 2*std::pow(mi,2)*s - 2*std::pow(mj,2)*s - 2*mi*ml*s + 2*mj*ml*s + 2*ml*mr*s)*lPScc[b][l][j]*lPScc[b][r][i]*lScc[a][l][i]*lScc[a][r][j];

 	if (n == 1)
		return (8*std::pow(mc,2) - 4*std::pow(Mc,2) - 2*mi*ml + 2*mj*ml + 2*std::pow(mpb,2) - 2*mi*mr + 2*mj*mr + 2*std::pow(msa,2) - 2*s)*lPScc[b][l][j]*lPScc[b][r][i]*lScc[a][l][i]*lScc[a][r][j];

 	if (n == 2)
		return -2*lPScc[b][l][j]*lPScc[b][r][i]*lScc[a][l][i]*lScc[a][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



// Generated with Mathematica - FeynCalc
int CrossSection_ccsp::Q_order(ParticlePair const&pp)
{
	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return -1;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::dm)
		return -1;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::dm)
		return -1;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	return -1;
}

// Generated with Mathematica - FeynCalc (does not include massless vector here)
dcomp CrossSection_ccsp::Q(ParticlePair const&pp, int n, double s, double t)
{
	if (pp.get_type()== INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_ssphipphip(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_ttchichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_uuchichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return Q_stphipchi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return Q_suphipchi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_tuchichi(n, pp.get_first(), pp.get_second(), s);

	return -1;
}
