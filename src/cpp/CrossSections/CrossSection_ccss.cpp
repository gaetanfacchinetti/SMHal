#include "../../headers/CrossSections/CrossSection_ccss.h"

CrossSection_ccss::CrossSection_ccss(DarkSectorModel const& DSMod, Particle const& chii, Particle const& chij, Particle const& phisa, Particle const& phisb)
    : CrossSection(chii, chij, phisa, phisb)
{

    // std::cout << chii.get_index() << std::endl;

    int n_chi = DSMod.get_n_DS_DM();
    int n_S = DSMod.get_n_DS_phis();

    lScc = DSMod.get_lambda_S_DM();
    csss = DSMod.get_c_sss();

    //std::cout << cSCC[0][0] << std::endl;

    // First cross section
    std::vector<Particle> as, at, au;

    int const ns = n_S;
    int const nt = n_chi;
    int nu = 0;

    // UPDATE -> what is this former comment : "// Avoid the u channel if the two incoming and outgoing particles are the same" ?
    if  ((&chii == &chij && chii.get_fermiontype() == Fermiontype::majorana) || &phisa == &phisb)
	{
		//std::cout << "here" << std::endl;
        nu = n_chi;
	}
    else
        nu = 0;

    as.resize(ns);
    at.resize(nt);
    au.resize(nu);

    for (int i = 0; i < ns; i++)
        as[i] = DSMod.get_phis(i);

    for (int i = 0; i < nt; i++)
        at[i] = DSMod.get_chi(i);

    if ((&chii == &chij && chii.get_fermiontype() == Fermiontype::majorana) || &phisa == &phisb)
        for (int i = 0; i < nu; i++)
            au[i] = DSMod.get_chi(i);

    mi = _m1;
    mj = _m2;
    msa = _m3;
    msb = _m4;

    mc = (mi + mj) / 2;
    Mc = sqrt(mi * mj);
    dc = (mj - mi) / 2;

    i = chii.get_specific_index();
    j = chij.get_specific_index();
    a = phisa.get_specific_index();
    b = phisb.get_specific_index();

    Initialise(as, at, au, "ccss");

    v = 246.;
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccss::Q_ssphisphis(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int c2 = p2.get_specific_index();

  	if (n == 0)
		return (-8*std::pow(mc,2) + 4*std::pow(Mc,2) - 4*mi*mj + 2*s)*std::pow(v,2)*csss[c1][a][b]*csss[c2][a][b]*lScc[c1][i][j]*lScc[c2][i][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccss::Q_ttchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (-4*std::pow(Mc,4) + 8*std::pow(mc,2)*std::pow(mi,2) - 4*std::pow(Mc,2)*std::pow(mi,2) - 2*std::pow(mi,4) + 8*std::pow(mc,2)*mi*ml - 4*std::pow(Mc,2)*mi*ml - 2*std::pow(mi,3)*ml - 2*std::pow(mi,2)*mj*ml - 4*mi*std::pow(mj,2)*ml + 8*std::pow(mc,2)*mi*mr - 4*std::pow(Mc,2)*mi*mr - 2*std::pow(mi,3)*mr - 2*std::pow(mi,2)*mj*mr - 4*mi*std::pow(mj,2)*mr - 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr - 8*std::pow(mc,2)*std::pow(msa,2) + 4*std::pow(Mc,2)*std::pow(msa,2) + 2*std::pow(mi,2)*std::pow(msa,2) + 4*std::pow(mj,2)*std::pow(msa,2) + 2*mj*ml*std::pow(msa,2) + 2*mj*mr*std::pow(msa,2) + 2*std::pow(mi,2)*std::pow(msb,2) + 2*mi*ml*std::pow(msb,2) + 2*mi*mr*std::pow(msb,2) - 2*std::pow(msa,2)*std::pow(msb,2) + 2*ml*mr*s)*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 1)
		return (8*std::pow(mc,2) - 8*std::pow(Mc,2) - 2*std::pow(mi,2) - 2*std::pow(mj,2) - 2*mi*ml - 2*mj*ml - 2*mi*mr - 2*mj*mr + 2*std::pow(msa,2) + 2*std::pow(msb,2) - 2*s)*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 2)
		return -2*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccss::Q_uuchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (-4*std::pow(Mc,4) - 4*std::pow(Mc,2)*std::pow(mi,2) - 2*std::pow(mi,4) - 8*std::pow(mc,2)*std::pow(mj,2) - 2*std::pow(mi,3)*ml - 8*std::pow(mc,2)*mj*ml + 4*std::pow(Mc,2)*mj*ml - 2*std::pow(mi,2)*mj*ml - 4*mi*std::pow(mj,2)*ml - 2*std::pow(mi,3)*mr - 8*std::pow(mc,2)*mj*mr + 4*std::pow(Mc,2)*mj*mr - 2*std::pow(mi,2)*mj*mr - 4*mi*std::pow(mj,2)*mr - 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr - 4*std::pow(Mc,2)*std::pow(msa,2) - 2*std::pow(mj,2)*std::pow(msa,2) - 2*mj*ml*std::pow(msa,2) - 2*mj*mr*std::pow(msa,2) - 8*std::pow(mc,2)*std::pow(msb,2) + 2*std::pow(mj,2)*std::pow(msb,2) - 2*mi*ml*std::pow(msb,2) - 2*mi*mr*std::pow(msb,2) - 2*std::pow(msa,2)*std::pow(msb,2) + 4*std::pow(Mc,2)*s + 2*std::pow(mi,2)*s + 2*std::pow(mj,2)*s + 2*mi*ml*s + 2*mj*ml*s + 2*mi*mr*s + 2*mj*mr*s + 2*ml*mr*s)*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 1)
		return (8*std::pow(mc,2) + 2*std::pow(mi,2) + 2*std::pow(mj,2) + 2*mi*ml + 2*mj*ml + 2*mi*mr + 2*mj*mr + 2*std::pow(msa,2) + 2*std::pow(msb,2) - 2*s)*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 2)
		return -2*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccss::Q_stphischi(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double mr = p2.get_mass();

  	if (n == 0)
		return (8*std::pow(mc,2)*mi - 4*std::pow(Mc,2)*mi - 2*std::pow(mi,3) - 2*std::pow(mi,2)*mj - 4*mi*std::pow(mj,2) - 4*std::pow(Mc,2)*mr - 2*std::pow(mi,2)*mr - 2*std::pow(mj,2)*mr + 2*mj*std::pow(msa,2) + 2*mi*std::pow(msb,2) + 2*mr*s)*v*csss[c1][a][b]*lScc[a][r][i]*lScc[b][r][j]*lScc[c1][i][j];

 	if (n == 1)
		return (-2*mi - 2*mj)*v*csss[c1][a][b]*lScc[a][r][i]*lScc[b][r][j]*lScc[c1][i][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccss::Q_suphischi(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double mr = p2.get_mass();

  	if (n == 0)
		return (-2*std::pow(mi,3) - 8*std::pow(mc,2)*mj + 4*std::pow(Mc,2)*mj - 2*std::pow(mi,2)*mj - 4*mi*std::pow(mj,2) - 4*std::pow(Mc,2)*mr - 2*std::pow(mi,2)*mr - 2*std::pow(mj,2)*mr - 2*mj*std::pow(msa,2) - 2*mi*std::pow(msb,2) + 2*mi*s + 2*mj*s + 2*mr*s)*v*csss[c1][a][b]*lScc[a][r][i]*lScc[b][r][j]*lScc[c1][i][j];

 	if (n == 1)
		return (2*mi + 2*mj)*v*csss[c1][a][b]*lScc[a][r][i]*lScc[b][r][j]*lScc[c1][i][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccss::Q_tuchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (16*std::pow(mc,4) - 16*std::pow(mc,2)*std::pow(Mc,2) + 4*std::pow(mc,2)*std::pow(mi,2) - 4*std::pow(Mc,2)*std::pow(mi,2) - 2*std::pow(mi,4) - 4*std::pow(mc,2)*std::pow(mj,2) - 2*std::pow(mi,3)*ml - 8*std::pow(mc,2)*mj*ml + 4*std::pow(Mc,2)*mj*ml - 2*std::pow(mi,2)*mj*ml - 4*mi*std::pow(mj,2)*ml + 8*std::pow(mc,2)*mi*mr - 4*std::pow(Mc,2)*mi*mr - 2*std::pow(mi,3)*mr - 2*std::pow(mi,2)*mj*mr - 4*mi*std::pow(mj,2)*mr - 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr + 4*std::pow(mc,2)*std::pow(msa,2) - std::pow(mi,2)*std::pow(msa,2) - std::pow(mj,2)*std::pow(msa,2) - 2*mj*ml*std::pow(msa,2) + 2*mj*mr*std::pow(msa,2) + 4*std::pow(mc,2)*std::pow(msb,2) - std::pow(mi,2)*std::pow(msb,2) - std::pow(mj,2)*std::pow(msb,2) - 2*mi*ml*std::pow(msb,2) + 2*mi*mr*std::pow(msb,2) + 2*std::pow(msa,2)*std::pow(msb,2) - 8*std::pow(mc,2)*s + 4*std::pow(Mc,2)*s + 2*std::pow(mi,2)*s + 2*std::pow(mj,2)*s + 2*mi*ml*s + 2*mj*ml*s + 2*ml*mr*s)*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 1)
		return (-8*std::pow(mc,2) + 4*std::pow(Mc,2) + 2*mi*ml + 2*mj*ml - 2*mi*mr - 2*mj*mr - 2*std::pow(msa,2) - 2*std::pow(msb,2) + 2*s)*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 2)
		return 2*lScc[a][l][i]*lScc[a][r][i]*lScc[b][l][j]*lScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



// Generated with Mathematica - FeynCalc
int CrossSection_ccss::Q_order(ParticlePair const&pp)
{
	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 0;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return 1;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return 1;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	return -1;
}

// Generated with Mathematica - FeynCalc (does not include massless vector here)
dcomp CrossSection_ccss::Q(ParticlePair const&pp, int n, double s, double t)
{
	if (pp.get_type()== INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_ssphisphis(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_ttchichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_uuchichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return Q_stphischi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return Q_suphischi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_tuchichi(n, pp.get_first(), pp.get_second(), s);

	return -1;
}


/*
dcomp CrossSection_ccss::Q(ParticlePair const& pp, int n, double s)
{
    if (pp.get_type() == INT_TYPE::SS)
        return Q_ssphisphis(n, pp.get_first(), pp.get_second(), s);
    if (pp.get_type() == INT_TYPE::TT)
        return Q_ttchilchir(n, pp.get_first(), pp.get_second(), s);
    if (pp.get_type() == INT_TYPE::UU)
        return Q_uuchilchir(n, pp.get_first(), pp.get_second(), s);
    if (pp.get_type() == INT_TYPE::ST)
        return Q_stphischil(n, pp.get_first(), pp.get_second(), s);
    if (pp.get_type() == INT_TYPE::SU)
        return Q_suphischil(n, pp.get_first(), pp.get_second(), s);
    if (pp.get_type() == INT_TYPE::TU)
        return Q_tuchilchir(n, pp.get_first(), pp.get_second(), s);

    std::cout << "FATAL ERROR : " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}

int CrossSection_ccss::Q_order(ParticlePair const &pp)
{
    if (pp.get_type() == INT_TYPE::SS)
        return 0;
    if (pp.get_type() == INT_TYPE::TT)
        return 2;
    if (pp.get_type() == INT_TYPE::UU)
        return 2;
    if (pp.get_type() == INT_TYPE::ST)
        return 1;
    if (pp.get_type() == INT_TYPE::SU)
        return 1;
    if (pp.get_type() == INT_TYPE::TU)
        return 2;

    std::cout << "FATAL ERROR : " << __PRETTY_FUNCTION__ << std::endl;
    exit(0);
}*/
