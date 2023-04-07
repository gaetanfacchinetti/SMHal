#include "../../headers/CrossSections/CrossSection_cscs.h"

CrossSection_cscs::CrossSection_cscs(DarkSectorModel const& DSMod, Particle const& chii, Particle const& phisa, Particle const& chij, Particle const& phisb)
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

    int const ns = n_chi;
    int const nt = n_S;
    int nu = 0;

    // We get a u-channel if the two fermions or the two scalars are the same
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
        as[i] = DSMod.get_chi(i);

    for (int i = 0; i < nt; i++)
        at[i] = DSMod.get_phis(i);

    if ((&chii == &chij && chii.get_fermiontype() == Fermiontype::majorana) || &phisa == &phisb)
        for (int i = 0; i < nu; i++)
            au[i] = DSMod.get_chi(i);

    mi = _m1;
    mj = _m3;
    msa = _m2;
    msb = _m4;

    i = chii.get_specific_index();
    j = chij.get_specific_index();
    a = phisa.get_specific_index();
    b = phisb.get_specific_index();

    Initialise(as, at, au, "cscs");

    w = 246.;
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cscs::Q_sschichi(int n, Particle const &c1, Particle const &c2, double s)
{
	int l = c1.get_specific_index(); 
	int r = c2.get_specific_index();

  	if (n == 0)
		return 2*std::pow(mi + mj,2)*lScc[a][i][l]*lScc[a][i][r]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 1)
		return -2*lScc[a][i][l]*lScc[a][i][r]*lScc[b][l][j]*lScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cscs::Q_ttphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	int l = s1.get_specific_index(); 
	int r = s2.get_specific_index();

  	if (n == 0)
		return 2*std::pow(mi + mj,2)*std::pow(w,2)*csss[l][a][b]*csss[r][a][b]*lScc[l][i][j]*lScc[r][i][j];

 	if (n == 1)
		return -2*std::pow(w,2)*csss[l][a][b]*csss[r][a][b]*lScc[l][i][j]*lScc[r][i][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cscs::Q_uuchichi(int n, Particle const &c1, Particle const &c2, double s)
{
	int l = c1.get_specific_index(); 
	int r = c2.get_specific_index();

  	if (n == 0)
		return 2*std::pow(mi + mj,2)*lScc[a][i][l]*lScc[a][i][r]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 1)
		return -2*lScc[a][i][l]*lScc[a][i][r]*lScc[b][l][j]*lScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cscs::Q_stchiphis(int n, Particle const &c1, Particle const &s2, double s)
{
	int l = c1.get_specific_index(); 
	int r = s2.get_specific_index();

  	if (n == 0)
		return 2*std::pow(mi + mj,2)*w*csss[r][a][b]*lScc[a][i][l]*lScc[b][l][j]*lScc[r][i][j];

 	if (n == 1)
		return -2*w*csss[r][a][b]*lScc[a][i][l]*lScc[b][l][j]*lScc[r][i][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cscs::Q_suchichi(int n, Particle const &c1, Particle const &c2, double s)
{
	int l = c1.get_specific_index(); 
	int r = c2.get_specific_index();

  	if (n == 0)
		return 2*std::pow(mi + mj,2)*lScc[a][i][l]*lScc[a][i][r]*lScc[b][l][j]*lScc[b][r][j];

 	if (n == 1)
		return -2*lScc[a][i][l]*lScc[a][i][r]*lScc[b][l][j]*lScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cscs::Q_tuphischi(int n, Particle const &s1, Particle const &c2, double s)
{
	int l = s1.get_specific_index(); 
	int r = c2.get_specific_index();

  	if (n == 0)
		return 2*std::pow(mi + mj,2)*w*csss[l][a][b]*lScc[a][i][r]*lScc[b][r][j]*lScc[l][i][j];

 	if (n == 1)
		return -2*w*csss[l][a][b]*lScc[a][i][r]*lScc[b][r][j]*lScc[l][i][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



// Generated with Mathematica - FeynCalc
int CrossSection_cscs::Q_order(ParticlePair const &pp)
{
	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 1;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 1;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 1;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::scalar)
		return 1;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 1;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return 1;

	return -1;
}

// Generated with Mathematica - FeynCalc (does not include massless vector here)
dcomp CrossSection_cscs::Q(ParticlePair const &pp, int n, double s, double t)
{
	if (pp.get_type()== INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_sschichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_ttphisphis(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_uuchichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_stchiphis(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_suchichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return Q_tuphischi(n, pp.get_first(), pp.get_second(), s);

	return -1;
}
