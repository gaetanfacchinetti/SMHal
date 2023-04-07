#include "../../headers/CrossSections/CrossSection_ccpp.h"

CrossSection_ccpp::CrossSection_ccpp(DarkSectorModel const& DSMod, Particle const& chii, Particle const& chij, Particle const& phipa, Particle const& phipb)
    : CrossSection(chii, chij, phipa, phipb)
{

    // std::cout << chii.get_index() << std::endl;

    int n_chi = DSMod.get_n_DS_DM();
    int n_S = DSMod.get_n_DS_phis();
    int n_PS = DSMod.get_n_DS_phip();

    lPScc = DSMod.get_lambda_PS_DM();
    lScc = DSMod.get_lambda_S_DM();
    dspp = DSMod.get_d_spp();

    //std::cout << cSCC[0][0] << std::endl;

    // First cross section
    std::vector<Particle> as, at, au;

    int const ns = n_S;
    int const nt = n_chi;
    int nu = 0;

    //std::cout << &chii << " " << &chij << " " <<  &phipa << " " <<  &phipb << std::endl;

    // Avoid the u channel if the two ingoing particles are the same
    if ((&chii == &chij && chii.get_fermiontype() == Fermiontype::majorana) || &phipa == &phipb)
    {
        //std::cout << "here" << std::endl;
        nu = n_chi;
    }
    else
        nu = 0;
    
    as.resize(ns);
    at.resize(nt);
    au.resize(nu);

    for (int ii = 0; ii < ns; ii++)
        as[ii] = DSMod.get_phis(ii);

    for (int ii = 0; ii < nt; ii++)
        at[ii] = DSMod.get_chi(ii);

    if ((&chii == &chij && chii.get_fermiontype() == Fermiontype::majorana) || &phipa == &phipb)
        for (int ii = 0; ii < nu; ii++)
            au[ii] = DSMod.get_chi(ii);


    mi = _m1;
    mj = _m2;
    mpa = _m3;
    mpb = _m4;

    mc = (mi + mj) / 2;
    Mc = sqrt(mi * mj);
    dc = (mj - mi) / 2;

    i = chii.get_specific_index();
    j = chij.get_specific_index();
    a = phipa.get_specific_index();
    b = phipb.get_specific_index();

    Initialise(as, at, au, "ccpp");

    v = 246.;
}


// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccpp::Q_ssphisphis(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int c2 = p2.get_specific_index();

  	if (n == 0)
		return (-8*std::pow(mc,2) + 4*std::pow(Mc,2) - 4*mi*mj + 2*s)*std::pow(v,2)*dspp[c1][a][b]*dspp[c2][a][b]*lScc[c1][i][j]*lScc[c2][i][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccpp::Q_ttchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (-4*std::pow(Mc,4) + 8*std::pow(mc,2)*std::pow(mi,2) - 4*std::pow(Mc,2)*std::pow(mi,2) - 2*std::pow(mi,4) - 8*std::pow(mc,2)*mi*ml + 4*std::pow(Mc,2)*mi*ml + 2*std::pow(mi,3)*ml - 2*std::pow(mi,2)*mj*ml + 4*mi*std::pow(mj,2)*ml - 8*std::pow(mc,2)*std::pow(mpa,2) + 4*std::pow(Mc,2)*std::pow(mpa,2) + 2*std::pow(mi,2)*std::pow(mpa,2) + 4*std::pow(mj,2)*std::pow(mpa,2) + 2*mj*ml*std::pow(mpa,2) + 2*std::pow(mi,2)*std::pow(mpb,2) - 2*mi*ml*std::pow(mpb,2) - 2*std::pow(mpa,2)*std::pow(mpb,2) - 8*std::pow(mc,2)*mi*mr + 4*std::pow(Mc,2)*mi*mr + 2*std::pow(mi,3)*mr - 2*std::pow(mi,2)*mj*mr + 4*mi*std::pow(mj,2)*mr + 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr + 2*mj*std::pow(mpa,2)*mr - 2*mi*std::pow(mpb,2)*mr + 2*ml*mr*s)*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

 	if (n == 1)
		return (8*std::pow(mc,2) - 2*std::pow(mi,2) - 2*std::pow(mj,2) + 2*mi*ml - 2*mj*ml + 2*std::pow(mpa,2) + 2*std::pow(mpb,2) + 2*mi*mr - 2*mj*mr - 2*s)*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

 	if (n == 2)
		return -2*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccpp::Q_uuchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (-4*std::pow(Mc,4) + 4*std::pow(Mc,2)*std::pow(mi,2) - 2*std::pow(mi,4) - 8*std::pow(mc,2)*std::pow(mj,2) + 8*std::pow(Mc,2)*std::pow(mj,2) + 2*std::pow(mi,3)*ml - 8*std::pow(mc,2)*mj*ml + 4*std::pow(Mc,2)*mj*ml - 2*std::pow(mi,2)*mj*ml + 4*mi*std::pow(mj,2)*ml + 4*std::pow(Mc,2)*std::pow(mpa,2) - 2*std::pow(mj,2)*std::pow(mpa,2) - 2*mj*ml*std::pow(mpa,2) - 8*std::pow(mc,2)*std::pow(mpb,2) + 8*std::pow(Mc,2)*std::pow(mpb,2) + 2*std::pow(mj,2)*std::pow(mpb,2) + 2*mi*ml*std::pow(mpb,2) - 2*std::pow(mpa,2)*std::pow(mpb,2) + 2*std::pow(mi,3)*mr - 8*std::pow(mc,2)*mj*mr + 4*std::pow(Mc,2)*mj*mr - 2*std::pow(mi,2)*mj*mr + 4*mi*std::pow(mj,2)*mr + 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr - 2*mj*std::pow(mpa,2)*mr + 2*mi*std::pow(mpb,2)*mr - 4*std::pow(Mc,2)*s + 2*std::pow(mi,2)*s + 2*std::pow(mj,2)*s - 2*mi*ml*s + 2*mj*ml*s - 2*mi*mr*s + 2*mj*mr*s + 2*ml*mr*s)*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

 	if (n == 1)
		return (8*std::pow(mc,2) - 8*std::pow(Mc,2) + 2*std::pow(mi,2) + 2*std::pow(mj,2) - 2*mi*ml + 2*mj*ml + 2*std::pow(mpa,2) + 2*std::pow(mpb,2) - 2*mi*mr + 2*mj*mr - 2*s)*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

 	if (n == 2)
		return -2*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccpp::Q_stphischi(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double mr = p2.get_mass();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccpp::Q_suphischi(int n, Particle const &p1, Particle const &p2, double s)
{
	int c1 = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double mr = p2.get_mass();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccpp::Q_tuchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_mass(); 
	double mr = p2.get_mass();

  	if (n == 0)
		return (16*std::pow(mc,4) - 16*std::pow(mc,2)*std::pow(Mc,2) + 4*std::pow(mc,2)*std::pow(mi,2) - 2*std::pow(mi,4) - 4*std::pow(mc,2)*std::pow(mj,2) + 4*std::pow(Mc,2)*std::pow(mj,2) + 2*std::pow(mi,3)*ml - 8*std::pow(mc,2)*mj*ml + 4*std::pow(Mc,2)*mj*ml - 2*std::pow(mi,2)*mj*ml + 4*mi*std::pow(mj,2)*ml + 4*std::pow(mc,2)*std::pow(mpa,2) - 4*std::pow(Mc,2)*std::pow(mpa,2) - std::pow(mi,2)*std::pow(mpa,2) - std::pow(mj,2)*std::pow(mpa,2) - 2*mj*ml*std::pow(mpa,2) + 4*std::pow(mc,2)*std::pow(mpb,2) - 4*std::pow(Mc,2)*std::pow(mpb,2) - std::pow(mi,2)*std::pow(mpb,2) - std::pow(mj,2)*std::pow(mpb,2) + 2*mi*ml*std::pow(mpb,2) + 2*std::pow(mpa,2)*std::pow(mpb,2) - 8*std::pow(mc,2)*mi*mr + 4*std::pow(Mc,2)*mi*mr + 2*std::pow(mi,3)*mr - 2*std::pow(mi,2)*mj*mr + 4*mi*std::pow(mj,2)*mr + 4*std::pow(Mc,2)*ml*mr - 2*std::pow(mi,2)*ml*mr - 2*std::pow(mj,2)*ml*mr + 2*mj*std::pow(mpa,2)*mr - 2*mi*std::pow(mpb,2)*mr - 8*std::pow(mc,2)*s + 4*std::pow(Mc,2)*s + 2*std::pow(mi,2)*s + 2*std::pow(mj,2)*s - 2*mi*ml*s + 2*mj*ml*s + 2*ml*mr*s)*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

 	if (n == 1)
		return (-8*std::pow(mc,2) + 4*std::pow(Mc,2) - 2*mi*ml + 2*mj*ml - 2*std::pow(mpa,2) - 2*std::pow(mpb,2) + 2*mi*mr - 2*mj*mr + 2*s)*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

 	if (n == 2)
		return 2*lPScc[a][l][i]*lPScc[a][r][i]*lPScc[b][l][j]*lPScc[b][r][j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



// Generated with Mathematica - FeynCalc
int CrossSection_ccpp::Q_order(ParticlePair const&pp)
{
	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 0;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return -1;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::dm)
		return -1;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	return -1;
}

// Generated with Mathematica - FeynCalc (does not include massless vector here)
dcomp CrossSection_ccpp::Q(ParticlePair const&pp, int n, double s, double t)
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
