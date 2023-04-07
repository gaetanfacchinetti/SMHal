#include "../../headers/CrossSections/CrossSection_ccXp.h"

CrossSection_ccXp::CrossSection_ccXp(DarkSectorModel const &DSMod, Particle const &chii, Particle const &chij, Particle const &Xa, Particle const &phipb)
    : CrossSection(chii, chij, Xa, phipb)
{
    // std::cout << chii.get_index() << std::endl;
    int n_chi = DSMod.get_n_DS_DM();
    int n_X = DSMod.get_n_DS_X();

    aXcc = DSMod.get_a_VEC_DM();
    bXcc = DSMod.get_b_VEC_DM();
    lPScc = DSMod.get_lambda_PS_DM();

    // First cross section
    std::vector<Particle> as, at, au;

    int const nX = n_X;
    int const nt = n_chi;
    int nu = 0;
    int ns = 0;

    // Avoid the u channel if the two incoming and outgoing particles are the same
    if ((n_chi > 1 && &chii != &chij))
        nu = n_chi;
    else
        nu = 0;

    as.resize(ns);
    at.resize(nt);
    au.resize(nu);

    for (int i = 0; i < nt; i++)
        at[i] = DSMod.get_chi(i);

    if ((n_chi > 1 && &chii != &chij))
        for (int i = 0; i < nu; i++)
            au[i] = DSMod.get_chi(i);

    mi = _m1;
    mj = _m2;
    mXa = _m3;
    mpb = _m4;

    mc = (mi + mj) / 2;
    Mc = sqrt(mi * mj);
    dc = (mj - mi) / 2;

    Initialise(as, at, au, "ccss");

    i = chii.get_specific_index();
    j = chij.get_specific_index();
    a = Xa.get_specific_index();
    b = phipb.get_specific_index();
}


// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccXp::Q_ttchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_specific_index(); 
	double mr = p2.get_specific_index();

  	if (n == 0)
		return -12*std::pow(Mc,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*mj*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*std::pow(mj,2)*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*ml*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*mj*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*std::pow(mj,2)*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*std::pow(mpb,2)*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,2)*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*ml*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,2)*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*mr*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*ml*mr*(4*std::pow(mc,2) + 2*std::pow(Mc,2) - s)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(Mc,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*mj*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*std::pow(mj,2)*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*ml*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*mj*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*std::pow(mj,2)*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*std::pow(mpb,2)*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,2)*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*ml*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,2)*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*mr*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*ml*mr*(4*std::pow(mc,2) - 6*std::pow(Mc,2) - s)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

 	if (n == 1)
		return -48*std::pow(Mc,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*std::pow(Mc,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

 	if (n == 2)
		return -12*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccXp::Q_tuchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_specific_index(); 
	double mr = p2.get_specific_index();

  	if (n == 0)
		return 12*std::pow(Mc,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mi,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*std::pow(mj,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mi,3)*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mi,2)*mj*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*mi*std::pow(mj,2)*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,3)*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,2)*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*ml*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*std::pow(mpb,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*mj*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*std::pow(mj,2)*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*std::pow(mpb,2)*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,2)*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*ml*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mpb,2)*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*mr*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*std::pow(mXa,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(Mc,2)*(4*std::pow(mc,2) - 2*std::pow(Mc,2) - std::pow(mpb,2) - std::pow(mXa,2))*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mi,2)*(std::pow(mj,2) + std::pow(mXa,2))*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mj,2)*(std::pow(mj,2) + std::pow(mXa,2))*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mpb,2)*(std::pow(mj,2) + std::pow(mXa,2))*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mXa,2)*(std::pow(mj,2) + std::pow(mXa,2))*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mpb,2)*(4*std::pow(mc,2) - 2*std::pow(Mc,2) + std::pow(mpb,2) - std::pow(mXa,2) - 2*s)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mi,2)*(-std::pow(mi,2) + 3*std::pow(mj,2) + std::pow(mpb,2) + std::pow(mXa,2) - 2*s)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*ml*mr*(4*std::pow(mc,2) + 2*std::pow(Mc,2) - s)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mj,2)*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*ml*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*ml*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,2)*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mXa,2)*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*(std::pow(mj,2) + std::pow(mXa,2))*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(Mc,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mi,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*std::pow(mj,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mi,3)*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mi,2)*mj*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*mi*std::pow(mj,2)*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,3)*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,2)*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*ml*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*std::pow(mpb,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*mj*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*std::pow(mj,2)*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*std::pow(mpb,2)*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,2)*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*ml*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mpb,2)*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*mr*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*std::pow(mXa,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mi,2)*(std::pow(mj,2) + std::pow(mXa,2))*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mj,2)*(std::pow(mj,2) + std::pow(mXa,2))*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mpb,2)*(std::pow(mj,2) + std::pow(mXa,2))*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mXa,2)*(std::pow(mj,2) + std::pow(mXa,2))*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(Mc,2)*(-4*std::pow(mc,2) + 2*std::pow(Mc,2) + std::pow(mpb,2) + std::pow(mXa,2))*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mpb,2)*(4*std::pow(mc,2) - 2*std::pow(Mc,2) + std::pow(mpb,2) - std::pow(mXa,2) - 2*s)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mi,2)*(-std::pow(mi,2) + 3*std::pow(mj,2) + std::pow(mpb,2) + std::pow(mXa,2) - 2*s)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*ml*mr*(4*std::pow(mc,2) - 6*std::pow(Mc,2) - s)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mj,2)*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*ml*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*ml*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,2)*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 6*std::pow(mXa,2)*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*(std::pow(mj,2) + std::pow(mXa,2))*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

 	if (n == 1)
		return -12*std::pow(mi,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 18*std::pow(mj,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 18*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*(std::pow(mj,2) + std::pow(mXa,2))*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 18*std::pow(mj,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 18*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 6*(std::pow(mj,2) + std::pow(mXa,2))*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

 	if (n == 2)
		return 12*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccXp::Q_uuchichi(int n, Particle const &p1, Particle const &p2, double s)
{
	int l = p1.get_specific_index(); 
	int r = p2.get_specific_index();
 	double ml = p1.get_specific_index(); 
	double mr = p2.get_specific_index();

  	if (n == 0)
		return -12*std::pow(Mc,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*std::pow(Mc,2)*std::pow(mi,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*std::pow(Mc,2)*std::pow(mj,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mj,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mi,3)*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mi,2)*mj*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*mi*std::pow(mj,2)*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mj,3)*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*std::pow(Mc,2)*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mj,2)*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*ml*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,4)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*(mi - mpb)*(mi + mpb)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mj,2)*(mi - mpb)*(mi + mpb)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*(mi - mpb)*std::pow(mpb,2)*(mi + mpb)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mi,3)*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mi,2)*mj*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*mi*std::pow(mj,2)*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mj,3)*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*std::pow(mpb,2)*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*std::pow(Mc,2)*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mj,2)*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*ml*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,2)*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*mr*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*ml*mr*(4*std::pow(mc,2) + 2*std::pow(Mc,2) - s)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*std::pow(Mc,2)*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,2)*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*ml*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*ml*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mpb,2)*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*(mi - mpb)*(mi + mpb)*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*mr*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*mr*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(Mc,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*std::pow(Mc,2)*std::pow(mi,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*std::pow(Mc,2)*std::pow(mj,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mj,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mi,3)*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mi,2)*mj*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*mi*std::pow(mj,2)*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mj,3)*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*std::pow(Mc,2)*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mj,2)*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*ml*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,4)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mi,2)*(mi - mpb)*(mi + mpb)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mj,2)*(mi - mpb)*(mi + mpb)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*(mi - mpb)*std::pow(mpb,2)*(mi + mpb)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mi,3)*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*std::pow(mi,2)*mj*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*mi*std::pow(mj,2)*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mj,3)*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*std::pow(mpb,2)*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 48*std::pow(Mc,2)*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mj,2)*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*ml*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*std::pow(mpb,2)*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*mj*mr*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*ml*mr*(4*std::pow(mc,2) - 6*std::pow(Mc,2) - s)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*std::pow(Mc,2)*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mj,2)*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*ml*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*ml*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mpb,2)*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*(mi - mpb)*(mi + mpb)*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*mr*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*mr*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

 	if (n == 1)
		return 48*std::pow(Mc,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mj,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*ml*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mpb,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*(mi - mpb)*(mi + mpb)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*mi*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*mr*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mXa,2)*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*s*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 48*std::pow(Mc,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mi,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mj,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*ml*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 24*std::pow(mpb,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*(mi - mpb)*(mi + mpb)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 24*mi*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*mj*mr*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) + 12*std::pow(mXa,2)*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*s*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

 	if (n == 2)
		return -12*aXcc[a][l][i]*aXcc[a][r][i]*std::pow(lPScc[b][i][j],2) - 12*bXcc[a][l][i]*bXcc[a][r][i]*std::pow(lPScc[b][i][j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



// Generated with Mathematica - FeynCalc
int CrossSection_ccXp::Q_order(ParticlePair const&pp)
{
	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return 2;

	return -1;
}

// Generated with Mathematica - FeynCalc (does not include massless vector here)
dcomp CrossSection_ccXp::Q(ParticlePair const&pp, int n, double s, double t)
{
	if (pp.get_type()== INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_ttchichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()== INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_tuchichi(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()== INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::dm && pp.get_second().get_proptype()==Proptype::dm)
		return Q_uuchichi(n, pp.get_first(), pp.get_second(), s);

	return -1;
}
