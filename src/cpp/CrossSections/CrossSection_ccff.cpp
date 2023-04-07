#include "../../headers/CrossSections/CrossSection_ccff.h"

CrossSection_ccff::CrossSection_ccff(DarkSectorModel const &DSMod, Particle const &chii, Particle const &chij, ParticlePair const &couple_SM_ferm)
    : CrossSection(chii, chij, couple_SM_ferm.get_particle(0), couple_SM_ferm.get_particle(1))
{

    //std::cout << chii.get_index() << std::endl;

    int n_S = DSMod.get_n_DS_phis();
    int n_PS = DSMod.get_n_DS_phip();
    int n_X = DSMod.get_n_DS_X();

    lSff.resize(n_S);
    lPSff.resize(n_PS);
    aXff.resize(n_X);
    bXff.resize(n_X);
    lScc.resize(n_S);
    lPScc.resize(n_PS);
    aXcc.resize(n_X);
    bXcc.resize(n_X);

    for (int i = 0; i < n_S; i++)
    {
        lSff[i] = DSMod.get_lambda_S_SM(i, couple_SM_ferm.ind());
        lScc[i] = DSMod.get_lambda_S_DM(i, chii.get_specific_index(), chij.get_specific_index());

        //std::cout << i << " " << lSff[i] << std::endl;
    }
    for (int i = 0; i < n_PS; i++)
    {
        lPSff[i] = DSMod.get_lambda_PS_SM(i, couple_SM_ferm.ind());
        lPScc[i] = DSMod.get_lambda_PS_DM(i, chii.get_specific_index(), chij.get_specific_index());

        //std::cout << i << " " << lSff[i] << " " << lPScc[i] << std::endl;
    }
    for (int i = 0; i < n_X; i++)
    {
        aXff[i] = DSMod.get_a_VEC_SM(i, couple_SM_ferm.ind());
        aXcc[i] = DSMod.get_a_VEC_DM(i, chii.get_specific_index(), chij.get_specific_index());
        bXff[i] = DSMod.get_b_VEC_SM(i, couple_SM_ferm.ind());
        bXcc[i] = DSMod.get_b_VEC_DM(i, chii.get_specific_index(), chij.get_specific_index());
    }

    // First cross section
    std::vector<Particle> as, at, au;

    int const ns = n_S + n_PS + n_X;
    int const nt = 0;
    int const nu = 0;

    as.resize(ns);
    at.resize(nt);
    au.resize(nu);

    int compt = 0;

    for (int i = 0; i < n_S; i++)
    {
        as[compt] = DSMod.get_phis(i);
        compt++;
    }
    for (int i = 0; i < n_PS; i++)
    {
        as[compt] = DSMod.get_phip(i);
        compt++;
    }
    for (int i = 0; i < n_X; i++)
    {
        as[compt] = DSMod.get_X(i);
        compt++;
    }

    mi = _m1;
    mj = _m2;
    mf = _m3;

    mc = (mi + mj) / 2;
    Mc = sqrt(mi * mj);
    dc = (mj - mi) / 2;

    // Remark the code does not take into account the case where mcouple_SM_ferm.get_particle(0) != mcouple_SM_ferm.get_particle(1)

    Initialise(as, at, au, "ccff");

    mX.resize(n_X);

    for (int i = 0; i < n_X; i++)
        mX[i] = DSMod.get_X(i).get_mass();

    //std::cout << mc << std::endl;
    //std::cout << Mc << std::endl;
    //std::cout << dc << std::endl;

    //std::vector<double> mS, mPS, wS, wPS;
    /*mS = phis.get_mass();
    mPS = phip.get_mass();

    wS = phis.get_width();
    wPS = phip.get_width();*/
}


// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccff::Q_ssphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return 4*(4*std::pow(mc,2) - s)*(4*std::pow(mf,2) - s)*lScc[i]*lScc[j]*lSff[i]*lSff[j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccff::Q_ssphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return -4*(4*std::pow(dc,2) - s)*s*lPScc[i]*lPScc[j]*lPSff[i]*lPSff[j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccff::Q_ssXX(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();

  	if (n == 0)
		return 32*std::pow(Mc,2)*std::pow(mf,2)*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] + 16*std::pow(mf,4)*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] + 16*std::pow(mi,2)*std::pow(mj,2)*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] + 16*std::pow(Mc,2)*s*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] - 8*std::pow(mi,2)*s*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] - 8*std::pow(mj,2)*s*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] + 8*std::pow(s,2)*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] - 32*std::pow(Mc,2)*std::pow(mf,2)*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] + 16*std::pow(mf,4)*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] + 16*std::pow(mi,2)*std::pow(mj,2)*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] - 16*std::pow(Mc,2)*s*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] - 8*std::pow(mi,2)*s*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] - 8*std::pow(mj,2)*s*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] + 8*std::pow(s,2)*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] + 128*std::pow(mc,2)*std::pow(mf,2)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] - 64*std::pow(Mc,2)*std::pow(mf,2)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] + 16*std::pow(mf,4)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] - 96*std::pow(mf,2)*mi*mj*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] + 16*std::pow(mi,2)*std::pow(mj,2)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] + 16*std::pow(Mc,2)*s*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] - 32*std::pow(mf,2)*s*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] - 8*std::pow(mi,2)*s*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] - 8*std::pow(mj,2)*s*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] + 8*std::pow(s,2)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] + 128*std::pow(mc,2)*std::pow(mf,2)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] - 64*std::pow(Mc,2)*std::pow(mf,2)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] + 16*std::pow(mf,4)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] + 96*std::pow(mf,2)*mi*mj*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] + 16*std::pow(mi,2)*std::pow(mj,2)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] - 16*std::pow(Mc,2)*s*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] - 32*std::pow(mf,2)*s*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] - 8*std::pow(mi,2)*s*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] - 8*std::pow(mj,2)*s*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] + 8*std::pow(s,2)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] - 32*std::pow(mc,2)*s*aXcc[j]*bXcc[i]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) + 16*std::pow(Mc,2)*s*aXcc[j]*bXcc[i]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) - 16*std::pow(mf,2)*s*aXcc[j]*bXcc[i]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) + 8*std::pow(s,2)*aXcc[j]*bXcc[i]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) - 32*std::pow(mc,2)*s*aXcc[i]*bXcc[j]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) + 16*std::pow(Mc,2)*s*aXcc[i]*bXcc[j]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) - 16*std::pow(mf,2)*s*aXcc[i]*bXcc[j]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) + 8*std::pow(s,2)*aXcc[i]*bXcc[j]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) + (64*std::pow(dc,2)*std::pow(mf,2)*(4*std::pow(mc,2) - s)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j])/std::pow(mX[i],2) + (64*std::pow(mc,2)*std::pow(mf,2)*(4*std::pow(dc,2) - s)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j])/std::pow(mX[i],2) - (64*std::pow(dc,2)*std::pow(mf,2)*(4*std::pow(mc,2) - s)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j]*(s - std::pow(mX[i],2)))/(std::pow(mX[i],2)*std::pow(mX[j],2)) - (64*std::pow(mc,2)*std::pow(mf,2)*(4*std::pow(dc,2) - s)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j]*(s - std::pow(mX[i],2)))/(std::pow(mX[i],2)*std::pow(mX[j],2));

 	if (n == 1)
		return -32*std::pow(mf,2)*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] - 16*std::pow(mi,2)*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] - 16*std::pow(mj,2)*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] + 16*s*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] - 32*std::pow(mf,2)*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] - 16*std::pow(mi,2)*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] - 16*std::pow(mj,2)*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] + 16*s*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] - 32*std::pow(mf,2)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] - 16*std::pow(mi,2)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] - 16*std::pow(mj,2)*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] + 16*s*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] - 32*std::pow(mf,2)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] - 16*std::pow(mi,2)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] - 16*std::pow(mj,2)*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] + 16*s*bXcc[i]*bXcc[j]*bXff[i]*bXff[j] + 16*s*aXcc[j]*bXcc[i]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]) + 16*s*aXcc[i]*bXcc[j]*(aXff[j]*bXff[i] + aXff[i]*bXff[j]);

 	if (n == 2)
		return 16*aXcc[i]*aXcc[j]*aXff[i]*aXff[j] + 16*aXff[i]*aXff[j]*bXcc[i]*bXcc[j] + 16*aXcc[i]*aXcc[j]*bXff[i]*bXff[j] + 16*bXcc[i]*bXcc[j]*bXff[i]*bXff[j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccff::Q_ssphisphip(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccff::Q_ssphisX(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();

  	if (n == 0)
		return -64*std::pow(mc,3)*mf*aXcc[j]*aXff[j]*lScc[i]*lSff[i] + 32*mc*std::pow(Mc,2)*mf*aXcc[j]*aXff[j]*lScc[i]*lSff[i] - 32*mc*std::pow(mf,3)*aXcc[j]*aXff[j]*lScc[i]*lSff[i] + 16*mc*mf*s*aXcc[j]*aXff[j]*lScc[i]*lSff[i];

 	if (n == 1)
		return 32*mc*mf*aXcc[j]*aXff[j]*lScc[i]*lSff[i];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_ccff::Q_ssphipX(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();

  	if (n == 0)
		return (-16*mc*mf*(4*std::pow(dc,2) - s)*bXcc[j]*bXff[j]*lPScc[i]*lPSff[i]*(s - std::pow(mX[j],2)))/std::pow(mX[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



// Generated with Mathematica - FeynCalc
int CrossSection_ccff::Q_order(ParticlePair const&pp)
{
	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 0;

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 0;

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
		return 2;

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return -1;

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
		return 1;

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
		return 0;

	return -1;
}

// Generated with Mathematica - FeynCalc (does not include massless vector here)
dcomp CrossSection_ccff::Q(ParticlePair const&pp, int n, double s, double t)
{
	if (pp.get_type()== INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_ssphisphis(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_ssphipphip(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
		return Q_ssXX(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_ssphisphip(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
		return Q_ssphisX(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
		return Q_ssphipX(n, pp.get_first(), pp.get_second(), s);

	return -1;
}



/*


*/
// Additional function

/*
std::vector<double> CrossSection_ccff::AppSigmaVChanSPropSPS()
{
    double mDM, mF;

    if (_m1 != _m2 || _m3 != _m4)
    {
        std::cout << "FATAL ERROR : CrossSection_ccff::AppSigmaVChanSPropSPS() . only works for one DM species and couples (ferm, bar{ferm})" << std::endl;
        exit(0);
    }
    else
    {
        mDM = _m1;
        mF = _m3;
    }

    std::vector<double> coeff;
    coeff.resize(2); // Compute only two coeff

    double mPS = 100, mS = 100, wS = 10, wPS = 10;

    double rFDM = pow(mF / mDM, 2);
    double rPSDM = pow(mPS / (2 * mDM), 2);
    double rSDM = pow(mS / (2 * mDM), 2);

    double uFDM = 1 - pow(mF / mDM, 2);
    double uPSDM = 1 - pow(mPS / (2 * mDM), 2);
    double uSDM = 1 - pow(mS / (2 * mDM), 2);
    double mDM2 = mDM * mDM;

    double rgs = mS * wS / (4 * mDM2);
    double rgp = mPS * wPS / (4 * mDM2);

    double couplPS2 = pow(lPScc[0] * lPSff[0], 2);
    double couplS2 = pow(lScc[0] * lSff[0], 2);

    coeff[0] = (couplPS2 / (32 * PI * mDM2)) * pow(uFDM, 0.5) * pow(pow(uPSDM, 2) + rgp * rgp, -1);
    coeff[1] = (-3 / (2 * mDM2)) * ((couplPS2 / (16 * PI)) * pow(uFDM, 0.5) * pow(pow(uPSDM, 2) + rgp * rgp, -1) -
                                    (couplS2 / (32 * PI)) * pow(uFDM, 1.5) * pow(pow(uSDM, 2) + rgs * rgs, -1) -
                                    (couplPS2 / (64 * PI)) * pow(uFDM, -0.5) * rFDM * pow(pow(uPSDM, 2) + rgp * rgp, -1) +
                                    (couplPS2 / (16 * PI)) * pow(uFDM, 0.5) * (uPSDM * rPSDM - rgp * rgp) * pow(pow(uPSDM, 2) + rgp * rgp, -2));

    return coeff;
}
*/