#include "../../headers/CrossSections/CrossSection_cfcf.h"

CrossSection_cfcf::CrossSection_cfcf(DarkSectorModel const &DSMod, Particle const &chii, ParticlePair const &couple_SM_ferm) : CrossSection(chii, couple_SM_ferm.get_particle(0), chii, couple_SM_ferm.get_particle(0))
{

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
		lScc[i] = DSMod.get_lambda_S_DM(i, chii.get_specific_index(), chii.get_specific_index());

		//std::cout << i << " " << lSff[i] << std::endl;
	}
	for (int i = 0; i < n_PS; i++)
	{
		lPSff[i] = DSMod.get_lambda_PS_SM(i, couple_SM_ferm.ind());
		lPScc[i] = DSMod.get_lambda_PS_DM(i, chii.get_specific_index(), chii.get_specific_index());

		//std::cout << i << " " << lSff[i] << " " << lPScc[i] << std::endl;
	}
	for (int i = 0; i < n_X; i++)
	{
		aXff[i] = DSMod.get_a_VEC_SM(i, couple_SM_ferm.ind());
		aXcc[i] = DSMod.get_a_VEC_DM(i, chii.get_specific_index(), chii.get_specific_index());
		bXff[i] = DSMod.get_b_VEC_SM(i, couple_SM_ferm.ind());
		bXcc[i] = DSMod.get_b_VEC_DM(i, chii.get_specific_index(), chii.get_specific_index());
	}

	// First cross section
	std::vector<Particle> as, at, au;

	int const ns = 0;
	int const nt = n_S + n_PS + n_X;
	int const nu = 0;

	as.resize(ns);
	at.resize(nt);
	au.resize(nu);

	int compt = 0;

	for (int i = 0; i < n_S; i++)
	{
		at[compt] = DSMod.get_phis(i);
		compt++;
	}
	for (int i = 0; i < n_PS; i++)
	{
		at[compt] = DSMod.get_phip(i);
		compt++;
	}
	for (int i = 0; i < n_X; i++)
	{
		at[compt] = DSMod.get_X(i);
		compt++;
	}

	mi = _m1;
	mf = _m2;

	// Remark the code does not take into account the case where mferma != mfermb

	//std::cout << "couplings : " << lSff[0] << std::endl;

	Initialise(as, at, au, "cfcf");
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cfcf::Q_ttphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	int i = s1.get_specific_index();
	int j = s2.get_specific_index();

	//std::cout << mf << " " << mi << " " << lScc[i]*lScc[j]*lSff[i]*lSff[j] << std::endl;

	if (n == 0)
		return 64 * std::pow(mf, 2) * std::pow(mi, 2) * lScc[i] * lScc[j] * lSff[i] * lSff[j];

	if (n == 1)
		return (-16 * std::pow(mf, 2) - 16 * std::pow(mi, 2)) * lScc[i] * lScc[j] * lSff[i] * lSff[j];

	if (n == 2)
		return 4 * lScc[i] * lScc[j] * lSff[i] * lSff[j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0);
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cfcf::Q_ttphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index();
	int j = ps2.get_specific_index();

	if (n == 0)
		return 0;

	if (n == 1)
		return 0;

	if (n == 2)
		return 4 * lPScc[i] * lPScc[j] * lPSff[i] * lPSff[j];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0);
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cfcf::Q_ttXX(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index();
	int j = X2.get_specific_index();
	double mXi = X1.get_mass();
	double mXj = X2.get_mass();

	if (n == 0)
		return 16 * std::pow(std::pow(mf, 2) + std::pow(mi, 2) - s, 2) * aXcc[i] * aXcc[j] * aXff[i] * aXff[j] + 16 * std::pow(mf, 4) * aXff[i] * aXff[j] * bXcc[i] * bXcc[j] + 16 * std::pow(std::pow(mi, 2) - s, 2) * aXff[i] * aXff[j] * bXcc[i] * bXcc[j] - 32 * std::pow(mf, 2) * (std::pow(mi, 2) + s) * aXff[i] * aXff[j] * bXcc[i] * bXcc[j] + 16 * std::pow(mf, 4) * aXcc[i] * aXcc[j] * bXff[i] * bXff[j] - 32 * std::pow(mf, 2) * std::pow(mi, 2) * aXcc[i] * aXcc[j] * bXff[i] * bXff[j] + 16 * std::pow(std::pow(mi, 2) - s, 2) * aXcc[i] * aXcc[j] * bXff[i] * bXff[j] - 32 * std::pow(mf, 2) * s * aXcc[i] * aXcc[j] * bXff[i] * bXff[j] + 16 * (std::pow(mf, 4) + std::pow(std::pow(mi, 2) - s, 2) + 2 * std::pow(mf, 2) * (5 * std::pow(mi, 2) - s)) * bXcc[i] * bXcc[j] * bXff[i] * bXff[j];

	if (n == 1)
		return 16 * s * aXcc[i] * aXcc[j] * aXff[i] * aXff[j] + 16 * (-2 * std::pow(mi, 2) + s) * aXff[i] * aXff[j] * bXcc[i] * bXcc[j] - 32 * std::pow(mf, 2) * aXcc[i] * aXcc[j] * bXff[i] * bXff[j] + 16 * s * aXcc[i] * aXcc[j] * bXff[i] * bXff[j] + 16 * (-2 * (std::pow(mf, 2) + std::pow(mi, 2)) + s) * bXcc[i] * bXcc[j] * bXff[i] * bXff[j] + 16 * (std::pow(mf, 2) + std::pow(mi, 2) - s) * aXcc[j] * bXcc[i] * (aXff[j] * bXff[i] + aXff[i] * bXff[j]) + 16 * (std::pow(mf, 2) + std::pow(mi, 2) - s) * aXcc[i] * bXcc[j] * (aXff[j] * bXff[i] + aXff[i] * bXff[j]) - (128 * std::pow(mf, 2) * std::pow(mi, 2) * bXcc[i] * bXcc[j] * bXff[i] * bXff[j]) / std::pow(mX[i], 2) - (64 * std::pow(mf, 2) * std::pow(mi, 2) * bXcc[i] * bXcc[j] * bXff[i] * bXff[j]) / std::pow(mX[j], 2);

	if (n == 2)
		return 8 * aXcc[i] * aXcc[j] * aXff[i] * aXff[j] + 8 * aXff[i] * aXff[j] * bXcc[i] * bXcc[j] + 8 * aXcc[i] * aXcc[j] * bXff[i] * bXff[j] + 8 * bXcc[i] * bXcc[j] * bXff[i] * bXff[j] - 8 * aXcc[j] * bXcc[i] * (aXff[j] * bXff[i] + aXff[i] * bXff[j]) - 8 * aXcc[i] * bXcc[j] * (aXff[j] * bXff[i] + aXff[i] * bXff[j]) + (64 * std::pow(mf, 2) * std::pow(mi, 2) * bXcc[i] * bXcc[j] * bXff[i] * bXff[j]) / (std::pow(mX[i], 2) * std::pow(mX[j], 2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0);
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cfcf::Q_ttphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index();
	int j = ps2.get_specific_index();

	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0);
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cfcf::Q_ttphisX(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index();
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (n == 0)
		return -((-32 * std::pow(mf, 3) * mi * aXcc[j] * aXff[j] - 32 * mf * std::pow(mi, 3) * aXcc[j] * aXff[j] + 32 * mf * mi * s * aXcc[j] * aXff[j]) * lScc[i] * lSff[i]);

	if (n == 1)
		return -16 * mf * mi * aXcc[j] * aXff[j] * lScc[i] * lSff[i];

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0);
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cfcf::Q_ttphipX(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index();
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (n == 0)
		return 0;

	if (n == 1)
		return 16 * mf * mi * bXcc[j] * bXff[j] * lScc[i] * lSff[i];

	if (n == 2)
		return (-16 * mf * mi * bXcc[j] * bXff[j] * lScc[i] * lSff[i]) / std::pow(mX[j], 2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0);
}

// Generated with Mathematica - FeynCalc
int CrossSection_cfcf::Q_order(ParticlePair const &pp)
{
	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::scalar && pp.get_second().get_proptype() == Proptype::scalar)
		return 2;

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::pseudoscalar && pp.get_second().get_proptype() == Proptype::pseudoscalar)
		return 2;

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::vector && pp.get_second().get_proptype() == Proptype::vector)
		return 2;

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::scalar && pp.get_second().get_proptype() == Proptype::pseudoscalar)
		return -1;

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::scalar && pp.get_second().get_proptype() == Proptype::vector)
		return 1;

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::pseudoscalar && pp.get_second().get_proptype() == Proptype::vector)
		return 2;

	return -1;
}

// Generated with Mathematica - FeynCalc (does not include massless vector here)
dcomp CrossSection_cfcf::Q(ParticlePair const &pp, int n, double s, double t)
{
	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::scalar && pp.get_second().get_proptype() == Proptype::scalar)
		return Q_ttphisphis(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::pseudoscalar && pp.get_second().get_proptype() == Proptype::pseudoscalar)
		return Q_ttphipphip(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::vector && pp.get_second().get_proptype() == Proptype::vector)
		return Q_ttXX(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::scalar && pp.get_second().get_proptype() == Proptype::pseudoscalar)
		return Q_ttphisphip(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::scalar && pp.get_second().get_proptype() == Proptype::vector)
		return Q_ttphisX(n, pp.get_first(), pp.get_second(), s);

	if (pp.get_type() == INT_TYPE::TT && pp.get_first().get_proptype() == Proptype::pseudoscalar && pp.get_second().get_proptype() == Proptype::vector)
		return Q_ttphipX(n, pp.get_first(), pp.get_second(), s);

	return -1;
}
