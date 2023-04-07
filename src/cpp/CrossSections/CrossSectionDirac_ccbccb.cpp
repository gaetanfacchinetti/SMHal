#include "../../headers/CrossSections/CrossSectionDirac_ccbccb.h"

CrossSectionDirac_ccbccb::CrossSectionDirac_ccbccb(DarkSectorModel const &DSMod, Particle const &chi, Particle const &chib)
    : CrossSection(chi, chib, chi, chib)
{

    //std::cout << chii.get_index() << std::endl;

    if (chi.get_fermiontype() != Fermiontype::dirac)
    {
        std::cout << "FATAL ERROR : Trying to evaluate Dirac self-interaction cross-section of majorana fermion in " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    int n_S = DSMod.get_n_DS_phis();
    int n_PS = DSMod.get_n_DS_phip();
    int n_X = DSMod.get_n_DS_X();

    lScc.resize(n_S);
    lPScc.resize(n_PS);
    aXcc.resize(n_X);
    bXcc.resize(n_X);

    chi_index = chi.get_specific_index();
    mc = chi.get_mass();

    for (int i = 0; i < n_S; i++)
        lScc[i] = DSMod.get_lambda_S_DM(i, chi_index, chi_index);

    for (int i = 0; i < n_PS; i++)
        lPScc[i] = DSMod.get_lambda_PS_DM(i, chi_index, chi_index);

    for (int i = 0; i < n_X; i++)
    {
        if (chi.get_fermiontype() == Fermiontype::majorana)
            aXcc[i] = 0;
        else
            aXcc[i] = DSMod.get_a_VEC_DM(i, chi_index, chi_index);

        bXcc[i] = DSMod.get_b_VEC_DM(i, chi_index, chi_index);
    }

    // First cross section
    std::vector<Particle> as, at, au;

    int const ns = n_S + n_PS + n_X;
    int const nt = n_S + n_PS + n_X;
    int nu = 0;


    as.resize(ns);
    at.resize(nt);
    au.resize(nu);

    int compt = 0;

    for (int i = 0; i < n_S; i++)
    {
        at[compt] = DSMod.get_phis(i);
        as[compt] = DSMod.get_phis(i);

        compt++;
    }

    for (int i = 0; i < n_PS; i++)
    {
        at[compt] = DSMod.get_phip(i);
        as[compt] = DSMod.get_phip(i);

        compt++;
    }
    for (int i = 0; i < n_X; i++)
    {

        at[compt] = DSMod.get_X(i);
        as[compt] = DSMod.get_X(i);

        compt++;
    }

    // Remark the code does not take into account the case where mcouple_SM_ferm.get_particle(0) != mcouple_SM_ferm.get_particle(1)

    Initialise(as, at, au, "ccbccb");

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
dcomp CrossSectionDirac_ccbccb::Q_ssphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return 4*std::pow(-4*std::pow(mc,2) + s,2)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return 4*std::pow(s,2)*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
 	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

  	if (n == 0)
		return 64*std::pow(mc,4)*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) + 8*std::pow(s,2)*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) - 32*std::pow(mc,2)*s*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) + 8*std::pow(s,2)*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) - 128*std::pow(mc,2)*s*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 32*std::pow(s,2)*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] - 32*std::pow(mc,2)*s*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(s,2)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) + 192*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) - 64*std::pow(mc,2)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 1)
		return -64*std::pow(mc,2)*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) + 16*s*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) + 16*s*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) - 128*std::pow(mc,2)*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 64*s*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 16*s*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) - 64*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 16*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 2)
		return 16*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) + 32*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 16*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssXXgppg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
 	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

	if (n == 1)
		return 0;

	if (n == 2)
		return 0;

 	if (n == 0)
		return (-64*std::pow(mc,4)*(std::pow(mXi,2) + std::pow(mXj,2))*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssXXpg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXi == 0) 
		return 0;

	if (n == 1)
		return 0;

	if (n == 2)
		return 0;

 	if (n == 0)
		return (-64*std::pow(mc,4)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssXXpp(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

	if (n == 1)
		return 0;

	if (n == 2)
		return 0;

 	if (n == 0)
		return (64*std::pow(mc,4)*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (n == 0)
		return -((64*std::pow(mc,4)*std::pow(aXcc[j],2) - 16*std::pow(mc,2)*s*std::pow(aXcc[j],2))*std::pow(lScc[i],2));

 	if (n == 1)
		return 32*std::pow(mc,2)*std::pow(aXcc[j],2)*std::pow(lScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssphisXp(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

	return 0;

	if (n == 0)
		return 0;

	if (n == 1)
		return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssphipXg(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (n == 0)
		return -16*std::pow(mc,2)*s*std::pow(bXcc[j],2)*std::pow(lPScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ssphipXp(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

	if (n == 0)
		return 0;

 	if (n == 0)
		return (16*std::pow(mc,2)*std::pow(s,2)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return 64*std::pow(mc,4)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

 	if (n == 1)
		return -32*std::pow(mc,2)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

 	if (n == 2)
		return 4*std::pow(lScc[i],2)*std::pow(lScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return 0;

 	if (n == 1)
		return 0;

 	if (n == 2)
		return 4*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

  	if (n == 0)
		return 16*std::pow(-2*std::pow(mc,2) + s,2)*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) - 128*std::pow(mc,2)*s*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 32*std::pow(s,2)*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 192*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) - 64*std::pow(mc,2)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 16*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 1)
		return 16*s*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) - 32*std::pow(mc,2)*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) + 16*s*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) - 128*std::pow(mc,2)*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 64*s*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] - 32*std::pow(mc,2)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) + 16*s*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) - 64*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 16*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 2)
		return 8*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) + 8*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) + 32*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 8*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttXXgppg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

	if (n == 2)
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return (-64*std::pow(mc,4)*(std::pow(mXi,2) + std::pow(mXj,2))*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttXXpg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXi == 0) 
		return 0;

	if (n == 2)
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return (-64*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttXXpp(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return 0;

 	if (n == 2)
		return (64*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

  	if (n == 0)
		return -((64*std::pow(mc,4)*std::pow(aXcc[j],2) - 32*std::pow(mc,2)*s*std::pow(aXcc[j],2))*std::pow(lScc[i],2));

 	if (n == 1)
		return 16*std::pow(mc,2)*std::pow(aXcc[j],2)*std::pow(lScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttphisXp(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

	return 0;

	if (n == 0)
		return 0;

	if (n == 1)
		return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttphipXg(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (n == 2)
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return -16*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_ttphipXp(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return 0;

 	if (n == 2)
		return (16*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return (-32*std::pow(mc,4) + 8*std::pow(mc,2)*s)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

 	if (n == 1)
		return (8*std::pow(mc,2) + 2*s)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return 0;

 	if (n == 1)
		return 2*s*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

  	if (n == 0)
		return -32*std::pow(mc,4)*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) + 8*std::pow(s,2)*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) + 96*std::pow(mc,4)*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) - 64*std::pow(mc,2)*s*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) + 8*std::pow(s,2)*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) - 128*std::pow(mc,2)*s*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 32*std::pow(s,2)*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 96*std::pow(mc,4)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) - 32*std::pow(mc,2)*s*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(s,2)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) + 96*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) - 32*std::pow(mc,2)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 1)
		return 16*s*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) - 32*std::pow(mc,2)*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) + 16*s*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) - 128*std::pow(mc,2)*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 64*s*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] - 64*std::pow(mc,2)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) + 16*s*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) - 32*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 16*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 2)
		return 8*std::pow(aXcc[i],2)*std::pow(aXcc[j],2) + 8*std::pow(aXcc[j],2)*std::pow(bXcc[i],2) + 32*aXcc[i]*aXcc[j]*bXcc[i]*bXcc[j] + 8*std::pow(aXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stXXgppg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

 	if (n == 0)
		return (16*std::pow(mc,2)*s*std::pow(bXcc[i],2)*((-2*std::pow(mc,2) + s)*std::pow(aXcc[j],2) - (2*std::pow(mc,2) + s)*std::pow(bXcc[j],2)))/std::pow(mXi,2);

 	if (n == 1)
		return (-32*std::pow(mc,4)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (32*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

 	if (n == 2)
		return (16*std::pow(mc,2)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (16*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stXXpg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXi == 0) 
		return 0;

 	if (n == 0)
		return (16*std::pow(mc,2)*s*std::pow(bXcc[i],2)*((-2*std::pow(mc,2) + s)*std::pow(aXcc[j],2) - (2*std::pow(mc,2) + s)*std::pow(bXcc[j],2)))/std::pow(mXi,2);

 	if (n == 1)
		return (-32*std::pow(mc,4)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (32*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

 	if (n == 2)
		return (16*std::pow(mc,2)*std::pow(aXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (16*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stXXpp(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

	if (n == 2)
			return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return (32*std::pow(mc,4)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return 0;

 	if (n == 1)
		return -2*(-4*std::pow(mc,2) + s)*std::pow(lPScc[j],2)*std::pow(lScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

  	if (n == 0)
		return -((-32*std::pow(mc,4)*std::pow(aXcc[j],2) + 24*std::pow(mc,2)*s*std::pow(aXcc[j],2) - 4*std::pow(s,2)*std::pow(aXcc[j],2) + 96*std::pow(mc,4)*std::pow(bXcc[j],2) - 40*std::pow(mc,2)*s*std::pow(bXcc[j],2) + 4*std::pow(s,2)*std::pow(bXcc[j],2))*std::pow(lScc[i],2));

 	if (n == 1)
		return -((-16*std::pow(mc,2)*std::pow(aXcc[j],2) - 16*std::pow(mc,2)*std::pow(bXcc[j],2))*std::pow(lScc[i],2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stphisXp(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return (-8*std::pow(mc,2)*(-4*std::pow(mc,2) + s)*std::pow(bXcc[j],2)*std::pow(lScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stphipXg(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (n == 1)
		return 0;

 	if (n == 0)
		return (4*s*(-2*std::pow(mc,2) + s)*std::pow(aXcc[j],2) - 4*s*(2*std::pow(mc,2) + s)*std::pow(bXcc[j],2))*std::pow(lPScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q_stphipXp(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (mXj == 0) 
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return (8*std::pow(mc,2)*s*std::pow(bXcc[j],2)*std::pow(lPScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



int CrossSectionDirac_ccbccb::Q_order(ParticlePair const&pp)
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


	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 2;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 2;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
		return 2;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return -1;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
		return 1;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
		return 2;


	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 1;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 1;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
		return 2;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 1;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
		return 1;

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
		return 1;


	return -1;

}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionDirac_ccbccb::Q(ParticlePair const & pp,int n, double s, double t)
{

	double res=0;

	// --------------------------------------
	// SS terms
	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_ssphisphis(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_ssphipphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()>0)
			return Q_ssXXgg(n,pp.get_first(),pp.get_second(),s)+Q_ssXXgppg(n,pp.get_first(),pp.get_second(),s)+Q_ssXXpp(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()==0)
			return  Q_ssXXgg(n,pp.get_second(),pp.get_first(),s)+Q_ssXXpg(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()==0 && pp.get_second().get_mass()>0)
			return Q_ssXXgg(n,pp.get_first(),pp.get_second(),s)+Q_ssXXpg(n,pp.get_second(),pp.get_first(),s);
		else
			return Q_ssXXgg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_ssphisphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_ssphisXg(n,pp.get_first(),pp.get_second(),s)+Q_ssphisXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_ssphisXg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::SS && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_ssphipXg(n,pp.get_first(),pp.get_second(),s)+Q_ssphipXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_ssphisXg(n,pp.get_first(),pp.get_second(),s);
	}
	// --------------------------------------

	// --------------------------------------
	// TT terms
	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_ttphisphis(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_ttphipphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()>0)
			return Q_ttXXgg(n,pp.get_first(),pp.get_second(),s)+Q_ttXXgppg(n,pp.get_first(),pp.get_second(),s)+Q_ttXXpp(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()==0)
			return  Q_ttXXgg(n,pp.get_second(),pp.get_first(),s)+Q_ttXXpg(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()==0 && pp.get_second().get_mass()>0)
			return Q_ttXXgg(n,pp.get_first(),pp.get_second(),s)+Q_ttXXpg(n,pp.get_second(),pp.get_first(),s);
		else
			return Q_ttXXgg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_ttphisphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_ttphisXg(n,pp.get_first(),pp.get_second(),s)+Q_ttphisXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_ttphisXg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_ttphipXg(n,pp.get_first(),pp.get_second(),s)+Q_ttphipXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_ttphisXg(n,pp.get_first(),pp.get_second(),s);
	}
	// --------------------------------------

	// --------------------------------------
	// ST terms
	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_stphisphis(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_stphipphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()>0)
			return Q_stXXgg(n,pp.get_first(),pp.get_second(),s)+Q_stXXgppg(n,pp.get_first(),pp.get_second(),s)+Q_stXXpp(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()==0)
			return  Q_stXXgg(n,pp.get_second(),pp.get_first(),s)+Q_stXXpg(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()==0 && pp.get_second().get_mass()>0)
			return Q_stXXgg(n,pp.get_first(),pp.get_second(),s)+Q_stXXpg(n,pp.get_second(),pp.get_first(),s);
		else
			return Q_stXXgg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_stphisphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_stphisXg(n,pp.get_first(),pp.get_second(),s)+Q_stphisXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_stphisXg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::ST && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_stphipXg(n,pp.get_first(),pp.get_second(),s)+Q_stphipXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_stphisXg(n,pp.get_first(),pp.get_second(),s);
	}

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

