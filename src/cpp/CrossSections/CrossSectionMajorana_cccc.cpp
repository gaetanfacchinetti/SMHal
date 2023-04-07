#include "../../headers/CrossSections/CrossSectionMajorana_cccc.h"

CrossSectionMajorana_cccc::CrossSectionMajorana_cccc(DarkSectorModel const &DSMod, Particle const &chi)
    : CrossSection(chi, chi, chi, chi)
{

    //std::cout << chii.get_index() << std::endl;

    if (chi.get_fermiontype() != Fermiontype::majorana)
    {
        std::cout << "FATAL ERROR : Trying to evaluate Majorana self-interaction cross-section of dirac fermion in " << __PRETTY_FUNCTION__ << std::endl;
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

    if (chi.get_fermiontype() == Fermiontype::majorana)
        nu = n_S + n_PS + n_X; // If Majorana !!
    else
        nu = 0;

    as.resize(ns);
    at.resize(nt);
    au.resize(nu);

    int compt = 0;

    for (int i = 0; i < n_S; i++)
    {
        as[compt] = DSMod.get_phis(i);
        at[compt] = DSMod.get_phis(i);

        if (chi.get_fermiontype() == Fermiontype::majorana)
            au[compt] = DSMod.get_phis(i);

        compt++;
    }

    for (int i = 0; i < n_PS; i++)
    {
        as[compt] = DSMod.get_phip(i);
        at[compt] = DSMod.get_phip(i);

        if (chi.get_fermiontype() == Fermiontype::majorana)
            au[compt] = DSMod.get_phip(i);

        compt++;
    }
    for (int i = 0; i < n_X; i++)
    {
        as[compt] = DSMod.get_X(i);
        at[compt] = DSMod.get_X(i);

        if (chi.get_fermiontype() == Fermiontype::majorana)
            au[compt] = DSMod.get_X(i);

        compt++;
    }

    // Remark the code does not take into account the case where mcouple_SM_ferm.get_particle(0) != mcouple_SM_ferm.get_particle(1)

    Initialise(as, at, au, "cccc");

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
dcomp CrossSectionMajorana_cccc::Q_ssphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	//return 0;
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return 4*std::pow(-4*std::pow(mc,2) + s,2)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ssphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return 4*std::pow(s,2)*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ssXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
 	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

  	if (n == 0)
		return 192*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) - 64*std::pow(mc,2)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 1)
		return -64*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 16*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 2)
		return 16*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ssXXgppg(int n, Particle const &X1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ssXXpg(int n, Particle const &X1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ssXXpp(int n, Particle const &X1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ssphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ssphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ssphisXp(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ssphipXg(int n, Particle const &ps1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ssphipXp(int n, Particle const &ps1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ttphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	//return 0;
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
dcomp CrossSectionMajorana_cccc::Q_ttphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ttXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

  	if (n == 0)
		return 192*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) - 64*std::pow(mc,2)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 16*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 1)
		return -64*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 16*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 2)
		return 8*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ttXXgppg(int n, Particle const &X1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ttXXpg(int n, Particle const &X1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ttXXpp(int n, Particle const &X1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ttphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ttphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ttphisXp(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_ttphipXg(int n, Particle const &ps1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_ttphipXp(int n, Particle const &ps1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_uuphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	//return 0;
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return 4*std::pow(s,2)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

 	if (n == 1)
		return 8*s*std::pow(lScc[i],2)*std::pow(lScc[j],2);

 	if (n == 2)
		return 4*std::pow(lScc[i],2)*std::pow(lScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return (64*std::pow(mc,4) - 32*std::pow(mc,2)*s + 4*std::pow(s,2))*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

 	if (n == 1)
		return (-32*std::pow(mc,2) + 8*s)*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

 	if (n == 2)
		return 4*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

  	if (n == 0)
		return 64*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 1)
		return 0;

 	if (n == 2)
		return 8*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuXXgppg(int n, Particle const &X1, Particle const &X2, double s)
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
		return (-256*std::pow(mc,6)*(std::pow(mXi,2) + std::pow(mXj,2))*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2)) + (64*std::pow(mc,4)*(std::pow(mXi,2) + std::pow(mXj,2))*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

 	if (n == 1)
		return (64*std::pow(mc,4)*(std::pow(mXi,2) + std::pow(mXj,2))*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuXXpg(int n, Particle const &X1, Particle const &X2, double s)
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
		return (-256*std::pow(mc,6)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2) + (64*std::pow(mc,4)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

 	if (n == 1)
		return (64*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuXXpp(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

 	if (n == 0)
		return (1024*std::pow(mc,8)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2)) - (512*std::pow(mc,6)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2)) + (64*std::pow(mc,4)*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

 	if (n == 1)
		return (-512*std::pow(mc,6)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2)) + (128*std::pow(mc,4)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

 	if (n == 2)
		return (64*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuphisXp(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuphipXg(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (n == 2)
		return 0;

 	if (n == 0)
		return (-64*std::pow(mc,4)*std::pow(bXcc[j],2) + 16*std::pow(mc,2)*s*std::pow(bXcc[j],2))*std::pow(lPScc[i],2);

 	if (n == 1)
		return 16*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_uuphipXp(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (mXj == 0) 
		return 0;

 	if (n == 0)
		return ((256*std::pow(mc,6)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (128*std::pow(mc,4)*s*std::pow(bXcc[j],2))/std::pow(mXj,2) + (16*std::pow(mc,2)*std::pow(s,2)*std::pow(bXcc[j],2))/std::pow(mXj,2))*std::pow(lPScc[i],2);

 	if (n == 1)
		return ((-128*std::pow(mc,4)*std::pow(bXcc[j],2))/std::pow(mXj,2) + (32*std::pow(mc,2)*s*std::pow(bXcc[j],2))/std::pow(mXj,2))*std::pow(lPScc[i],2);

 	if (n == 2)
		return (16*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_stphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	//return 0;
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
dcomp CrossSectionMajorana_cccc::Q_stphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_stXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

  	if (n == 0)
		return 96*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) - 32*std::pow(mc,2)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 8*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 1)
		return -32*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2) + 16*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 2)
		return 8*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_stXXgppg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

 	if (n == 0)
		return (-16*std::pow(mc,2)*s*(2*std::pow(mc,2) + s)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

 	if (n == 1)
		return (-32*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

 	if (n == 2)
		return (-16*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_stXXpg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXi == 0) 
		return 0;

 	if (n == 0)
		return (-16*std::pow(mc,2)*s*(2*std::pow(mc,2) + s)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

 	if (n == 1)
		return (-32*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

 	if (n == 2)
		return (-16*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_stXXpp(int n, Particle const &X1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_stphisphip(int n, Particle const &s1, Particle const &ps2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_stphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

  	if (n == 0)
		return -((96*std::pow(mc,4)*std::pow(bXcc[j],2) - 40*std::pow(mc,2)*s*std::pow(bXcc[j],2) + 4*std::pow(s,2)*std::pow(bXcc[j],2))*std::pow(lScc[i],2));

 	if (n == 1)
		return 16*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_stphisXp(int n, Particle const &s1, Particle const &X2, double s)
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
dcomp CrossSectionMajorana_cccc::Q_stphipXg(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (n == 1)
		return 0;

 	if (n == 0)
		return -4*s*(2*std::pow(mc,2) + s)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_stphipXp(int n, Particle const &ps1, Particle const &X2, double s)
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

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return (8*std::pow(mc,2)*s - 2*std::pow(s,2))*std::pow(lScc[i],2)*std::pow(lScc[j],2);

 	if (n == 1)
		return (-8*std::pow(mc,2) - 2*s)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return (8*std::pow(mc,2)*s - 2*std::pow(s,2))*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

 	if (n == 1)
		return -2*s*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

  	if (n == 0)
		return 96*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 1)
		return -32*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

 	if (n == 2)
		return 8*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suXXgppg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

 	if (n == 0)
		return (-384*std::pow(mc,6)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) + (160*std::pow(mc,4)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (16*std::pow(mc,2)*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (16*std::pow(mc,2)*s*(2*std::pow(mc,2) + s)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

 	if (n == 1)
		return (160*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (32*std::pow(mc,2)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

 	if (n == 2)
		return (-16*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suXXpg(int n, Particle const &X1, Particle const &X2, double s)
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
		return (-16*std::pow(mc,2)*s*(2*std::pow(mc,2) + s)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suXXpp(int n, Particle const &X1, Particle const &X2, double s)
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
		return (128*std::pow(mc,6)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2)) - (32*std::pow(mc,4)*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

 	if (n == 1)
		return (-32*std::pow(mc,4)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return -((-8*std::pow(mc,2)*(4*std::pow(mc,2) - s) + 2*(4*std::pow(mc,2) - s)*s)*std::pow(lPScc[j],2)*std::pow(lScc[i],2));

 	if (n == 1)
		return -2*(4*std::pow(mc,2) - s)*std::pow(lPScc[j],2)*std::pow(lScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (n == 0)
		return -((32*std::pow(mc,4)*std::pow(bXcc[j],2) - 24*std::pow(mc,2)*s*std::pow(bXcc[j],2) + 4*std::pow(s,2)*std::pow(bXcc[j],2))*std::pow(lScc[i],2));

 	if (n == 1)
		return -16*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suphisXp(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (mXj == 0) 
		return 0;

 	if (n == 0)
		return -(((-32*std::pow(mc,4)*(4*std::pow(mc,2) - s)*std::pow(bXcc[j],2))/std::pow(mXj,2) + (8*std::pow(mc,2)*(4*std::pow(mc,2) - s)*s*std::pow(bXcc[j],2))/std::pow(mXj,2))*std::pow(lScc[i],2));

 	if (n == 1)
		return (-8*std::pow(mc,2)*(4*std::pow(mc,2) - s)*std::pow(bXcc[j],2)*std::pow(lScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suphipXg(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (n == 1)
		return 0;

 	if (n == 0)
		return -4*s*(2*std::pow(mc,2) + s)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_suphipXp(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (mXj == 0) 
		return 0;

 	if (n == 0)
		return ((32*std::pow(mc,4)*s*std::pow(bXcc[j],2))/std::pow(mXj,2) - (8*std::pow(mc,2)*std::pow(s,2)*std::pow(bXcc[j],2))/std::pow(mXj,2))*std::pow(lPScc[i],2);

 	if (n == 1)
		return (-8*std::pow(mc,2)*s*std::pow(bXcc[j],2)*std::pow(lPScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuphisphis(int n, Particle const &s1, Particle const &s2, double s)
{
	//return 0;
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return -8*std::pow(mc,2)*s*std::pow(lScc[i],2)*std::pow(lScc[j],2);

 	if (n == 1)
		return (8*std::pow(mc,2) - 2*s)*std::pow(lScc[i],2)*std::pow(lScc[j],2);

 	if (n == 2)
		return -2*std::pow(lScc[i],2)*std::pow(lScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuphipphip(int n, Particle const &ps1, Particle const &ps2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return 0;

 	if (n == 1)
		return (8*std::pow(mc,2) - 2*s)*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

 	if (n == 2)
		return -2*std::pow(lPScc[i],2)*std::pow(lPScc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuXXgg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (n == 1)
		return 0;

	if (n == 2)
		return 0;

 	if (n == 0)
		return 8*(12*std::pow(mc,4) - 4*std::pow(mc,2)*s + std::pow(s,2))*std::pow(bXcc[i],2)*std::pow(bXcc[j],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuXXgppg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXj == 0 || mXi == 0) 
		return 0;

 	if (n == 0)
		return (-384*std::pow(mc,6)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) + (160*std::pow(mc,4)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (16*std::pow(mc,2)*std::pow(s,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

 	if (n == 1)
		return (-32*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2) + (160*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (32*std::pow(mc,2)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

 	if (n == 2)
		return (-16*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2) - (16*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuXXpg(int n, Particle const &X1, Particle const &X2, double s)
{
	int i = X1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXi = X1.get_mass(); 
	double mXj = X2.get_mass();

 	if (mXi == 0) 
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return (-32*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

 	if (n == 2)
		return (-16*std::pow(mc,2)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/std::pow(mXi,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuXXpp(int n, Particle const &X1, Particle const &X2, double s)
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
		return (128*std::pow(mc,6)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2)) - (32*std::pow(mc,4)*s*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

 	if (n == 2)
		return (-32*std::pow(mc,4)*std::pow(bXcc[i],2)*std::pow(bXcc[j],2))/(std::pow(mXi,2)*std::pow(mXj,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuphisphip(int n, Particle const &s1, Particle const &ps2, double s)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return -((-32*std::pow(mc,4) + 8*std::pow(mc,2)*s)*std::pow(lPScc[j],2)*std::pow(lScc[i],2));

 	if (n == 1)
		return -((16*std::pow(mc,2) - 2*s)*std::pow(lPScc[j],2)*std::pow(lScc[i],2));

 	if (n == 2)
		return 2*std::pow(lPScc[j],2)*std::pow(lScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuphisXg(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (n == 0)
		return -((32*std::pow(mc,4)*std::pow(bXcc[j],2) + 16*std::pow(mc,2)*s*std::pow(bXcc[j],2))*std::pow(lScc[i],2));

 	if (n == 1)
		return 24*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lScc[i],2);

 	if (n == 2)
		return -4*std::pow(bXcc[j],2)*std::pow(lScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuphisXp(int n, Particle const &s1, Particle const &X2, double s)
{
	int i = s1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (mXj == 0) 
		return 0;

 	if (n == 0)
		return -(((-128*std::pow(mc,6)*std::pow(bXcc[j],2))/std::pow(mXj,2) + (32*std::pow(mc,4)*s*std::pow(bXcc[j],2))/std::pow(mXj,2))*std::pow(lScc[i],2));

 	if (n == 1)
		return -(((64*std::pow(mc,4)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (8*std::pow(mc,2)*s*std::pow(bXcc[j],2))/std::pow(mXj,2))*std::pow(lScc[i],2));

 	if (n == 2)
		return (8*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuphipXg(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

 	if (n == 0)
		return 0;

 	if (n == 1)
		return -8*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2);

 	if (n == 2)
		return -4*std::pow(bXcc[j],2)*std::pow(lPScc[i],2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q_tuphipXp(int n, Particle const &ps1, Particle const &X2, double s)
{
	int i = ps1.get_specific_index(); 
	int j = X2.get_specific_index();
	double mXj = X2.get_mass();

	if (mXj == 0) 
		return 0;

 	if (n == 0)
		return 0;

 	if (n == 1)
		return ((32*std::pow(mc,4)*std::pow(bXcc[j],2))/std::pow(mXj,2) - (8*std::pow(mc,2)*s*std::pow(bXcc[j],2))/std::pow(mXj,2))*std::pow(lPScc[i],2);

 	if (n == 2)
		return (-8*std::pow(mc,2)*std::pow(bXcc[j],2)*std::pow(lPScc[i],2))/std::pow(mXj,2);

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



int CrossSectionMajorana_cccc::Q_order(ParticlePair const&pp)
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
		return -1;

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
		return -1;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
		return 2;


	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 2;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 2;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
		return 2;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return -1;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
		return -1;

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
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


	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 1;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 1;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
		return 2;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 1;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
		return 1;

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
		return 1;


	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 2;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 2;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
		return 2;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 2;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
		return 2;

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
		return 2;

	return -1;

}

// Generated with Mathematica - FeynCalc
dcomp CrossSectionMajorana_cccc::Q(ParticlePair const & pp,int n, double s, double t)
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
	// UU terms
	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_uuphisphis(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_uuphipphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()>0)
			return Q_uuXXgg(n,pp.get_first(),pp.get_second(),s)+Q_uuXXgppg(n,pp.get_first(),pp.get_second(),s)+Q_uuXXpp(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()==0)
			return  Q_uuXXgg(n,pp.get_second(),pp.get_first(),s)+Q_uuXXpg(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()==0 && pp.get_second().get_mass()>0)
			return Q_uuXXgg(n,pp.get_first(),pp.get_second(),s)+Q_uuXXpg(n,pp.get_second(),pp.get_first(),s);
		else
			return Q_uuXXgg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_uuphisphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_uuphisXg(n,pp.get_first(),pp.get_second(),s)+Q_uuphisXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_uuphisXg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::UU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_uuphipXg(n,pp.get_first(),pp.get_second(),s)+Q_uuphipXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_uuphisXg(n,pp.get_first(),pp.get_second(),s);
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
	// --------------------------------------

	// --------------------------------------
	// SU terms
	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_suphisphis(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_suphipphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()>0)
			return Q_suXXgg(n,pp.get_first(),pp.get_second(),s)+Q_suXXgppg(n,pp.get_first(),pp.get_second(),s)+Q_suXXpp(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()==0)
			return  Q_suXXgg(n,pp.get_second(),pp.get_first(),s)+Q_suXXpg(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()==0 && pp.get_second().get_mass()>0)
			return Q_suXXgg(n,pp.get_first(),pp.get_second(),s)+Q_suXXpg(n,pp.get_second(),pp.get_first(),s);
		else
			return Q_suXXgg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_suphisphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_suphisXg(n,pp.get_first(),pp.get_second(),s)+Q_suphisXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_suphisXg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::SU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_suphipXg(n,pp.get_first(),pp.get_second(),s)+Q_suphipXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_suphisXg(n,pp.get_first(),pp.get_second(),s);
	}
	// --------------------------------------

	// --------------------------------------
	// TU terms
	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_tuphisphis(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_tuphipphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::vector && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()>0)
			return Q_tuXXgg(n,pp.get_first(),pp.get_second(),s)+Q_tuXXgppg(n,pp.get_first(),pp.get_second(),s)+Q_tuXXpp(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()>0 && pp.get_second().get_mass()==0)
			return  Q_tuXXgg(n,pp.get_second(),pp.get_first(),s)+Q_tuXXpg(n,pp.get_first(),pp.get_second(),s);
		else if (pp.get_first().get_mass()==0 && pp.get_second().get_mass()>0)
			return Q_tuXXgg(n,pp.get_first(),pp.get_second(),s)+Q_tuXXpg(n,pp.get_second(),pp.get_first(),s);
		else
			return Q_tuXXgg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_tuphisphip(n,pp.get_first(),pp.get_second(),s);

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_tuphisXg(n,pp.get_first(),pp.get_second(),s)+Q_tuphisXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_tuphisXg(n,pp.get_first(),pp.get_second(),s);
	}

	if (pp.get_type()==INT_TYPE::TU && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::vector)
	{
		if (pp.get_second().get_mass()>0)
			return Q_tuphipXg(n,pp.get_first(),pp.get_second(),s)+Q_tuphipXp(n,pp.get_first(),pp.get_second(),s);
		else
			return Q_tuphisXg(n,pp.get_first(),pp.get_second(),s);
	}

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

