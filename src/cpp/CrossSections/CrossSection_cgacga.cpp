#include "../../headers/CrossSections/CrossSection_cgacga.h"

CrossSection_cgacga::CrossSection_cgacga(DarkSectorModel const &DSMod, Particle const &chii) : CrossSection(chii, *(StandardModel::getInstance()->get_photon_ptr()), chii,  *(StandardModel::getInstance()->get_photon_ptr()))
{

	int n_S = DSMod.get_n_DS_phis();
	int n_PS = DSMod.get_n_DS_phip();

	lSff.resize(n_S);
	lPSff.resize(n_PS);
	lScc.resize(n_S);
	lPScc.resize(n_PS);

    

	for (int i = 0; i < n_S; i++)
	{
        lSff[i].resize(StandardModel::getInstance()->get_n_couples_SM_ferm());
        for(int j = 0; j < lSff[i].size(); j++)
		    lSff[i][j] = DSMod.get_lambda_S_SM(i, j);

		lScc[i] = DSMod.get_lambda_S_DM(i, chii.get_specific_index(), chii.get_specific_index());
        // Specific index is 
	}
	for (int i = 0; i < n_PS; i++)
	{
        lPSff[i].resize(StandardModel::getInstance()->get_n_couples_SM_ferm());
        for(int j = 0; j < lPSff[i].size(); j++)
		    lPSff[i][j] = DSMod.get_lambda_PS_SM(i, j);
		
        lPScc[i] = DSMod.get_lambda_PS_DM(i, chii.get_specific_index(), chii.get_specific_index());

		//std::cout << i << " " << lSff[i] << " " << lPScc[i] << std::endl;
	}


	// First cross section
	std::vector<Particle> as, at, au;

	int const ns = 0;
	int const nt = n_S + n_PS;
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

	mi = _m1;

	Initialise(as, at, au, "cgacga");
}

dcomp CrossSection_cgacga::Fphis(int i, double q2)
// Here i represents the index on the propagator
// if q2 <= 0 we assume energy echanged to be arbitrarily small and we sum over every contrubution
{
    dcomp res = 0, floop(0,0);
    double mf = 0, qf2 = 0;

    for(int j = 0; j < StandardModel::getInstance()->get_n_couples_SM_ferm(); j++)
    {
        mf = StandardModel::getInstance()->get_couples_SM_ferm(j).get_first().get_mass();
		qf2 = pow(StandardModel::getInstance()->get_couples_SM_ferm(j).get_first().get_electric_charge(), 2);
		
		if(q2 > 0)
		{
			if(mf > 0)
				floop = fOneLoopScalarVectors(4 * mf * mf / q2) / mf;
			else
				floop = 0; // limit of fOneLoopScalaVector/mf in mf->0
		}
		else
		{
			if(mf > 0)
				floop = 2./3. / mf; // limit of fOneLoopScalarVectors/mf in q2->0
			else
				floop = 0;// limit of fOneLoopScalaVector/mf in mf->0
		}


		//std::cout << j << " " << 4 * mf * mf / q2 << " " << floop << std::endl;

		if (StandardModel::getInstance()->get_couples_SM_ferm(j).get_first().get_QCDType() == QCDTYPE::QUARK)
            res += 3 * qf2 * lSff[i][j] * floop;
        else
            res += qf2 * lSff[i][j] * floop;
    }
	
	return one_i * res;
}

dcomp CrossSection_cgacga::Fphip(int i, double q2)
{
    dcomp res = 0, floop(0, 0);
    double mf = 0, qf2 = 0;

    for(int j = 0; j < StandardModel::getInstance()->get_n_couples_SM_ferm(); j++)
    {
        mf = StandardModel::getInstance()->get_couples_SM_ferm(j).get_first().get_mass();
		qf2 = pow(StandardModel::getInstance()->get_couples_SM_ferm(j).get_first().get_electric_charge(), 2);
		
		if(q2 > 0)
		{
			if(mf > 0)
				floop = fOneLoopPseudoScalarVectors(4 * mf * mf / q2) / mf;
			else
				floop = 0; // limit of fOneLoopPseudoScalaVector/mf in mf->0
		}
		else
		{
			if(mf > 0)
				floop = 2./3. / mf; // limit of fOneLoopPseudoScalaVector/mf in q2->0
			else
				floop = 0;// limit of fOneLoopPseudoScalaVector/mf in mf->0
		}


		if (StandardModel::getInstance()->get_couples_SM_ferm(j).get_first().get_QCDType() == QCDTYPE::QUARK)
            res += 3 * qf2 * lPSff[i][j] * floop;
        else
            res += qf2 * lPSff[i][j] * floop;
    }
	
	return -res;

}


// Generated with Mathematica - FeynCalc
dcomp CrossSection_cgacga::Q_ttphisphis(int n, Particle const &s1, Particle const &s2, double s, double t)
{
	int i = s1.get_specific_index(); 
	int j = s2.get_specific_index();

  	if (n == 0)
		return 0;

 	if (n == 1)
		return 0;

 	if (n == 2)
		return (4*std::pow(mi,2)*std::pow(aem(t),2)*conj(Fphis(i,t))*Fphis(j,t)*lScc[i]*lScc[j])/std::pow(PI,2);

 	if (n == 3)
		return -((std::pow(aem(t),2)*conj(Fphis(i,t))*Fphis(j,t)*lScc[i]*lScc[j])/std::pow(PI,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cgacga::Q_ttphipphip(int n, Particle const &ps1, Particle const &ps2, double s, double t)
{
	int i = ps1.get_specific_index(); 
	int j = ps2.get_specific_index();

  	if (n == 0)
		return 0;

 	if (n == 1)
		return 0;

 	if (n == 2)
		return 0;

 	if (n == 3)
		return -((std::pow(aem(t),2)*conj(Fphip(i,t))*Fphip(j,t)*lPScc[i]*lPScc[j])/std::pow(PI,2));

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cgacga::Q_ttphisphip(int n, Particle const &s1, Particle const &ps2, double s, double t)
{
	int i = s1.get_specific_index(); 
	int j = ps2.get_specific_index();

 	return 0;

	printErrorMessage(__PRETTY_FUNCTION__);
	exit(0); 
}



// Generated with Mathematica - FeynCalc
int CrossSection_cgacga::Q_order(ParticlePair const&pp)
{
	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return 3;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return 3;

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return -1;

	return -1;
}

// Generated with Mathematica - FeynCalc
dcomp CrossSection_cgacga::Q(ParticlePair const&pp, int n, double s, double t)
{
	if (pp.get_type()== INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::scalar)
		return Q_ttphisphis(n, pp.get_first(), pp.get_second(), s, t);

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::pseudoscalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_ttphipphip(n, pp.get_first(), pp.get_second(), s, t);

	if (pp.get_type()==INT_TYPE::TT && pp.get_first().get_proptype()==Proptype::scalar && pp.get_second().get_proptype()==Proptype::pseudoscalar)
		return Q_ttphisphip(n, pp.get_first(), pp.get_second(), s, t);

	return -1;
}
