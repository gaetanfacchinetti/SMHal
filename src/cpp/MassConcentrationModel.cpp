#include "../headers/MassConcentrationModel.h"

double MassConcentrationModel::RedshiftOfCollapse_vs_Mass(double M)
{
    std::vector<double> zM;
    zM.resize(2);
    zM[1] = M; 

    double S = _power_spectrum.SigmaM2(M, 0);
    double dc = _mass_function->get_delta_c();
    double D0 = _cosmo.growth_factor_D1_Carroll(0);
    return _cosmo.Inverse_growth_factor_Caroll(dc*D0/sqrt(S));
    //return sqrt(_power_spectrum.SigmaM2(M, 0)/_mass_function->get_delta_c()/_mass_function->get_delta_c()) - 1;
    //return DichotomieLn_Static(0, 2, zM, this, CallBack_fToDichotomieFor_RedshiftOfCollapse_vs_Mass, 1e-1, 1e+5, 1e-3, 0);
}


double MassConcentrationModel::Mass_vs_RedshiftOfCollapse(double zc)
{
    std::vector<double> zM;
    zM.resize(2);
    zM[0] = zc; 

    return DichotomieLn_Static(1, 2, zM, this, CallBack_fToDichotomieFor_RedshiftOfCollapse_vs_Mass, 1e-12, 1e+16, 1e-3, 0);
}

double MassConcentrationModel::fToDichotomieFor_RedshiftOfCollapse_vs_Mass(std::vector<double> zM)
{

    //double SF = _power_spectrum.SigmaM2(0.05*zM[1], zM[0]);
    double S  = _power_spectrum.SigmaM2(zM[1], zM[0]);
    double dc = _mass_function->get_delta_c();
    
    return  S-dc*dc; 
}

double MassConcentrationModel::ProbaRedshiftOfCollapse(double zc, double M, bool norm)
{
    double D1_zc = _cosmo.growth_factor_exact(zc);
    double dc = _mass_function->get_delta_c();
    double sigma_M = _power_spectrum.SigmaM(M, zc);
    double sigma_M_0 = _power_spectrum.SigmaM(M,0);
    double derLnD1 = fabs(_cosmo.derLn_growth_factor_exact(zc)); 

    double N = 1;

    if(norm == true)
        N = (1-erf(dc/(sqrt(2)*sigma_M_0)));

    return sqrt(2./PI)/N*(dc/(sigma_M))*derLnD1*exp(-0.5*dc*dc/(sigma_M*sigma_M)) ;
}

double MassConcentrationModel::CumulativeRedshiftOfCollapse(double zc, double M)
{
    std::vector<double> zcM{0, M, 0};
    return GaussLegendre_IntegralLn_Static(0, 3, gl100, zc, 1e+7, zcM, this, CallBack_fToIntFor_MomentProbazc);
}



double MassConcentrationModel::MeanRedshiftOfCollapse(double M)
{
    std::vector<double> zcM{0, M, 1};
    return GaussLegendre_IntegralLn_Static(0, 3, gl1000, 1e-5, 1000, zcM, this, CallBack_fToIntFor_MomentProbazc);
}

double MassConcentrationModel::fToIntOnM2_RedshiftOfFormation_LaceyCole(double St, double omFt, double M1, double S1, double Sh)
// M1 is the mass of the final halo, M2, that is integrated on, is the mass of the progenitor
// z1 is the collapse time of halo 1 of mass M1 and zF is the formation time of halo 1 of mass M1
{
    double S2 = (Sh - S1)*St + S1;
    double M2 = _power_spectrum.Mass_vs_Sigma2(S2, 0);

    return (M1/M2)*omFt/pow(St, 3./2.) * (1-pow(omFt, 2)/St)* exp(-pow(omFt, 2)/ 2./ St);
} 



double MassConcentrationModel::fToIntOnM2_IntegratedProbaRedshiftOfFormation_LaceyCole(double St, double omFt, double M1, double S1, double Sh)
// M1 is the mass of the final halo, M2, that is integrated on, is the mass of the progenitor
// z1 is the collapse time of halo 1 of mas M1 and zF is the formation time of halo 1 of mass M1
{
    double S2 = (Sh - S1)*St + S1;
    double M2 = _power_spectrum.Mass_vs_Sigma2(S2, 0);

    return (M1/M2)*omFt/pow(St, 3./2.) * exp(-pow(omFt, 2)/ 2./ St);
} 


void MassConcentrationModel::plot_fToIntOnM2_RedshiftOfFormation_LaceyCole(double omFt, double M1)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/fToIntOnM2_RedshiftOfFormation_LaceyCole_vs_M2.out");
    outfile << "# M2 [Msol] | fToIntOnM2_RedshiftOfFormation_LaceyCole " << std::endl;


    double S1 = _power_spectrum.SigmaM2(M1, 0);
    double Sh = _power_spectrum.SigmaM2(M1/2., 0);

    int N = 1000;
    double dSt = log10(1./0.001)/N, St;


    for(int i = 0; i < N+1; i++)
    {
        St = 0.001*pow(10, i*dSt);
        outfile << St << "\t" << fToIntOnM2_RedshiftOfFormation_LaceyCole(St, omFt, M1, S1, Sh) << std::endl;
    }

    outfile.close();
}


double MassConcentrationModel::ProbaRedshiftOfFormation_LaceyCole_vs_omFt(double omFt, double M1)
{

    double S1 = _power_spectrum.SigmaM2(M1, 0);
    double Sh = _power_spectrum.SigmaM2(M1/2., 0);

    std::vector<double> xx = {0, omFt, M1, S1, Sh};
    double integral = - GaussLegendre_Integral_Static(0, 5, gl100, 1e-3, 1, xx, this, CallBack_fToIntOnM2_RedshiftOfFormation_LaceyCole);
    
    return integral/sqrt(2*PI);
}

double MassConcentrationModel::IntegratedProbaRedshiftOfFormation_LaceyCole_vs_omFt(double omFt, double M1)
{

    double S1 = _power_spectrum.SigmaM2(M1, 0);
    double Sh = _power_spectrum.SigmaM2(M1/2., 0);
    
    std::vector<double> M2M1z2z1 = {0, omFt, M1, S1, Sh};
    double integral = GaussLegendre_IntegralLn_Static(0, 5, gl100, 1e-3, 1, M2M1z2z1, this, CallBack_fToIntOnM2_IntegratedProbaRedshiftOfFormation_LaceyCole);

    return integral/sqrt(2*PI);
}

double MassConcentrationModel::ProbaRedshiftOfFormation_LaceyCole_vs_zF(double zF, double M1, double z1)
{

    double D1_z0 = _cosmo.growth_factor_D1_Carroll(0);
    double D1_z1 = _cosmo.growth_factor_D1_Carroll(z1);
    double D1_zF = _cosmo.growth_factor_D1_Carroll(zF);
    double om1 = _mass_function->get_delta_c()/D1_z1*D1_z0;
    double omF = _mass_function->get_delta_c()/D1_zF*D1_z0;

    double S1 = _power_spectrum.SigmaM2(M1, 0);
    double Sh = _power_spectrum.SigmaM2(M1/2., 0);

    double omFt = (omF - om1)/sqrt(Sh-S1);
    
    double der_omFt_zF = omF/sqrt(Sh-S1)*fabs(_cosmo.derLn_growth_factor_exact(zF));

    double res = ProbaRedshiftOfFormation_LaceyCole_vs_omFt(omFt, M1)*der_omFt_zF;

    if (res < 0)
        return 0;
    
    return res;
}



void MassConcentrationModel::plot_ProbaRedshiftOfFormation_LaceyCole_vs_omFt(double M1)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/ProbaRedshiftOfFormation_LaceyCole_vs_omegaFtilde.out");
    outfile << "# zF | p(z_f) " << std::endl;

    int N = 50;
    double d_omt = 4./N, omt;

    for(int i = 0; i < N+1; i++)
    {
        omt = i*d_omt;
        outfile << omt << "\t" << ProbaRedshiftOfFormation_LaceyCole_vs_omFt(omt, M1) << std::endl;
    }

    outfile.close();
}

void MassConcentrationModel::plot_ProbaRedshiftOfFormation_LaceyCole_vs_zF(double M1, double z1)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/ProbaRedshiftOfFormation_LaceyCole_vs_zF.out");
    outfile << "# zF | p(z_f) " << std::endl;

    int N = 100;
    double dzF = log10(std::max(100., 100.*z1)/z1)/N, zF;

    for(int i = 0; i < N+1; i++)
    {
        zF = z1*pow(10, i*dzF);
        outfile << zF << "\t" << ProbaRedshiftOfFormation_LaceyCole_vs_zF(zF, M1, z1) << std::endl;
    }

    outfile.close();
}




void MassConcentrationModel::plot_IntegratedProbaRedshiftOfFormation_LaceyCole_vs_omFt(double M1)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/IntegratedProbaRedshiftOfFormation_LaceyCole_vs_omegaFtilde.out");
    outfile << "# z2 | p(zf < z2) " << std::endl;

    int N = 100;
    double d_omt = 4./N, omt;

    for(int i = 0; i < N+1; i++)
    {
        omt = i*d_omt;
        outfile << omt << "\t" << IntegratedProbaRedshiftOfFormation_LaceyCole_vs_omFt(omt, M1) << std::endl;
    }

    outfile.close();
}







//
// Maccio and Dutton model
//

double MassConcentrationModel::ConcentrationMaccio(double zc, double z)
{
    double K = 3.8;
    return K*pow(_cosmo.Hubble_parameter(zc)/_cosmo.Hubble_parameter(z), 2./3.);
}

double MassConcentrationModel::ConcentrationMaccioMass(double z, double M)
{
    double zc = RedshiftOfCollapse_vs_Mass(M/100.);
    return ConcentrationMaccio(zc, z);
}

double MassConcentrationModel::MeanConcentrationMaccio(double z, double M)
{
    std::vector<double> zcMz{0, M, z};
    return GaussLegendre_IntegralLn_Static(0, 3, gl100, 1e-5, 1000, zcMz, this, CallBack_fToIntFor_MeanConcentrationMaccio);
}

double MassConcentrationModel::SigmaConcentrationMaccio(double z, double M)
// Attention gives the result as sigma_log10c (in dex)
{
    std::vector<double> zcMz1{0, M, z, 1};
    std::vector<double> zcMz2{0, M, z, 2};
    double meanmean = pow(GaussLegendre_IntegralLn_Static(0, 4, gl1000, 1e-5, 1000, zcMz1, this, CallBack_fToIntFor_Log10ConcentrationNMaccio), 2);
    double mean2 = GaussLegendre_IntegralLn_Static(0, 4, gl1000, 1e-5, 1000, zcMz2, this, CallBack_fToIntFor_Log10ConcentrationNMaccio);
    
    return sqrt(mean2-meanmean); 
}


//
//
double MassConcentrationModel::fToIntOnM_for_dndzc(double zc, double M, double z)
{
    std::cout << zc << " " << M << " " << ProbaRedshiftOfCollapse(zc, M) << " " << _mass_function->NumberDensity(M, zc) << std::endl;
    double res = ProbaRedshiftOfCollapse(zc, M)*_mass_function->NumberDensity(M, z);


    if(res != res)
        return 0;

    return res;
}

double MassConcentrationModel::dndzc(double zc, double z)
{
    std::vector<double> xx = {zc, 0, z};
    return GaussLegendre_IntegralLn_Static(1, 3, gl200, 1e-14, 1e+12, xx, this, CallBack_fToIntOnM_for_dndzc);
}

void MassConcentrationModel::plot_dndzc(double z)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/dndzc.out");
    outfile.precision(8);

    int Npts = 200;
    double zc_min = 10, zc_max = 500, d_zc = log10(zc_max/zc_min)/Npts, zc;


    for(int i = 0; i< Npts+1; i++)
    {
        zc = zc_min*pow(10, i*d_zc);
        outfile << zc << "\t" << dndzc(zc, z) <<std::endl;
        exit(0);
    }

    outfile.close();
}   


//
//
//
//

double MassConcentrationModel::ConcentrationFacchinetti(double zc, double z)
{
    std::vector<double> xx = {0, z, zc};
    return DichotomieLn_Static(0, 3, xx, this, CallBack_fForBissection_ConcentrationFacchinetti, 1e-5, 1e+6, 1e-3, 0);
}

double MassConcentrationModel::fForBissection_ConcentrationFacchinetti(double c, double z, double zc)
{
    double res = 0;
    // For the moment we just look at the NFW case
    double mu_c = 0;
    if(c > 1e-2)
        mu_c = log(1+c) - c/(1+c);
    else
        mu_c = pow(c, 2)/2. - 2*pow(c, 3)/3 + 3*pow(c, 4)/4.;

    res =  pow(c, 3.) - pow(4.2, 3) * _cosmo.critical_density(zc)/_cosmo.critical_density(z);

    //std::cout << c << " " << res << std::endl;
    return res;
}

double MassConcentrationModel::RedshiftOfCollapse_Facchinetti(double c, double m, double z)
{
    std::vector<double> xx = {0, c, m, z};
    //return RedshiftOfCollapse_vs_Mass(m/100);
    return DichotomieLn_Static(0, 4, xx, this, CallBack_fForBissection_RedshiftOfCollapse_Facchinetti, 1e-2, 1e+5, 1e-3, 0);
}

double MassConcentrationModel::fForBissection_RedshiftOfCollapse_Facchinetti(double zc, double c, double m, double z)
{
    double c_zc = ConcentrationFacchinetti(zc, zc);
    double m_zc = Mass_vs_RedshiftOfCollapse(zc);


    double mu_c_zc, mu_c;
    
    if(c_zc > 1e-2)
        mu_c_zc = log(1+c_zc) - c_zc/(1+c_zc);
    else
        mu_c_zc = pow(c_zc, 2)/2. - 2*pow(c_zc, 3)/3. + 3*pow(c_zc, 4)/4.;
    
    if(c > 1e-2)
        mu_c = log(1+c) - c/(1+c);
    else
        mu_c = pow(c, 2)/2. - 2*pow(c, 3)/3. + 3*pow(c, 4)/4.;
    
    //std::cout << zc << "\t" << m_zc - m*mu_c/mu_c_zc << "\t" << mu_c << "\t" << mu_c_zc << "\t" << c_zc << std::endl;

    return m_zc - m*mu_c_zc/mu_c;
}

double MassConcentrationModel::ConcentrationFacchinetti_vs_m(double m, double z)
{
    std::vector<double> xx = {0, m, z};
    return DichotomieLn_Static(0, 3, xx, this, CallBack_fForBissection_ConcentrationFacchinetti_vs_m, 1e-3, 1e+3, 1e-3, 0);
}

double MassConcentrationModel::fForBissection_ConcentrationFacchinetti_vs_m(double c, double m, double z)
{
    return ConcentrationFacchinetti(RedshiftOfCollapse_Facchinetti(c, m, z), z) - c;
}


//
//
// Okoli and Afshordi model
//

double MassConcentrationModel::ConcentrationOkoli(double z, double m200)
{

    if(_power_spectrum.get_window() != PSWindowType::real_space_top_hat)
    {
        std::cout << "FATAL ERROR: cannot use ConcentrationOkoli with a window that is not real_space_top_hat" << std::endl;
        exit(0);
    }

    double zc = RedshiftOfCollapse_vs_Mass(m200);
    double nu = _mass_function->get_delta_c()/_power_spectrum.SigmaM(m200, zc);
    double Ht = _cosmo.CosmicTime_times_Hubble_parameter(zc)*_cosmo.Hubble_parameter(z)/_cosmo.Hubble_parameter(zc);
    double y = (0.42 + 0.20*pow(nu, -1.23))/pow(Ht, 2./3.);
    std::cout << m200 << " " << z << " " << zc << std::endl;
    return pow(10, 0.78*log(y) + 1.09);
}

double MassConcentrationModel::A2_Okoli(double z, double R)
// Use the usual definition of the variance dor the real space top hat
{
    return pow(4*PI/3.,2)*_power_spectrum.SigmaR2(R, z);
}

double MassConcentrationModel::B2_Okoli(double z, double R)
{
    double lnk_min = -6, lnk_max = 3.0 - log(R);

    if (lnk_max <= lnk_min)
        return 0;

    std::vector<double> lnkRz;
    lnkRz.resize(3);
    lnkRz[1] = R;
    lnkRz[2] = z;

    return GaussLegendre_Integral_Static(0, 3, gl1000, lnk_min, lnk_max, lnkRz, this, CallBack_fToIntFor_B2_Okoli);
}

double MassConcentrationModel::AB_Okoli(double z, double R)
// Need to use the new B window (not defined in the Press-Schechter class) for this computation
{
    double lnk_min = -6, lnk_max = 3.0 - log(R);

    if (lnk_max <= lnk_min)
        return 0;

    std::vector<double> lnkRz;
    lnkRz.resize(3);
    lnkRz[1] = R;
    lnkRz[2] = z;

    return GaussLegendre_Integral_Static(0, 3, gl1000, lnk_min, lnk_max, lnkRz, this, CallBack_fToIntFor_AB_Okoli);

}

double MassConcentrationModel::fToIntFor_AB_Okoli(std::vector<double> lnkRz)
{
    double lnk = lnkRz[0];
    double R = lnkRz[1];
    double z = lnkRz[2];

    double k = exp(lnk);
    double power_spectrum = pow(k, 3) * _power_spectrum.matter_power_spectrum(k, z) / (2 * PI * PI);

    return power_spectrum * (4*PI/3.) * _power_spectrum.Window(k * R) * Window_B_Okoli(k*R);
}

double MassConcentrationModel::fToIntFor_B2_Okoli(std::vector<double> lnkRz)
{
    double lnk = lnkRz[0];
    double R = lnkRz[1];
    double z = lnkRz[2];

    double k = exp(lnk);
    double power_spectrum = pow(k, 3) * _power_spectrum.matter_power_spectrum(k, z) / (2 * PI * PI);

    return power_spectrum *pow(Window_B_Okoli(k*R), 2);
}

double MassConcentrationModel::Window_B_Okoli(double kR)
{
    return (4*PI*pow(kR, -5))*(6*sin(kR)- 6*kR*cos(kR) - 2*kR*kR*sin(kR));
}



//
//
// Sanchez-Conde and Prada model 
//


double MassConcentrationModel::ConcentrationSanchezCondePrada(double m200)
{
    double cbar = 0;

    double cn[6] = {37.5153, -1.5093, 1.636e-2, 3.66e-4, -2.89237e-5, 5.32e-7};

    double m200_bis = m200;

    if (m200 <= 7.24e-10)
        m200_bis = 7.24e-10;

    for (int i = 0; i < 6; i++)
        cbar += cn[i] * pow(log(m200_bis * _cosmo.get_h()), i);

    return cbar;
}



//
//
//
//

// Plot functions

void MassConcentrationModel::plot_ProbaRedshiftOfCollapse_vs_zc(double M)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/ProbaRedshiftOfCollapse_vs_zc.out");
    outfile.precision(8);

    int Npts = 1000;
    double zmin = 1e-2, zmax = 1e+4, d_logz = log10(zmax/zmin)/Npts, z;
    //double zmin = 0, zmax = 100, dz = (zmax - zmin)/Npts, z;

    for(int i = 0; i< Npts+1; i++)
    {
        z = zmin*pow(10, i*d_logz);
        //z = zmin +i*dz;
        outfile << z << "\t" << ProbaRedshiftOfCollapse(z, M)  << std::endl;
    }

    outfile.close();
}   


void MassConcentrationModel::plot_ProbaRedshiftOfCollapse_vs_M(double zc)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/ProbaRedshiftOfCollapse_vs_M.out");
    outfile.precision(8);

    int Npts = 1000;
    double Mmin = 1e-12, Mmax = 1e+14, d_logM = log10(Mmax/Mmin)/Npts, M;
    //double zmin = 0, zmax = 100, dz = (zmax - zmin)/Npts, z;

    for(int i = 0; i< Npts+1; i++)
    {
        M = Mmin*pow(10, i*d_logM);
        //z = zmin +i*dz;
        outfile << M << "\t" << ProbaRedshiftOfCollapse(zc, M)  << std::endl;
    }

    outfile.close();
}   


void MassConcentrationModel::plot_CumulativeRedshiftOfCollapse_vs_zc(double M)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/CumulativeRedshiftOfCollapse_vs_zc.out");
    outfile.precision(8);

    int Npts = 100;
    double zmin = 1e-2, zmax = 1e+4, d_logz = log10(zmax/zmin)/Npts, z;
    //double zmin = 0, zmax = 100, dz = (zmax - zmin)/Npts, z;

    for(int i = 0; i< Npts+1; i++)
    {
        z = zmin*pow(10, i*d_logz);
        //z = zmin +i*dz;
        outfile << z << "\t" << CumulativeRedshiftOfCollapse(z, M) << std::endl;
    }

    outfile.close();
}   


void MassConcentrationModel::plot_B_over_A_Okoli(double z)
{
    std::ofstream outfile;
    outfile.open("../output/A_over_B_Okoli.out");
    outfile.precision(8);

    int Npts = 1000;
    double Rhmin = 1e-1, Rhmax = 50, dRh = (Rhmax-Rhmin)/Npts, Rh;
    //double zmin = 0, zmax = 100, dz = (zmax - zmin)/Npts, z;
    double A2, B2, AB;
    double mean, sigma;

    for(int i = 0; i< Npts+1; i++)
    {
        Rh = Rhmin + i*dRh;
        //z = zmin +i*dz;
        A2 = A2_Okoli(z, Rh/_cosmo.get_h());
        B2 = B2_Okoli(z, Rh/_cosmo.get_h());
        AB = AB_Okoli(z, Rh/_cosmo.get_h());
        mean = AB/A2;
        sigma = sqrt((AB*AB-A2*B2)/(A2*A2));
        outfile << Rh << "\t" << mean << "\t" << mean - sigma << "\t" << mean+sigma << std::endl;
    }

    outfile.close();
}   


void MassConcentrationModel::plot_concentration(double z)
{
    std::ofstream outfile;
    outfile.open("../output/UCMHs/Concentrations.out");
    outfile.precision(8);

    int Npts = 20;
    double Mhmin = 1e-10, Mhmax = 1e+15, d_logMh = log10(Mhmax/Mhmin)/Npts, Mh;
    double c_Facchinetti = 0, zc_Facchinetti;

    for(int i = 0; i< Npts+1; i++)
    {
        Mh = Mhmin*pow(10, i*d_logMh);
        //z = zmin +i*dz;
        c_Facchinetti = ConcentrationFacchinetti_vs_m(Mh, z);
        zc_Facchinetti = RedshiftOfCollapse_Facchinetti(c_Facchinetti, Mh, z);
        outfile << Mh << "\t" << ConcentrationMaccioMass(z, Mh) << "\t" << c_Facchinetti;
         
        if (z == 0)
            outfile << "\t" << ConcentrationSanchezCondePrada(Mh);
        else
            outfile << "\t" << 0 ;

        outfile << "\t" << zc_Facchinetti << "\t" << RedshiftOfCollapse_vs_Mass(Mh/100) << "\t" << Mass_vs_RedshiftOfCollapse(zc_Facchinetti)/Mh << std::endl;
    }

    outfile.close();
}   


void MassConcentrationModel::plot_fToIntFor_AB_Okoli(double z, double R)
{
    std::ofstream outfile;
    outfile.open("../output/fToIntFor_AB_Okoli.out");
    outfile.precision(8);

    int Npts = 200;
    double kmin = exp(-6), kmax = exp(3)/R, d_logk = log10(kmax/kmin)/Npts, k;

    std::vector<double> lnkRz{0, R, z};

    for(int i = 0; i< Npts+1; i++)
    {
        k = kmin*pow(10, i*d_logk);
        lnkRz[0] = log(k);
        //z = zmin +i*dz;
        outfile << k << "\t" << fToIntFor_AB_Okoli(lnkRz) << std::endl;
         
    }

    outfile.close();
}   


/*
void MassConcentrationModel::plot_ProbaRedshiftOfCollapse_vs_k(double zc)
{
    std::ofstream outfile;
    outfile.open("../output/ProbaRedshiftOfCollapse_vs_k.out");
    outfile.precision(8);

    int Npts = 1000;
    double kmin = 1e+10, kmax = 1e+30, d_logk = log10(kmax/kmin)/Npts, k;

    for(int i = 0; i< Npts+1; i++)
    {
        k = kmin*pow(10, i*d_logk);
        outfile << k << "\t" << ProbaRedshiftOfCollapse(zc, k) << std::endl;
    }

    outfile.close();
}   */