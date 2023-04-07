#include "../headers/MinimalMass.h"



MinimalMass::MinimalMass(DegreeFreedom *_degFree, Cosmology Cosmo)
{

    Tstar = 1e-6; //(GeV) i.e. Tstar = 1 keV

    spline_gSmooth = _degFree->get_gSmoothInterp();
    spline_hSmooth = _degFree->get_hSmoothInterp();
    spline_derLogHSmooth = _degFree->get_derLogHSmoothInterp();
    spline_gStarHalf = _degFree->get_gStarHalfInterp();

    zeq = Cosmo.Zeq_rad_mat();
    a0OutOfaeq = zeq + 1;
    aeq = 1./(1+zeq);
    Heq = Cosmo.Hubble_parameter(zeq)/kpc_to_m/s_to_GeVm1; // GeV
    Teq_GeV = T0_CMB_GeV * a0OutOfaeq;
    rhoEq = Cosmo.cosmological_density(Cosmo.Zeq_rad_mat(), Species::MATTER)*Msun_to_kg*kg_to_GeV*pow(kpc_to_cm*cm_to_GeVm1, -3); // GeV^4
    rhoM0 = Cosmo.cosmological_density(0, Species::MATTER)*Msun_to_kg*kg_to_GeV*pow(kpc_to_cm*cm_to_GeVm1, -3); // GeV^4

    //double T = 1;
    //double One_aH = 1/((T_CMB_GeV/T)*sqrt(4.*PI*PI*PI/45.*spline_gSmooth(log10(T)))*T*T/M_PLANCK*GeV_to_cmm1*kpc_to_cm*1e+3);
    //std::cout << One_aH << std::endl; 

    //std::cout << "Intro : " << aeq << " " << zeq << " " << Teq_GeV << std::endl;
    std::clog << std::endl;
    std::clog << "===============================================" << std::endl;
    std::clog << "MINIMAL MASS  INITIALIZATION" << std::endl;
    std::clog << std::endl;

    std::clog.precision(4);
    std::clog << std::scientific;



    std::clog << " -- Cosmological quantities ------- " << std::endl;
    std::clog << "| Teq     : " << Teq_GeV << " (GeV)       |" << std::endl;
    std::clog << "| rho_eq  : " << rhoEq << " (GeV^4)     |" << std::endl;
    std::clog << "| a0/aeq  : " << a0OutOfaeq << "             |" << std::endl;
    //std::clog << "| akd/aeq : " << akdOutOfaeq << "             |" << std::endl;
    //std::clog << "| rhoC    : " << rhoC_GeV4 << " (GeV^4)     |" << std::endl;
    std::clog << "| rhoM0  : " << rhoM0 << " (GeV^4)     |" << std::endl;
    std::clog << " ---------------------------------- " << std::endl;
    std::clog << std::endl;
    std::clog << sqrt(8 * PI * rhoEq / 3) / M_PLANCK << std::endl;

    std::cout << "--> MinimalMass initialised -> see log file" << std::endl;

    // std::cout << "| Input - Kinetic decoupling temperature : " << std::endl;
}

/*
void MinimalMass::set_Tkd(double const &TkdNew)
{
    Tkd = TkdNew;
    akdOutOfaeq = (Teq_GeV / Tkd) * pow(spline_hSmooth(log10(Teq_GeV)), 1. / 3.) / pow(spline_hSmooth(log10(Tkd)), 1. / 3.);
    //std::cout << "|| INFO.   : (Minimal Mass Comp.) Tkd = " << Tkd << " (GeV)" << std::endl;
}
*/

double MinimalMass::FreeStreamingLength(double Tkd, double mDM)
{

    double akdOutOfaeq = (Teq_GeV / Tkd) * pow(spline_hSmooth(log10(Teq_GeV)), 1. / 3.) / pow(spline_hSmooth(log10(Tkd)), 1. / 3.);
    
    std::vector<double> variables;
    variables.resize(1);

    double vkd = sqrt(6 * Tkd / (5 * mDM));

    double C = (1 / a0OutOfaeq) * (sqrt(PI / 3.) / M_PLANCK) * sqrt(rhoEq);

    double Ctau0 = sqrt(1 + Teq_GeV / T0_CMB_GeV) - 1;
    double Ctaukd = sqrt(1 + Teq_GeV / Tkd) - 1;
    double factor = log((Ctau0 / Ctaukd) * ((2 + Ctaukd) / (2 + Ctau0)));

    int nPointsInt = 5000;
    double err;

    double Kfs = Simpson_Integral1_Static(0, 1, 0, nPointsInt, log(Tstar), log(Tkd), variables, this, CallBackFToIntToComputeKfs, err);
    //double Kfs = Integram

    double lambdafs = (vkd / (2 * C)) * akdOutOfaeq * (factor + Kfs);
    return lambdafs;
}


double MinimalMass::ComovingFreeStreamingLength(double Tkd, double mDM, double z) // computation with tespect to the scale factor
// The z dependance simply indicates when do we stop evaluating the free-streaming length
// The true non comoving length is given by this function multplied by 1/(1+z)
{
    /*
    std::vector<double> variables;
    variables.resize(1);

    double akdOutOfaeq = (Teq_GeV / Tkd) * pow(spline_hSmooth(log10(Teq_GeV)), 1. / 3.) / pow(spline_hSmooth(log10(Tkd)), 1. / 3.);

    double vkd = VelocityAtKineticDecoupling(Tkd, mDM);

    double C = (1 / (a * a0OutOfaeq)) * (sqrt(PI / 3.) / M_PLANCK) * sqrt(rhoEq);

    double Ctau0 = sqrt(1 + Teq_GeV / (T_CMB_GeV / a)) - 1;
    double Ctaukd = sqrt(1 + Teq_GeV / Tkd) - 1;
    double factor = log((Ctau0 / Ctaukd) * ((2 + Ctaukd) / (2 + Ctau0)));

    int nPointsInt = 5000;
    double err;

    double Kfs = Simpson_Integral1_NonStatic(0, 1, 0, nPointsInt, log(Tstar), log(Tkd), variables, this, CallBackFToIntToComputeKfs, err);

    double lambdafs = (vkd / (2 * C)) * akdOutOfaeq * (factor + Kfs);
    return lambdafs;*/
    double vkd = VelocityAtKineticDecoupling(Tkd, mDM);
    double kfs = ComovingFreeStreamingScale(Tkd, mDM, z);

    return vkd/kfs*sqrt(2*mDM/Tkd); 
}


double MinimalMass::VelocityAtKineticDecoupling(double Tkd, double mDM)
{
    return sqrt(6 * Tkd /(5 *mDM));
}


double MinimalMass::ComovingCollisionalDampingScale(double Tkd, double mDM) 
// Approximation from http://stacks.iop.org/1475-7516/2005/i=08/a=003?key=crossref.b746e23fee2ea665a08baead49ee77d9
{
    double Hkd = sqrt(4.*PI*PI*PI/45.*spline_gSmooth(log10(Tkd)))*Tkd*Tkd/M_PLANCK; // GeV
    double akdOutOfaeq = (Teq_GeV / Tkd) * pow(spline_hSmooth(log10(Tstar)), 1. / 3.) / pow(spline_hSmooth(log10(Tkd)), 1. / 3.);
    
    return 1.8*sqrt(mDM/Tkd)*akdOutOfaeq*aeq*Hkd;
}

double MinimalMass::ComovingFreeStreamingScale(double Tkd, double mDM, double z) // computation with tespect to the scale factor
{

    std::vector<double> variables;
    variables.resize(1);

    double akdOutOfaeq = (Teq_GeV / Tkd) * pow(spline_hSmooth(log10(Tstar)), 1. / 3.) / pow(spline_hSmooth(log10(Tkd)), 1. / 3.);

    double a = 1./(1+z);

    //std::cout << akdOutOfaeq << " " << spline_hSmooth(log10(Tstar)) << " " << spline_hSmooth(log10(Tkd)) << " " << Teq_GeV  << " " << Tkd  << std::endl;
    //double vkd = VelocityAtKineticDecoupling(Tkd, mDM);

    //std::cout << spline_hSmooth(log10(Tstar)) << std::endl;
    //double C = (1 / a0OutOfaeq) * (sqrt(PI / 3.) / M_PLANCK) * sqrt(rhoEq);

    double Ctau0 = sqrt(1 + Teq_GeV / (T0_CMB_GeV / a) ) - 1;
    double Ctaukd = sqrt(1 + Teq_GeV / Tkd) - 1;
    double factor = log((Ctau0 / Ctaukd) * ((2 + Ctaukd) / (2 + Ctau0)));

    //int nPointsInt = 5000;
    //double err;

    //double Kfs = Simpson_Integral1_NonStatic(0, 1, 0, nPointsInt, log(Tstar), log(Tkd), variables, this, CallBackFToIntToComputeKfs, err);
    double Kfs = GaussLegendre_Integral_Static(0, 1, gl200, log(Tstar), log(Tkd), variables, this, CallBackFToIntToComputeKfs);

    double akd_int_over_time = sqrt(2)*akdOutOfaeq/aeq/Heq*(factor + Kfs);
    //double _kfs =  2*C/(factor + Kfs)/akdOutOfaeq*sqrt(2*mDM/Tkd);// (vkd / (2 * C)) * akdOutOfaeq * (factor + Kfs);
    double _kfs = sqrt(2*mDM/Tkd)/akd_int_over_time;

    return _kfs;
}



double MinimalMass::FToIntToComputeKfs(std::vector<double> variables)
{

    double result;
    double T = exp(variables[0]);

    double prefactor = sqrt(spline_gSmooth(log10(Tstar))) * pow(spline_hSmooth(log10(Tstar)), -2. / 3.) *
                       pow(spline_hSmooth(log10(T)), -1. / 3.);

    result = (prefactor * spline_gStarHalf(log10(T)) - 1);

    return result;
}

std::vector<double> MinimalMass::MassesComputation(double Tkd, double mDM)
{
    double rhoMkd, Hkd;

    double masskd, massfs;

    double kfs_eq = ComovingFreeStreamingScale(Tkd, mDM, 0);

    // Free streaming minimal mass
    //massfs = (4 * PI / 3) * sqrt(PI) * pow(ComovingFreeStreamingLength(Tkd, mDM, 0) / 4., 3) * rhoM0;
    massfs = (4 * PI / 3) * pow(PI/kfs_eq, 3) * rhoM0;


    // Acoustic minimal mass
    rhoMkd = rhoM0 * spline_hSmooth(log10(Tkd)) * pow(Tkd / T0_CMB_GeV, 3) / spline_hSmooth(log10(T0_CMB_GeV));
    Hkd = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(Tkd))) * pow(Tkd, 2) / M_PLANCK;
    masskd = (4 * PI / 3) * (rhoMkd / pow(Hkd, 3)) * pow(3, -3. / 2.);


    std::vector<double> res;
    res.resize(0);
    res.push_back(masskd);
    res.push_back(massfs);

    return res;
}

void MinimalMass::PlotMinimalMass(double mDM)
{
    std::ofstream outfile;
    std::string name_file = "MassHalo_vsTkd_massDM";

    std::ostringstream strs;
    strs << mDM;
    name_file += strs.str();

    outfile.open("../output/MinimalMassFreeStreaming/" + name_file + ".out");

    outfile << "# Temperature [GeV] | Acoustic mass [Msol] | Free streaming mass [Msol]" << std::endl;

    double logTkdMin = log10(0.0001);
    double logTkdMax = log10(10);
    int numPoints = 500;

    double Tkd;

    double delta = (logTkdMax - logTkdMin) / (1.0 * numPoints);

    for (int i = 0; i < numPoints + 1; i++)
    {

        Tkd = pow(10, logTkdMin + delta * i);
        //set_Tkd(Tkd);
        //std::cout << Tkd <<  " " <<  logTkdMin << " " << delta << std::endl;
        outfile << Tkd << "\t" << MassesComputation(Tkd, mDM)[0] * GeV_to_MSOL << " \t " << MassesComputation(Tkd, mDM)[1] * GeV_to_MSOL << std::endl;
        // std::cout << MassesComputation()[0] << " \t " << MassesComputation()[1] << std::endl;
    }

    outfile.close();
}

void MinimalMass::PlotFreeStreamingLength_vs_z(double Tkd, double mDM)
{
    std::ofstream outfile;
    std::string namefile = "../output/MinimalMassFreeStreaming/Free_streaming_length_vs_z_massDM100_Tkd0001.out";

    outfile.open(namefile);
    outfile.precision(8);

    outfile << "# mchi = " << mDM << " GeV | T_kd " << Tkd << " GeV";
    outfile << "# z | lambda/a [pc] | kfs [pc^{-1}]" << std::endl;

    int numPoints = 500;

    double akdOutOfaeq = (Teq_GeV / Tkd) * pow(spline_hSmooth(log10(Teq_GeV)), 1. / 3.) / pow(spline_hSmooth(log10(Tkd)), 1. / 3.);
    double zeq = a0OutOfaeq - 1;
    double zkd = a0OutOfaeq / akdOutOfaeq - 1;
    double zmin = 1, zmax = zkd;
    double dz = log10(zmax / zmin) / (1.0 * numPoints), z = 0; //a = 0;

    std::cout << zeq << " " << zkd << " " << akdOutOfaeq << std::endl;
    double kfs = 0, lfs = 0;

    for (int i = 0; i < numPoints + 1; i++)
    {
        z = zmax * pow(10, -dz * i);
        lfs = ComovingFreeStreamingLength(Tkd, mDM, z);
        kfs = ComovingFreeStreamingScale(Tkd, mDM, z);
        outfile << z << " \t " << lfs / (GeV_to_mm1 * kpc_to_cm * 1e-2) * 1e+3 << "\t" << kfs *(GeV_to_mm1 * kpc_to_cm * 1e-2) / 1e+3  << std::endl;
    }

    outfile.close();
}

void MinimalMass::PlotFreeStreamingLength_vs_Tkd(double z, double mDM)
{
    std::ofstream outfile;
    std::string namefile = "../output/MinimalMassFreeStreaming/Free_streaming_length_vs_Tkd_massDM1000_z0.out";

    outfile.open(namefile);
    outfile.precision(8);

    outfile << "# mchi = " << mDM << " GeV | z=" << z << std::endl;
    outfile << "# Tkd [GeV] | lambda/a [pc] | kfs [pc^{-1}]" << std::endl;

    int numPoints = 500;

    //double akdOutOfaeq = (Teq_GeV / Tkd) * pow(spline_hSmooth(log10(Teq_GeV)), 1. / 3.) / pow(spline_hSmooth(log10(Tkd)), 1. / 3.);
    double Tmin = 1e-4 , Tmax = 10;
    double dT = log10(Tmax / Tmin) / (1.0 * numPoints), Tkd = 0;// a = 0;

    //std::cout << zeq << " " << zkd << " " << akdOutOfaeq << std::endl;
    double kfs = 0, lfs = 0;

    for (int i = 0; i < numPoints + 1; i++)
    {
        Tkd = Tmin * pow(10, i * dT);
        lfs = ComovingFreeStreamingLength(Tkd, mDM, z);
        kfs = ComovingFreeStreamingScale(Tkd, mDM, z);
        outfile << Tkd << " \t " << lfs / (GeV_to_mm1 * kpc_to_cm * 1e-2) * 1e+3 << "\t" << kfs *(GeV_to_mm1 * kpc_to_cm * 1e-2) / 1e+3  << std::endl;
    }

    outfile.close();
}
