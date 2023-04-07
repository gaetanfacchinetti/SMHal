#include "../headers/FSLModel.h"



FSLModel::FSLModel(double n_m200_min, double n_epst, DarkHalo n_hostHalo) : hostHalo(n_hostHalo), subModel(DarkHalo(n_hostHalo.get_cosmo(), 0, 1, 3, 1))
{
    //hostHalo = n_hostHalo;
    cosmo =  hostHalo.get_cosmo();
    R_max = hostHalo.virial_radius(200, true);

    _N_sub_un = 0;
    is_NsubUnevolved_initialised = false;
    m200_min = n_m200_min;

    //subModel = DarkHalo(cosmo, 0, 1, 3, 1); // NFW model

    epst = n_epst;

    m200_max = hostHalo.virial_mass(200, true);
    c200_max = 500;

    is_rhoSub_initialized = false;
    is_Kt_Tot_initialized = false;
    is_N_sub_initialized = false;
}



double FSLModel::rt_over_rs_DMO(double R, double c200)
{
    std::vector<double> xtRc200;
    xtRc200.resize(3);
    xtRc200[1] = R;
    xtRc200[2] = c200;

    return Dichotomie_Static(0, 3, xtRc200, this, CallBack_funcToSolve_tidal_radius_DMO, 1e-10, 1e+10, 1e-8, 0);
}


double FSLModel::funcToSolve_tidal_radius_DMO(std::vector<double> xtRc200)
{
    double xt = xtRc200[0];
    double R = xtRc200[1];
    double c200 = xtRc200[2];

    double M = hostHalo.dm_mass(R);
    double eta = subModel.rm2_rs();
    double rho_s = (200. * cosmo.critical_density(0) / 3.) * pow(eta * c200, 3) / subModel.mass_profile(eta * c200);

    return subModel.mass_profile(xt) * pow(xt, -3) - ((3 * M) / (4 * PI * rho_s)) * pow(R, -3) * (1 - (4 * PI / 3.) * pow(R, 3) * hostHalo.dm_density(R) / M);
}


void FSLModel::plot_dimensionless_tidal_radius_DMO()
{
    std::ofstream outfile;
    outfile.open("../input/rt_over_rs_core_DMO.in");

    outfile << "0"
            << "    ";

    int N = 500;
    double dlogR = log10(R_max / 0.001) / N;
    double dlogc = log10(500 / 1.) / N;
    double R, c;

    for (int i = 0; i < N + 1; i++)
    {
        c = pow(10, i * dlogc);
        outfile << c;
        if (i != N)
            outfile << "    ";
    }

    outfile << std::endl;

    for (int j = 0; j < N + 1; j++)
    {
        R = 0.001 * pow(10, j * dlogR);
        outfile << R << "    ";

        for (int i = 0; i < N + 1; i++)
        {
            c = pow(10, i * dlogc);
            outfile << rt_over_rs_DMO(R, c);
            if (i != N)
                outfile << "    ";
        }
        outfile << std::endl;
    }
}


void FSLModel::ReadFormTable_rt_over_rs_DMO()
{
    //std::ifstream in_rt_over_rs("../input/rt_over_rs_nfw_DMO.in");
    std::ifstream in_rt_over_rs("../input/rt_over_rs_core_DMO.in");

    std::string line;
    std::string R_str;

    int place_line = 0;
    int line_number = 0;

    if (in_rt_over_rs.is_open())
    {
        while (getline(in_rt_over_rs, line))
        {
            //std::cout << line << '\n';

            if (line[0] != '#')
            {
                place_line = 0;

                std::string delimiter = "    ";

                size_t pos = 0;
                std::string token;

                if (line_number > 0)
                    vec_rt_over_rs_DMO.push_back(std::vector<double>());

                while ((pos = line.find(delimiter)) != std::string::npos)
                {
                    //std::cout << "here oui " << std::endl;
                    token = line.substr(0, pos);

                    //std::cout << pos << std::endl;

                    if (place_line == 0 && line_number > 0)
                        vec_R_DMO.push_back(stod(token));

                    if (place_line > 0 && line_number == 0)
                        vec_c200_DMO.push_back(stod(token));

                    if (place_line > 0 && line_number > 0)
                        vec_rt_over_rs_DMO[line_number - 1].push_back(stod(token));

                    line.erase(0, pos + delimiter.length());

                    place_line++;
                }

                if (place_line > 0 && line_number == 0)
                    vec_c200_DMO.push_back(stod(line));

                if (place_line > 0 && line_number > 0)
                    vec_rt_over_rs_DMO[line_number - 1].push_back(stod(line));

                //std::cout << line << std::endl;

                line_number++;
            }
        }

        in_rt_over_rs.close();
    }
    else
        std::cout << "FATAL ERROR : Impossible to read input file in : " << __PRETTY_FUNCTION__ << std::endl;
}




double FSLModel::Interpolation_rt_over_rs_DMO(double R, double c200)
{

    //std::cout << "Ok ici" << std::endl;
    int pos_R = 0;
    int pos_c = 0;

    while (vec_R_DMO[pos_R] < R)
        pos_R++;


    if (R > vec_R_DMO[vec_R_DMO.size() - 1])
        std::cout << "R to large " << vec_R_DMO[vec_R_DMO.size() - 1] << " " << R << std::endl;

    while (vec_c200_DMO[pos_c] < c200)
        pos_c++;

    if (pos_R > 0)
        pos_R--;
    if (pos_c > 0)
        pos_c--;

    double fQ11 = vec_rt_over_rs_DMO[pos_R][pos_c];
    double fQ12 = vec_rt_over_rs_DMO[pos_R][pos_c + 1];
    double fQ21 = vec_rt_over_rs_DMO[pos_R + 1][pos_c];
    double fQ22 = vec_rt_over_rs_DMO[pos_R + 1][pos_c + 1];

    //std::cout << fQ11 << " " << fQ12 << " " << fQ21 << " " << fQ22 << std::endl;
    double x1 = log10(vec_R_DMO[pos_R]);
    double x2 = log10(vec_R_DMO[pos_R + 1]);
    double y1 = log10(vec_c200_DMO[pos_c]);
    double y2 = log10(vec_c200_DMO[pos_c + 1]);

    //std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << std::endl;

    double Dx = x2 - x1;
    double Dy = y2 - y1;

    double dx = log10(R) - x1;
    double dy = log10(c200) - y1;

    double Delfx = fQ21 - fQ11;
    double Delfy = fQ12 - fQ11;
    double Delfxy = fQ11 + fQ22 - fQ21 - fQ12;

    return Delfx * dx / Dx + Delfy * dy / Dy + Delfxy * dx * dy / (Dx * Dy) + fQ11;
}



double FSLModel::c200_min_DMO(double R, double eps_t)
{
    std::vector<double> Rc200eps;
    Rc200eps.resize(3);
    Rc200eps[0] = R;
    Rc200eps[2] = eps_t;

    double c200_min = Dichotomie_Static(1, 3, Rc200eps, this, CallBack_funcToSolvec200min_DMO, 1, c200_max, 1e-5, 1);

    if (c200_min == 1 && rt_over_rs_DMO(R, c200_max) < eps_t)
        c200_min = c200_max;

    return c200_min;
}

double FSLModel::funcToSolvec200min_DMO(std::vector<double> Rc200eps)
{
    return rt_over_rs_DMO(Rc200eps[0], Rc200eps[1]) - Rc200eps[2];
}

//
//
//
//
//
//

double FSLModel::ProbDistributionSpace(double R)
{
    double M200 =  m200_max;
    double rho = hostHalo.dm_density(R);

    //std::ofstream outfile;
    //outfile.open("../output/ProbDistributionSpace_out", std::ios_base::app);
    //outfile << R << " " << rho/M200 << std::endl;

    return rho / M200; // en kpc^(-3)
}



double FSLModel::MassFunction(double m200, double M200)
{
    /*
    double gamma1 = 0.014;
    double alpha1 = -0.965;
    double gamma2 = 0.41;
    double alpha2 = -0.57;
    double beta = 20;
    double zeta = 3.4;*/
    
    
    double gamma1 = 0.019;
    double alpha1 = -0.94;
    double gamma2 = 0.464;
    double alpha2 = -0.58;
    double beta = 24;
    double zeta = 3.4;

    /*double gamma1 = 0.12090712;
    double alpha1 = -0.99257559;
    double gamma2 = 0.12197409;
    double alpha2 = -0.82846965;
    double beta = 94.75905245;
    double zeta =  5.40502835;*/

    double m200_M200 = m200/M200;

    return (gamma1*pow(m200_M200, alpha1) + gamma2*pow(m200_M200, alpha2))*exp(-beta*pow(m200_M200, zeta))/m200;
}

double FSLModel::TotalNumberOfSubhalosUnevolved()
{
    /*
    double gamma1 = 0.014;
    double alpha1 = -0.965;
    double gamma2 = 0.41;
    double alpha2 = -0.57;
    double beta = 20;
    double zeta = 3.4;*/


    double gamma1 = 0.019;
    double alpha1 = -0.94;
    double gamma2 = 0.464;
    double alpha2 = -0.58;
    double beta = 24;
    double zeta = 3.4;

    
    /*double gamma1 = 0.12090712;
    double alpha1 = -0.99257559;
    double gamma2 = 0.12197409;
    double alpha2 = -0.82846965;
    double beta = 94.75905245;
    double zeta =  5.40502835;*/
    
    double M200 = m200_max;

    double xlim_1 = 1e-3*pow(gamma1/gamma2, 1./(alpha2 - alpha1));
    double xlim_2 = pow(1e-3/beta, 1/zeta);
    double xlim = myMin(xlim_1, xlim_2);
    double res = gamma1/(-alpha1)*(pow(m200_min/M200, alpha1) - pow(xlim, alpha1));
    
    //std::cout << "xlim : " << xlim << " " << xlim_1 << " " << xlim_2 << " " << res << " " << m200_min/M200 << std::endl;

    std::vector<double> xx ={0, M200};
    
    if(M200/100. > M200*xlim)
        res += GaussLegendre_IntegralLn_Static(0, 2, gl100, M200*xlim, M200/100., xx, this, CallBack_MassFunction);
    
    res += GaussLegendre_IntegralLn_Static(0, 2, gl100, M200/100., M200, xx, this, CallBack_MassFunction);
    return res;
}


double FSLModel::ProbDistributionMass(double m200)
{
    double M200 = m200_max;
    return MassFunction(m200, M200)/NsubUnevolved();
}


double FSLModel::ProbDistributionConcentration(double c200, double m200)
// m200 in Msol
{
    // Definition of properties from Sanchez-Conde & Prada relation
    double sigma_c = 0.14 * log(10);
    double cbar = c_bar(m200);

    double Kc = 0.5 * erfc(-log(cbar) / (sqrt(2) * sigma_c));

    return 1. / Kc / c200 / sqrt(2 * PI) / sigma_c * exp(-pow((log(c200) - log(cbar)) / sqrt(2) / sigma_c, 2));
}

double FSLModel::c_bar(double m200)
// m200 in Msol
{
    double cbar = 0;

    double cn[6] = {37.5153, -1.5093, 1.636e-2, 3.66e-4, -2.89237e-5, 5.32e-7};

    double m200_bis = m200;

    if (m200 <= 7.24e-10)
        m200_bis = 7.24e-10;

    for (int i = 0; i < 6; i++)
        cbar += cn[i] * pow(log(m200_bis * cosmo.get_h()), i);

    return cbar;
}



double FSLModel::fToIntOnm200_MassFraction(std::vector<double> yy)
{
    double m200 = yy[1];
    double R = yy[2];

    double c200min = c200_min_DMO(R, 1);

    double sigma_c = 0.14 * log(10);
    double cbar = c_bar(m200);

    double Kc = 0.5 * erfc(-log(cbar) / (sqrt(2) * sigma_c));

    //c200min = 1;

    if (c200min < c200_max)
        //return m200 * ProbDistributionMass(m200) * (-0.5 * erf(-(log(c200_max) - log(cbar)) / (sqrt(2) * sigma_c)) + 0.5 * erf(-(log(c200min) - log(cbar)) / (sqrt(2) * sigma_c))) / Kc;
        return m200 * ProbDistributionMass(m200) * (0.5 * my_erf_diff((log(c200_max) - log(cbar)) / (sqrt(2) * sigma_c), (log(c200min) - log(cbar)) / (sqrt(2) * sigma_c))) / Kc;
    else
        return 0;
}

double FSLModel::fToIntOnR_MassFraction(std::vector<double> yy)
{
    double R = yy[2];

    double M200 =  m200_max;
    double m1 = 2.2e-6 * M200;
    double m2 = 8.8e-4 * M200;
    return R * R * ProbDistributionSpace(R) * GaussLegendre_IntegralLn_Static(1, 3, gl1000, m1, m2, yy, this, CallBack_fToIntOnm200_MassFraction);
}

double FSLModel::MassFraction()
{
    std::vector<double> yy;
    yy.resize(3);

    double M200 =   m200_max;

    double integral = GaussLegendre_IntegralLn_Static(2, 3, gl100, 0, R_max, yy, this, CallBack_fToIntOnR_MassFraction);
    return 4*PI*NsubUnevolved()*integral/M200;
}




double FSLModel::fToIntOnc200_Norm_Tot(std::vector<double> varNorm)
{
    double c200 = varNorm[0];
    double m200 = varNorm[1];

    return ProbDistributionConcentration(c200, m200);
}

double FSLModel::fToIntOnm200_Norm_Tot(std::vector<double> varNorm)
{
    double R = varNorm[2];
    double m200 = varNorm[1];

    double c200min = c200_min_DMO(R, epst);
    //double cmax = c_bar(m200) * exp(-pow(0.14 * log(10), 2));

    double sigma_c = 0.14 * log(10);
    double cbar = c_bar(m200);

    double Kc = 0.5 * erfc(-log(cbar) / (sqrt(2) * sigma_c));

    //return ProbDistributionMass(m200)*GaussLegendre_IntegralLn_Static(0, 3, gl200, c200min, c200_max, varNorm, this, CallBack_fToIntOnc200_Norm_Tot);
    //return ProbDistributionMass(m200) * (-0.5 * erf(-(log(c200_max) - log(cbar)) / (sqrt(2) * sigma_c)) + 0.5 * erf(-(log(c200min) - log(cbar)) / (sqrt(2) * sigma_c))) / Kc;

    //return ProbDistributionMass(m200);
    //std::cout << "Coucou here " << m200 << std::endl;

    if (c200min < c200_max)
        return ProbDistributionMass(m200) * (0.5 * my_erf_diff((log(c200_max) - log(cbar)) / (sqrt(2) * sigma_c), (log(c200min) - log(cbar)) / (sqrt(2) * sigma_c))) / Kc;
    else
        return 0;
}

double FSLModel::fToIntOnR_Norm_Tot(std::vector<double> varNorm)
{
    //std::cout << GaussLegendre_IntegralLn_Static(1, 3, gl200, m200_min, m200_max, varNorm, this, CallBack_fToIntOnm200_Norm_Tot) << std::endl;
  
    double R = varNorm[2];
    //return R*R*hostHalo.dm_density(R);
    //return R * R * ProbDistributionSpace(R);
    //std::cout << "Coucou here " << R << std::endl;
    return R * R * ProbDistributionSpace(R) * GaussLegendre_IntegralLn_Static(1, 3, gl200, m200_min, m200_max, varNorm, this, CallBack_fToIntOnm200_Norm_Tot);
    //return GaussLegendre_IntegralLn_Static(1, 3, gl1000, m200_min, m200_max, varNorm, this, CallBack_fToIntOnm200_Norm_Tot);
}

double FSLModel::Normalisation_Kt_Tot()
{
    std::vector<double> varNorm;
    varNorm.resize(3);

    //std::cout << hostHalo.virial_mass(200, true) << std::endl;

    //std::cout << hostHalo.virial_radius(200, true) << std::endl;
    //std::cout << "here " << R_max << " " << 4 * PI * GaussLegendre_IntegralLn_Static(2, 3, gl100, 0, R_max, varNorm, this, CallBack_fToIntOnR_Norm_Tot) << std::endl;

    double res =  4 * PI * GaussLegendre_IntegralLn_Static(2, 3, gl200, 1e-3, R_max, varNorm, this, CallBack_fToIntOnR_Norm_Tot);

    return res;
}



double FSLModel::NumberDensity_subhalos(double R)
{
    std::vector<double> varNorm={0,0,R};
    double Nsub = get_Nsub();
    double myKt = get_Kt_Tot();
    //return fToIntOnR_Norm_Tot(varNorm);
    return Nsub/myKt*fToIntOnR_Norm_Tot(varNorm)/R/R;
}

void FSLModel::plot_NumberDensity_subhalos(int id)
{
    std::ofstream outfile;
    outfile.open("../output/NumberClumps/NumberDensity_subhalos" + std::to_string(id) + "_v2.out");
    outfile << "# Number density of subhalos inside at distance R from the center of the galaxy" << std::endl;
    outfile << "# alpha_m = " << alpha_m << " | epst = " << epst << " | m200_min = " << m200_min << " [Msol]" << std::endl;
    outfile << "# R [kpc] | Number_density [pc^3] | Number_density without tides [pc^{-3}]" << std::endl;

    int N = 100;

    double Rmin = 1e-1, Rmax = R_max;
    double dR = log10(Rmax / Rmin) / N, R;
    double Num, Nsub = get_Nsub();

    for (int i = 0; i < N + 1; i++)
    {
        R = Rmin * pow(10, i * dR);
        Num = NumberDensity_subhalos(R);
        outfile << R << "\t" << Num*1e-9 << "\t" <<  Nsub*ProbDistributionSpace(R)*1e-9 << std::endl;
    }

    outfile.close();
}




std::vector<double> FSLModel::Nsub_calibrated()
{
    double Kt = Normalisation_Kt_Tot();
    double NsimKtsim = Nsim_Ktsim();

    //std::cout << NsimKtsim << std::endl;

    std::vector<double> Nsub_Kt;
    Nsub_Kt.resize(2);
    Nsub_Kt[0] = Kt * NsimKtsim;
    Nsub_Kt[1] = Kt;

    _Kt_Tot = Kt;
    _N_sub = Kt * NsimKtsim;

    is_Kt_Tot_initialized = true;
    is_N_sub_initialized = true;

    return Nsub_Kt;
}

double FSLModel::fToIntOnc200AverageRhoSub(std::vector<double> yy)
{
    double c200 = yy[0];
    double m200 = yy[1];
    double R = yy[2];

    double eta = subModel.rm2_rs();
    double xt = rt_over_rs_DMO(R, c200);
    //double xt = 1;
    //double rs = subModel.scale_radius(c200, m200, 200, true);

    double mass = m200 * subModel.mass_profile(xt) * pow(subModel.mass_profile(eta * c200), -1);
    return mass * ProbDistributionConcentration(c200, m200);
}

double FSLModel::fToIntOnm200AverageRhoSub(std::vector<double> yy)
{
    double m200 = yy[1];
    double R = yy[2];

    double c200min = c200_min_DMO(R, epst);
    double cmax = c_bar(m200) * exp(-pow(0.14 * log(10), 2));

    //std::cout << "here" << std::endl;

    if (c200min < c200_max && c200min < cmax)
        return ProbDistributionMass(m200) * (GaussLegendre_Integral_Static(0, 3, gl200, c200min, cmax, yy, this, CallBack_fToIntOnc200AverageRhoSub) + GaussLegendre_IntegralLn_Static(0, 3, gl200, cmax, c200_max, yy, this, CallBack_fToIntOnc200AverageRhoSub));
    else if (c200min < c200_max && c200min >= cmax)
        return ProbDistributionMass(m200) * GaussLegendre_IntegralLn_Static(0, 3, gl200, c200min, c200_max, yy, this, CallBack_fToIntOnc200AverageRhoSub);
    else
        return 0;
}

double FSLModel::UnNormalizedAverageRhoSub(double R)
{
    std::vector<double> varNorm;
    varNorm.resize(3);

    varNorm[2] = R;

    return ProbDistributionSpace(R) * GaussLegendre_IntegralLn_Static(1, 3, gl200, m200_min, m200_max, varNorm, this, CallBack_fToIntOnm200AverageRhoSub);
}

int FSLModel::InitialiseAverageRhoSub()
{
    int N_R = 200;
    std::vector<double> R;
    std::vector<double> Rho_R;

    std::ofstream out_RhoSub;
    out_RhoSub.open("../output/RhoSub_vs_R.out");
    out_RhoSub.precision(8);
    out_RhoSub << "# R [kpc] | rho [pc^{—3}] | rho_tot gal [pc^{-3}]" << std::endl;

    R.resize(N_R + 1);
    Rho_R.resize(N_R + 1);

    double R_min_bis = 1e-2;

    double d_logR = (log10(R_max) - log10(R_min_bis)) / N_R;

    double NsubOverKt = NsubUnevolved();

    for (int i = 0; i < N_R + 1; i++)
    {
        
        R[i] = R_min_bis * pow(10, i * d_logR);
        //std::cout << "In there " << R[i] << std::endl;
        
        Rho_R[i] = NsubOverKt * UnNormalizedAverageRhoSub(R[i]);
        //out_RhoSub << R[i] << "\t" << Rho_R[i]*1e-9 << "\t" << hostHalo.dm_density(R[i])*1e-9 << std::endl;
    }

    spline_AverageRhoSub.set_points(R, Rho_R, true);

    out_RhoSub.close();

    is_rhoSub_initialized = true;

    return 0;
}





double FSLModel::fToIntOnc200_EvolvedMassFunction(double c200, double mt, double R, double c200_m)
{

    double xt = rt_over_rs_DMO(R, c200);
    //double xt = Interpolation_rt_over_rs(R, c200);
    double eta = subModel.rm2_rs();

    double m200 = mt * subModel.mass_profile(eta * c200) / subModel.mass_profile(xt);

    if (m200 > m200_max) // Means that the subhalo is destoryed
        return 0;

    // The integral is done originally over m200 and the dirac distribution select a tidal mass mt
    // therefore we should not forget the conversion factor between mt and m200 (that comes from a change of variable
    // to recover an integral over mt from the integral over m200)

    if (c200 < c200_m)
        return 0;

    return ProbDistributionConcentration(c200, m200) * ProbDistributionMass(m200) * ProbDistributionSpace(R) * subModel.mass_profile(eta * c200) / subModel.mass_profile(xt);
}

double FSLModel::EvolvedMassFunction_Rfixed(double mt, double R, double c200_m)
{
    std::vector<double> xx = {0, mt, R, c200_m};
    return NsubUnevolved()*GaussLegendre_IntegralLn_Static(0, 4, gl100, c200_m, c200_max, xx, this, CallBack_fToIntOnc200_EvolvedMassFunction);
}

double FSLModel::fToIntOnR_EvolvedMassFunction(double mt, double R)
{
    double c200_m = c200_min_DMO(R, epst);
    return 4*PI*R*R*EvolvedMassFunction_Rfixed(mt, R, c200_m);
}

double FSLModel::EvolvedMassFunction(double mt)
{
    std::vector<double> xx = {mt, 0};
    return GaussLegendre_IntegralLn_Static(1, 2, gl100, 1e-3, R_max, xx, this, CallBack_fToIntOnR_EvolvedMassFunction);
}


void FSLModel::plot_EvolveMassFunction(int version)
{
    int Npts = 200;

    std::ofstream outfile;
    outfile.open("../output/EvolvedMassFunction/EvolvedMassFunction_AndoComp_v" + std::to_string(version) + ".out");
    outfile << "# m/M200 | m*EvolvedMassFunction | m*MassFunction" << std::endl;

    double M200 =  m200_max;
    double dm = log10(1./(m200_min/M200))/Npts, m;

    for (int i = 0; i < Npts + 1; i++)
    {
        m = m200_min*pow(10, i*dm);
        outfile << m << "\t" << m*m*EvolvedMassFunction(m) << "\t" << m*m*MassFunction(m, M200) << std::endl;
    }

    outfile.close();

}
