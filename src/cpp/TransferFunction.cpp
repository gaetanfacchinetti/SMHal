#include "../headers/TransferFunction.h"


TransferFunction_EH98::TransferFunction_EH98(Cosmology cosmology, bool with_baryons) : TransferFunction()
{
    _cosmology   = cosmology;
    _with_baryons = with_baryons;
    
    Omega_m_h2 = _cosmology.get_Omega_m_h2();
    Omega_b_h2 = _cosmology.get_Omega_b_h2();
    Omega_c_h2 = _cosmology.get_Omega_c_h2();


    // z_eq
    //z_eq = 2.5e4*Omega_m_h2*pow(Theta27,-4)*Z_EQUIVALENCE/3427.91; // to have the Planck normalization
    //z_eq = 2.5e4*Omega_m_h2*pow(Theta27,-4); // Changed here
    z_eq = _cosmology.Zeq_rad_mat();
    //z_eq = 2.5e4*Omega_m_h2;

    //k_eq
    //k_eq = sqrt(2*Omega_m_h2*z_eq)* 1e5/LIGHT_SPEED; // [Mpc^-1]
    k_eq = _cosmology.keq_rad_mat();

    //std::cout << "Transfer function : " << k_eq << " " << z_eq << std::endl;
    
    // z_drag
    double b1 = 0.313*pow(Omega_m_h2,-0.419)*(1+0.607*pow(Omega_m_h2,0.674));
    double b2 = 0.238*pow(Omega_m_h2,0.223);
    //z_drag =  1291.*pow(Omega_m_h2,0.251)/(1+0.659*pow(Omega_m_h2,0.828))*(1+b1*pow(Omega_b_h2,b2))* 1059.94/1020.66;
    z_drag =  1291.*pow(Omega_m_h2,0.251)/(1+0.659*pow(Omega_m_h2,0.828))*(1+b1*pow(Omega_b_h2,b2)); // Changed here

    // sound horizon
    double Rd = R_ratio(z_drag);
    double Req = R_ratio(z_eq);
    sound_horizon = 2./(3*k_eq)*sqrt(6/Req)*log((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req))); //[Mpc]

    // alpha_c
    double a1 = pow(46.9*Omega_m_h2,0.67)*(1+pow(32.1*Omega_m_h2,-0.532));
    double a2 = pow(12*Omega_m_h2,0.424)*(1+pow(45*Omega_m_h2,-0.582));
    alpha_c = pow(a1,-Omega_b_h2/Omega_m_h2)*pow(a2,-pow(Omega_b_h2/Omega_m_h2,3));

    // beta_c
    double b1_2 = 0.944/(1+pow(458*Omega_m_h2,-0.708));
    double b2_2 = pow(0.395*Omega_m_h2,-0.0266);
    beta_c = 1./(1+b1_2*(pow(Omega_c_h2/Omega_m_h2,b2_2)-1));

    // k_silk
    k_silk = 1.6*pow(Omega_b_h2,0.52)*pow(Omega_m_h2,0.73)*( 1+pow(10.4*Omega_m_h2,-0.95) ); // [Mpc^-1]

    // alpha_b
    alpha_b = 2.07*k_eq*sound_horizon*pow(1+Rd,-0.75)*G_func((1+z_eq)/(1+z_drag));

    // beta_b
    beta_b = 0.5+Omega_b_h2/Omega_m_h2+(3-2*Omega_b_h2/Omega_m_h2)*sqrt(1+pow(17.2*Omega_m_h2,2));

}

double TransferFunction_EH98::R_ratio(double z)
{
    return  31.5*Omega_b_h2*pow(Theta27,-4)*1e+3/z;
}


double TransferFunction_EH98::T0_tilde(double q, double alphac, double betac)
{
    double C = 14.2/alphac+386./(1.0+69.9*pow(q,1.08));
    return log(exp(1)+1.8*betac*q)/(log(exp(1)+1.8*betac*q)+C*q*q);

}
double TransferFunction_EH98::shape_parameter(double k)
// k in Mpc^-1
{
     return k*pow(Theta27,2)/Omega_m_h2;
}

double TransferFunction_EH98::T_cdm(double k)
{
    double q = shape_parameter(k);
    double T0_1 = T0_tilde(q,1,beta_c);
    double T0_2 = T0_tilde(q,alpha_c,beta_c);
    double f = 1./(1+pow(k*sound_horizon/5.4,4));

    return f*T0_1+(1-f)*T0_2;
}

/////////////////////////////

double TransferFunction_EH98::G_func(double y)
{
    double res = -6*sqrt(1+y)+(2+3*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1));
    res *= y;
    return res;
}

double TransferFunction_EH98::s_tilde(double k)
// [Mpc]
{
    double beta_node = 8.41*pow(Omega_m_h2,0.435);
    double res = pow(1+pow(beta_node/(k*sound_horizon),3.),-1./3.);
    res *= sound_horizon;
    return res;
}
double TransferFunction_EH98::T_baryons(double k)
{
    double q = shape_parameter(k);
    double sh_t = s_tilde(k);
    double j0 = sinc(k*sh_t);
    double res = T0_tilde(q,1,1)/(1+pow(k*sound_horizon/5.2,2)) + alpha_b/(1+pow(beta_b/(k*sound_horizon),3))*exp(-pow(k/k_silk,1.4));
    res *= j0;
    return res;
}

////////////////////////////

double TransferFunction_EH98::evaluate(double k)
{    
    //std::cout << "Evaluate EH98" << std::endl;

    double Tc = T_cdm(k);
    double Tb = T_baryons(k);
    double res =  (_with_baryons == true) ? Omega_b_h2/Omega_m_h2 * Tb + Omega_c_h2/Omega_m_h2 * Tc : Omega_c_h2/Omega_m_h2 * Tc;

    //double res = Tc;
    return fabs(res);
}

/////////////////////////////

double TransferFunction_EH98::effective_CDM_transfer_function(double k)
// Einsenstein and Hu 1998
// k in Mpc^-1
{
    double q = shape_parameter(k);
    double L0 = log(2*exp(1)+1.8*q);
    double C0 = 14.2+731./(1+62.5*q);
    double res = L0/(L0+C0*q*q);
    return res;
}



void TransferFunction::plot()
{
    std::ofstream outfile;
    outfile.open("../output/TransferFunction.out");

    double kmin = 1e-4;
    double kmax = 1e+10;

    double Npoints = 1000;

    double k, d_logk = log10(kmax / kmin) / Npoints;

    outfile << "# k [Mpc^-1] | T(k) " << std::endl;

    for (int i = 0; i < Npoints + 1; i++)
    {
        k = kmin * pow(10, i * d_logk);
        outfile << k << " " << evaluate(k) << std::endl;
    }

    outfile.close();
}


TransferFunction_WIMP_GHS05::TransferFunction_WIMP_GHS05(Cosmology n_cosmo, double x, double y, bool Tm) : TransferFunction(n_cosmo)
// Here we connect mchi and Tkf to kfs and kd through the slightly modified relation of http://arxiv.org/abs/1604.02457
{
    if (Tm == true)
    {
        double h = _cosmology.get_h();
        double m_adim = x/100.; //  with x = mchi in GeV
        double T_adim = y/30e-3; // with y = Tkd in GeV

        kfs = 2.4e+6*sqrt(m_adim*T_adim)/(1+log(T_adim)/19.2)/h;
        kd = 5.4e+7*sqrt(m_adim*T_adim)/h; 
    }
    else
    {
        kfs = x;
        kd = y;
    }
}

double TransferFunction_WIMP_GHS05::evaluate(double k)
{
    //std::cout << "Evaluate GHS05" << std::endl;
    return (1+(2./3.)*pow(k/kfs, 2))*exp(-pow(k/kfs, 2)-pow(k/kd, 2));
}
