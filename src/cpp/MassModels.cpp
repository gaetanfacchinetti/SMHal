#include "../headers/MassModels.h"

//#include "dark_halo.h"
//#include "baryons.h"

using namespace std;


MassModel::MassModel(){}

MassModel::MassModel(int model_dm, int model_baryons)
{
  m_model_dm = model_dm;
  m_model_baryons = model_baryons;
}


double MassModel::dm_alpha() const
{
  double res = 0;
  switch(m_model_dm)
    {
    case 0:
      res = 1; break;
    case 1:
      res = 1; break;
    case 2:
      res = 1; break;
    case 3:
      res = 1; break;
    case 4:
      res = 1; break;
    case 5:
      res = 2.639; break;
    case 6:
      res = 2.165; break;
    case 7:
      res = 1.768; break;
    case 8:
      res = 1.407; break;
    case 9:
      res = 1; break;
    case 10:
      res = 1; break;
    case 11:
      res = 1; break;
    case 12:
      res = 1; break;
    }
  return res;
}
double MassModel::dm_beta() const
{
  double res = 0;
  switch(m_model_dm)
    {
    case 0:
      res = 3; break;
    case 1:
      res = 3; break;
    case 2:
      res = 3; break;
    case 3:
      res = 3; break;
    case 4:
      res = 3; break;
    case 5:
      res = 2.621; break;
    case 6:
      res = 2.504; break;
    case 7:
      res = 2.633; break;
    case 8:
      res = 2.756; break;
    case 9:
      res = 2.785; break;
    case 10:
      res = 2.805; break;
    case 11:
      res = 2.721; break;
    case 12:
      res = 2.744; break;
    }
  return res;
}
double MassModel::dm_gamma() const
{
  double res = 0;
  switch(m_model_dm)
    {
    case 0:
      res = 1; break;
    case 1:
      res = 0.79; break;
    case 2:
      res = 0.5; break;
    case 3:
      res = 0.25; break;
    case 4:
      res = 0; break;
    case 5:
      res = 0.127; break;
    case 6:
      res = 0.189; break;
    case 7:
      res = 0.468; break;
    case 8:
      res = 0.090; break;
    case 9:
      res = 1.101; break;
    case 10:
      res = 1.065; break;
    case 11:
      res = 1.066; break;
    case 12:
      res = 0.975; break;
    }
  return res;
}
double MassModel::dm_scale_density() const
{
  double res = 0;
  switch(m_model_dm)
    {
    case 0:
      res = 8.51723e+06; break;
    case 1:
      res = 1.57691e+07; break;
    case 2:
      res = 3.19009e+07; break;
    case 3:
      res = 5.26136e+07; break;
    case 4:
      res = 9.08606e+07; break;
    case 5:
      res = pow(10,7.622); break;
    case 6:
      res = pow(10,7.746); break;
    case 7:
      res = pow(10,7.414); break;
    case 8:
      res = pow(10,7.667); break;
    case 9:
      res = pow(10,7.101); break;
    case 10:
      res = pow(10,6.848); break;
    case 11:
      res = pow(10,6.963); break;
    case 12:
      res = pow(10,7.181); break;
    }
  return res;
}

double MassModel::dm_scale_radius() const
{
  double res = 0;
  switch(m_model_dm)
    {
    case 0:
      res = 18.6; break;
    case 1:
      res = 15.4; break;
    case 2:
      res = 11.7; break;
    case 3:
      res = 9.6; break;
    case 4:
      res = 7.7; break;
    case 5:
      res = 4.885; break;
    case 6:
      res = 4.210; break;
    case 7:
      res = 7.364; break;
    case 8:
      res = 6.671; break;
    case 9:
      res = 10.338; break;
    case 10:
      res = 14.291; break;
    case 11:
      res = 13.786; break;
    case 12:
      res = 11.288; break;
    }
  return res;
}

// thick stellar disc

double MassModel::thick_disc_Sigma_d() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 183e6; break;
    case 1:
      res = pow(10,7.885)*2*1.219; break;
    case 2:
      res = pow(10,7.860)*2*0.674; break;
    case 3:
      res = pow(10,8.892)*2*1.067; break;
    case 4:
      res = pow(10,8.493)*2*0.900; break;
    case 5: // mass of MN disc
      res = 1.6e11; break;
    }
  return res;
}
double MassModel::thick_disc_Rd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 3.02; break;
    case 1:
      res = 6.159; break;
    case 2:
      res = 8.671; break;
    case 3:
      res = 3.017; break;
    case 4:
      res = 4.864; break;
    case 5:
      res = 6.5; break;
    }
  return res;
}
double MassModel::thick_disc_zd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 0.9; break;
    case 1:
      res = 1.219; break;
    case 2:
      res = 0.674; break;
    case 3:
      res = 1.067; break;
    case 4:
      res = 0.900; break;
    case 5:
      res = 0.26; break;
    }
  return res;
}

// thin stellar disc

double MassModel::thin_disc_Sigma_d() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 896e6; break;
    case 1:
      res = pow(10,8.634)*2*0.608; break;
    case 2:
      res = pow(10,9.080)*2*0.021; break;
    case 3:
      res = pow(10,7.384)*2*1.049; break;
    case 4:
      res = pow(10,9.303)*2*0.787; break;
    }
  return res;
}
double MassModel::thin_disc_Rd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 2.5; break;
    case 1:
      res = 3.104; break;
    case 2:
      res = 1.622; break;
    case 3:
      res = 2.098; break;
    case 4:
      res = 0.892; break;
    }
  return res;
}
double MassModel::thin_disc_zd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 0.3; break;
    case 1:
      res = 0.608; break;
    case 2:
      res = 0.021; break;
    case 3:
      res = 1.049; break;
    case 4:
      res = 0.787; break;
    }
  return res;
}

// H1 disc McMillan 2017

double MassModel::h1_disc_Sigma_d() const
// [M_sun kpc^{-2}]
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 53.1e6; break;
    }
  return res;
}
double MassModel::h1_disc_Rm() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 4; break;
    }
  return res;
}
double MassModel::h1_disc_Rd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 7; break;
    }
  return res;
}
double MassModel::h1_disc_zd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 8.5e-2; break;
    }
  return res;
}

// H2 disc McMillan 2017

double MassModel::h2_disc_Sigma_d() const
// [M_sun kpc^{-2}]
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 2180e6; break;
    }
  return res;
}
double MassModel::h2_disc_Rm() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 12; break;
    }
  return res;
}
double MassModel::h2_disc_Rd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 1.5; break;
    }
  return res;
}
double MassModel::h2_disc_zd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 4.5e-2; break;
    }
  return res;
}

// gas disc simulations

double MassModel::gas_disc_Sigma_d() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 1:
      res = pow(10,9.037)*2*0.031; break;
    case 2:
      res = pow(10,9.327)*2*0.135; break;
    case 3:
      res = pow(10,8.497)*2*0.119; break;
    case 4:
      res = pow(10,8.415)*2*0.131; break;
    }
  return res;
}

double MassModel::gas_disc_Rd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 1:
      res = 6.116; break;
    case 2:
      res = 1.223; break;
    case 3:
      res = 5.971; break;
    case 4:
      res = 7.510; break;
    }
  return res;
}

double MassModel::gas_disc_zd() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 1:
      res = 0.031; break;
    case 2:
      res = 0.135; break;
    case 3:
      res = 0.119; break;
    case 4:
      res = 0.131; break;
    }
  return res;
}

// Bulge

double MassModel::bulge_rho_b() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 98.4e9; break;
    case 1:
      res = pow(10,8.095); break;
    case 2:
      res = pow(10,9.581); break;
    case 3:
      res = pow(10,9.381); break;
    case 4:
      res = pow(10,9.559); break;
    }
  return res;
}
double MassModel::bulge_q() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 0.5; break;
    case 1:
      res = 0.540; break;
    case 2:
      res = 0.330; break;
    case 3:
      res = 0.677; break;
    case 4:
      res = 0.984; break;
    }
  return res;
}
double MassModel::bulge_r0() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 0.075; break;
    case 1:
      res = 0.340; break;
    case 2:
      res = 0.129; break;
    case 3:
      res = 0.678; break;
    case 4:
      res = 0.641; break;
    }
  return res;
}
double MassModel::bulge_r_cut() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 2.1; break;
    case 1:
      res = 2.297; break;
    case 2:
      res = 4.896; break;
    case 3:
      res = 1.063; break;
    case 4:
      res = 0.618; break;
    }
  return res;
}
double MassModel::bulge_alpha_b() const
{
  double res = 0;
  switch(m_model_baryons)
    {
    case 0:
      res = 1.8; break;
    case 1:
      res = 0.260; break;
    case 2:
      res = 0.818; break;
    case 3:
      res = 0.381; break;
    case 4:
      res = 3.714; break;
    }
  return res;
}

// radial boundary

double MassModel::R_max() const
{
  double res = 0;
  if(m_model_dm==0 || m_model_dm==1 || m_model_dm==2 || m_model_dm==3 || m_model_dm==4)
    {res = 500.;}
  else if(m_model_dm==5)
    {res = 881.61;}
  else if(m_model_dm==6)
    {res = 1725.09;}
  else if(m_model_dm==7)
    {res = 2908.43;}
  else if(m_model_dm==8)
    {res = 1453.96;}
  else if(m_model_dm==9)
    {res = 785.94;}
  else if(m_model_dm==10)
    {res = 1688.46;}
  else if(m_model_dm==11)
    {res = 2792.98;}
  else if(m_model_dm==12)
    {res = 1418.67;}
  else
    {cout << "Rmax is not defined for this mass model!" << endl;}
  return res;
}

double MassModel::bulge_density_mcmillan(double R, double z, double q, double r0,
				       double r_cut, double rho_b, double alpha_b)
{
  double r = sqrt(R*R+z*z/(q*q));
  return rho_b*pow(1+r/r0,-alpha_b)*exp(-r*r/(r_cut*r_cut));
}

double MassModel::exponential_disc_density(double R, double z, double Rd, double zd, double Sigma_d)
{return 0.5*Sigma_d/zd*exp(-R/Rd-fabs(z)/zd);}

double MassModel::gas_disc_density_mcmillan(double R, double z, double Rd, double Rm, double zd, double Sigma_d)
{return 0.25*Sigma_d/zd*exp(-Rm/R-R/Rd)*pow(cosh(0.5*z/zd),-2);}

double MassModel::exponential_disc_surface_density(double R, double Rd, double Sigma_d)
{return Sigma_d*exp(-R/Rd);}

double MassModel::gas_disc_surface_density_mcmillan(double R, double Rd, double Rm, double Sigma_d)
{return Sigma_d*exp(-Rm/R-R/Rd);}

double MassModel::miyamoto_nagai_disc_density(double R, double z, double a, double b, double M)
{
  double z0 = sqrt(z*z+b*b);
  double res = a*R*R+(a+3*z0)*pow(a+z0,2);
  res *= pow(R*R+pow(a+z0,2),-2.5)*pow(z0,-3);
  res *= b*b*M/(4*pi);
  return res;
}

double MassModel::miyamoto_nagai_disc_potential(double R, double z, double a, double b, double M) // [(km/s)^]
{
  double z0 = sqrt(z*z+b*b);
  double res = -G_NEWTON*M/sqrt(R*R+pow(a+z0,2));
  return res*pow(speed_in_kms,2);
}
