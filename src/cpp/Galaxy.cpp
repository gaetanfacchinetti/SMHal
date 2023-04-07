#include "../headers/Galaxy.h"

using namespace std;

Galaxy::Galaxy(Cosmology cosmo) : DarkHalo(cosmo)
{
  m_model_dm = 0; // NFW profile, (alpha, beta, gamma) = (1, 3, 1)
  m_model_baryons = 0;
  m_profile = 0; // NFW-like profile
  model = MassModel(m_model_dm, m_model_baryons);
  update_rhos(model.dm_scale_density());
  update_rs(model.dm_scale_radius());
  update_alpha(model.dm_alpha());
  update_beta(model.dm_beta());
  update_gamma(model.dm_gamma());
}

Galaxy::Galaxy(Cosmology cosmo, int model_dm, int model_baryons) : DarkHalo(cosmo)
{
  m_model_dm = model_dm;
  m_model_baryons = model_baryons;
  m_profile = 0; // NFW-like profile
  model = MassModel(m_model_dm, m_model_baryons);
  update_rhos(model.dm_scale_density());
  update_rs(model.dm_scale_radius());
  update_alpha(model.dm_alpha());
  update_beta(model.dm_beta());
  update_gamma(model.dm_gamma());
}

double Galaxy::get_rmax() const
{
  return model.R_max();
}

double Galaxy::total_baryon_density(double R, double z)
{
  double Rd_thick, zd_thick, Sigmad_thick, thick_disc;
  double Rd_thin, zd_thin, Sigmad_thin, thin_disc;
  double Rd_h1, Rm_h1, zd_h1, Sigma_h1, h1_disc;
  double Rd_h2, Rm_h2, zd_h2, Sigma_h2, h2_disc;
  double q, r0, r_cut, rho_b, alpha_b, bulge;
  double Sigma_gas, Rd_gas, zd_gas, gas_disc;
  double res = 0;
  if (m_model_baryons == 0) // McMillan 2017
  {
    Rd_thick = model.thick_disc_Rd();
    zd_thick = model.thick_disc_zd();
    Sigmad_thick = model.thick_disc_Sigma_d();
    thick_disc = model.exponential_disc_density(R, z, Rd_thick, zd_thick, Sigmad_thick);
    Rd_thin = model.thin_disc_Rd();
    zd_thin = model.thin_disc_zd();
    Sigmad_thin = model.thin_disc_Sigma_d();
    thin_disc = model.exponential_disc_density(R, z, Rd_thin, zd_thin, Sigmad_thin);
    Rd_h1 = model.h1_disc_Rd();
    Rm_h1 = model.h1_disc_Rm();
    zd_h1 = model.h1_disc_zd();
    Sigma_h1 = model.h1_disc_Sigma_d();
    h1_disc = model.gas_disc_density_mcmillan(R, z, Rd_h1, Rm_h1, zd_h1, Sigma_h1);
    Rd_h2 = model.h2_disc_Rd();
    Rm_h2 = model.h2_disc_Rm();
    zd_h2 = model.h2_disc_zd();
    Sigma_h2 = model.h2_disc_Sigma_d();
    h2_disc = model.gas_disc_density_mcmillan(R, z, Rd_h2, Rm_h2, zd_h2, Sigma_h2);
    q = model.bulge_q();
    r0 = model.bulge_r0();
    r_cut = model.bulge_r_cut();
    rho_b = model.bulge_rho_b();
    alpha_b = model.bulge_alpha_b();
    bulge = model.bulge_density_mcmillan(R, z, q, r0, r_cut, rho_b, alpha_b);
    res = thick_disc + thin_disc + h1_disc + h2_disc + bulge;
  }
  else if (m_model_baryons == 1 || m_model_baryons == 2 || m_model_baryons == 3 || m_model_baryons == 4) // Simulations
  {
    Rd_thick = model.thick_disc_Rd();
    zd_thick = model.thick_disc_zd();
    Sigmad_thick = model.thick_disc_Sigma_d();
    thick_disc = model.exponential_disc_density(R, z, Rd_thick, zd_thick, Sigmad_thick);
    Rd_thin = model.thin_disc_Rd();
    zd_thin = model.thin_disc_zd();
    Sigmad_thin = model.thin_disc_Sigma_d();
    thin_disc = model.exponential_disc_density(R, z, Rd_thin, zd_thin, Sigmad_thin);
    q = model.bulge_q();
    r0 = model.bulge_r0();
    r_cut = model.bulge_r_cut();
    rho_b = model.bulge_rho_b();
    alpha_b = model.bulge_alpha_b();
    bulge = model.bulge_density_mcmillan(R, z, q, r0, r_cut, rho_b, alpha_b);
    Sigma_gas = model.gas_disc_Sigma_d();
    Rd_gas = model.gas_disc_Rd();
    zd_gas = model.gas_disc_zd();
    gas_disc = model.exponential_disc_density(R, z, Rd_gas, zd_gas, Sigma_gas);
    res = thick_disc + thin_disc + gas_disc + bulge;
  }
  else if (m_model_baryons == 5) // Miyamoto-Nagai disc
  {
    Rd_thick = model.thick_disc_Rd();
    zd_thick = model.thick_disc_zd();
    Sigmad_thick = model.thick_disc_Sigma_d();
    thick_disc = model.miyamoto_nagai_disc_density(R, z, Rd_thick, zd_thick, Sigmad_thick);
    res = thick_disc;
  }
  return res;
}

double Galaxy::total_ISM_density(double R, double z)
{
  double Rd_h1, Rm_h1, zd_h1, Sigma_h1, h1_disc;
  double Rd_h2, Rm_h2, zd_h2, Sigma_h2, h2_disc;

  double res = 0;

  if (m_model_baryons == 0) // McMillan 2017
  {
    Rd_h1 = model.h1_disc_Rd();
    Rm_h1 = model.h1_disc_Rm();
    zd_h1 = model.h1_disc_zd();
    Sigma_h1 = model.h1_disc_Sigma_d();
    h1_disc = model.gas_disc_density_mcmillan(R, z, Rd_h1, Rm_h1, zd_h1, Sigma_h1);
    Rd_h2 = model.h2_disc_Rd();
    Rm_h2 = model.h2_disc_Rm();
    zd_h2 = model.h2_disc_zd();
    Sigma_h2 = model.h2_disc_Sigma_d();
    h2_disc = model.gas_disc_density_mcmillan(R, z, Rd_h2, Rm_h2, zd_h2, Sigma_h2);

    res = h1_disc + h2_disc;
  }
  else
  {
    std::cout << "Model not implemented in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return res;
}

double Galaxy::HI_density(double R, double z)
{
  double Rd_h1, Rm_h1, zd_h1, Sigma_h1, h1_disc;
  double res = 0;

  if (m_model_baryons == 0) // McMillan 2017
  {
    Rd_h1 = model.h1_disc_Rd();
    Rm_h1 = model.h1_disc_Rm();
    zd_h1 = model.h1_disc_zd();
    Sigma_h1 = model.h1_disc_Sigma_d();
    h1_disc = model.gas_disc_density_mcmillan(R, z, Rd_h1, Rm_h1, zd_h1, Sigma_h1);
    res = h1_disc;
  }
  else
  {
    std::cout << "Model not implemented in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return res;
}

double Galaxy::der_R_HI_density(double R, double z)
{
  double Rd_h1, Rm_h1;
  double res = 0;

  if (m_model_baryons == 0) // McMillan 2017
  {
    Rd_h1 = model.h1_disc_Rd();
    Rm_h1 = model.h1_disc_Rm();
    res = HI_density(R, z) * (-1 / Rd_h1 + Rm_h1 / (R * R));
  }
  else
  {
    std::cout << "Model not implemented in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return res;
}

double Galaxy::der_z_HI_density(double R, double z)
{
  double zd_h1;
  double res = 0;

  if (m_model_baryons == 0) // McMillan 2017
  {
    zd_h1 = model.h1_disc_zd();
    res = -HI_density(R, z) * tanh(z / 2 / zd_h1) / zd_h1;
  }
  else
  {
    std::cout << "Model not implemented in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return res;
}

double Galaxy::H2_density(double R, double z)
{
  double Rd_h2, Rm_h2, zd_h2, Sigma_h2, h2_disc;
  double res = 0;

  if (m_model_baryons == 0) // McMillan 2017
  {
    Rd_h2 = model.h2_disc_Rd();
    Rm_h2 = model.h2_disc_Rm();
    zd_h2 = model.h2_disc_zd();
    Sigma_h2 = model.h2_disc_Sigma_d();
    h2_disc = model.gas_disc_density_mcmillan(R, z, Rd_h2, Rm_h2, zd_h2, Sigma_h2);
    res = h2_disc;
  }
  else
  {
    std::cout << "Model not implemented in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return res;
}

double Galaxy::der_R_H2_density(double R, double z)
{
  double Rd_h2, Rm_h2;
  double res = 0;

  if (m_model_baryons == 0) // McMillan 2017
  {
    Rd_h2 = model.h2_disc_Rd();
    Rm_h2 = model.h2_disc_Rm();
    res = H2_density(R, z) * (-1 / Rd_h2 + Rm_h2 / (R * R));
  }
  else
  {
    std::cout << "Model not implemented in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return res;
}

double Galaxy::der_z_H2_density(double R, double z)
{
  double zd_h2;
  double res = 0;

  if (m_model_baryons == 0) // McMillan 2017
  {
    zd_h2 = model.h2_disc_zd();
    res = -H2_density(R, z) * tanh(z / 2 / zd_h2) / zd_h2;
  }
  else
  {
    std::cout << "Model not implemented in " << __PRETTY_FUNCTION__ << std::endl;
  }
  return res;
}



double Galaxy::disc_surface_density(double R)
{
  double Rd_thick, Sigmad_thick, thick_disc;
  double Rd_thin, Sigmad_thin, thin_disc;
  double Rd_h1, Rm_h1, Sigma_h1, h1_disc;
  double Rd_h2, Rm_h2, Sigma_h2, h2_disc;
  double res = 0;
  if (m_model_baryons == 0) // McMillan 2017
  {
    Rd_thick = model.thick_disc_Rd();
    Sigmad_thick = model.thick_disc_Sigma_d();
    thick_disc = model.exponential_disc_surface_density(R, Rd_thick, Sigmad_thick);
    Rd_thin = model.thin_disc_Rd();
    Sigmad_thin = model.thin_disc_Sigma_d();
    thin_disc = model.exponential_disc_surface_density(R, Rd_thin, Sigmad_thin);
    Rd_h1 = model.h1_disc_Rd();
    Rm_h1 = model.h1_disc_Rm();
    Sigma_h1 = model.h1_disc_Sigma_d();
    h1_disc = model.gas_disc_surface_density_mcmillan(R, Rd_h1, Rm_h1, Sigma_h1);
    Rd_h2 = model.h2_disc_Rd();
    Rm_h2 = model.h2_disc_Rm();
    Sigma_h2 = model.h2_disc_Sigma_d();
    h2_disc = model.gas_disc_surface_density_mcmillan(R, Rd_h2, Rm_h2, Sigma_h2);
    res = thick_disc + thin_disc + h1_disc + h2_disc;
  }
  return res;
}

double Galaxy::integrand_spherical_density(vector<double> z_r)
{
  double z = z_r[0], r = z_r[1];
  double R = (r > fabs(z)) ? sqrt(r * r - z * z) : 0;
  return total_baryon_density(R, z);
}

static double CallBack_integrand_spherical_density(void *pt2object, vector<double> z_r)
{
  return ((Galaxy *)pt2object)->integrand_spherical_density(z_r);
}

void Galaxy::compute_spherical_density()
{
  static int _baryons = -1;
  ostringstream smodel;
  string file;
  double res, rmin = 1e-4, rmax = 500;
  int ivar = 0, isize = 2, ilog = 1, N = 1000;
  double dr = log(rmax / rmin) / (double)N, err;
  vector<double> z_r(2);
  ofstream fout;
  if (_baryons != m_model_baryons)
  {
    smodel << m_model_baryons;
    file = "../input/Galaxy_models/Files/baryon_spherical_density_model" + smodel.str() + ".dat";
    fout.open(file.c_str());
    for (int i = 0; i <= N; ++i)
    {
      z_r[1] = rmin * exp(i * dr);
      res = Simpson_Integral1_Static(ivar, isize, ilog, N, 0.1 * rmin, z_r[1], z_r, this,
                                     CallBack_integrand_spherical_density, err);
      fout << z_r[1] << "   " << res / z_r[1] << endl;
    }
    fout.close();
    _baryons = m_model_baryons;
  }
  return;
}

double Galaxy::interpolate_spherical_density(double r)
{
  static my_spline::spline ss;
  static int _baryons = -1;
  ostringstream smodel;
  string file, line;
  double rr, dens, rmin = 1e-4, rmax = 500;
  vector<double> RR, DENS;
  if (_baryons != m_model_baryons)
  {
    compute_spherical_density();
    smodel << m_model_baryons;
    file = "../input/Galaxy_models/Files/baryon_spherical_density_model" + smodel.str() + ".dat";
    ifstream inputFile(file.c_str());
    while (getline(inputFile, line))
    {
      if (!line.length() || line[0] == '#')
        continue;
      sscanf(line.c_str(), "%lf %lf", &rr, &dens);
      RR.push_back(log(rr)), DENS.push_back(log(dens));
    }
    ss.set_points(RR, DENS);
    RR.clear(), DENS.clear();
    _baryons = m_model_baryons;
  }
  double res = exp(ss(log(r)));
  if (r > rmax)
  {
    res = exp(ss(log(rmax)));
  }
  else if (r < rmin)
  {
    res = exp(ss(log(rmin)));
  }
  return res;
}

double Galaxy::density(double r)
{
  double res = 0;
  if (m_model_dm >= 0)
  {
    res += dm_density(r);
  }
  if (m_model_baryons >= 0)
  {
    res += interpolate_spherical_density(r);
  }
  return res;
}

double Galaxy::integrand_baryon_mass(vector<double> rr)
{
  double r = rr[0];
  return r * r * interpolate_spherical_density(r);
}

static double CallBack_integrand_baryon_mass(void *pt2object, vector<double> rr)
{
  return ((Galaxy *)pt2object)->integrand_baryon_mass(rr);
}

void Galaxy::compute_baryon_mass()
{
  static int _baryons = -1;
  ostringstream smodel;
  string file;
  double res = 0, rmin = 1e-4, rmax = 500, r, err;
  int N = 1000, ivar = 0, isize = 1, ilog = 1;
  double dr = log(rmax / rmin) / (double)N;
  vector<double> rr(1);
  ofstream fout;
  if (_baryons != m_model_baryons)
  {
    smodel << m_model_baryons;
    file = "../input/Galaxy_models/Files/baryon_mass_model" + smodel.str() + ".dat";
    fout.open(file.c_str());
    for (int i = 0; i <= N; ++i)
    {
      r = rmin * exp(i * dr);
      res = Simpson_Integral1_Static(ivar, isize, ilog, N, 0.1 * rmin, r, rr, this,
                                     CallBack_integrand_baryon_mass, err);
      fout << r << "   " << res << endl; // the factor 4*pi is put in the interpolate function
    }
    fout.close();
    _baryons = m_model_baryons;
  }
  return;
}

double Galaxy::interpolate_baryon_mass(double r)
{
  static my_spline::spline ss;
  static int _baryons = -1;
  ostringstream smodel;
  string file, line;
  double rr, dens, rmin = 1e-4, rmax = 500;
  vector<double> RR, DENS;
  if (_baryons != m_model_baryons)
  {
    compute_baryon_mass();
    smodel << m_model_baryons;
    file = "../input/Galaxy_models/Files/baryon_mass_model" + smodel.str() + ".dat";
    ifstream inputFile(file.c_str());
    while (getline(inputFile, line))
    {
      if (!line.length() || line[0] == '#')
        continue;
      sscanf(line.c_str(), "%lf %lf", &rr, &dens);
      RR.push_back(log(rr)), DENS.push_back(log(dens));
    }
    ss.set_points(RR, DENS);
    RR.clear(), DENS.clear();
    _baryons = m_model_baryons;
  }
  double res = 4 * pi * exp(ss(log(r)));
  if (r > rmax)
  {
    res = 4 * pi * exp(ss(log(rmax)));
  }
  else if (r < rmin)
  {
    res = 4 * pi * exp(ss(log(rmin)));
  }
  return res;
}

double Galaxy::mass(double r)
{
  double res = 0;
  if (m_model_dm >= 0)
  {
    res += dm_mass(r);
  }
  if (m_model_baryons >= 0)
  {
    res += interpolate_baryon_mass(r);
  }
  return res;
}

double Galaxy::integrand_baryon_potential(vector<double> rr)
{
  double r = rr[0];
  return interpolate_baryon_mass(r) / (r * r);
}

static double CallBack_integrand_baryon_potential(void *pt2object, vector<double> rr)
{
  return ((Galaxy *)pt2object)->integrand_baryon_potential(rr);
}

void Galaxy::compute_baryon_potential()
{
  static int _baryons = -1;
  ostringstream smodel;
  string file;
  double res = 0, rmin = 1e-4, rmax = 500, r, err;
  int N = 1000, ivar = 0, isize = 1, ilog = 1;
  double dr = log(rmax / rmin) / (double)N;
  vector<double> rr(1);
  ofstream fout;
  if (_baryons != m_model_baryons)
  {
    smodel << m_model_baryons;
    file = "../input/Galaxy_models/Files/baryon_potential_model" + smodel.str() + ".dat";
    fout.open(file.c_str());
    for (int i = 0; i <= N; ++i)
    {
      r = rmin * exp(i * dr);
      res = Simpson_Integral1_Static(ivar, isize, ilog, N, 0.1 * rmin, r, rr, this,
                                     CallBack_integrand_baryon_potential, err);
      fout << r << "   " << res << endl;
    }
    fout.close();
    _baryons = m_model_baryons;
  }
  return;
}

double Galaxy::interpolate_baryon_potential(double r) // [(km/s)^2]
{
  static my_spline::spline ss;
  static int _baryons = -1;
  ostringstream smodel;
  string file, line;
  double rr, dens, rmin = 1e-4, rmax = 500;
  vector<double> RR, DENS;
  if (_baryons != m_model_baryons)
  {
    compute_baryon_potential();
    smodel << m_model_baryons;
    file = "../input/Galaxy_models/Files/baryon_potential_model" + smodel.str() + ".dat";
    ifstream inputFile(file.c_str());
    while (getline(inputFile, line))
    {
      if (!line.length() || line[0] == '#')
        continue;
      sscanf(line.c_str(), "%lf %lf", &rr, &dens);
      RR.push_back(log(rr)), DENS.push_back(log(dens));
    }
    ss.set_points(RR, DENS);
    RR.clear(), DENS.clear();
    _baryons = m_model_baryons;
  }
  double res = G_NEWTON * exp(ss(log(r)));
  if (r > rmax)
  {
    res = G_NEWTON * exp(ss(log(rmax)));
  }
  else if (r < rmin)
  {
    res = G_NEWTON * exp(ss(log(rmin)));
  }
  return res * 1e-6 * Msun_to_kg * m_to_kpc;
}

double Galaxy::gravitational_potential(double r) // [(km/s)^2]
{
  double res = 0;
  if (m_model_dm >= 0)
  {
    res += dm_grav_potential(r);
  }
  if (m_model_baryons >= 0)
  {
    res += interpolate_baryon_potential(r);
  }
  return res;
}

double Galaxy::circular_velocity(double r) // [km/s]
{
  double res = sqrt(G_NEWTON * mass(r) / r);
  res *= 1e-3 * sqrt(m_to_kpc * Msun_to_kg);
  return res;
}


double Galaxy::circular_velocity_baryons_only(double r) // [km/s]
{
  double res = sqrt(G_NEWTON * (mass(r) - dm_mass(r)) / r);
  res *= 1e-3 * sqrt(m_to_kpc * Msun_to_kg);
  return res;
}

void Galaxy::plot_circular_velocity()
{
  std::ofstream outfile;
  outfile.open("../output/RotationCurves.out");

  double R = 0;

  outfile << "# R [kpc] | v (tot) [km/s] | v (bar) [km/s] | v (DM) [km/s]" << std::endl;
  for(int i = 0; i < 201; i++)
  {
    R = 0.5 + i*(20 - 0.5)/200;
    outfile << R << "\t" << circular_velocity(R) << "\t" << circular_velocity_baryons_only(R) << "\t" << circular_velocity_DM_only(R) << std::endl;
  } 

  outfile.close();
}

double Galaxy::circular_period(double r) // [Myr]
{
  double res = 2 * pi * r / circular_velocity(r);
  return res * 1e-3 * kpc_to_m * s_to_Myr;
}

double Galaxy::escape_speed(double r) // [km/s]
{
  return sqrt(2 * fabs(gravitational_potential(r) - gravitational_potential(get_rmax())));
}

double Galaxy::circular_Ncross(double r) 
// r in kpc 
{
  double age_disk = 1e+4; // Myr <- Maybe this should be galaxy (and redshift) dependant
  return 2 * 1e+4 / circular_period(r);
}


// double miyamoto_nagai_potential(double R, double z, double a, double b, double M)
// {
//   double res = G_NEWTON*M/sqrt(R*R+pow(a+sqrt(z*z+b*b),2));
//   return res;
// }
// static double integrand1(vector<double> mu_a_b_M_r)
// {
//   double mu = mu_a_b_M_r[0], a = mu_a_b_M_r[1], b = mu_a_b_M_r[2];
//   double M = mu_a_b_M_r[3], r = mu_a_b_M_r[4];
//   return miyamoto_nagai_potential(r*sqrt(1-mu*mu),r*mu,a,b,M);
// }
// double MN_spherical_potential(double r, double a, double b, double M)
// {
//   vector<double> mu_a_b_M_r(5);
//   mu_a_b_M_r[1] = a, mu_a_b_M_r[2] = b, mu_a_b_M_r[3] = M;
//   mu_a_b_M_r[4] = r;
//   int ivar = 0, isize = 5, ilog = 0, Nsteps = 100;
//   double Rmin = 0, Rmax = 1, err;
//   double res = Simpson_Integral1(ivar,isize,ilog,Nsteps,
// 				 Rmin,Rmax,mu_a_b_M_r,
// 				 &integrand1,err);
//   return res*pow(speed_in_kms,2);
// }
