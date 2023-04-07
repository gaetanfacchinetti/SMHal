#include "../headers/Cosmology.h"


Cosmology::Cosmology(CosmoModel n_model) : model(n_model)
{
  if(model == CosmoModel::PLANCKONLY)
  {
    _Omega_c_h2 = 0.1200;
    _Omega_b_h2 = 0.02237;
    _ns = 0.9649;
    _As = 1e-10*exp(3.044);
    _h = 0.6736;
    _sigma8 = 0.8111;
    // From https://arxiv.org/abs/1807.06209 (Table 2 - 5th column)
    // Note that here we consider no neutrino species
  }
  else if (model == CosmoModel::EdSPLANCK)
  {
    _Omega_c_h2 = 0.3*0.6736*0.6736;
    _Omega_b_h2 = 0;
    _ns = 0.9649;
    _As = 1e-10*exp(3.044);
    _h = 0.6736;
    _sigma8 = 0.8111;
  }
  else
  {   
    _Omega_c_h2 = 0.1200;
    _Omega_b_h2 = 0.02237;
    _ns = 0.9649;
    _As = 1e-10*exp(3.044);
    _h = 0.6736;
    _sigma8 = 0.8111;
    // Evaluated to PLANCKONLY by default
   } 

    double H0_GeV = ((1e+4*_h)/(kpc_to_cm))*(H_BAR/GeV_to_J);
    double rho_crit_0_GeV = (3.0/(8.0*PI))*pow(H0_GeV,2)*pow(PLANCK_MASS,2);
    double Omega_gamma = 2*((PI*PI)/30)*pow(T0_CMB_GeV,4)/rho_crit_0_GeV;
    double Omega_neutr = Omega_gamma * N_EFF * (7./8.) * pow(4./11.,4./3.);

    _Omega_m_h2 = _Omega_c_h2 + _Omega_b_h2;
    _Omega_r_h2 = (Omega_gamma + Omega_neutr)*_h*_h;
    _Omega_l_h2 = pow(_h,2) - _Omega_m_h2 - _Omega_r_h2;


    // Initialisation of the equivalence parameters
    _aeq_rm = _Omega_r_h2 / _Omega_m_h2; // Assuming a0 =1
    _zeq_rm = 1./_aeq_rm - 1;
    _keq_Mpcm1 = _aeq_rm*Hubble_parameter(_zeq_rm)/(LIGHT_SPEED*1e-3);

    // Set the interpolation of the inverse growth to false
    is_inverse_growth_Caroll_interpolated = false;
}

Cosmology::Cosmology(double n_Omega_c_h2, double n_Omega_b_h2, double n_ns, double n_As, double n_h) : _Omega_c_h2(n_Omega_c_h2), _Omega_b_h2(n_Omega_b_h2), _ns(n_ns), _As(n_As)
{
    _h = n_h;
    
    double h0_GeV = ((1e+4*_h)/(kpc_to_cm))*(H_BAR/GeV_to_J);
    double rho_crit_0_GeV = (3.0/(8.0*PI))*pow(h0_GeV,2)*pow(PLANCK_MASS,2);
    double Omega_gamma = (2*((PI*PI)/30)*pow(T0_CMB_GeV,4)/rho_crit_0_GeV);
    double Omega_neutr = Omega_gamma * N_EFF * (7./8.) * pow(4./11.,4./3.);

    _Omega_m_h2 = _Omega_c_h2 + _Omega_b_h2;
    _Omega_r_h2 = (Omega_gamma + Omega_neutr)*_h*_h;
    _Omega_l_h2 = pow(_h,2) - _Omega_m_h2 - _Omega_r_h2;

    _sigma8 = -1;
    model = CosmoModel::UNKNOWN; // A priori unknown
  
}

double Cosmology::cosmic_abundance(double z, Species spe)
// Omega_{i} with i=0 : CDM
//                i=1 : pressureless baryons
//                i=2 : radiation (all neutrinos massless)
//                i=3 : Lambda
{

  //std::cout << " aa " << _Omega_l_h2 << std::endl;
  double res = 0;
  double E2_h2 = _Omega_m_h2*pow(1+z,3)+ _Omega_r_h2*pow(1+z,4)+_Omega_l_h2;
  
  if (spe == Species::CDM)
      res = _Omega_c_h2*pow(1+z,3); 
  
  if (spe == Species::BARYONS)
      res = _Omega_b_h2*pow(1+z,3); 
  
  if (spe == Species::RADIATION)
      res = _Omega_r_h2*pow(1+z,4); 
  
  if (spe == Species::LAMBDA)
      res = _Omega_l_h2; 

  if (spe == Species::MATTER)
      res = _Omega_m_h2*pow(1+z,3); 
  
  return res/E2_h2;
}

double Cosmology::Ez(double z)
{
  return sqrt((_Omega_m_h2*pow(1+z,3) + _Omega_r_h2*pow(1+z,4) + _Omega_l_h2)/(_h*_h));
}


double Cosmology::Hubble_parameter(double z) // [km/s/Mpc]
{
  double E2_h2 = _Omega_m_h2*pow(1+z,3) + _Omega_r_h2*pow(1+z,4) + _Omega_l_h2;
  double H_z = 100*sqrt(E2_h2);
  return H_z;
}


double Cosmology::critical_density(double z) // [Msun/kpc^3]
{
  double H2_z = pow(Hubble_parameter(z),2);
  double res = 3*H2_z/(8*PI*G_NEWTON);
  return res*kg_to_Msun*kpc_to_m;
}

double Cosmology::cosmological_density(double z, Species spe) 
// [Msun/kpc^3] with i=0 : CDM
//                   i=1 : pressureless baryons
//                   i=2 : radiation (all neutrinos massless)
//                   i=3 : Lambda
{
  return cosmic_abundance(z, spe)*critical_density(z);  // rho_i
}


double Cosmology::growth_factor_D1_Carroll(double z)
// Carroll+ 1992  // Mo and White p. 172
// corresponds to D1(a=1) with the definition of Dodelson 2003 -> To be checked
{
  //double Omega_c_z = cosmic_abundance(z, Species::CDM);
  //double Omega_b_z = cosmic_abundance(z, Species::BARYONS);
  //double Omega_l_z = cosmic_abundance(z, Species::LAMBDA);
  //double Omega_m_z = cosmic_abundance(z, Species::MATTER);

  // Abundances in a Universe with no radiation
  double E2_h2_bis = _Omega_m_h2*pow(1+z,3) + _Omega_l_h2;
  double Om = _Omega_m_h2*pow(1+z,3)/E2_h2_bis;
  double Ol = _Omega_l_h2/E2_h2_bis;

  double res = 2.5*Om/(pow(Om,4./7)-Ol+(1+0.5*Om)*(1+1./70.*Ol))/(1+z);
  return res;
}


double Cosmology::growth_factor_D1_Belloso(double z)
// Belloso+01
// corresponds to D1(a=1) with the definition of Dodelson 2003 -> To be checked
{
  double Omega_c_z = cosmic_abundance(z, Species::CDM);
  double Omega_b_z = cosmic_abundance(z, Species::BARYONS);
  double Omega_l_z = cosmic_abundance(z, Species::LAMBDA);
  double Omega_m_z = cosmic_abundance(z, Species::MATTER);

  //std::cout << Omega_m_z << std::endl;

  double res = pow(Omega_m_z, 1./3.)*gsl_sf_hyperg_2F1(1./3., 5./6., 1 + 5. / 6., 1 -  Omega_m_z)/(1+z);
  return res;
}


double Cosmology::ftoIntOnz_growth_factor_exact(double z)
{
  double Ez_bis = sqrt(_Omega_m_h2*pow(1+z,3) + _Omega_l_h2)/_h;
  return (1+z)/pow(Ez_bis, 3);
}


double Cosmology::growth_factor_exact(double z)
{
  std::vector<double> xx = {0};
  double E2_h2_bis = _Omega_m_h2*pow(1+z,3) + _Omega_l_h2;
  double norm = (5./2.)*_Omega_m_h2/(_h*_h)*sqrt(E2_h2_bis)/(_h);
  return norm*GaussLegendre_IntegralLn_Static(0, 1, gl100, z, 1e+8, xx, this, CallBack_ftoIntOnz_growth_factor_exact);
}


double Cosmology::derLn_growth_factor_exact(double z)
{
  double Ez_bis = sqrt(_Omega_m_h2*pow(1+z,3) + _Omega_l_h2)/_h;
  double res = 3*_Omega_m_h2/(_h*_h)*pow(1+z, 2)/pow(Ez_bis, 2)/2. - 5./2.*_Omega_m_h2/(_h*_h)*(1+z)/pow(Ez_bis, 2)/growth_factor_exact(z);
  return res;
}


double Cosmology::Der_growth_factor_exact(double z)
{
  return derLn_growth_factor_exact(z)*growth_factor_exact(z);
}


double Cosmology::derLn_growth_factor_D1_Carroll(double z)
// Carroll+ 1992  // Mo and White p. 172
// dlnD/dz -> WARNING it is a negative quantity
{
  // Abundances in a Universe with no radiation
  double E2_h2_bis = _Omega_m_h2*pow(1+z,3) + _Omega_l_h2;
  double Om = _Omega_m_h2*pow(1+z,3)/E2_h2_bis;
  double Ol = _Omega_l_h2/E2_h2_bis;

  double res = -pow((Om*pow(1+z,3))/(Om*pow(1+z,3) - (Om+Ol-1)*pow(1+z,2) + Ol), 4./7.)/(1+z);
  
  return res;
}

double Cosmology::Inverse_growth_factor_Caroll(double D1)
{
  if(D1 > growth_factor_D1_Carroll(0))
  {
    // By default we return -1 ... probably not a goot prescrption
    return -1.;
  }

  // If not interpolated we interpolate first
  if(is_inverse_growth_Caroll_interpolated == false)
    Interpolate_inverse_growth_factor_Caroll();

  return pow(10, spline_inverse_log10z_Caroll_log10D1(log10(D1)));
}

double Cosmology::_Inverse_growth_factor_Caroll(double y)
{
  std::vector<double> xx = {0, y};
  return DichotomieLn_Static(0, 2, xx, this, CallBack_f_forBissection_Inverse_growth_factor_Caroll, 1e-2, 1e+5, 1e-6, 0);
}

double Cosmology::f_forBissection_Inverse_growth_factor_Caroll(double z, double y)
{
  //std::cout << growth_factor_D1_Carroll(z) << " " << growth_factor_D1_Carroll(z) - y << std::endl;
  return growth_factor_D1_Carroll(z) - y;
}

void Cosmology::Interpolate_inverse_growth_factor_Caroll()
{
  int Npts = 500;
  double Dmin = 1.1/(1+1e+5), Dmax = growth_factor_D1_Carroll(0), dD = log10(Dmax/Dmin)/(1.0*Npts), D;
  std::vector<double> log10_D_vec, log10_z_vec;
  log10_D_vec.resize(0);
  log10_z_vec.resize(0);

  for(int i = 0; i<= Npts; i++)
  {
    D = Dmin*pow(10, i*dD);
    log10_D_vec.push_back(log10(D));
    log10_z_vec.push_back(log10(_Inverse_growth_factor_Caroll(D)));
  }

  spline_inverse_log10z_Caroll_log10D1.set_points(log10_D_vec, log10_z_vec, false);

  is_inverse_growth_Caroll_interpolated = true;

}


//
//

double Cosmology::Der_growth_factor_D1_Carroll(double z)
// DOES NOT WORK
{
  double Omm = cosmic_abundance(0, Species::MATTER);
  double Oml = cosmic_abundance(0, Species::LAMBDA);

  double E = sqrt(Omm*pow(1+z,3) + Oml);
  double der_E =( 3*Omm*pow(1+z,2) )/(2*E);
  
  // need to investigate on why a factor h here
  return (der_E/E)*growth_factor_D1_Carroll(z) - (1+z)*pow(E, -2);
}

double Cosmology::Der_growth_factor_D1_Carroll_numerical(double z)
{
  std::vector<double> zz{z};
  return abs(myDerivative_Static(1, 0, zz, this, CallBack_growth_factor_D1_Carroll, 1e-3));
}





void Cosmology::plot_cosmic_abundance(std::string cosmo)
{

  std::ofstream outfile;
  outfile.open("../output/cosmic_abundances_" + cosmo + ".out");
  outfile.precision(8);

  outfile << "# Omega_c_h2 = " << _Omega_c_h2 << std::endl;
  outfile << "# Omega_v_h2 = " << _Omega_b_h2 << std::endl;
  outfile << "# Spectral index : ns = " <<  _ns << std::endl;
  outfile <<  "# Spectral amplitude = " << _As << std::endl; 
  outfile << "# Hubble parameter (reduced) = "<< _h << std::endl;

  double z, zmin = 1e-3, zmax = 1e+7;

  int Npts = 200;

  double d_logz = log10(zmax/zmin)/Npts;

  outfile << "# z | Omega_cdm | Omega_b | Omega_r | Omega_l" << std::endl;

  for(int i = 0; i<Npts+1; i++)
  {
    z = zmin*pow(10, i*d_logz);
    outfile << z << "\t" << cosmic_abundance(z, Species::CDM) << "\t" << cosmic_abundance(z, Species::BARYONS) << "\t" << cosmic_abundance(z, Species::RADIATION) << "\t" << cosmic_abundance(z,Species::LAMBDA) << std::endl;
  }

  outfile.close();
  
}

void Cosmology::plot_Der_growth_factor()
{

  std::ofstream outfile;
  outfile.open("../output/UCMHS/Der_growth_factor.out");
  outfile.precision(8);

  double z, zmin = 1e-2, zmax = 1e+5;

  int Npts = 5000;

  double d_logz = log10(zmax/zmin)/Npts;

  outfile << "# z | Der_growth_factor  | Der_growth_factor_approximate" << std::endl;

  double z0;

  for(int i = 0; i<Npts+1; i++)
  {
    z = zmin*pow(10, i*d_logz);
    outfile << z << "\t" << abs(Der_growth_factor_D1_Carroll(z)) << "\t" << abs((growth_factor_D1_Carroll(z)-growth_factor_D1_Carroll(z0))/abs(z-z0)) 
           << "\t" << Der_growth_factor_D1_Carroll_numerical(z)<< "\t" << abs(Der_growth_factor_exact(z)) << std::endl;
    z0 = z;
  }

  outfile.close();
  
}

void Cosmology::plot_growth_factor()
{

  std::ofstream outfile;
  outfile.open("../output/UCMHs/growth_factor.out");
  outfile.precision(8);

  double z, zmin = 1e-2, zmax = 1e+5;

  int Npts = 5000;

  double d_logz = log10(zmax/zmin)/Npts;

  outfile << "# z | growth_factor" << std::endl;

  double z0;

  for(int i = 0; i<Npts+1; i++)
  {
    z = zmin*pow(10, i*d_logz);
    outfile << z << "\t" << growth_factor_D1_Carroll(z) << "\t" << growth_factor_D1_Belloso(z) << "\t" << growth_factor_exact(z) << std::endl;
    z0 = z;
  }

  outfile.close();
  
}


void Cosmology::plot_Hubble_Parameter(std::string cosmo)
{

  std::ofstream outfile;
  outfile.open("../output/Hubble_Parameter_" + cosmo +".out");
  outfile.precision(8);

  outfile << "# Omega_c_h2 = " << _Omega_c_h2 << std::endl;
  outfile << "# Omega_v_h2 = " << _Omega_b_h2 << std::endl;
  outfile << "# Spectral index : ns = " <<  _ns << std::endl;
  outfile <<  "# Spectral amplitude = " << _As << std::endl; 
  outfile << "# Hubble parameter (reduced) = "<< _h << std::endl;

  double z, zmin = 1e-3, zmax = 1e+7;

  int Npts = 500;

  double d_logz = log10(zmax/zmin)/Npts;

  outfile << "# z | H" << std::endl;

    for(int i = 0; i<Npts+1; i++)
  {
    z = zmin*pow(10, i*d_logz);
    outfile << z << "\t" << Hubble_parameter(z) << "\t" << CosmicTime_in_yr(z) << std::endl;
  }

  outfile.close();

}

double Cosmology::redshift_horizon_entry(double k)
{
  std::vector<double> xx = {0, k};
  return DichotomieLn_Static(0, 2, xx, this, CallBack_fToSolve_redshift_horizon_entry, 1e-2, 1e+10, 1e-6, 0);
}

double Cosmology::fToSolve_redshift_horizon_entry(double z, double k)
{
  return Hubble_parameter(z)/(1+z)/(LIGHT_SPEED*1e-3) - k;
}

//
//

double Cosmology::CosmicTime_in_yr(double z)
// Result in [yr]
{
  std::vector<double> zz {0};
  GaussLegendre gl200(200);
  return GaussLegendre_IntegralLn_Static(0, 1, gl200, z, 1e+6, zz, this, CallBack_fToIntForCosmicTime)*1e-2/cm_to_kpc/365.25/24./3600.;
  //return 0;
}

double Cosmology::CosmicTime_times_Hubble_parameter(double z)
// Result dimensionless
{
  std::vector<double> zz {0};
  GaussLegendre gl200(200);
  return Hubble_parameter(z)*GaussLegendre_IntegralLn_Static(0, 1, gl200, z, 1e+6, zz, this, CallBack_fToIntForCosmicTime);
  //return 0;
}

double Cosmology::fToIntForCosmicTime(double z)
{
  return 1./(1+z)/Hubble_parameter(z); 
}



// Overcharging the operators

bool operator== (const Cosmology& c1, const Cosmology& c2)
{
    return (c1._Omega_c_h2 == c2._Omega_c_h2 && c1._Omega_b_h2 == c2._Omega_b_h2 && c1._ns == c2._ns && c1._As == c2._As && c1._h == c2._h);
}

bool operator!= (const Cosmology& c1, const Cosmology& c2)
{
    return (c1._Omega_c_h2 != c2._Omega_c_h2 || c1._Omega_b_h2 != c2._Omega_b_h2 || c1._ns != c2._ns || c1._As != c2._As || c1._h != c2._h);
}
