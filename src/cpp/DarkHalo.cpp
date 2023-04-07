#include "../headers/DarkHalo.h"


using namespace std;


const double XMIN = 1e-6;
const double XMAX = 1e3;


DarkHalo::DarkHalo(Cosmology cosmo)
{
  m_profile = 0;
  m_alpha = 1;
  m_beta = 3;
  m_gamma = 1;
  m_rhos = 1e7;
  m_rs = 20;
  m_cosmo = cosmo;
}

DarkHalo::DarkHalo(Cosmology cosmo, DensityProfile profile, double rhos, double rs)
{
  if(profile == DensityProfile::NFW)
  {
    m_profile = 0;
    m_alpha = 1;
    m_beta = 3;
    m_gamma = 1;
    m_rhos = rhos;
    m_rs = rs;
    m_cosmo = cosmo;
  }

}

DarkHalo::DarkHalo(Cosmology cosmo, int profile, double alpha, double beta, double gamma, double rhos, double rs)
{
  m_profile = profile;
  m_alpha = alpha;
  m_beta = beta;
  m_gamma = gamma;
  m_rhos = rhos;
  m_rs = rs;
  m_cosmo = cosmo;
}

DarkHalo::DarkHalo(Cosmology cosmo, double rhos, double rs)
{
  m_profile = 0;
  m_alpha = 1;
  m_beta = 3;
  m_gamma = 1;
  m_rhos = rhos;
  m_rs = rs;
  m_cosmo = cosmo;
}

void DarkHalo::Initialise_from_virial_parameters(int profile, double alpha, double beta, double gamma, double mvir, double cvir, double delta, bool is_critical)
{
  m_profile = profile;
  m_alpha = alpha;
  m_beta = beta;
  m_gamma = gamma;

  double n_cvir = cvir;

  if (cvir < 0 && delta == 200)
    n_cvir = c_bar(mvir);
  else if (cvir < 0 && delta != 200)
  {
    std::cout << "FATAL ERROR : cvir must be a real value greater or equal to 1 if delta != 200 in " << __PRETTY_FUNCTION__ << std::endl; 
    exit(0);
  }

  m_rhos = scale_density(n_cvir, delta, is_critical);
  m_rs = scale_radius(n_cvir, mvir, delta, is_critical);
}


void DarkHalo::Initialise_from_virial_parameters(double mvir, double cvir, double delta, bool is_critical)
{
   double n_cvir = cvir;

  if (cvir < 0 && delta == 200)
    n_cvir = c_bar(mvir);
  else if (cvir < 0 && delta != 200)
  {
    std::cout << "FATAL ERROR : cvir must be a real value greater or equal to 1 if delta != 200 in " << __PRETTY_FUNCTION__ << std::endl; 
    exit(0);
  }

  m_rhos = scale_density(n_cvir, delta, is_critical);
  m_rs = scale_radius(n_cvir, mvir, delta, is_critical);
}


double DarkHalo::c_bar(double m200)
{
  double cbar = 0;

  double cn[6] = {37.5153, -1.5093, 1.636e-2, 3.66e-4, -2.89237e-5, 5.32e-7};

  double m200_bis = m200;

  if (m200 <= 7.24e-10)
    m200_bis = 7.24e-10;

  for (int i = 0; i < 6; i++)
    cbar += cn[i] * pow(log(m200_bis * m_cosmo.get_h()), i);

  return cbar;
}

void DarkHalo::update_rhos(double rhos){m_rhos = rhos;}
void DarkHalo::update_rs(double rs){m_rs = rs;}
void DarkHalo::update_profile(int profile){m_profile = profile;}
void DarkHalo::update_alpha(double alpha){m_alpha = alpha;}
void DarkHalo::update_beta(double beta){m_beta = beta;}
void DarkHalo::update_gamma(double gamma){m_gamma = gamma;}

int DarkHalo::get_profile() const {return m_profile;}
double DarkHalo::get_alpha() const {return m_alpha;}
double DarkHalo::get_beta() const {return m_beta;}
double DarkHalo::get_gamma() const {return m_gamma;}
double DarkHalo::get_rhos() const {return m_rhos;}
double DarkHalo::get_rs() const {return m_rs;}


double DarkHalo::virial_radius(double delta, bool is_critical)
{
  double rho_ref = (is_critical) ? m_cosmo.critical_density(0) : m_cosmo.cosmological_density(0,Species::MATTER);
  rho_ref *= delta;
  double rmin = 0.1*m_rs, rmax = 1e3*m_rs;
  double rint = sqrt(rmin*rmax);
  double rhomax = dm_mass(rmin)/(4./3*pi)*pow(rmin,-3);
  double rhomin = dm_mass(rmax)/(4./3*pi)*pow(rmax,-3);
  double rhoint = dm_mass(rint)/(4./3*pi)*pow(rint,-3);
  int i = 0;
  if(rho_ref>=rhomax){rint = rmin;}
  else if(rho_ref<=rhomin){rint = rmax;}
  else
    {
      while(fabs(rho_ref-rhoint)/rho_ref>1e-8 && i<100)
	{
	  (rhoint>rho_ref) ? rmin = rint : rmax = rint;
	  rint  = sqrt(rmin*rmax);
	  rhoint = dm_mass(rint)/(4./3*pi)*pow(rint,-3);
	  i += 1;
	}
    }
  return rint;
}

double DarkHalo::virial_mass(double delta, bool is_critical)
{
  double r_vir = virial_radius(delta,is_critical);
  double rho_ref = (is_critical) ? m_cosmo.critical_density(0) : m_cosmo.cosmological_density(0,Species::MATTER);
  rho_ref *= delta;
  double m_vir = 4./3*pi*pow(r_vir,3)*rho_ref;
  return m_vir;
}

double DarkHalo::virial_concentration(double delta, bool is_critical)
{
  double r_2 = m_rs;
  if(m_profile==0){r_2 *= pow((2-m_gamma)/(m_beta-2),1./m_alpha);}
  double r_vir = virial_radius(delta,is_critical);
  return r_vir/r_2;
}

double DarkHalo::scale_density(double c_vir, double delta, bool is_critical)
{
  double rho_ref = (is_critical) ? m_cosmo.critical_density(0) : m_cosmo.cosmological_density(0,Species::MATTER);
  rho_ref *= delta;
  double r2_rs = 1;
  if(m_profile==0){r2_rs *= pow((2-m_gamma)/(m_beta-2),1./m_alpha);}
  double c_eff = r2_rs*c_vir; 
  return 1./3*pow(c_eff,3)*rho_ref/mass_profile(c_eff);
}

double DarkHalo::scale_radius(double c_vir, double m_vir, double delta, bool is_critical) // [kpc]
{
  double rho_ref = (is_critical) ? m_cosmo.critical_density(0) : m_cosmo.cosmological_density(0,Species::MATTER);
  rho_ref *= delta;
  double r2_rs = 1;
  if(m_profile==0){r2_rs *= pow((2-m_gamma)/(m_beta-2),1./m_alpha);}
  double c_eff = r2_rs*c_vir; 
  return pow(3./(4*pi)*m_vir/rho_ref,1./3)/c_eff;
}

double DarkHalo::scale_radius_bis(double rhos, double m_vir, double delta, bool is_critical)
{
  double rho_ref = (is_critical) ? m_cosmo.critical_density(0) : m_cosmo.cosmological_density(0,Species::MATTER);
  rho_ref *= delta;
  double ref = rho_ref/(3*rhos);
  double xmin = 0.1, xmax = 1e4, xmed = sqrt(xmin*xmax);
  double my_max = mass_profile(xmin)*pow(xmin,-3);
  double my_min = mass_profile(xmax)*pow(xmax,-3);
  double my_med = mass_profile(xmed)*pow(xmed,-3);
  int i = 0;
  if(ref>=my_max){xmed = xmin;}
  else if(ref<=my_min){xmed = xmax;}
  else
    {
      while(fabs(ref-my_med)/ref>1e-6 && i<100)
	{
	  (my_med>ref) ? xmin = xmed : xmax = xmed;
	  xmed = sqrt(xmin*xmax);
	  my_med = mass_profile(xmed)*pow(xmed,-3);
	  i += 1;
	}
    }
  double r_vir = pow(3*m_vir/(4*pi*rho_ref),1./3);
  return r_vir/xmed;
}


double DarkHalo::density_profile(double x)
{
  double res = 0;//, x_sat = m_rsat/m_rs;
  switch(m_profile)
    {
    case 0: // alpha-beta-gamma
      if(m_gamma==0)
	{res = pow(1+pow(x,m_alpha),-m_beta/m_alpha);}
      else
	{res = pow(x,-m_gamma)*pow(1+pow(x,m_alpha),(m_gamma-m_beta)/m_alpha);}
      break;
    case 1: // einasto
      res = exp(-2./m_alpha*(pow(x,m_alpha)-1)); break;
    case 2: // burkert
      res = 1./(1+x)*1./(1+pow(x,2)); break;
    }
  return res;
}

double DarkHalo::dm_density(double r)
{return m_rhos*density_profile(r/m_rs);}

double DarkHalo::integrand_mass_profile(vector<double> xx)
{double x = xx[0];
  return x*x*density_profile(x);}

static double CallBack_integrand_mass_profile(void *pt2object, vector<double> xx)
{return ((DarkHalo *)pt2object)->integrand_mass_profile(xx);}

void DarkHalo::compute_mass_profile()
{
  static int _profile = -1;
  static double _alpha = -1, _beta = -1, _gamma = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string file;
  double res = 0, x;
  int N = 1000, ivar = 0, isize = 1, ilog = 1;
  double dx = log(XMAX/XMIN)/(double)N, err;
  vector<double> xx(1);
  ofstream fout;
  if(_profile!=m_profile || _alpha!=m_alpha || _beta!=m_beta || _gamma!=m_gamma)
    {
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/mass_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      fout.open(file.c_str());
      for(int i=0;i<=N;++i)
	{
	  x = XMIN*exp(i*dx);
	  res = Simpson_Integral1_Static(ivar,isize,ilog,N,0.1*XMIN,x,xx,this,
					 CallBack_integrand_mass_profile,err);
	  fout << x << "   " << res << endl;
	}
      fout.close();
      _profile = m_profile, _alpha = m_alpha;
      _beta = m_beta, _gamma = m_gamma;
    }
  return;
}

double DarkHalo::interpolate_mass_profile(double x)
{
  static my_spline::spline ss;
  static int _profile2 = -1;
  static double _alpha2 = -1, _beta2 = -1, _gamma2 = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string line,file;
  double xx, mass;
  vector<double> XX,MASS;
  if(_profile2!=m_profile || _alpha2!=m_alpha || _beta2!=m_beta || _gamma2!=m_gamma)
    {
      compute_mass_profile();
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/mass_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      ifstream inputFile(file.c_str());
      while(getline(inputFile,line))
	{
	  if (!line.length() || line[0] == '#')
	    continue;
	  sscanf(line.c_str(),"%lf %lf",&xx, &mass);
	  XX.push_back(log(xx)), MASS.push_back(log(mass));
	}
      ss.set_points(XX,MASS);
      XX.clear(), MASS.clear();
      _profile2 = m_profile, _alpha2 = m_alpha;
      _beta2 = m_beta, _gamma2 = m_gamma;      
    }
  double res = 0;
  if(x>XMAX){res = exp(ss(log(XMAX)));}
  else if(XMIN<=x && x<=XMAX){res = exp(ss(log(x)));}
  return res;
}

double DarkHalo::mass_profile(double x)
{
  double res = 0;
  switch(m_profile)
    {
    default:
      if(m_alpha==1 && m_beta==3 && (m_gamma==1 || m_gamma==0)) // generalized NFW
	{
	  if(m_gamma==1) // NFW
	    {res = (x>1.e-4) ? log(1+x)-x/(1+x) : 
	      x*x*(0.5-2./3.*x+0.75*x*x-0.8*pow(x,3));}
	  else if(m_gamma==0) // cored NFW
	    {res = (x>1.e-4) ? log(1+x)-0.5*x*(3*x+2)*pow(1+x,-2)
		: pow(x,3)/3.-0.75*pow(x,4)+6./5.*pow(x,5);}
	}
      else if(m_beta==2 && m_gamma==2){res = x;} // SIS
      else if(m_alpha==2 && m_gamma==0 && (m_beta==2 || m_beta==5))
	{
	  if(m_beta==2) // cored IS
	    {res = (x>1.e-4) ? x-atan(x) : pow(x,3)/3.-pow(x,5)/5.;}
	  else if(m_beta==5) // Plummer
	    {res = (x>1.e-4) ? pow(x,3)/3*pow(1+x*x,-1.5)
		: pow(x,3)/3-0.5*pow(x,5);}
	}
      else
	{res = interpolate_mass_profile(x);}
      break;
    case 1: // einasto
      res = interpolate_mass_profile(x);
      break;
    case 2: // burkert
      res = (x>1.e-4) ? -0.5*atan(x)+0.5*log(1+x)+0.25*log(1+x*x)
      : pow(x,3)/3.-pow(x,4)/4.;
      break;
    }
  return res;
}

double DarkHalo::dm_mass(double r)
{return 4*pi*m_rhos*pow(m_rs,3)*mass_profile(r/m_rs);}

double DarkHalo::integrand_grav_potential(vector<double> xx)
{
  double x = xx[0];
  double res = mass_profile(x)/(x*x);
  return res;
}

static double CallBack_integrand_grav_potential(void *pt2object, vector<double> xx)
{return ((DarkHalo *)pt2object)->integrand_grav_potential(xx);}

void DarkHalo::compute_grav_potential_profile()
{
  static int _profile = -1;
  static double _alpha = -1, _beta = -1, _gamma = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string file;
  double res = 0, x;
  int N = 1000, ivar = 0, isize = 1, ilog = 1;
  double dx = log(XMAX/XMIN)/(double)N, err;
  vector<double> xx(1);
  ofstream fout;
  if(_profile!=m_profile || _alpha!=m_alpha || _beta!=m_beta || _gamma!=m_gamma)
    {
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/grav_potential_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      fout.open(file.c_str());
      for(int i=0;i<=N;++i)
	{
	  x = XMIN*exp(i*dx);
	  res = Simpson_Integral1_Static(ivar,isize,ilog,N,0.1*XMIN,x,xx,this,
					 CallBack_integrand_grav_potential,err);
	  fout << x << "   " << res << endl;
	}
      fout.close();
      _profile = m_profile, _alpha = m_alpha;
      _beta = m_beta, _gamma = m_gamma;
    }
  return;
}

double DarkHalo::interpolate_grav_potential_profile(double x)
{
  static my_spline::spline ss;
  static int _profile2 = -1;
  static double _alpha2 = -1, _beta2 = -1, _gamma2 = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string line,file;
  double xx, pot;
  vector<double> XX,POT;
  if(_profile2!=m_profile || _alpha2!=m_alpha || _beta2!=m_beta || _gamma2!=m_gamma)
    {
      compute_grav_potential_profile();
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/grav_potential_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      ifstream inputFile(file.c_str());
      while(getline(inputFile,line))
	{
	  if (!line.length() || line[0] == '#')
	    continue;
	  sscanf(line.c_str(),"%lf %lf",&xx, &pot);
	  XX.push_back(log(xx)), POT.push_back(log(pot));
	}
      ss.set_points(XX,POT);
      XX.clear(), POT.clear();
      _profile2 = m_profile, _alpha2 = m_alpha;
      _beta2 = m_beta, _gamma2 = m_gamma;      
    }
  double res = 0;
  if(x>XMAX){res = exp(ss(log(XMAX)));}
  else if(XMIN<=x && x<=XMAX){res = exp(ss(log(x)));}
  return res;
}

double DarkHalo::grav_potential_profile(double x)
{
  double res = 0;
  switch(m_profile)
    {
    default:
      if(m_alpha==1 && m_beta==3 && (m_gamma==1 || m_gamma==0)) // generalized NFW
	{
	  if(m_gamma==1) // NFW
	    {res = (x>1.e-4) ? 1-log(1+x)/x 
		: x*(0.5-x/3.+0.25*x*x-pow(x,3)/5.+pow(x,4)/6.);}
	  else if(m_gamma==0) // cored NFW
	    {res = (x>1.e-4) ? 0.5*(2+x)/(1+x)-log(1+x)/x
		: x*x/6.-0.25*pow(x,3)+3./10*pow(x,4)-pow(x,5)/3.;}
	}
      else if(m_beta==2 && m_gamma==2){res = log(x);} // SIS
      else if(m_alpha==2 && m_gamma==0 && (m_beta==2 || m_beta==5))
	{
	  if(m_beta==2) // cored IS
	    {res = (atan(x)-x)/x+0.5*log(1+x*x);}
	  else if(m_beta==5) // Plummer
	    {res = (x>1.e-4) ? 1./3.*(1-1/sqrt(1+x*x)) : x*x/6-pow(x,4)/8;}
	}
      else
	{res = interpolate_grav_potential_profile(x);}
      break;
    case 1: // einasto
      res = interpolate_grav_potential_profile(x);
      break;
    case 2: // burkert
      res = (x>1.e-4) ? 0.25/x*(2*(1+x)*atan(x)-2*(1+x)*log(1+x)+(x-1)*log(1+x*x))
	: pow(x,2)/6.-pow(x,3)/12.;
      break;
    }
  return res;
}

double DarkHalo::dm_grav_potential(double r) // [(km/s)^2]
{
  double res = 4*pi*G_NEWTON*m_rhos*m_rs*m_rs*grav_potential_profile(r/m_rs);
  return res*1e-6*Msun_to_kg*m_to_kpc;
}

double DarkHalo::integrand_potential_energy(vector<double> xx)
{
  double x = xx[0];
  double grav_pot = grav_potential_profile(x);
  double rho = density_profile(x);
  return x*x*rho*grav_pot;
}

static double CallBack_integrand_potential_energy(void *pt2object, vector<double> xx)
{return ((DarkHalo *)pt2object)->integrand_potential_energy(xx);}

void DarkHalo::compute_potential_energy_profile()
{
  static int _profile = -1;
  static double _alpha = -1, _beta = -1, _gamma = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string file;
  double res = 0, x;
  int N = 1000, ivar = 0, isize = 1, ilog = 1;
  double dx = log(XMAX/XMIN)/(double)N, err;
  vector<double> xx(1);
  ofstream fout;
  if(_profile!=m_profile || _alpha!=m_alpha || _beta!=m_beta || _gamma!=m_gamma)
    {
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/potential_energy_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      fout.open(file.c_str());
      for(int i=0;i<=N;++i)
	{
	  x = XMIN*exp(i*dx);
	  res = Simpson_Integral1_Static(ivar,isize,ilog,N,0.1*XMIN,x,xx,this,
					 CallBack_integrand_potential_energy,err);
	  fout << x << "   " << res << endl;
	}
      fout.close();
      _profile = m_profile, _alpha = m_alpha;
      _beta = m_beta, _gamma = m_gamma;
    }
  return;
}

double DarkHalo::interpolate_potential_energy_profile(double x)
{
  static my_spline::spline ss;
  static int _profile2 = -1;
  static double _alpha2 = -1, _beta2 = -1, _gamma2 = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string line,file;
  double xx, pot;
  vector<double> XX,POT;
  if(_profile2!=m_profile || _alpha2!=m_alpha || _beta2!=m_beta || _gamma2!=m_gamma)
    {
      compute_potential_energy_profile();
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../SRC/Galaxy_models/Files/potential_energy_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      ifstream inputFile(file.c_str());
      while(getline(inputFile,line))
	{
	  if (!line.length() || line[0] == '#')
	    continue;
	  sscanf(line.c_str(),"%lf %lf",&xx, &pot);
	  XX.push_back(log(xx)), POT.push_back(log(pot));
	}
      ss.set_points(XX,POT);
      XX.clear(), POT.clear();
      _profile2 = m_profile, _alpha2 = m_alpha;
      _beta2 = m_beta, _gamma2 = m_gamma;      
    }
  double res = 0;
  if(x>XMAX){res = exp(ss(log(XMAX)));}
  else if(XMIN<=x && x<=XMAX){res = exp(ss(log(x)));}
  return res;
}

double DarkHalo::potential_energy_profile(double x) // negative
{
  double res = 0;
  switch(m_profile)
    {
    default:
      if(m_alpha==1 && m_beta==3 && m_gamma==1) // NFW
	{res = (x>1e-4) ? 0.5*((2+x)*log(1+x)-2*x)/(1+x) : pow(x,3)/12.-pow(x,4)/6.+29./120*pow(x,5);}
      else
	{res = interpolate_potential_energy_profile(x);}
      break;
    case 1: // einasto
      res = interpolate_potential_energy_profile(x);
      break;
    case 2: // burkert
      res = interpolate_potential_energy_profile(x);
      break;
    }
  return res-0.5*grav_potential_profile(x)*mass_profile(x);
}

double DarkHalo::dm_potential_energy(double r) // [Msun*(km/s)^2], negative
{
  double res = pow(4*pi*m_rhos,2)*G_NEWTON*pow(m_rs,5)*potential_energy_profile(r/m_rs);
  return res*1e-6*m_to_kpc*Msun_to_kg;
}

double DarkHalo::integrand_binding_energy(vector<double> xx)
{
  double x = xx[0];
  double rho = density_profile(x);
  double mass = mass_profile(x);
  return x*rho*mass;
}

static double CallBack_integrand_binding_energy(void *pt2object, std::vector<double> xx)
{return ((DarkHalo *)pt2object)->integrand_binding_energy(xx);}

void DarkHalo::compute_binding_energy_profile()
{
  static int _profile = -1;
  static double _alpha = -1, _beta = -1, _gamma = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string file;
  double res = 0, x;
  int N = 1000, ivar = 0, isize = 1, ilog = 1;
  double dx = log(XMAX/XMIN)/(double)N, err;
  vector<double> xx(1);
  ofstream fout;
  if(_profile!=m_profile || _alpha!=m_alpha || _beta!=m_beta || _gamma!=m_gamma)
    {
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/binding_energy_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      fout.open(file.c_str());
      for(int i=0;i<=N;++i)
	{
	  x = XMIN*exp(i*dx);
	  res = Simpson_Integral1_Static(ivar,isize,ilog,N,0.1*XMIN,x,xx,this,
					 CallBack_integrand_binding_energy,err);
	  fout << x << "   " << res << endl;
	}
      fout.close();
      _profile = m_profile, _alpha = m_alpha;
      _beta = m_beta, _gamma = m_gamma;
    }
  return;
}

double DarkHalo::interpolate_binding_energy_profile(double x)
{
  static my_spline::spline ss;
  static int _profile2 = -1;
  static double _alpha2 = -1, _beta2 = -1, _gamma2 = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string line,file;
  double xx, pot;
  vector<double> XX,POT;
  if(_profile2!=m_profile || _alpha2!=m_alpha || _beta2!=m_beta || _gamma2!=m_gamma)
    {
      compute_binding_energy_profile();
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/binding_energy_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      ifstream inputFile(file.c_str());
      while(getline(inputFile,line))
	{
	  if (!line.length() || line[0] == '#')
	    continue;
	  sscanf(line.c_str(),"%lf %lf",&xx, &pot);
	  XX.push_back(log(xx)), POT.push_back(log(pot));
	}
      ss.set_points(XX,POT);
      XX.clear(), POT.clear();
      _profile2 = m_profile, _alpha2 = m_alpha;
      _beta2 = m_beta, _gamma2 = m_gamma;      
    }
  double res = 0;
  if(x>XMAX){res = exp(ss(log(XMAX)));}
  else if(XMIN<=x && x<=XMAX){res = exp(ss(log(x)));}
  return res;
}

double DarkHalo::binding_energy_profile(double x)
{
  double res = 0;
  switch(m_profile)
    {
    default:
      if(m_alpha==1 && m_beta==3 && m_gamma==1) // NFW
	{res = (x>1e-4) ? 0.5*(2*x+x*x-2*(1+x)*log(1+x))*pow(1+x,-2) : pow(x,3)/6.-5./12*pow(x,4)+43./60*pow(x,5);}
      else
	{res = interpolate_binding_energy_profile(x);}
      break;
    case 1: // einasto
      res = interpolate_binding_energy_profile(x);
      break;
    case 2: // burkert
      res = interpolate_binding_energy_profile(x);
      break;
    }
  return res;
}

double DarkHalo::dm_binding_energy(double r) // [Msun*(km/s)^2]
{
  double res = pow(4*pi*m_rhos,2)*G_NEWTON*pow(m_rs,5)*binding_energy_profile(r/m_rs);
  return res*1e-6*m_to_kpc*Msun_to_kg;
}

double DarkHalo::integrand_luminosity_profile(vector<double> xx)
{
  double x = xx[0];
  return x*x*density_profile(x)*density_profile(x);
}

static double CallBack_integrand_luminosity_profile(void *pt2object, vector<double> xx)
{return ((DarkHalo *)pt2object)->integrand_luminosity_profile(xx);}

void DarkHalo::compute_luminosity_profile()
{
  static int _profile = -1;
  static double _alpha = -1, _beta = -1, _gamma = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string file;
  double res = 0, x;
  int N = 1000, ivar = 0, isize = 1, ilog = 1;
  double dx = log(XMAX/XMIN)/(double)N, err;
  vector<double> xx(1);
  ofstream fout;
  if(_profile!=m_profile || _alpha!=m_alpha || _beta!=m_beta || _gamma!=m_gamma)
    {
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/luminosity_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      fout.open(file.c_str());
      for(int i=0;i<=N;++i)
	{
	  x = XMIN*exp(i*dx);
	  res = Simpson_Integral1_Static(ivar,isize,ilog,N,0.1*XMIN,x,xx,this,
					 CallBack_integrand_luminosity_profile,err);
	  fout << x << "   " << res << endl;
	}
      fout.close();
      _profile = m_profile, _alpha = m_alpha;
      _beta = m_beta, _gamma = m_gamma;
    }
  return;
}

double DarkHalo::interpolate_luminosity_profile(double x)
{
  static my_spline::spline ss;
  static int _profile2 = -1;
  static double _alpha2 = -1, _beta2 = -1, _gamma2 = -1;
  ostringstream sprof,salpha,sbeta,sgamma;
  string line,file;
  double xx, lum;
  vector<double> XX,LUM;
  if(_profile2!=m_profile || _alpha2!=m_alpha || _beta2!=m_beta || _gamma2!=m_gamma)
    {
      compute_luminosity_profile();
      sprof << m_profile, salpha << m_alpha, sbeta << m_beta, sgamma << m_gamma;
      file = "../input/Galaxy_models/Files/luminosity_profile"+sprof.str()+"_alpha"+salpha.str()+"_beta"+sbeta.str()+"_gamma"+sgamma.str()+".dat";
      ifstream inputFile(file.c_str());
      while(getline(inputFile,line))
	{
	  if (!line.length() || line[0] == '#')
	    continue;
	  sscanf(line.c_str(),"%lf %lf",&xx, &lum);
	  XX.push_back(log(xx)), LUM.push_back(log(lum));
	}
      ss.set_points(XX,LUM);
      XX.clear(), LUM.clear();
      _profile2 = m_profile, _alpha2 = m_alpha;
      _beta2 = m_beta, _gamma2 = m_gamma;      
    }
  double res = 0;
  if(x>XMAX){res = exp(ss(log(XMAX)));}
  else if(XMIN<=x && x<=XMAX){res = exp(ss(log(x)));}
  return res;
}

double DarkHalo::luminosity_profile(double x)
{
  double res = 0;
  switch(m_profile)
    {
    default:
      if(m_alpha==1 && m_beta==3 && m_gamma==1) // NFW
	{res = (x>1e-3) ? 1./3.*(1-pow(1+x,-3.)) : 
	    x-2*x*x+10./3*pow(x,3.)-5*pow(x,4.)+7*pow(x,5.);}
      else if(m_alpha==2 && m_gamma==0 && m_beta==2) // cored SIS
	{res = 0.5*(x/(1+x*x)+atan(x));}
      else
	{res = interpolate_luminosity_profile(x);}
      break;
    case 1: // einasto
      res = interpolate_luminosity_profile(x);
      break;
    case 2: // burkert
      res = interpolate_luminosity_profile(x);
      break;
    }
  return res;
}

double DarkHalo::dm_luminosity(double r) // [Msun/kpc^3]
{return 4*pi*pow(m_rhos,2)*pow(m_rs,3)*luminosity_profile(r/m_rs);}



double DarkHalo::rm2_rs()
{
  if (m_profile == 0)
    return pow((2 - m_gamma) / (m_beta - 2), 1 / m_alpha);
  else
  {
    std::cout << "FATAL ERROR : No eta defined in " << __PRETTY_FUNCTION__ << " for this model" << std::endl; 
    exit(0);
  }
}



// Velocity dispersion, depends on the profile and on the mass scale
double DarkHalo::velocity_dispersion_profile(double x, double x0)
// Gives the result dimensionless
{
  if (m_alpha == 1 && m_beta == 3 && m_gamma == 1)
  {
    double first_term = 0;
    double J = 0, res = 0;

    if (abs(x - x0) / abs(x0) > 1e-2)
    {
      first_term = f_for_velocity_disperison(x0) - f_for_velocity_disperison(x);
      J = 3 * (polyLog(2, 1 / (1 + x)) - polyLog(2, 1 / (1 + x0)));
      //std::cout << sqrt(x * pow(1 + x, 2.) * (first_term + J)) << std::endl;
      return sqrt(x * pow(1 + x, 2.) * (first_term + J));
    }
    else
    {
      //std::cout << "here" << std::endl;
      res = sqrt((-1. + 1. / (1 + x) + log(1. + x)) * (x0 - x) / (x * x) + (3 * x * (1 + 2 * x) - (1 + x) * (3 + 5 * x) * log(1 + x)) * (x0 - x) * (x0 - x) / (2 * x * x * x * (1 + x) * (1 + x)));
      //std::cout << res << std::endl;
      return res;
    }

    //std::cout << J << " " << first_term << " " << x * pow(1 + x, 2.) << std::endl;

    //return sqrt(x * pow(1 + x, 2.) * (first_term + J));
    //return first_term;
  }
  if (m_alpha == 2 && m_beta == 5 && m_gamma == 0)
    return sqrt(pow(1 + x * x, 5. / 2.) * (pow(1 + x * x, -3) - pow(1 + x0 * x0, -3)) / 18);

  exit(0);
}

double DarkHalo::velocity_dispersion(double r, double r0)
// Gives the result en km/s
{
  if ((m_profile == 0 && m_alpha == 1 && m_beta == 3 && m_gamma == 1) || (m_profile == 0 && m_alpha == 2 && m_beta == 5 && m_gamma == 0))
  {
    double E0 = 4 * PI * G_NEWTON * m_rhos * pow(m_rs, 2) * 1e-6 * m_to_kpc * MSOL_to_kg;
    return sqrt(E0) * velocity_dispersion_profile(r / m_rs, r0 / m_rs);
  }

  // If it is not a tabulated profile we use the
  return velocity_dispersion_integral(r, r0);
}

double DarkHalo::orbital_frequency_profile(double x, double x0)
{
  return velocity_dispersion_profile(x, x0) / x;
}

double DarkHalo::orbital_frequency(double r, double r0)
// Result in s^{-1}
{
  return velocity_dispersion(r, r0) / (r * kpc_to_m * 1e-3);
}

double DarkHalo::orbital_frequency_integral(double r, double r0)
// Result in s^{-1}
{
  return velocity_dispersion_integral(r, r0) / (r * kpc_to_m * 1e-3);
}

double DarkHalo::f_for_velocity_disperison(double y)
{
  return (1 + 9 * y + 7 * y * y) / (2 * y * (1 + y) * (1 + y)) - 0.5 * log((1 + y) / y) + ((6 * y * y + 3 * y - 1) / (2 * y * y * (1 + y)) - 3 * log((1 + y) / y)) * log(1 + y);
}

double DarkHalo::fToInt_velocity_dispersion_integral(double r)
{
  return dm_mass(r) * dm_density(r) / (r * r);
}

double DarkHalo::velocity_dispersion_integral(double r, double r0)
// result in km/s
{
  std::vector<double> rvec = {0};
  return sqrt(G_NEWTON * MSOL_to_kg * pow(m_to_kpc, 3) / dm_density(r) * GaussLegendre_IntegralLn_Static(0, 1, gl100, r, r0, rvec, this, CallBack_fToInt_velocity_dispersion_integral)) * kpc_to_m * 1e-3;
}


double DarkHalo::circular_velocity_DM_only(double r) // r [kpc] -> [km/s] 
{
  double res = sqrt(G_NEWTON * dm_mass(r) / r);
  res *= 1e-3 * sqrt(m_to_kpc * Msun_to_kg);
  return res;
}