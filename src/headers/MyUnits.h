#ifndef DEF_MYUNITS
#define DEF_MYUNITS
///
/// Definitions of unit conversions, fundamental constants, standard masses  
/// and couplings.                                                           
///


#include "CppLib.h"
#include "mymath.h"


/*

  Martin definition

*/

static const double G_NEWTON = 6.67408e-11; // [m^3.kg^-1.s^-2]
static const double LIGHT_SPEED = 299792458; // [m/s]
static const double K_BOLTZMANN = 1.3806488e-23; // Boltzmann constant [J.K^-1]
static const double H_BAR = 1.054571726e-34; // reduced Planck constant [J.s]
static const double E_CHARGE = 1.602176565e-19;

// conversion factors 

static const double Msun_to_kg = 1.32712442099e20/G_NEWTON; 
static const double kg_to_Msun = 1./Msun_to_kg;
static const double kpc_to_m = 3.08567758149e19;
static const double kpc_to_cm = 3.08567758149e+21;
static const double cm_to_kpc = 1./kpc_to_cm;
static const double m_to_kpc = 1./kpc_to_m; 
static const double ev_to_kg = 1.782661845e-36;
static const double kg_to_ev = 1./ev_to_kg;
static const double Mskpc3_to_gevcm3 = 1e-15*Msun_to_kg*kg_to_ev*pow(m_to_kpc,3.);
static const double gevcm3_to_Mskpc3 = 1./Mskpc3_to_gevcm3;
static const double Myr_to_s = 1e6*365.25*24*3600;
static const double s_to_Myr = 1/Myr_to_s;
static const double m_eV = 1e9/(197.3269718); // 1 meter times 1 eV
static const double eV_to_J = E_CHARGE; // Proton charge in [C]
static const double GeV_to_J = 1.e9*eV_to_J;
static const double J_to_eV = 1./eV_to_J;
static const double eV_to_K = eV_to_J/K_BOLTZMANN;
static const double GeV_to_K = 1.e9*eV_to_K;
static const double K_to_eV = 1./eV_to_K;
static const double K_to_GeV = 1./GeV_to_K;
static const double s_ev = m_eV*LIGHT_SPEED; // 1 second times 1 eV
static const double speed_in_kms = 1e-3*sqrt(m_to_kpc*Msun_to_kg); // conversion factor: (G*Msun/kpc)^1/2 ---> km/s

static const double PLANCK_MASS = sqrt(H_BAR*LIGHT_SPEED/G_NEWTON)*kg_to_ev*1e-9; // Planck mass [GeV]

// cosmological parameters (Planck 2018) - Everything except BAO


static const double T0_CMB_K = 2.72548; // From Fixsen arXiv:0911.1955
static const double T0_CMB_GeV = T0_CMB_K * K_to_GeV;
static const double N_EFF = 3.046; // From arXiv:1606.06986

static const double OMEGA_C_h2_WMAP5 = 0.1081;
static const double OMEGA_B_h2_WMAP5 = 2.268/100;
static const double SPECTRAL_INDEX_WMAP5 = 0.961;
static const double SPECTRAL_AMPLITUDE_WMAP5 = 2.41e-9;


/*

  Julien definition


*/

///
/// Fundamental constants, masses and couplings.
/// (from CODATA 2010)
///

/// Newton constant [m^3/kg/s^2].
static const double G_N = 6.67408e-11;
/// Planck constant = h/2PI [J.s].
static const double HBAR = 1.054571726e-34;
/// Boltzmann constant [J/K].
static const double K_BOLTZ = 1.38066488e-23;
/// speed of light in vacuum [m/s].
static const double C_LIGHT = 2.99792458e+8;
/// speed of light in vacuum [cm/s].
static const double C_light = 1.e2*C_LIGHT;
/// permeability of free space [H/m = T m / A = m kg / (s2 A2) = m kg / C2].
static const double MU_0 = 4.e-7 * PI;
/// permittivity of free space [A^2.s^2.N^-1.m^-2].
static const double EPSILON_0 = 1./(MU_0 * C_LIGHT * C_LIGHT);
/// Fine structure constant.
static const double ALPHA_EM = pow(E_CHARGE,2)/(4.*PI*EPSILON_0*HBAR*C_LIGHT);
/// Strong interaction constant (at MZ).
static const double ALPHA_S = 0.1184;
/// Fermi constant GF/(hbar c)^3 in [GeV-2]
static const double FERMI_CST = 1.1663787e-5;
/// Weak interaction mixing angle.
static const double SIN2_THETA_W = 0.23146;
static const double SIN_THETA_W = sqrt(SIN2_THETA_W);
static const double THETA_W = asin(SIN_THETA_W);
/// Planck mass in international units [J.kg.s2/m2]^1/2.
static const double M_PLANCK_IU = sqrt(HBAR * C_LIGHT / G_N);
/// Avogadro constant
static const double N_AVOGADRO = 6.02214129e23;





///
/// Other useful quantities
///
/// Sun's mass : [kg]
static const double M_SUN = 1.3271244e20/G_N;
/// Sun's average radius [m]
static const double R_SUN = 6.955e8;
/// Earth's mass [kg]
static const double M_EARTH = 5.9736e24;
/// Earth's average radius [m]
static const double R_EARTH = 6.371e6;
/// Temperature of QCD Phase transition
static const double T_QCD_GeV = 0.150;






///
/// Reminder: magnetic fields are usually expressed in units of muG. Note that
/// 1 G = 1.e-4 T, and 1 T = 1 kg/(A s2) = 1 kg / (C s). This might be helpful
/// for further conversions.
///


///
/// Lepton/antilepton masses [GeV/c^2].
///

/// \f$ \nu_{e} \f$.
static const double M_nue = 0;
static const double W_nue = 0;
/// \f$ \nu_{\mu} \f$.
static const double M_numu = 0;
static const double W_numu = 0;
/// \f$ \nu_{\tau} \f$.
static const double M_nutau = 0;
static const double W_nutau = 0;
/// electron mass.
static const double M_e = 0.510998928e-3;
static const double W_e = 0;
/// muon mass.
static const double M_MU = 1.05658369e-01;
static const double W_MU = 2.996e-19;
/// \f$ \tau \f$ mass.
static const double M_TAU = 1.77699;
static const double W_TAU = 2.267e-12;

/// Charged leptons
static const double M_CLEP[3] = {M_e,M_MU,M_TAU};
/// Neutral leptons
static const double M_NLEP[3] = {M_nue,M_numu,M_nutau};
//static const double M_NLEP[3] = {0,0,0};

///
/// Quark pole masses [GeV/c^2].
///

/// up.
static const double M_QUARK_UP = 2.3e-03;
//static const double M_qu = 3.3e-1;
/// down.
static const double M_QUARK_DOWN = 4.8e-03;
//static const double M_qd = 3.3e-1;
/// strange.
//static const double M_qs = 1.0e-01;
static const double M_QUARK_STRANGE = 5.0e-01;
/// charm.
static const double M_QUARK_CHARM = 1.5;
/// bottom.
static const double M_QUARK_BOTTOM = 4.18;
/// top.
static const double M_QUARK_TOP = 172.9;

/// Up quarks
static const double M_QU[3] = {M_QUARK_UP,M_QUARK_CHARM,M_QUARK_TOP};
/// Down quarks
static const double M_QD[3] = {M_QUARK_DOWN,M_QUARK_STRANGE,M_QUARK_BOTTOM};


///
/// Massive boson masses and widths (standard model) [GeV/c^2].
///

/// Z.
static const double M_Z = 91.1876;
static const double W_Z = 2.4952;

/// W.
static const double M_W = 80.385;
static const double W_W = 2.085;

/// Higgs.
static const double M_H = 126.0;
static const double W_H = 1.0;


/// Weak interaction coupling.
static const double G_WEAK = sqrt(8*M_W*M_W*FERMI_CST/sqrt(2.0));


///
/// Hadron masses [GeV/c^2].
///

/// Proton/anti-proton mass.
static const double M_PROTON = 0.938272046;
/// Neutron mass.
static const double M_NEUTRON = 0.93956533;
/// Deuteron mass.
static const double M_deut = 1.87561282202674984;
/// charged pion mass [GeV/c^2].
static const double M_PI_C = 0.13957018;
static const double W_PI_C = 2.528e-17;
/// neutral pion mass [GeV/c^2].
static const double M_PI_0 = 0.1349766;
static const double W_PI_0 = 7.725e-9;
/// kaon charged mass [GeV/c^2].
static const double M_K_C = 0.493677;
static const double W_K_C = 5.317e-17;
/// kaon neutral mass [GeV/c^2].
static const double M_K_0 = 0.497614;
static const double W_K_0_S = 7.351e-15;
static const double M_K_0_L = 1.287e-17;
/// eta (meson) mass [GeV/c^2].
static const double M_ETA = 0.547853;
static const double W_ETA = 1.31e-6;
/// eta prime (meson) mass [GeV/c^2].
static const double M_ETA_P = 0.95778;
static const double W_ETA_P = 0.198e-3;
/// sigma or f0(600) (meson) mass [GeV/c^2].
static const double M_SIGMA = 0.6;
/// rho(770) (meson) mass [GeV/c^2].
static const double M_RHO770 = 0.77549;
/// omega(782) (meson) mass [GeV/c^2].
static const double M_OMEGA782 = 0.78265;
/// eta'(958) (meson) mass [GeV/c^2].
static const double M_ETAP958 = 0.95778;
/// f0(980) (meson) mass [GeV/c^2].
static const double M_F0980 = 0.980;
/// a0(980) (meson) mass [GeV/c^2].
static const double M_A0980 = 0.980;
/// phi(1020) (meson) mass [GeV/c^2].
static const double M_PHI1020 = 1.019455;
/// h1(1170) (meson) mass [GeV/c^2].
static const double M_H11170 = 1.170;
/// Delta+ mass [GeV/c^2].
static const double M_D1232 = 1.232;
/// Delta+ width (Breit-Wigner).
static const double G_D1232 = 0.120/2.;
/// b1(1235) mass [GeV/c^2].
static const double M_B11235 = 1.2295;


/// Unified atomic mass
static const double UAM = 1.660538921e-27/1.782661845e-27;
/// Atomic number to mass in GeV/c2
static const double A_to_GeV = UAM;

//static const double ATOMIC_WEIGHT[200]={
//  1,2,3,4,5,6,7,8,9,10,11,12,
//};

///---------------------------------------------------------------------
///
/// Unit conversions
///
///---------------------------------------------------------------------

static const double DEG_TO_RAD = PI/180.;
static const double RAD_TO_DEG = 180./PI;
static const double DEG_TO_ARCMIN = 60.;
static const double ARCMIN_TO_DEG = 1./DEG_TO_ARCMIN;
static const double ARCMIN_TO_RAD = ARCMIN_TO_DEG * DEG_TO_RAD;
static const double ARCMIN_TO_ARCSEC = 60.;
static const double ARCSEC_TO_RAD = 1./(RAD_TO_DEG * DEG_TO_ARCMIN * ARCMIN_TO_ARCSEC);


/// mb to cm^2 (1b = 10^-28 m^2).
static const double mb_to_cm2 = 1.e-27;
static const double cm2_to_mb = 1.e+27;
/// pb to cm^2 (1b = 10^-28 m^2).
static const double pb_to_cm2 = 1.e-36;
static const double cm2_to_pb = 1.e36;

/// eV/c^2 to kg.
static const double eV_to_kg = 1.782661845e-36;
static const double GeV_to_kg = 1.782661845e-27;
static const double kg_to_eV = 1./eV_to_kg;
static const double kg_to_GeV = 1./GeV_to_kg;
/// solar mass to kg.
static const double MSOL_to_kg = M_SUN;
static const double kg_to_MSOL = 1.0/MSOL_to_kg;
/// solar mass to GeV/c^2.
static const double MSOL_to_GeV = (MSOL_to_kg*kg_to_GeV);
/// GeV/c^2 to solar mass.
static const double GeV_to_MSOL = 1./MSOL_to_GeV;

/// mass/energy density conversion
static const double GeVpercm3_to_MSOLperkpc3 = GeV_to_MSOL * pow(kpc_to_cm,3);

/// Convert GM/r in km2/s2 when M is in Msun and r in kpc.
static const double GN_Msunkpcm2_to_km2sm2 = G_N*1.e-6*M_SUN*100*cm_to_kpc;


/// Earth mass
static const double MEARTH_to_kg = M_EARTH;

/// time to time.
static const double yr_to_sec = 365.25 * 24. * 60. * 60.;
static const double Myr_to_sec = 1.e6 * yr_to_sec;
static const double sec_to_yr = 1./yr_to_sec;
static const double sec_to_Myr = 1./Myr_to_sec;

/// eV to Joules.
//static const double eV_to_J = E_CHARGE;

//static const double J_to_eV = 1./eV_to_J;
static const double J_to_GeV = 1./GeV_to_J;

/// erg to GeV (1 erg = 1.e-7 J).
static const double erg_to_J = 1.0e-7;
static const double erg_to_GeV = 1.0e-7 * J_to_GeV;
static const double erg_to_eV = 1.0e-7 * J_to_eV;
static const double GeV_to_erg = 1.0/erg_to_GeV;
static const double eV_to_erg = 1.0/erg_to_eV;

/// Kelvin to eV (at 300 K).


/// eV to Hz
static const double Hz_to_eV = (2.*PI*HBAR) * J_to_eV;
static const double Hz_to_GeV = (2.*PI*HBAR) * J_to_GeV;
static const double eV_to_Hz = 1./Hz_to_eV;
static const double GeV_to_Hz = 1./Hz_to_GeV;

/// Magnetic field.
static const double Gauss_to_Tesla = 1.e-4;
static const double muGauss_to_Tesla = 1.e-10;
static const double Tesla_to_Gauss = 1.e4;
static const double Tesla_to_muGauss = 1.e10;

/// Jansky to Joules, erg, GeV --- [Jy] = [10^-26 J/s/m^2/Hz]
/// I guess I can define this as a (nu Fnu).
static const double Jansky_to_Jperm2sHz = 1.e-26;
static const double Jansky_to_ergpercm2sGeV = 
  1.e-4 * Jansky_to_Jperm2sHz * J_to_GeV * GeV_to_erg * GeV_to_Hz;
static const double Jansky_to_GeVpercm2sGeV = Jansky_to_ergpercm2sGeV * erg_to_GeV;


/// eV to 1/time [s^-1].
static const double eV_to_sm1 = eV_to_J/HBAR;
/// GeV to 1/time [s^-1].
static const double GeV_to_sm1 = GeV_to_J/HBAR;
/// time [s] to 1/eV [eV^-1].
static const double s_to_eVm1 = eV_to_sm1;
/// time [s] to 1/GeV [GeV^-1].
static const double s_to_GeVm1 = GeV_to_sm1;
/// eV to 1/distance [m^-1].
static const double eV_to_mm1 = eV_to_sm1/C_LIGHT;
/// GeV to 1/distance [m^-1].
static const double GeV_to_mm1 = GeV_to_sm1/C_LIGHT;
/// eV to 1/distance [cm^-1].
static const double eV_to_cmm1 = eV_to_sm1/C_light;
/// GeV to 1/distance [cm^-1].
static const double GeV_to_cmm1 = GeV_to_sm1/C_light;
/// distance [m] to 1/eV [eV^-1].
static const double m_to_eVm1 = eV_to_mm1;
/// distance [m] to 1/GeV [GeV^-1].
static const double m_to_GeVm1 = GeV_to_mm1;
/// distance [cm] to 1/eV [eV^-1].
static const double cm_to_eVm1 = eV_to_cmm1;
/// distance [cm] to 1/GeV [GeV^-1].
static const double cm_to_GeVm1 = GeV_to_cmm1;

/// Cross-section conversion [GeV^-2]->[pb].
static const double GeVm2_to_pb = pow(HBAR*J_to_GeV * C_light,2) * cm2_to_pb;
static const double pb_to_GeVm2 = 1./GeVm2_to_pb;
// (hbar c)^3 / hbar = 1.167e-17 [cm^3/s][GeV^2]
// => 1 GeV^-2 = 1.167e-17 [cm^3/s]
static const double GeVm2_to_cm3persec = 
  pow(HBAR*J_to_GeV * C_light,3)/ (HBAR * J_to_GeV);
static const double cm3persec_to_GeVm2 = 1./GeVm2_to_cm3persec;



///------------------------------------------------------------------------
///
/// Cross sections.
///

/// Classical electron radius [m].
static const double ELECTRON_RADIUS = ALPHA_EM * HBAR * J_to_GeV * C_LIGHT / M_e;
/// Thomson cross section [barn].
static const double SIGMA_THOMSON = 8.e28*PI*pow(ELECTRON_RADIUS,2)/3.;


/// Magnetic energy density
static const double muGauss2_to_GeVpercm3 = 1./(4./3. * SIGMA_THOMSON * 
					 1.e3 * mb_to_cm2 * C_light /
					 ( pow(8.*PI*M_e,2) * 2.53e-18));

///------------------------------------------------------------------------
///
/// Cosmological and galactic quantities.
///

/// Planck mass [J.kg.s2/m2]^1/2 -> GeV.
static const double M_PLANCK = M_PLANCK_IU * 
  sqrt(J_to_GeV * kg_to_GeV) * s_to_GeVm1 / m_to_GeVm1;


/// Hubble parameters (VL2)
//static const double HUBBLE_h = 0.73; // dimensionless
/// Hubble parameters (WMAP7)
//static const double HUBBLE_h = 0.71; // dimensionless
//static const double HUBBLE_H = 100 * HUBBLE_h; // [km/s/Mpc]
/// Hubble parameters (Planck 2015)
//static const double HUBBLE_h = 0.6781; // dimensionless
//static const double HUBBLE_h = 0.6774; // dimensionless

//static const double HUBBLE_H = 100 * HUBBLE_h; // [km/s/Mpc]



/// Closure density of the Universe \f$ \rho = \frac{3H^2}{8\pi G} \f$
/// [Msol/kpc^3].
//static const double rhoC = 3 * pow(HUBBLE_H * 1.e5*cm_to_kpc*1.e-3,2)/(8*PI*G_N) * 
 // pow(1.e-2*kpc_to_cm,3) * kg_to_MSOL;
//static const double rhoC_GeV_per_cm3 = rhoC * MSOL_to_GeV * pow(cm_to_kpc,3);
//static const double rhoC_GeV4 = rhoC_GeV_per_cm3 / pow(cm_to_GeVm1,3);

//static const double H0GeV = ((1e+2*HUBBLE_H) / (kpc_to_cm)) * (HBAR / GeV_to_J);
//static const double rhoC_GeV4Bis = (3.0/(8.0*PI))*pow(H0GeV,2)*pow(M_PLANCK,2);

/// CMB energy density in [GeV.cm^-3].
static const double U_CMB = PI*PI * pow(K_BOLTZ * T0_CMB_K * J_to_GeV,4) / 
  (15. * pow(HBAR * J_to_GeV * C_light,3));


/// Reduced cosmological matter density.
/// (from WMAP7)
//static const double Omega_M = 0.14170/pow(HUBBLE_h,2);   // Planck all included / previous value : 0.1334/pow(HUBBLE_h,2)
//static const double Omega_B = 2.230/pow(10*HUBBLE_h,2); // PLANCK all included / otherwise :  2.258/pow(10*HUBBLE_h,2);
// static const double Omega_NU = (M_nue+M_numu+M_nutau)/(93.14e-9*pow(HUBBLE_h,2));
// static const double Omega_CDM = Omega_M-Omega_B-Omega_NU;
//static const double Omega_Lambda = 1.0 - Omega_M;
/// (from Planck)
//static const double Omega_CDM = 0.1188/pow(HUBBLE_h,2);


// static const double N_EFF = 3.000;

//static const double Omega_Gamma = (2*((PI*PI)/30)*pow(T0_CMB_GeV,4)/rhoC_GeV4Bis);


/// Local magnetic field in muG.
//static const double B_MW = 50.;
//static const double B_MW = 10.;
//static const double B_MW = 3.;
static const double B_MW = 1.; // when averaged over a 1 kpc sphere
//static const double B_MW = 0.;
// test with Timur and Roberto
//static const double B_MW = 2.631;

/// Local magnetic fluctuation [muG]
static const double dB_MW = 0.;

/// Local magnetic energy density in [GeV/cm^3]: 
/// \f[ U_B = \frac{B^2}{2\mu_0}\;({\rm J/m^3}) \f].
/// \f[ 1\mu{\rm G} = 10^{-10}{\rm T}\f]
static const double U_BMW = 1.e-20*B_MW*B_MW/(2*MU_0) * J_to_GeV * 1.e-6;


/// Local IR + UV + starlight temperatures and energy densities [GeV/cm^3].
/// 3 different models.
const static int N_ISRF_MODEL = 3;

static int ISRF_MODEL_ID = 1;
// IR component.
//static const double T_IR_GeV[N_ISRF_MODEL] = {3.45e-12,2.85e-12,2.9e-12};
static const double T_IR_GeV[N_ISRF_MODEL] = {2.85e-12,2.85e-12,2.9e-12};
static const double T_IR_K[N_ISRF_MODEL] = {T_IR_GeV[0] * GeV_to_K , 
				     T_IR_GeV[1] * GeV_to_K ,
				     T_IR_GeV[2] * GeV_to_K};
static const double U_IR[N_ISRF_MODEL] = {1.e-10,2.54e-10,3.2117e-10};
// Optical Component
static const double T_SL_GeV[N_ISRF_MODEL] = {3.e-10,2.7e-11,2.7e-11};
static const double T_SL_K[N_ISRF_MODEL] = {T_SL_GeV[0] * GeV_to_K,
				     T_SL_GeV[1] * GeV_to_K,
				     T_SL_GeV[2] * GeV_to_K};
static const double U_SL[N_ISRF_MODEL] = {3.98e-10,5.47e-11,6.1992e-11};
// UV component
const static int N_UV_COMP = 3;
static const double T_UV_GeV[N_ISRF_MODEL][N_UV_COMP] = {{0,0,0},
						  {2.8e-10,5.3e-10,2.e-9},
						  {2.5e-10,4.8e-10,1.9e-9}};
static const double T_UV_K[N_ISRF_MODEL][N_UV_COMP] = 
  {{0,0,0},
   {T_UV_GeV[1][0]*GeV_to_K,T_UV_GeV[1][1]* GeV_to_K,T_UV_GeV[1][2]*GeV_to_K},
   {T_UV_GeV[2][0]*GeV_to_K,T_UV_GeV[2][1]* GeV_to_K,T_UV_GeV[2][2]*GeV_to_K}};
static const double U_UV[N_ISRF_MODEL][N_UV_COMP] = 
  {{0,0,0},
   {3.7e-10,2.29e-10,1.1885192e-10},
   {3.376475e-10,2.593e-10,1.02645e-10}};


/// If one wants to define his own ISRF model
/// Vector of BB temperatures
static std::vector<double> T_ISRF_GeV;
static std::vector<double> T_ISRF_K;
/// Vector of energy densities
static std::vector<double> U_ISRF;
/// Vector of normalisations assuming bb distributions. 
static std::vector<double> CORR_ISRF_BB;


/// Average dispersion velocity in the Galaxy
static const double BETA_MW = 220.e3/C_LIGHT;


///
/// Energy loss timescales for different processes in the Galaxy or elsewhere.
///

///
/// Electrons and positrons.
///

/// Inverse Compton on CMB.
/// \f$ dE/dt = - b_{0}^{\rm IC} E^2/E_{0} \f$ , where \f$ E_0 = 1 \f$ GeV.
static const double B0_IC_CMB = 4./3. * SIGMA_THOMSON * 1.e3 * mb_to_cm2 * 
  U_CMB * C_light / (M_e*M_e);
/// Associated characteristic time [s].
static const double TAU_IC_CMB = 1./B0_IC_CMB;

/// Inverse Compton on IR field.
/// \f$ dE/dt = - b_{0}^{\rm IR} E^2/E_{0} \f$ , where \f$ E_0 = 1 \f$ GeV.
static const double B0_IC_IR = 4./3. * SIGMA_THOMSON * 1.e3 * mb_to_cm2 * 
  U_IR[0] * C_light / (M_e*M_e);
/// Associated characteristic time [s].
static const double TAU_IC_IR = 1./B0_IC_IR;

/// Synchrotron on magnetic field.
/// \f$ dE/dt = - b_{0}^{\rm sync} E^2/E_{0} \f$ , where \f$ E_0 = 1 \f$ GeV.
static const double B0_SYNC = 2.53e-18 * B_MW * B_MW;
/// Associated characteristic time [s].
static const double TAU_SYNC = (B_MW > 0.0) ? 1./B0_SYNC : 1.e30;

/// Total timescale for E^-2 losses [s].
static const double TAU_TOT = 1./(B0_IC_CMB + B0_IC_IR + B0_SYNC);

#endif
