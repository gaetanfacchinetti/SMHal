#include "../headers/CPOddHiggsLowMass.h"

CPOddHiggsLowMass::CPOddHiggsLowMass(std::string const &n_name_file)
{
    V = 1. / (348 * sqrt(2));
    vv = 174;

    // Parameters for the mass matrix
    del = 0;
    del1 = 0.0083759;
    del2 = 0.0036189;  // AEP3Pi *)
    del3 = 0.00913525; // AEPi3PiC *)
    del4 = 0.00197181; // AEPPi3PiC *)
    del5 = 6.05444;    // AEPEPi3 *)
    del6 = 6.31207;    // AEPEPiC *)
    del7 = 6.05446;    // AEPEPPi3 *)
    del8 = 6.31257;    // AEPEPPiC *)
    del9 = 0.005064;   // mKz ^ 2 - mK ^ 2 + mpip ^ 2 - mpi3 ^ 2 *)
    Delta = -(1369. / 10000.);
    fpi = 93. / 1000.;
    mpi = (27. / 200.);
    mpim = 139. / 1000.;
    mpip = 139. / 1000.;
    mpi8 = (23. / 40.);
    mpi9 = (47. / 50.);
    meta = 548. / 1000.;
    metp = 958. / 1000.;
    mK = 494. / 1000.;
    mKz = 498. / 1000.;

    // Wess-Zumino-Witten couplings
    Cpi3 = -10.75;
    Ceta = -10.8;
    Cetp = 13.6;

    // Charges
    Qu = 2. / 3.;
    Qd = -1. / 3.;
    Qe = -1.;
    Qch = -1.;

    // Number of colors
    Nc = 3;

    // Coupling ... ?
    g = 0.3;

    // Theta_eta transformation angle between (pi8,pi9) and (eta, eta_prime)
    te = -0.226893;

    Npoints = 0;

    // ...
    Yepi3 = 3e-7;
    Ymueta = 2e-5;

    // Input file name
    name_file = n_name_file;


    GA1N1N1.resize(0);
    

    Read();

    OAA.resize(Npoints);
    OA3.resize(Npoints);
    OA8.resize(Npoints);
    OA9.resize(Npoints);
    OAeta.resize(Npoints);
    OAetp.resize(Npoints);

    AmpA1.resize(Npoints);
    AmpPi3.resize(Npoints);
    AmpEta.resize(Npoints);
    AmpEtap.resize(Npoints);


    // indices of mesons
    _pi3 = 0;
    _pip = 1;
    _pim = 2;
    _eta = 3;
    _etp = 4;
    _Kp = 5;
    _Km = 6;
    _Kz = 7;

    MM.resize(8);
    MM[_pi3] = mpi;
    MM[_pip] = mpip;
    MM[_pim] = mpim;
    MM[_eta] = meta;
    MM[_etp] = metp;
    MM[_Kp] = mK;
    MM[_Km] = mK;
    MM[_Kz] = mKz;

    //AmplitudeInitialise();

    // Diagonalisation of the mass matrix
    for (int i = 0; i < Npoints; i++)
        Diagonalise(i);

    // Initialise coupling DM - A1
    //for (int i = 0; i < Npoints; i++)
    //    ach0A[i] = OAA[i] * GA1N1N1[i];


    //int i = 3, j = 2, k = 1;
    //Order(i,j,k);
    //std::cout << i << " " << j << " " << k << std::endl;
}

void CPOddHiggsLowMass::Diagonalise(int i)
{

    double dm3 = -fpi * P11[i] * V * (mpi * mpi * ((1. / Tbeta[i]) - Tbeta[i]) + del * ((1. / Tbeta[i]) + Tbeta[i]));
    double dm8 = -fpi * P11[i] * V * (-sqrt(3) * Tbeta[i] * mpi8 * mpi8 + (mpi * mpi / sqrt(3)) * ((1. / Tbeta[i]) + 2 * Tbeta[i]) + (del / sqrt(3)) * ((1. / Tbeta[i]) - Tbeta[i]));
    dcomp dm9 = -fpi * P11[i] * V * (sqrt(3. / 2.) * Tbeta[i] * mpi8 * mpi8 + (mpi * mpi / sqrt(6)) * (2 * (1. / Tbeta[i]) + Tbeta[i]) + sqrt(2. / 3.) * del * ((1. / Tbeta[i]) - Tbeta[i])) - sqrt(2. / 3.) * fpi * Cg(i) * (mpi9 * mpi9 - (1. / 3.) * (mpi8 * mpi8 + mpi * mpi));

    //std::cout << Cg(i) << std::endl;
    //double xt = pow(M_QUARK_TOP / mA1[i], 2);
    //std::cout << FF(xt) << std::endl;

    // Initialisation of the mass matrix
    Eigen::Matrix4cd M = Eigen::Matrix4cd::Constant(4, 4, 0);
    M(0, 0) = mpi * mpi;
    M(0, 1) = sqrt(1. / 3.) * del;
    M(0, 2) = sqrt(2. / .3) * del;
    M(0, 3) = dm3;
    M(1, 0) = M(0, 1);
    M(1, 1) = mpi8 * mpi8;
    M(1, 2) = Delta;
    M(1, 3) = dm8;
    M(2, 0) = M(0, 2);
    M(2, 1) = M(1, 2);
    M(2, 2) = mpi9 * mpi9;
    M(2, 3) = dm9;
    M(3, 0) = M(0, 3);
    M(3, 1) = M(1, 3);
    M(3, 2) = M(2, 3);
    M(3, 3) = mA1[i] * mA1[i];
    //std::cout << M << std::endl;


    // Diagonalisation of the mass matrix
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4cd> eigensolver(M);
    if (eigensolver.info() != Eigen::Success)
        abort();

    //std::cout << "The eigenvalues of M are:\n"
    //          << eigensolver.eigenvalues() << std::endl;
    //          << eigensolver.eigenvectors().transpose() << std::endl;

    // Eigenvectors are stacked in column of the eigen matrix
    //Eigen::Vector4d eigval = eigensolver.eigenvalues();
    Eigen::Matrix4cd eigMat = eigensolver.eigenvectors();
    Eigen::Matrix4d eigMat_T; // = eigensolver.eigenvectors().transpose();

    for(int i = 0; i<4; i++)
        for(int j = 0; j<4; j++)
            eigMat_T(i,j) =  (eigensolver.eigenvectors().transpose()(i,j)).real();



    /*Eigen::Matrix4d diagMat = Eigen::Matrix4d::Constant(4, 4, 0);
    diagMat(0, 0) = eigval(0);
    diagMat(1, 1) = eigval(1);
    diagMat(2, 2) = eigval(2);
    diagMat(3, 3) = eigval(3);
    Eigen::Vector4d u(0, 0, 0, 1);
    Eigen::Vector4d OA = eigMat.transpose() * u;
    Eigen::Matrix4d A = OA * (OA.transpose());

    Eigen::Matrix2d N = Eigen::Matrix2d::Constant(2, 2, 0);
    N(0, 0) = mpi8 * mpi8;
    N(0, 1) = Delta;
    N(1, 0) = Delta;
    N(1, 1) = mpi9 * mpi9;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolve(N);
    if (eigensolve.info() != Eigen::Success)
        abort();
    
    std::cout << std::endl <<  "The eigenvalues of N are:\n"
              << eigensolve.eigenvalues() << std::endl
              << eigensolve.eigenvectors()<< std::endl << std::endl;*/

    // Associate the good line of the matrix
    int line_pi3 = 0;
    int line_A1 = 0;
    int line_eta = 0;
    int line_etap = 0;

    for (int i = 0; i < 4; i++)
    {
        if (abs(eigMat_T(i, 0)) > abs(eigMat_T(line_pi3, 0)))
            line_pi3 = i;
    }

    for (int i = 0; i < 4; i++)
    {
        if (abs(eigMat_T(i, 3)) > abs(eigMat_T(line_A1, 3)))
            line_A1 = i;
    }

    for (int i = 0; i < 4; i++)
    {
        if (abs(abs(eigMat_T(i, 1)) - cos(te)) < abs(abs(eigMat_T(line_eta, 1)) - cos(te)))
            line_eta = i;
    }

    for (int i = 0; i < 4; i++)
    {
        if (abs(abs(eigMat_T(i, 2)) - cos(te)) < abs(abs(eigMat_T(line_etap, 2)) - cos(te)))
            line_etap = i;
    }

    // Change the signs
    if (eigMat_T(line_pi3, 0) < 0)
        for (int j = 0; j < 4; j++)
            eigMat_T(line_pi3, j) = -eigMat_T(line_pi3, j);

    if (eigMat_T(line_eta, 1) < 0)
        for (int j = 0; j < 4; j++)
            eigMat_T(line_eta, j) = -eigMat_T(line_eta, j);

    if (eigMat_T(line_etap, 2) < 0)
        for (int j = 0; j < 4; j++)
            eigMat_T(line_etap, j) = -eigMat_T(line_etap, j);

    if (eigMat_T(line_A1, 3) < 0)
        for (int j = 0; j < 4; j++)
            eigMat_T(line_A1, j) = -eigMat_T(line_A1, j);

    Eigen::Matrix4d eigMat_T_ordered;

    for (int j = 0; j < 4; j++)
    {
        eigMat_T_ordered(0, j) = eigMat_T(line_pi3, j);
        eigMat_T_ordered(1, j) = eigMat_T(line_eta, j);
        eigMat_T_ordered(2, j) = eigMat_T(line_etap, j);
        eigMat_T_ordered(3, j) = eigMat_T(line_A1, j);
    }

    //std::cout << "ici " << std::endl << eigMat_T_ordered << std::endl;

    OAA[i] = eigMat_T_ordered(3, 3);
    OA3[i] = eigMat_T_ordered(3, 0);
    OAeta[i] = cos(te) * eigMat_T_ordered(3, 1) - sin(te) * eigMat_T_ordered(3, 2);
    OAetp[i] = sin(te) * eigMat_T_ordered(3, 1) + cos(te) * eigMat_T_ordered(3, 2);
    OA8[i] = eigMat_T_ordered(3, 1);
    OA9[i] = eigMat_T_ordered(3, 2);

    //std::cout << "everything is fine here " << OAA[0] << " " << GA1N1N1[0] << std::endl;

    //std::cout << OAA << " " << OA3 << " " << OAeta << " " << OAetp << " " << OA8 << " " << OA9 << std::endl;

    // ATTENTION the eigenvectors may not be ordered as they should in the matrix

    //std::cout << A(3,3) << std::endl;
    //std::cout  << std::endl << eigMat*diagMat*eigMat.transpose() << std::endl;
}

dcomp CPOddHiggsLowMass::FF(double x)
{
    dcomp _x = x;
    dcomp res = 2.0 * _x * pow(log((sqrt(1.0 - 4.0 * _x) - 1.0) / (sqrt(1.0 - 4.0 * _x) + 1.0)), 2);

    return res;
}

dcomp CPOddHiggsLowMass::Cg(int i)
{
    double xt = pow(M_QUARK_TOP / mA1[i], 2);
    double xc = pow(M_QUARK_CHARM / mA1[i], 2);
    double xb = pow(M_QUARK_BOTTOM / mA1[i], 2);

    dcomp res = -0.5 * P11[i] * V * ((1. / Tbeta[i]) * (FF(xt) + FF(xc)) + Tbeta[i] * FF(xb));

    return res;
}

// Coupling between A1_tilde and photons
dcomp CPOddHiggsLowMass::sFg(int i, double s, double T_over_TQCD)
{
    // SHOULD WE REALLY SWITCH OFF THE TEMPERATURE DEPENDANCE ? 
    //if (T_over_TQCD <= 1)
    //std::cout << abs(sCgam(i, s) + sdCgam(i, s)) << " " << Cpi3 << std::endl;
    return OAA[i] * (sCgam(i, s) + sdCgam(i, s)) + OA3[i] * Cpi3 + OAeta[i] * Ceta + OAetp[i] * Cetp;

    //else
    //    return sCgam(i, s) + sdCgam(i, s);
}

// Coupling between A1_tilde and photons
dcomp CPOddHiggsLowMass::sFg_test(int i, double s, double T_over_TQCD)
{
        //return OA3[i] * Cpi3;
        //std::cout << OA3[i] << "\t" << Cpi3 << std::endl;
        return sCgam(i, s) + sdCgam(i, s);
        exit(0);
}


// Coupling between A1 and photons through quarks and non SM particles
dcomp CPOddHiggsLowMass::sCgam(int i, double s)
{

    double xt = pow(M_QUARK_TOP, 2) / s;
    double xc = pow(M_QUARK_CHARM, 2) / s;
    double xb = pow(M_QUARK_BOTTOM, 2) / s;
    double xta = pow(M_TAU, 2) / s;
    double xch1 = pow(mch1[i], 2) / s;
    double xch2 = pow(mch2[i], 2) / s;

    return -P11[i] * V * (Nc * Qu * Qu * Cotb[i] * (FF(xt) + FF(xc)) + Nc * Qd * Qd * Tbeta[i] * FF(xb) + Qe * Qe * Tbeta[i] * FF(xta)) - ((Qch * Qch) / (2 * sqrt(2) * mch1[i])) * FF(xch1) * (lamb[i] * P12[i] * U12[i] * V12[i] - g * P11[i] * (Sbeta[i] * U12[i] * V11[i] + Cbeta[i] * U11[i] * V12[i])) - ((Qch * Qch) / (2 * sqrt(2) * mch2[i])) * FF(xch2) * (lamb[i] * P12[i] * U22[i] * V22[i] - g * P11[i] * (Sbeta[i] * U22[i] * V21[i] + Cbeta[i] * U21[i] * V22[i]));
}

// Coupling between A1 and photons through leptons
dcomp CPOddHiggsLowMass::sdCgam(int i, double s)
{

    double xe = pow(M_e, 2) / s;
    double xmu = pow(M_MU, 2) / s;

    return -P11[i] * V * Qe * Qe * Tbeta[i] * (FF(xe) + FF(xmu));
}

// coupling of A1 (tilde) to dark matter -- through A1 (not tilde)
double CPOddHiggsLowMass::ach0A(int i, double T_over_TQCD)
{
    //std::cout << "[in ach0A] " << GA1N1N1.size() << " " << OAA.size() <<  std::endl;
    //std::cout << "GA1N1N1 " << GA1N1N1[i] << std::endl;
    //std::cout << "OAA " << OAA[i] << std::endl;

    //if(T_over_TQCD <=1.)
        return OAA[i] * GA1N1N1[i];
    //else
    //    return GA1N1N1[i];
}


double CPOddHiggsLowMass::GammaTAS(int i, double T_over_TQCD)
{
    if (mA1[i] <= 1.2)
        return GamaL(i, T_over_TQCD) + GamaAch(i, T_over_TQCD) + GamaAt1S(i, T_over_TQCD) + Gamach0(i, T_over_TQCD);
    else
        return 0;
}

double CPOddHiggsLowMass::GamaAt1S(int i, double T_over_TQCD)
{
    // NOTE : In Gilbert's code mixing coefficients were the approximated one i.e. OAA->OAA0, ...
    //if (T_over_TQCD <= 1)
        return (ALPHA_EM * ALPHA_EM * pow(mA1[i], 3) * pow(abs(OAA[i] * (sCgam(i, mA1[i] * mA1[i]) + sdCgam(i, mA1[i] * mA1[i])) + OA3[i] * Cpi3 + OAeta[i] * Ceta + OAetp[i] * Cetp), 2)) / (64 * PI * PI * PI);
    //else
    //    return (ALPHA_EM * ALPHA_EM * pow(mA1[i], 3) * pow(abs(sCgam(i, mA1[i] * mA1[i]) + sdCgam(i, mA1[i] * mA1[i])), 2)) / (64 * PI * PI * PI);
}

double CPOddHiggsLowMass::GamaL(int i, double T_over_TQCD)
{
    if (mA1[i] <= 0.21)
        return GamaLe(i, T_over_TQCD);
    else
        return GamaLe(i, T_over_TQCD) + GamaLmu(i, T_over_TQCD);
}

double CPOddHiggsLowMass::GamaLe(int i, double T_over_TQCD)
{
    //if (T_over_TQCD <= 1)
        return pow(abs(OAA[i] * YllA1e[i] + OA3[i] * Yepi3), 2) * mA1[i] * sqrt(1 - 4 * M_e * M_e / (mA1[i] * mA1[i])) / (8 * PI);
    //else
     //   return pow(abs(YllA1e[i]), 2) * mA1[i] * sqrt(1 - 4 * M_e * M_e / (mA1[i] * mA1[i])) / (8 * PI);
}
double CPOddHiggsLowMass::GamaLmu(int i, double T_over_TQCD)
{
    //if (T_over_TQCD <= 1)
        return pow(abs(OAA[i] * YllA1mu[i] + OAeta[i] * Ymueta), 2) * mA1[i] * sqrt(1 - 4 * M_MU * M_MU / (mA1[i] * mA1[i])) / (8 * PI);
    //else
    //    return pow(abs(OAA[i] * YllA1mu[i]), 2) * mA1[i] * sqrt(1 - 4 * M_MU * M_MU / (mA1[i] * mA1[i])) / (8 * PI);
}

double CPOddHiggsLowMass::Gamach0(int i, double T_over_TQCD)
{
    if (mA1[i] <= 2 * mch0[i])
        return 0;
    else
        return pow(abs(ach0A(i, T_over_TQCD)), 2) * mA1[i] * sqrt(1 - 4 * mch0[i] * mch0[i] / (mA1[i] * mA1[i])) / (8 * PI) / 2;
 
    // The last division by two comes from Majorana case
}

double CPOddHiggsLowMass::GamaAch(int i, double T_over_TQCD)
{
    //if(T_over_TQCD <=1.)
        return GamaA(i, _pip, _pi3, _pim) + GamaA(i, _pi3, _pi3, _pi3) + GamaA(i, _eta, _pip, _pim) + GamaA(i, _eta, _pi3, _pi3) + GamaA(i, _etp, _pi3, _pi3) + GamaA(i, _etp, _pip, _pim) + GamaA(i, _pi3, _eta, _eta) + GamaA(i, _pi3, _eta, _etp) + GamaA(i, _pi3, _etp, _etp) + GamaA(i, _pi3, _Km, _Kp) + GamaA(i, _eta, _Kp, _Km) + GamaA(i, _etp, _Kp, _Km) + GamaA(i, _pip, _Km, _Kz) + GamaA(i, _pim, _Kp, _Kz);
    //else 
    //    return 0;
}

double CPOddHiggsLowMass::GamaA(int i, int mes1, int mes2, int mes3)
{

    if (mA1[i] < MM[mes1] + MM[mes2] + MM[mes3])
        return 0;

    std::vector<double> ismes1mes2mes3;
    ismes1mes2mes3.resize(5);
    ismes1mes2mes3[0] = (double)i;
    ismes1mes2mes3[2] = (double)mes1;
    ismes1mes2mes3[3] = (double)mes2;
    ismes1mes2mes3[4] = (double)mes3;

    double s_min = pow(MM[mes2] + MM[mes3], 2);
    double s_max = pow(mA1[i] - MM[mes1], 2);

    double res_int = GaussLegendre_Integral_Static(1, 5, gl100, s_min, s_max, ismes1mes2mes3, this, CallBack_fToIntForGamaA);
    double res = 1. / (256. * Sfactor(mes1, mes2, mes3) * PI * PI * PI * mA1[i]) * pow(AmpAf(i, mes1, mes2, mes3), 2) * res_int;

    return res;
}

double CPOddHiggsLowMass::fToIntForGamaA(std::vector<double> ismes1mes2mes3)
{
    int i = (int)ismes1mes2mes3[0];
    double s = ismes1mes2mes3[1];
    int mes1 = (int)ismes1mes2mes3[2];
    int mes2 = (int)ismes1mes2mes3[3];
    int mes3 = (int)ismes1mes2mes3[4];
    double mP = mA1[i];

    return sqrt(1 - 2 * (MM[mes2] * MM[mes2] + MM[mes3] * MM[mes3]) / s + pow(MM[mes2] * MM[mes2] - MM[mes3] * MM[mes3], 2) / (s * s)) * sqrt(pow(1 + (s - MM[mes1] * MM[mes1]) / (mP * mP), 2) - 4 * s / (mP * mP));
}

double CPOddHiggsLowMass::AmpAf(int i, int mes1, int mes2, int mes3)
{
    return OAA[i] * AmpA1f(i, mes1, mes2, mes3) + OA3[i] * AmpPi3f(i, mes1, mes2, mes3) + OAeta[i] * AmpEtaf(i, mes1, mes2, mes3) + OAetp[i] * AmpEtapf(i, mes1, mes2, mes3);
}

// Amplitudes for the decay of meson into other mesons
void CPOddHiggsLowMass::AmplitudeInitialise()
{

    // std::cout << AmpA1.size() << std::endl;
    for (int i = 0; i < Npoints; i++)
    {
        AmpA1[i].resize(8);
        AmpPi3[i].resize(8);
        AmpEta[i].resize(8);
        AmpEtap[i].resize(8);
        for (int j = 0; j < 8; j++)
        {
            AmpA1[i][j].resize(8);
            AmpPi3[i][j].resize(8);
            AmpEta[i][j].resize(8);
            AmpEtap[i][j].resize(8);
            for (int k = 0; k < 8; k++)
            {
                AmpA1[i][j][k].resize(8);
                AmpPi3[i][j][k].resize(8);
                AmpEta[i][j][k].resize(8);
                AmpEtap[i][j][k].resize(8);
            }
        }
    }

    // 0 -> pi3
    // 1 -> pi+
    // 2 -> pi-
    // 3 -> eta
    // 4 -> etap
    // 5 -> K+
    // 6 -> K-
    // 7 -> K0 / (bar{K0})

    // Amplitudes for the desintegration of the A1
    for (int i = 0; i < Npoints; i++)
    {
        AmpA1[i][_pi3][_pi3][_pi3] = -P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i]));
        AmpA1[i][_pi3][_pi3][_eta] = -P11[i] / (12 * sqrt(2) * vv * fpi) * (sqrt(3) * mpi * mpi * (Cotb[i] + Tbeta[i]) + del * (Cotb[i] - Tbeta[i])) * (cos(te) - sqrt(2) * sin(te));
        AmpA1[i][_pi3][_pi3][_etp] = -P11[i] / (12 * sqrt(2) * vv * fpi) * (sqrt(3) * mpi * mpi * (Cotb[i] + Tbeta[i]) + del * (Cotb[i] - Tbeta[i])) * (sin(te) + sqrt(2) * cos(te));
        AmpA1[i][_pi3][_pip][_pim] = 2 * AmpA1[i][_pi3][_pi3][_pi3];
        AmpA1[i][_pip][_pim][_eta] = 2 * AmpA1[i][_pi3][_pi3][_eta];
        AmpA1[i][_pip][_pim][_etp] = 2 * AmpA1[i][_pi3][_pi3][_etp];
        AmpA1[i][_pi3][_Kp][_Km] = -P11[i] / (6 * sqrt(2) * vv * fpi) * (mK * mK * (2 * Cotb[i] + Tbeta[i]) + (mpi * mpi - mK * mK) * (2 * Cotb[i] - Tbeta[i]));
        AmpA1[i][_pi3][_Kz][_Kz] = -P11[i] / (6 * sqrt(2) * vv * fpi) * (mK * mK - mpi * mpi - 3 * mKz * mKz) * Tbeta[i];
        AmpA1[i][_pip][_Km][_Kz] = -P11[i] / (6 * vv * fpi) * ((mK * mK + mpi * mpi) * Cotb[i] - mKz * mKz * (Cotb[i] - 2 * Tbeta[i]));
        AmpA1[i][_pim][_Kp][_Kz] = -P11[i] / (6 * vv * fpi) * ((mK * mK + mpi * mpi) * Cotb[i] - mKz * mKz * (Cotb[i] - 2 * Tbeta[i]));
        AmpA1[i][_pi3][_eta][_eta] = -P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i])) * pow(cos(te) - sqrt(2) * sin(te), 2);
        AmpA1[i][_pi3][_eta][_etp] = -P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i])) * 2 * (cos(te) - sqrt(2) * sin(te)) * (sin(te) + sqrt(2) * cos(te));
        AmpA1[i][_pi3][_etp][_etp] = -P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i])) * pow(sin(te) + sqrt(2) * cos(te), 2);
    }

    // Amplitudes for the decay of the pion 3 pi0
    for (int i = 0; i < Npoints; i++)
    {
        AmpPi3[i][_pi3][_pi3][_pi3] = mpi * mpi / (24 * fpi * fpi);
        AmpPi3[i][_pi3][_pip][_pim] = 4 * AmpPi3[i][_pi3][_pi3][_pi3];
        AmpPi3[i][_pi3][_pi3][_eta] = del / (6 * sqrt(3) * fpi * fpi) * (cos(te) - sqrt(2) * sin(te));
        AmpPi3[i][_pi3][_pi3][_etp] = del / (6 * sqrt(3) * fpi * fpi) * (sin(te) + sqrt(2) * cos(te));
        AmpPi3[i][_pip][_pim][_eta] = 2 * AmpPi3[i][_pi3][_pi3][_eta];
        AmpPi3[i][_pip][_pim][_etp] = 2 * AmpPi3[i][_pi3][_pi3][_etp];
        AmpPi3[i][_pi3][_eta][_eta] = mpi * mpi / (12 * fpi * fpi) * pow(cos(te) - sqrt(2) * sin(te), 2);
        AmpPi3[i][_pi3][_eta][_etp] = mpi * mpi / (12 * fpi * fpi) * 2 * (cos(te) - sqrt(2) * sin(te)) * (sin(te) + sqrt(2) * cos(te));
        AmpPi3[i][_pi3][_etp][_etp] = mpi * mpi / (12 * fpi * fpi) * pow(sin(te) + sqrt(2) * cos(te), 2);
        AmpPi3[i][_pi3][_Kp][_Km] = 1. / (12 * fpi * fpi) * (2 * mK * mK - mKz * mKz + mpi * mpi);
        AmpPi3[i][_pi3][_Kz][_Kz] = 1. / (12 * fpi * fpi) * (2 * mKz * mKz - mK * mK + mpi * mpi);
        AmpPi3[i][_pip][_Km][_Kz] = 1. / (12 * fpi * fpi) * sqrt(2) * del;
        AmpPi3[i][_pim][_Kp][_Kz] = 1. / (12 * fpi * fpi) * sqrt(2) * del;
        AmpPi3[i][_eta][_Kp][_Km] = 1. / (6 * sqrt(3) * fpi * fpi) * (-3 * sqrt(2) * mK * mK * sin(te) + (mpi * mpi - mKz * mKz) * (cos(te) - sqrt(2) * sin(te)));
        AmpPi3[i][_etp][_Kp][_Km] = 1. / (6 * sqrt(3) * fpi * fpi) * (3 * sqrt(2) * mK * mK * cos(te) + (mpi * mpi - mKz * mKz) * (sin(te) + sqrt(2) * cos(te)));
        AmpPi3[i][_eta][_Kz][_Kz] = 1. / (6 * sqrt(3) * fpi * fpi) * (3 * sqrt(2) * mKz * mKz * sin(te) - (mpi * mpi - mK * mK) * (cos(te) - sqrt(2) * sin(te)));
        AmpPi3[i][_etp][_Kz][_Kz] = 1. / (6 * sqrt(3) * fpi * fpi) * (-3 * sqrt(2) * mKz * mKz * cos(te) + (mpi * mpi - mK * mK) * (sin(te) + sqrt(2) * cos(te)));
        AmpPi3[i][_eta][_eta][_eta] = del / (18 * sqrt(3) * fpi * fpi) * pow(cos(te) - sqrt(2) * sin(te), 3);
        AmpPi3[i][_eta][_eta][_etp] = del / (18 * sqrt(3) * fpi * fpi) * 3 * pow(cos(te) - sqrt(2) * sin(te), 2) * pow(sin(te) + sqrt(2) * cos(te), 1);
        AmpPi3[i][_eta][_etp][_etp] = del / (18 * sqrt(3) * fpi * fpi) * 3 * pow(cos(te) - sqrt(2) * sin(te), 1) * pow(sin(te) + sqrt(2) * cos(te), 2);
        AmpPi3[i][_etp][_etp][_etp] = del / (18 * sqrt(3) * fpi * fpi) * pow(sin(te) + sqrt(2) * cos(te), 3);
    }

    for (int i = 0; i < Npoints; i++)
    {
        AmpEta[i][_pi3][_pi3][_pi3] = AmpPi3[i][_pi3][_pi3][_eta];
        AmpEta[i][_pi3][_pip][_pim] = AmpPi3[i][_pip][_pim][_eta];
        AmpEta[i][_pi3][_pi3][_eta] = AmpPi3[i][_pi3][_eta][_eta];
        AmpEta[i][_pi3][_pi3][_etp] = AmpPi3[i][_pi3][_eta][_etp];
        AmpEta[i][_pip][_pim][_eta] = 2 * AmpEta[i][_pi3][_pi3][_eta];
        AmpEta[i][_pip][_pim][_etp] = 2 * AmpEta[i][_pi3][_pi3][_etp];
        AmpEta[i][_pi3][_Kp][_Km] = AmpPi3[i][_eta][_Kp][_Km];
        AmpEta[i][_pi3][_Kz][_Kz] = AmpPi3[i][_eta][_Kz][_Kz];
        AmpEta[i][_pip][_Km][_Kz] = 1. / (6 * sqrt(6) * fpi * fpi) * ((2 * mpi * mpi - mK * mK - mKz * mKz) * cos(te) - 2 * sqrt(2) * (mKz * mKz + mK * mK + mpi * mpi) * sin(te));
        AmpEta[i][_pim][_Kp][_Kz] = AmpEta[i][_pip][_Km][_Kz];
        AmpEta[i][_pi3][_eta][_eta] = AmpPi3[i][_eta][_eta][_eta];
        AmpEta[i][_pi3][_eta][_etp] = AmpPi3[i][_eta][_eta][_etp];
        AmpEta[i][_pi3][_etp][_etp] = AmpPi3[i][_eta][_etp][_etp];
    }

    for (int i = 0; i < Npoints; i++)
    {
        AmpEtap[i][_pi3][_pi3][_pi3] = AmpPi3[i][_pi3][_pi3][_etp];
        AmpEtap[i][_pi3][_pip][_pim] = AmpPi3[i][_pip][_pim][_etp];
        AmpEtap[i][_pi3][_pi3][_eta] = AmpPi3[i][_pi3][_eta][_etp];
        AmpEtap[i][_pi3][_pi3][_etp] = AmpPi3[i][_pi3][_etp][_etp];
        AmpEtap[i][_pip][_pim][_eta] = 2 * AmpEtap[i][_pi3][_pi3][_eta];
        AmpEtap[i][_pip][_pim][_etp] = 2 * AmpEtap[i][_pi3][_pi3][_etp];
        AmpEtap[i][_pi3][_Kp][_Km] = AmpPi3[i][_etp][_Kp][_Km];
        AmpEtap[i][_pi3][_Kz][_Kz] = AmpPi3[i][_etp][_Kz][_Kz];
        AmpEtap[i][_pip][_Km][_Kz] = 1. / (6 * sqrt(6) * fpi * fpi) * ((2 * mpi * mpi - mK * mK - mKz * mKz) * sin(te) + 2 * sqrt(2) * (mKz * mKz + mK * mK + mpi * mpi) * cos(te));
        AmpEtap[i][_pim][_Kp][_Kz] = AmpEtap[i][_pip][_Km][_Kz];
        AmpEtap[i][_pi3][_eta][_eta] = AmpPi3[i][_eta][_eta][_etp];
        AmpEtap[i][_pi3][_eta][_etp] = AmpPi3[i][_eta][_etp][_etp];
        AmpEtap[i][_pi3][_etp][_etp] = AmpPi3[i][_etp][_etp][_etp];
    }

    //std::cout << Npoints << std::endl;
    //AmpA1f(0, _pim,0, 1);
}

double CPOddHiggsLowMass::AmpA1f(int i, int mes1, int mes2, int mes3)
{
    Order(mes1, mes2, mes3);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _pi3)
        return -3 * 2 * P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i]));

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _eta)
        return -2 * P11[i] / (12 * sqrt(2) * vv * fpi) * (sqrt(3) * mpi * mpi * (Cotb[i] + Tbeta[i]) + del * (Cotb[i] - Tbeta[i])) * (cos(te) - sqrt(2) * sin(te));

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _etp)
        return -2 * P11[i] / (12 * sqrt(2) * vv * fpi) * (sqrt(3) * mpi * mpi * (Cotb[i] + Tbeta[i]) + del * (Cotb[i] - Tbeta[i])) * (sin(te) + sqrt(2) * cos(te));

    if (mes1 == _pi3 && mes2 == _pip && mes3 == _pim)
        return -2 * P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i]));

    if (mes1 == _pip && mes2 == _pim && mes3 == _eta)
        return -2 * P11[i] / (12 * sqrt(2) * vv * fpi) * (sqrt(3) * mpi * mpi * (Cotb[i] + Tbeta[i]) + del * (Cotb[i] - Tbeta[i])) * (cos(te) - sqrt(2) * sin(te));

    if (mes1 == _pip && mes2 == _pim && mes3 == _etp)
        return -2 * P11[i] / (12 * sqrt(2) * vv * fpi) * (sqrt(3) * mpi * mpi * (Cotb[i] + Tbeta[i]) + del * (Cotb[i] - Tbeta[i])) * (sin(te) + sqrt(2) * cos(te));

    if (mes1 == _pi3 && mes2 == _Kp && mes3 == _Km)
        return -P11[i] / (6 * sqrt(2) * vv * fpi) * (mK * mK * (2 * Cotb[i] + Tbeta[i]) + (mpi * mpi - mK * mK) * (2 * Cotb[i] - Tbeta[i]));

    if (mes1 == _pi3 && mes2 == _Kz && mes3 == _Kz)
        return -P11[i] / (6 * sqrt(2) * vv * fpi) * (mK * mK - mpi * mpi - 3 * mKz * mKz) * Tbeta[i];

    if (mes1 == _pip && mes2 == _Km && mes3 == _Kz)
        return -P11[i] / (6 * vv * fpi) * ((mK * mK + mpi * mpi) * Cotb[i] - mKz * mKz * (Cotb[i] - 2 * Tbeta[i]));

    if (mes1 == _pim && mes2 == _Kp && mes3 == _Kz)
        return -P11[i] / (6 * vv * fpi) * ((mK * mK + mpi * mpi) * Cotb[i] - mKz * mKz * (Cotb[i] - 2 * Tbeta[i]));

    if (mes1 == _pi3 && mes2 == _eta && mes3 == _eta)
        return -2 * P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i])) * pow(cos(te) - sqrt(2) * sin(te), 2);

    if (mes1 == _pi3 && mes2 == _eta && mes3 == _etp)
        return -P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i])) * 2 * (cos(te) - sqrt(2) * sin(te)) * (sin(te) + sqrt(2) * cos(te));

    if (mes1 == _pi3 && mes2 == _etp && mes3 == _etp)
        return -2 * P11[i] / (12 * sqrt(2) * vv * fpi) * (mpi * mpi * (Cotb[i] - Tbeta[i]) + del * (Cotb[i] + Tbeta[i])) * pow(sin(te) + sqrt(2) * cos(te), 2);

    return 0;
}

double CPOddHiggsLowMass::AmpPi3f(int i, int mes1, int mes2, int mes3)
{
    Order(mes1, mes2, mes3);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _pi3)
        return 4 * 3 * 2 * mpi * mpi / (24 * fpi * fpi);

    if (mes1 == _pi3 && mes2 == _pip && mes3 == _pim)
        return 2 * 4 * mpi * mpi / (24 * fpi * fpi);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _eta)
    {
        //std::cout << "also here " << 3 * 2 * del / (6 * sqrt(3) * fpi * fpi) * (cos(te) - sqrt(2) * sin(te)) << std::endl;
        return 3 * 2 * del1 / (6 * sqrt(3) * fpi * fpi) * (cos(te) - sqrt(2) * sin(te));
    }

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _etp)
        return 3 * 2 * del2 / (6 * sqrt(3) * fpi * fpi) * (sin(te) + sqrt(2) * cos(te));

    if (mes1 == _pip && mes2 == _pim && mes3 == _eta)
        return 2 * del3 / (6 * sqrt(3) * fpi * fpi) * (cos(te) - sqrt(2) * sin(te));

    if (mes1 == _pip && mes2 == _pim && mes3 == _etp)
        return 2 * del4 / (6 * sqrt(3) * fpi * fpi) * (sin(te) + sqrt(2) * cos(te));

    if (mes1 == _pi3 && mes2 == _eta && mes3 == _eta)
        return 2 * 2 * mpi * mpi / (12 * fpi * fpi) * pow(cos(te) - sqrt(2) * sin(te), 2);

    if (mes1 == _pi3 && mes2 == _eta && mes3 == _etp)
        return 2 * del5 * mpi * mpi / (12 * fpi * fpi) * 2 * (cos(te) - sqrt(2) * sin(te)) * (sin(te) + sqrt(2) * cos(te));

    if (mes1 == _pi3 && mes2 == _etp && mes3 == _etp)
        return 2 * 2 * del7 * mpi * mpi / (12 * fpi * fpi) * pow(sin(te) + sqrt(2) * cos(te), 2);

    if (mes1 == _pi3 && mes2 == _Kp && mes3 == _Km)
        return 2 * 1. / (12 * fpi * fpi) * (2 * mK * mK - mKz * mKz + mpi * mpi);

    if (mes1 == _pi3 && mes2 == _Kz && mes3 == _Kz)
        return 2 * 1. / (12 * fpi * fpi) * (2 * mKz * mKz - mK * mK + mpi * mpi);

    if (mes1 == _pip && mes2 == _Km && mes3 == _Kz)
        return 1. / (12 * fpi * fpi) * sqrt(2) * del;

    if (mes1 == _pim && mes2 == _Kp && mes3 == _Kz)
        return 1. / (12 * fpi * fpi) * sqrt(2) * del;

    if (mes1 == _eta && mes2 == _Kp && mes3 == _Km)
        return 1. / (6 * sqrt(3) * fpi * fpi) * (-3 * sqrt(2) * mK * mK * sin(te) + (mpi * mpi - mKz * mKz) * (cos(te) - sqrt(2) * sin(te)));

    if (mes1 == _etp && mes2 == _Kp && mes3 == _Km)
        return 1. / (6 * sqrt(3) * fpi * fpi) * (3 * sqrt(2) * mK * mK * cos(te) + (mpi * mpi - mKz * mKz) * (sin(te) + sqrt(2) * cos(te)));

    if (mes1 == _eta && mes2 == _Kz && mes3 == _Kz)
        return 1. / (6 * sqrt(3) * fpi * fpi) * (3 * sqrt(2) * mKz * mKz * sin(te) - (mpi * mpi - mK * mK) * (cos(te) - sqrt(2) * sin(te)));

    if (mes1 == _etp && mes2 == _Kz && mes3 == _Kz)
        return 1. / (6 * sqrt(3) * fpi * fpi) * (-3 * sqrt(2) * mKz * mKz * cos(te) - (mpi * mpi + mK * mK) * (sin(te) + sqrt(2) * cos(te)));

    if (mes1 == _eta && mes2 == _eta && mes3 == _eta)
        return 3 * 2 * del / (18 * sqrt(3) * fpi * fpi) * pow(cos(te) - sqrt(2) * sin(te), 3);

    if (mes1 == _eta && mes2 == _eta && mes3 == _etp)
        return 2 * del / (18 * sqrt(3) * fpi * fpi) * 3 * pow(cos(te) - sqrt(2) * sin(te), 2) * pow(sin(te) + sqrt(2) * cos(te), 1);

    if (mes1 == _eta && mes2 == _etp && mes3 == _etp)
        return 2 * del / (18 * sqrt(3) * fpi * fpi) * 3 * pow(cos(te) - sqrt(2) * sin(te), 1) * pow(sin(te) + sqrt(2) * cos(te), 2);

    if (mes1 == _etp && mes2 == _etp && mes3 == _etp)
        return 3 * 2 * del / (18 * sqrt(3) * fpi * fpi) * pow(sin(te) + sqrt(2) * cos(te), 3);

    return 0;
}

double CPOddHiggsLowMass::AmpEtaf(int i, int mes1, int mes2, int mes3)
{
    Order(mes1, mes2, mes3);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _pi3)
    {
        //std::cout << "here " << AmpPi3f(i, _pi3, _pi3, _eta) << std::endl;
        return AmpPi3f(i, _pi3, _pi3, _eta);
    }

    if (mes1 == _pi3 && mes2 == _pip && mes3 == _pim)
        return AmpPi3f(i, _pip, _pim, _eta);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _eta)
        return AmpPi3f(i, _pi3, _eta, _eta);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _etp)
        return AmpPi3f(i, _pi3, _eta, _etp);

    if (mes1 == _pip && mes2 == _pim && mes3 == _eta)
        return AmpPi3f(i, _pi3, _eta, _eta);

    if (mes1 == _pip && mes2 == _pim && mes3 == _etp)
        return 2 * del6 * mpi * mpi / (12 * fpi * fpi) * 2 * (cos(te) - sqrt(2) * sin(te)) * (sin(te) + sqrt(2) * cos(te));

    if (mes1 == _pi3 && mes2 == _Kp && mes3 == _Km)
        return AmpPi3f(i, _eta, _Kp, _Km);

    if (mes1 == _pi3 && mes2 == _Kz && mes3 == _Kz)
        return AmpPi3f(i, _eta, _Kz, _Kz);

    if (mes1 == _pip && mes2 == _Km && mes3 == _Kz)
        return 1. / (6 * sqrt(6) * fpi * fpi) * ((2 * mpi * mpi - mK * mK - mKz * mKz) * cos(te) - 2 * sqrt(2) * (mKz * mKz + mK * mK + mpi * mpi) * sin(te));

    if (mes1 == _pim && mes2 == _Kp && mes3 == _Kz)
        return 1. / (6 * sqrt(6) * fpi * fpi) * ((2 * mpi * mpi - mK * mK - mKz * mKz) * cos(te) - 2 * sqrt(2) * (mKz * mKz + mK * mK + mpi * mpi) * sin(te));

    if (mes1 == _pi3 && mes2 == _eta && mes3 == _eta)
        return AmpPi3f(i, _eta, _eta, _eta);

    if (mes1 == _pi3 && mes2 == _eta && mes3 == _etp)
        return AmpPi3f(i, _eta, _eta, _etp);

    if (mes1 == _pi3 && mes2 == _etp && mes3 == _etp)
        return AmpPi3f(i, _eta, _etp, _etp);

    return 0;
}

double CPOddHiggsLowMass::AmpEtapf(int i, int mes1, int mes2, int mes3)
{
    Order(mes1, mes2, mes3);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _pi3)
        return AmpPi3f(i, _pi3, _pi3, _etp);

    if (mes1 == _pi3 && mes2 == _pip && mes3 == _pim)
        return AmpPi3f(i, _pip, _pim, _etp);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _eta)
        return AmpPi3f(i, _pi3, _eta, _etp);

    if (mes1 == _pi3 && mes2 == _pi3 && mes3 == _etp)
        return AmpPi3f(i, _pi3, _etp, _etp);

    if (mes1 == _pip && mes2 == _pim && mes3 == _eta)
        return AmpPi3f(i, _pi3, _pi3, _etp);

    if (mes1 == _pip && mes2 == _pim && mes3 == _etp)
        return 2 * 2 * del8 * mpi * mpi / (12 * fpi * fpi) * pow(sin(te) + sqrt(2) * cos(te), 2);

    if (mes1 == _pi3 && mes2 == _Kp && mes3 == _Km)
        return AmpPi3f(i, _etp, _Kp, _Km);

    if (mes1 == _pi3 && mes2 == _Kz && mes3 == _Kz)
        return AmpPi3f(i, _etp, _Kz, _Kz);

    if (mes1 == _pip && mes2 == _Km && mes3 == _Kz)
        return 1. / (6 * sqrt(6) * fpi * fpi) * ((2 * mpi * mpi - mK * mK - mKz * mKz) * sin(te) + 2 * sqrt(2) * (mKz * mKz + mK * mK + mpi * mpi) * cos(te));

    if (mes1 == _pim && mes2 == _Kp && mes3 == _Kz)
        return 1. / (6 * sqrt(6) * fpi * fpi) * ((2 * mpi * mpi - mK * mK - mKz * mKz) * sin(te) + 2 * sqrt(2) * (mKz * mKz + mK * mK + mpi * mpi) * cos(te));

    if (mes1 == _pi3 && mes2 == _eta && mes3 == _eta)
        return AmpPi3f(i, _eta, _eta, _etp);

    if (mes1 == _pi3 && mes2 == _eta && mes3 == _etp)
        return AmpPi3f(i, _eta, _etp, _etp);

    if (mes1 == _pi3 && mes2 == _etp && mes3 == _etp)
        return AmpPi3f(i, _etp, _etp, _etp);

    return 0;
}

double CPOddHiggsLowMass::Sfactor(int mes1, int mes2, int mes3)
{
    if ((mes1 == mes2 && mes1 != mes3) || (mes1 == mes3 && mes1 != mes2) || (mes2 == mes3 && mes2 != mes1))
        return 2;
    if (mes1 == mes2 && mes2 == mes3)
        return 6;

    return 1;
}

void CPOddHiggsLowMass::Order(int &i, int &j, int &k)
{
    int _i = i, _j = j, _k = k;

    if (_i != _j && _j != _k && _i != _k)
    {
        if (_k < _i && _i < _j)
        {
            i = _k;
            j = _i;
            k = _j;
        }
        if (_j < _k && _k < _i)
        {
            i = _j;
            j = _k;
            k = _i;
        }
        if (_j < _i && _i < _k)
        {
            i = _j;
            j = _i;
        }
        if (_k < _j && _j < _i)
        {
            i = _k;
            k = _i;
        }
        if (_i < _k && _k < _j)
        {
            j = _k;
            k = _j;
        }
    }
    if (_i == _j && _i != _k && _k < _i)
    {
        i = _k;
        k = _i;
    }
    if (_j == _k && _i != _j && _j < _i)
    {
        i = _k;
        k = _i;
    }
    if (_i == _k && _i != _j)
    {
        if (_j < _i)
        {
            i = _j;
            j = _k;
        }
        if (_i < _j)
        {
            j = _k;
            k = _j;
        }
    }
}

//
//
//
//

void CPOddHiggsLowMass::Read()
{
    std::string line, token;
    size_t pos0, pos1;

    lamb.resize(0);
    kappa.resize(0);
    Tbeta.resize(0);
    mH1.resize(0);
    mA1.resize(0);
    mch0.resize(0);
    mch1.resize(0);
    mch1.resize(0);
    GA1N1N1.resize(0);
    U11.resize(0);
    U12.resize(0);
    U21.resize(0);
    U22.resize(0);
    V11.resize(0);
    V12.resize(0);
    V21.resize(0);
    V22.resize(0);
    P11.resize(0);
    P12.resize(0);

    std::string delimiter = " ";

    int j = 0;
    bool end_reached = false;

    std::fstream infile(name_file);
    if (infile.is_open())
    {
        while (getline(infile, line))
        {
            pos1 = line.find("#");

            //Check that the line is not a comment
            if (pos1 == std::string::npos)
            {
                do
                {

                    pos0 = line.find(delimiter);
                    lamb.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    kappa.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    Tbeta.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    mH1.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    mA1.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    mch0.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    mch1.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    mch2.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    GA1N1N1.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    U11.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    U12.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    U21.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    U22.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    V11.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    V12.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    V21.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    V22.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    P11.push_back(stod(line.substr(0, pos0)));
                    //std::cout << line << " | " << stod(line.substr(0, pos0)) << //std::endl;
                    line.erase(0, pos0 + delimiter.length());

                    pos0 = line.find(delimiter);
                    P12.push_back(stod(line.substr(0, pos0)));

                    end_reached = true;

                } while (!end_reached);
            }

            //std::cout << lamb.size() - 1 << " " << lamb[lamb.size() - 1] << std::endl;
        }
    }

    Npoints = lamb.size();
    Cotb.resize(Npoints);
    Sbeta.resize(Npoints);
    Cbeta.resize(Npoints);
    YllA1e.resize(Npoints);
    YllA1mu.resize(Npoints);

    for (int i = 0; i < Npoints; i++)
    {
        Cotb[i] = 1. / Tbeta[i];
        Sbeta[i] = sqrt(Tbeta[i] * Tbeta[i] / (1 + Tbeta[i] * Tbeta[i]));
        Cbeta[i] = sqrt(1. / (Tbeta[i] * Tbeta[i] + 1));

        // Initialisation of the couplings A1 - e and mu
        YllA1e[i] = M_e / (sqrt(2) * vv) * P11[i] * Tbeta[i];
        YllA1mu[i] = M_MU / (sqrt(2) * vv) * P11[i] * Tbeta[i];
    }
}

//
//

void CPOddHiggsLowMass::Write(std::string const &name_file_out)
{
    std::ofstream out_file;
    out_file.open("../input/" + name_file_out);

    time_t now = time(0);

    // convert now to string form
    char *dt = ctime(&now);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);

    out_file << "####################################################" << std::endl;
    out_file << "#" << std::endl;
    out_file << "# CrossSection input file" << std::endl;
    out_file << "# Automatically generated from CPOddHiggsLowMass class" << std::endl;
    out_file << "#" << std::endl;
    out_file << "# Gaetan Facchinetti" << std::endl;
    out_file << "# gaetan.facchinetti@umontpellier.fr" << std::endl;
    out_file << "# Laboratoire Univers et Particules Montpellier" << std::endl;
    out_file << "# Created on : (UTC) " << dt;
    out_file << "#" << std::endl;
    out_file << "####################################################" << std::endl;
    out_file << std::endl;

    out_file << "## Number of points to evaluate ##" << std::endl;
    out_file << "Npoints :: " << lamb.size() << std::endl;
    out_file << std::endl;
    out_file << std::endl;

    out_file << "## Dark sector particles ##" << std::endl;
    out_file << "DSCONT :: NDM :: 1" << std::endl;
    out_file << "DSCONT :: DSPart :: (pseudoscalar, 1)" << std::endl;
    out_file << std::endl;
    out_file << std::endl;

    out_file << "## Masses dark sector [GeV] ##" << std::endl;

    for (int i = 0; i < lamb.size(); i++)
        out_file << "MDMDS :: " << i << " :: " << mch0[i] << ", " << mA1[i] << std::endl;

    out_file << std::endl;
    out_file << std::endl;

    out_file << "## Interactions DS-SM [GeV] ##" << std::endl;

    for (int i = 0; i < lamb.size(); i++)
        out_file << "IntDSSM :: " << i << " :: 0 :: "
                 << "l_e=" << P11[i] * 2 * V * M_e * Tbeta[i]
                 << ", l_mu=" << P11[i] * 2 * V * M_MU * Tbeta[i]
                 << ", l_d=" << P11[i] * 2 * V * M_QUARK_DOWN * Tbeta[i]
                 << ", l_u=" << P11[i] * 2 * V * M_QUARK_UP * (1. / Tbeta[i])
                 << ", l_s=" << P11[i] * 2 * V * M_QUARK_STRANGE * Tbeta[i] << std::endl;

    out_file << std::endl;

    out_file << "## Interactions DS-DM [GeV] ##" << std::endl;

     for (int i = 0; i < lamb.size(); i++)
        out_file << "IntDSDM :: " << i << " :: 0 :: "
                 << "l_0_0=" << GA1N1N1[i] << std::endl;

}

//
//
//
//

void CPOddHiggsLowMass::plot_MixingParameters()
{
    std::ofstream outfile;
    outfile.open("../output/CPOddHiggsLowMass/MixingParameters.out");

    //std::cout << "../output/MixingParameters_" + name_file + ".out" << std::endl;

    outfile << "# From file : " << name_file << std::endl;
    outfile << "# mA1 [GeV] | abs(OA3) | abs(OAeta) | abs(OAetp) | abs(OAA) | abs(OA8) | abs(OA9) " << std::endl;
    outfile << "#" << std::endl;

    for (int i = 0; i < Npoints; i++)
    {
        //std::cout << mA1[i] << std::endl;
        outfile << mA1[i] << "\t" << abs(OA3[i]) << "\t" << abs(OAeta[i]) << "\t" << abs(OAetp[i]) << "\t"
                << abs(OAA[i]) << "\t" << abs(OA8[i]) << "\t" << abs(OA9[i]) << std::endl;
    }

    outfile.close();
}

void CPOddHiggsLowMass::plot_GammaA()
{

    std::ofstream outfile;
    outfile.open("../output/GammaA.out");

    //std::cout << "../output/MixingParameters_" + name_file + ".out" << std::endl;

    outfile << "# From file : " << name_file << std::endl;
    outfile << "# mA1 [GeV] | GamamTAS [GeV] | GamaL [GeV] | GamaAch [GeV] | GamaAt1S [GeV] | Gamach0 [GeV]" << std::endl;
    outfile << std::endl;

    for (int i = 0; i < Npoints; i++)
        outfile << mA1[i] << "\t" << GammaTAS(i, 0.1) << "\t" << GamaL(i, 0.1) << "\t" << GamaAch(i, 0.1) << "\t" << GamaAt1S(i, 0.1) << "\t" << Gamach0(i, 0.1) << std::endl;

    outfile.close();
}

void CPOddHiggsLowMass::plot_GamaA(int mes1, int mes2, int mes3)
{

    std::ofstream outfile;
    outfile.open("../output/GamaA.out");

    //std::cout << "../output/MixingParameters_" + name_file + ".out" << std::endl;

    outfile << "# From file : " << name_file << std::endl;
    outfile << "# mA1 [GeV] | GamaA [GeV]" << std::endl;
    outfile << std::endl;

    for (int i = 0; i < Npoints; i++)
        outfile << mA1[i] << "\t" << GamaA(i, mes2, mes2, mes3) << std::endl;

    outfile.close();
}

void CPOddHiggsLowMass::plot_Amp()
{

    std::ofstream outfile;
    outfile.open("../output/Amp.out");

    //std::cout << "../output/MixingParameters_" + name_file + ".out" << std::endl;

    outfile << "# From file : " << name_file << std::endl;
    outfile << "# mA1 [GeV] | GamaA [GeV]" << std::endl;
    outfile << std::endl;

    for (int mes1 = 0; mes1 < 8; mes1++)
        for (int mes2 = mes1; mes2 < 8; mes2++)
            for (int mes3 = mes2; mes3 < 8; mes3++)
                outfile << mes1 << "\t" << mes2 << "\t" << mes3 << "\t" << AmpA1f(0, mes1, mes2, mes3) << std::endl;

    outfile.close();
}

void CPOddHiggsLowMass::plot_sFg(int i)
{

    std::ostringstream ival;
    ival << i;
    std::string ival_str = ival.str();

    std::ofstream outfile;
    outfile.open("../output/sFg_" + ival_str + ".out");

    //std::cout << "../output/MixingParameters_" + name_file + ".out" << std::endl;

    outfile << "# From file : " << name_file << std::endl;
    outfile << "# s [GeV^2] | sFg(s) | mA1 [GeV]" << std::endl;
    outfile << std::endl;

    int Npts = 500;

    double s, smin = 1e-2 * mA1[i] * mA1[i], smax = 100 * mA1[i] * mA1[i];
    double d_logs = log10(smax / smin) / Npts;

    for (int j = 0; j < Npts; j++)
    {
        s = smin * pow(10, j * d_logs);
        outfile << s << "\t" << abs(sFg(i,  s, 0.1)) << "\t" << mA1[i] << std::endl;
        //std::cout << "here" << std::endl;
    }

    outfile.close();
}


void CPOddHiggsLowMass::plot_sCgam(int i)
{

    std::ostringstream ival;
    ival << i;
    std::string ival_str = ival.str();

    std::ofstream outfile;
    outfile.open("../output/sCgam_" + ival_str + ".out");

    //std::cout << "../output/MixingParameters_" + name_file + ".out" << std::endl;

    outfile << "# From file : " << name_file << std::endl;
    outfile << "# s [GeV^2] | sCgam(s) | mA1 [GeV]" << std::endl;
    outfile << std::endl;

    int Npts = 500;

    double s, smin = 1e-2 * mA1[i] * mA1[i], smax = 100 * mA1[i] * mA1[i];
    double d_logs = log10(smax / smin) / Npts;

    for (int j = 0; j < Npts; j++)
    {
        s = smin * pow(10, j * d_logs);
        outfile << s << "\t" << abs(sdCgam(i,  s)) << "\t" << mA1[i] << std::endl;
        //std::cout << "here" << std::endl;
    }

    outfile.close();
}

