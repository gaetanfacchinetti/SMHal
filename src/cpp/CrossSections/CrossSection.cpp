#include "../../headers/CrossSections/CrossSection.h"

// !!!!!!
// The 26/04/2019 I changed L1 for LP in the evaluation
// of the cross section for ST and SU interference terms
//
// The 29/03/21 I combined CrossSectionScatt and CrossSection together

// ATTENTION : in this version the particles are passed as references and not as pointers
// therefore if in the code their properties are changed, a new cross section object must
// be defined in order to take these changes into account.

CrossSection::CrossSection(Particle const &p_1, Particle const &p_2, Particle const &p_3, Particle const &p_4)
{

    _m1 = p_1.get_mass();
    _m2 = p_2.get_mass();
    _m3 = p_3.get_mass();
    _m4 = p_4.get_mass();

    _g1 = p_1.get_degFree();
    _g2 = p_2.get_degFree();

    _p_1 = p_1;
    _p_2 = p_2;
    _p_3 = p_3;
    _p_4 = p_4;

    // ------ Symmetry factor for outgoing particles -------- //
    // WARNING : Be carefull that with this definition in the code the exact same Particle objects
    // ... have to be passes to the constructor for identical outgoing particles
    if (&p_3 == &p_4)
        _sym = 2.;
    else
        _sym = 1.;

    //std::cout << "The symmetry factor is " << _sym << std::endl;

    mandSum = _m1 * _m1 + _m2 * _m2 + _m3 * _m3 + _m4 * _m4;

    ns = 0;
    nt = 0;
    nu = 0;

    ParticlePairsSSSq.resize(0);
    ParticlePairsSSInt.resize(0);
    ParticlePairsUUSq.resize(0);
    ParticlePairsTTSq.resize(0);
    ParticlePairsTTInt.resize(0);
    ParticlePairsUUInt.resize(0);
    ParticlePairsSUInt.resize(0);
    ParticlePairsSTInt.resize(0);
    ParticlePairsTUInt.resize(0);

    speInt = SpecialIntegrations();
}

void CrossSection::Initialise(std::vector<Particle> &as, std::vector<Particle> &at, std::vector<Particle> &au, std::string name)
{
    ns = as.size();
    nt = at.size();
    nu = au.size();

    ParticlePairsSSSq.resize(ns);
    ParticlePairsSSInt.resize((int)(ns * (ns - 1.) / 2.));
    ParticlePairsTTSq.resize(nt);
    ParticlePairsTTInt.resize((int)(nt * (nt - 1.) / 2.));
    ParticlePairsUUSq.resize(nu);
    ParticlePairsUUInt.resize((int)(nu * (nu - 1.) / 2.));
    ParticlePairsSUInt.resize(ns * nu);
    ParticlePairsSTInt.resize(ns * nt);
    ParticlePairsTUInt.resize(nt * nu);

    int compt = 0;
    int pairIndex = 0;

    if (ParticlePairsSSSq.size() != 0)
        for (int i = 0; i < ns; i++)
        {
            ParticlePairsSSSq[i] = ParticlePair(as[i], as[i], pairIndex, INT_TYPE::SS);
            pairIndex++;
        }

    if (ParticlePairsSSInt.size() != 0)
        for (int i = 0; i < ns; i++)
            for (int j = i + 1; j < ns; j++)
            {
                ParticlePairsSSInt[compt] = ParticlePair(as[i], as[j], pairIndex, INT_TYPE::SS);
                compt++;
                pairIndex++;
            }

    compt = 0;

    if (ParticlePairsTTSq.size() != 0)
        for (int i = 0; i < nt; i++)
        {
            ParticlePairsTTSq[i] = ParticlePair(at[i], at[i], pairIndex, INT_TYPE::TT);
            pairIndex++;
        }

    if (ParticlePairsTTInt.size() != 0)
        for (int i = 0; i < nt; i++)
            for (int j = i + 1; j < nt; j++)
            {
                ParticlePairsTTInt[compt] = ParticlePair(at[i], at[j], pairIndex, INT_TYPE::TT);
                compt++;
                pairIndex++;
            }

    compt = 0;

    if (ParticlePairsUUSq.size() != 0)
        for (int i = 0; i < nu; i++)
        {
            ParticlePairsUUSq[i] = ParticlePair(au[i], au[i], pairIndex, INT_TYPE::UU);
            pairIndex++;
        }

    if (ParticlePairsUUInt.size() != 0)
        for (int i = 0; i < nu; i++)
            for (int j = i + 1; j < nu; j++)
            {
                ParticlePairsUUInt[compt] = ParticlePair(au[i], au[j], pairIndex, INT_TYPE::UU);
                compt++;
                pairIndex++;
            }

    compt = 0;

    if (ParticlePairsSTInt.size() != 0)
        for (int i = 0; i < ns; i++)
            for (int j = 0; j < nt; j++)
            {
                ParticlePairsSTInt[compt] = ParticlePair(as[i], at[j], pairIndex, INT_TYPE::ST);
                compt++;
                pairIndex++;
            }

    compt = 0;

    if (ParticlePairsSUInt.size() != 0)
        for (int i = 0; i < ns; i++)
            for (int j = 0; j < nu; j++)
            {
                ParticlePairsSUInt[compt] = ParticlePair(as[i], au[j], pairIndex, INT_TYPE::SU);
                compt++;
                pairIndex++;
            }

    compt = 0;

    if (ParticlePairsTUInt.size() != 0)
        for (int i = 0; i < nt; i++)
            for (int j = 0; j < nu; j++)
            {
                ParticlePairsTUInt[compt] = ParticlePair(at[i], au[j], pairIndex, INT_TYPE::TU);
                compt++;
                pairIndex++;
            }

    //std::cout << "aa " << Verbose << std::endl;

    if (Verbose == true)
    {
        std::clog << std::endl;
        std::clog << "## Particle pairs - CrossSection : " << name << std::endl;
        std::clog << "SS Squared : " << ParticlePairsSSSq << std::endl;
        std::clog << "SS Interf. : " << ParticlePairsSSInt << std::endl;
        std::clog << "TT Squared : " << ParticlePairsTTSq << std::endl;
        std::clog << "TT Interf. : " << ParticlePairsTTInt << std::endl;
        std::clog << "UU Squared : " << ParticlePairsUUSq << std::endl;
        std::clog << "UU Interf. : " << ParticlePairsUUInt << std::endl;
        std::clog << "ST Interf. : " << ParticlePairsSTInt << std::endl;
        std::clog << "SU Interf. : " << ParticlePairsSUInt << std::endl;
        std::clog << "TU Interf. : " << ParticlePairsTUInt << std::endl;
        std::cout << "--> CrossSection " << name << " initialised : see log file" << std::endl;
    }

    //coeffSE.resize(pairIndex - 1);
    //InitializeCoeffSE_TTUU(ParticlePairsTTSq);
    //InitializeCoeffSE_TTUU(ParticlePairsTTInt);
    //InitializeCoeffSE_TTUU(ParticlePairsUUSq);
    //InitializeCoeffSE_TTUU(ParticlePairsUUInt);
    // std::cout << (this->*polyQ[0])(0, 4 * 100 * 100 + 1) << std::endl;
}

// =============
// Classical cross section
double CrossSection::EvalSigmaP1cm2(double s)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2) || s < pow(_m1 + _m2, 2))
        return 0;

    //std::cout << "ooooooaaaa " << std::endl;
    //std::cout << ParticlePairsSSSq[0] << std::endl;

    for (int i = 0; i < ParticlePairsSSSq.size(); i++)
    {
        m0 = ParticlePairsSSSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSSq[i].get_particle(0).get_width();

        res += 1. / (pow(s - m0 * m0, 2) + m0 * m0 * w0 * w0) *
               real(IntL1(ParticlePairsSSSq[i], s));
    }

    //std::cout << "oooooo " << std::endl;
    //std::cout << "res here " << res << std::endl;

    for (int i = 0; i < ParticlePairsTTSq.size(); i++)
        res += real(IntLPtPt(ParticlePairsTTSq[i], s));

    //std::cout << "res here " << res << std::endl;

    for (int i = 0; i < ParticlePairsUUSq.size(); i++)
        res += real(IntLPuPu(ParticlePairsUUSq[i], s));

    for (int i = 0; i < ParticlePairsSSInt.size(); i++)
    {
        m0 = ParticlePairsSSInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSInt[i].get_particle(0).get_width();

        m1 = ParticlePairsSSInt[i].get_particle(1).get_mass();
        w1 = ParticlePairsSSInt[i].get_particle(1).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0) * ((s - m1 * m1) + one_i * m1 * w1)) * IntL1(ParticlePairsSSInt[i], s));
    }

    for (int i = 0; i < ParticlePairsTTInt.size(); i++)
        res += 2 * real(IntLPtPt(ParticlePairsTTInt[i], s));

    for (int i = 0; i < ParticlePairsUUInt.size(); i++)
        res += 2 * real(IntLPuPu(ParticlePairsUUInt[i], s));

    for (int i = 0; i < ParticlePairsSTInt.size(); i++)
    {
        m0 = ParticlePairsSTInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSTInt[i].get_particle(0).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntLP(ParticlePairsSTInt[i], s));
    }

    for (int i = 0; i < ParticlePairsSUInt.size(); i++)
    {
        m0 = ParticlePairsSUInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSUInt[i].get_particle(0).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntLP(ParticlePairsSUInt[i], s));
    }

    for (int i = 0; i < ParticlePairsTUInt.size(); i++)
        res += 2 * real(IntLPtPu(ParticlePairsTUInt[i], s));

    // double p1cm = (1. / (2. * sqrt(s))) * sqrt((s - pow(_m1 - _m2, 2)) * (s - pow(_m1 + _m2, 2)));
    //std::cout << "res here " << res << std::endl;
    // res /= (64 * PI * s * p1cm * p1cm * _g1 * _g2);
    //std::cout << (64 * PI * s * f_p1cm(s) * f_p1cm(s) * _g1 * _g2) << " " << f_p1cm(s) << " " << res << " " << _g1  << " " << _g2 << std::endl;
    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::EvalSigmaP1cm2_fromIntegral(double s)
{
    double res  = 0;
    if (s < pow(_m3 + _m4, 2) || s < pow(_m1 + _m2, 2))
        return 0;

    std::vector<double> ts = {0, s};
    res =  GaussLegendre_IntegralLn_Static(0, 2, gl100, tmin(s), tmax(s), ts, this, CallBack_fToIntOnt_for_EvalSigmaP1cm2_fromIntegral);
    res /= (64 * PI * s * _g1 * _g2) * _sym;
    
    return res;
}

double CrossSection::fToIntOnt_for_EvalSigmaP1cm2_fromIntegral(double t, double s)
{
    double res = 0;
    double m0, w0, m1, w1;

    for (int i = 0; i < ParticlePairsSSSq.size(); i++)
    {
        m0 = ParticlePairsSSSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSSq[i].get_particle(0).get_width();

        for (int n = 0; n < Q_order(ParticlePairsSSSq[i]) + 1; n++)
            res += real(Q(ParticlePairsSSSq[i], n, s, t)) * pow(t, n) / (pow(s - m0 * m0, 2) + m0 * m0 * w0 * w0);
    }

    for (int i = 0; i < ParticlePairsTTSq.size(); i++)
    {
        m0 = ParticlePairsTTSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsTTSq[i].get_particle(0).get_width();
        
        for (int n = 0; n < Q_order(ParticlePairsTTSq[i]) + 1; n++)
            res += real(Q(ParticlePairsTTSq[i], n, s, t)) * pow(t, n) / (pow(t- m0 * m0, 2) + m0 * m0 * w0 * w0 );
    }

    double u = 0;
    for (int i = 0; i < ParticlePairsUUSq.size(); i++)
    {
        m0 = ParticlePairsUUSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsUUSq[i].get_particle(0).get_width();
        
        u = _m1*_m1 + _m2*_m2 + _m3*_m3 + _m4*_m4 - s -t;

        for (int n = 0; n < Q_order(ParticlePairsUUSq[i]) + 1; n++)
            res += real(Q(ParticlePairsUUSq[i], n, s, t)) * pow(t, n) / (pow(u - m0 * m0, 2) + m0 * m0 * w0 * w0 );
    }

    if(ParticlePairsSSInt.size() > 0 || ParticlePairsTTInt.size() > 0 || ParticlePairsUUInt.size() > 0 ||
       ParticlePairsSTInt.size() > 0 || ParticlePairsTUInt.size() > 0 || ParticlePairsSUInt.size() > 0)
       {
           std::cout << "ERROR : Cannot use " << __PRETTY_FUNCTION__ << " in these configurations yet" << std::endl;
           exit(0);
       }

    return res;
}



// Transfer CrossSection
double CrossSection::EvalSigmaTransferP1cm4(double s)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2) || s < pow(_m1 + _m2, 2))
        return 0;

    for (int i = 0; i < ParticlePairsSSSq.size(); i++)
    {
        m0 = ParticlePairsSSSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSSq[i].get_particle(0).get_width();

        res += 1. / (pow(s - m0 * m0, 2) + m0 * m0 * w0 * w0) * real(IntR1(ParticlePairsSSSq[i], s));
    }

    for (int i = 0; i < ParticlePairsTTSq.size(); i++)
        res += real(IntRPtPt(ParticlePairsTTSq[i], s));

    for (int i = 0; i < ParticlePairsUUSq.size(); i++)
        res += real(IntRPuPu(ParticlePairsUUSq[i], s));

    for (int i = 0; i < ParticlePairsSSInt.size(); i++)
    {
        m0 = ParticlePairsSSInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSInt[i].get_particle(0).get_width();

        m1 = ParticlePairsSSInt[i].get_particle(1).get_mass();
        w1 = ParticlePairsSSInt[i].get_particle(1).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0) * ((s - m1 * m1) + one_i * m1 * w1)) * IntR1(ParticlePairsSSInt[i], s));
    }

    for (int i = 0; i < ParticlePairsTTInt.size(); i++)
        res += 2 * real(IntRPtPt(ParticlePairsTTInt[i], s));

    for (int i = 0; i < ParticlePairsUUInt.size(); i++)
        res += 2 * real(IntRPuPu(ParticlePairsUUInt[i], s));

    for (int i = 0; i < ParticlePairsSTInt.size(); i++)
    {
        m0 = ParticlePairsSTInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSTInt[i].get_particle(0).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntRP(ParticlePairsSTInt[i], s));
    }

    for (int i = 0; i < ParticlePairsSUInt.size(); i++)
    {
        m0 = ParticlePairsSUInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSUInt[i].get_particle(0).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntRP(ParticlePairsSUInt[i], s));
    }

    for (int i = 0; i < ParticlePairsTUInt.size(); i++)
        res += 2 * real(IntRPtPu(ParticlePairsTUInt[i], s));

    //std::cout << (16 * PI * s * pow(f_p1cm(s), 4) * 8 * _g1 * _g2) << " " << f_p1cm(s) << " " << res << " " << PI << " " << _g1  << " " << _g2 << std::endl;
    return (-res) / (16 * PI * s * _g1 * _g2 * 8 * _sym);
    //return (-res);
}

// Transfer CrossSection according to Bringmann
double CrossSection::EvalSigmaTransferBringmannP1cm4(double omega)
{
    double res = 0;

    double s = pow(_m1, 2) + 2 * omega * _m1 + pow(_m2, 2);
    std::cout.precision(8);

    if (s < pow(_m3 + _m4, 2) || s < pow(_m1 + _m2, 2))
        return 0;

    double m0, w0, m1, w1;

    for (int i = 0; i < ParticlePairsSSSq.size(); i++)
    {
        m0 = ParticlePairsSSSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSSq[i].get_particle(0).get_width();

        res += 1. / (pow(s - m0 * m0, 2) + m0 * m0 * w0 * w0) * real(IntR1(ParticlePairsSSSq[i], s));
    }

    for (int i = 0; i < ParticlePairsTTSq.size(); i++)
        res += real(IntRPtPt(ParticlePairsTTSq[i], s));

    for (int i = 0; i < ParticlePairsUUSq.size(); i++)
        res += real(IntRPuPu(ParticlePairsUUSq[i], s));

    for (int i = 0; i < ParticlePairsSSInt.size(); i++)
    {
        m0 = ParticlePairsSSInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSInt[i].get_particle(0).get_width();

        m1 = ParticlePairsSSInt[i].get_particle(1).get_mass();
        w1 = ParticlePairsSSInt[i].get_particle(1).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0) * ((s - m1 * m1) + one_i * m1 * w1)) * IntR1(ParticlePairsSSInt[i], s));
    }

    for (int i = 0; i < ParticlePairsTTInt.size(); i++)
        res += 2 * real(IntRPtPt(ParticlePairsTTInt[i], s));

    for (int i = 0; i < ParticlePairsUUInt.size(); i++)
        res += 2 * real(IntRPuPu(ParticlePairsUUInt[i], s));

    for (int i = 0; i < ParticlePairsSTInt.size(); i++)
    {
        m0 = ParticlePairsSTInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSTInt[i].get_particle(0).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntRP(ParticlePairsSTInt[i], s));
    }

    for (int i = 0; i < ParticlePairsSUInt.size(); i++)
    {
        m0 = ParticlePairsSUInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSUInt[i].get_particle(0).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntRP(ParticlePairsSUInt[i], s));
    }

    for (int i = 0; i < ParticlePairsTUInt.size(); i++)
        res += 2 * real(IntRPtPu(ParticlePairsTUInt[i], s));

    //std::cout << s << " " << omega << " " << res << " " << tmin(s) << " " << _m1 << " " << _m2 << std::endl;

    return (-res) / (16 * PI * _m1 * _m1 * 8 * _sym);
    // BE CAREFULL : here we do not want to devide by the number of degree of freedom of p_1 and p_2
    // because this is included in the definition of gamma (the scattering rate) ..
    // This implementation should be modified to be more consistent
    // ATTENTION : needs to impose _m1 as the DM candidate !!!
}

double CrossSection::EvalSigmaTransferBringmannP1cm4_fromIntegral(double omega)
{
    double res = 0;
    double s = pow(_m1, 2) + 2 * omega * _m1 + pow(_m2, 2);
    std::cout.precision(8);

    if (s < pow(_m3 + _m4, 2) || s < pow(_m1 + _m2, 2))
        return 0;

    std::vector<double> ts = {0, s};
    res =  GaussLegendre_IntegralLn_Static(0, 2, gl100, tmin(s), tmax(s), ts, this, CallBack_fToIntOnt_for_EvalSigmaTransferBringmannP1cm4_fromIntegral);
    
    return (-res) / (16 * PI * _m1 * _m1 * 8 * _sym);
}

double CrossSection::fToIntOnt_for_EvalSigmaTransferBringmannP1cm4_fromIntegral(double t, double s)
{
    double res = 0;
    double m0, w0, m1, w1;

    for (int i = 0; i < ParticlePairsSSSq.size(); i++)
    {
        m0 = ParticlePairsSSSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSSq[i].get_particle(0).get_width();

        for (int n = 0; n < Q_order(ParticlePairsSSSq[i]) + 1; n++)
            res += real(Q(ParticlePairsSSSq[i], n, s, t)) * pow(t, n+1) / (pow(s - m0 * m0, 2) + m0 * m0 * w0 * w0);
    }

    for (int i = 0; i < ParticlePairsTTSq.size(); i++)
    {
        m0 = ParticlePairsTTSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsTTSq[i].get_particle(0).get_width();
        
        for (int n = 0; n < Q_order(ParticlePairsTTSq[i]) + 1; n++)
            res += real(Q(ParticlePairsTTSq[i], n, s, t)) * pow(t, n+1) / (pow(t- m0 * m0, 2) + m0 * m0 * w0 * w0 );
    }

    double u = 0;
    for (int i = 0; i < ParticlePairsUUSq.size(); i++)
    {
        m0 = ParticlePairsUUSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsUUSq[i].get_particle(0).get_width();
        
        u = _m1*_m1 + _m2*_m2 + _m3*_m3 + _m4*_m4 - s -t;

        for (int n = 0; n < Q_order(ParticlePairsUUSq[i]) + 1; n++)
            res += real(Q(ParticlePairsUUSq[i], n, s, t)) * pow(t, n+1) / (pow(u - m0 * m0, 2) + m0 * m0 * w0 * w0 );
    }

    if(ParticlePairsSSInt.size() > 0 || ParticlePairsTTInt.size() > 0 || ParticlePairsUUInt.size() > 0 ||
       ParticlePairsSTInt.size() > 0 || ParticlePairsTUInt.size() > 0 || ParticlePairsSUInt.size() > 0)
       {
           std::cout << "ERROR : Cannot use " << __PRETTY_FUNCTION__ << " in these configurations yet" << std::endl;
           exit(0);
       }

    return res;
}

//
//
//
//
//
//
//
//
// =======================
// Intgrals (functionals)
// =======================

dcomp CrossSection::IntLPtPt(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
        res += Q(partPair, i, s) * speInt.IntJttComp(i, 1, 1, m0 * m0, w0 * m0, m1 * m1, w1 * m1, tmin(s), tmax(s));

    //std::cout << m0 << " " << w0 << " " << m1 << " " << w1 << " " << tmin(s) << " " << tmax(s) << " " << f_p1cm(s) << std::endl;
    //res = speInt.IntJttComp(0, 1, 1, m0 * m0, w0 * m0, m1 * m1, w1 * m1, tmin(s), 0);

    return res;
}

dcomp CrossSection::IntLPtPu(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    double r1 = m1 * m1 + s - mandSum;

    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
        res += Q(partPair, i, s) * speInt.IntJtuComp(i, 1, 1, m0 * m0, w0 * m0, r1, w1 * m1, tmin(s), tmax(s));

    return res;
}

dcomp CrossSection::IntLPuPu(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    double r0 = m0 * m0 + s - mandSum;
    double r1 = m1 * m1 + s - mandSum;

    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
        res += Q(partPair, i, s) * speInt.IntJuuComp(i, 1, 1, r0, w0 * m0, r1, w1 * m1, tmin(s), tmax(s));

    return res;
}

dcomp CrossSection::IntLP(ParticlePair const &partPair, double s)
{

    double m = partPair.get_second().get_mass();
    double w = partPair.get_second().get_width();

    INT_TYPE type = partPair.get_type();

    //std::cout << partPair << std::endl;

    if (type != INT_TYPE::ST && type != INT_TYPE::SU)
    {
        std::cout << "FATAL ERROR : CrossSection::IntLP should not be used in this context " << std::endl;
        //std::cout << partPair.ind() << std::endl;
        exit(0);
    }

    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
    {
        if (type == INT_TYPE::ST)
            res += Q(partPair, i, s) * speInt.IntJtComp(i, 1, m * m, w * m, tmin(s), tmax(s));
        if (type == INT_TYPE::SU)
        {
            double r = m * m + s - mandSum;
            res += Q(partPair, i, s) * speInt.IntJuComp(i, 1, r, w * m, tmin(s), tmax(s));
        }
    }

    return res;
}

dcomp CrossSection::IntL1(ParticlePair const &partPair, double s)
{
    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
    {
        res += Q(partPair, i, s) * speInt.IntK(i, tmin(s), tmax(s));
    }

    //std::cout << "on a dans : " << speInt.IntK(1, tmin(s), tmax(s)) << " " << Q(partPair, 1, s)  << " " << tmin(s) << " " <<  tmax(s) << std::endl;
    //exit(0);

    //std::cout << "res ici : " << 2*real(res) << std::endl;

    return res;
}

//
//
// Function for the evaluation of the s-wave term
//

double CrossSection::EvalSWaveTerm()
{
    double res = 0;

    double m0, w0, m1, w1;

    if (pow(_m1 + _m2, 2) <= pow(_m3 + _m4, 2))
        return 0;

    // s is fixed at its minimal value
    double s = pow(_m1 + _m2, 2);

    for (int i = 0; i < ParticlePairsSSSq.size(); i++)
    {
        m0 = ParticlePairsSSSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSSq[i].get_particle(0).get_width();

        res += 1. / (pow(s - m0 * m0, 2) + m0 * m0 * w0 * w0) *
               real(IntLH1(ParticlePairsSSSq[i], s));
    }

    for (int i = 0; i < ParticlePairsTTSq.size(); i++)
        res += real(IntLHPtPt(ParticlePairsTTSq[i], s));

    //std::cout << "houra 1 " << res  << " " <<  real(IntLHPtPt(*ParticlePairsTTSq[0], s)) << std::endl;

    for (int i = 0; i < ParticlePairsUUSq.size(); i++)
        res += real(IntLHPuPu(ParticlePairsUUSq[i], s));

    //std::cout << "houra 2 " << res << " " << real(IntLHPtPt(*ParticlePairsUUSq[0], s)) << std::endl;

    for (int i = 0; i < ParticlePairsSSInt.size(); i++)
    {
        m0 = ParticlePairsSSInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSInt[i].get_particle(0).get_width();

        m1 = ParticlePairsSSInt[i].get_particle(1).get_mass();
        w1 = ParticlePairsSSInt[i].get_particle(1).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0) * ((s - m1 * m1) + one_i * m1 * w1)) * IntLH1(ParticlePairsSSInt[i], s));
    }

    for (int i = 0; i < ParticlePairsTTInt.size(); i++)
        res += 2 * real(IntLHPtPt(ParticlePairsTTInt[i], s));

    for (int i = 0; i < ParticlePairsUUInt.size(); i++)
        res += 2 * real(IntLHPuPu(ParticlePairsUUInt[i], s));

    for (int i = 0; i < ParticlePairsSTInt.size(); i++)
    {
        m0 = ParticlePairsSTInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSTInt[i].get_particle(0).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntLHP(ParticlePairsSTInt[i], s));
    }

    for (int i = 0; i < ParticlePairsSUInt.size(); i++)
    {
        m0 = ParticlePairsSUInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSUInt[i].get_particle(0).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntLHP(ParticlePairsSUInt[i], s));
    }

    for (int i = 0; i < ParticlePairsTUInt.size(); i++)
        res += 2 * real(IntLHPtPu(ParticlePairsTUInt[i], s));

    //std::cout << "houra 3 " << res << " " << 2 * real(IntLHPtPu(*ParticlePairsTUInt[0], s)) << std::endl;

    // double p1cm = (1. / (2. * sqrt(s))) * sqrt((s - pow(_m1 - _m2, 2)) * (s - pow(_m1 + _m2, 2)));

    // res /= (64 * PI * s * p1cm * p1cm * _g1 * _g2);
    res /= (64 * PI * sqrt(s) * _g1 * _g2) * _sym * _m1 * _m2;
    //    std::cout << s << " " << p1cm  << " " << res << std::endl;
    return res;
}

dcomp CrossSection::IntLHPtPt(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
        res += Q(partPair, i, s) * speInt.IntHttComp(i, 1, 1, m0 * m0, w0 * m0, m1 * m1, w1 * m1, tmax(s), f_p3cm(s));

    return res;
}

dcomp CrossSection::IntLHPtPu(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    double r1 = m1 * m1 + s - mandSum;

    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
        res += Q(partPair, i, s) * speInt.IntHtuComp(i, 1, 1, m0 * m0, w0 * m0, r1, w1 * m1, tmax(s), f_p3cm(s));

    return res;
}

dcomp CrossSection::IntLHPuPu(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    double r0 = m0 * m0 + s - mandSum;
    double r1 = m1 * m1 + s - mandSum;

    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
        res += Q(partPair, i, s) * speInt.IntHuuComp(i, 1, 1, r0, w0 * m0, r1, w1 * m1, tmax(s), f_p3cm(s));

    return res;
}

dcomp CrossSection::IntLHP(ParticlePair const &partPair, double s)
{

    double m = partPair.get_second().get_mass();
    double w = partPair.get_second().get_width();

    INT_TYPE type = partPair.get_type();

    if (type != INT_TYPE::ST && type != INT_TYPE::SU)
    {
        std::cout << "FATAL ERROR : CrossSection::IntLP should not be used in this context " << std::endl;
        exit(0);
    }

    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
    {
        if (type == INT_TYPE::ST)
            res += Q(partPair, i, s) * speInt.IntHtComp(i, 1, m * m, w * m, tmax(s), f_p3cm(s));
        if (type == INT_TYPE::SU)
        {
            double r = m * m + s - mandSum;
            res += Q(partPair, i, s) * speInt.IntHuComp(i, 1, r, w * m, tmax(s), f_p3cm(s));
        }
    }

    return res;
}

dcomp CrossSection::IntLH1(ParticlePair const &partPair, double s)
{
    dcomp res = 0;

    for (int i = 0; i < Q_order(partPair) + 1; i++)
    {
        res += Q(partPair, i, s) * speInt.IntHtComp(i, 0, 0, 0, tmax(s), f_p3cm(s));
        //std::cout << "Ici : " << partPair.get_first().get_index()  << "|" << partPair.get_second().get_index()  << " " << i << " " << s << " " << Q(partPair, i, s) << std::endl;
    }

    return res;
}

//
//
// Integrals
//
//

dcomp CrossSection::IntRPtPt(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    //std::cout << "Value of the width :"  << w0 << " " << w1 << " " << m0 << " " << m1 << std::endl;

    dcomp res = 0;

    // here we compute the amplitude multiplied by t -> explains the i+1 as first argument
    for (int i = 0; i < Q_order(partPair) + 1; i++)
    {
        res += Q(partPair, i, s) * speInt.IntJttComp(i + 1, 1, 1, m0 * m0, w0 * m0, m1 * m1, w1 * m1, tmin(s), 0);
        //std::cout << Q(partPair, i, s) << " " << speInt.IntJttComp(i + 1, 1, 1, m0 * m0, w0 * m0, m1 * m1, w1 * m1, tmin(s), 0) << std::endl;
    }

    //std::cout << "ici : "  << Q(partPair, 0, s) << " " << Q(partPair, 1, s) << " " << Q(partPair, 2, s) << " " <<  tmin(s) << " " << res << std::endl;

    return res;
}

dcomp CrossSection::IntRPtPu(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    double r1 = m1 * m1 + s - mandSum;

    dcomp res = 0;

    // here we compute the amplitude multiplied by t -> explains the i+1 as first argument
    for (int i = 0; i < Q_order(partPair) + 1; i++)
        res += Q(partPair, i, s) * speInt.IntJtuComp(i + 1, 1, 1, m0 * m0, w0 * m0, r1, w1 * m1, tmin(s), 0);

    return res;
}

dcomp CrossSection::IntRPuPu(ParticlePair const &partPair, double s)
{

    double m0 = partPair.get_particle(0).get_mass();
    double m1 = partPair.get_particle(1).get_mass();

    double w0 = partPair.get_particle(0).get_width();
    double w1 = partPair.get_particle(1).get_width();

    double r0 = m0 * m0 + s - mandSum;
    double r1 = m1 * m1 + s - mandSum;

    dcomp res = 0;

    // here we compute the amplitude multiplied by t -> explains the i+1 as first argument
    for (int i = 0; i < Q_order(partPair) + 1; i++)
        res += Q(partPair, i, s) * speInt.IntJuuComp(i + 1, 1, 1, r0, w0 * m0, r1, w1 * m1, tmin(s), 0);

    return res;
}

dcomp CrossSection::IntRP(ParticlePair const &partPair, double s)
{

    double m = partPair.get_second().get_mass();
    double w = partPair.get_second().get_width();

    INT_TYPE type = partPair.get_type();

    if (type != INT_TYPE::ST && type != INT_TYPE::SU)
    {
        std::cout << "FATAL ERROR : CrossSectionScatt::IntRP should not be used in this context " << std::endl;
        exit(0);
    }

    dcomp res = 0;

    // here we compute the amplitude multiplied by t -> explains the i+1 as first argument
    for (int i = 0; i < Q_order(partPair) + 1; i++)
    {
        if (type == INT_TYPE::ST)
            res += Q(partPair, i, s) * speInt.IntJtComp(i + 1, 1, m * m, w * m, tmin(s), 0);
        if (type == INT_TYPE::SU)
        {
            double r = m * m + s - mandSum;
            res += Q(partPair, i, s) * speInt.IntJuComp(i + 1, 1, r, w * m, tmin(s), 0);
        }
    }

    return res;
}

dcomp CrossSection::IntR1(ParticlePair const &partPair, double s)
{
    dcomp res = 0;

    // here we compute the amplitude multiplied by t -> explan the i+1 as first argument
    for (int i = 0; i < Q_order(partPair) + 1; i++)
    {
        res += Q(partPair, i, s) * speInt.IntK(i + 1, tmin(s), 0);
    }

    return res;
}

//
//
//
//
// Other functions
//

// Definition of tmin and tmax
double CrossSection::tmin(double s)
{

    double p1cm = f_p1cm(s);
    double p3cm = (1. / (2. * sqrt(s))) * sqrt((s - pow(_m3 - _m4, 2)) * (s - pow(_m3 + _m4, 2)));
    double E1cm = sqrt(_m1 * _m1 + p1cm * p1cm);
    double E3cm = sqrt(_m3 * _m3 + p3cm * p3cm);

    return _m1 * _m1 + _m3 * _m3 - 2 * (E1cm * E3cm + p1cm * p3cm);
}

double CrossSection::tmax(double s)
{
    double p1cm = f_p1cm(s);
    double p3cm = (1. / (2. * sqrt(s))) * sqrt((s - pow(_m3 - _m4, 2)) * (s - pow(_m3 + _m4, 2)));
    double E1cm = sqrt(_m1 * _m1 + p1cm * p1cm);
    double E3cm = sqrt(_m3 * _m3 + p3cm * p3cm);

    return _m1 * _m1 + _m3 * _m3 - 2 * (E1cm * E3cm - p1cm * p3cm);
}

double CrossSection::f_p1cm(double s)
{
    double pow_m1m22 = pow(_m1 + _m2, 2);
    double p1cm;

    if (s / pow_m1m22 < 1 && s / pow_m1m22 > 1 - 1e-8 / pow_m1m22)
        return 0;
    else
        return (1. / (2. * sqrt(s))) * sqrt((s - pow(_m1 - _m2, 2)) * (s - pow_m1m22));
}

double CrossSection::f_p3cm(double s)
{
    return (1. / (2. * sqrt(s))) * sqrt((s - pow(_m3 - _m4, 2)) * (s - pow(_m3 + _m4, 2)));
}

// Plotting function to plot the cross section with s
void CrossSection::plot(double var_max, std::string name, std::string var_type)
{
    std::ofstream outfile;
    std::string outfile_name = "../output/SimplifiedModels/CrossSection/CrossSection_" + name + ".out";
    outfile.open(outfile_name);
    outfile.precision(8);

    int Npts = 1000;
    double var_min, dvar, var, om;

    if (var_type == "vrel")
    {
        var_min = 0;
        dvar = (var_max - var_min) / (1.0 * Npts);
        outfile << "# vrel | vrel*sigma(vrel) [GeV^{-2}] | s-wave term [GeV^{-2}] " << std::endl;
    }
    else // if not velocity we use the Mandelstam variable s
    {
        var_min = (_m1 + _m2) * (_m1 + _m2);
        dvar = (log10(var_max) - log10(var_min)) / (1.0 * Npts);
        outfile << "# s [GeV^2] | sigma(s) [GeV^{-2}] | s-wave term [GeV^{-2}]" << std::endl;
    }

    outfile << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        if (var_type == "vrel")
        {
            var = var_min + i * dvar;
            outfile << var << "\t" << var * EvalSigma_fVrel(var) << "\t" << EvalSWaveTerm() << std::endl;
        }
        else
        {
            var = var_min * pow(10, i * dvar);
            outfile << var << "\t" << EvalSigmaP1cm2(var) << "\t" << EvalSWaveTerm() << std::endl;
        }
    }
    //*(var-4*_m1*_m1)*sqrt(var)*BesselK1(sqrt(var)/(1.62378*1e-1)) 
   
    outfile.close();
}

// Plotting function to plot the cross section with s
void CrossSection::plot_fromIntegral(double var_max, std::string name, std::string var_type)
{
    std::ofstream outfile;
    std::string outfile_name = "../output/CrossSection/CrossSection_fromIntegral_" + name + ".out";
    outfile.open(outfile_name);
    outfile.precision(8);

    int Npts = 1000;
    double var_min, dvar, var, om;

    if (var_type == "vrel")
    {
        var_min = 0;
        dvar = (var_max - var_min) / (1.0 * Npts);
        outfile << "# vrel | vrel*sigma(vrel) [GeV^{-2}] | s-wave term [GeV^{-2}] " << std::endl;
    }
    else // if not velocity we use the Mandelstam variable s
    {
        var_min = (_m1 + _m2) * (_m1 + _m2);
        dvar = (log10(var_max) - log10(var_min)) / (1.0 * Npts);
        outfile << "# s [GeV^2] | sigma(s) [GeV^{-2}] | s-wave term [GeV^{-2}]" << std::endl;
    }

    outfile << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        if (var_type == "vrel")
        {
            var = var_min + i * dvar;
            outfile << var << "\t" << var * EvalSigma_fromIntegral_fVrel(var) << "\t" << std::endl;
        }
        else
        {
            var = var_min * pow(10, i * dvar);
            outfile << var << "\t" << EvalSigma_fromIntegral(var) << std::endl;
        }
    }

    outfile.close();
}


void CrossSection::plot_sigmaTransfer(double var_max, std::string name, std::string var_type)
{
    std::ofstream outfile;
    std::string outfile_name = "../output/CrossSection/CrossSectionTransfer_" + name + ".out";
    outfile.open(outfile_name);
    outfile.precision(8);

    int Npts = 1000;
    double var_min, dvar, var, om;

    if (var_type == "vrel")
    {
        var_min = 0;
        dvar = (var_max - var_min) / (1.0 * Npts);
        outfile << "# vrel | sigma_T(vrel) [GeV^{-2}]" << std::endl;
    }
    else // if not velocity we use the Mandelstam variable s
    {
        var_min = (_m1 + _m2) * (_m1 + _m2);
        dvar = (log10(var_max) - log10(var_min)) / (1.0 * Npts);
        outfile << "# s [GeV^2] | sigma_T(s) [GeV^{-2}]" << std::endl;
    }

    outfile << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        if (var_type == "vrel")
        {
            var = var_min + i * dvar;
            om = (s_vs_vrel(var) - _m1 * _m1 - _m2 * _m2) / (2 * _m1);
            outfile << var << "\t" << EvalSigmaTransferBringmannP1cm4(om) / pow(f_p1cm(s_vs_vrel(var)), 4) << "\t" << EvalSigmaTransferBringmannP1cm4(om) << "\t" << f_p1cm(s_vs_vrel(var)) << std::endl;
        }
        else
        {
            var = var_min * pow(10, i * dvar);
            om = (var - _m1 * _m1 - _m2 * _m2) / (2 * _m1);
            outfile << var << "\t" << EvalSigmaTransferBringmannP1cm4(om) / pow(f_p1cm(var), 4) << "\t" << EvalSigmaTransferBringmannP1cm4(om) << "\t" << f_p1cm(var) << std::endl;
        }
    }

    outfile.close();
}


void CrossSection::plot_sigmaTransfer_fromIntegral(double var_max, std::string name, std::string var_type)
{
    std::ofstream outfile;
    std::string outfile_name = "../output/CrossSection/CrossSectionTransfer_fromIntegral_" + name + ".out";
    outfile.open(outfile_name);
    outfile.precision(8);

    int Npts = 1000;
    double var_min, dvar, var, om;

    if (var_type == "vrel")
    {
        var_min = 0;
        dvar = (var_max - var_min) / (1.0 * Npts);
        outfile << "# vrel | sigma_T(vrel) [GeV^{-2}]" << std::endl;
    }
    else // if not velocity we use the Mandelstam variable s
    {
        var_min = (_m1 + _m2) * (_m1 + _m2);
        dvar = (log10(var_max) - log10(var_min)) / (1.0 * Npts);
        outfile << "# s [GeV^2] | sigma_T(s) [GeV^{-2}]" << std::endl;
    }

    outfile << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        if (var_type == "vrel")
        {
            var = var_min + i * dvar;
            om = (s_vs_vrel(var) - _m1 * _m1 - _m2 * _m2) / (2 * _m1);
            outfile << var << "\t" << EvalSigmaTransferBringmannP1cm4_fromIntegral(om) / pow(f_p1cm(s_vs_vrel(var)), 4) << "\t" << EvalSigmaTransferBringmannP1cm4_fromIntegral(om) << "\t" << f_p1cm(s_vs_vrel(var)) << std::endl;
        }
        else
        {
            var = var_min * pow(10, i * dvar);
            om = (var - _m1 * _m1 - _m2 * _m2) / (2 * _m1);
            outfile << var << "\t" << EvalSigmaTransferBringmannP1cm4_fromIntegral(om) / pow(f_p1cm(var), 4) << "\t" << EvalSigmaTransferBringmannP1cm4_fromIntegral(om) << "\t" << f_p1cm(var) << std::endl;
        }
    }

    outfile.close();
}


// Plotting function to plot the tmin and tmax with s
void CrossSection::plot_tminAndtmax(double s_max)
{
    std::ofstream f_t;
    f_t.open("../output/tminAndtmax.out");
    f_t.precision(8);

    double s_min = (_m1 + _m2) * (_m1 + _m2);

    int Npts = 10000;
    double logdel_s = (log10(s_max) - log10(s_min)) / (1.0 * Npts);

    double s;

    f_t << "# s [GeV^2] \t tmin(s) [GeV^2}] \t tmax(s) [GeV^2]" << std::endl;
    f_t << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        s = s_min * pow(10, i * logdel_s);
        f_t << s << "\t" << tmin(s) << "\t" << tmax(s) << std::endl;
    }

    f_t.close();
}

// Plotting function to plot L1
void CrossSection::plot_L1(ParticlePair const &partPair, double s_max)
{
    std::ofstream f_L;
    f_L.open("../output/L1.out");
    f_L.precision(8);

    double s_min = (_m1 + _m2) * (_m1 + _m2);

    int Npts = 10000;
    double logdel_s = (log10(s_max) - log10(s_min)) / (1.0 * Npts);

    double s;

    f_L << "# s [GeV^2] \t L1(s) [...]" << std::endl;
    f_L << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        s = s_min * pow(10, i * logdel_s);
        f_L << s << "\t" << real(IntL1(partPair, s)) << std::endl;
    }

    f_L.close();
}

// Plotting function to plot L1
void CrossSection::plot_RealQ(ParticlePair const &partPair, int order, double s_max)
{
    std::ofstream f_Q;
    f_Q.open("../output/Q.out");
    f_Q.precision(8);

    double s_min = (_m1 + _m2) * (_m1 + _m2);

    int Npts = 10000;
    double logdel_s = (log10(s_max) - log10(s_min)) / (1.0 * Npts);

    double s;

    f_Q << "# s [GeV^2] \t Q(s)" << std::endl;
    f_Q << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        s = s_min * pow(10, i * logdel_s);
        f_Q << s << "\t" << real(Q(partPair, order, s)) << std::endl;
    }

    f_Q.close();
}

//
//
//
//
//
//
//
//
// ==================================================================
// Function that are used to test all the possible interference terms
// ==================================================================

double CrossSection::EvalSigmaInterfChannelSP1cm2(double s)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsSSInt.size(); i++)
    {
        m0 = ParticlePairsSSInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSInt[i].get_particle(0).get_width();

        m1 = ParticlePairsSSInt[i].get_particle(1).get_mass();
        w1 = ParticlePairsSSInt[i].get_particle(1).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0) * ((s - m1 * m1) + one_i * m1 * w1)) * IntL1(ParticlePairsSSInt[i], s));
    }

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    //std::cout << "res en bas : " << res << std::endl;

    return res;
}

double CrossSection::EvalSWaveTermInterfChannelS()
{
    double res = 0;

    double m0, w0, m1, w1;

    // s is fixed at its minimal value
    double s = pow(_m1 + _m2, 2);

    if (s <= pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsSSInt.size(); i++)
    {
        m0 = ParticlePairsSSInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSInt[i].get_particle(0).get_width();

        m1 = ParticlePairsSSInt[i].get_particle(1).get_mass();
        w1 = ParticlePairsSSInt[i].get_particle(1).get_width();

        res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0) * ((s - m1 * m1) + one_i * m1 * w1)) * IntLH1(ParticlePairsSSInt[i], s));
    }

    res /= (64 * PI * sqrt(s) * _g1 * _g2) * _sym * (pow(_m1 + _m2, 2) / 4);

    return res;
}

double CrossSection::EvalSigmaInterfChannelTP1cm2(double s)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsTTInt.size(); i++)
        res += 2 * real(IntLPtPt(ParticlePairsTTInt[i], s));

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::EvalSigmaInterfChannelUP1cm2(double s)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsUUInt.size(); i++)
        res += 2 * real(IntLPuPu(ParticlePairsUUInt[i], s));

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::EvalSigmaInterfChannelSTP1cm2(double s, bool turnOnForSameParticle)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsSTInt.size(); i++)
    {
        m0 = ParticlePairsSTInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSTInt[i].get_particle(0).get_width();

        if (ParticlePairsSTInt[i].get_particle(0).get_index() != ParticlePairsSTInt[i].get_particle(1).get_index() || turnOnForSameParticle)
            res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntLP(ParticlePairsSTInt[i], s));
        //else
        //    std::cout << "on est bien dans ce cas parfois" << std::endl;
    }

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::EvalSigmaInterfChannelSUP1cm2(double s, bool turnOnForSameParticle)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsSUInt.size(); i++)
    {
        m0 = ParticlePairsSUInt[i].get_particle(0).get_mass();
        w0 = ParticlePairsSUInt[i].get_particle(0).get_width();

        if (ParticlePairsSUInt[i].get_particle(0).get_index() != ParticlePairsSUInt[i].get_particle(1).get_index() || turnOnForSameParticle)
            res += 2 * real(1. / (((s - m0 * m0) - one_i * m0 * w0)) * IntLP(ParticlePairsSUInt[i], s));
    }

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::EvalSigmaInterfChannelTUP1cm2(double s, bool turnOnForSameParticle)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsTUInt.size(); i++)
        if (ParticlePairsTUInt[i].get_particle(0).get_index() != ParticlePairsTUInt[i].get_particle(1).get_index() || turnOnForSameParticle)
            res += 2 * real(IntLPtPu(ParticlePairsTUInt[i], s));

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::EvalSigmaSquareChannelSP1cm2(double s)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsSSSq.size(); i++)
    {
        m0 = ParticlePairsSSSq[i].get_particle(0).get_mass();
        w0 = ParticlePairsSSSq[i].get_particle(0).get_width();

        res += 1. / (pow(s - m0 * m0, 2) + m0 * m0 * w0 * w0) *
               real(IntL1(ParticlePairsSSSq[i], s));
    }

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::EvalSigmaSquareChannelTP1cm2(double s)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsTTSq.size(); i++)
        res += real(IntLPtPt(ParticlePairsTTSq[i], s));

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::EvalSigmaSquareChannelUP1cm2(double s)
{
    double res = 0;

    double m0, w0, m1, w1;

    if (s < pow(_m3 + _m4, 2))
        return 0;

    for (int i = 0; i < ParticlePairsUUSq.size(); i++)
        res += real(IntLPuPu(ParticlePairsUUSq[i], s));

    res /= (64 * PI * s * _g1 * _g2) * _sym;

    return res;
}

double CrossSection::s_vs_vrel(double vrel)
{

    //    if(vrel < 1e-3)
    //        return pow(_m1 +_m2, 2) + _m1*_m2*vrel*vrel;

    double gamma_rel = 1 / sqrt(1 - vrel * vrel);
    return pow(_m1 - _m2, 2) + 2 * _m1 * _m2 * (1 + gamma_rel);
}
