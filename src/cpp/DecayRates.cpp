#include "../headers/DecayRates.h"

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

/** Decay rate of the prop_s scalar mediator
* For the quarks we only consider them if the mass of the decaying particle is higher than T_QCD
* Otherwise we should consider the decay into mesons and nucleons
* Note that as nucleons and meson masses are at least of the order of T_QCD these terms
* will be neglected in a first approximation */
double DecayRates::DecayRate_sff(DarkSectorModel const &_DS, int prop_s)
{
    double res = 0;

    //std::cout << "ici" << std::endl;

    int n_SM_couples_Fermions = StandardModel::getInstance()->get_n_couples_SM_ferm();
    std::vector<double> lambda_S_SM = _DS.get_lambda_S_SM(prop_s);

    double ms = _DS.get_phis(prop_s).get_mass();
    double mf = 0;

    QCDTYPE type;

    for (int i = 0; i < n_SM_couples_Fermions; i++)
    {

        mf = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_mass();
        type = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_QCDType();

        if (2 * mf < ms)
        {
            if (type != QCDTYPE::QUARK)
                res += pow(lambda_S_SM[i], 2) * ms * pow(1 - 4 * mf * mf / (ms * ms), 3 / 2);
            else
            {
                //std::cout << mf << std::endl;
                res += 3 * regulatorTQCD(log10(ms), 1.0 / 0.3, log10(T_transQCD), 1, 0) * pow(lambda_S_SM[i], 2) * ms * pow(1 - 4 * mf * mf / (ms * ms), 3. / 2.);
            }
            // The mass of the scalar needs to be higher than the QCD phase transition ... otherwise decay is not into quarks directly
            // The factor 3 here takes into account the 3 colors of the quarks
        }
    }

    //std::cout << " " << std::endl;
    res = res / (8 * PI);

    return res;
}

// Decay-rate of scalar in WIMPs
double DecayRates::DecayRate_scc(DarkSectorModel const &_DS, int prop_s)
{
    double res = 0;

    std::vector<std::vector<double> > lambda_S_DM = _DS.get_lambda_S_DM(prop_s); // coupling to the mediator
    double ms = _DS.get_phis(prop_s).get_mass();                                 // mass of the mediator

    double n_chi = _DS.get_n_DS_DM(); // number of dark matter particles

    double mi, mj;

    for (int i = 0; i < n_chi; i++)
    {
        mi = _DS.get_chi(i).get_mass();

        // We devide by two if WIMPs are majorana when they are the same : particle = anti-particle
        if (2 * mi < ms && _DS.get_chi(i).get_fermiontype() == Fermiontype::majorana)
            res += pow(lambda_S_DM[i][i], 2) * ms * pow(1 - 4 * mi * mi / (ms * ms), 3 / 2) / 4;

        if (2 * mi < ms && _DS.get_chi(i).get_fermiontype() == Fermiontype::dirac)
            res += pow(lambda_S_DM[i][i], 2) * ms * pow(1 - 4 * mi * mi / (ms * ms), 3 / 2) / 2;

        for (int j = 0; j < i; j++)
        {
            mj = _DS.get_chi(j).get_mass();

            if (mi + mj < ms)
                res += pow(lambda_S_DM[i][j], 2) * (ms * ms - (mi + mj) * (mi + mj)) * pf(ms, mi, mj);
        }
    }

    res = res / (4 * PI);

    return res;
}

// Decay-rate of scalar in scalar + scalar
double DecayRates::DecayRate_sss(DarkSectorModel const &_DS, int prop_s)
{
    double res = 0;

    std::vector<std::vector<double> > c_sss = _DS.get_c_sss(prop_s); // coupling to the mediator
    double ms_0 = _DS.get_phis(prop_s).get_mass();                   // mass of the mediator

    double n_s = _DS.get_n_DS_phis(); // number of scalar propagators

    double ms_1, ms_2;

    for (int i = 0; i < n_s; i++)
    {
        ms_1 = _DS.get_phis(i).get_mass();

        if (2 * ms_1 < ms_0)
            res += pow(c_sss[i][i] * v, 2) * pf(ms_0, ms_1, ms_1) / (2 * ms_0 * ms_0);

        for (int j = 0; j < i; j++)
        {
            ms_2 = _DS.get_phis(j).get_mass();

            if (ms_1 + ms_2 < ms_0)
                res += pow(c_sss[i][j] * v, 2) * pf(ms_0, ms_1, ms_2) / (ms_0 * ms_0);
        }
    }

    res = res / (8 * PI);
    return res;
}

// Decay-rate of scalar in pseudoscalar + pseudoscalar
double DecayRates::DecayRate_spp(DarkSectorModel const &_DS, int prop_s)
{
    double res = 0;

    std::vector<std::vector<double> > d_spp = _DS.get_d_spp(prop_s); // coupling to the mediator
    double ms_0 = _DS.get_phis(prop_s).get_mass();                   // mass of the mediator

    double n_p = _DS.get_n_DS_phip(); // number of scalar propagators

    double mp_1, mp_2;

    for (int i = 0; i < n_p; i++)
    {
        mp_1 = _DS.get_phip(i).get_mass();

        if (2 * mp_1 < ms_0)
            res += pow(d_spp[i][i] * v, 2) * pf(ms_0, mp_1, mp_1) / (2 * ms_0 * ms_0);

        for (int j = 0; j < i; j++)
        {
            mp_2 = _DS.get_phip(j).get_mass();

            if (mp_1 + mp_2 < ms_0)
                res += pow(d_spp[i][j] * v, 2) * pf(ms_0, mp_1, mp_2) / (ms_0 * ms_0);
        }
    }

    res = res / (8 * PI);

    return res;
}

// Decay-rate of scalar in vector + vector
double DecayRates::DecayRate_sXX(DarkSectorModel const &_DS, int prop_s)
{
    double res = 0;

    std::vector<std::vector<double> > g_sXX = _DS.get_g_sXX(prop_s); // coupling to the mediator
    double ms = _DS.get_phis(prop_s).get_mass();                     // mass of the mediator

    double n_X = _DS.get_n_DS_X(); // number of scalar propagators

    double mX = 0;

    for (int i = 0; i < n_X; i++)
    {
        mX = _DS.get_X(i).get_mass();

        if (2 * mX < ms && mX != 0)
            res += pow(g_sXX[i][i], 2) * (v * v / (ms * ms)) * (3 + pow(ms / mX, 2) * (0.25 * pow(ms / mX, 2) - 1)) * sqrt(1 - 4 * mX * mX / (ms * ms)) * ms;
        else if (2 * mX < ms && mX == 0)
            res += pow(g_sXX[i][i], 2) * (v * v / (ms * ms)) * 3 * ms;
    }

    res = res / (32 * PI);

    return res;
}

// Decay rate into two gluons
// Here the temperature matters as above a given threshold we know longer have decay into gluons
double DecayRates::DecayRate_sgluongluon(DarkSectorModel const &_DS, int prop_s)
// Attention here we use the value tabulated at M_Z
{
    double res = 0;
    //std::cout << "ici" << std::endl;

    int n_SM_couples_Fermions = StandardModel::getInstance()->get_n_couples_SM_ferm();
    std::vector<double> lambda_S_SM = _DS.get_lambda_S_SM(prop_s);

    double ms = _DS.get_phis(prop_s).get_mass();
    double mf = 0;

    dcomp amp = 0;
    double xf = 0;
    std::string name = "";
    QCDTYPE type;

    for (int i = 0; i < n_SM_couples_Fermions; i++)
    {
        mf = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_mass();
        name = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_name();
        type = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_QCDType();
        xf = 4 * mf * mf / (ms * ms);

        //std::cout << name << std::endl;

        if (xf > 0 && name != "quark_up" && name != "quark_down" && name != "quark_str" && type == QCDTYPE::QUARK)
        {
            //std::cout << name << " " << lambda_S_SM[i]/mf*M_QUARK_TOP  << std::endl;
            amp += lambda_S_SM[i] * fOneLoopScalarVectors(xf) / mf;
        }
    }

    //std::cout << pow(abs(amp) * ALPHA_S/PI* M_QUARK_TOP*10, 2)  << std::endl;
 
    // std::cout << ms << " " << amp << std::endl;
    res = regulatorTQCD(log10(ms), 1.0 / 0.3, log10(T_transQCD), 1, 0) * pow(abs(amp) * ALPHA_S, 2) * pow(ms, 3) / (32 * PI * PI * PI);

    return res;
}

// Decay rate into two gluons
double DecayRates::DecayRate_sphotonphoton(DarkSectorModel const &_DS, int prop_s)
{
    double res = 0;
    //std::cout << "ici" << std::endl;

    int n_SM_couples_Fermions = StandardModel::getInstance()->get_n_couples_SM_ferm();
    std::vector<double> lambda_S_SM = _DS.get_lambda_S_SM(prop_s);

    double ms = _DS.get_phis(prop_s).get_mass();
    double mf = 0;

    dcomp amp = 0;
    double xf = 0;
    double qf2 = 0;
    std::string name = "";

    for (int i = 0; i < n_SM_couples_Fermions; i++)
    {
        mf = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_mass();
        qf2 = pow(StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_electric_charge(), 2);
        name = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_name();
        xf = 4 * mf * mf / (ms * ms);

        if (xf > 0 && name != "quark_up" && name != "quark_down" && name != "quark_str")
        {
            if (StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_QCDType() == QCDTYPE::QUARK)
                amp += 3 * qf2 * lambda_S_SM[i] * fOneLoopScalarVectors(xf) / mf;
            else
                amp += qf2 * lambda_S_SM[i] * fOneLoopScalarVectors(xf) / mf;
        }
    }

    res = pow(abs(amp) * ALPHA_EM, 2) * pow(ms, 3) / (64 * PI * PI * PI);
    return res;
}

//
//
//
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//

// decay rate of the prop_s scalar mediator
double DecayRates::DecayRate_pff(DarkSectorModel const &_DS, int prop_p)
{
    double res = 0;

    //std::cout << "ici" << std::endl;

    int n_SM_couples_Fermions = StandardModel::getInstance()->get_n_couples_SM_ferm();
    std::vector<double> lambda_PS_SM = _DS.get_lambda_PS_SM(prop_p);

    double mp = _DS.get_phip(prop_p).get_mass();
    double mf = 0;
    QCDTYPE type;

    for (int i = 0; i < n_SM_couples_Fermions; i++)
    {
        mf = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_mass();
        type = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_QCDType();

        if (2 * mf < mp)
        {
            if (type != QCDTYPE::QUARK)
                res += pow(lambda_PS_SM[i], 2) * mp * pow(1 - 4 * mf * mf / (mp * mp), 1 / 2);
            else
                res += 3 * regulatorTQCD(log10(mp), 1.0 / 0.3, log10(T_transQCD), 1, 0) * pow(lambda_PS_SM[i], 2) * mp * pow(1 - 4 * mf * mf / (mp * mp), 1 / 2);
            // The factor 3 here takes into account the 3 colors of the quarks
        }
    }

    res = res / (8 * PI);

    return res;
}

// Attention that is done for Majorana particles
// Decay-rate of scalar in WIMPs
double DecayRates::DecayRate_pcc(DarkSectorModel const &_DS, int prop_p)
{
    double res = 0;

    std::vector<std::vector<double> > lambda_PS_DM = _DS.get_lambda_PS_DM(prop_p); // coupling to the mediator
    double mp = _DS.get_phip(prop_p).get_mass();                                   // mass of the mediator

    double n_chi = _DS.get_n_DS_DM(); // number of dark matter particles

    double mi, mj;

    for (int i = 0; i < n_chi; i++)
    {
        mi = _DS.get_chi(i).get_mass();

        // We devide by two if WIMPs are majorana when they are the same : particle = anti-particle
        if (2 * mi < mp && _DS.get_chi(i).get_fermiontype() == Fermiontype::majorana)
            res += pow(lambda_PS_DM[i][i], 2) * mp * pow(1 - 4 * mi * mi / (mp * mp), 1 / 2) / 4;

        if (2 * mi < mp && _DS.get_chi(i).get_fermiontype() == Fermiontype::dirac)
            res += pow(lambda_PS_DM[i][i], 2) * mp * pow(1 - 4 * mi * mi / (mp * mp), 1 / 2) / 2;

        for (int j = 0; j < i; j++)
        {
            mj = _DS.get_chi(j).get_mass();

            if (mi + mj < mp)
                res += pow(lambda_PS_DM[i][j], 2) * (mp * mp - (mi - mj) * (mi - mj)) * pf(mp, mi, mj);
        }
    }

    res = res / (4 * PI);

    return res;
}

// Attention that is done for Majorana particles
// Decay-rate of scalar in scalar +  scalar
double DecayRates::DecayRate_psp(DarkSectorModel const &_DS, int prop_p)
{
    double res = 0;

    std::vector<std::vector<double> > d_psp = _DS.get_d_psp(prop_p); // coupling to the mediator
    double mp_0 = _DS.get_phip(prop_p).get_mass();                   // mass of the mediator

    double n_s = _DS.get_n_DS_phis(); // number of scalar propagators
    double n_p = _DS.get_n_DS_phip(); // number of pseudoscalar propagators

    double ms_1, mp_2;

    for (int i = 0; i < n_s; i++)
    {
        ms_1 = _DS.get_phis(i).get_mass();

        for (int j = 0; j < n_p; j++)
        {
            mp_2 = _DS.get_phip(j).get_mass();

            if (ms_1 + mp_2 < mp_0)
                res += pow(d_psp[i][j], 2) * v * v * pf(mp_0, ms_1, mp_2) / (mp_0 * mp_0);
        }
    }

    res = res / (8 * PI);

    return res;
}

// Decay rate into two gluons
// Here the temperature matters as abovea given threshold we know longer have decay into gluons
double DecayRates::DecayRate_pgluongluon(DarkSectorModel const &_DS, int prop_p)
// Attention here we use the value tabulated at M_Z
{
    double res = 0;
    //std::cout << "ici" << std::endl;

    int n_SM_couples_Fermions = StandardModel::getInstance()->get_n_couples_SM_ferm();
    std::vector<double> lambda_PS_SM = _DS.get_lambda_PS_SM(prop_p);

    double mp = _DS.get_phip(prop_p).get_mass();
    double mf = 0;

    dcomp amp = 0;
    double xf = 0;
    std::string name = "";
    QCDTYPE type;

    for (int i = 0; i < n_SM_couples_Fermions; i++)
    {
        mf = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_mass();
        name = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_name();
        type = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_QCDType();
        xf = 4 * mf * mf / (mp * mp);

        //std::cout << name << std::endl;

        if (xf > 0 && name != "quark_up" && name != "quark_down" && name != "quark_str" && type == QCDTYPE::QUARK)
        {
            //std::cout << name << std::endl;
            amp += lambda_PS_SM[i] * fOneLoopPseudoScalarVectors(xf) / mf;
        }
    }

    // std::cout << ms << " " << amp << std::endl;
    res = regulatorTQCD(log10(mp), 1.0 / 0.3, log10(T_transQCD), 1, 0) * pow(abs(amp) * ALPHA_S, 2) * pow(mp, 3) / (32 * PI * PI * PI);

    return res;
}

// Decay rate into two photons
double DecayRates::DecayRate_pphotonphoton(DarkSectorModel const &_DS, int prop_p)
{
    double res = 0;
    //std::cout << "ici" << std::endl;

    int n_SM_couples_Fermions = StandardModel::getInstance()->get_n_couples_SM_ferm();
    std::vector<double> lambda_PS_SM = _DS.get_lambda_PS_SM(prop_p);

    double mp = _DS.get_phip(prop_p).get_mass();
    double mf = 0;

    dcomp amp = 0;
    double xf = 0;
    double qf2 = 0;
    std::string name = "";

    for (int i = 0; i < n_SM_couples_Fermions; i++)
    {
        mf = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_mass();
        qf2 = pow(StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_electric_charge(), 2);
        name = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_name();
        xf = 4 * mf * mf / (mp * mp);

        if (xf > 0 && name != "quark_up" && name != "quark_down" && name != "quark_str")
        {
            if (StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_QCDType() == QCDTYPE::QUARK)
                amp += 3 * qf2 * lambda_PS_SM[i] * fOneLoopPseudoScalarVectors(xf) / mf;
            else
                amp += qf2 * lambda_PS_SM[i] * fOneLoopPseudoScalarVectors(xf) / mf;
        }
    }

    res = pow(abs(amp) * ALPHA_EM, 2) * pow(mp, 3) / (64 * PI * PI * PI);

    return res;
}

//
//
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//

double DecayRates::DecayRate_Xff(DarkSectorModel const &_DS, int prop_X)
{
    double res = 0;

    //std::cout << "ici" << std::endl;

    int n_SM_couples_Fermions = StandardModel::getInstance()->get_n_couples_SM_ferm();
    std::vector<double> a_VEC_SM = _DS.get_a_VEC_SM(prop_X);
    std::vector<double> b_VEC_SM = _DS.get_b_VEC_SM(prop_X);

    double mX = _DS.get_X(prop_X).get_mass();
    double mf = 0;

    for (int i = 0; i < n_SM_couples_Fermions; i++)
    {
        mf = StandardModel::getInstance()->get_couples_SM_ferm(i).get_first().get_mass();

        if (2 * mf < mX)
        {
            res += (pow(a_VEC_SM[i], 2) * (1 + 2 * (mf * mf) / (mX * mX)) + pow(b_VEC_SM[i], 2) * (1 - 4 * (mf * mf) / (mX * mX))) * mX * sqrt(1 - 4 * (mf * mf) / (mX * mX));
            //std::cout << i << " " << lambda_S_SM[i] << std::endl;
        }
    }

    res = res / (12 * PI);

    return res;
}

double DecayRates::DecayRate_Xcc(DarkSectorModel const &_DS, int prop_X)
{
    double res = 0;

    std::vector<std::vector<double> > a_VEC_DM = _DS.get_a_VEC_DM(prop_X); // couplings to the mediator
    std::vector<std::vector<double> > b_VEC_DM = _DS.get_b_VEC_DM(prop_X);

    double mX = _DS.get_X(prop_X).get_mass(); // mass of the mediator

    double n_chi = _DS.get_n_DS_DM(); // number of dark matter particles

    double a2 = 0;
    double b2 = 0;

    double mi, mj;

    for (int i = 0; i < n_chi; i++)
    {
        mi = _DS.get_chi(i).get_mass();

        a2 = pow(a_VEC_DM[i][i], 2);
        b2 = pow(b_VEC_DM[i][i], 2);

        // We devide by two if WIMPs are majorana when they are the same : particle = anti-particle
        if (2 * mi < mX && _DS.get_chi(i).get_fermiontype() == Fermiontype::majorana)
            res += (a2 * (1 + 2 * mi * mi / (mX * mX)) + b2 * (1 - 4 * mi * mi / (mX * mX))) * mX * sqrt(1 - 4 * mi * mi / (mX * mX));

        if (2 * mi < mX && _DS.get_chi(i).get_fermiontype() == Fermiontype::dirac)
            res += (a2 * (1 + 2 * mi * mi / (mX * mX)) + b2 * (1 - 4 * mi * mi / (mX * mX))) * mX * sqrt(1 - 4 * mi * mi / (mX * mX)) * 2;

        for (int j = 0; j < i; j++)
        {
            mj = _DS.get_chi(j).get_mass();

            a2 = pow(a_VEC_DM[i][j], 2);
            b2 = pow(b_VEC_DM[i][j], 2);

            if (mi + mj < mX)
                res += (4 * a2 * (mX * mX + 2 * mi * mj - (mi - mj) * (mi - mj)) + 4 * b2 * (mX * mX - 2 * mi * mj - (mi + mj) * (mi + mj)) + 2 * a2 * ((mi - mj) * (mi - mj) / (mX * mX)) * (mX * mX - (mi + mj) * (mi + mj)) + 2 * b2 * ((mi + mj) * (mi + mj) / (mX * mX)) * (mX * mX - (mi - mj) * (mi - mj))) * pf(mX, mi, mj) / (mX * mX);
        }
    }

    res = res / (24 * PI);

    return res;
}

//
//
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//

double DecayRates::DecayRate_csc(DarkSectorModel const &_DS, int i_chi)
{
    double res = 0;

    int n_chi = _DS.get_n_DS_DM();
    int n_s = _DS.get_n_DS_phis();

    double mi = 0, mj = 0, ms = 0;
    double lambda = 0;

    mi = _DS.get_chi(i_chi).get_mass();

    for (int j_chi = 0; j_chi < n_chi; j_chi++)
    {
        mj = _DS.get_chi(j_chi).get_mass();

        if (j_chi != i_chi)
            for (int prop_s = 0; prop_s < n_s; prop_s++)
            {
                ms = _DS.get_phis(prop_s).get_mass();
                lambda = _DS.get_lambda_S_DM(prop_s, i_chi, j_chi);

                if (mj + ms < mi)
                    res += pow(lambda, 2) * (1 + mj / mi - mj * mj / (mi * mi) + ms * ms / (mi * mi)) * pf(mi, ms, mj);
            }
    }

    res = res / (8 * pi);

    return res;
}

double DecayRates::DecayRate_cpc(DarkSectorModel const &_DS, int i_chi)
{
    double res = 0;

    int n_chi = _DS.get_n_DS_DM();
    int n_p = _DS.get_n_DS_phip();

    double mi = 0, mj = 0, mp = 0;
    double lambda = 0;

    mi = _DS.get_chi(i_chi).get_mass();

    for (int j_chi = 0; j_chi < n_chi; j_chi++)
    {
        mj = _DS.get_chi(j_chi).get_mass();

        if (j_chi != i_chi)
            for (int prop_p = 0; prop_p < n_p; prop_p++)
            {
                mp = _DS.get_phip(prop_p).get_mass();
                lambda = _DS.get_lambda_PS_DM(prop_p, i_chi, j_chi);

                if (mj + mp < mi)
                    res += pow(lambda, 2) * (1 - mj / mi - mj * mj / (mi * mi) + mp * mp / (mi * mi)) * pf(mi, mp, mj);
            }
    }

    res = res / (8 * pi);

    return res;
}

double DecayRates::DecayRate_cXc(DarkSectorModel const &_DS, int i_chi)
{
    double res = 0;

    int n_chi = _DS.get_n_DS_DM();
    int n_X = _DS.get_n_DS_X();

    double mi = 0, mj = 0, mX = 0;
    double a = 0, b = 0;

    mi = _DS.get_chi(i_chi).get_mass();

    for (int j_chi = 0; j_chi < n_chi; j_chi++)
    {
        mj = _DS.get_chi(j_chi).get_mass();

        if (j_chi != i_chi)
            for (int prop_X = 0; prop_X < n_X; prop_X++)
            {
                mX = _DS.get_X(prop_X).get_mass();
                a = _DS.get_a_VEC_DM(prop_X, i_chi, j_chi);
                b = _DS.get_b_VEC_DM(prop_X, i_chi, j_chi);

                if (mj + mX < mi)
                    res += ((a * a + b * b) * (1 - mj * mj / (mi * mi)) * (1 + (mi - mj) * (mi - mj) / (mX * mX)) + 6 * (b * b - a * a) * mj / mi) * pf(mi, mX, mj);
            }
    }

    res = res / (8 * pi);

    return res;
}

//
//
//
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//
//
double DecayRates::DecayRate_s(DarkSectorModel const &_DS, int prop_s)
{
    double Gammas = DecayRate_sff(_DS, prop_s) + DecayRate_scc(_DS, prop_s) + DecayRate_sss(_DS, prop_s) + DecayRate_spp(_DS, prop_s) + DecayRate_sXX(_DS, prop_s);
    double ms = _DS.get_phis(prop_s).get_mass();

    return Gammas;
}

double DecayRates::DecayRate_p(DarkSectorModel const &_DS, int prop_p)
{
    double Gammap = DecayRate_pff(_DS, prop_p) + DecayRate_pcc(_DS, prop_p) + DecayRate_psp(_DS, prop_p);
    double mp = _DS.get_phip(prop_p).get_mass();

    return Gammap;
}

double DecayRates::DecayRate_X(DarkSectorModel const &_DS, int prop_X)
{
    double GammaX = DecayRate_Xff(_DS, prop_X) + DecayRate_Xcc(_DS, prop_X);
    double mX = _DS.get_X(prop_X).get_mass();

    return GammaX;
}

double DecayRates::DecayRate_chi(DarkSectorModel const &_DS, int prop_chi)
{
    double GammaChi = DecayRate_csc(_DS, prop_chi) + DecayRate_cpc(_DS, prop_chi) + DecayRate_cXc(_DS, prop_chi);
    double mc = _DS.get_chi(prop_chi).get_mass();

    return GammaChi;
}

void DecayRates::set_width_phis(DarkSectorModel &_DS)
{
    double ns = _DS.get_n_DS_phis();
    std::vector<double> width;
    width.resize(ns);

    //std::cout << "We enter here : " << std::endl;
    //std::cout << ns << std::endl;

    for (int i = 0; i < ns; i++)
    {
        //std::cout << "there first" << std::endl;
        width[i] = DecayRate_s(_DS, i);
        //std::cout << "Pourquoi" << std::endl;
        _DS.get_phis_ptr(i)->set_width(width[i]);
        //std::cout << "there end" << std::endl;
    }

    for (int i = 0; i < ns; i++)
    {
        if (width[i] > 0.1 * _DS.get_phis(i).get_mass())
        {
            std::string message = "WARNING :: Width of scalar too large exception : ratio width/mass = " + std::to_string(width[i] / 0.1 * _DS.get_phis(i).get_mass());
            char msg[message.size()];
            strcpy(msg, message.c_str());
            ExceptionHandler e(msg, __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::WidthTooLarge);
            throw e;
        }
    }
    //std::cout << "here" << std::endl;
}

void DecayRates::set_width_phip(DarkSectorModel &_DS)
{
    double np = _DS.get_n_DS_phip();
    std::vector<double> width;
    width.resize(np);

    for (int i = 0; i < np; i++)
    {
        width[i] = DecayRate_p(_DS, i);
        _DS.get_phip_ptr(i)->set_width(width[i]);
    }

    for (int i = 0; i < np; i++)
    {
        if (width[i] > 0.1 * _DS.get_phip(i).get_mass())
        {
            std::string message = "WARNING :: Width of psuedoscalar too large exception : ratio width/mass = " + std::to_string(width[i] / 0.1 * _DS.get_phip(i).get_mass());
            char msg[message.size()];
            strcpy(msg, message.c_str());
            ExceptionHandler e(msg, __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::WidthTooLarge);
            throw e;
        }
    }
}

void DecayRates::set_width_X(DarkSectorModel &_DS)
{
    double nX = _DS.get_n_DS_X();
    std::vector<double> width;
    width.resize(nX);

    for (int i = 0; i < nX; i++)
    {
        width[i] = DecayRate_X(_DS, i);
        _DS.get_X_ptr(i)->set_width(width[i]);
    }

    for (int i = 0; i < nX; i++)
    {
        if (width[i] > 0.1 * _DS.get_X(i).get_mass())
        {
            std::string message = "WARNING :: Width of vector too large exception : ratio width/mass = " + std::to_string(width[i] / 0.1 * _DS.get_X(i).get_mass());
            char msg[message.size()];
            strcpy(msg, message.c_str());
            ExceptionHandler e(msg, __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::WidthTooLarge);
            throw e;
        }
    }
}

void DecayRates::set_width_chi(DarkSectorModel &_DS)
{
    double nc = _DS.get_n_DS_DM();
    double width = 0;

    for (int i = 0; i < nc; i++)
    {
        width = DecayRate_chi(_DS, i);
        _DS.get_chi_ptr(i)->set_width(width);
    }
}

double DecayRates::pf(double m0, double m1, double m2)
{
    return (1. / (2 * m0)) * sqrt((m0 * m0 - (m1 + m2) * (m1 + m2)) * (m0 * m0 - (m1 - m2) * (m1 - m2)));
}