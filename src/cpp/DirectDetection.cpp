#include "../headers/DirectDetection.h"

DirectDetection::DirectDetection(DarkSectorModel const &DS)
{

    // Factors for protons
    f_u_p = 0.0153;
    f_d_p = 0.0191;
    f_s_p = 0.0447;

    DD_u_p = 0.842;
    DD_d_p = -0.427;
    DD_s_p = -0.085;

    d_u_p = 0.84;
    d_d_p = -0.23;
    d_s_p = -0.046;

    // Factors for netrons
    f_u_n = 0.011;
    f_d_n = 0.0273;
    f_s_n = 0.0447;

    DD_u_n = -0.427;
    DD_d_n = 0.842;
    DD_s_n = -0.085;

    d_u_n = -0.23;
    d_d_n = 0.84;
    d_s_n = -0.046;

    /** 
    * We resize this table to the total number of possible
    * quark and gluon operators.
    * However to avoid any confusion we start the table at row 1
    * as operators are only counted from 1 to 15. */

    cq.resize(15);
    cg.resize(5, 0); // Note that in simplified models cg = 0, cg != 0 only in effective field theories
    cq[0].resize(6, std::nan(""));
    for (int i = 1; i < cq.size(); i++)
        cq[i].resize(6, 0); // for the six quarks playing a role in the order (u, d, s, c, b, t)

    /** *****************************************************
     * Initialisation of the coefficients coupling DM-to quarks
     * !! In this version of the code, couplings are not renormalised !!
     * This assumes that only one DM species makes all of DM today
     * ***************************************************** */
    int idm = (DS.get_chi(0).get_fermiontype() == Fermiontype::dirac) ? 1 : 0;

    for (int i = 0; i < DS.get_n_DS_phis(); i++)
        for (int q = 0; q < 6; q++)
            cq[1][q] += DS.get_lambda_S_DM(i, 0, idm) * DS.get_lambda_S_SM(i, q + 6) * pow(DS.get_phis(i).get_mass(), -2);

    for (int i = 0; i < DS.get_n_DS_phip(); i++)
        for (int q = 0; q < 6; q++)
            cq[4][q] += DS.get_lambda_PS_DM(i, 0, idm) * DS.get_lambda_PS_SM(i, q + 6) * pow(DS.get_phip(i).get_mass(), -2);

    for (int i = 0; i < DS.get_n_DS_X(); i++)
        for (int q = 0; q < 6; q++)
        {
            cq[5][q] += DS.get_a_VEC_DM(i, 0, idm) * DS.get_a_VEC_SM(i, q + 6) * pow(DS.get_X(i).get_mass(), -2); // vector/vector

            if (idm == 1)
            {
                cq[6][q] += DS.get_b_VEC_DM(i, 0, idm) * DS.get_a_VEC_SM(i, q + 6) * pow(DS.get_X(i).get_mass(), -2); // axial-vector/vector
                cq[7][q] += DS.get_a_VEC_DM(i, 0, idm) * DS.get_b_VEC_SM(i, q + 6) * pow(DS.get_X(i).get_mass(), -2); // vector/axial-vector
                cq[8][q] += DS.get_b_VEC_DM(i, 0, idm) * DS.get_b_VEC_SM(i, q + 6) * pow(DS.get_X(i).get_mass(), -2); // axial-vector/axial-vector
            }
        }

    //-> Still need to define cg[] properly (partly done in CPOddHigssLowMass.cpp)

    /** Note that here we start at q+6 for lambda as quarks are put after the leptons 
     *  in the table of SM couples (but also in the order u, d, s, c, b, t) */
    // ******************************************************

    // Nucleon mass in GeV defined in my units.h (assuming M_NEUTRON \simeq M_PROTON)
    mN = M_NEUTRON;
    // Dark Matter mass in GeV
    mchi = DS.get_chi(0).get_mass();
    // Quark masses in GeV (ordered)
    mq = {M_QUARK_UP, M_QUARK_DOWN, M_QUARK_STRANGE, M_QUARK_CHARM, M_QUARK_BOTTOM, M_QUARK_TOP};
    // Ordered list of the light quarks
    vec_Qtypes = {QuarkType::UP, QuarkType::DOWN, QuarkType::STRANGE};

    /** ****************************************
     * Initialisation of the parameters C3 and C4
     * ******************************************** */
    C3 = 0;
    C4 = 0;

    m_bar = 0;
    for (int q = 0; q < 3; q++)
    {
        C3 += cq[3][q] / mq[q];
        C4 += cq[4][q] / mq[q];
    }
    for (int q = 0; q < 3; q++)
        m_bar += 1 / mq[q];

    m_bar = 1. / m_bar;
    C3 = C3 * m_bar;
    C4 = C4 * m_bar;

    /** *****************************************************
     * Initialisation of the coefficients coupling DM-to quarks
     * ***************************************************** */
    vec_ON = {EffectiveOperator::N0, EffectiveOperator::N1, EffectiveOperator::N2, EffectiveOperator::N3,
              EffectiveOperator::N4, EffectiveOperator::N5, EffectiveOperator::N6, EffectiveOperator::N7,
              EffectiveOperator::N8, EffectiveOperator::N9, EffectiveOperator::N10};
    vec_ONR = {NROperator::ONR0, NROperator::ONR1, NROperator::ONR2, NROperator::ONR3, NROperator::ONR4,
               NROperator::ONR5, NROperator::ONR6, NROperator::ONR7, NROperator::ONR8,
               NROperator::ONR9, NROperator::ONR10, NROperator::ONR11, NROperator::ONR12,
               NROperator::ONR13, NROperator::ONR14, NROperator::ONR15};

    cN.resize(11);
    cN[0].resize(2, std::nan(""));
    for (int i = 1; i < cN.size(); i++)
    {
        cN[i].resize(2, 0);
        cN[i][0] = DMNucleonCoefficient(vec_ON[i], NucleonType::PROTON);
        cN[i][1] = DMNucleonCoefficient(vec_ON[i], NucleonType::NEUTRON);

        //std::cout << i << " " << cN[i][0] << " " << cN[i][1] << std::endl;
    }

    //std::cout << " ----------- " << std::endl;

    vec_ONtoONR.resize(vec_ON.size());
    vec_ONtoONR[0].resize(vec_ONR.size(), std::nan(""));
    for (int i = 1; i < vec_ONtoONR.size(); i++)
    {
        vec_ONtoONR[i].resize(vec_ONR.size());
        vec_ONtoONR[i][0] = std::nan("");

        for (int j = 1; j < vec_ONR.size(); j++)
        {
            vec_ONtoONR[i][j] = ONtoONR(vec_ON[i], vec_ONR[j]);
            //std::cout << i << " " << j << " " << vec_ONtoONR[i][j] << std::endl;
        }
    }

    //std::cout << " ----------- " << std::endl;

    cNR.resize(vec_ONR.size());
    cNR[0].resize(2, std::nan(""));
    for (int i = 1; i < cNR.size(); i++)
    {
        cNR[i].resize(2, 0);
        for (int n = 0; n < 2; n++)
        {
            for (int j = 1; j < vec_ON.size(); j++)
                cNR[i][n] += vec_ONtoONR[j][i] * cN[j][n];
        }

        //std::cout << cNR[i][0] << " " << cNR[i][1] << std::endl;
    }

    cNR_iso.resize(vec_ONR.size());
    cNR_iso[0].resize(2, std::nan(""));
    for (int i = 1; i < cNR_iso.size(); i++)
    {
        cNR_iso[i].resize(2, 0);
        cNR_iso[i][0] = 0.5 * (cNR[i][0] + cNR[i][1]);
        cNR_iso[i][1] = 0.5 * (cNR[i][0] - cNR[i][1]);
    }
}

//

//

double DirectDetection::NROperatorCoefficient_ISO(int index_op, IsospinType isospin)
{
    int index_iso = 0;

    if (isospin == IsospinType::ISOSCALAR)
        index_iso = 0;

    if (isospin == IsospinType::ISOVECTOR)
        index_iso = 1;

    //std::cout << "Index_iso : " << index_iso << std::endl;

    return cNR_iso[index_op][index_iso];
}

double DirectDetection::NROperatorCoefficient_ISO(NROperator op, IsospinType isospin)
{
    int index_op = 0;
    for (int i = 1; i < vec_ONR.size(); i++)
        if (vec_ONR[i] == op)
            index_op = i;

    //std::cout << "Index_op : " << index_op << std::endl;

    return NROperatorCoefficient_ISO(index_op, isospin);
}

double DirectDetection::NROperatorCoefficient_ISO(NROperator op, int index_iso)
{
    int index_op = 0;
    for (int i = 1; i < vec_ONR.size(); i++)
        if (vec_ONR[i] == op)
            index_op = i;

    return NROperatorCoefficient_ISO(index_op, index_iso);
}

double DirectDetection::NROperatorCoefficient(int index_op, NucleonType nucleon)
{
    int index_nucleon = 0;

    if (nucleon == NucleonType::PROTON)
        index_nucleon = 0;

    if (nucleon == NucleonType::NEUTRON)
        index_nucleon = 1;

    return cNR[index_op][index_nucleon];
}

double DirectDetection::NROperatorCoefficient(NROperator op, NucleonType nucleon)
{
    int index_op = 0;
    for (int i = 1; i < vec_ONR.size(); i++)
        if (vec_ONR[i] == op)
            index_op = i;

    return NROperatorCoefficient(index_op, nucleon);
}

double DirectDetection::NROperatorCoefficient(NROperator op, int index_nucleon)
{
    int index_op = 0;
    for (int i = 1; i < vec_ONR.size(); i++)
        if (vec_ONR[i] == op)
            index_op = i;

    return NROperatorCoefficient(index_op, index_nucleon);
}

//

//

double DirectDetection::ONtoONR(EffectiveOperator opeff, NROperator opnr)
{
    double res = 0;

    if (opeff == EffectiveOperator::N1 && opnr == NROperator::ONR1)
        res = 1;
    if (opeff == EffectiveOperator::N2 && opnr == NROperator::ONR11)
        res = -mN / mchi;
    if (opeff == EffectiveOperator::N3 && opnr == NROperator::ONR10)
        res = 1;
    if (opeff == EffectiveOperator::N4 && opnr == NROperator::ONR6)
        res = mN / mchi;
    if (opeff == EffectiveOperator::N5 && opnr == NROperator::ONR1)
        res = 1;
    if (opeff == EffectiveOperator::N6 && opnr == NROperator::ONR8)
        res = 2;
    if (opeff == EffectiveOperator::N6 && opnr == NROperator::ONR9)
        res = 2;
    if (opeff == EffectiveOperator::N7 && opnr == NROperator::ONR8)
        res = -2;
    if (opeff == EffectiveOperator::N7 && opnr == NROperator::ONR9)
        res = 2;
    if (opeff == EffectiveOperator::N8 && opnr == NROperator::ONR4)
        res = -4;
    if (opeff == EffectiveOperator::N9 && opnr == NROperator::ONR4)
        res = 8;
    if (opeff == EffectiveOperator::N10 && opnr == NROperator::ONR11)
        res = 2;
    if (opeff == EffectiveOperator::N10 && opnr == NROperator::ONR10)
        res = -2 * mN / mchi;
    if (opeff == EffectiveOperator::N10 && opnr == NROperator::ONR12)
        res = -1;

    return res;
}

double DirectDetection::DMNucleonCoefficient(EffectiveOperator op, NucleonType nucleon)
{
    double res = 0;
    if (op == EffectiveOperator::N1)
    {
        for (int i = 0; i < 3; i++)
            res += cq[1][i] * mN / mq[i] * fCoeff(vec_Qtypes[i], nucleon);

        for (int j = 3; j < 6; j++)
            res += 2. / 27. * fGcoeff(nucleon) * (cq[1][j] * mN / mq[j] - 1 / (3 * PI) * cg[1] * mN);
    }
    if (op == EffectiveOperator::N2)
    {
        for (int i = 0; i < 3; i++)
            res += cq[2][i] * mN / mq[i] * fCoeff(vec_Qtypes[i], nucleon);

        for (int j = 3; j < 6; j++)
            res += 2. / 27. * fGcoeff(nucleon) * (cq[2][j] * mN / mq[j] - 1 / (3 * PI) * cg[2] * mN);
    }
    if (op == EffectiveOperator::N3)
    {
        for (int i = 0; i < 3; i++)
            res += mN / mq[i] * ((cq[3][i] - C3) + 1. / (2 * PI) * cg[3] * m_bar) * DDCoeff(vec_Qtypes[i], nucleon);
    }
    if (op == EffectiveOperator::N4)
    {
        for (int i = 0; i < 3; i++)
            res += mN / mq[i] * ((cq[4][i] - C4) + 1. / (2 * PI) * cg[4] * m_bar) * DDCoeff(vec_Qtypes[i], nucleon);
    }
    if (op == EffectiveOperator::N5)
    {
        if (nucleon == NucleonType::PROTON)
            res = 2 * cq[5][0] + cq[5][1];

        if (nucleon == NucleonType::NEUTRON)
            res = cq[5][0] + 2 * cq[5][1];
    }
    if (op == EffectiveOperator::N6)
    {
        if (nucleon == NucleonType::PROTON)
            res = 2 * cq[6][0] + cq[6][1];

        if (nucleon == NucleonType::NEUTRON)
            res = cq[6][0] + 2 * cq[6][1];
    }
    if (op == EffectiveOperator::N7)
    {
        for (int q = 0; q < 3; q++)
        {
            res += cq[7][q] * DDCoeff(vec_Qtypes[q], nucleon);
            //std::cout << q << " " << res << std::endl;
        }
    }
    if (op == EffectiveOperator::N8)
    {
        for (int q = 0; q < 3; q++)
            res += cq[8][q] * DDCoeff(vec_Qtypes[q], nucleon);
    }
    if (op == EffectiveOperator::N9)
    {
        for (int q = 0; q < 3; q++)
            res += cq[9][q] * dCoeff(vec_Qtypes[q], nucleon);
    }
    if (op == EffectiveOperator::N10)
    {
        for (int q = 0; q < 3; q++)
            res += cq[10][q] * dCoeff(vec_Qtypes[q], nucleon);
    }

    return res;
}

double DirectDetection::fCoeff(QuarkType quark, NucleonType nucleon)
{
    if (nucleon == NucleonType::PROTON)
    {
        if (quark == QuarkType::UP)
            return f_u_p;
        if (quark == QuarkType::DOWN)
            return f_d_p;
        if (quark == QuarkType::STRANGE)
            return f_s_p;
    }

    if (nucleon == NucleonType::NEUTRON)
    {
        if (quark == QuarkType::UP)
            return f_u_n;
        if (quark == QuarkType::DOWN)
            return f_d_n;
        if (quark == QuarkType::STRANGE)
            return f_s_n;
    }

    return std::nan("");
}

double DirectDetection::DDCoeff(QuarkType quark, NucleonType nucleon)
{
    if (nucleon == NucleonType::PROTON)
    {
        if (quark == QuarkType::UP)
            return DD_u_p;
        if (quark == QuarkType::DOWN)
            return DD_d_p;
        if (quark == QuarkType::STRANGE)
            return DD_s_p;
    }

    if (nucleon == NucleonType::NEUTRON)
    {
        if (quark == QuarkType::UP)
            return DD_u_n;
        if (quark == QuarkType::DOWN)
            return DD_d_n;
        if (quark == QuarkType::STRANGE)
            return DD_s_n;
    }

    return std::nan("");
}

double DirectDetection::dCoeff(QuarkType quark, NucleonType nucleon)
{
    if (nucleon == NucleonType::PROTON)
    {
        if (quark == QuarkType::UP)
            return d_u_p;
        if (quark == QuarkType::DOWN)
            return d_d_p;
        if (quark == QuarkType::STRANGE)
            return d_s_p;
    }

    if (nucleon == NucleonType::NEUTRON)
    {
        if (quark == QuarkType::UP)
            return d_u_n;
        if (quark == QuarkType::DOWN)
            return d_d_n;
        if (quark == QuarkType::STRANGE)
            return d_s_n;
    }

    return std::nan("");
}

double DirectDetection::fGcoeff(NucleonType nucleon)
{
    double res = 0;
    for (int i = 0; i < 3; i++)
        res += fCoeff(vec_Qtypes[i], nucleon);
    return 1 - res;
}
