#ifndef DEF_DIRECTDETECTION
#define DEF_DIRECTDETECTION

/**
 * 
 * In this class we convert the relativistic lagrangian 
 * in the non relativistic version by evaluating the coefficients
 * of the non relativistic operators.
 * 
 * The conversion is made with the normalisation of the spinors as in [1203.3542]
 * This convention is used to match DDCalc relying on [1607.04418] based on [1308.6288], [1504.06554], and [1505.03117]
 * 
 * We start from \f$\mathcal{L}_{\rm int} = \sum_{i} c_i^q \mathcal{O}_i\f^q$ 
 * with \f$ \mathcal{O}_i = \chi\chi \bar q q, \, \chi i \gamma^5 \chi \bar q i\gamma^5 q\f\, \dots$.
 * and we evaluate \f$\mathcal{L}_{\rm int}^{\rm NR} = \sum_{i} c_i^N \mathcal{O}_i^{N} \f^q$
 * with  \f$ \mathcal{O}_i = 1, \, i {\bf s}_N \cdot {\bf s}_\chi\dots \f$.
 * In the end, we can evaluate \f$ c_i^0 \f$ and \f$ c_i^1 \f$ in the isospin base. 
 * Indeed, the code DDCalc v2 [arXiv:1705.07920] relies on a decomposition in the isospin base
 * 
 */

#include "CppLib.h"
#include "MyUnits.h"
#include "mymath.h"
#include "DarkSectorModel.h"

enum class EffectiveOperator
{
    N0, /// No meaning
    N1, /// chi chi N N
    N2, /// chi i gamma_5 chi N N
    N3, /// chi chi N i gamma_5 N
    N4, /// chi i gamma_5 chi N i gamma_5 N
    N5, /// chi gamma_mu chi N gamma_mu N
    N6, /// chi gamma_mu gamma_5 chi N gamma_mu N
    N7, /// chi gamma_mu chi N gamma_mu gamma_5 N
    N8, /// chi gamma_mu gamma_5 chi N gamma_mu gamma_5 N
    N9, /// chi sigma_mu_nu chi N sigma_mu_nu N
    N10 /// chi i sigma_mu_nu gamma_5 chi N sigma_mu_nu N
};

/** Note that conventionnaly ONR2 = (v^\perp)^2 is ommitted as it does not arise from any relativistic theory
 *  Here all the NRoperators are dimensionless */
enum class NROperator
{
    ONR0,  /// No meaning
    ONR1,  /// 1
    ONR2,  /// (v^\perp)^2
    ONR3,  /// i s_N . (q/m_N \times v^\perp)
    ONR4,  /// s_\chi s_N
    ONR5,  /// i s_chi . (q/m_N \times v^\perp)
    ONR6,  /// (s_chi . q/m_N)(s_N . q/m_N)
    ONR7,  /// s_N . v^\perp
    ONR8,  /// s_\chi . v^\perp
    ONR9,  /// i s_chi . (s_N \times q/m_N)
    ONR10, /// i s_N . q/m_N
    ONR11, /// i s_/chi . q/m_N
    ONR12, /// s_\chi . (s_N \times v^\perp) = v^\perp . (s_\chi \times s_N)
    ONR13, /// i (s_\chi . v^\perp) (s_N . q/m_N)
    ONR14, /// i (s_\chi . q/m_N) (s_N . v^\perp)
    ONR15  /// -(s_\chi . q/m_N) [(s_N \times v^\perp) . q/m_N]
};

enum class QuarkType
{
    UP,
    DOWN,
    STRANGE
};

enum class NucleonType
{
    PROTON,
    NEUTRON
};

enum class IsospinType
{
    ISOSCALAR,
    ISOVECTOR
};

/** 
* Note that it is esay to rotate from the Nucleon base to the Isospin base
* by a simple transformation of the coupling coefficients.
* If c_N, c_P are the coefficients in the Nucleon base, then
* c_0(Isoscalar) = 1/2(c_P+c_N) and c_1(Isovector) = 1/2 (c_P - c_N)
* are the coefficients in the Isospin base
*/

class DirectDetection
{

public:
    DirectDetection(DarkSectorModel const &DS);
    ~DirectDetection(){};

    /// Coupling to the NR operator : index_nucleon=0 (proton), index_nucleon=1 (neutron)
    double NROperatorCoefficient(int index_op, int index_nucleon) { return cNR[index_op][index_nucleon]; };
    /// Overriden coupling to the NR operator : index_nucleon=0 (proton), index_nucleon=1 (neutron)
    double NROperatorCoefficient(NROperator op, int index_nucleon);
    /// Overriden coupling to the NR operator
    double NROperatorCoefficient(NROperator op, NucleonType nucleon);
    /// Overriden coupling to the NR operator
    double NROperatorCoefficient(int i, NucleonType nucleon);

    /// Coupling to the NR operator in the isospin basis : index_iso=0 (isoscalar), index_iso = 1 (isovector)
    double NROperatorCoefficient_ISO(int index_op, int index_iso) { return cNR_iso[index_op][index_iso]; };
    /// Overriden coupling to the NR operator in the isospin basis : index_iso=0 (isoscalar), index_iso = 1 (isovector)
    double NROperatorCoefficient_ISO(NROperator op, int index_nucleon);
    /// Overriden coupling to the NR operator in the isospin basis
    double NROperatorCoefficient_ISO(NROperator op, IsospinType isospin);
    /// Overriden coupling to the NR operator in the isospin basis
    double NROperatorCoefficient_ISO(int i, IsospinType isospin);

    /// Makes the conversion between quark couplings and nucleon couplings
    double DMNucleonCoefficient(EffectiveOperator op, NucleonType nucleon);

private:
    /**
     * Members
     */

    /// Coefficient f for the conversion between quark couplings and nucleon couplings
    double fCoeff(QuarkType quark, NucleonType nucleon);
    /// Coefficient Delta for the conversion between quark couplings and nucleon couplings
    double DDCoeff(QuarkType quark, NucleonType nucleon);
    /// Coefficient delta for the conversion between quark couplings and nucleon couplings
    double dCoeff(QuarkType quark, NucleonType nucleon);
    /// Coefficient f_G for the conversion between quark couplings and nucleon couplings
    double fGcoeff(NucleonType nucleon);
    /// Coefficients linking effective operators to the non-relativistic operators <ON_i> = sum_i ONtoONR_(i, j) ONR_j
    double ONtoONR(EffectiveOperator opeff, NROperator opnr);

    /**
     * Attributes
     */

    double mN;   ///< nucleon mass
    double mchi; ///< DM mass

    std::vector<EffectiveOperator> vec_ON; ///< Ordered list of effective operators
    std::vector<NROperator> vec_ONR;       ///< Ordered list of non-relativistic operators
    std::vector<QuarkType> vec_Qtypes;     ///< Ordered list of the light quarks (u, d, s)

    std::vector<std::vector<double> > cNR_iso; ///< coefficient non-relativistic operators (isospin basis)
    std::vector<std::vector<double> > cNR;     ///< coefficient non-relativistic operators
    std::vector<std::vector<double> > cN;      ///< coupling DM-nucleons
    std::vector<std::vector<double> > cq;      ///< coupling DM-quarks
    std::vector<double> cg;                    ///< coupling DM-gluons (at one loop order in our case)

    std::vector<double> mq; ///< quark masses (up, down, strange, charm, bottom, top)

    std::vector<std::vector<double> > vec_ONtoONR;

    /// Parameters for the computation
    double C3, C4, m_bar;
    double f_u_p, f_d_p, f_s_p;
    double f_u_n, f_d_n, f_s_n;
    double DD_u_p, DD_d_p, DD_s_p;
    double DD_u_n, DD_d_n, DD_s_n;
    double d_u_p, d_d_p, d_s_p;
    double d_u_n, d_d_n, d_s_n;
};

#endif