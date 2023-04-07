#include "../headers/mymath.h"

const double PI = 3.1415926535897932384626433832795028841971693993751;
const double pi = PI;
const double PI2 = PI * PI;
const double PI4 = PI2 * PI2;
const double PI8 = PI4 * PI4;
const double PI16 = PI8 * PI8;

const GaussLegendre gl100 = GaussLegendre(100);
const GaussLegendre gl200 = GaussLegendre(200);
const GaussLegendre gl500 = GaussLegendre(500);
const GaussLegendre gl1000 = GaussLegendre(1000);
const GaussLegendre gl2000 = GaussLegendre(2000);
const GaussLegendre gl5000 = GaussLegendre(5000);

/// Compute the modified Bessel function I_0(x) for any real x.
///
/// --- NvE 12-mar-2000 UU-SAP Utrecht
double BesselI0(double x)
{

    // Parameters of the polynomial approximation
    const double p1 = 1.0, p2 = 3.5156229, p3 = 3.0899424, p4 = 1.2067492, p5 = 0.2659732,
                 p6 = 3.60768e-2, p7 = 4.5813e-3;

    const double q1 = 0.39894228, q2 = 1.328592e-2, q3 = 2.25319e-3, q4 = -1.57565e-3,
                 q5 = 9.16281e-3, q6 = -2.057706e-2, q7 = 2.635537e-2, q8 = -1.647633e-2,
                 q9 = 3.92377e-3;

    const double k1 = 3.75;
    double ax = fabs(x);

    double y = 0, result = 0;

    if (ax < k1)
    {
        double xx = x / k1;
        y = xx * xx;
        result = p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7)))));
    }
    else
    {
        y = k1 / ax;
        result = (exp(ax) / sqrt(ax)) * (q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * (q7 + y * (q8 + y * q9))))))));
    }
    return result;
}

double BesselI1(double x)
{
    double ax, ans;
    double y;

    if ((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        ans = (ax * (0.5 + y * (0.87890594 +
                                y * (0.51498869 +
                                     y * (0.15084934 +
                                          y * (0.2658733e-1 +
                                               y * (0.301532e-2 + y * 0.32411e-3)))))));
    }
    else
    {
        y = 3.75 / ax;
        ans = (0.2282967e-1 +
               y * (-0.2895312e-1 +
                    y * (0.1787654e-1 -
                         y * 0.420059e-2)));
        ans = (0.39894228 +
               y * (-0.3988024e-1 +
                    y * (-0.362018e-2 +
                         y * (0.163801e-2 +
                              y * (-0.1031555e-1 + y * ans)))));
        ans *= ((exp(ax) / sqrt(ax)));
    }
    return x < 0.0 ? -ans : ans;
}

///
/// Compute the modified Bessel function K_0(x) for positive real x.
///  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
///     Applied Mathematics Series vol. 55 (1964), Washington.
/// --- NvE 12-mar-2000 UU-SAP Utrecht
///
double BesselK0(double x)
{

    // Parameters of the polynomial approximation
    const double p1 = -0.57721566, p2 = 0.42278420, p3 = 0.23069756,
                 p4 = 3.488590e-2, p5 = 2.62698e-3, p6 = 1.0750e-4, p7 = 7.4e-6;

    const double q1 = 1.25331414, q2 = -7.832358e-2, q3 = 2.189568e-2,
                 q4 = -1.062446e-2, q5 = 5.87872e-3, q6 = -2.51540e-3, q7 = 5.3208e-4;

    if (x <= 0)
    {
        //printf("BesselK0 Invalid argument x = %g\n",x);
        return 0;
    }

    double y = 0, result = 0;

    if (x <= 2)
    {
        y = x * x / 4;
        result = (-log(x / 2.) * BesselI0(x)) +
                 (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7))))));
    }
    else
    {
        y = 2 / x;
        result = (exp(-x) / sqrt(x)) * (q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * q7))))));
    }
    return result;
}

double BesselK1(double x)
{
    double y, ans;

    if (x <= 2.0)
    {
        y = x * x / 4.0;
        ans = (log(x / 2.0) * BesselI1(x)) +
              (1.0 / x) * (1.0 + y * (0.15443144 +
                                      y * (-0.67278579 +
                                           y * (-0.18156897 +
                                                y * (-0.1919402e-1 +
                                                     y * (-0.110404e-2 +
                                                          y * (-0.4686e-4)))))));
    }
    else
    {
        y = 2.0 / x;
        ans = (exp(-x) / sqrt(x)) * (1.25331414 +
                                     y * (0.23498619 +
                                          y * (-0.3655620e-1 +
                                               y * (0.1504268e-1 +
                                                    y * (-0.780353e-2 +
                                                         y * (0.325614e-2 + y * (-0.68245e-3)))))));
    }
    return ans;
}

double BesselK2(double x)
{
    return BesselKn(2, x);
}

double BesselK3(double x)
{
    return BesselKn(3, x);
}

///
/// Compute the Integer Order Modified Bessel function K_n(x)
/// for n=0,1,2,... and positive real x.
///
double BesselKn(int n, double x)
{
    if (x <= 0 || n < 0)
        return 0;

    //if(x > 299)
    //   return exp(-x)*sqrt(PI/2)*pow(x, -0.5)*(1 + 1./8.*(-1+4*n*n)/x);

    if (n == 0)
        return BesselK0(x);
    if (n == 1)
        return BesselK1(x);

    // Perform upward recurrence for all x
    double tox = 2 / x;
    double bkm = BesselK0(x);
    double bk = BesselK1(x);
    double bkp = 0;
    for (int j = 1; j < n; j++)
    {
        bkp = bkm + (double)j * tox * bk;
        bkm = bk;
        bk = bkp;
    }
    return bk;
}

///
/// First derivative of Kn => dKn(x)/dx.
///
double DBesselKn(int n, double x)
{
    if (x <= 0 || n < 1)
        return 0;
    return 0.5 * (BesselKn(n + 1, x) - BesselKn(n - 1, x));
}

///
/// Modif Bessel function K2 => multiplied by exp(x) to avoid
/// crashes when x gets high.
/// Compute the modified Bessel function K_0(x) for positive real x.
///  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
///     Applied Mathematics Series vol. 55 (1964), Washington.
/// --- NvE 12-mar-2000 UU-SAP Utrecht.
///
double expBesselK0(double x)
{

    // Parameters of the polynomial approximation
    const double p1 = -0.57721566, p2 = 0.42278420, p3 = 0.23069756,
                 p4 = 3.488590e-2, p5 = 2.62698e-3, p6 = 1.0750e-4, p7 = 7.4e-6;

    const double q1 = 1.25331414, q2 = -7.832358e-2, q3 = 2.189568e-2,
                 q4 = -1.062446e-2, q5 = 5.87872e-3, q6 = -2.51540e-3, q7 = 5.3208e-4;

    if (x <= 0)
    {
        printf("BesselK0 Invalid argument x = %g\n", x);
        return 0;
    }

    double y = 0, result = 0;

    if (x <= 2)
    {
        y = x * x / 4;
        result = exp(x) * ((-log(x / 2.) * BesselI0(x)) +
                           (p1 + y * (p2 + y * (p3 + y * (p4 + y * (p5 + y * (p6 + y * p7)))))));
    }
    else
    {
        y = 2 / x;
        result = (1. / sqrt(x)) * (q1 + y * (q2 + y * (q3 + y * (q4 + y * (q5 + y * (q6 + y * q7))))));
    }
    return result;
}

///
/// Modif Bessel function K2 => multiplied by exp(x) to avoid
/// crashes when x gets high.
///
double expBesselK1(double x)
{
    double y, ans;

    if (x <= 2.0)
    {
        y = x * x / 4.0;
        ans = exp(x) * ((log(x / 2.0) * BesselI1(x)) +
                        (1.0 / x) * (1.0 +
                                     y * (0.15443144 +
                                          y * (-0.67278579 +
                                               y * (-0.18156897 +
                                                    y * (-0.1919402e-1 +
                                                         y * (-0.110404e-2 +
                                                              y * (-0.4686e-4))))))));
    }
    else
    {
        y = 2.0 / x;
        ans = (1. / sqrt(x)) * (1.25331414 +
                                y * (0.23498619 +
                                     y * (-0.3655620e-1 +
                                          y * (0.1504268e-1 +
                                               y * (-0.780353e-2 +
                                                    y * (0.325614e-2 + y * (-0.68245e-3)))))));
    }
    return ans;
}

///
/// Modif Bessel function K2 => multiplied by exp(x) to avoid
/// crashes when x gets high.
///
double expBesselK2(double x)
{
    double bk, bkm, bkp, tox;

    tox = 2.0 / x;
    bkm = expBesselK0(x);
    bk = expBesselK1(x);
    bkp = bkm + tox * bk;

    bkm = bk;
    bk = bkp;
    return bk;
}

double polyLog(unsigned int n, double x)
{
    if (abs(x) > 1)
    {
        std::cout << "FATAL ERROR : the code cannot compute polylog functions for |x| > 1" << std::endl;
        exit(0);
    }

    if(x == 1)
        return PI*PI/6.; 

    double err, old_term, new_term, sum = 0;
    int k = 1;

    old_term = x / k;
    sum = old_term;

    //std::cout << 1 << " " << old_term << std::endl;
    do
    {
        k++;
        new_term = pow(x, k) * pow(1. * k, -1. * n);
        err = abs(new_term) / abs(sum);
        sum += new_term;
        old_term = new_term;
        //std::cout << err << " " << k  << " " << sum << std::endl;
        //std::cout << k << " " << old_term << " " << std::pow(1.*k, -1.*n) <<  std::endl;
    } while (err > 1e-7 && k < 5000);

    if (err > 1e-7)
    {
        std::cout << "FATAL ERROR : polylog function has not reached wanted precision" << std::endl;
        exit(0);
    }

    return sum;
}

// Regulator function for the QCD phase transition 
double regulatorTQCD(double x, double alpha, double xshift, double max, double min)
{
  return ((tanh(alpha * (x - xshift)) + 1) / 2) * (max - min) + min;
}

dcomp fOneLoopPseudoScalarVectors(double x)
{
    dcomp _x = x;
    dcomp res = _x / 4.0 * pow(log((sqrt(1.0 - _x) - 1.0) / (sqrt(1.0 - _x) + 1.0)), 2);

    return res;
}

dcomp fOneLoopScalarVectors(double x)
{
    dcomp _x = x;
    dcomp res = _x * ((_x - 1.0) / 4.0 * pow(log((sqrt(1.0 - _x) - 1.0) / (sqrt(1.0 - _x) + 1.0)), 2) + 1.0);

    return res;
}

///
/// Simpson method that computes the integral of myfunc(x) between xmin
/// and xmax. This time, x is a vector and the integral is performed over the
/// ivar_th coordinate. The steps are linear/log/exp/log2 according to ilog
/// (0/1/2/3). Nsteps is the number of steps, and err the approximate relative
/// error.
/// For the lin steps, we integrate along:
/// \f[ {\cal I} = \int_{a}^{b} dx\, f(x)\;; \f]
/// for the log steps:
/// \f[ {\cal I} = \int_{\ln(a)}^{\ln(b)} d\ln(x)\,(x\times f(x))\;; \f]
/// for the exp steps:
/// \f[ {\cal I} = \int_{e^{a/b}}^{e^{1}} d\exp(x)\,(b\,e^{-x/b}\, f(x))\;; \f]
/// for the double-log steps:
/// \f[ {\cal I} = \int_{\ln(\ln(a))}^{\ln(\ln(b))} d\ln(\ln(x))\, x\ln(x)\times f(x)\;.\f]
///

double Simpson_Integral1(int ivar, int isize, int ilog, int Nsteps,
                         double xmin, double xmax, std::vector<double> xx,
                         double (*myfunc)(std::vector<double>), double &err)
{
    using namespace std;
    string funname = "Simpson_Integral1";
    ostringstream oss;

    int i;

    int ndim = xx.size();
    if (isize != ndim)
    {
        cout << "Check vector length!!!\n";
        exit(1);
    }
    if (ivar > (ndim - 1))
    {
        cout << "Check vector length!!!\n";
        exit(1);
    }

    // We do not modify the initial vector
    vector<double> x = xx;

    int Nmax = 100000, Nmin = 10;

    oss << "Nsteps = " << Nsteps << " is too large: it has been decreased down to N =" << Nmax;
    if (Nsteps > Nmax)
    {
        cout << (oss.str()).c_str() << endl;
        exit(1);
    }

    double epsilon = 1.e-3;

    Nsteps = myMax(Nmin, Nsteps);
    Nsteps = myMin(Nmax, Nsteps);

    double epsx = 1.e-1 / (double)Nsteps;

    // Nsteps has to be even for the Simpson method to be used.
    Nsteps += (Nsteps % 2);

    // linear steps
    double dx = (xmax - xmin) / (double)Nsteps;
    // log steps
    double dlx = 0, dlxcorr = 0;
    if (myAbs(xmin) > 0 && myAbs(xmax) > 0 && (xmax / xmin) > 0)
    {
        dlx = log(xmax / xmin) / (double)Nsteps;
    }
    else if (ilog == 1)
    {
        // if xmin = 0 or xmax = 0 or xmax/xmin<0, we make a change y = x + dlxcorr such that ymin>0
        // so that we can define log steps for y.
        //if(xmin<=0 && xmax>xmin)
        if (xmin <= 0 || xmax <= 0)
        {
            dlxcorr = (xmax > xmin) ? -xmin + 0.1 : -xmax + 0.1;
            dlx = log((xmax + dlxcorr) / (xmin + dlxcorr)) / (double)Nsteps;
        }
        else
        {
            cout << "Log steps ---> Lin steps FORCED!!!!" << endl;
            ilog = 0;
        }
    }
    // exp steps
    double ex = 0, dexmin = exp(xmin / xmax), dexmax = exp(1.);
    double dex = (dexmax - dexmin) / (double)Nsteps;
    // inverted double-log steps
    double delta = log10(xmin) - 1e-10, l2x = 0;
    double l2xmin = log10(log10(xmin) - delta), l2xmax = log10(log10(xmax) - delta);
    double dl2x = (l2xmax - l2xmin) / (double)Nsteps;
    // double-exp steps
    double e2x = 0, de2xmin = exp(exp(xmin / xmax)), de2xmax = exp(exp(1.));
    double de2x = (de2xmax - de2xmin) / (double)Nsteps;

    double res = 0, resref = res, cx = 0;
    switch (ilog)
    {
    case 0: // linear steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            x[ivar] = xmin + i * dx;
            res += cx * dx * (*myfunc)(x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
        break;
    case 1: // log steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            x[ivar] = (xmin + dlxcorr) * exp(i * dlx) - dlxcorr;
            dx = dlx * (x[ivar] + dlxcorr);
            res += cx * dx * (*myfunc)(x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
        break;
    case 2: // exp steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            ex = dexmin + i * dex;
            x[ivar] = xmax * log(ex);
            dx = xmax * dex / ex;
            res += cx * dx * (*myfunc)(x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
        break;
    case 3: // inverted double-log steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            l2x = l2xmin + i * dl2x;
            x[ivar] = xmin + xmax - pow(10, delta + pow(10, l2xmin + (Nsteps - i) * dl2x));
            dx = log(10) * log(10) * dl2x * pow(10, l2x) * pow(10, delta + pow(10, l2x));
            res += cx * dx * (*myfunc)(x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
    case 4: // double-exp steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            e2x = de2xmin + i * de2x;
            x[ivar] = xmax * log(log(e2x));
            dx = xmax * de2x / (e2x * log(e2x));
            res += cx * dx * (*myfunc)(x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
        break;
    default:
        cout << "undefined integration method!" << endl;
        break;
    }
    err = myAbs(err);

    return res;
}

double Simpson_Integral1_Static(int ivar, int isize, int ilog, int Nsteps,
                                   double xmin, double xmax, std::vector<double> xx, void *pt2Object,
                                   double (*myfunc)(void *pt2Object, std::vector<double>), double &err)
{
    using namespace std;
    string funname = "Simpson_Integral1";
    ostringstream oss;

    int i;

    int ndim = xx.size();
    if (isize != ndim)
    {
        cout << "Check vector length!!!\n";
        exit(1);
    }
    if (ivar > (ndim - 1))
    {
        cout << "Check vector length!!!\n";
        exit(1);
    }

    // We do not modify the initial vector
    vector<double> x = xx;

    int Nmax = 100000, Nmin = 10;

    oss << "Nsteps = " << Nsteps << " is too large: it has been decreased down to N =" << Nmax;
    if (Nsteps > Nmax)
    {
        cout << (oss.str()).c_str() << endl;
        exit(1);
    }

    double epsilon = 1.e-3;

    Nsteps = myMax(Nmin, Nsteps);
    Nsteps = myMin(Nmax, Nsteps);

    double epsx = 1.e-1 / (double)Nsteps;

    // Nsteps has to be even for the Simpson method to be used.
    Nsteps += (Nsteps % 2);

    // linear steps
    double dx = (xmax - xmin) / (double)Nsteps;
    // log steps
    double dlx = 0, dlxcorr = 0;
    if (myAbs(xmin) > 0 && myAbs(xmax) > 0 && (xmax / xmin) > 0)
    {
        dlx = log(xmax / xmin) / (double)Nsteps;
    }
    else if (ilog == 1)
    {
        // if xmin = 0 or xmax = 0 or xmax/xmin<0, we make a change y = x + dlxcorr such that ymin>0
        // so that we can define log steps for y.
        //if(xmin<=0 && xmax>xmin)
        if (xmin <= 0 || xmax <= 0)
        {
            dlxcorr = (xmax > xmin) ? -xmin + 0.1 : -xmax + 0.1;
            dlx = log((xmax + dlxcorr) / (xmin + dlxcorr)) / (double)Nsteps;
        }
        else
        {
            cout << "Log steps ---> Lin steps FORCED!!!!" << endl;
            ilog = 0;
        }
    }
    // exp steps
    double ex = 0, dexmin = exp(xmin / xmax), dexmax = exp(1.);
    double dex = (dexmax - dexmin) / (double)Nsteps;
    // inverted double-log steps
    double delta = log10(xmin) - 1e-10, l2x = 0;
    double l2xmin = log10(log10(xmin) - delta), l2xmax = log10(log10(xmax) - delta);
    double dl2x = (l2xmax - l2xmin) / (double)Nsteps;
    // double-exp steps
    double e2x = 0, de2xmin = exp(exp(xmin / xmax)), de2xmax = exp(exp(1.));
    double de2x = (de2xmax - de2xmin) / (double)Nsteps;

    double res = 0, resref = res, cx = 0;
    switch (ilog)
    {
    case 0: // linear steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            x[ivar] = xmin + i * dx;
            res += cx * dx * (*myfunc)(pt2Object, x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
        break;
    case 1: // log steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            x[ivar] = (xmin + dlxcorr) * exp(i * dlx) - dlxcorr;
            dx = dlx * (x[ivar] + dlxcorr);
            res += cx * dx * (*myfunc)(pt2Object, x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
        break;
    case 2: // exp steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            ex = dexmin + i * dex;
            x[ivar] = xmax * log(ex);
            dx = xmax * dex / ex;
            res += cx * dx * (*myfunc)(pt2Object, x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
        break;
    case 3: // inverted double-log steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            l2x = l2xmin + i * dl2x;
            x[ivar] = xmin + xmax - pow(10, delta + pow(10, l2xmin + (Nsteps - i) * dl2x));
            dx = log(10) * log(10) * dl2x * pow(10, l2x) * pow(10, delta + pow(10, l2x));
            res += cx * dx * (*myfunc)(pt2Object, x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
    case 4: // double-exp steps
        for (i = 0; i <= Nsteps; ++i)
        {
            resref = res;
            cx = simpson_coef(i, Nsteps);
            e2x = de2xmin + i * de2x;
            x[ivar] = xmax * log(log(e2x));
            dx = xmax * de2x / (e2x * log(e2x));
            res += cx * dx * (*myfunc)(pt2Object, x);
            if ((i > (Nsteps / 4)) && res != 0.)
                err = 1. - resref / res;
        }
        break;
    default:
        cout << "undefined integration method!" << endl;
        break;
    }
    err = myAbs(err);

    return res;
}

double GaussLegendre_Integral_Static(int const ivar, int const isize, GaussLegendre const &gl,
                                     double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
                                     double (*myfunc)(void *pt2Object, std::vector<double>))
{

    if (isize != xx.size())
    {
        std::cout << "FATAL ERROR : isize does not agree -> " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    int ngl = gl.get_ngl();
    double res = 0;

    std::vector<double> _xx = xx;

    double A = (xmax - xmin) / 2;
    double B = (xmax + xmin) / 2;

    for (int i = 1; i < ngl + 1; i++)
    {
        _xx[ivar] = A * gl.get_x(i) + B;
        res += gl.get_w(i) * myfunc(pt2Object, _xx);
    }

    res *= A;

    return res;
}

double GaussLegendre_IntegralExp_Static(int const ivar, int const isize, GaussLegendre const &gl,
                                        double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
                                        double (*myfunc)(void *pt2Object, std::vector<double>))
{
    if (isize != xx.size())
    {
        std::cout << "FATAL ERROR : isize does not agree -> " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    int ngl = gl.get_ngl();
    double res = 0;

    std::vector<double> _xx = xx;

    double A = (exp(xmax) - exp(xmin)) / 2;
    double B = (exp(xmax) + exp(xmin)) / 2;
    double C = 0;

    for (int i = 1; i < ngl + 1; i++)
    {
        C = A * gl.get_x(i) + B;
        _xx[ivar] = log(C);
        res += gl.get_w(i) * myfunc(pt2Object, _xx) / C;
    }

    res *= A;

    return res;
}

double GaussLegendre_IntegralLn_Static(int const ivar, int const isize, GaussLegendre const &gl,
                                       double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
                                       double (*myfunc)(void *pt2Object, std::vector<double>))
{
    if (isize != xx.size())
    {
        std::cout << "FATAL ERROR : isize does not agree -> " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    int ngl = gl.get_ngl();
    double res = 0;

    double dxcorr = 0;
    double eps = 1;

    std::vector<double> _xx = xx;

    if (xmin * xmax <= 0 || (xmin * xmax > 0 && xmin < 0))
    {
        dxcorr = (xmax > xmin) ? -xmin + 0.1 : -xmax + 0.1;
        xmax += dxcorr;
        xmin += dxcorr;
    }
    /*else if()
    {
        //eps = -1;
        dxcorr = (xmax > xmin) ? -xmin + 0.1 : -xmax + 0.1;
        xmax += dxcorr;
        xmin += dxcorr;
    }*/

    double A = eps * log(xmax / xmin) / 2;
    double B = log(xmax * xmin) / 2;
    double C = 0;

    for (int i = 1; i < ngl + 1; i++)
    {
        C = A * gl.get_x(i) + B;
        _xx[ivar] = eps * exp(C) - dxcorr;
        res += gl.get_w(i) * myfunc(pt2Object, _xx) * exp(C);
    }

    res *= A;

    return res;
}

double GaussLegendre_IntegralInverseLnLn_Static(int const ivar, int const isize, GaussLegendre const &gl,
                                                double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
                                                double (*myfunc)(void *pt2Object, std::vector<double>))
{

    if (isize != xx.size())
    {
        std::cout << "FATAL ERROR : isize does not agree -> " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    int ngl = gl.get_ngl();
    double res = 0;

    double dxcorr = 0;

    std::vector<double> _xx = xx;

    if (xmin * xmax <= 0 || (xmin * xmax > 0 && xmin < 1))
    {
        dxcorr = (xmax > xmin) ? -xmin + 1 + 1e-10 : -xmax + 1 + 1e-10;
        xmax += dxcorr;
        xmin += dxcorr;
    }
    /*else if()
    {
        //eps = -1;
        dxcorr = (xmax > xmin) ? -xmin + 0.1 : -xmax + 0.1;
        xmax += dxcorr;
        xmin += dxcorr;
    }*/

    double A = log(log(xmax) / log(xmin)) / 2;
    double B = log(log(xmax) * log(xmin)) / 2;
    double C = 0;

    for (int i = 1; i < ngl + 1; i++)
    {
        C = A * gl.get_x(i) + B;
        _xx[ivar] = -exp(exp(C)) - dxcorr + xmin + xmax;
        res += gl.get_w(i) * myfunc(pt2Object, _xx) * exp(C) * exp(exp(C));
    }

    res *= A;

    return res;
}

double GaussLegendre_IntegralInverseLn_Static(int const ivar, int const isize, GaussLegendre const &gl,
                                              double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
                                              double (*myfunc)(void *pt2Object, std::vector<double>))
{

    if (isize != xx.size())
    {
        std::cout << "FATAL ERROR : isize does not agree -> " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    int ngl = gl.get_ngl();
    double res = 0;

    double dxcorr = 0;

    std::vector<double> _xx = xx;

    if (xmin * xmax <= 0 || (xmin * xmax > 0 && xmin < 1))
    {
        dxcorr = (xmax > xmin) ? -xmin + 1 + 1e-4 : -xmax + 1 + 1e-4;
        xmax += dxcorr;
        xmin += dxcorr;
    }
    /*else if()
    {
        //eps = -1;
        dxcorr = (xmax > xmin) ? -xmin + 0.1 : -xmax + 0.1;
        xmax += dxcorr;
        xmin += dxcorr;
    }*/

    double A = log(xmax / xmin) / 2;
    double B = log(xmax * xmin) / 2;
    double C = 0;

    for (int i = 1; i < ngl + 1; i++)
    {
        C = A * gl.get_x(i) + B;
        _xx[ivar] = -exp(C) - dxcorr + xmin + xmax;
        res += gl.get_w(i) * myfunc(pt2Object, _xx) * exp(C);
    }

    res *= A;

    return res;
}

double GaussLegendre_IntegralLnLn_Static(int const ivar, int const isize, GaussLegendre const &gl,
                                         double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
                                         double (*myfunc)(void *pt2Object, std::vector<double>))
{

    if (isize != xx.size())
    {
        std::cout << "FATAL ERROR : isize does not agree -> " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    int ngl = gl.get_ngl();
    double res = 0;

    double dxcorr = 0;

    std::vector<double> _xx = xx;

    if (xmin * xmax <= 0 || (xmin * xmax > 0 && xmin < 1))
    {
        dxcorr = (xmax > xmin) ? -xmin + 1 + 1e-10 : -xmax + 1 + 1e-10;
        xmax += dxcorr;
        xmin += dxcorr;
    }
    /*else if()
    {
        //eps = -1;
        dxcorr = (xmax > xmin) ? -xmin + 0.1 : -xmax + 0.1;
        xmax += dxcorr;
        xmin += dxcorr;
    }*/

    double A = log(log(xmax) / log(xmin)) / 2;
    double B = log(log(xmax) * log(xmin)) / 2;
    double C = 0;

    for (int i = 1; i < ngl + 1; i++)
    {
        C = A * gl.get_x(i) + B;
        _xx[ivar] = exp(exp(C)) - dxcorr;
        res += gl.get_w(i) * myfunc(pt2Object, _xx) * exp(C) * exp(exp(C));
    }

    res *= A;

    return res;
}

double GSL_Integral_Static_QAGS(void *params, gsl_integration_workspace *w, size_t limit, double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
                                double (*myfunc)(double, void *), double errabs, double errrel, double &err)
{
    gsl_function F;
    F.params = params;
    F.function = myfunc;

    double integral;

    gsl_integration_qags(&F, xmin, xmax, errabs, errrel, limit, w, &integral, &err);

    return integral;
}

// ///
// /// Simpson method that computes the integral of myfunc(x) between xmin
// /// and xmax. This time, x is a vector and the integral is performed over the
// /// ivar_th coordinate. The steps are linear/log/exp according to ilog
// /// (0/1/2). Nsteps is the number of steps, and err the approximate relative
// /// error.
// /// For the lin steps, we integrate along:
// /// \f[ {\cal I} = \int_{a}^{b} dx\, f(x)\;; \f]
// /// for the log steps:
// /// \f[ {\cal I} = \int_{\ln(a)}^{\ln(b)} d\ln(x)\,(x\times f(x))\;; \f]
// /// and finally, for the exp steps:
// /// \f[ {\cal I} = \int_{e^{a/b}}^{e^{1}} d\exp(x)\,(b\,e^{-x/b}\, f(x))\;. \f]
// ///
// double Simpson_Integral3(int ivar,int isize,int ilog,int Nsteps,
// 			 double xmin,double xmax,std::vector<double> xx,
// 			 double (*myfunc)(std::vector<double>),double &err)
// {
//   using namespace std;
//   vector<double> xchange,res{0};
//   int ichange = -1;
//   Simpson_Integral3(ivar,isize,ilog,Nsteps,
// 		    xmin,xmax,xx,ichange,xchange,
// 		    myfunc,res,err);
//   return res[0];
// }

// ///
// /// Simpson method that computes the integral of myfunc(x) between xmin
// /// and xmax. This time, x is a vector and the integral is performed over the
// /// ivar_th coordinate. The steps are linear/log/exp according to ilog
// /// (0/1/2). Nsteps is the number of steps, and err the approximate relative
// /// error.
// /// For the lin steps, we integrate along:
// /// \f[ {\cal I} = \int_{a}^{b} dx\, f(x)\;; \f]
// /// for the log steps:
// /// \f[ {\cal I} = \int_{\ln(a)}^{\ln(b)} d\ln(x)\,(x\times f(x))\;; \f]
// /// and finally, for the exp steps:
// /// \f[ {\cal I} = \int_{e^{a/b}}^{e^{1}} d\exp(x)\,(b\,e^{-x/b}\, f(x))\;. \f]
// ///
// void Simpson_Integral3(int ivar,// variable id
// 		       int isize,// size of vector of input parameters
// 		       int ilog,// defines method (lin/log/exp steps)
// 		       int Nsteps,
// 		       double xmin,double xmax,// boundaries
// 		       const std::vector<double> &xx,// contains input parameters + variable slot
// 		       int ichange,// parameter that may change -> indicate position
// 		       const std::vector<double> &xchange,// contains the changes
// 		       double (*myfunc)(std::vector<double>),// function to integrate over
// 		       std::vector<double> &res,// vector of results, same size as xchange
// 		       double &err)// numerical error estimate
// {
//   using namespace std;
//   string funname = "Simpson_Integral1";
//   ostringstream oss;

//   int i,j,Nchange=xchange.size();
//   if(ichange<0) Nchange = 1;
//   res.clear(); res.resize(Nchange,0.0);

//   int ndim = xx.size();
//   if((isize!=ndim) || (ivar > (ndim-1))) Error_mess(2,funname,"Check vector length!!!\n");
//   else if(ivar == ichange) Error_mess(2,funname,"Parameter id = variable id!!!\n");

//   // We do not modify the initial vector
//   vector<double> x = xx;

//   double cx,resref=0;
//   double dx,ex,alpha,x0,_xmax,dexmin,dexmax,Jacobian;
//   int nexp = 0;

//   int Nmax = 50000;
//   int Nmin = 10;

//   oss << "Nsteps = " << Nsteps << " is too large: it has been decreased down to N =" << Nmax;
//   if(Nsteps<Nmin) Nsteps = Nmin;
//   else if(Nsteps>Nmax){ Nsteps = Nmax; Error_mess(0,funname,oss.str()); }

//   // Nsteps has to be even for the Simpson method to be used.
//   Nsteps += (Nsteps%2);

//   // linear steps => check whether the interval is well ordered.
//   dx = (xmax-xmin)/(double)Nsteps;
//   if(dx<=0.0) return;

//   switch(ilog){

//   // inverted super-log steps (Thomas Lacroix)
// case -4:
//   // This method is aimed at providing an alternative to the super-exponential method which
//   // may lead to numerical overflow.
//   // We make two changes of variables:
//   // (i) xbar = xmin+xmax-x
//   // (ii) y = log(log(alpha*(xbar+x0))), such that x0 = max(1,-xmin+1)
//   // and alpha = 2/(xmin+x0) ensures that ymin = log(2).
//   // The Jacobian reads dxbar/dy = (xbar+x0) log(alpha*(xbar+x0))
//   x0 = myMax(1.0,1.-xmin);
//   alpha = 2./(xmin+x0);
//   dexmin = log(log(alpha*(xmin+x0)));
//   dexmax = log(log(alpha*(xmax+x0)));
//   dx = log(log(alpha*(xmax+x0))/log(alpha*(xmin+x0)))/Nsteps;

//   for(i=0;i<=Nsteps;++i)
//     {
// 	resref = res[0];
// 	cx = simpson_coef(i,Nsteps);
// 	ex = dexmin+i*dx;
// 	ex = exp(ex);
// 	Jacobian = ex*exp(ex)/alpha;
// 	x[ivar] = xmin+xmax - (Jacobian/ex - x0);

// 	switch(ichange)
// 	  {
// 	  case -1:
// 	    res[0] += cx * Jacobian * (* myfunc)(x);
// 	    break;
// 	  default:
// 	    for(j=0;j<Nchange;++j)
// 	      {
// 		x[ichange] = xchange[j];
// 		res[j] += cx * Jacobian * (* myfunc)(x);
// 	      }
// 	    break;
// 	  }
// 	err = 1. - resref/res[0] ;
//     }
//   break;

//     // inverse log steps (not really accurate)
//   case -3:
//     // Here we use the change of variable y=log(alpha/(x+x0)), such that x0=xmin gives an order
//     // of magnitude for the considered numbers, and alpha = 2(xmax+x0) allows y>log(2). This is
//     // similar to the exp change of variable in spirit (see case 2).
//     // We then have x = alpha*exp(-y)-x0 and the Jacobian dx/dy = - (x+x0).
//     // We can then scan the integration range from ymin = y(xmax) to ymax = y(xmin)
//     // => we change the sign of the Jacobian accordingly.
//     x0 = myMax(1.e-2,xmin);
//     alpha = 2*(xmax+x0);
//     dexmin = log(alpha/(xmax+x0));
//     dexmax = log(alpha/(xmin+x0));
//     dx = log((xmax+x0)/(xmin+x0))/Nsteps;
//     // Above, dx is defined such that dx>0. We integrate from ymin=y(xmax) up to ymax=y(xmin).

//     for(i=0;i<=Nsteps;++i){
//       resref = res[0];
//       cx = simpson_coef(i,Nsteps);
//       ex = dexmin+i*dx;
//       // increment on the new variable y = log(1/(x+x0)). Note that exp(-y) = x+x0.
//       Jacobian = alpha*exp(-ex);
//       x[ivar] = Jacobian - x0;

//       switch(ichange){
//       case -1:
// 	res[0] += cx * Jacobian * (* myfunc)(x);
// 	break;
//       default:
// 	for(j=0;j<Nchange;++j){
// 	  x[ichange] = xchange[j];
// 	  res[j] += cx * Jacobian * (* myfunc)(x);
// 	}
// 	break;
//       }
//       //if(ichange<0) res[0] += cx * Jacobian * (* myfunc)(x);
//       //else for(j=0;j<Nchange;++j){
//       //	  x[ichange] = xchange[j];
//       //	  res[j] += cx * Jacobian * (* myfunc)(x);
//       //	}
//       //cout << xmin << " " << x[ivar] << " " << xmax << " " << dx << " " << res[0]*dx << endl;
//       err = 1. - resref/res[0];
//     }
//     break;

//     // inverse log steps (not really accurate)
//   case -2:
//     // Here we use the change of variable y=log(xmax-x+x0) which may apply to any range such that
//     // x0>x-xmax. This is similar to the exp change of variable in spirit (see case 2).
//     // We then have x = xmax+x0-exp(y) and the Jacobian dx/dy = - (xmax+x0-x).
//     // We can then scan the integration range from ymin = log(x0) to
//     // ymax = log(xmax-xmin+x0) => we change the sign of the Jacobian accordingly.
//     x0 = myMax(10.,xmin);
//     dexmin = log(x0);
//     dexmax = log(xmax-xmin+x0);
//     dx = log((xmax-xmin+x0)/(x0))/Nsteps;
//     // Above, dx is defined such that dx>0. We integrate from ymin=y(xmax) up to ymax=y(xmin).

//     for(i=0;i<=Nsteps;++i){
//       resref = res[0];
//       cx = simpson_coef(i,Nsteps);
//       ex = dexmin+i*dx;// we start from ymin = y(xmax)
//       // increment on the new variable y = log(xmax-x+x0). Note that exp(y) = xmax-x+x0.
//       Jacobian = exp(ex);
//       x[ivar] = xmax + x0 - Jacobian;

//       switch(ichange){
//       case -1:
// 	res[0] += cx * Jacobian * (* myfunc)(x);
// 	break;
//       default:
// 	for(j=0;j<Nchange;++j){
// 	  x[ichange] = xchange[j];
// 	  res[j] += cx * Jacobian * (* myfunc)(x);
// 	}
// 	break;
//       }
//       //if(ichange<0)
//       //	res[0] += cx * Jacobian * (* myfunc)(x);
//       //else for(j=0;j<Nchange;++j){
//       //	  x[ichange] = xchange[j];
//       //	  res[j] += cx * Jacobian * (* myfunc)(x);
//       //	}
//       ////cout << xmin << " " << x[ivar] << " " << xmax << " " << dx << " " << res[0]*dx << endl;
//       err = 1. - resref/res[0];
//     }
//     break;

//     // super-log steps
//   case -1:
//     // Here we use the change of variable y=log(log(x+x0)) which may apply to any x>-x0+exp(delta).
//     // We then have x = exp(exp(y))-x0 and the Jacobian dx/dy = (x+x0) log(x+x0).
//     // We can then scan the integration range from ymin = log(1/(xmax+x0)) to
//     // ymax = log(1/(xmin+x0)).
//     x0 = myMax(exp(1.e-5)-xmin,exp(1.e-5));
//     dexmax = log(log(xmax+x0));
//     dexmin = log(log(xmin+x0));
//     dx = log(log(xmax+x0)/log(xmin+x0))/(double)Nsteps;

//     for(i=0;i<=Nsteps;++i){
//       resref = res[0];
//       cx = simpson_coef(i,Nsteps);
//       ex = dexmin+i*dx;
//       Jacobian = exp(exp(ex));// x+x0 = exp(exp(y))
//       x[ivar] = Jacobian - x0;
//       Jacobian = Jacobian * log(Jacobian);// Jacobian = (x+x0) log(x+x0)

//       switch(ichange){
//       case -1:
// 	res[0] += cx * Jacobian * (* myfunc)(x);
// 	break;
//       default:
// 	for(j=0;j<Nchange;++j){
// 	  x[ichange] = xchange[j];
// 	  res[j] += cx * Jacobian * (* myfunc)(x);
// 	}
// 	break;
//       }
//       // Just above, we have considered a change of variable y=log(log(x+corr)). The Jacobian for
//       // this change is dx/dy = (x+corr) log(x+corr), so that dx f(x) -> dy |dx/dy| f(x).
//       err = 1. - resref/res[0] ;
//     }
//     break;

//     // linear steps (accurate if relevant)
//   case 0:
//     dx = (xmax-xmin)/(double)Nsteps;
//     for(i=0;i<=Nsteps;++i){
//       resref = res[0];
//       cx = simpson_coef(i,Nsteps);
//       x[ivar] = xmin + i*dx;

//       switch(ichange){
//       case -1:
// 	res[0] += cx * (* myfunc)(x);
// 	break;
//       default:
// 	for(j=0;j<Nchange;++j){
// 	  x[ichange] = xchange[j];
// 	  res[j] += cx * (* myfunc)(x);
// 	}
// 	break;
//       }
//       err = 1. - resref/res[0] ;
//     }
//     break;

//     // log steps (accurate if relevant)
//   case 1:
//     // We make the change of variable y=log(x+x0) such that x0>-x over the full range.
//     // if xmin = 0 or xmax = 0 or xmax/xmin<0, we make a change y = x + x0 such that ymin>0
//     // so that we can define log steps for y.
//     x0 = myMax(0.1-xmin,0.1);
//     dexmin = log(xmin+x0);
//     dexmax = log(xmax+x0);
//     dx = log((xmax+x0)/(xmin+x0))/(double)Nsteps;

//     for(i=0;i<=Nsteps;++i){
//       resref = res[0];
//       cx = simpson_coef(i,Nsteps);
//       ex = dexmin + i*dx;// y = log(x+x0)
//       Jacobian = exp(ex);// Jacobian = x+x0 = exp(y)
//       x[ivar] = Jacobian - x0;

//       switch(ichange){
//       case -1:
// 	res[0] += cx * Jacobian * (* myfunc)(x);
// 	break;
//       default:
// 	for(j=0;j<Nchange;++j){
// 	  x[ichange] = xchange[j];
// 	  res[j] += cx * Jacobian * (* myfunc)(x);
// 	}
// 	break;
//       }
//       // Just above, we have considered a change of variable y=log(x+corr). The Jacobian
//       // for this change is dx/dy = x+corr, so that dx f(x) -> dy |dx/dy| f(x) = dy (x+corr) f(x).
//       err = 1. - resref/res[0] ;
//     }
//     break;

//   // Exponential or super-exponential steps: x -> y = exp((x/xmax)^nexp), where nexp an odd
//   // positive integer (to the correct maping of real numbers - if even, then y(-x) = y(x) would
//   // lead to errors).
//   // The inverse transformation is x = xmax*sign(log(y))*pow(abs(log(y)), 1/nexp).
//   // Note that we can recover negative signs from the above inversion.
//   // The Jacobian is dx/dy = (xmax/nexp) * pow(x/xmax,1-n) * exp(-pow(x/xmax,nexp))
//   //                       = (xmax/nexp) * pow(x/xmax,1-n)/ y.
//   //
//   case 2: // (accurate if relevant)
//     //x0 = myMax(1.e-5*xmax,xmin);// this induces instabilities (to be checked)
//     x0 = 0.0;
//     dexmin = exp((xmin+x0)/xmax); dexmax = exp(1.+x0/xmax);
//     dx = (dexmax-dexmin)/(double)Nsteps;
//     for(i=0;i<=Nsteps;++i){
//       resref = res[0];
//       cx = simpson_coef(i,Nsteps);
//       ex = dexmin + i*dx; // evaluate y = exp(x/xmax) (y_i+1 = y_i + dy)
//       Jacobian = xmax/ex;
//       x[ivar] = xmax*log(ex)-x0;

//       switch(ichange){
//       case -1:
// 	res[0] += cx * Jacobian * (* myfunc)(x);
// 	break;
//       default:
// 	for(j=0;j<Nchange;++j){
// 	  x[ichange] = xchange[j];
// 	  res[j] += cx * Jacobian * (* myfunc)(x);
// 	}
// 	break;
//       }
//       // Just above, we have considered a change of variable y=exp((x/xmax)^n), where n=2m+1 to be
//       // able to recover changes in signs. The Jacobian for this change is
//       // dx/dy = (xmax/(n y)) (x/xmax)^(1-n), so that dx f(x) -> dy |dx/dy| f(x).
//       err = 1. - resref/res[0] ;
//     }
//     break;

//   case 4:
//     nexp = (ilog-2)*2+1;// the power needs to be odd to be able to recover negative signs.
//     if(xmin<=0.0) Error_mess(2,funname,"The lower bound must be > 0 for ilog>=2");
//     _xmax = xmin*pow(1.e-15,-1./nexp);// to avoid numerical crash
//     if(_xmax<xmax){
//       Error_mess(0,funname,"WARNING: numerical integral in too extended a range -> split");
//       vector<double> res1=res,res2=res;
//       Simpson_Integral1(ivar,isize,ilog,Nsteps/2,
//     			xmin,_xmax,xx,ichange,xchange,myfunc,res1,err);
//       Simpson_Integral1(ivar,isize,ilog,Nsteps/2,
//     			_xmax,xmax,xx,ichange,xchange,myfunc,res2,err);
//       for(j=0;j<Nchange;++j) res[j] = res1[j]+res2[j];
//       return;
//     }
//     dexmin = exp(pow(xmin/xmax,nexp)); dexmax = exp(1.0);
//     dx = (dexmax-dexmin)/(double)Nsteps;
//     for(i=0;i<=Nsteps;++i){
//       resref = res[0];
//       cx = simpson_coef(i,Nsteps);
//       ex = dexmin + i*dx; // evaluate y = exp( (x/xmax)^nexp) (y_i+1 = y_i + dy)
//       x[ivar] = xmax*mySign(log(ex))*pow(myAbs(log(ex)),1./nexp);
//       Jacobian = (xmax/(nexp*ex)) * pow(x[ivar]/xmax,1-nexp);

//       switch(ichange){
//       case -1:
// 	res[0] += cx * Jacobian * (* myfunc)(x);
// 	break;
//       default:
// 	for(j=0;j<Nchange;++j){
// 	  x[ichange] = xchange[j];
// 	  res[j] += cx * Jacobian * (* myfunc)(x);
// 	}
// 	break;
//       }
//       // Just above, we have considered a change of variable y=exp((x/xmax)^n), where n=2m+1 to be
//       // able to recover changes in signs. The Jacobian for this change is
//       // dx/dy = (xmax/(n y)) (x/xmax)^(1-n), so that dx f(x) -> dy |dx/dy| f(x).
//       err = 1. - resref/res[0] ;
//     }
//     break;

//   default:
//     Error_mess(2,funname,"Integration method undefined");
//     break;
//   }

//   err = myAbs(err);

//   // Multiply once for all by dx.
//   for(j=0;j<Nchange;++j) res[j] = res[j]*dx;

//   return;
// }

// ///
// /// Function that discretizes a given range as a vector of length N according to a selected method:
// /// linearly, etc. Similar to the type of sampling in the function Simpson 1 (see below).
// ///
// void my_sample_vector(std::vector<double> par,std::vector<double> &xx)
// {
//   using namespace std;
//   string funcname = "my_sample_vector";
//   if(par.size()<4) Error_mess(2,funcname,"param vector size incorrect");

//   int ilog = (int)par[0];// method to discretize the range
//   int i,nexp,Nsteps = (int)par[1];// number of points
//   double x,ex,xmin,xmax,_xmax,dx,dexmin,dexmax,alpha,x0,Jacobian;

//   // Boundaries
//   xmin = par[2];
//   xmax = par[3];
//   if(xmin>xmax) Error_mess(2,funcname,"the range limits should be ordered");

//   xx.resize(0);

//   switch(ilog){
//     // inverted super-log steps (inspired by Thomas Lacroix)
//   case -4:
//     // This method is aimed at providing an alternative to the super-exponential method which
//     // may lead to numerical overflow.
//     // We make two changes of variables:
//     // (i) xbar = xmin+xmax-x
//     // (ii) y = log(log(alpha(xbar+x0))), such that x0 = max(1,-xmin+1)
//     // and alpha = 2/(xmin+x0) ensures that ymin = log(2).
//     // The Jacobian reads dxbar/dy = (xbar+x0) log(alpha(xbar+x0))
//     x0 = myMax(1.0,1.0-xmin);
//     alpha = 2./(xmin+x0);
//     dexmin = log(log(alpha*(xmin+x0)));
//     dexmax = log(log(alpha*(xmax+x0)));
//     dx = log(log(alpha*(xmax+x0))/log(alpha*(xmin+x0)))/Nsteps;
//     for(i=0;i<=Nsteps;++i){
//       ex = dexmin+i*dx;
//       ex = exp(ex);
//       Jacobian = ex*exp(ex)/alpha;
//       xx.push_back(xmin+xmax - (Jacobian/ex - x0));
//     }
//     break;

//     // inverse log steps
//   case -3:
//     // Here we use the change of variable y=log(alpha/(x+x0)), such that x0=xmin gives an order
//     // of magnitude for the considered numbers, and alpha = 2(xmax+x0) allows y>log(2). This is
//     // similar to the exp change of variable in spirit (see case 2).
//     // We then have x = alpha*exp(-y)-x0 and the Jacobian dx/dy = - (x+x0).
//     // We can then scan the integration range from ymin = y(xmax) to ymax = y(xmin)
//     // => we change the sign of the Jacobian accordingly.
//     x0 = myMax(1.e-2,xmin);
//     alpha = 2*(xmax+x0);
//     dexmin = log(alpha/(xmax+x0));
//     dexmax = log(alpha/(xmin+x0));
//     dx = log((xmax+x0)/(xmin+x0))/Nsteps;
//     // Above, dx is defined such that dx>0. We integrate from ymin=y(xmax) up to ymax=y(xmin).

//     for(i=0;i<=Nsteps;++i){
//       ex = dexmin+i*dx;
//       // increment on the new variable y = log(1/(x+x0)). Note that exp(-y) = x+x0.
//       Jacobian = alpha*exp(-ex);
//       xx.push_back(Jacobian - x0);
//     }

//     break;
//     // inverse log steps

//   case -2:
//     // Here we use the change of variable y=log(xmax-x+x0) which may apply to any range such that
//     // x0>x-xmax. This is similar to the exp change of variable in spirit (see case 2).
//     // We then have x = xmax+x0-exp(y) and the Jacobian dx/dy = - (xmax+x0-x).
//     // We can then scan the integration range from ymin = log(x0) to
//     // ymax = log(xmax-xmin+x0) => we change the sign of the Jacobian accordingly.
//     x0 = myMax(10.,xmin);
//     dexmin = log(x0);
//     dexmax = log(xmax-xmin+x0);
//     dx = log((xmax-xmin+x0)/(x0))/Nsteps;
//     // Above, dx is defined such that dx>0. We integrate from ymin=y(xmax) up to ymax=y(xmin).

//     for(i=0;i<=Nsteps;++i){
//       ex = dexmin+i*dx;// we start from ymin = y(xmax)
//       // increment on the new variable y = log(xmax-x+x0). Note that exp(y) = xmax-x+x0.
//       Jacobian = exp(ex);
//       xx.push_back(xmax + x0 - Jacobian);
//     }
//     break;
//     // super-log steps
//   case -1:
//     // Here we use the change of variable y=log(log(x+x0)) which may apply to any x>-x0+exp(delta).
//     // We then have x = exp(exp(y))-x0 and the Jacobian dx/dy = (x+x0) log(x+x0).
//     // We can then scan the integration range from ymin = log(1/(xmax+x0)) to
//     // ymax = log(1/(xmin+x0)).
//     x0 = myMax(exp(1.e-5)-xmin,exp(1.e-5));
//     dexmax = log(log(xmax+x0));
//     dexmin = log(log(xmin+x0));
//     dx = log(log(xmax+x0)/log(xmin+x0))/(double)Nsteps;

//     for(i=0;i<=Nsteps;++i){
//       ex = dexmin+i*dx;
//       Jacobian = exp(exp(ex));// x+x0 = exp(exp(y))
//       xx.push_back(Jacobian - x0);
//     }
//     break;

//     // linear steps
//   case 0:
//     dx = (xmax-xmin)/(double)Nsteps;
//     for(i=0;i<=Nsteps;++i){
//       xx.push_back(xmin + i*dx);
//     }
//     break;

//     // log steps
//   case 1:
//     // We make the change of variable y=log(x+x0) such that x0>-x over the full range.
//     // if xmin = 0 or xmax = 0 or xmax/xmin<0, we make a change y = x + x0 such that ymin>0
//     // so that we can define log steps for y.
//     x0 = myMax(0.1-xmin,0.1);
//     dexmin = log(xmin+x0);
//     dexmax = log(xmax+x0);
//     dx = log((xmax+x0)/(xmin+x0))/(double)Nsteps;

//     for(i=0;i<=Nsteps;++i){
//       ex = dexmin + i*dx;// y = log(x+x0)
//       Jacobian = exp(ex);// Jacobian = x+x0 = exp(y)
//       xx.push_back(Jacobian - x0);
//     }
//     break;

//   // Exponential or super-exponential steps: x -> y = exp((x/xmax)^nexp), where nexp an odd
//   // positive integer (to the correct maping of real numbers - if even, then y(-x) = y(x) would
//   // lead to errors).
//   // The inverse transformation is x = xmax*sign(log(y))*pow(abs(log(y)), 1/nexp).
//   // Note that we can recover negative signs from the above inversion.
//   // The Jacobian is dx/dy = (xmax/nexp) * pow(x/xmax,1-n) * exp(-pow(x/xmax,nexp))
//   //                       = (xmax/nexp) * pow(x/xmax,1-n)/ y.
//   //
//   case 2:
//     //x0 = myMax(1.e-5*xmax,xmin);// this induces instabilities (to be checked)
//     x0 = 0.0;
//     dexmin = exp((xmin+x0)/xmax); dexmax = exp(1.+x0/xmax);
//     dx = (dexmax-dexmin)/(double)Nsteps;
//     for(i=0;i<=Nsteps;++i){
//       ex = dexmin + i*dx; // evaluate y = exp(x/xmax) (y_i+1 = y_i + dy)
//       Jacobian = xmax/ex;
//       xx.push_back(xmax*log(ex)-x0);
//     }
//     break;
//   case 4:
//     nexp = (ilog-2)*2+1;// the power needs to be odd to apply to negative signs.
//     if(xmin<=0.0) Error_mess(2,funcname,"The lower bound must be > 0 for ilog>=2");
//     _xmax = xmin*pow(1.e-15,-1./nexp);// to avoid numerical crash
//     if(_xmax<xmax){
//       Error_mess(0,funcname,"WARNING: range too extended wrt the sampling method -> split");
//       vector<double> x1,x2,parr=par;
//       parr[1] = par[1]/2; parr[3] = _xmax; my_sample_vector(parr,x1);
//       parr[2] = _xmax; parr[3] = xmax; my_sample_vector(parr,x2);
//       x1.insert(x1.end(),x2.begin(),x2.end());
//       xx.insert(xx.end(),x1.begin(),x1.end());
//       return;
//     }
//     dexmin = exp(pow(xmin/xmax,nexp)); dexmax = exp(1.0);
//     dx = (dexmax-dexmin)/(double)Nsteps;
//     for(i=0;i<=Nsteps;++i){
//       ex = dexmin + i*dx; // evaluate y = exp( (x/xmax)^nexp) (y_i+1 = y_i + dy)
//       xx.push_back(xmax*mySign(log(ex))*pow(myAbs(log(ex)),1./nexp));
//     }
//     break;

//   default:
//     Error_mess(2,funcname,"Sampling method undefined");
//     break;
//   }
//   return;
// }

double trapezoid_integral(int ivar, int isize,
                          std::vector<double> xx, std::vector<double> interval,
                          double (*myfunc)(std::vector<double>), double &err)
{
    using namespace std;
    string funname = "trapezoid_integral";
    ostringstream oss;

    int i, ndim = xx.size(); // size of the vector in argument of myfunc
    if (isize != ndim)
    {
        cout << "Check vector length!!!\n";
        exit(1);
    }
    if (ivar > (ndim - 1))
    {
        cout << "Check vector length!!!\n";
        exit(1);
    }

    int Nsteps = interval.size(); // number of integration steps

    vector<double> x0 = xx, x1 = xx; // xx is not modified

    double dx, resref, res = 0;
    for (i = 0; i < Nsteps - 1; ++i)
    {
        resref = res;
        x0[ivar] = interval[i], x1[ivar] = interval[i + 1];
        dx = x1[ivar] - x0[ivar];
        res += 0.5 * dx * ((*myfunc)(x0) + (*myfunc)(x1));

        if ((i > (Nsteps / 4)) && res != 0.)
            err = 1. - resref / res;
    }
    err = myAbs(err);

    return res;
}

// copy of the Python "simps" function

double Simpson_Integral2(int ivar, int isize,
                         std::vector<double> xx, std::vector<double> interval,
                         double (*myfunc)(std::vector<double>), double &err)
{
    using namespace std;
    string funname = "Simpson_Integral2";
    ostringstream oss;

    int i, ndim = xx.size(); // size of the vector in argument of myfunc
    if (isize != ndim)
    {
        cout << "Check vector length!!!\n";
        exit(1);
    }
    if (ivar > (ndim - 1))
    {
        cout << "Check vector length!!!\n";
        exit(1);
    }

    int Nsteps = interval.size(); // number of integration steps

    vector<double> h(Nsteps - 1), xx0 = xx, xx1 = xx, xx2 = xx;

    for (i = 0; i < Nsteps - 1; ++i) // filling vector of differences
    {
        h[i] = interval[i + 1] - interval[i];
    }

    double h0, h1, hsum, hprod, h0divh1, yi0, yi1, yi2, res = 0, resref;
    for (i = 0; i < Nsteps - 2; i += 2)
    {
        resref = res;
        xx0[ivar] = interval[i], xx1[ivar] = interval[i + 1], xx2[ivar] = interval[i + 2];
        h0 = h[i], h1 = h[i + 1];
        hsum = h0 + h1, hprod = h0 * h1, h0divh1 = h0 / h1;
        yi0 = (*myfunc)(xx0), yi1 = (*myfunc)(xx1), yi2 = (*myfunc)(xx2);
        res += hsum / 6. * (yi0 * (2. - 1. / h0divh1) + yi1 * hsum * hsum / hprod + yi2 * (2. - h0divh1));

        if ((i > (Nsteps / 4)) && res != 0.)
            err = 1. - resref / res;
    }
    err = myAbs(err);
    return res;
}
// double Simpson_Integral2(int ivar, int isize,
// 			 std::vector<double> xx, std::vector<double> interval,
// 			 double (*myfunc)(std::vector<double>), double &err)
// {
//   using namespace std;
//   string funname = "Simpson_Integral2";
//   ostringstream oss;

//   int i, ndim = xx.size(); // size of the vector in argument of myfunc
//   if(isize!=ndim) { cout << "Check vector length!!!\n"; exit(1); }
//   if(ivar > (ndim-1)) {cout << "Check vector length!!!\n"; exit(1); }

//   int Nsteps = interval.size(); // number of integration steps

//   vector<double> h(Nsteps-1), xx0 = xx, xx1 = xx, xx2 = xx;

//   for(i=0;i<Nsteps-1;++i) // filling vector of differences
//     {h[i] = interval[i+1]-interval[i];}

//   double h0,h1,hsum,hprod,h0onh1,yi,yi1,yi2,res = 0,resref;
//   for(i=0;i<Nsteps-2;i+=2)
//     {
//       resref = res;
//       xx0[ivar] = interval[i], xx1[ivar] = interval[i+1], xx2[ivar] = interval[i+2];
//       h0 = h[i], h1 = h[i+1];
//       hsum = h0+h1, hprod = h0*h1, h0onh1 = h0/h1;
//       yi = (*myfunc)(xx0), yi1 = (*myfunc)(xx1), yi2 = (*myfunc)(xx2);
//       res += hsum/6.*( (2-1./h0onh1)*yi
//   		       + pow(hsum,2)/hprod*yi1
//   		       + (2-h0onh1)*yi2 );

//       if((i > (Nsteps/4)) && res!=0.) err = 1. - resref/res;
//     }
//   err = myAbs(err);
//   return res;
// }

///
/// Numerical derivative of order n, for a given function.
/// Richardson extrapolation?
///
/*
double myDerivative(int n, int ivar, std::vector<double> xx,
                    double (*myfunc)(std::vector<double>), double eps)
{
    using namespace std;
    string funname = "myDerivative";
    double res = 0;
    int isize = xx.size();
    if (isize <= ivar)
        cout << "Check myderiv parameters" << endl;

    vector<double> _xx = xx;
    double x = xx[ivar];
    double delta = eps * x;

    int i, nn = n;
    double fp[8], fm[8];

    const double ap[8] = {6.4064e5, -2.24224e5, 8.1536e4, -2.548e4, 6.272e3, -1.12e3, 1.28e2, -7.};

    if (n >= 1)
    {
        nn--;
        res = 0;
        for (i = 0; i < 8; ++i)
        {
            _xx[ivar] = x + (i + 1) * delta;
            fp[i] = myDerivative(nn, ivar, _xx, myfunc, eps);
            _xx[ivar] = x - (i + 1) * delta;
            fm[i] = myDerivative(nn, ivar, _xx, myfunc, eps);
            res += ap[i] * fp[i] - ap[i] * fm[i];
        }
        res = res / (720720. * delta);
    }
    else if (n == 0)
        res = myfunc(xx);
    return res;
}*/


double myDerivative_Static(int n, int ivar, std::vector<double> xx, void *pt2Object, double (*myfunc)(void *pt2Object, std::vector<double>), double eps)
{
    using namespace std;
    string funname = "myDerivative";
    double res = 0;
    int isize = xx.size();
    if (isize <= ivar)
        cout << "Check myderiv parameters" << endl;

    vector<double> _xx = xx;
    double x = xx[ivar];
    double delta = eps * x;

    int i, nn = n;
    double fp[8], fm[8];

    const double ap[8] = {6.4064e5, -2.24224e5, 8.1536e4, -2.548e4, 6.272e3, -1.12e3, 1.28e2, -7.};

    if (n >= 1)
    {
        nn--;
        res = 0;
        for (i = 0; i < 8; ++i)
        {
            _xx[ivar] = x + (i + 1) * delta;
            fp[i] = myDerivative_Static(nn, ivar, _xx, pt2Object, myfunc, eps);
            _xx[ivar] = x - (i + 1) * delta;
            fm[i] = myDerivative_Static(nn, ivar, _xx, pt2Object, myfunc, eps);
            res += ap[i] * fp[i] - ap[i] * fm[i];
        }
        res = res / (720720. * delta);
    }
    else if (n == 0)
        res = myfunc(pt2Object, xx);
    return res;
}

///
/// Runge-Kutta 4 method. We seek for y(x+dx) given y(x) = y_in and
/// dy/dx = f(x,y).
///
double myRungeKutta4(int ivar, double y_in, std::vector<double> xx, double &dx,
                     double (*dytodx)(std::vector<double>, double))
{
    using namespace std;
    string funname = "myRungeKutta4";
    double y = y_in;
    vector<double> x1, x2;
    x1 = x2 = xx;
    double k1 = 0, k2 = 1, k3, k4;
    const double eps = 1.e-5;
    double err = myAbs(k1 - k2);
    dx *= 2;
    int icount = 0;
    int imax = 50;
    while (myAbs(k1 - k2) > eps && icount < imax)
    {
        dx *= 0.5;
        k1 = dx * dytodx(xx, y);
        x1[ivar] = xx[ivar] + 0.5 * dx;
        x2[ivar] = xx[ivar] + dx;
        k2 = dx * dytodx(x1, y + 0.5 * k1);
        k3 = dx * dytodx(x1, y + 0.5 * k2);
        //k4 = dx * dytodx(x2,y+0.5*k3);
        k4 = dx * dytodx(x2, y + k3); // modif 08/06/2016
        icount++;
    }

    if (icount == imax)
        printf("dx too large\n");

    double ynext = y + (k1 + 2 * (k2 + k3) + k4) / 6.0;
    return ynext;
}

double myRungeKutta4_Static(int ivar, double y_in, std::vector<double> xx, double &dx,
                               void *pt2Object, double (*dytodx)(void *pt2Object, std::vector<double>, double))
{
    using namespace std;
    string funname = "myRungeKutta4";
    double y = y_in;
    vector<double> x1, x2;
    x1 = x2 = xx;
    double k1 = 0, k2 = 1, k3, k4;
    const double eps = 1.e-5;
    double err = myAbs(k1 - k2);
    dx *= 2;
    int icount = 0;
    int imax = 50;
    while (myAbs(k1 - k2) > eps && icount < imax)
    {
        dx *= 0.5;
        k1 = dx * dytodx(pt2Object, xx, y);
        x1[ivar] = xx[ivar] + 0.5 * dx;
        x2[ivar] = xx[ivar] + dx;
        k2 = dx * dytodx(pt2Object, x1, y + 0.5 * k1);
        k3 = dx * dytodx(pt2Object, x1, y + 0.5 * k2);
        //k4 = dx * dytodx(x2,y+0.5*k3);
        k4 = dx * dytodx(pt2Object, x2, y + k3); // modif 08/06/2016
        icount++;
    }

    if (icount == imax)
        printf("dx too large\n");

    double ynext = y + (k1 + 2 * (k2 + k3) + k4) / 6.0;
    return ynext;
}

// example:
//
// int N = 500, ivar = 0;
// double err, dx = 1e-3;
// vector<double> xx(1);
// xx[0] = 0;
// double y_in = 1; // y_in = y(xx[0])
// for(int i=0;i<=N;++i)
//   {
//    xx[0] += dx;
//    y_in = myRungeKutta4(ivar,y_in,xx,dx,&deriv);
//    cout << xx[0] << "   " << y_in << endl;
//   }
// return;

///
/// Based on the Fehlberg method (cf theory in Numerical Recipes).
/// Evaluates the error and provides adaptative steps.
///
double myRungeKutta6(int ivar, std::vector<double> x, double y0,
                     void (*derivs)(std::vector<double>, double, double &),
                     bool varstep, double &dx, double dxmin,
                     double prec, bool &okprec)
{
    using namespace std;

    // Cash-Karp parameters
    const double aa[6] = {0, 0.2, 0.3, 3. / 5., 1, 7. / 8.};
    const double bb[6][5] = {{0, 0, 0, 0, 0},
                             {0.2, 0, 0, 0, 0},
                             {0.3 / 4., 0.9 / 4., 0, 0, 0},
                             {0.3, -0.9, 1.2, 0, 0},
                             {-11. / 54., 2.5, -70. / 27., 35. / 27., 0},
                             {16231. / 55296, 175. / 512., 575. / 13284., 44275. / 110592.}};
    const double cc[6] = {37. / 378., 0, 250. / 621., 125. / 594., 0, 512. / 1771.};
    const double _cc[6] = {2825. / 27648., 0, 18575. / 48384., 13525. / 55296., 277. / 14336., 0.25};

    int i, j, k, N;
    double dxx, err;
    err = dxx = 0;
    double kk[6] = {0};
    double xx[6] = {0};

    vector<double> _x = x;
    double ynew, dytodx;
    ynew = y0;

    for (i = 0; i < 6; ++i)
        xx[i] = _x[ivar] + aa[i] * dx;

    _x[ivar] = xx[0];
    derivs(_x, y0, dytodx); // get dy/dx=f(x,y)
    kk[0] = dx * dytodx;    // k0 = dx*f(x,y)

    for (i = 1; i < 6; ++i)
    {
        kk[i] = 0;
        for (j = 0; j < i; ++j)
            kk[i] += bb[i][j] * kk[j]; // ex: k3 = b30*k0+b31*k1+b32*k2
        _x[ivar] = xx[i];
        derivs(_x, y0 + kk[i], dytodx);
        kk[i] = dx * dytodx;
    }

    ynew = y0;
    for (i = 0; i < 6; ++i)
    {
        ynew += cc[i] * kk[i];
        err += (cc[i] - _cc[i]) * kk[i]; // error estimate
    }

    // checking the step size.
    dxx = dx * pow(myAbs(err / prec), 0.2);
    if (varstep)
        dx = myMax(dxx, dxmin); // adapt the step if required
    else if (dxx < dx)
        okprec = false;
    return ynew;
}

double myRungeKutta6_Static(int ivar, std::vector<double> x, double y0, void *pt2Object,
                               void (*derivs)(void *pt2Object, std::vector<double>, double, double &),
                               bool varstep, double &dx, double dxmin,
                               double prec, bool &okprec)
{
    using namespace std;

    // Cash-Karp parameters
    const double aa[6] = {0, 0.2, 0.3, 3. / 5., 1, 7. / 8.};
    const double bb[6][5] = {{0, 0, 0, 0, 0},
                             {0.2, 0, 0, 0, 0},
                             {0.3 / 4., 0.9 / 4., 0, 0, 0},
                             {0.3, -0.9, 1.2, 0, 0},
                             {-11. / 54., 2.5, -70. / 27., 35. / 27., 0},
                             {16231. / 55296, 175. / 512., 575. / 13284., 44275. / 110592.}};
    const double cc[6] = {37. / 378., 0, 250. / 621., 125. / 594., 0, 512. / 1771.};
    const double _cc[6] = {2825. / 27648., 0, 18575. / 48384., 13525. / 55296., 277. / 14336., 0.25};

    int i, j, k, N;
    double dxx, err;
    err = dxx = 0;
    double kk[6] = {0};
    double xx[6] = {0};

    vector<double> _x = x;
    double ynew, dytodx;
    ynew = y0;

    for (i = 0; i < 6; ++i)
        xx[i] = _x[ivar] + aa[i] * dx;

    _x[ivar] = xx[0];
    derivs(pt2Object, _x, y0, dytodx); // get dy/dx=f(x,y)
    kk[0] = dx * dytodx;               // k0 = dx*f(x,y)

    for (i = 1; i < 6; ++i)
    {
        kk[i] = 0;
        for (j = 0; j < i; ++j)
            kk[i] += bb[i][j] * kk[j]; // ex: k3 = b30*k0+b31*k1+b32*k2
        _x[ivar] = xx[i];
        derivs(pt2Object, _x, y0 + kk[i], dytodx);
        kk[i] = dx * dytodx;
    }

    ynew = y0;
    for (i = 0; i < 6; ++i)
    {
        ynew += cc[i] * kk[i];
        err += (cc[i] - _cc[i]) * kk[i]; // error estimate
    }

    // checking the step size.
    dxx = dx * pow(myAbs(err / prec), 0.2);
    if (varstep)
        dx = myMax(dxx, dxmin); // adapt the step if required
    else if (dxx < dx)
        okprec = false;
    return ynew;
}

// example:
//
// int ivar = 0, N = 500;
// vector<double> xx(1);
// xx[0] = 0;
// double y0 = 1, err;
// bool varstep = false;
// double dx = 1e-2, dxmin = dx, prec = 1e-3;
// bool okprec;
// for(int i=0;i<=N;++i)
//   {
//     xx[0] += dx;
//     y0 = myRungeKutta6(ivar,xx,y0,derivs,
// 			 varstep,dx,dxmin,
// 			 prec,okprec);
//     cout << xx[0] << "   " << y0 << endl;
//   }
// return;

//
// Implicit Euler solver for specific equations dY(x)/dx = A(x)(Z(x)^2-Y(x)^2)
double ImplicitEulerSolverSquares_Static(int i, std::vector<double> x, double y0, void *pt2Object,
                                            double (*func_A)(void *pt2Object, std::vector<double>),
                                            double (*func_Z)(void *pt2Object, std::vector<double>),
                                            double const &dx)
{
    std::vector<double> _x = x;

    double A = func_A(pt2Object, _x);
    double Z = func_Z(pt2Object, _x);

    double y1 = (sqrt(1 + 4 * A * dx * (A * dx * pow(Z, 2) + y0)) - 1) / (2 * A * dx);

    return y1;
}

// Implicit Euler solver for specific equations dY(x)/dx = (A(x)*Z(x)^2-B(x)*Y(x)^2)
// Overriden version of a previous solver where B = A
double ImplicitEulerSolverSquares_Static(int i, std::vector<double> x, double y0, void *pt2Object,
                                            double (*func_A)(void *pt2Object, std::vector<double>),
                                            double (*func_B)(void *pt2Object, std::vector<double>),
                                            double (*func_Z)(void *pt2Object, std::vector<double>),
                                            double const &dx)
{
    std::vector<double> _x = x;

    double A = func_A(pt2Object, _x);
    double B = func_B(pt2Object, _x);
    double Z = func_Z(pt2Object, _x);

    double y1 = (sqrt(1 + 4 * B * dx * (A * dx * pow(Z, 2) + y0)) - 1) / (2 * B * dx);

    return y1;
}

//
// Implicit Euler solver for specific equations dy(x)/dx = A(x)(Z(x)-y(x))
double ImplicitEulerSolverLinear_Static(int i, std::vector<double> x, double y0, void *pt2Object,
                                           double (*func_A)(void *pt2Object, std::vector<double>),
                                           double (*func_Z)(void *pt2Object, std::vector<double>),
                                           double const &dx)
{
    std::vector<double> _x = x;

    double A = func_A(pt2Object, _x);
    double Z = func_Z(pt2Object, _x);

    double y1 = (y0 + A * Z * dx) / (1 + A * dx);

    return y1;
}

// Implicit Backward differentiation solver of order 5 for specific equations dy(x)/dx = A(x)(Z(x)-y(x))
double ImplicitBackwardDiffSolver5Linear_Static(int i, std::vector<double> x, double y0, double ym1, double ym2, double ym3, double ym4, void *pt2Object,
                                                   double (*func_A)(void *pt2Object, std::vector<double>),
                                                   double (*func_Z)(void *pt2Object, std::vector<double>),
                                                   double const &dx)
{
    std::vector<double> _x = x;

    double A = func_A(pt2Object, _x);
    double Z = func_Z(pt2Object, _x);

    double y1 = (300. / 137. * y0 - 300. / 137. * ym1 + 200. / 137. * ym2 - 75. / 137. * ym3 + 12. / 137. * ym4 + 60. / 137. * A * Z * dx) / (1 + 60. / 137. * A * dx);

    return y1;
}

// Implicit Adam-Moulton solver of order 5 for specific equations dy(x)/dx = A(x)(Z(x)-y(x))
double ImplicitAdamsMoultonSolver5Linear_Static(int i, std::vector<double> x, std::vector<double> x0, std::vector<double> xm1, std::vector<double> xm2, std::vector<double> xm3,
                                                   double y0, double ym1, double ym2, double ym3, void *pt2Object,
                                                   double (*func_A)(void *pt2Object, std::vector<double>),
                                                   double (*func_Z)(void *pt2Object, std::vector<double>),
                                                   double const &dx)
{
    std::vector<double> _x = x;
    std::vector<double> _x0 = x0;
    std::vector<double> _xm1 = xm1;
    std::vector<double> _xm2 = xm2;
    std::vector<double> _xm3 = xm3;

    double A = func_A(pt2Object, _x);
    double Z = func_Z(pt2Object, _x);

    double func0 = func_A(pt2Object, _x0) * (func_Z(pt2Object, _x0) - y0);
    double funcm1 = func_A(pt2Object, _xm1) * (func_Z(pt2Object, _xm1) - ym1);
    double funcm2 = func_A(pt2Object, _xm2) * (func_Z(pt2Object, _xm2) - ym2);
    double funcm3 = func_A(pt2Object, _xm3) * (func_Z(pt2Object, _xm3) - ym3);

    double y1 = (y0 + dx * (251. / 720. * A * Z + 646. / 720. * func0 - 264. / 720. * funcm1 + 106. / 720. * funcm2 - 19. / 720. * funcm3)) / (1 + 251. / 720. * dx * A);

    return y1;
}

//

double Dichotomie_Static(void *pt2Object, double (*myfunc)(void *pt2Object, double),
                            double const &xMin, double const &xMax, double const &prec, double dflt)
{

    double _y1, _y2, _zz;
    int compt, comptMax = 100000;

    compt = 0;

    // Initialisation
    _y1 = xMin;
    _y2 = xMax;

    double _myfunc_y1, _myfunc_y2, _myfunc_zz;
    _myfunc_y1 = myfunc(pt2Object, _y1);
    _myfunc_y2 = myfunc(pt2Object, _y2);

    if (_myfunc_y1 == 0)
        return xMin;
    if (_myfunc_y2 == 0)
        return xMax;

    _zz = (_y1 + _y2) / 2.0;

    double _zzm1;

    do
    {
        _zzm1 = _zz;

        // Update the value of _my_func_zz
        _myfunc_zz = myfunc(pt2Object, _zz);

        double prodzy1 = _myfunc_zz * _myfunc_y1;
        double prodzy2 = _myfunc_zz * _myfunc_y2;

        //std::cout <<  _zz[ivar] << " " << myfunc(pt2Object, _zz) << " " << myfunc(pt2Object, _y1) << " " << myfunc(pt2Object, _y2) <<  std::endl;

         if (prodzy1 <= 0 && prodzy2 > 0)
        {
            _y2 = _zz;
            _myfunc_y2 = _myfunc_zz;
        }
        else if (prodzy1 >= 0 && prodzy2 < 0)
        {
            _y1 = _zz;
            _myfunc_y1 = _myfunc_zz;
        }
        else
        {
            //std::cout << "WARNING : Dichotomie not possible" << std::endl;
            return dflt; // Exiting function
        }

        _zz = (_y1 + _y2) / 2.0;

        compt++;

    } while (fabs((_zzm1 - _zz) / (_zzm1 + _zz)) > prec && compt < comptMax);

    if (compt == comptMax)
        std::cout << "WARNING : Dichotomie max precision may not be obtained : " << fabs(myfunc(pt2Object, _zz)) << std::endl;

    return _zz;
}


double Dichotomie_Static(int ivar, int isize, std::vector<double> xx, void *pt2Object,
                         double (*myfunc)(void *pt2Object, std::vector<double> xx),
                         double const &xMin, double const &xMax, double const &prec, double dflt)
{

    if (isize != xx.size())
    {
        std::cout << "FATAL ERROR : Wrong vector length in " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    std::vector<double> _y1 = xx, _y2 = xx, _zz = xx;
    int compt, comptMax = 100000;

    compt = 0;

    // Initialisation
    _y1[ivar] = xMin;
    _y2[ivar] = xMax;

    double _myfunc_y1, _myfunc_y2, _myfunc_zz;
    _myfunc_y1 = myfunc(pt2Object, _y1);
    _myfunc_y2 = myfunc(pt2Object, _y2);

    if (_myfunc_y1 == 0)
        return xMin;
    if (_myfunc_y2 == 0)
        return xMax;

    _zz[ivar] = (_y1[ivar] + _y2[ivar]) / 2.0;

    double _zzm1;

    do
    {
        _zzm1 = _zz[ivar];

        // Update the value of _my_func_zz
        _myfunc_zz = myfunc(pt2Object, _zz);

        double prodzy1 = _myfunc_zz * _myfunc_y1;
        double prodzy2 = _myfunc_zz * _myfunc_y2;

        //std::cout <<  _zz[ivar] << " " << myfunc(pt2Object, _zz) << " " << myfunc(pt2Object, _y1) << " " << myfunc(pt2Object, _y2) <<  std::endl;

         if (prodzy1 <= 0 && prodzy2 > 0)
        {
            _y2[ivar] = _zz[ivar];
            _myfunc_y2 = _myfunc_zz;
        }
        else if (prodzy1 >= 0 && prodzy2 < 0)
        {
            _y1[ivar] = _zz[ivar];
            _myfunc_y1 = _myfunc_zz;
        }
        else
        {
            //std::cout << "WARNING : Dichotomie not possible" << std::endl;
            return dflt; // Exiting function
        }

        _zz[ivar] = (_y1[ivar] + _y2[ivar]) / 2.0;

        compt++;

    } while (fabs((_zzm1 - _zz[ivar]) / (_zzm1 + _zz[ivar])) > prec && compt < comptMax);

    if (compt == comptMax)
        std::cout << "WARNING : Dichotomie max precision may not be obtained : " << fabs(myfunc(pt2Object, _zz)) << std::endl;

    return _zz[ivar];
}

double DichotomieLn_Static(int ivar, int isize, std::vector<double> xx, void *pt2Object,
                           double (*myfunc)(void *pt2Object, std::vector<double> xx),
                           double const &xMin, double const &xMax, double const &prec, double dflt)
{

    if (isize != xx.size())
    {
        std::cout << "FATAL ERROR : Wrong vector length in " << __PRETTY_FUNCTION__ << std::endl;
        exit(0);
    }

    std::vector<double> _y1 = xx, _y2 = xx, _zz = xx;
    int compt, comptMax = 100000;

    compt = 0;

    // Initialisation
    _y1[ivar] = xMin;
    _y2[ivar] = xMax;

    double _myfunc_y1, _myfunc_y2, _myfunc_zz;
    _myfunc_y1 = myfunc(pt2Object, _y1);
    _myfunc_y2 = myfunc(pt2Object, _y2);

    if (_myfunc_y1 == 0)
        return xMin;
    if (_myfunc_y2 == 0)
        return xMax;

    _zz[ivar] = pow(10, (log10(_y1[ivar]) + log10(_y2[ivar])) / 2.0);

    double _zzm1;

    do
    {
        _zzm1 = _zz[ivar];

        // Update the value of _my_func_zz
        _myfunc_zz = myfunc(pt2Object, _zz);

        double prodzy1 = _myfunc_zz * _myfunc_y1;
        double prodzy2 = _myfunc_zz * _myfunc_y2;

        //std::cout <<  _zz[ivar] << " " << myfunc(pt2Object, _zz) << " " << myfunc(pt2Object, _y1) << " " << myfunc(pt2Object, _y2) <<  std::endl;

        if (prodzy1 <= 0 && prodzy2 > 0)
        {
            _y2[ivar] = _zz[ivar];
            _myfunc_y2 = _myfunc_zz;
        }
        else if (prodzy1 >= 0 && prodzy2 < 0)
        {
            _y1[ivar] = _zz[ivar];
            _myfunc_y1 = _myfunc_zz;
        }
        else
        {
            //std::cout << "WARNING : Dichotomie not possible" << std::endl;
            return dflt; // Exiting function
        }

        _zz[ivar] = pow(10, (log10(_y1[ivar]) + log10(_y2[ivar])) / 2.0);

        compt++;

    } while (fabs((_zzm1 - _zz[ivar]) / (_zzm1 + _zz[ivar])) > prec && compt < comptMax);

    if (compt == comptMax)
        std::cout << "WARNING : Dichotomie max precision may not be obtained : " << fabs(myfunc(pt2Object, _zz)) << std::endl;

    return _zz[ivar];
}



//
//
//
//
//
// cardinal sinus: sinc(x) = sin(x)/x
//

double sinc(double x)
{
    using namespace std;
    double sc = 0;
    if (fabs(x) < 1e-3)
    {
        sc = 1 - x * x / 6. + pow(x, 4.) / 120. - pow(x, 6.) / 5040.;
    }
    else
    {
        sc = sin(x) / x;
    }
    return sc;
}

//
// exp(x^2)*erfc(x)
//

double erfcx(double x)
{
    using namespace std;
    double res = 0;
    if (fabs(x) > 5)
    {
        res = 1 / (x * sqrt(pi)) * (1 - 0.5 * pow(x, -2) + 0.75 * pow(x, -4) - 0.625 * pow(x, -6));
    }
    else
    {
        res = exp(x * x) * erfc(x);
    }
    return res;
}

double my_erfc(double x)
{
    if (x < 20)
        return erfc(x);
    else
        return exp(-x*x)*(1./sqrt(PI)/x - 1./2./sqrt(PI)/(x*x*x));
}

double my_erf_diff(double const &b, double const &a)
// Return erf(b) - erf(a)
{
    int eps = 1;

    double _b = b;
    double _a = a;

    if (b == a)
        return 0;

    if (b < a)
    {
        _b = a;
        _a = b;
        eps = -1;
    }

    //std::cout << _b << " " << _a << " " << eps << std::endl;

    if (_a >= 0)
        return eps * (erfc(_a) - erfc(_b));
    if (_b <= 0 && _a < 0)
        return eps * (erfc(-_b) - erfc(-_a));
    else
        return eps * (2 - erfc(-_a) - erfc(_b));

    /*
    if (_a >= 0)
        return eps*(gsl_sf_erfc(_a) - gsl_sf_erfc(_b));
    if (_b <= 0 && _a < 0)
        return eps*(gsl_sf_erfc(-_b) - gsl_sf_erfc(-_a));
    else
        return eps*(2 - gsl_sf_erfc(-_a) - gsl_sf_erfc(_b));*/
}

///