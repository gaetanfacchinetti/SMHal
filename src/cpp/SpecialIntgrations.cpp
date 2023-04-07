#include "../headers/SpecialIntegrations.h"

dcomp SpecialIntegrations::IntJt(int n, double a, double b, double tmin, double tmax)
{
    dcomp res = 0;

    double ymax;
    double ymin;

    if (b != 0)
    {
        ymax = (tmax - a) / b;
        ymin = (tmin - a) / b;
    }
    else
    {
        ymax = (tmax - a);
        ymin = (tmin - a);

        if (ymin * ymax < 0 || ymin == 0 || ymax == 0)
        {
            std::cout << "FATAL ERROR : IntJt -> Integral is not well defined" << std::endl;
            exit(0);
        }
    }

    // Initialisaiton of jk at k = 0
    dcomp j0 = 0.5 * (log(1 + ymax * ymax) - log(1 + ymin * ymin)) - one_i * (atan(ymax) - atan(ymin));
    dcomp jk, sumxi = 0;

    if (b != 0)
    {
        for (int k = 0; k < n + 1; k++)
        {
            if (k > 0)
                sumxi += pow(one_i, k) * (pow(ymax, k) - pow(ymin, k)) / (1.0 * k);
            jk = pow(-one_i, k) * (sumxi + j0);
            res += Binom(k, n) * pow(b, k) * pow(a, n - k) * jk;
        }
    }
    else
    {

        res = pow(a, n) * log(ymax / ymin);
        for (int k = 1; k < n + 1; k++)
            res += Binom(k, n) * pow(a, n - k) * (pow(ymax, k) - pow(ymin, k)) / (1.0 * k);
    }

    return res;
}

dcomp SpecialIntegrations::IntJu(int n, double r, double b, double tmin, double tmax)
{
    return pow(-1, 1) * IntJt(n, -r, b, tmin, tmax);
}

dcomp SpecialIntegrations::IntJuComp(int n, int m, double r, double b, double tmin, double tmax)
{
    return pow(-1, m) * IntJtComp(n, m, -r, -b, tmin, tmax);
}

double SpecialIntegrations::IntK(int n, double tmin, double tmax)
{
    return (1. / (n + 1)) * (pow(tmax, n + 1) - pow(tmin, n + 1));
}

int SpecialIntegrations::Binom(int k, int n)
{
    double res = 1;

    if (k > n)
    {
        std::cout << "ERROR   : SpecialIntegration::Binom() -> k>n" << std::endl;
        return -1;
    }

    for (int i = n; i > 0; i--)
    {

        if (i < k + 1 && i < (n - k + 1))
            res /= (1.0 * i);
        else if (i > k && i > (n - k))
            res *= i;
    }

    return (int)res;
}

int SpecialIntegrations::BinomGen(int k, int n)
{
    int num = 1, denom = 1;

    int compt = 0;
    int i = n;
    int j = 1;

    //std::cout << compt << " " << num << " " << denom << " " << i << " " << j << " " << std::endl;

    while (compt < k)
    {
        num *= i;
        i = i - 1;

        denom *= j;
        j = j + 1;

        compt++;

        //std::cout << compt << " " << num << " " << denom << " " << i << " " << j << " " << std::endl;
    }

    return (int)(num / denom);
}

int SpecialIntegrations::GG(int n, int m, int p, int q)
{
    int res = 0;

    for (int k = 0; k < n + 1; k++)
    {
        res += BinomGen(k, n) * BinomGen(p - 1, k - m) * BinomGen(q, k - m + 1) * pow(-1, k);
    }

    return (int)res;
}

dcomp SpecialIntegrations::FF(int n, int m, int p, dcomp x)
{
    dcomp res = 0;

    for (int k = 0; k < n + 1; k++)
    {
        res += (dcomp)(BinomGen(k, n) * BinomGen(p - 1, k - m)) * pow(x, k);
    }

    return res;
}

dcomp SpecialIntegrations::IntJtCompTest(int n, int m, double a, double b, double tmin, double tmax)
{
    dcomp res = 0;

    double ymin, ymax;
    dcomp zmax, zmin;

    dcomp z0 = a - one_i * b;

    if (b != 0)
    {

        zmax = (tmax - a + one_i * b) / (1. * a);
        zmin = (tmin - a + one_i * b) / (1. * a);
    }
    else
    {
        ymax = tmax / a - 1;
        ymin = tmin / a - 1;

        if (ymin * ymax < 0 || ymin == 0 || ymax == 0)
        {
            std::cout << "FATAL ERROR : IntJt -> Integral is not well defined" << std::endl;
            exit(0);
        }
    }

    dcomp u = zmin / zmax - (dcomp)1.;

    //std::cout << "oui " << abs(tmax / a) << " " << abs(u) << std::endl;

    if (abs(u) < 0.2 && abs(tmax / a) < 0.1)
    {
        dcomp res1 = 0;

        for (int p = 1; p < 20; p++)
        {
            res1 = (dcomp)GG(n, m, p, 0);

            for (int q = 1; q < 20; q++)
            {
                res1 += pow(-((dcomp)tmax) / z0, q) * ((dcomp)GG(n, m, p, q));
                //std::cout << "ici " << p << " " << q << " " << res1 << " " << pow(-tmax/z0, q) << " " << tmax << " " << z0 << std::endl;
            }

            res1 *= pow(u, p) / ((dcomp)p);
            res += res1;
        }

        res *= pow(-1, m) * pow(z0, n - m + 1);
    }

    if (abs(u) < 0.2 && abs(tmax / a) >= 0.1)
    {
        for (int p = 1; p < 20; p++)
        {
            res += pow(u, p) / ((dcomp)p) * (FF(n, m, p, a * zmax / z0)) * pow(a * zmax, -m + 1) * pow(z0, n);
            //std::cout << "ici " << p  << " " << res << " " << FF(n, m, p, a*zmax/z0) << " " << tmax << " " << z0 << std::endl;
        }

        res = -res;
    }

    return res;
}

dcomp SpecialIntegrations::IntJtComp(int n, int m, double a, double b, double tmin, double tmax)
// Integral of (t^n)/((t-a+ib)^m)
{

    /* Remark 
        In the package std::complex, log(std::complex<T> z) is defined with a
        branch cut on the negative axis. Therefore as soon as b is different
        from 0 there are no problem and any integration path will never go
        through the cut. However, when b = 0 one needs to compute specifically 
        the integral in order to be sure to have no problems. When m = 1 one 
        can use IntJt and IntJtComp indifferently however this method is more
        general and allows cases where m != 0 that we may need
    */

    dcomp res = 0;

    //double ymin, ymax;
    dcomp zmax, zmin;

    dcomp z0 = a - one_i * b;

    zmax = (tmax - a + one_i * b) / (1. * a);
    zmin = (tmin - a + one_i * b) / (1. * a);

    dcomp u = zmin / zmax - (dcomp)1.;

    if (abs(u) >= 0.2)
    {
        for (int k = 0; k < n + 1; k++)
        {
            if (k - m == -1)
            {
                res += ((dcomp)Binom(k, n)) * pow(z0, n - k) * pow(a, k - m + 1) * (log(zmax / zmin));
                //std::cout << "rrrr : " << zmax << " " << zmin << " " <<  ((dcomp)Binom(k,n))*pow(z0, n-k)*pow(a,k-m+1)*(log(zmax)-log(zmin)) << " " <<  ((dcomp)Binom(k,n))*pow(z0, n-k)*pow(a,k-m+1)*(log(zmax/zmin)) << std::endl;
            }
            else
            {
                res += ((dcomp)Binom(k, n)) * pow(z0, n - k) * pow(a, k - m + 1) * (1. / (k - m + 1)) * (pow(zmax, k - m + 1) - pow(zmin, k - m + 1));
                //std::cout << "ssss : " << (pow(zmax,k-m+1)-pow(zmin,k-m+1)) << " " << ((dcomp)Binom(k,n))*pow(z0, n-k)*pow(a,k-m+1)*(1./(k-m+1))*(pow(zmax,k-m+1)-pow(zmin,k-m+1)) << std::endl;
            }
        }
    }

    else if (abs(u) < 0.2 && abs(tmax / a) < 0.1)
    {
        dcomp res1 = 0;

        for (int p = 1; p < 10; p++)
        {
            res1 = (dcomp)GG(n, m, p, 0);

            if (tmax != 0)
            {
                for (int q = 1; q < 10; q++)
                {
                    res1 += pow(-((dcomp)tmax) / z0, q) * ((dcomp)GG(n, m, p, q));
                    //std::cout << "ici " << p << " " << q << " " << res1 << " " << pow(-tmax/z0, q) << " " << tmax << " " << z0 << std::endl;
                }
            }

            res1 *= pow(u, p) / ((dcomp)p);
           //std::cout << p << " " << pow(u, p) << std::endl;
            res += res1;
        }

        res *= pow(-1, m) * pow(z0, n - m + 1);
    }


    else if (abs(u) < 0.2 && abs(tmax / a) >= 0.1)
    {
        for (int p = 1; p < 20; p++)
        {
            res += pow(u, p) / ((dcomp)p) * (FF(n, m, p, a * zmax / z0)) * pow(a * zmax, -m + 1) * pow(z0, n);
            //std::cout << "ici " << p  << " " << res << " " << FF(n, m, p, a*zmax/z0) << " " << tmax << " " << z0 << std::endl;
        }

        res = -res;
    }

    return res;
}

/*
dcomp SpecialIntegrations::IntJtComp(int n, int m, double a, double b, double tmin, double tmax)
// Integral of (t^n)/((t-a+ib)^m)
{

    */
/* Remark 
        In the package std::complex, log(std::complex<T> z) is defined with a
        branch cut on the negative axis. Therefore as soon as b is different
        from 0 there are no problem and any integration path will never go
        through the cut. However, when b = 0 one needs to compute specifically 
        the integral in order to be sure to have no problems. When m = 1 one 
        can use IntJt and IntJtComp indifferently however this method is more
        general and allows cases where m != 0 that we may need
    */
/*
    dcomp res = 0;

    double ymin, ymax;
    dcomp zmax, zmin;

    dcomp z0 = a - one_i * b;

    if (b != 0)
    {

        zmax = (tmax - a + one_i * b) / (1. * a);
        zmin = (tmin - a + one_i * b) / (1. * a);
    }
    else
    {
        ymax = tmax / a - 1;
        ymin = tmin / a - 1;

        if (ymin * ymax < 0 || ymin == 0 || ymax == 0)
        {
            std::cout << "FATAL ERROR : IntJt -> Integral is not well defined" << std::endl;
            exit(0);
        }
    }

    //std::cout << "bizarre : " << log(zmax)-log(zmin) << " " << log(zmax/zmin) << std::endl;
    std::cout.precision(15);

    if (b != 0)
    {

        for (int k = 0; k < n + 1; k++)
        {
            if (k - m == -1)
            {
                res += ((dcomp)Binom(k, n)) * pow(z0, n - k) * pow(a, k - m + 1) * (log(zmax/zmin));
                //std::cout << "rrrr : " << zmax << " " << zmin << " " <<  ((dcomp)Binom(k,n))*pow(z0, n-k)*pow(a,k-m+1)*(log(zmax)-log(zmin)) << " " <<  ((dcomp)Binom(k,n))*pow(z0, n-k)*pow(a,k-m+1)*(log(zmax/zmin)) << std::endl;
            }
            else
            {
                res += ((dcomp)Binom(k, n)) * pow(z0, n - k) * pow(a, k - m + 1) * (1. / (k - m + 1)) * (pow(zmax, k - m + 1) - pow(zmin, k - m + 1));
                //std::cout << "ssss : " << (pow(zmax,k-m+1)-pow(zmin,k-m+1)) << " " << ((dcomp)Binom(k,n))*pow(z0, n-k)*pow(a,k-m+1)*(1./(k-m+1))*(pow(zmax,k-m+1)-pow(zmin,k-m+1)) << std::endl;
            }
        }
    }
    else
    {
        //std::cout << ymin << " " << ymax << std::endl;
        for (int k = 0; k < n + 1; k++)
        {
            if (k - m == -1)
                res += ((dcomp)Binom(k, n)) * pow(a, n - m + 1) * (log(ymax / ymin));
            else
                res += ((dcomp)Binom(k, n)) * pow(a, n - m + 1) * (1. / (k - m + 1)) * (pow(ymax, k - m + 1) - pow(ymin, k - m + 1));
        }
    }

    return res;
}*/

dcomp SpecialIntegrations::IntJttComp(int n, int p, int q, double a, double b, double c, double d, double tmin, double tmax)
// Integral of (t^n)/((t-a-ib)^p*(t-c+id)^q)
{
    dcomp res = 0;
    dcomp A, B;

    if ((b != 0 || d != 0) || ((b == 0 && d == 0) && a != c))
    {
        if (q == 1 && p == 1)
        {
            A = pow((a - c) + one_i * (b + d), -1);
            B = pow((c - a) - one_i * (b + d), -1);

            res = A * IntJtComp(n, p, a, -b, tmin, tmax) + B * IntJtComp(n, q, c, d, tmin, tmax);
        }
        else if (q == 0 && p == 1)
            res = IntJtComp(n, p, a, -b, tmin, tmax);
        else if (p == 0 && q == 1)
            res = IntJtComp(n, q, c, d, tmin, tmax);
        else if (p == 0 && q == 0)
            res = IntK(n, tmin, tmax);
        else
        {
            std::cout << "FATAL ERROR : SpecialIntegration::IntJttComp() -> p and q must be equal or lower than 1" << std::endl;
            exit(0);
        }
    }

    else
        res = IntJtComp(n, p + q, a, 0, tmin, tmax);

    return res;
}

// Integral of (t^n)/((t-a-ib)^p*(u-c+id)^q) = (-1)^q*(t^n)/((t-a-ib)^p*(t+r-id)^q)
dcomp SpecialIntegrations::IntJtuComp(int n, int p, int q, double a, double b, double r, double d, double tmin, double tmax)
{
    return pow(-1, q) * IntJttComp(n, p, q, a, b, -r, -d, tmin, tmax);
}

// Integral of (t^n)/((u-a-ib)^p*(u-c+id)^q) = (-1)^(p+q)*(t^n)/((t+r1+ib)^p*(t+r2-id)^q)
dcomp SpecialIntegrations::IntJuuComp(int n, int p, int q, double r1, double b, double r2, double d, double tmin, double tmax)
{
    return pow(-1, p + q) * IntJttComp(n, p, q, -r1, -b, -r2, -d, tmin, tmax);
}

//
//
// Special integrations for the computation of the s-wave and p-wave term
//

dcomp SpecialIntegrations::IntHtComp(int n, int m, double a, double b, double t, double pkl)
{
    return 4 * pkl * pow(t, n) * pow(t - a + one_i * b, -m);
}

dcomp SpecialIntegrations::IntHuComp(int n, int m, double r, double b, double t, double pkl)
{
    return pow(-1, n + 1) * IntHtComp(n, m, r, b, t, pkl);
}

dcomp SpecialIntegrations::IntHttComp(int n, int p, int q, double a, double b, double c, double d, double t, double pkl)
// lim s->(m_i+m_j)^2 of Integral of (1/pij_cm(s))*(t^n)/((t-a-ib)^p*(t-c+id)^q)
{
    dcomp res = 0;
    dcomp A, B;

    if ((b != 0 || d != 0) || ((b == 0 && d == 0) && a != c))
    {
        if (q == 1 && p == 1)
        {
            A = pow((a - c) + one_i * (b + d), -1);
            B = pow((c - a) - one_i * (b + d), -1);

            res = A * IntHtComp(n, p, a, -b, t, pkl) + B * IntHtComp(n, q, c, d, t, pkl);
        }
        else if (q == 0 && p == 1)
            res = IntHtComp(n, p, a, -b, t, pkl);
        else if (p == 0 && q == 1)
            res = IntHtComp(n, q, c, d, t, pkl);
        else if (p == 0 && q == 0)
            res = IntHtComp(n, 0, 1, 1, t, pkl);
        else
        {
            std::cout << "FATAL ERROR : SpecialIntegration::IntJttComp() -> p and q must be equal or lower than 1" << std::endl;
            exit(0);
        }
    }

    else
        res = IntHtComp(n, p + q, a, 0, t, pkl);

    return res;
}

// lim s->(m_i+m_j)^2 of Integral of (t^n)/((t-a-ib)^p*(u-c+id)^q) = (-1)^q*(t^n)/((t-a-ib)^p*(t+r-id)^q)
dcomp SpecialIntegrations::IntHtuComp(int n, int p, int q, double a, double b, double r, double d, double t, double pkl)
{
    return pow(-1, q) * IntHttComp(n, p, q, a, b, -r, -d, t, pkl);
}

// lim s->(m_i+m_j)^2 of Integral of (t^n)/((u-a-ib)^p*(u-c+id)^q) = (-1)^(p+q)*(t^n)/((t+r1+ib)^p*(t+r2-id)^q)
dcomp SpecialIntegrations::IntHuuComp(int n, int p, int q, double r1, double b, double r2, double d, double t, double pkl)
{
    return pow(-1, p + q) * IntHttComp(n, p, q, -r1, -b, -r2, -d, t, pkl);
}
