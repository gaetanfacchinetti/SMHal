#include "../headers/Derivateur.h"


//1st Derivative
double Derivateur::der1bulk(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return (fm2/4. - 2*fm1 + 2*fp1 - fp2/4.)/3.;
}

double Derivateur::der1FstPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return -25*fm2/12. + 4*fm1 - 3*f0 + 4*fp1/3. - fp2/4.;
}

double Derivateur::der1ScdPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return (-fm2/2. - 5*fm1/3. + 3*f0 - fp1 + fp2/6.)/2.;
}
double Derivateur::der1LastPtm1(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return (fp2/2. + 5*fp1/3. - 3*f0 + fm1 - fm2/6.)/2.;
}
double Derivateur::der1LastPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return 25*fp2/12. - 4*fp1 + 3*f0 - 4*fm1/3. + fm2/4. ;
}



//2nd Derivative
double Derivateur::der2bulk(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return (-fm2/4. + 4*fm1 + 4*fp1 - fp2/4.)/3. - 5*f0/2.;
}
double Derivateur::der2FstPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return 35*fm2/12. - 26*fm1/3. + 19*f0/2. - 14*fp1/3. + 11*fp2/12.;
}
double Derivateur::der2ScdPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return (11*fm2)/12. - (5*fm1)/3. + f0/2. + fp1/3. - fp2/12.;
}
double Derivateur::der2LastPtm1(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return -fm2/12. + fm1/3. + f0/2. - (5*fp1)/3. + (11*fp2)/12.;
}
double Derivateur::der2LastPt(double fm2, double fm1, double f0, double fp1, double fp2)
{
    return 11*fm2/12. - 14*fm1/3. + 19*f0/2. - 26*fp1/3. + 35*fp2/12.;
}



// First derivative calculation
double Derivateur::derivePrm(std::vector<double> const& fff, int iii, int nbPtRho)
{
    double fm2, fm1, f0, fp1, fp2;
    if(iii>=2&&iii<(nbPtRho-1))
    {
        fm2 = fff[iii-2];
        fm1 = fff[iii-1];
        f0 = fff[iii];
        fp1 = fff[iii+1];
        fp2 = fff[iii+2];
        return der1bulk(fm2,fm1,f0,fp1,fp2);
    }
    else if(iii==0)
    {
        fm2 = fff[0];
        fm1 = fff[1];
        f0 = fff[2];
        fp1 = fff[3];
        fp2 = fff[4];
        return der1FstPt(fm2,fm1,f0,fp1,fp2);
    }
    else if(iii==1)
    {
        fm2 = fff[0];
        fm1 = fff[1];
        f0 = fff[2];
        fp1 = fff[3];
        fp2 = fff[4];
        return der1ScdPt(fm2,fm1,f0,fp1,fp2);
    }
    else if(iii==(nbPtRho-1))
    {
        fm2 = fff[nbPtRho-4];
        fm1 = fff[nbPtRho-3];
        f0 = fff[nbPtRho-2];
        fp1 = fff[nbPtRho-1];
        fp2 = fff[nbPtRho];
        return der1LastPtm1(fm2,fm1,f0,fp1,fp2);
    }
    else if(iii==nbPtRho)
    {
        fm2 = fff[nbPtRho-4];
        fm1 = fff[nbPtRho-3];
        f0 = fff[nbPtRho-2];
        fp1 = fff[nbPtRho-1];
        fp2 = fff[nbPtRho];
        return der1LastPt(fm2,fm1,f0,fp1,fp2);
    }
    exit(0);
}




// Second derivative calculation
double Derivateur::deriveScd(std::vector<double> const& fff, int iii, int nbPtRho)
{
    double fm2, fm1, f0, fp1, fp2;
    if(iii>=2&&iii<(nbPtRho-1))
    {
        fm2 = fff[iii-2];
        fm1 = fff[iii-1];
        f0 = fff[iii];
        fp1 = fff[iii+1];
        fp2 = fff[iii+2];
        return der2bulk(fm2,fm1,f0,fp1,fp2);
    }
    else if(iii==0)
    {
        fm2 = fff[0];
        fm1 = fff[1];
        f0 = fff[2];
        fp1 = fff[3];
        fp2 = fff[4];
        return der2FstPt(fm2,fm1,f0,fp1,fp2);
    }
    else if(iii==1)
    {
        fm2 = fff[0];
        fm1 = fff[1];
        f0 = fff[2];
        fp1 = fff[3];
        fp2 = fff[4];
        return der2ScdPt(fm2,fm1,f0,fp1,fp2);
    }
    else if(iii==(nbPtRho-1))
    {
        fm2 = fff[nbPtRho-4];
        fm1 = fff[nbPtRho-3];
        f0 = fff[nbPtRho-2];
        fp1 = fff[nbPtRho-1];
        fp2 = fff[nbPtRho];
        return der2LastPtm1(fm2,fm1,f0,fp1,fp2);
    }
    else if(iii==nbPtRho)
    {
        fm2 = fff[nbPtRho-4];
        fm1 = fff[nbPtRho-3];
        f0 = fff[nbPtRho-2];
        fp1 = fff[nbPtRho-1];
        fp2 = fff[nbPtRho];
        return der2LastPt(fm2,fm1,f0,fp1,fp2);
    }
    exit(0);
}
