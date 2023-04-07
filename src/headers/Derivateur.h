#ifndef DEF_DERIVATEUR
#define DEF_DERIVATEUR


#include "CppLib.h"


class Derivateur
{
    
public:
    
    // Constructeur et destructeurs
    Derivateur(){};
    ~Derivateur(){};
    
    // Fonction de classe
    double der1bulk(double, double, double,double, double);
    double der1FstPt(double, double, double,double, double);
    double der1ScdPt(double, double, double,double, double);
    double der1LastPtm1(double, double, double,double, double);
    double der1LastPt(double, double, double,double, double);
    
    double der2bulk(double, double, double,double, double);
    double der2FstPt(double, double, double,double, double);
    double der2ScdPt(double, double, double,double, double);
    double der2LastPtm1(double, double, double,double, double);
    double der2LastPt(double, double, double,double, double);
    
    double derivePrm(std::vector<double> const& fff, int iii, int nbPtRho);
    double deriveScd(std::vector<double> const& fff, int iii, int nbPtRho);
    
    
};

#endif
