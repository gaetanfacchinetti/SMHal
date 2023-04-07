#ifndef DEF_GAUSSLEGENDRE
#define DEF_GAUSSLEGENDRE

#include "CppLib.h"

class GaussLegendre
{
  
 public:
  
  // Constructeur et destructeur
  GaussLegendre(){};
  GaussLegendre(int n_ngl);
  ~GaussLegendre(){};

  // Accesseur
  int const& get_ngl() const {return ngl;};
  double const& get_w(int iii) const {return w[iii];};
  double const& get_x(int iii) const {return x[iii];};
  

 private:
  int ngl;
  double EPS;
  std::vector<double> w,x;

};


#endif