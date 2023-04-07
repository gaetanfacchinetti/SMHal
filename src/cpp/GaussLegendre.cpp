#include "../headers/GaussLegendre.h"


 
// Calculation of Gauss Legendre weight and abscissas
GaussLegendre::GaussLegendre(int n_ngl)
{

  ngl = n_ngl;

  EPS = 3.0e-14;

  int m,i,j;
  double z1,z,pp,p3,p2,p1; //High precision is a good idea for this routine.

  x.resize(ngl+2);
  w.resize(ngl+2);

    // The roots are symmetric in the interval, so we only have to find half of them
  m=(ngl+1)/2;

  
  for(i=1;i<=m;i++)
    // Loop over the desired roots
    {
      // Starting with the above approximation the ith root,
      // we enter the main loop of refinement by Newtow's method
      z = cos(M_PI*(i-0.25)/(ngl+0.5));
      do
        {
	  p1=1.0;
	  p2=0.0;

	  // Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
	  for(j=1;j<=ngl;j++)
            {
	      p3=p2;
	      p2=p1;
	      p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
	  //p1 is now the desired Legendre polynomial.
	  //We next compute pp, its derivative, by a standard relation
	  //involving also p2, the polynomial of one lower order.
	  pp=ngl*(z*p1-p2)/(z*z-1.0);
	  z1=z;
	  //Newthon's method.
	  z=z1-p1/pp;       
        }while(fabs(z-z1)>EPS);
      
      //Scale the root to the desired interval, and put in its symmetric counterpart.
      x[i]=-z;                       
      x[ngl+1-i]=z;
      
      //Compute the weight and its symmetric counterpart
      w[i]=2./((1.0-z*z)*pp*pp);      
      w[ngl+1-i]=w[i];
    }

/*
  for(int i = 0; i < ngl+2; i++)
  {
    std::cout << i << " " << x[i] << " " << w[i] << std::endl;
  } */
 
  //std::cout << " | Initialisation GaussLegendre complete " << std::endl;
  //std::cout << " | Nombre de Points : " << ngl << std::endl;
  //std::cout << std::endl; 
 
}

