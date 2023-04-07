#ifndef MY_MATH_H
#define MY_MATH_H

#include "../headers/CppLib.h"
#include "../headers/GaussLegendre.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf.h"


extern const double PI2;
extern const double PI4;
extern const double PI8;
extern const double PI16;
extern const double PI;
extern const double pi;

extern const GaussLegendre gl100, gl200, gl500, gl1000, gl2000, gl5000;

// complex variables
typedef std::complex<double> dcomp;
static const dcomp one_i(0.0, 1.0);
static const dcomp one(1.0, 0.0);

/// Define a new type for function pointers (in order to be able of declaring
/// vectors of function pointers).
typedef double (*Function_Pointer)(std::vector<double>);
typedef dcomp (*fcomp_ptr)(double);

// Inline functions ----------------------------------------------------

template <typename T>
T myMax(T a, T b) { return a < b ? b : a; }
template <typename T>
T myMin(T a, T b) { return a > b ? b : a; }
template <typename T>
T myAbs(T a) { return a >= 0 ? a : -a; }
template <typename T>
T mySign(T a)
{
	if (myAbs(a) > 0)
		return a / myAbs(a);
	else
		printf("sign of 0 undef!!! => return 0\n");
	return 0;
}
//else return 1;} // We consider 0 as positive ...
template <typename T>
void myOrderDown(T &a, T &b)
{
	T c;
	if (a > b)
		return;
	else
	{
		c = b;
		b = a;
		a = c;
		return;
	}
}
template <typename T>
void myOrderUp(T &a, T &b)
{
	T c;
	if (a < b)
		return;
	else
	{
		c = b;
		b = a;
		a = c;
		return;
	}
}

template <typename T>
std::vector<T> sorter(std::vector<T> vec)
{
	std::vector<T> res;
	std::vector<T> vec2;

	int jj = 0;

	T val;

	int size;
	bool add_one = true;

	vec2.resize(0);

	for (int i = 0; i < vec.size(); i++)
	{
		for (int j = 0; j < vec2.size(); j++)
		{
			if (vec2[j] == vec[i])
				add_one = false;
		}

		if (add_one == true)
			vec2.push_back(vec[i]);

		add_one = true;
	}

	res.resize(vec2.size());
	size = vec2.size();

	for (int i = 0; i < size; i++)
	{
		res[i] = vec2[0];
		jj = 0;

		for (int j = 0; j < vec2.size(); j++)
		{

			if (vec2[j] < res[i])
			{
				res[i] = vec2[j];
				jj = j;
			}
			//std::cout << "here : " << i << " " << j << " " << jj << " " << vec[j] << " " << res[i] << std::endl;
		}

		// erase the last position that correspond to the min
		/// std::cout << "erase : " << vec[jj] << std::endl;
		vec2.erase(vec2.begin() + jj);
	}

	return res;
}

inline double simpson_coef(int i, int N)
{
	if (i == 0 || i == N)
		return 1. / 3.;
	else
		return 2. * (1. + (i % 2)) / 3.; // gives 4/3,2/3,4/3 ...
}
//-----------------------------------------------------------------------

// Prototypes of normal functions ---------------------------------------

double BesselI0(double x);
double BesselI1(double x);
double BesselK0(double x);
double BesselK1(double x);
double BesselK2(double x);
double BesselK3(double x);
double BesselKn(int n, double x);
double BesselKx(double x, double xx);
double DBesselKn(int n, double x);
double expBesselK0(double x);
double expBesselK1(double x);
double expBesselK2(double x);

double polyLog(unsigned int n, double x);

/// Function that are used to evaluate the one loop-correction of scalar into two vectors
dcomp fOneLoopScalarVectors(double x);
dcomp fOneLoopPseudoScalarVectors(double x);

// Special regulator for the QCD phase transition
double regulatorTQCD(double x, double alpha, double xshift, double max, double min);

double Simpson_Integral1(int ivar, int isize, int ilog, int Nsteps,
						 double xmin, double xmax, std::vector<double> xx,
						 double (*myfunc)(std::vector<double>), double &err);
/* double Simpson_Integral3(int ivar,int isize,int ilog,int Nsteps, */
/* 			 double xmin,double xmax,std::vector<double> xx, */
/* 			 double (*myfunc)(std::vector<double>),double &err); */
double trapezoid_integral(int ivar, int isize,
						  std::vector<double> xx, std::vector<double> interval,
						  double (*myfunc)(std::vector<double>), double &err);
double Simpson_Integral2(int ivar, int isize,
						 std::vector<double> xx, std::vector<double> interval,
						 double (*myfunc)(std::vector<double>), double &err);
double myDerivative_Static(int n, int ivar, std::vector<double> xx, void *pt2Object,
					double (*myfunc)(void * pt2Object, std::vector<double>), double eps = 1.e-3);

// Method too integrate non static (member) functions through wrapper
double Simpson_Integral1_Static(int ivar, int isize, int ilog, int Nsteps,
								   double xmin, double xmax, std::vector<double> xx, void *pt2Object,
								   double (*myfunc)(void *pt2Object, std::vector<double>), double &err);

double GaussLegendre_Integral_Static(int const ivar, int const isize, GaussLegendre const &gl,
									 double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
									 double (*myfunc)(void *pt2Object, std::vector<double>));
double GaussLegendre_IntegralLn_Static(int const ivar, int const isize, GaussLegendre const &gl,
									   double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
									   double (*myfunc)(void *pt2Object, std::vector<double>));
double GaussLegendre_IntegralInverseLn_Static(int const ivar, int const isize, GaussLegendre const &gl,
											  double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
											  double (*myfunc)(void *pt2Object, std::vector<double>));
double GaussLegendre_IntegralLnLn_Static(int const ivar, int const isize, GaussLegendre &gl,
										 double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
										 double (*myfunc)(void *pt2Object, std::vector<double>));
double GaussLegendre_IntegralInverseLnLn_Static(int const ivar, int const isize, GaussLegendre const &gl,
												double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
												double (*myfunc)(void *pt2Object, std::vector<double>));
double GaussLegendre_IntegralExp_Static(int const ivar, int const isize, GaussLegendre const &gl,
										double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
										double (*myfunc)(void *pt2Object, std::vector<double>));
double GSL_Integral_Static_QAGS(void *params, gsl_integration_workspace *w, size_t limit, double xmin, double xmax, std::vector<double> const &xx, void *pt2Object,
								double (*myfunc)(double, void *), double errabs, double errrel, double &err);

double Dichotomie_Static(int ivar, int isize, std::vector<double> xx, void *pt2Object,
						 double (*myfunc)(void *pt2Object, std::vector<double> xx),
						 double const &xMin, double const &xMax, double const &prec, double dflt);
double Dichotomie_Static(void *pt2Object, double (*myfunc)(void *pt2Object, double),
							double const &xMin, double const &xMax, double const &prec, double dlft); // Overriden function for a single parameter
double DichotomieLn_Static(int ivar, int isize, std::vector<double> xx, void *pt2Object,
						   double (*myfunc)(void *pt2Object, std::vector<double> xx),
						   double const &xMin, double const &xMax, double const &prec, double dflt);
double myRungeKutta4_Static(int ivar, double y_in, std::vector<double> xx, double &dx,
							   void *pt2Object, double (*dytodx)(void *pt2Objcet, std::vector<double>, double));

double myRungeKutta4(int ivar, double y_in, std::vector<double> xx, double &dx,
					 double (*dytodx)(std::vector<double>, double));
double myRungeKutta6(int ivar, std::vector<double> x, double y0,
					 void (*derivs)(std::vector<double>, double, double &),
					 bool varstep, double &dx, double dxmin,
					 double prec, bool &okprec);
double myRungeKutta6_Static(int ivar, std::vector<double> x, double y0, void *pt2Object,
							   void (*derivs)(void *pt2Object, std::vector<double>, double, double &),
							   bool varstep, double &dx, double dxmin,
							   double prec, bool &okprec);

// Implicit solver
double ImplicitEulerSolverSquares_Static(int i, std::vector<double> x, double y0, void *pt2Object,
											double (*func_A)(void *pt2Object, std::vector<double>),
											double (*func_Z)(void *pt2Object, std::vector<double>),
											double const &dx);
double ImplicitEulerSolverSquares_Static(int i, std::vector<double> x, double y0, void *pt2Object,
											double (*func_A)(void *pt2Object, std::vector<double>),
											double (*func_B)(void *pt2Object, std::vector<double>),
											double (*func_Z)(void *pt2Object, std::vector<double>),
											double const &dx);
double ImplicitEulerSolverLinear_Static(int i, std::vector<double> x, double y0, void *pt2Object,
										   double (*func_A)(void *pt2Object, std::vector<double>),
										   double (*func_Z)(void *pt2Object, std::vector<double>),
										   double const &dx);
double ImplicitBackwardDiffSolver5Linear_Static(int i, std::vector<double> x, double y0, double ym1, double ym2, double ym3, double ym4, void *pt2Object,
												   double (*func_A)(void *pt2Object, std::vector<double>),
												   double (*func_Z)(void *pt2Object, std::vector<double>),
												   double const &dx);
double ImplicitAdamsMoultonSolver5Linear_Static(int i, std::vector<double> x, std::vector<double> x0, std::vector<double> xm1, std::vector<double> xm2, std::vector<double> xm3,
												   double y0, double ym1, double ym2, double ym3, void *pt2Object,
												   double (*func_A)(void *pt2Object, std::vector<double>),
												   double (*func_Z)(void *pt2Object, std::vector<double>),
												   double const &dx);

double sinc(double x);
double erfcx(double x);
double my_erfc(double x);
double my_erf_diff(double const &b, double const &a);

//----------------------------------------------------------------------

#endif
