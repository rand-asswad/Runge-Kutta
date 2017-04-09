/* 
 * =========================================================================
 * 		Resolution of Ordinary Differential Equations
 *		    Runge-Kutta Methods : RK12 et RK24 
 * 				Rand ASSWAD
 * ========================================================================
 *	The program solves an ODE of the form:	y'(t) = f(t,y(t))
 * ========================================================================
 */

/*** USED LIBRARIES ***/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "interface_utilisateur.h"

/*** TYPE DEFINITIONS ***/
typedef double(*eqn_fct)(double,double);	// function type f(t,y(t))
typedef double(*function)(double,point);	// fonction type y(t) passing by the point M0=(t0,y(t0))
typedef double *tab;	// a pointer array of double (dynamic array)
typedef struct res_struct {tab t; tab y; int nb_pts;} points;	// structure of points (t_n,y_n) and number of points

/******************************* FUNCTION SIGNATURES *******************************/

/*** function examples ***/
double f1 (double t, double y);
double sol1 (double t, point M0);
double f2 (double t, double y);
double sol2 (double t, point M0);

/*** Runge-Kutta functions ***/
points RK12(eqn_fct f, point M0, double interval_end, int nb_points, double error_tolerance);
points RK24(eqn_fct f, point M0, double interval_end, int nb_points, double error_tolerance);
#define MAX_pts 100000		// The maximal number of points for solving the ODE
#define beta	0.8		// The step correction factor (if erreur>tolérance)

/*** Result display functions ***/
void displaySolution (points solution_points, function y, point M0);
#define nb_digits 6		// Number of digits to display (after the point)

/*** Result manipulation functions ***/
void extractPoints (char fileName[], points solution_points, function exact_solution, point M0);	// Extract points in a file to compare results with exact solution
void gnuplot(char fileName[]);	// Plot the points extracted in a file

/*************************************** LES FONCTIONS ***************************************/

// Dynamic array functions
tab init_tab (int n) { return malloc(sizeof(double)*n); }	// initialize array
tab copy_tab (int n, double old_tab[])
{
	tab t = init_tab (n);
	for (int i=0; i<n; i++) t[i] = old_tab[i];
	return t;
}

void displaySolution (points c, function y, point M0)
{
	printf("k \t t_k \t\t y_k \t\t y(t_k) \t error\n");
	for (int k=0; k<c.dim; k++) printf("%d \t %5.*f \t %5.*f \t %5.*f \t %e\n",k,nb_digits,c.t[k],nb_digits,c.y[k],nb_digits,y(c.t[k],M0),fabs(c.y[k]-y(c.t[k],M0)));
}

points RK12(eqn_fct f, point M0, double b, int n, double tol)
{
	double t[MAX_pts], y[MAX_pts], tk,yk,k1,k2,h,eps;
	int k=0;
	
	h = (b-M0.t)/n;	// initialize h
	t[0] = M0.t;	// t0 = a
	y[0] = M0.y;	// y0 = y(x0) = y(a)
	
	while (t[k]<=b)
	{
		if (h<DBL_EPSILON)
		/* If the error tolerance is too small, it might need an h smaller that the double precision
		 * in that case, we resolve the problem by increasing the error tolerance by a factor of 10 */
		{
			h = (b-M0.t)/n;		// re-initalize h
			eps *= 10;		// increase the error tolerance
			k = 0;			// restart algorithme
		}
		
		tk = t[k] + h;
		
		k1 = f(t[k],y[k]);		// 1st order
		k2 = f(tk, y[k]+h*k1);		// 2nd order
		yk = y[k] + h*(k1 + k2)/2.0;	// RK12 formula
		
		eps = h*fabs(k2-k1)/2.0;	// consistance error(k) = h*|phi_2 - phi_1|
		
		// Step correction
		if (eps >= tol) h *= beta*sqrt(tol/eps);
		else
		{
			h *= beta*(tol/eps);

			k++;
			x[k] = xk;
			y[k] = yk;
		}
	}

	points curve;
	curve.t = copy_tab(k,t);
	curve.y = copy_tab(k,y);
	curve.nb_pts = k;
	
	return curve;
}

points RK24(eqn_fct f, point M0, double b, int n, double tol)
{
	double t[MAX_pts], y[MAX_pts], tk,yk,k1,k2,k3,k4,h,mid_h,eps;
	int k=0;
	
	h = (b-M0.t)/n;	// initialiser h
	t[0] = M0.t;	// t0 = a
	y[0] = M0.y;	// y0 = y(t0) = y(a)
	
	while (t[k]<=b)
	{
		if (h<DBL_EPSILON)
		/* If the error tolerance is too small, it might need an h smaller that the double precision
		 * in that case, we resolve the problem by increasing the error tolerance by a factor of 10 */
		{
			h = (b-M0.t)/n;		// re-initalize h
			eps *= 10;		// increase the error tolerance
			k = 0;			// restart algorithme
		}
		
		mid_h = h/2.0;
		tk = t[k] + mid_h;
		
		// RK2
		k1 = f(t[k],y[k]);
		k2 = f(tk, y[k] + mid_h*k1);
		
		// RK4
		k3 = f(tk, y[k] + mid_h*k2);
		tk += mid_h;
		k4 = f(tk, y[k] + h*k3);
		
		// RK24
		yk = y[k] + h*(k1 + 2*(k2 + k3) + k4)/6.0;
		
		eps = h*fabs(k1 - k2 + 2*k3 + k4)/6.0;		// consistance error(k) = h*|phi_2 - phi_1|
		
		// Step correction
		if (eps >= tol) h *= beta*pow(tol/eps,0.25);
		else
		{
			h *= beta*sqrt(tol/eps);

			k++;
			x[k] = xk;
			y[k] = yk;
		}
	}

	points curve;
	curve.t = copy_tab(k,t);
	curve.y = copy_tab(k,y);
	curve.nb_pts = k;
	
	return curve;
}

#define STR_MAX 100
void extractPoints (char fileName[], points sol, function exact, point M0)
{
	char path[STR_MAX] = "./files/";
	strcat(path,fileName);
	strcat(path,".dat");
	
	FILE* f = NULL;
	
	f = fopen(path, "w");
	for (int i=0; i<sol.nb_pts; i++)
		fprintf(f,"%f %f %f\n", sol.t[i], sol.y[i], exact(points.t[i],M0));
	
	fclose(f);
}

void gnuplot(char fileName[])
{
	char path[STR_MAX] = "./files/";
	strcat(path,fileName);
	strcat(path,".dat");
	
	FILE * gnuplotPipe = popen("gnuplot -persistent","w");
	fprintf(gnuplotPipe, "set title \"Runge-Kutta Methods\"\n");
	fprintf(gnuplotPipe, "plot '%s' using 1:2 title \"Runge-Kutta\" with points, '' using 1:3 title \"Exact Solution\"\n",path);
}

int main()
{
	puts("***** Resolution of Ordinary Differential Equations *****");
	puts("*****	 Runge-Kutta Methodes : RK12 & RK24	  *****");
	puts("*****		     Rand ASSWAD		  *****");
	puts("");
	
	int methode, equation, parameters, nb_pts = 100;
	function sol;
	eqn_fct f;
	point M;
	double tol=0.001, b=20.0;
	points results;
	char fileName[STR_MAX], c[]="0";
	
	methode = chooseMethode();
		
	equation = chooseEquation();
	if (equation==1) {
		f = &f1;	sol = &sol1;
		M.t = 0.0;	M.y = 0.5;
	}
	else {
		f = &f2;	sol = &sol2;
		M.t = 0.0;	M.y = 0.0;
	}
	
	parameters = YES_NO("Would you like to enter the methode parameters ? (if not, we assign parameters by default)");
	
	if (parameters==1) {
		M = input_point(equation);
		puts("Please enter the end point of the interval (a value greater than t0)");
		do scanf("%lf",&b); while (b<=M.t);
		puts("Please enter the maximal tolerated error"); scanf("%lf",&tol); tol=fabs(tol);
		puts("Please enter the initial number of points you would like to use"); scanf("%d",&nb_pts); nb_pts=abs(nb_pts);
	}
	
	switch (methode) {
		case 1 :
			results = RK12(f,M,b,nb_pts,tol);
			strcpy(fileName,"RK12-eq");
			break;
		case 2 :
			results = RK24(f,M,b,nb_pts,tol);
			strcpy(fileName,"RK24-eq");
	}
	
	c[0] += equation;
	strcat(fileName,c);
	
	extractPoints (fileName,results,sol,M);
	afficher_solution(results,sol,M);
	gnuplot(fileName);
	
	return EXIT_SUCCESS;
}





// first example : y'(t) = y(1-y) ==> y(t) = 1/(1 + c*exp(-t))
double f1 (double t, double y)
{
	return y*(1-y);
}
double sol1 (double t, point M0)
// ATTENTION!!! y0 = M0.y cannot be equal to (0)
{
	double c = (1.0/M0.y - 1.0)*exp(M0.t);
	return 1.0/(1.0+exp(-t)*c);
}

// second example : y'(t) = t*exp(-y) ==> y(t) = ln(t²/2 + c)
double f2 (double t, double y)
{
	return t*exp(-y);
}
double sol2 (double t, point M0)
{
	double c = exp(M0.y) - M0.x*M0.x/2.0;
	return log(t*t/2.0+c);
}
