/* ==================================================================
 *	ODE solver using Runge-Kutta methods  -  by Rand ASSWAD
 * ==================================================================
 * 			    USER INTERFACE
 * ==================================================================
 */

/*** USED LIBRARIES ***/
#include <stdio.h>
#include <math.h>

typedef struct point_struct { double t; double y; } point;

/******************************* LES SIGNATURES DES FONCTIONS *******************************/
int chooseMethode ();
int chooseEquation ();
int YES_NO (char message[]);
point input_point(int equation);


/*************************************** LES FONCTIONS ***************************************/

int chooseMethode ()
{
	int n;
	puts("Which Runge-Kutta methode would you like to use?");
	puts("1) RK12");
	puts("2) RK24");
	scanf("%d",&n);
	while ((n>2)||(n<1))
	{
		puts("Your choice is not possible, please enter 1 or 2");
		scanf("%d",&n);
	}
	return n;
}

int chooseEquation ()
{
	int n;
	printf("Which ODE would you like to solve?\n\n");
	printf("1) Verhulst's model (bacteria population evolution model) - simplified\n");
	printf("\t y'(t) = y(1-y) \t for y(t) in ]0;1[\n\n");
	printf("2) Example of a non-linear first order ODE\n");
	printf("\t y'(t) = t*exp(-y)\n\n");
	scanf("%d",&n);
	while ((n>2)||(n<1))
	{
		puts("Your choice is not possible, please enter 1 or 2");
		scanf("%d",&n);
	}
	return n;
}

int YES_NO (char message[])
{
	int n;
	puts(message);
	puts("1) YES");
	puts("2) NO");
	scanf("%d",&n);
	while ((n>2)||(n<1))
	{
		puts("Your choice is not possible, please enter 1 or 2");
		scanf("%d",&n);
	}
	return n;
}

point input_point(int equation)
{
	point M;
	puts("Please enter a point M0=(t0,y(t0)) to find y(t)");
	printf("t0 = "); scanf("%lf",&M.t);
	printf("y0 = "); scanf("%lf",&M.y);
	if (equation==1) {
		while (!((M.y>0)&&(M.y<1))) {
		puts("The equation is defined for y(t) in ]0;1[");
		scanf("%lf",&M.y);
		}
	}
	else { if (equation==2) {
		while ((exp(M.y)<(M.t*M.t/2.0))&&(M.y<1)) {
			puts("To obtain a solution for all real t the initial point M0 has to verify the following inequality");
			printf("\t exp(y0) > (t0^2)/2 \t \t (hint: check for t0 = 0)\n");
			printf("\t also y0>1\n");
			printf("t0 = "); scanf("%lf",&M.t);
			printf("y0 = "); scanf("%lf",&M.y);
		}
	}}
	return M;
}
