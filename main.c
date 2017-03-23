/* 
 * =========================================================================
 * 		Projet MMSN : Résolution d'équations differentielles ordinaires
 *				   Méthode de Runge-Kutta : RK12 et RK24 
 * 				 		Fait par : Rand ASSWAD
 * ========================================================================
 */

/*** les bibliothèques utilisées ***/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "interface_utilisateur.h"

/*** définitions des types ***/
typedef double(*eqn_fct)(double,double);	// type correspondant à la fonction f(x,y(x))
typedef double(*fonction)(double,point);	// type correspondant à la fonction y(x) pour un point M0=(x0,y(x0))
typedef double *tab;	// un pointeur tableau (tableau dynamique) de types double
typedef struct res_struct {tab x; tab y; int dim;} res;		// enregistrement les xk,yk, et le nombre total de points

/******************************* LES SIGNATURES DES FONCTIONS *******************************/

/*** les fonctions exemples ***/
double f1 (double t, double y);
double sol1 (double t, point M0);
double f2 (double t, double y);
double sol2 (double t, point M0);

/*** fonctions de Runge-Kutta ***/
res RK12(eqn_fct f, point M0, double fin_intervale, int nb_points, double tolerance);
res RK24(eqn_fct f, point M0, double fin_intervale, int nb_points, double tolerance);
#define MAX_pts 100000		// le nombre de points maximal à prendre pour la résolution de l'EDO
#define beta 0.8			// coéficient de correction du pas (si erreur>tolérance)

/*** fonction d'affichage de résultats ***/
void afficher_solution (res points_courbe, fonction y, point M0);
#define nb_chiffres	6		// nombre de chiffres à afficher après la virgule

/*** fonctions de traitement de résultats***/
void extraire_points (char nomFichier[], res points, fonction sol_exacte, point M0);	// extraire les résultats dans un fichier, avec la solution exacte (pour comparer)
void gnuplot(char nomFichier[]);	// dessiner les points sur gnuplot à partir d'un fichier

/*************************************** LES FONCTIONS ***************************************/

// fonctions de traitement de tableaux dynamiques
tab init_tab (int n) { return malloc(sizeof(double)*n); }
tab copy_tab (int n, double old_tab[])
{
	tab t = init_tab (n);
	for (int i=0; i<n; i++) t[i] = old_tab[i];
	return t;
}

void afficher_solution (res c, fonction y, point M0)
{
	printf("k \t x_k \t\t y_k \t\t y(x_k) \t erreur\n");
	for (int k=0; k<c.dim; k++) printf("%d \t %5.*f \t %5.*f \t %5.*f \t %e\n",k,nb_chiffres,c.x[k],nb_chiffres,c.y[k],nb_chiffres,y(c.x[k],M0),fabs(c.y[k]-y(c.x[k],M0)));
}

res RK12(eqn_fct f, point M0, double b, int n, double tol)
{
	double x[MAX_pts], y[MAX_pts], xk,yk,k1,k2,h,eps;
	int k=0;
	
	h = (b-M0.x)/n;	// initialiser h
	x[0] = M0.x;	// x0 = a
	y[0] = M0.y;	// y0 = y(x0) = y(a)
	
	while (x[k]<=b)
	{
		if (h<2E-50) {
			h = (b-M0.x)/n;	// réinitialiser h
			eps *= 10;		// augmenter la tolérance
			k = 0;			// recommencer
		}
		
		xk = x[k] + h;
		
		k1 = f(x[k],y[k]);		// ordre 1
		k2 = f(xk, y[k]+h*k1);	// ordre 2
		yk = y[k] + h*(k1 + k2)/2.0;	// RK12
		
		eps = h*fabs(k2-k1)/2.0;		// erreur(k) = h*|phi_2 - phi_1|
		
		if (eps >= tol) h *= beta*sqrt(tol/eps);
		else
		{
			h *= beta*(tol/eps);

			k++;
			x[k] = xk;
			y[k] = yk;
		}
	}

	res courbe;
	courbe.x = copy_tab(k,x);
	courbe.y = copy_tab(k,y);
	courbe.dim = k;
	
	return courbe;
}

res RK24(eqn_fct f, point M0, double b, int n, double tol)
{
	double x[MAX_pts], y[MAX_pts], xk,yk,k1,k2,k3,k4,h,mid_h,eps;
	int k=0;
	
	h = (b-M0.x)/n;	// initialiser h
	x[0] = M0.x;	// x0 = a
	y[0] = M0.y;	// y0 = y(x0) = y(a)
	
	while (x[k]<=b)
	{
		if (h<2E-50) {
			h = (b-M0.x)/n;	// réinitialiser h
			eps *= 10;		// augmenter la tolérance
			k = 0;			// recommencer
		}
		
		mid_h = h/2.0;
		xk = x[k] + mid_h;
		
		// Ordre 2
		k1 = f(x[k],y[k]);
		k2 = f(xk, y[k] + mid_h*k1);
		
		// Ordre 4
		k3 = f(xk, y[k] + mid_h*k2);
		xk += mid_h;
		k4 = f(xk, y[k] + h*k3);
		
		// RK24
		yk = y[k] + h*(k1 + 2*(k2 + k3) + k4)/6.0;
		
		eps = h*fabs(k1 - k2 + 2*k3 + k4)/6.0;
		
		if (eps >= tol) h *= beta*pow(tol/eps,0.25);
		else
		{
			h *= beta*sqrt(tol/eps);

			k++;
			x[k] = xk;
			y[k] = yk;
		}
	}

	res courbe;
	courbe.x = copy_tab(k,x);
	courbe.y = copy_tab(k,y);
	courbe.dim = k;
	
	return courbe;
}

#define CHAINE_MAX 100
void extraire_points (char nomFichier[], res points, fonction sol_theorique, point M0)
{
	char path[CHAINE_MAX] = "./fichiers/";
	strcat(path,nomFichier);
	strcat(path,".dat");
	
	FILE* fichier = NULL;
	
	fichier = fopen(path, "w");
	for (int i=0; i<points.dim; i++)
		fprintf(fichier,"%f %f %f\n", points.x[i], points.y[i], sol_theorique(points.x[i],M0));
	
	fclose(fichier);
}

void gnuplot(char nomFichier[])
{
	char path[CHAINE_MAX] = "./fichiers/";
	strcat(path,nomFichier);
	strcat(path,".dat");
	
	FILE * gnuplotPipe = popen("gnuplot -persistent","w");
	fprintf(gnuplotPipe, "set title \"Méthode de Runge-Kutta\"\n");
	fprintf(gnuplotPipe, "plot '%s' using 1:2 title \"Runge-Kutta\" with points, '' using 1:3 title \"Solution Exacte\"\n",path);
}

int main()
{
	puts("***** Projet MMSN : Résolution d'équations differentielles ordinaires *****");
	puts("*****             Méthode de Runge-Kutta : RK12 et RK24               *****");
	puts("*****        Projet realise par Rand ASSWAD et Cedric FRAGNAUD        *****");
	puts("");
	
	int methode, equation, parametres, nb_pts = 100;
	fonction sol;
	eqn_fct f;
	point M;
	double tol=0.001, b=20.0;
	res resultats;
	char nomFichier[CHAINE_MAX], c[]="0";
	
	methode = choisir_methode();
		
	equation = choisir_equation();
	if (equation==1) {
		f = &f1;	sol = &sol1;
		M.x = 0.0;	M.y = 0.5;
	}
	else {
		f = &f2;	sol = &sol2;
		M.x = 0.0;	M.y = 0.0;
	}
	
	parametres = oui_ou_non("Voulez-vous entrez les paramètres de la méthode ? (sinon on vous propose des valeurs par défaut)");
	
	if (parametres==1) {
		M = input_point(equation);
		puts("Entrez la borne supérieure de votre intervale (une valeure supérieure à x0)");
		do scanf("%lf",&b); while (b<=M.x);
		puts("Entrez l'erreur maximale tolérée"); scanf("%lf",&tol); tol=fabs(tol);
		puts("Entrez le nombre de points à utiliser au départ"); scanf("%d",&nb_pts); nb_pts=abs(nb_pts);
	}
	
	switch (methode) {
		case 1 :
			resultats = RK12(f,M,b,nb_pts,tol);
			strcpy(nomFichier,"RK12-eq");
			break;
		case 2 :
			resultats = RK24(f,M,b,nb_pts,tol);
			strcpy(nomFichier,"RK24-eq");
	}
	
	c[0] += equation;
	strcat(nomFichier,c);
	
	extraire_points (nomFichier,resultats,sol,M);
	afficher_solution(resultats,sol,M);
	gnuplot(nomFichier);
	
	return EXIT_SUCCESS;
}





// premier exemple
double f1 (double t, double y)
{
	return y*(1-y);
}
double sol1 (double t, point M0)
// ATTENTION!!! y0 = M0.y ne peux pas être égale à (0)
{
	double c = (1.0/M0.y - 1.0)*exp(M0.x);
	return 1.0/(1.0+exp(-t)*c);
}

// deuxième exemple
double f2 (double t, double y)
{
	return t*exp(-y);
}
double sol2 (double t, point M0)
{
	double c = exp(M0.y) - M0.x*M0.x/2.0;
	return log(t*t/2.0+c);
}
