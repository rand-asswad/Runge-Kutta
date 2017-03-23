/* ========================================================================
 *					   Projet MMSN - Rand ASSWAD
 * ========================================================================
 * 						L'INTERFACE UTILISATEUR
 * ========================================================================
 */

/*** les bibliothèques utilisées ***/
#include <stdio.h>
#include <math.h>

typedef struct point_struct { double x; double y; } point;

/******************************* LES SIGNATURES DES FONCTIONS *******************************/
int choisir_methode ();
int choisir_equation ();
int oui_ou_non (char message[]);
point input_point(int equation);


/*************************************** LES FONCTIONS ***************************************/

int choisir_methode ()
{
	int n;
	puts("Quelle méthode de Runge-Kutta voudriez-vous utiliser ?");
	puts("1) RK12");
	puts("2) RK24");
	scanf("%d",&n);
	while ((n>2)||(n<1))
	{
		puts("Votre choix n\'est pas valable, veuillez entrer 1 ou 2");
		scanf("%d",&n);
	}
	return n;
}

int choisir_equation ()
{
	int n;
	printf("Quelle équation voudriez-vous résoudre ?\n\n");
	printf("1) Modèle de Verhulst (modèle d\'évolution d\'une population de bactérie) - simplifié\n");
	printf("\t y'(t) = y(1-y) \t pour y(t) dans ]0;1[\n\n");
	printf("2) Exemple d\'une EDO du 1er ordre non-linéaire\n");
	printf("\t y'(t) = t*exp(-y)\n\n");
	scanf("%d",&n);
	while ((n>2)||(n<1))
	{
		puts("Votre choix n\'est pas valable, veuillez entrer 1 ou 2");
		scanf("%d",&n);
	}
	return n;
}

int oui_ou_non (char message[])
{
	int n;
	puts(message);
	puts("1) OUI");
	puts("2) NON");
	scanf("%d",&n);
	while ((n>2)||(n<1))
	{
		puts("Votre choix n\'est pas valable, veuillez entrer 1 ou 2");
		scanf("%d",&n);
	}
	return n;
}

point input_point(int equation)
{
	point M;
	puts("Entrez un point M0=(x0,y(x0)) pour trouver la solution");
	printf("x0 = "); scanf("%lf",&M.x);
	printf("y0 = "); scanf("%lf",&M.y);
	if (equation==1) {
		while (!((M.y>0)&&(M.y<1))) {
		puts("L\'equation est définie pour y(t) dans ]0;1[");
		scanf("%lf",&M.y);
		}
	}
	else { if (equation==2) {
		while ((exp(M.y)<(M.x*M.x/2.0))&&(M.y<1)) {
			puts("Pour avoir une solution pour tout t réel il faut que ce point vérifie l'inégalité suivante");
			printf("\t exp(y0) > (t0^2)/2 \t \t (conseil: vérifier si t0 = 0)\n");
			printf("\t Il faut aussi que y0>1\n");
			printf("t0 = "); scanf("%lf",&M.x);
			printf("y0 = "); scanf("%lf",&M.y);
		}
	}}
	return M;
}
