#include <iostream>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

/* - Retourne L tel que L*L' = A (factorisation de cholesky) par l algorithme classique des boucles imbriquees 
 - L matrice triangulaire inferieur
- A doit etre une matrice carre de taille n, symetrique et definie positive, la fonction ne procede a aucune verification 
- La matrice retournee est allouee dynamiquement (tableau 1D) 
donc ne pas oublier de liberer la memoire apres utilisation  */
double * choleskyLoopsFacto(const double *A, const int &n);

void choleskyTileFacto(double *A, const int &n);
