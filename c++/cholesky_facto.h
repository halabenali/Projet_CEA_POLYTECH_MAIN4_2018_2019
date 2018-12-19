#include <iostream>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

/* - Retourne L telle que L*L' = A (factorisation de cholesky) par l algorithme classique des boucles imbriquees 
 - L matrice triangulaire inferieur
- A doit etre une matrice carre de taille n, symetrique et definie positive, la fonction ne procede a aucune verification 
- La matrice retournee est allouee dynamiquement (tableau 1D) 
donc ne pas oublier de liberer la memoire apres utilisation ! */
double * choleskyLoopsFacto(const double *A, const int &n);

/* Tile cholesky factorization algorithm
 Modifie la matrice carre A (de taille n) tel que la partie triangulaire inferieur de A est egale a la partie triangulaire inferieur de L 
 avec L la matrice telle que L*L' = A (factorisation de cholesky) par le "tile cholesky factorization algorithm" 
- A doit etre une matrice carre de taille n, symetrique et definie positive, la fonction ne procede a aucune verification
- Retourne 1 si l'algorithme a pu etre execute sans encombre et que L a bien ete calcule, 0 sinon 
 */
int choleskyTileFacto(double *A, const int &n);
