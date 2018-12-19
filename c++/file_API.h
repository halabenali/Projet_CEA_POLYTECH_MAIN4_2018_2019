#include <iostream>
#include <sstream>

// retourne la premiere matrice lu dans ss_file_content a partir de la position courante du curseur dans la stringstream
// la taille de la matrice est inscrite dans n
// retourne la matrice sous forme de tableau 1D alloue dynamiquement donc ne pas oublier de liberer la memoire apres utilisation
// le format des matrices dans la stringstream doit suivre le format specifie dans le fichier README.md  
double * getMatrixFromFile(std::stringstream &ss_file_content, int &n);
