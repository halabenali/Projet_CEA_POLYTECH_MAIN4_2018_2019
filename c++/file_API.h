#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// retourne la premiere matrice lu dans ss_file_content a partir de la position courante du curseur dans la stringstream
// la taille de la matrice est inscrite dans n
// retourne la matrice sous forme de tableau 1D alloue dynamiquement donc ne pas oublier de liberer la memoire apres utilisation
// le format des matrices dans la stringstream doit suivre le format specifie dans le fichier README.md  
double * getMatrixFromSS(std::stringstream &ss_file_content, int &n);

/*
  Ecrit le contenu du fichier file_path dans la stringstream ss_file_content
  Le contenu de ss_file_content n'est pas efface ! Le contenu du fichier est ecrit a la suite du contenu de la stringstream ! 
  Retourne 1 si le fichier a pu etre ouvert et que l'operation a ete un succes, sinon 0  
 */
int fileToSS(std::stringstream &ss_file_content, const std::string &file_path);

/*
  Ecrit les temps d'execution mesures des algorithmes de factorisation de cholesky dans le fichier dont le chemin est file_path
  Le contenu du fichier est ecrase si il existe deja. 
  Ecrit de la maniere suivant ligne par ligne : valeur_de_n temps_double_boucle temps_tile_algorithme
 */
void writeResToFile(std::string file_path, const std::vector<int> & n_vector, const std::vector<float> &temps_double_boucle, const std::vector<float> &temps_tile_algorithm);
