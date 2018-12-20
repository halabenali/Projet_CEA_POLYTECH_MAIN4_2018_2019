#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include "file_API.h"
#include "cholesky_facto.h"


//Affiche une matrice carre de taille n stockee en tableau 1D M (matrice de doubles)
void printMatrice(double *M, const int &n){
  for(int i=0; i<n*n; i++){
    std::cout << M[i] << "\t";
    if((i+1)%n == 0) std::cout << std::endl;
  }
}

// le nom du fichier contenant les matrices au format specifie doit etre passse comme unique argument au programme (argv[1])
// et doit se trouver dans le repertoire data/matrices

int main (int argc, char *argv[]) {

  if(argc > 1) {
    
    std:: stringstream ss_file_content{}; // pour stocker le contenu du fichier (les matrices)
    std::string matrices_file_name{argv[1]};
    std::string matrices_file_path{"../data/matrices/" + matrices_file_name};
    std::string temps_res_file_path{"../data/temps/"+matrices_file_name};
      
    if(fileToSS(ss_file_content,matrices_file_path)){

      std::vector<int> n_vector{}; // pour stocker l'ensemble des tailles de matrice existants dans le fichier
      std::vector<float> temps_double_boucle{}; // pour stocker l'ensemble des temps d'execution de choleskyLoopsFacto
      std::vector<float> temps_tile_algorithm{}; // pour stocker l'ensemble des temps d'execution de choleskyTileFacto
      int n{}; // stockage de la taille des matrices
      double * M{}; // stockage des matrices
      double *L{}; // stockage des matrices L de la factorisation de cholesky (L*L'= M), L triangulaire inferieure

      // Pour la mesure des temps d'execution
      std::clock_t begin{};
      std::clock_t end{};
      
      while(not ss_file_content.eof()){ 
	M = getMatrixFromSS(ss_file_content,n); // recuperation de la matrice
	n_vector.push_back(n);
	
	begin = clock();
	L = choleskyLoopsFacto(M, n); // calcul de L par algorithme des doubles boucles
	end = clock();
	temps_double_boucle.push_back( float(end - begin) / CLOCKS_PER_SEC);
	
	begin = clock();
	if(choleskyTileFacto(M, n)){ // calcul de L par algorithme des blocs
	  end = clock();
	  temps_tile_algorithm.push_back( float(end - begin) / CLOCKS_PER_SEC);
	}
	else {
	  temps_tile_algorithm.push_back( -1.0); // si probleme avec choleskyTileFacto(M, n), signalement avec la valeur -1 
	}
	
	delete [] L;
	L = nullptr;
	delete[] M;
	M = nullptr;
      }

      writeResToFile(temps_res_file_path,n_vector,temps_double_boucle,temps_tile_algorithm);
    }
    
  }

  else{
    std::cerr << "Aucun argument n'a ete passe ! Le chemin du fichier contenant les matrices doit etre passe en argument."<< std::endl;
  }
  
  return 0;
}
