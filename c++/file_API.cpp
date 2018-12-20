#include "file_API.h"

double * getMatrixFromSS(std::stringstream &ss_file_content, int &n) {
  // On procede en suivant le format d'ecriture des matrices definit
  char c{};
  do {
    ss_file_content >> c;
  }
  while(c!='M' and (not ss_file_content.eof()));

  ss_file_content >> n;

  double *M = new double[n*n];
  for(int i=0; i<n*n; i++){
    ss_file_content >> M[i];
  }
  
  return M;
}

int fileToSS(std::stringstream &ss_file_content, const std::string &file_path){
  std::ifstream if_file{file_path};
  if(if_file){
    while(not if_file.eof()){
      char c;
      if_file.get(c);
      ss_file_content << c;
    }
    if_file.close();
    return 1;
  }
  else std::cerr << "ERREUR ! Impossible d'ouvrir le fichier " << file_path << std::endl;
  return 0;
}

void writeResToFile(std::string file_path, const std::vector<int> & n_vector, const std::vector<float> &temps_double_boucle, const std::vector<float> &temps_tile_algorithm){
  
  std::ofstream of_file{file_path,std::ios::trunc};
  if(of_file){
    for(int i=0; i<n_vector.size(); i++)
      of_file << n_vector.at(i) << " " << temps_double_boucle.at(i) << " " << temps_tile_algorithm.at(i) << std::endl;
    
    of_file.close();
  }
  return;
}

