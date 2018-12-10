#include "file_API.h"

double * getMatrixFromFile(std::stringstream &ss_file_content, int &n) {
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
