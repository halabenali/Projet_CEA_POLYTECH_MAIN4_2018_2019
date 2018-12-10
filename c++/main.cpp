#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "file_API.h"
#include "cholesky_facto.h"

int main(){

  /*float A[9]={1,0,0,
	      0,1,0,
	      0,0,3};
  float B[9]={1,0,0,
	      0,1,0,
	      0,0,1};
  
	      float C[9]={0,0,0,0,0,0,0,0,0};

  for(int i=0; i<9; i++) {
    std::cout << C[i] << "\t";
    if (i == 2 or i==5) std::cout << std::endl;
  }

  std::cout << std::endl << std::endl;

  cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,3,3,3,1,A,3,B,3,1,C,3);

  for(int i=0; i<9; i++) {
    std::cout << C[i] << "\t";
    if (i == 2 or i==5) std::cout << std::endl;
    }
  */

  std:: stringstream ss_file_content{};

  std::ifstream if_file{"../data/matrices/test"};
  if(if_file){
    std::cout << "good" << std::endl;
    while(not if_file.eof()){
      char c;
      if_file.get(c);
      ss_file_content << c;
    }
    if_file.close();
  }

  int n{};

  double * M = getMatrixFromFile(ss_file_content,n);

  std::cout << "Matrice M : " << std::endl;
  for(int i=0; i<n*n; i++){
    std::cout << M[i] << '\t';
    if((i+1)%n == 0 ) std::cout << std::endl;
  }

  double *L = choleskyLoopsFacto(M, n);

  std::cout << std::endl << std::endl << "Matrice L : " << std::endl;
  for(int i=0; i<n*n; i++){
    std::cout << L[i] << '\t';
    if((i+1)%n == 0 ) std::cout << std::endl;
    }

  choleskyTileFacto(M, n);

  std::cout<< std::endl << "Matrice M cholesky_facto PORTF : " << std::endl;
  for(int i=0; i<n*n; i++){
    std::cout << M[i] << '\t';
    if((i+1)%n == 0 ) std::cout << std::endl;
  }

  delete [] L;
  delete[] M;
  std::cout << std::endl;
  std::cout << std::endl;


  return 0;
}
