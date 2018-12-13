#include "cholesky_facto.h"

double * choleskyLoopsFacto(const double *A, const int &n){

  double *L = new double[n*n];

  // Initialisation de la partie triangulaire superieure de la matrice L a 0
  // Parcourt ligne par ligne 
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      L[i*n+j] = 0.0;
    }
  }

  
   // Remplissage de la premiere colonne comme l indique l algorithme

  // Initialisation de L00 (premier coefficient de la matrice, premiere ligne et premiere colonne)
  L[0] = sqrt(A[0]);
  // Le reste de la colonne 
  for(int i=1; i<n; i++){
    L[i*n] = A[i]/L[0];
  }
  
  //  Les coefficients restants, colonne par colonne necessairement en suivant l algorithme 
  for(int j=1; j<n; j++){
    // calcul du coefficient Ljj
    double sum{};
    for(int k=0; k<j; k++){
      sum += L[j*n+k]*L[j*n+k];
    }
    L[j*n+j] = sqrt(A[j*n+j] - sum); // Ljj
    
    // calcul des coefficients Lij avec j < i < n 
    for(int i=j+1; i<n; i++){
      sum = 0.0;
      for(int k=0; k<j; k++){
	sum += L[j*n+k]*L[i*n+k];
      }
      L[i*n+j] = (A[j*n+i]-sum)/L[j*n +j]; // Lij
    }
  }

  return L;
}


void choleskyTileFacto(double *A, const int &n){

  lapack_int info;
  double *D = new double[n*n];
  double *C = new double[n*n];
  
  for(int i=0; i<n; i++){

    info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'L',1,&A[i*n+i],1);
    
    if(info == 0) std::cerr << " : info=0, the execution is successful" << std::endl;
    else if(info < 0) std::cerr << " : INFO = -i, the i-th argument had an illegal value\tinfo = " << info << std::endl;
    else  std::cerr << " : INFO = i, the leading minor of order i is not positive definite, and the factorization could not be completed\tinfo = " << info << std::endl;

    for(int t=0; t<n-i-1; t++) D[t] = A[(i+1+t)*n+i];

    if(i<n-1)
      cblas_dtrsm ( CblasRowMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, n-(i+1), 1, 1,A+i*n+i, 1, D, 1 );

    for(int t=0; t<n-i-1; t++) A[(i+1+t)*n+i]= D[t];

    cblas_dsyrk (CblasRowMajor,CblasLower,CblasNoTrans,1,((i+1)),-1,A+((i+1))*n,((i+1)),1,&A[((i+1))*n+((i+1))],1);

    for(int t=0; t<n-((i+1))-1; t++) D[t]= A[(((i+1))+1+t)*n+((i+1))];

    for(int t=0; t<n-((i+1))-1; t++)
      for(int e=0; e<((i+1)); e++)
	C[t*((i+1))+e] = A[(((i+1))+1+t)*n+e];

    if(i<n-1)
      cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasTrans, n-(((i+1))+1), 1, ((i+1)),-1,C, ((i+1)), A+((i+1))*n , ((i+1)), 1, D, 1);
      
    for(int t=0; t<n-((i+1))-1; t++) A[(((i+1))+1+t)*n+((i+1))]= D[t];

   }

  delete[] D;
  D = nullptr;
  delete[] C;
  C = nullptr;

    return;
}
