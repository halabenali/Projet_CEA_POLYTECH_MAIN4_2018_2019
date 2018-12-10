#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip> 
#include <stdlib.h> 
#include <stdio.h>
#include <time.h>
using namespace std;


//////////////////////////Fonctions//////////////////////////////////////////////////////                                      

float ** cholesky_double_boucle(int n, float** A){
	float ** L= new float*[n];
 	for(int i=0 ; i< n ; ++ i)
	{
			L[i]=new float[n];
	}
	L[0][0]=sqrt(A[0][0]);

	//remplir la première colonne:
	for(int i=1;i<n;i++){
		L[i][0]=A[0][i]/L[0][0];		
	}
	
	//remplir la diagonale et le reste:
	for(int i=1;i<n;i++){
		L[i][i]=A[i][i];
		for(int k=0;k<i-1;k++){
			L[i][i]-=pow(L[i][k],2);
		}
		L[i][i]=sqrt(L[i][i]);
		if(i<n-1){
		
			for(int j=i+1;j<n;j++){
				L[j][i]=A[i][j];
				for(int k=0;k<i-1;k++){
					L[j][i]-=L[i][k]*L[j][k];
					}
				L[j][i]/=L[i][i];
				}
			}
		}
	
	for(int j=1;j<n;j++){
		for(int i=0;i<j;i++)	
			{L[i][j]=0;}
		}
return(L);
}


///////////////////////////////////////////main()//////////////////////////////////////////////////

int main (){
	
	ifstream flux;
	flux.open("MATRIX.txt",fstream::in);
	
	ofstream flux_ecriture;
	flux_ecriture.open("Temps_execution.txt",fstream::out);
	flux_ecriture<<"taille n  temps éxecution/s"<<endl;
	if(flux){
		
	flux.seekg(0, ios::end);  //On se déplace à la fin du fichier et on regarde l'indice du curseur pour connaitre la taille du 						ficher
	int taille_fichier;
	 taille_fichier = flux.tellg();
	flux.seekg(0, ios::beg); // on remet le curseur au début pour commencer à le lire au début 
	
	while(flux.tellg()<=taille_fichier){
	
		if(flux.tellg()==0) {flux.seekg(1, ios::cur);}
		else {flux.seekg(2, ios::cur);}
	
		//récupération de la taille de la matrice sous forme d'un string après la convertir en entier
		string taille;
		flux>>taille;
		char ch[taille.size()+1];
			for(int i=0;i<=taille.size();i++){
				ch[i]=taille[i];
			}
			int n=atoi(ch); //conversion en entier

	
		//allocation de la mémoire pour la matrice A
		float ** A= new float*[n];
	 	for(int i=0 ; i< n ; ++ i)
		{
				A[i]=new float[n];
		}
		//remplissage de la matrice A
		float a;
	
		for (int i=0;i<n;i++)
			{for(int j=0;j<n;j++)
				{flux>>a;
				A[i][j]=a;
				//cout<<A[i][j]<<" ";
				}
			}
	
		//allocation de la mémoire pour la matrice L
		float ** L= new float*[n];
		 	for(int i=0 ; i< n ; ++ i)
			{
					L[i]=new float[n];
			}
		//calcul temps d'éxécution de la factorisation
	 	float temps;
	   	 clock_t t1, t2;
	   	 t1 = clock();
		L=cholesky_double_boucle(n, A);
		t2 = clock();
	   	temps = (float)(t2-t1)/CLOCKS_PER_SEC;

	   	flux_ecriture<<n<<"             "<<temps<<"s"<<endl;

		//libération de la mémoire 
		for(int i=0 ; i< n ; ++ i){
			delete[] A[i] ;
			delete[] L[i] ;
		}
		delete[] A;
		delete[] L;
	}

	flux_ecriture.close();
	flux.close();	
}

	else
	{cout<<"ERREUR:Impossible d'ouvrir le fichier."<<endl;}

	return 0;
}











