/*
 * matutils.c
 *
 *  Created on: 1 abr 2021
 *      Author: Administrador
 */


#include "matutils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

double** zeros(int n, int m) {
	double **mat = (double**) calloc(m, sizeof(double*));

	if (mat == NULL){
		printf("No se pudo reservar memoria.");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < m; i++) {
		mat[i] = (double*) calloc(n, sizeof(double));
		if (mat[i] == NULL){
			printf("No se pudo reservar memoria.");
			exit(EXIT_FAILURE);
		}
	}

	// calloc inicializa a ceros, así que nos ahorramos ese código.

	return mat;
}

int compare(double* A, double* B, int n){
	int res = 1;

	for(int i = 0; i<n; i++)
		if(A[i] != B[i])
			res = 0;

	return res;
}

int compare_mat(double** A, double** B, int n, int m){
	int res = 1;

	for(int i = 0; i<n; i++)
		for(int j = 0; j<m; j++)
			if(A[j][i] != B[j][i])
				res = 0;

	return res;
}

double norm(double* v, int n){
	double res = 0.0;

	for(int i = 0; i<n; i++)
		res += v[i]*v[i];

	return sqrt(res);
}

double dot(double* u, double* v, int n){
	double res = 0.0;

	for(int i = 0; i < n; i++)
		res += u[i]*v[i];

	return res;
}

double* unit(double* v, int n){
	double* res = (double*)calloc(n, sizeof(double));

	double norma = norm(v,n);

	if(norma > DBL_EPSILON) {
		for(int i = 0; i<n; i++)
			res[i] = v[i]/norma;
	}
	else {
		for(int i = 0; i<n; i++)
			res[i] = 0.0;
	}

	return res;
}

double** producto(double** A, double** B, int n, int m, int p){
	double** prod = zeros(n,p);

	for(int i = 0; i<n; i++){
		for(int j = 0; j<p; j++){
			for(int k = 0; k<m;k++){
				prod[j][i] += A[k][i]*B[j][k];
			}
		}
	}

	return prod;
}

double** producto_escalar_mat(double** A, int n, double k){
	double** prod = zeros(n,n);

	for(int i = 0; i<n; i++){
		for(int j = 0; j<n; j++){
			prod[j][i] = k*A[j][i];
		}
	}

	return prod;
}

double* producto_escalar(double *v, int n, double k) {
	double* prod = (double*)calloc(n,sizeof(double));

	for (int i = 0; i < n; i++)
		prod[i] = k*v[i];

	return prod;
}

double** vect_a_matr(double* v){
	double** res = (double**) calloc(1, sizeof(double*));

	if(res == NULL)
		exit(EXIT_FAILURE);

	res[0] = v;

	return res;
}
double** transpose(double** A, int n, int m){
	double** T = zeros(m,n);

	for(int i = 0;i<n;i++)
		for(int j = 0;j<m;j++)
			T[i][j] = A[j][i];

	return T;
}

double** suma(double** A, double** B, int n, int m){
	double** S = zeros(n,m);

	for(int i = 0;i<m;i++)
			for(int j = 0;j<n;j++)
				S[i][j] = A[i][j]+B[i][j];

	return S;
}

double* matrizFila(double** A, int nFila, int nColumnas){
	double* f = calloc(nColumnas, sizeof(double));

	for(int i = 0;i<nColumnas;i++)
		f[i] = A[i][nFila];

	return f;
}

double* Cheb3D(double t, int N, double Ta, double Tb, double* Cx, double* Cy, double* Cz){
	// Comprobamos fronteras.
	if ( (t<Ta) || (Tb<t) )
	    printf("ERROR: Time out of range in Cheb3D::Value\n");

	// Algoritmo de Clenshaw
	double tau = (2*t-Ta-Tb)/(Tb-Ta);

	double* f1 = calloc(3,sizeof(double));
	double* f2 = calloc(3,sizeof(double));
	double* old_f1 = calloc(3,sizeof(double));

	for(int i = N; i>=2; i--){
		for(int j = 0; j < 3; j++){
			old_f1[j] = f1[j];
			f1[j] = 2*tau*f1[j]-f2[j];
			f2[j] = old_f1[j];
		}
		f1[0] += Cx[i-1];
		f1[1] += Cy[i-1];
		f1[2] += Cz[i-1];
	}

	double* ChebRet = calloc(3,sizeof(double));

	ChebRet[0]=tau*f1[0]-f2[0]+Cx[0];
	ChebRet[1]=tau*f1[1]-f2[1]+Cy[0];
	ChebRet[2]=tau*f1[2]-f2[2]+Cz[0];

	free(f1);
	free(f2);
	free(old_f1);

	return ChebRet;
}

void liberar(double** A, int m){
	for(int i = 0; i<m; i++)
		free(A[i]);
	free(A);
}
