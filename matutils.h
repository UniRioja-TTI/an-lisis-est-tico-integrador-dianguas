/*
 * matutils.h
 *
 *  Created on: 1 abr 2021
 *      Author: Administrador
 */

#ifndef MATUTILS_H_
#define MATUTILS_H_

double** zeros(int n, int m); // Devuelve una matriz nxm inicializada a 0.
int compare(double* A, double* B, int n); // Compara los vectores A y B, de tamaño n.
int compare_mat(double** A, double** B, int n, int m); // Compara las matrices A y B, de tamaño nxm.
double norm(double* v, int n); // Devuelve la norma de v. (n es el tamaño de v)
double dot(double* u, double* v, int n); // Devuelve el producto escalar u·v. (u y v son de tamaño n)
double* unit(double* v, int n); // Devuelve el vector v/|v|. Si v=0, devuelve el vector 0. (n es el tamaño de v)
double** producto(double** A, double** B, int n, int m, int p); // Devuelve el producto A*B, con A matriz nxm y B matriz mxp.
double** producto_escalar_mat(double** A, int n, double k); // Devuelve el producto k*A. (A es una matriz nxn)
//double** producto_escalar(double** A, int n, int m, double k); // Devuelve el producto k*A. (A es una matriz nxm)
double* producto_escalar(double* v, int n, double k); // Devuelve el producto k*v.
double** vect_a_matr(double* v); // Devuelve un puntero que apunta a v.
double** transpose(double** A, int n, int m); // Devuelve la matriz transpuesta de A. (A es de tamaño nxm)
double** suma(double** A, double** B, int n, int m); // Devuelve la matriz A+B. (A y B son de tamaño nxm)
double* matrizFila(double** A, int nFila, int nColumnas); // Devuelve el vector fila en nFila de la matriz A.
double* Cheb3D(double t, int N, double Ta, double Tb, double* Cx, double* Cy, double* Cz); // Aproximación de Chebyshev de vectores 3-dimensionales.
void liberar(double** A, int m); // Libera la memoria ocupada por A. (A tiene m columnas)

#endif /* MATUTILS_H_ */
