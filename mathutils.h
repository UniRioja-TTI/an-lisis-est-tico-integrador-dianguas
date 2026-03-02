/*
 * mathutils.h
 *
 *  Created on: 1 abr 2021
 *      Author: Administrador
 */

#ifndef MATHUTILS_H_
#define MATHUTILS_H_

void Legendre(int n, int m, double fi, /* Salidas */ double*** pnm, double*** dpnm); // TODO descr
double Frac(double x); // Devuelve la parte fraccionaria de x.
double sign(double a); // Devuelve 1.0 si a es positivo y -1.0 si a es negativo.
//double sign(double a, double b); // Devuelve el valor absoluto de a con el signo de b.
double angl(double* u, double* v, int n); // Devuelve el ángulo entre los vectores u y v. (n es la dimensión de u y v)
double* VarEqn(double x, double* yPhi); // Computa las ecuaciones variacionales.
void VarEqnWrapper(double x, double* yPhi, double** yPhip);

#endif /* MATHUTILS_H_ */
