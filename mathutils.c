/*
 * mathutils.c
 *
 *  Created on: 1 abr 2021
 *      Author: Administrador
 */

#include "globales.h"
#include "rotations.h"
#include "matutils.h"
#include "mathutils.h"
#include "fechas.h"
#include "logica.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

void Legendre(int n, int m, double fi, double*** pnm, double*** dpnm){ // TODO descr
	*pnm = zeros(n+1, m+1);
	*dpnm = zeros(n+1, m+1);

	(*pnm)[0][0] = 1.0;
	(*dpnm)[0][0] = 0.0;

	(*pnm)[1][1]=sqrt(3)*cos(fi);
	(*dpnm)[1][1]=-sqrt(3)*sin(fi);

	// diagonal coefficients
	for(int i=2;i<n+1;i++)
	    (*pnm)[i][i] = sqrt((double)(2*i+1)/(2*i))*cos(fi)*(*pnm)[i-1][i-1];

	for(int i=2;i<n+1;i++)
	    (*dpnm)[i][i]= sqrt((double)(2*i+1)/(2*i))*((cos(fi)*(*dpnm)[i-1][i-1])-(sin(fi)*(*pnm)[i-1][i-1]));

	// horizontal first step coefficients
	for(int i=1;i<n+1;i++)
	    (*pnm)[i-1][i]= sqrt((double)(2*i+1))*sin(fi)*((*pnm)[i-1][i-1]);

	for(int i=1;i<n+1;i++)
	    (*dpnm)[i-1][i]= sqrt((double)(2*i+1))*((cos(fi)*((*pnm)[i-1][i-1]))+(sin(fi)*((*dpnm)[i-1][i-1])));

	// horizontal second step coefficients
	for(int j = 0, k = 2;1==1;){
	    for(int i = k;i<n+1;i++)
	        (*pnm)[j][i]=sqrt((double)(2*i+1)/((i-j)*(i+j)))*((sqrt((double)(2*i-1))*sin(fi)*(*pnm)[j][i-1])-(sqrt((double)((i+j-1)*(i-j-1))/(2*i-3))* (*pnm)[j][i-2]));
	    j++;
	    k++;
	    if(j > m)
	    	break;
	}

	for(int j = 0, k = 2;1==1;){
	    for(int i = k;i<n+1;i++)
	        (*dpnm)[j][i]=sqrt((double)(2*i+1)/((i-j)*(i+j)))*((sqrt((double)(2*i-1))*sin(fi)*(*dpnm)[j][i-1])+(sqrt((double)(2*i-1))*cos(fi)* (*pnm)[j][i-1])-(sqrt((double)((i+j-1)*(i-j-1))/(2*i-3))* (*dpnm)[j][i-2]));
	    j++;
	    k++;
	    if(j > m)
	    	break;
	}
}

double Frac(double x){
	return x-floor(x);
}

double sign(double a){
	double result;
	if (a>=0.0)
	    result = 1.0;
	else
	    result = -1.0;
	return result;
}

//double sign(double a, double b){
//	double result;
//	if (b>=0.0)
//	    result = abs(a);
//	else
//	    result = - abs(a);
//	return result;
//}

double angl(double* u, double* v, int n){
	double undefined = 999999.1;

	double Nu = norm(u,n);
	double Nv = norm(v,n);

	double theta;

	if (Nu*Nv > DBL_EPSILON){
	    double temp = dot(u,v,n) / (Nu*Nv);
	    if (abs(temp) > 1.0)
	        temp = sign(temp) * 1.0;
	    theta = acos(temp);
	}
	else
	    theta = undefined;

	return theta;
}

double* VarEqn(double x, double* yPhi){

	//extern double** eopdata;
	int nAux =      20;
	int mAux =      20;
	double MJD_J2000 = 51544.5;

	double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

	IERS(eopdata, Mjd_UTCAux + x / 86400.0, 'l', &x_pole, &y_pole, &UT1_UTC, &LOD,
			&dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

	double Mjd_UT1 = Mjd_TTAux + (UT1_UTC - TT_UTC) / 86400.0;

	double **P = PrecMatrix(MJD_J2000, Mjd_TTAux + x / 86400.0);
	double **N = NutMatrix(Mjd_TTAux + x / 86400.0);

	double **T = producto(N, P, 3,3,3);
	double **E = producto(producto(PoleMatrix(x_pole, y_pole), GHAMatrix(Mjd_UT1), 3,3,3), T, 3,3,3);

	double *r = calloc(3, sizeof(double));
	double *v = calloc(3, sizeof(double));

	for (int j = 0; j < 3; j++) {
		r[j] = yPhi[j];
		v[j] = yPhi[j + 3];
	}

	double **Phi = zeros(6,6);

	for (int j = 1; j <= 6; j++)
		for (int jj = 0; jj < 6; jj++)
			Phi[j-1][jj] = yPhi[6 * j + jj];

	double *a = AccelHarmonic(r, E, nAux, mAux);
	double **G = G_AccelHarmonic(r, E, nAux, mAux);

	double *yPhip = calloc(42, sizeof(double));
	double **dfdy = zeros(6,6);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			dfdy[j][i] = 0.0;
			dfdy[j][i + 3] = G[j][i];
			if (i == j)
				dfdy[j + 3][i] = 1.0;
			else
				dfdy[j + 3][i] = 0.0;
			dfdy[j + 3][i + 3] = 0.0;
		}
	}

	double **Phip = producto(dfdy, Phi, 6,6,6);

	for (int j = 0; j < 3; j++) {
		yPhip[j] = v[j];
		yPhip[j+3] = a[j];
	}

	for (int j = 1; j <= 6; j++)
			for (int jj = 1; jj <= 6; jj++)
				yPhip[6*j+jj-1] = Phip[j-1][jj-1];

	liberar(P, 3);
	liberar(N, 3);
	liberar(T, 3);
	liberar(E, 3);
	liberar(G, 3);
	liberar(dfdy, 6);
	liberar(Phip, 6);
	liberar(Phi, 6);
	free(r);
	free(v);
	free(a);

	return yPhip;
}

void VarEqnWrapper(double x, double* yPhi, double** yPhip){
	*yPhip = VarEqn(x, yPhi);
}
