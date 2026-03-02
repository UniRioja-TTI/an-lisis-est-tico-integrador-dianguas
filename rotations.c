/*
 * rotations.c
 *
 *  Created on: 1 abr 2021
 *      Author: Administrador
 */

#include <math.h>
#include "matutils.h"
#include "rotations.h"

double** R_x(double ang){ // Devuelve una matriz rotacional 3x3 para girar un vector en el eje x ang radianes
	double** matrizRotacional = zeros(3,3);

	double S = sin(ang);
	double C = cos(ang);

	//rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
	//rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
	//rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;
	// En C transpondremos las matrices para poder acceder fácilmente a los vectores columna

	matrizRotacional[0][0] = 1.0;
	matrizRotacional[0][1] = 0.0;
	matrizRotacional[0][2] = 0.0;
	matrizRotacional[1][0] = 0.0;
	matrizRotacional[1][1] = C;
	matrizRotacional[1][2] = -1.0*S;
	matrizRotacional[2][0] = 0.0;
	matrizRotacional[2][1] = S;
	matrizRotacional[2][2] = C;

	return matrizRotacional;
}

double** R_y(double ang){ // Devuelve una matriz rotacional 3x3 para girar un vector en el eje y ang radianes
	double** matrizRotacional = zeros(3,3);

	double S = sin(ang);
	double C = cos(ang);

	//rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
	//rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
	//rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;
	// En C transpondremos las matrices para poder acceder fácilmente a los vectores columna

	matrizRotacional[0][0] = C;
	matrizRotacional[0][1] = 0.0;
	matrizRotacional[0][2] = S;
	matrizRotacional[1][0] = 0.0;
	matrizRotacional[1][1] = 1.0;
	matrizRotacional[1][2] = 0.0;
	matrizRotacional[2][0] = -1.0*S;
	matrizRotacional[2][1] = 0.0;
	matrizRotacional[2][2] = C;

	return matrizRotacional;
}

double** R_z(double ang){ // Devuelve una matriz rotacional 3x3 para girar un vector en el eje z ang radianes
	double** matrizRotacional = zeros(3,3);

	double S = sin(ang);
	double C = cos(ang);

	//	rotmat(1,1) =      C;  rotmat(1,2) =   S;  rotmat(1,3) = 0.0;
	//	rotmat(2,1) = -1.0*S;  rotmat(2,2) =   C;  rotmat(2,3) = 0.0;
	//	rotmat(3,1) =    0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;
	// En C transpondremos las matrices para poder acceder fácilmente a los vectores columna

	matrizRotacional[0][0] = C;
	matrizRotacional[0][1] = -1.0*S;
	matrizRotacional[0][2] = 0.0;
	matrizRotacional[1][0] = S;
	matrizRotacional[1][1] = C;
	matrizRotacional[1][2] = 0.0;
	matrizRotacional[2][0] = 0.0;
	matrizRotacional[2][1] = 0.0;
	matrizRotacional[2][2] = 1.0;

	return matrizRotacional;
}
