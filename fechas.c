/*
 * fechas.c
 *
 *  Created on: 1 abr 2021
 *      Author: Administrador
 */

#include "rotations.h"
#include "matutils.h"
#include "mathutils.h"
#include "fechas.h"
#include "logica.h"
#include "globales.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double gmst(double Mjd_UT1) {
	double MJD_J2000 = 51544.5;

	double Mjd_0 = floor(Mjd_UT1);
	double UT1 = 86400.0 * (Mjd_UT1 - Mjd_0);
	double T_0 = (Mjd_0 -MJD_J2000)/36525.0;
	double T = (Mjd_UT1 - MJD_J2000) / 36525.0;

	double gmstv = 24110.54841 + 8640184.812866 * T_0 + 1.002737909350795 * UT1+ (0.093104 - 6.2e-6 * T) * T * T;

	return 2*M_PI*Frac(gmstv/86400.0);
}

double gast(double Mjd_UT1){
	return fmod(gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*M_PI);
}

double Mjday(int yr, int mon, int day, int hr, int min, int sec){
	double jd = 367.0 * yr
	    - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )
	    + floor( 275 * mon / 9.0 )
	    + day + 1721013.5
	    + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;

	return jd-2400000.5;
}

//double Mjday(int yr, int mon, int day){
//	double jd = 367.0 * yr
//    - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 )
//    + floor( 275 * mon / 9.0 )
//    + day + 1721013.5;
//
//	return jd-2400000.5;
//}

double Mjday_TDB(double Mjd_TT){
	// Computa los siglos Julianos para TT.
	double T_TT = (Mjd_TT - 51544.5)/36525;

	// Computa la fecha Juliana modificada de TDB.
	return Mjd_TT + ( 0.001658*sin(628.3076*T_TT + 6.2401)
	              +   0.000022*sin(575.3385*T_TT+4.2970)
	              +   0.000014*sin(1256.6152*T_TT + 6.1969)
	              +   0.000005*sin(606.9777*T_TT+4.0212)
	              +   0.000005*sin(52.9691*T_TT+0.4444)
	              +   0.000002*sin(21.3299*T_TT+5.5431)
	              +   0.000010*sin(628.3076*T_TT+4.2490) )/86400;
}

void TimeUpdate(double*** P, double** Phi, int n){
	double** Phit = transpose(Phi, n,n);
	double** PhiP = producto(Phi,*P,n,n,n);
	*P = producto(PhiP,Phit,n,n,n);

	liberar(Phit, n);
	liberar(PhiP, n);
	liberar(*P, n);
}
