/*
 * logica.c
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
#include <float.h>
#include <math.h>


// eps = MeanObliquity (Mjd_TT);


//extern double** producto(double** A, double** B, int n, int m, int p);
//extern double** suma(double** A, double** B, int n, int m);

double** LTC(double lon, double lat){
	double** M = producto(R_y(-1.0*lat), R_z(lon), 3,3,3);

	for(int j = 0; j<3; j++){
	    double AuxFuncLTC = M[j][0];
	    M[j][0]=M[j][1];
	    M[j][1]=M[j][2];
	    M[j][2]= AuxFuncLTC;
	}

	return M;
}

double** PoleMatrix(double xp, double yp){
	return producto(R_y(-xp), R_x(-yp), 3,3,3);
}

double** PrecMatrix(double Mjd_1, double Mjd_2){
	double MJD_J2000 = 51544.5;
	double Arcs      = 3600*180/M_PI;

	double T  = (Mjd_1-MJD_J2000)/36525;
	double dT = (Mjd_2-Mjd_1)/36525;

	// Ángulos de precesión
	double zeta  = ((2306.2181+(1.39656-0.000139*T)*T)+((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Arcs;
	double z     = zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Arcs;
	double theta = ((2004.3109-(0.85330+0.000217*T)*T)-((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Arcs;

	return producto(producto(R_z(-z), R_y(theta), 3,3,3), R_z(-zeta), 3,3,3);
}

double* Position(double lon, double lat, double h){
	double R_Earth   = 6378.1363e3;
	double f_Earth   = 1/298.257223563;

	double* r = (double*)calloc(3, sizeof(double));

	double e2     = f_Earth*(2.0-f_Earth);   // Cuadrado de la eccentricidad
	double CosLat = cos(lat);    // (Co)seno de la latitud geodética
	double SinLat = sin(lat);

	// Vector posición
	double N = R_Earth / sqrt(1.0-e2*SinLat*SinLat);

	r[0] =  (         N+h)*CosLat*cos(lon);
	r[1] =  (         N+h)*CosLat*sin(lon);
	r[2] =  ((1.0-e2)*N+h)*SinLat;

	return r;
}

void Geodetic(double* r,/* Salida */ double* lon, double* lat, double* h){
	double R_Earth   = 6378.1363e3;
	double f_Earth   = 1/298.257223563;

	double epsRequ = DBL_EPSILON*R_Earth;        // Criterio de convergencia
	double e2      = f_Earth*(2.0-f_Earth);        // Cuadrado de la eccentricidad

	double X = r[0]; // Coordenadas cartesianas
	double Y = r[1];
	double Z = r[2];
	double rho2 = X*X + Y*Y;           // Cuadrado de la distancia desde el eje z

	// Comprobamos validez de los datos de entrada.
	if (norm(r,3)==0.0){
	    printf( " entrada invalida en la funcion Geodetic\n" );
	    lon[0] = 0.0;
	    lat[0] = 0.0;
	    h[0]   = -R_Earth;
	}

	// Iteración
	double dZ = e2*Z;
	double ZdZ;
	double Nh;
	double N;
	while(1){
		ZdZ    =  Z + dZ;
		Nh     =  sqrt ( rho2 + ZdZ*ZdZ );
		double SinPhi =  ZdZ / Nh; // Seno de la latitud geodética.
		N      =  R_Earth / sqrt(1.0-e2*SinPhi*SinPhi);
		double dZ_new =  N*e2*SinPhi;
	    if (abs(dZ-dZ_new) < epsRequ)
	        break;
	    dZ = dZ_new;
	}

	// Longitud, latitud, altitud.
	lon[0] = atan2 ( Y, X );
	lat[0] = atan2 ( ZdZ, sqrt(rho2) );
	h[0]   = Nh - N;
}

void MeasUpdate(double** x, double z, double g, double s, double* G, double*** P, int n){
	double W = s*s;

	double** GMat = vect_a_matr(G);
	double** PT = transpose(*P, 6, 6);

	double** GP = producto(PT, GMat, 6,6,1);
	liberar(PT, 6);

	double multi = dot(GP[0],GMat[0],6)+W;
	liberar(GP,1);
	double** Kprev = producto(*P, GMat, 6, 6, 1);
	double* K = producto_escalar(Kprev[0], 6, 1.0/multi);
	liberar(GMat,1);
	liberar(Kprev,1);

	for(int j = 0; j<6; j++)
		(*x)[j] += K[j]*(z-g);

	double** GT = transpose(vect_a_matr(G), 6, 1);
	double** KG = producto(vect_a_matr(K),GT,6,1,6);
	liberar(GT, 6);
	double** eye = zeros(n,n);
	for(int jj = 0;jj < n; jj++)
		eye[jj][jj] = 1.0;

	double** mKG = producto_escalar_mat(KG, 6, -1.0);
	liberar(KG, 6);
	double** eyemKG = suma(eye, mKG, 6, 6);
	liberar(eye,6);
	liberar(mKG,6);

	*P = producto(eyemKG, *P, 6,6,6);
	free(K);
}

void NutAngles(double Mjd_TT, double* dpsi, double* deps){
	double MJD_J2000 = 51544.5;
	double Arcs      = 3600*180/M_PI;

	double T  = (Mjd_TT-MJD_J2000)/36525;
	double T2 = T*T;
	double T3 = T2*T;
	int rev = 360*3600;

	int N_coeff = 106;
	int C[954] = {
	  // l  l' F  D Om    dpsi    *T     deps     *T
	    0, 0, 0, 0, 1,-1719960,-1742,  920250,   89 ,   //   1
	    0, 0, 0, 0, 2,   20620,    2,   -8950,    5 ,   //   2
	   -2, 0, 2, 0, 1,     460,    0,    -240,    0 ,   //   3
	    2, 0,-2, 0, 0,     110,    0,       0,    0 ,   //   4
	   -2, 0, 2, 0, 2,     -30,    0,      10,    0 ,   //   5
	    1,-1, 0,-1, 0,     -30,    0,       0,    0 ,   //   6
	    0,-2, 2,-2, 1,     -20,    0,      10,    0 ,   //   7
	    2, 0,-2, 0, 1,      10,    0,       0,    0 ,   //   8
	    0, 0, 2,-2, 2, -131870,  -16,   57360,  -31 ,   //   9
	    0, 1, 0, 0, 0,   14260,  -34,     540,   -1 ,   //  10
	    0, 1, 2,-2, 2,   -5170,   12,    2240,   -6 ,   //  11
	    0,-1, 2,-2, 2,    2170,   -5,    -950,    3 ,   //  12
	    0, 0, 2,-2, 1,    1290,    1,    -700,    0 ,   //  13
	    2, 0, 0,-2, 0,     480,    0,      10,    0 ,   //  14
	    0, 0, 2,-2, 0,    -220,    0,       0,    0 ,   //  15
	    0, 2, 0, 0, 0,     170,   -1,       0,    0 ,   //  16
	    0, 1, 0, 0, 1,    -150,    0,      90,    0 ,   //  17
	    0, 2, 2,-2, 2,    -160,    1,      70,    0 ,   //  18
	    0,-1, 0, 0, 1,    -120,    0,      60,    0 ,   //  19
	   -2, 0, 0, 2, 1,     -60,    0,      30,    0 ,   //  20
	    0,-1, 2,-2, 1,     -50,    0,      30,    0 ,   //  21
	    2, 0, 0,-2, 1,      40,    0,     -20,    0 ,   //  22
	    0, 1, 2,-2, 1,      40,    0,     -20,    0 ,   //  23
	    1, 0, 0,-1, 0,     -40,    0,       0,    0 ,   //  24
	    2, 1, 0,-2, 0,      10,    0,       0,    0 ,   //  25
	    0, 0,-2, 2, 1,      10,    0,       0,    0 ,   //  26
	    0, 1,-2, 2, 0,     -10,    0,       0,    0 ,   //  27
	    0, 1, 0, 0, 2,      10,    0,       0,    0 ,   //  28
	   -1, 0, 0, 1, 1,      10,    0,       0,    0 ,   //  29
	    0, 1, 2,-2, 0,     -10,    0,       0,    0 ,   //  30
	    0, 0, 2, 0, 2,  -22740,   -2,    9770,   -5 ,   //  31
	    1, 0, 0, 0, 0,    7120,    1,     -70,    0 ,   //  32
	    0, 0, 2, 0, 1,   -3860,   -4,    2000,    0 ,   //  33
	    1, 0, 2, 0, 2,   -3010,    0,    1290,   -1 ,   //  34
	    1, 0, 0,-2, 0,   -1580,    0,     -10,    0 ,   //  35
	   -1, 0, 2, 0, 2,    1230,    0,    -530,    0 ,   //  36
	    0, 0, 0, 2, 0,     630,    0,     -20,    0 ,   //  37
	    1, 0, 0, 0, 1,     630,    1,    -330,    0 ,   //  38
	   -1, 0, 0, 0, 1,    -580,   -1,     320,    0 ,   //  39
	   -1, 0, 2, 2, 2,    -590,    0,     260,    0 ,   //  40
	    1, 0, 2, 0, 1,    -510,    0,     270,    0 ,   //  41
	    0, 0, 2, 2, 2,    -380,    0,     160,    0 ,   //  42
	    2, 0, 0, 0, 0,     290,    0,     -10,    0 ,   //  43
	    1, 0, 2,-2, 2,     290,    0,    -120,    0 ,   //  44
	    2, 0, 2, 0, 2,    -310,    0,     130,    0 ,   //  45
	    0, 0, 2, 0, 0,     260,    0,     -10,    0 ,   //  46
	   -1, 0, 2, 0, 1,     210,    0,    -100,    0 ,   //  47
	   -1, 0, 0, 2, 1,     160,    0,     -80,    0 ,   //  48
	    1, 0, 0,-2, 1,    -130,    0,      70,    0 ,   //  49
	   -1, 0, 2, 2, 1,    -100,    0,      50,    0 ,   //  50
	    1, 1, 0,-2, 0,     -70,    0,       0,    0 ,   //  51
	    0, 1, 2, 0, 2,      70,    0,     -30,    0 ,   //  52
	    0,-1, 2, 0, 2,     -70,    0,      30,    0 ,   //  53
	    1, 0, 2, 2, 2,     -80,    0,      30,    0 ,   //  54
	    1, 0, 0, 2, 0,      60,    0,       0,    0 ,   //  55
	    2, 0, 2,-2, 2,      60,    0,     -30,    0 ,   //  56
	    0, 0, 0, 2, 1,     -60,    0,      30,    0 ,   //  57
	    0, 0, 2, 2, 1,     -70,    0,      30,    0 ,   //  58
	    1, 0, 2,-2, 1,      60,    0,     -30,    0 ,   //  59
	    0, 0, 0,-2, 1,     -50,    0,      30,    0 ,   //  60
	    1,-1, 0, 0, 0,      50,    0,       0,    0 ,   //  61
	    2, 0, 2, 0, 1,     -50,    0,      30,    0 ,   //  62
	    0, 1, 0,-2, 0,     -40,    0,       0,    0 ,   //  63
	    1, 0,-2, 0, 0,      40,    0,       0,    0 ,   //  64
	    0, 0, 0, 1, 0,     -40,    0,       0,    0 ,   //  65
	    1, 1, 0, 0, 0,     -30,    0,       0,    0 ,   //  66
	    1, 0, 2, 0, 0,      30,    0,       0,    0 ,   //  67
	    1,-1, 2, 0, 2,     -30,    0,      10,    0 ,   //  68
	   -1,-1, 2, 2, 2,     -30,    0,      10,    0 ,   //  69
	   -2, 0, 0, 0, 1,     -20,    0,      10,    0 ,   //  70
	    3, 0, 2, 0, 2,     -30,    0,      10,    0 ,   //  71
	    0,-1, 2, 2, 2,     -30,    0,      10,    0 ,   //  72
	    1, 1, 2, 0, 2,      20,    0,     -10,    0 ,   //  73
	   -1, 0, 2,-2, 1,     -20,    0,      10,    0 ,   //  74
	    2, 0, 0, 0, 1,      20,    0,     -10,    0 ,   //  75
	    1, 0, 0, 0, 2,     -20,    0,      10,    0 ,   //  76
	    3, 0, 0, 0, 0,      20,    0,       0,    0 ,   //  77
	    0, 0, 2, 1, 2,      20,    0,     -10,    0 ,   //  78
	   -1, 0, 0, 0, 2,      10,    0,     -10,    0 ,   //  79
	    1, 0, 0,-4, 0,     -10,    0,       0,    0 ,   //  80
	   -2, 0, 2, 2, 2,      10,    0,     -10,    0 ,   //  81
	   -1, 0, 2, 4, 2,     -20,    0,      10,    0 ,   //  82
	    2, 0, 0,-4, 0,     -10,    0,       0,    0 ,   //  83
	    1, 1, 2,-2, 2,      10,    0,     -10,    0 ,   //  84
	    1, 0, 2, 2, 1,     -10,    0,      10,    0 ,   //  85
	   -2, 0, 2, 4, 2,     -10,    0,      10,    0 ,   //  86
	   -1, 0, 4, 0, 2,      10,    0,       0,    0 ,   //  87
	    1,-1, 0,-2, 0,      10,    0,       0,    0 ,   //  88
	    2, 0, 2,-2, 1,      10,    0,     -10,    0 ,   //  89
	    2, 0, 2, 2, 2,     -10,    0,       0,    0 ,   //  90
	    1, 0, 0, 2, 1,     -10,    0,       0,    0 ,   //  91
	    0, 0, 4,-2, 2,      10,    0,       0,    0 ,   //  92
	    3, 0, 2,-2, 2,      10,    0,       0,    0 ,   //  93
	    1, 0, 2,-2, 0,     -10,    0,       0,    0 ,   //  94
	    0, 1, 2, 0, 1,      10,    0,       0,    0 ,   //  95
	   -1,-1, 0, 2, 1,      10,    0,       0,    0 ,   //  96
	    0, 0,-2, 0, 1,     -10,    0,       0,    0 ,   //  97
	    0, 0, 2,-1, 2,     -10,    0,       0,    0 ,   //  98
	    0, 1, 0, 2, 0,     -10,    0,       0,    0 ,   //  99
	    1, 0,-2,-2, 0,     -10,    0,       0,    0 ,   // 100
	    0,-1, 2, 0, 1,     -10,    0,       0,    0 ,   // 101
	    1, 1, 0,-2, 1,     -10,    0,       0,    0 ,   // 102
	    1, 0,-2, 2, 0,     -10,    0,       0,    0 ,   // 103
	    2, 0, 0, 2, 0,      10,    0,       0,    0 ,   // 104
	    0, 0, 2, 4, 2,     -10,    0,       0,    0 ,   // 105
	    0, 1, 0, 1, 0,      10,    0,       0,    0     // 106
	};

	double l  = fmod (  485866.733 + (1325.0*rev +  715922.633)*T + 31.310*T2 + 0.064*T3, rev );
	double lp = fmod ( 1287099.804 + (  99.0*rev + 1292581.224)*T -  0.577*T2 - 0.012*T3, rev );
	double F  = fmod (  335778.877 + (1342.0*rev +  295263.137)*T - 13.257*T2 + 0.011*T3, rev );
	double D  = fmod ( 1072261.307 + (1236.0*rev + 1105601.328)*T -  6.891*T2 + 0.019*T3, rev );
	double Om = fmod (  450160.280 - (   5.0*rev +  482890.539)*T +  7.455*T2 + 0.008*T3, rev );

	// Nutación en longitud y oblicuidad [rad]

	dpsi[0] = 0.0;
	deps[0] = 0.0;

	for(int i=0;i<N_coeff;i++){
		  double arg  =  ( C[9*i+0]*l+C[9*i+1]*lp+C[9*i+2]*F+C[9*i+3]*D+C[9*i+4]*Om )/Arcs;
		  dpsi[0] += ( C[9*i+5]+C[9*i+6]*T ) * sin(arg);
		  deps[0] += ( C[9*i+7]+C[9*i+8]*T ) * cos(arg);
	}

	// free(C);

	dpsi[0] = 1.0e-5 * dpsi[0]/Arcs;
	deps[0] = 1.0e-5 * deps[0]/Arcs;
}

double MeanObliquity(double Mjd_TT){
	double MJD_J2000 = 51544.5;
	double Rad       = M_PI/180;

	double T = (Mjd_TT-MJD_J2000)/36525;
	return Rad*(84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600);
}

double** NutMatrix(double Mjd_TT){
	// Oblicuidad media de la eclíptica
	double eps = MeanObliquity(Mjd_TT);

	// Nutación eanObliqun longitud y oblicuidad
	double dpsi = 0.0;
	double deps = 0.0;
	NutAngles(Mjd_TT,&dpsi,&deps);

	// Transformación de ecuador y equinocio medios a reales
	double** res = producto(producto(R_x(-eps-deps), R_z(-dpsi), 3,3,3), R_x(+eps), 3,3,3);

	return res;
}

double EqnEquinox(double Mjd_TT){
	// Nutación en longitud y oblicuidad
	double dpsi = 0.0;
	double deps = 0.0;
	NutAngles(Mjd_TT,&dpsi,&deps);

	// Ecuación de los equinocios
	double res = dpsi*cos(MeanObliquity(Mjd_TT));

	return res;
}

void AzElPa(double* s, /* Salidas */ double* Az, double* El, double** dAds, double** dEds){
	double pi2 = 2.0*M_PI;

	double rho = sqrt(s[0]*s[0]+s[1]*s[1]);

	*Az = atan2(s[0],s[1]);

	if(*Az < 0.0)
		*Az = *Az + pi2;

	*El = atan(s[2]/rho);

	*dAds = calloc(3,sizeof(double));
	*dEds = calloc(3,sizeof(double));

	(*dAds)[0] = s[1]/(rho*rho);
	(*dAds)[1] = -s[0]/(rho*rho);
	(*dAds)[2] = 0.0;

	double dotv = dot(s,s,3);

	(*dEds)[0] = -s[0]*s[2]/(rho*dotv);
	(*dEds)[1] = -s[1]*s[2]/(rho*dotv);
	(*dEds)[2] = rho/dotv;
}

void IERS(double** eop, double Mjd_UTCParam, char interp, double *x_pole, double *y_pole,
		double *UT1_UTC, double *LOD, double *dpsi, double *deps,
		double *dx_pole, double *dy_pole, double *TAI_UTC){

	double Arcs      = 3600*180/M_PI;
	// Código anulado.
	/*
	 * if (nargin == 2)
   	 * interp = 'n';
	 * end
	 */

	if (interp =='l'){

		// linear interpolation
		double mjd = floor(Mjd_UTCParam);
		double* fila3 = matrizFila(eop,3,21413);
		int i = 0;

		for(;i<21413;i++){
			if(fila3[i]==mjd)
				break;
		}

		free(fila3);

		double *preeop = eop[i]; // No hacer free.
		double *nexteop = eop[i+1];

		double mfme = 1440.0*(Mjd_UTCParam-floor(Mjd_UTCParam));
		double fixf = mfme/1440.0;

		// Configuramos los parámetros IERS de rotación terrestre
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		*x_pole  = preeop[4]+(nexteop[4]-preeop[4])*fixf;
		*y_pole  = preeop[5]+(nexteop[5]-preeop[5])*fixf;
		*UT1_UTC = preeop[6]+(nexteop[6]-preeop[6])*fixf;
		*LOD     = preeop[7]+(nexteop[7]-preeop[7])*fixf;
		*dpsi    = preeop[8]+(nexteop[8]-preeop[8])*fixf;
		*deps    = preeop[9]+(nexteop[9]-preeop[9])*fixf;
		*dx_pole = preeop[10]+(nexteop[10]-preeop[10])*fixf;
		*dy_pole = preeop[11]+(nexteop[11]-preeop[11])*fixf;
		*TAI_UTC = preeop[12];

		*x_pole  = (*x_pole)/Arcs;  // Coordenada de polo [rad]
		*y_pole  = (*y_pole)/Arcs;  // Coordenada de polo [rad]
		*dpsi    = (*dpsi)/Arcs;
		*deps    = (*deps)/Arcs;
		*dx_pole = (*dx_pole)/Arcs; // Coordenada de polo [rad]
		*dy_pole = (*dy_pole)/Arcs; // Coordenada de polo [rad]
	}

	else if (interp =='n'){

		double mjd = (floor(Mjd_UTCParam));

		double* fila3 = matrizFila(eop,3,21413);
		int i = 0;

		for(;i<21413;i++){
			if(fila3[i]==mjd)
				break;
		}

		free(fila3);

	    double* eopV = eop[i]; // No hacer free.

	    // Configuramos los parámetros IERS de rotación terrestre
	    // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
	    *x_pole  = eopV[4]/Arcs;  // Coordenada de polo [rad]
	    *y_pole  = eopV[5]/Arcs;  // Coordenada de polo [rad]
		*UT1_UTC = eopV[6];             // Diferencia temporal UT1-UTC [s]
	    *LOD     = eopV[7];             // Longitud del día [s]
	    *dpsi    = eopV[8]/Arcs;
	    *deps    = eopV[9]/Arcs;
	    *dx_pole = eopV[10]/Arcs; // Coordenada de polo [rad]
	    *dy_pole = eopV[11]/Arcs; // Coordenada de polo [rad]
		*TAI_UTC = eopV[12];            // Diferencia temporal UT1-UTC [s]
	}
}

void timediff(double UT1_UTC, double TAI_UTC, double* UT1_TAI, double* UTC_GPS, double* UT1_GPS, double* TT_UTC, double* GPS_UTC){
	double TT_TAI  = +32.184;          // Diferencia temporal TT-TAI  [s]
	double GPS_TAI = -19.0;            // Diferencia temporal GPS-TAI [s]
//	double TT_GPS  =  TT_TAI-GPS_TAI;  // Diferencia temporal TT-GPS  [s]
//	double TAI_GPS = -GPS_TAI;         // Diferencia temporal TAI-GPS [s]
	double UTC_TAI = -TAI_UTC;         // Diferencia temporal UTC-TAI [s]

	*UT1_TAI = UT1_UTC-TAI_UTC;  // Diferencia temporal UT1-TAI [s]
	*UTC_GPS = UTC_TAI-GPS_TAI;  // Diferencia temporal UTC_GPS [s]
	*UT1_GPS = *UT1_TAI-GPS_TAI;  // Diferencia temporal UT1-GPS [s]
	*TT_UTC  = TT_TAI-UTC_TAI;   // Diferencia temporal TT-UTC  [s]
	*GPS_UTC = GPS_TAI-UTC_TAI;  // Diferencia temporal GPS-UTC [s]
}

void sacarCX(int inicio, int salto, int fin, int itBucle, double* PCtemp, /* Salidas */ double** Cx_Ret, double** Cy_Ret, double** Cz_Ret){
	int n = (fin-inicio)/salto+1;

	int temp[n];

	for (int j = 0; j < n; j++)
		temp[j] = inicio + salto * j - 1; // Reducimos en 1 el valor para pasar a índices de C.

	Cx_Ret[0] = calloc(salto*(itBucle+1), sizeof(double));
	if(Cx_Ret[0] == NULL)
		printf("Fallo de memoria\n");
	Cy_Ret[0] = calloc(salto*(itBucle+1), sizeof(double));
	if(Cy_Ret[0] == NULL)
		printf("Fallo de memoria\n");
	Cz_Ret[0] = calloc(salto*(itBucle+1), sizeof(double));
	if(Cz_Ret[0] == NULL)
		printf("Fallo de memoria\n");

	double* Cx = calloc(salto,sizeof(double));
	if(Cx == NULL)
		printf("Fallo de memoria\n");
	double* Cy = calloc(salto,sizeof(double));
	if(Cy == NULL)
		printf("Fallo de memoria\n");
	double* Cz = calloc(salto,sizeof(double));
	if(Cz == NULL)
		printf("Fallo de memoria\n");

	for (int j = 0; j < salto; j++) {
		(*Cx_Ret)[j] = PCtemp[temp[0] + j];
		(*Cy_Ret)[j] = PCtemp[temp[1] + j];
		(*Cz_Ret)[j] = PCtemp[temp[2] + j];
	}
	for(int i = 1; i <= itBucle; i++){
		for (int j = 0; j < n; j++)
			temp[j] += salto*3;

		for (int j = 0; j < salto; j++) {
			Cx[j] = PCtemp[temp[0] + j];
			Cy[j] = PCtemp[temp[1] + j];
			Cz[j] = PCtemp[temp[2] + j];
		}

		for (int j = 0; j < salto; j++) {
			(*Cx_Ret)[j+salto*i] = Cx[j];
			(*Cy_Ret)[j+salto*i] = Cy[j];
			(*Cz_Ret)[j+salto*i] = Cz[j];
		}
	}
	//free(temp); Crashea a veces sin ningún motivo. Lo pasamos a array y arreglado.

	free(Cx);
	free(Cy);
	free(Cz);
}

void JPL_Eph_DE430(double Mjd_TDB, /* Salidas */ double** r_Mercury, double** r_Venus, double** r_Earth, double** r_Mars, double** r_Jupiter, double** r_Saturn, double** r_Uranus, double** r_Neptune, double** r_Pluto, double** r_Moon, double** r_Sun){
	//extern double** PC;

	double JD = Mjd_TDB + 2400000.5;
	int i = 0;

	double* col0 = calloc(2285,sizeof(double));
	double* col1 = calloc(2285,sizeof(double));

	for(int j = 0;j < 2285; j++){
		col0[j] = PC[0][j];
		col1[j] = PC[1][j];
	}

	for(;i < 2285; i++){
		if(col0[i]<=JD && JD<=col1[i])
			break;
	}

	free(col0);
	free(col1);

	double* PCtemp = calloc(1020,sizeof(double));

	for(int j = 0; j < 1020; j++)
		PCtemp[j] = PC[j][i];

	double t1 = PCtemp[0]-2400000.5; // MJD al inicio del intervalo

	double dt = Mjd_TDB - t1;

	int j = 0;
	double Mjd0 = 0.0;

	double* Cx_Earth;
	double* Cy_Earth;
	double* Cz_Earth;

	sacarCX(231, 13, 270, 1, PCtemp, &Cx_Earth, &Cy_Earth, &Cz_Earth);

	if (0<=dt && dt<=16){
	    j=0;
	    Mjd0 = t1;
	}
	else if(16<dt && dt<=32){
	    j=1;
	    Mjd0 = t1+16*j;
	}

	double* r_EarthCheb = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, &Cx_Earth[13*j], &Cy_Earth[13*j], &Cz_Earth[13*j]);

	free(Cx_Earth);
	free(Cy_Earth);
	free(Cz_Earth);

	*r_Earth = producto_escalar(r_EarthCheb, 3, 1e3);
	free(r_EarthCheb);

	/* LUNA */

	double *Cx_Moon = calloc(104, sizeof(double));
	double *Cy_Moon = calloc(104, sizeof(double));
	double *Cz_Moon = calloc(104, sizeof(double));

	sacarCX(441, 13, 480, 7, PCtemp, &Cx_Moon, &Cy_Moon, &Cz_Moon);

	if (0<=dt && dt<=4){
	    j=0;
	    Mjd0 = t1;
	}
	else if(4<dt && dt<=8){
	    j=1;
	    Mjd0 = t1+4*j;
	}
	else if(8<dt && dt<=12){
	    j=2;
	    Mjd0 = t1+4*j;
	}
	else if(12<dt && dt<=16){
	    j=3;
	    Mjd0 = t1+4*j;
	}
	else if(16<dt && dt<=20){
	    j=4;
	    Mjd0 = t1+4*j;
	}
	else if(20<dt && dt<=24){
	    j=5;
	    Mjd0 = t1+4*j;
	}
	else if(24<dt && dt<=28){
	    j=6;
	    Mjd0 = t1+4*j;
	}
	else if(28<dt && dt<=32){
	    j=7;
	    Mjd0 = t1+4*j;
	}

	double* r_MoonCheb = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, &Cx_Moon[13*j], &Cy_Moon[13*j], &Cz_Moon[13*j]);

	free(Cx_Moon);
	free(Cy_Moon);
	free(Cz_Moon);

	*r_Moon = producto_escalar(r_MoonCheb, 3, 1e3);
	free(r_MoonCheb);

	/* SOL */

	double *Cx_Sun;
	double *Cy_Sun;
	double *Cz_Sun;

	sacarCX(753, 11, 786, 1, PCtemp, &Cx_Sun, &Cy_Sun, &Cz_Sun);

	if (0 <= dt && dt <= 16) {
		j = 0;
		Mjd0 = t1;
	} else if (16 < dt && dt <= 32) {
		j = 1;
		Mjd0 = t1 + 16 * j;
	}

	double *r_SunCheb = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16, &Cx_Sun[11 * j], &Cy_Sun[11 * j], &Cz_Sun[11 * j]);

	free(Cx_Sun);
	free(Cy_Sun);
	free(Cz_Sun);

	*r_Sun = producto_escalar(r_SunCheb, 3, 1e3);
	free(r_SunCheb);

	/* MERCURIO */

	double *Cx_Mercury;
	double *Cy_Mercury;
	double *Cz_Mercury;

	sacarCX(3, 14, 45, 3, PCtemp, &Cx_Mercury, &Cy_Mercury, &Cz_Mercury);
	if (0<=dt && dt<=8){
	    j=0;
	    Mjd0 = t1;
	}
	else if(8<dt && dt<=16){
	    j=1;
	    Mjd0 = t1+8*j;
	}
	else if(16<dt && dt<=24){
	    j=2;
	    Mjd0 = t1+8*j;
	}
	else if(24<dt && dt<=32){
	    j=3;
	    Mjd0 = t1+8*j;
	}

	double* r_MercuryCheb = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, &Cx_Mercury[14*j], &Cy_Mercury[14*j], &Cz_Mercury[14*j]);

	free(Cx_Mercury);
	free(Cy_Mercury);
	free(Cz_Mercury);

	*r_Mercury = producto_escalar(r_MercuryCheb, 3, 1e3);
	free(r_MercuryCheb);

	/* VENUS */

	double *Cx_Venus;
	double *Cy_Venus;
	double *Cz_Venus;

	sacarCX(171, 10, 201, 1, PCtemp, &Cx_Venus, &Cy_Venus, &Cz_Venus);
	if (0<=dt && dt<=16){
	    j=0;
	    Mjd0 = t1+16*j;
	}
	else if(16<dt && dt<=32){
	    j=1;
	    Mjd0 = t1+16*j;
	}

	double* r_VenusCheb = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, &Cx_Venus[10*j], &Cx_Venus[10*j], &Cx_Venus[10*j]);

	free(Cx_Venus);
	free(Cy_Venus);
	free(Cz_Venus);

	*r_Venus = producto_escalar(r_VenusCheb, 3, 1e3);
	free(r_VenusCheb);

	/* MARTE */

	double *Cx_Mars;
	double *Cy_Mars;
	double *Cz_Mars;

	sacarCX(309, 11, 342, 0, PCtemp, &Cx_Mars, &Cy_Mars, &Cz_Mars);

	j=0;
	Mjd0 = t1;

	double* r_MarsCheb = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, &Cx_Mars[11*j], &Cy_Mars[11*j], &Cz_Mars[11*j]);

	free(Cx_Mars);
	free(Cy_Mars);
	free(Cz_Mars);

	*r_Mars = producto_escalar(r_MarsCheb, 3, 1e3);
	free(r_MarsCheb);

	/* JUPITER */

	double *Cx_Jupiter;
	double *Cy_Jupiter;
	double *Cz_Jupiter;

	sacarCX(342, 8, 366, 0, PCtemp, &Cx_Jupiter, &Cy_Jupiter, &Cz_Jupiter);

	j=0;
	Mjd0 = t1;

	double* r_JupiterCheb = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, &Cx_Jupiter[8*j], &Cy_Jupiter[8*j], &Cz_Jupiter[8*j]);

	free(Cx_Jupiter);
	free(Cy_Jupiter);
	free(Cz_Jupiter);

	*r_Jupiter = producto_escalar(r_JupiterCheb, 3, 1e3);
	free(r_JupiterCheb);

	/* SATURNO */

	double *Cx_Saturn;
	double *Cy_Saturn;
	double *Cz_Saturn;

	sacarCX(366, 7, 387, 0, PCtemp, &Cx_Saturn, &Cy_Saturn, &Cz_Saturn);

	j=0;
	Mjd0 = t1;

	double* r_SaturnCheb = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, &Cx_Saturn[7*j], &Cy_Saturn[7*j], &Cz_Saturn[7*j]);

	free(Cx_Saturn);
	free(Cy_Saturn);
	free(Cz_Saturn);

	*r_Saturn = producto_escalar(r_SaturnCheb, 3, 1e3);
	free(r_SaturnCheb);

	/* URANO */

	double *Cx_Uranus;
	double *Cy_Uranus;
	double *Cz_Uranus;

	sacarCX(387, 6, 405, 0, PCtemp, &Cx_Uranus, &Cy_Uranus, &Cz_Uranus);

	j=0;
	Mjd0 = t1;

	double* r_UranusCheb = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, &Cx_Uranus[6*j], &Cy_Uranus[6*j], &Cz_Uranus[6*j]);

	free(Cx_Uranus);
	free(Cy_Uranus);
	free(Cz_Uranus);

	*r_Uranus = producto_escalar(r_UranusCheb, 3, 1e3);
	free(r_UranusCheb);

	/* NEPTUNO */

	double *Cx_Neptune;
	double *Cy_Neptune;
	double *Cz_Neptune;

	sacarCX(405, 6, 423, 0, PCtemp, &Cx_Neptune, &Cy_Neptune, &Cz_Neptune);

	j=0;
	Mjd0 = t1;

	double* r_NeptuneCheb = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, &Cx_Neptune[6*j], &Cy_Neptune[6*j], &Cz_Neptune[6*j]);

	free(Cx_Neptune);
	free(Cy_Neptune);
	free(Cz_Neptune);

	*r_Neptune = producto_escalar(r_NeptuneCheb, 3, 1e3);
	free(r_NeptuneCheb);

	/* PLUTON */

	double *Cx_Pluto;
	double *Cy_Pluto;
	double *Cz_Pluto;

	sacarCX(423, 6, 441, 0, PCtemp, &Cx_Pluto, &Cy_Pluto, &Cz_Pluto);

	j=0;
	Mjd0 = t1;

	double* r_PlutoCheb = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, &Cx_Pluto[6*j], &Cy_Pluto[6*j], &Cz_Pluto[6*j]);

	free(Cx_Pluto);
	free(Cy_Pluto);
	free(Cz_Pluto);

	*r_Pluto = producto_escalar(r_PlutoCheb, 3, 1e3);
	free(r_PlutoCheb);

	/* NUTACIONES */

//	double *Cx_Nutations = calloc(40, sizeof(double));
//	double *Cy_Nutations = calloc(40, sizeof(double));
//	double *Cz_Nutations = calloc(10, sizeof(double));
//
//	int n = 3;
//
//	int *temp = calloc(n, sizeof(int));
//
//	for (int j = 0; j < n; j++)
//		temp[j] = 818 + 10 * j; // Reducimos en 1 el valor para pasar a índices de C.
//
//
//	double *Cx = calloc(10, sizeof(double));
//	double *Cy = calloc(10, sizeof(double));
//
//	for (int j = 0; j < 10; j++) {
//		Cx_Nutations[j] = PCtemp[temp[0] + j];
//		Cy_Nutations[j] = PCtemp[temp[1] + j];
//	}
//	for (int i = 1; i <= 3; i++) {
//		for (int j = 0; j < 10; j++)
//			temp[j] += 20;
//
//		for (int j = 0; j < 10; j++) {
//			Cx[j] = PCtemp[temp[0] + j];
//			Cy[j] = PCtemp[temp[1] + j];
//		}
//
//		for (int j = 0; j < 10; j++) {
//			Cx_Nutations[j + 10 * i] = Cx[j];
//			Cy_Nutations[j + 10 * i] = Cy[j];
//		}
//	}
//
//	free(temp);
//	free(Cx);
//	free(Cy);
//
//	if (0 <= dt && dt <= 8) {
//		j = 0;
//	    Mjd0 = t1;
//	}
//	else if(8<dt && dt<=16){
//	    j=1;
//	    Mjd0 = t1+8*j;
//	}
//	else if(16<dt && dt<=24){
//	    j=2;
//	    Mjd0 = t1+8*j;
//	}
//	else if(24<dt && dt<=32){
//	    j=3;
//	    Mjd0 = t1+8*j;
//	}
//	double* Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, &Cx_Nutations[10*j], &Cy_Nutations[10*j], &Cz_Nutations);
//
//	free(Cx_Nutations);
//	free(Cy_Nutations);
//	free(Cz_Nutations);
//
//	/* LIBRACIONES */
//
//	double *Cx_Librations;
//	double *Cy_Librations;
//	double *Cz_Librations;
//
//	sacarCX(899, 10, 929, 3, PCtemp, &Cx_Librations, &Cy_Librations, &Cz_Librations);
//	if (0<=dt && dt<=8){
//	    j=0;
//	    Mjd0 = t1;
//	}
//	else if(8<dt && dt<=16){
//	    j=1;
//	    Mjd0 = t1+8*j;
//	}
//	else if(16<dt && dt<=24){
//	    j=2;
//	    Mjd0 = t1+8*j;
//	}
//	else if(24<dt && dt<=32){
//	    j=3;
//	    Mjd0 = t1+8*j;
//	}
//	double* Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, &Cx_Librations[10*j], &Cy_Librations[10*j], &Cz_Librations[10*j]);
//
//	free(Cx_Librations);
//	free(Cy_Librations);
//	free(Cz_Librations);

	double EMRAT = 81.30056907419062; // DE430
	double EMRAT1 = 1/(1+EMRAT);

	// Acabar
	for(int j = 0; j < 3; j++){
		(*r_Earth)[j] = (*r_Earth)[j]-EMRAT1*(*r_Moon)[j];
		(*r_Mercury)[j] = -(*r_Earth)[j]+(*r_Mercury)[j];
		(*r_Venus)[j] = -(*r_Earth)[j]+(*r_Venus)[j];
		(*r_Mars)[j] = -(*r_Earth)[j]+(*r_Mars)[j];
		(*r_Jupiter)[j] = -(*r_Earth)[j]+(*r_Jupiter)[j];
		(*r_Saturn)[j] = -(*r_Earth)[j]+(*r_Saturn)[j];
		(*r_Uranus)[j] = -(*r_Earth)[j]+(*r_Uranus)[j];
		(*r_Neptune)[j] = -(*r_Earth)[j]+(*r_Neptune)[j];
		(*r_Pluto)[j] = -(*r_Earth)[j]+(*r_Pluto)[j];
		(*r_Sun)[j] = -(*r_Earth)[j]+(*r_Sun)[j];
	}
}

double EccAnom(double M, double e){
	int maxit = 15;
	int i = 1;
	double E = 0;

	// Valor inicial
	M = fmod(M, 2.0*M_PI);

	if (e<0.8)
	    E = M;
	else
	    E = M_PI;

	double f = E - e*sin(E) - M;
	E = E - f / ( 1.0 - e*cos(E) );

	// Iteración
	while (abs(f) > 1e2 * DBL_EPSILON){
	    f = E - e*sin(E) - M;
	    E = E - f / (1.0 - e*cos(E));
	    i = i+1;
	    if (i == maxit)
	        printf(" problemas de convergencia en EccAnom");
	}

	return E;
}

double** GHAMatrix(double Mjd_UT1){
	return R_z(gast(Mjd_UT1));
}
double* AccelHarmonic(double* r, double** E, int n_max, int m_max){
	//extern double** Cnm;
	//extern double** Snm;

	double R_Earth   = 6378.1363e3;
	//double f_Earth   = 1/298.257223563;

	double GM_Earth    = 398600.4415e9;                  // [m^3/s^2]; DE430

	// Posición de cuerpos fijos.
	double* r_bf = producto(E,vect_a_matr(r),3,3,1)[0];

	// Cantidades auxiliares
	double d = norm(r_bf,3); // Distancia
	double latgc = asin((double)r_bf[2]/d);
	double lon = atan2(r_bf[1],r_bf[0]);

	double** pnm;
	double** dpnm;

	Legendre(n_max, m_max, latgc, &pnm, &dpnm);

	double dUdr = 0.0;
	double dUdlatgc = 0.0;
	double dUdlon = 0.0;
	double q3 = 0.0;
	double q2 = q3;
	double q1 = q2;
    double b1 = 0.0;
    double b2 = 0.0;
    double b3 = 0.0;
	for(int n=0;n < n_max+1;n++){
	    b1 = (-(double)GM_Earth/(d*d))*pow(((double)R_Earth/d),n)*(n+1);
	    b2 =  ((double)GM_Earth/d)*pow(((double)R_Earth/d),n);
	    b3 =  ((double)GM_Earth/d)*pow(((double)R_Earth/d),n);
	    for(int m=0;m < m_max+1;m++){
	        q1 = q1 + pnm[m][n]*(Cnm[m][n]*cos(m*lon)+Snm[m][n]*sin(m*lon));
	        q2 = q2 + dpnm[m][n]*(Cnm[m][n]*cos(m*lon)+Snm[m][n]*sin(m*lon));
	        q3 = q3 + m*pnm[m][n]*(Snm[m][n]*cos(m*lon)-Cnm[m][n]*sin(m*lon));
	    }
	    dUdr     += q1*b1;
	    dUdlatgc += q2*b2;
	    dUdlon   += q3*b3;
	    q3 = 0; q2 = q3; q1 = q2;
	}
	// Vamos liberando espacio.
	liberar(pnm, m_max+1);
	liberar(dpnm, m_max+1);

	// Aceleración de cuerpos fijos.
	double r2xy = r_bf[0]*r_bf[0]+r_bf[1]*r_bf[1];

	double ax = (1.0/d*dUdr-r_bf[2]/(d*d*sqrt((double)r2xy))*dUdlatgc)*r_bf[0]-(1.0/r2xy*dUdlon)*r_bf[1];
	double ay = (1.0/d*dUdr-r_bf[2]/(d*d*sqrt((double)r2xy))*dUdlatgc)*r_bf[1]+(1.0/r2xy*dUdlon)*r_bf[0];
	double az =  1.0/d*dUdr*r_bf[2]+sqrt((double)r2xy)/(d*d)*dUdlatgc;

	double* a_bf = (double*)calloc(3,sizeof(double));
	a_bf[0] = ax;
	a_bf[1] = ay;
	a_bf[2] = az;

//	printf("\n");
//	for(int jj = 0; jj<3; jj++)
//		printf("%3.13f\n",a_bf[jj]);
//	system("pause");

	// Aceleración por inercia
	double* res = producto(transpose(E,3,3),vect_a_matr(a_bf),3,3,1)[0];

	free(a_bf);
	return res;
}

double* AccelPointMass(double* r, double* s, double GM){
	// Vector de posición relativa del satélite con respecto del punto de masa
	double* d = *suma(vect_a_matr(r),vect_a_matr(producto_escalar(s, 3, -1.0)),3,1);

	// Aceleración
	double* res = producto_escalar(((double**)suma( vect_a_matr(producto_escalar(d, 3, 1.0/(norm(d,3)*norm(d,3)*norm(d,3)))), vect_a_matr(producto_escalar(s, 3, 1.0/(norm(s,3)*norm(s,3)*norm(s,3)))),3,1 ))[0], 3, -GM);


	free(d);

	return res;
}

double** G_AccelHarmonic(double* r, double** U, double n_max, double m_max){
	double d = 1.0;   // Incremento de posición [m]

	double** G = calloc(3,sizeof(double*));
	double* dr = calloc(3,sizeof(double));
	double* drsuma = calloc(3,sizeof(double));
	double* drresta = calloc(3,sizeof(double));
	double* a1 = calloc(3,sizeof(double));
	double* a2 = calloc(3,sizeof(double));
	double* da = calloc(3,sizeof(double));

	// Gradiente
	for(int i = 0; i<3; i++) {
		// Establecer offset en la componente i-ésima del vector posición r.
		dr[0] = 0.0;
		dr[1] = 0.0;
		dr[2] = 0.0;
	    dr[i] = d;

	    for(int j = 0; j < 3; j++){
	    	drsuma[j] = r[j]+dr[j]/2.0;
	    	drresta[j] = r[j]-dr[j]/2.0;
	    }

	    a1 = AccelHarmonic(drsuma, U, n_max, m_max);
	    a2 = AccelHarmonic(drresta, U, n_max, m_max);

	    for(int j = 0; j < 3; j++){
	    	da[j] = a1[j] - a2[j];
	    }

	    // Derivada con respecto al eje i-ésimo.
	    G[i] = producto_escalar(da, 3, 1.0/d);
	}

	free(dr);
	free(drsuma);
	free(drresta);
	free(a1);
	free(a2);

	return G;
}

double* Accel(double x, double* Y){
	double GM_Sun      = 132712440041.939400e9;            // [m^3/s^2]; DE430
	double GM_Moon     = 398600.435436e9/81.30056907419062;// [m^3/s^2]; DE430
	double GM_Mercury  = 22031.780000e9;                   // [m^3/s^2]; DE430
	double GM_Venus    = 324858.592000e9;                  // [m^3/s^2]; DE430
	double GM_Mars     = 42828.375214e9;                   // [m^3/s^2]; DE430
	double GM_Jupiter  = 126712764.800000e9;               // [m^3/s^2]; DE430
	double GM_Saturn   = 37940585.200000e9;                // [m^3/s^2]; DE430
	double GM_Uranus   = 5794548.600000e9;                 // [m^3/s^2]; DE430
	double GM_Neptune  = 6836527.100580e9;                 // [m^3/s^2]; DE430
	double GM_Pluto    = 977.0000000000009e9;              // [m^3/s^2]; DE430

	int nAux =      20;
	int mAux =      20;
	int sunAux =     1;
	int moonAux =    1;
	int planetsAux = 1;

	//extern double** eopdata;

	double MJD_J2000 = 51544.5;

	double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

	IERS(eopdata,Mjd_UTCAux+x/86400.0, 'l', &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

	double Mjd_UT1 = Mjd_UTCAux + x/86400.0 + UT1_UTC/86400.0;
	double Mjd_TT = Mjd_UTCAux + x/86400.0 + TT_UTC/86400.0;

	double** P = PrecMatrix(MJD_J2000, Mjd_TT);
	double** N = NutMatrix(Mjd_TT);

	double** T = producto(N,P,3,3,3);
	double **E = producto(producto(PoleMatrix(x_pole, y_pole),GHAMatrix(Mjd_UT1),3,3,3),T,3,3,3);

	double MJD_TDB = Mjday_TDB(Mjd_TT);

	double* r_Mercury;
	double* r_Venus;
	double* r_Earth;
	double* r_Mars;
	double* r_Jupiter;
	double* r_Saturn;
	double* r_Uranus;
	double* r_Neptune;
	double* r_Pluto;
	double* r_Moon;
	double* r_Sun;

	JPL_Eph_DE430(MJD_TDB, &r_Mercury, &r_Venus, &r_Earth, &r_Mars, &r_Jupiter, &r_Saturn, &r_Uranus, &r_Neptune, &r_Pluto, &r_Moon, &r_Sun);

	// Aceleración por campo gravitatorio armónico.
	double* a = AccelHarmonic(Y, E, nAux, mAux);

	// Perturbaciones luni-solares.
	if(sunAux){
		double* aHar = AccelPointMass(Y, r_Sun, GM_Sun);
		for(int j = 0; j < 3; j++)
			a[j] += aHar[j];
		free(aHar);
	}

	if (moonAux) {
		double *aHar = AccelPointMass(Y, r_Moon, GM_Moon);
		for (int j = 0; j < 3; j++)
			a[j] += aHar[j];
		free(aHar);
	}

	if(planetsAux){
		double *aHarMercury = AccelPointMass(Y, r_Mercury, GM_Mercury);
		double *aHarVenus = AccelPointMass(Y, r_Venus, GM_Venus);
		double *aHarMars = AccelPointMass(Y, r_Mars, GM_Mars);
		double *aHarJupiter = AccelPointMass(Y, r_Jupiter, GM_Jupiter);
		double *aHarSaturn = AccelPointMass(Y, r_Saturn, GM_Saturn);
		double *aHarUranus = AccelPointMass(Y, r_Uranus, GM_Uranus);
		double *aHarNeptune = AccelPointMass(Y, r_Neptune, GM_Neptune);
		double *aHarPluto = AccelPointMass(Y, r_Pluto, GM_Pluto);
		for (int j = 0; j < 3; j++){
			a[j] += aHarMercury[j];
			a[j] += aHarVenus[j];
			a[j] += aHarMars[j];
			a[j] += aHarJupiter[j];
			a[j] += aHarSaturn[j];
			a[j] += aHarUranus[j];
			a[j] += aHarNeptune[j];
			a[j] += aHarPluto[j];
		}
		free(aHarMercury);
		free(aHarVenus);
		free(aHarMars);
		free(aHarJupiter);
		free(aHarSaturn);
		free(aHarUranus);
		free(aHarNeptune);
		free(aHarPluto);
	}

	double* ret = calloc(6,sizeof(double));

	for(int j = 0; j < 3; j++){
		ret[j]=Y[j+3];
		ret[j+3]=a[j];
	}

	liberar(P, 3);
	liberar(N, 3);
	liberar(T, 3);
	liberar(E, 3);
	free(r_Mercury);
	free(r_Venus);
	free(r_Earth);
	free(r_Mars);
	free(r_Jupiter);
	free(r_Saturn);
	free(r_Uranus);
	free(r_Neptune);
	free(r_Pluto);
	free(r_Moon);
	free(r_Sun);
	free(a);

	return ret;
}

void AccelWrapper(double x, double* Y, /* Salidas */ double** Yp){
	*Yp = Accel(x,Y);
}
