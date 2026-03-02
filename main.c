/*
 ============================================================================
 Name        : Proyecto.c
 Author      : alcampo
 Version     :
 Copyright   : 
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include "rotations.h"
#include "matutils.h"
#include "mathutils.h"
#include "fechas.h"
#include "logica.h"
#include "ode.h"
#include "globales.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

double **Cnm;
double **Snm;
double **PC;
double **eopdata;

double Mjd_UTCAux;
double Mjd_TTAux;

int main(void) {

	//extern int nAux;
	//extern int mAux;
	//extern int sunAux;
	//extern int moonAux;
	//extern int planetsAux;

	//extern double MJD_J2000;
	//extern double Arcs;
	//extern double R_Earth;
	//extern double f_Earth;
	//extern double GM_Earth;
	double Rad       = M_PI/180;

	Cnm = zeros(181,181);
	Snm = zeros(181,181);
	PC = zeros(2285,1020);

	FILE *ptr;

	ptr = fopen("Data\\DE430Coeff.txt","r");

	if(ptr == NULL){
		printf("Error al abrir DE430Coeff.txt\n");
		exit(EXIT_FAILURE);
	}

	for(int n = 0; n < 2285; n++){
		for(int m = 0; m < 1020; m++){
			fscanf(ptr, "%lf", &PC[m][n]);
		}
	}

	fclose(ptr);

	printf("Leido DE430Coeff.txt\n");

	ptr = fopen("Data\\GGM03S.txt","r");

	if(ptr == NULL){
		printf("Error al abrir GGM03S.txt\n");
		exit(EXIT_FAILURE);
	}

	int n = 0, m = 0, f = 0, c = 0;
	double aux1 = 0.0, aux2 = 0.0;

	for(n = 0; n < 181; n++){
		for(m = 0; m <= n; m++){
			fscanf(ptr,"%d%d%lf%lf%lf%lf", &f, &c, &Cnm[m][n], &Snm[m][n], &aux1, &aux2);
		}
	}

	fclose(ptr);

	printf("Leido GGM03S.txt\n");

	eopdata = zeros(13,21413);

	ptr = fopen("Data\\eop19620101.txt","r");

	if(ptr == NULL){
		printf("Error al abrir eop19620101.txt\n");
		exit(EXIT_FAILURE);
	}

	for(f = 0; f < 21413; f++){
		fscanf(ptr,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &eopdata[f][0], &eopdata[f][1], &eopdata[f][2], &eopdata[f][3], &eopdata[f][4], &eopdata[f][5], &eopdata[f][6], &eopdata[f][7], &eopdata[f][8], &eopdata[f][9], &eopdata[f][10], &eopdata[f][11], &eopdata[f][12]);
	}

	fclose(ptr);

	printf("Leido eop19620101.txt\n");

	int nobs = 46;
	double** obs = zeros(nobs, 4);
	char line[54];
	int YY = 0, M = 0, D = 0, hh = 0, mm = 0;
	double ss = 0.0, az = 0.0, el = 0.0, Dist = 0.0;
	char y[5], mo[3], d[3], h[3], mi[3], s[7], a[9], e[9], di[10];

	ptr = fopen("Data\\GEOS3.txt","r");

	if(ptr == NULL){
		printf("Error al abrir GEOS3.txt\n");
		exit(EXIT_FAILURE);
	}

	for(c = 0; c < nobs; c++){
		fgets(line, sizeof(line)+2, ptr);

		strncpy(y, &line[0], 4);
		y[4] = '\0';
		YY = atoi(y);

		strncpy(mo, &line[5], 2);
		mo[2] = '\0';
		M = atoi(mo);

		strncpy(d, &line[8], 2);
		d[2] = '\0';
		D = atoi(d);

		strncpy(h, &line[12], 2);
		h[2] = '\0';
		hh = atoi(h);

		strncpy(mi, &line[15], 2);
		mi[2] = '\0';
		mm = atoi(mi);

		strncpy(s, &line[18], 6);
		s[6] = '\0';
		ss = atof(s);

		strncpy(a, &line[25], 8);
		a[8] = '\0';
		az = atof(a);

		strncpy(e, &line[35], 8);
		e[8] = '\0';
		el = atof(e);

		strncpy(di, &line[44], 9);
		di[9] = '\0';
		Dist = atof(di);

		obs[0][c] = Mjday(YY,M,D,hh,mm,ss);
		obs[1][c] = Rad * az;
		obs[2][c] = Rad * el;
		obs[3][c] = 1e3 * Dist;
	}

	fclose(ptr);

	printf("Leido GEOS3.txt\n");

	double sigma_range = 92.5, sigma_az = 0.0224 * Rad, sigma_el = 0.0139 * Rad;
	double lat = Rad * 21.5748, lon = Rad * (-158.2706), alt = 300.20;

	double* Rs = Position(lon, lat, alt); // Transpuesta

//	double Mjd1 = obs[0][0];
//	double Mjd2 = obs[0][8];
//	double Mjd3 = obs[0][17];

	double* Y0_apr = calloc(6,sizeof(double));
	Y0_apr[0] = 6221397.62857869;
	Y0_apr[1] =	2867713.77965741;
	Y0_apr[2] =	3006155.9850995;
	Y0_apr[3] =	4645.0472516175;
	Y0_apr[4] = -2752.21591588182;
	Y0_apr[5] = -7507.99940986939;

	double Mjd0 = Mjday(1995,1,29,02,38,0);

	double Mjd_UTC = obs[0][8];
	Mjd_UTCAux = Mjd_UTC;
	Mjd_TTAux = 0.0;
	double t = 0.0;
	int n_eqn = 6;
	double relerr = 1e-13;
	double abserr = 1e-6;
	int iflag = 1;
	double* work = malloc((100+21*n_eqn)*sizeof(double));
	double* workVar = malloc((100+21*42)*sizeof(double));
	int* iwork = malloc(5*sizeof(int));

	printf("Iniciando integracion previa...");

//	for(int j = 0; j < 6; j++)
//		printf("%f\n",Y0_apr[j]);

	ode(AccelWrapper, n_eqn, &Y0_apr, &t, -(obs[0][8]-Mjd0)*86400.0,relerr,abserr,&iflag,work,iwork);

	double* Y = calloc(6, sizeof(double));

	for(int j = 0; j < 6; j++)
		Y[j] = Y0_apr[j];

//	for(int j = 0; j < 6; j++)
//		printf("%f\n",Y[j]);

	iflag = 1;

	printf(" Hecho\n");

	double** P = zeros(6,6);
	for(int j = 0; j < 3; j++){
		P[j][j] = 1e8;
		P[j+3][j+3] = 1e3;
	}

	double** LT = LTC(lon, lat);

	double* yPhi = calloc(42,sizeof(double));
	double** Phi = zeros(6,6);

	double t_old = 0.0;
	double* Y_old = calloc(6,sizeof(6));

	t = 0.0;

	double tin = 0.0;
	double x_pole = 0.0;
	double y_pole = 0.0;
	double UT1_UTC = 0.0;
	double LOD = 0.0;
	double dpsi = 0.0;
	double deps = 0.0;
	double dx_pole = 0.0;
	double dy_pole = 0.0;
	double TAI_UTC = 0.0;

	double UT1_TAI = 0.0;
	double UTC_GPS = 0.0;
	double UT1_GPS = 0.0;
	double TT_UTC = 0.0;
	double GPS_UTC = 0.0;

	double Mjd_TT = 0.0;
	double Mjd_UT1 = 0.0;

	double theta = 0.0;
	double** U;
	double* r = calloc(3, sizeof(double));

	double* pos;
	double** R;
	double** Ur;

	double Azim = 0.0;
	double Elev = 0.0;
	double* dAds;
	double* dEds;
	double* dDds = calloc(3, sizeof(double));

	double* dAdY = calloc(6, sizeof(double));
	double* dEdY = calloc(6, sizeof(double));
	double* dDdY = calloc(6, sizeof(double));

	double** LTU;

	double** dAdsT;
	double** dEdsT;
	double** dDdsT;

	double** dAdsMult;
	double** dEdsMult;
	double** dDdsMult;

	printf("Empezando bucle principal\n");

	for(int i = 0;i < nobs;i++){
		printf("Iteracion: %d de %d\n",i+1,nobs);

		t_old = t;
		for(int j = 0; j < 6; j++)
			Y_old[j] = Y[j];

		Mjd_UTC = obs[0][i];			  // Fecha Juliana modificada.
		t		= (Mjd_UTC-Mjd0)*86400.0; // Tiempo desde la parte temporal. [s]

		IERS(eopdata, Mjd_UTC, 'l', &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);
		timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);

		Mjd_TT = Mjd_UTC + TT_UTC/86400;
		Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
		Mjd_UTCAux = Mjd_UTC;
		Mjd_TTAux = Mjd_TT;

		for(int ii = 1; ii <= 6; ii++){
			yPhi[ii-1] = Y_old[ii-1];
			for(int j = 1; j <= 6; j++){
				if(ii == j)
					yPhi[6*j+ii-1] = 1;
				else
					yPhi[6*j+ii-1] = 0;
			}
		}

		iflag = 1;

		n_eqn = 42;
		tin = 0.0;
		ode(VarEqnWrapper, n_eqn, &yPhi, &tin, t-t_old,relerr,abserr,&iflag,workVar,iwork);

		for(int j = 1; j <= 6; j++)
			for(int jj = 0; jj < 6; jj++)
				Phi[j-1][jj] = yPhi[6*j+jj];

		iflag = 1;

		n_eqn = 6;

		tin = 0.0;
		ode(AccelWrapper, n_eqn, &Y, &tin, t-t_old,relerr,abserr,&iflag,work,iwork);

		theta = gmst(Mjd_UT1);
		U = R_z(theta);

		for(int j = 0; j < 3; j++)
			r[j] = Y[j];
		R = vect_a_matr(r);
		Ur = producto(U, R, 3, 3, 1);
		liberar(Ur,1);
		liberar(R,1);
		for(int j = 0; j < 3; j++)
			Ur[0][j] = Ur[0][j]-Rs[j];
		pos = producto(LT,Ur,3,3,1)[0];

		TimeUpdate(&P, Phi, 6);

		AzElPa(pos, &Azim, &Elev, &dAds, &dEds);

		LTU = producto(LT,U,3,3,3);
		dAdsT = transpose(vect_a_matr(dAds), 3, 1);
		dAdsMult = producto(dAdsT,LTU, 1, 3, 3);
		liberar(dAdsT,3);
		for(int j = 0; j < 3; j++){
			dAdY[j] = dAdsMult[j][0];
			dAdY[j+3] = 0.0;
		}
		liberar(dAdsMult,3);
		MeasUpdate(&Y, obs[1][i], Azim, sigma_az, dAdY, &P, 6);

		for(int j = 0; j < 3; j++)
			r[j] = Y[j];
		R = vect_a_matr(r);
		Ur = producto(U, R, 3, 3, 1);
		liberar(R,1);
		for(int j = 0; j < 3; j++)
			Ur[0][j] = Ur[0][j]-Rs[j];
		free(pos);
		pos = producto(LT,Ur,3,3,1)[0];
		liberar(Ur,1);

		AzElPa(pos, &Azim, &Elev, &dAds, &dEds);

		liberar(LTU,3);
		LTU = producto(LT,U,3,3,3);
		dEdsT = transpose(vect_a_matr(dEds), 3, 1);
		dEdsMult = producto(dEdsT,LTU, 1, 3, 3);
		liberar(dEdsT,3);
		for(int j = 0; j < 3; j++){
			dEdY[j] = dEdsMult[j][0];
			dEdY[j+3] = 0.0;
		}
		liberar(dEdsMult,3);

		MeasUpdate(&Y, obs[2][i], Elev, sigma_el, dEdY, &P, 6);

		for(int j = 0; j < 3; j++)
			r[j] = Y[j];
		R = vect_a_matr(r);
		liberar(Ur,1);
		Ur = producto(U, R, 3, 3, 1);
		liberar(R,1);
		for(int j = 0; j < 3; j++)
			Ur[0][j] = Ur[0][j]-Rs[j];
		free(pos);
		pos = producto(LT,Ur,3,3,1)[0];

		Dist = norm(pos, 3);
		for(int j = 0; j < 3; j++)
			dDds[j] = pos[j]/Dist;
		liberar(LTU,3);
		LTU = producto(LT,U,3,3,3);
		dDdsT = transpose(vect_a_matr(dDds), 3, 1);
		dDdsMult = producto(dDdsT, LTU, 1, 3, 3);
		liberar(dDdsT, 3);
		for (int j = 0; j < 3; j++) {
			dDdY[j] = dDdsMult[j][0];
			dDdY[j + 3] = 0.0;
		}
		liberar(dDdsMult, 3);

		MeasUpdate(&Y, obs[3][i], Dist, sigma_range, dDdY, &P, 6); // Culpable

		liberar(LTU,3);
		liberar(U,3);
		free(pos);
	}

	printf("Acabando...\n");

	IERS(eopdata,obs[0][45], 'l', &x_pole, &y_pole, &UT1_UTC, &LOD, &dpsi, &deps, &dx_pole, &dy_pole, &TAI_UTC);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	Mjd_UTCAux = Mjd_UTC;
	Mjd_TTAux = Mjd_TT;

	iflag = 1;

	n_eqn = 6;
	tin = 0.0;
	ode(AccelWrapper, n_eqn, &Y, &tin,-(obs[0][45]-obs[0][0])*86400.0,relerr,abserr,&iflag,work,iwork);

	double* Y_true = calloc(6,sizeof(double));
	Y_true[0] = 5753.173e3;
	Y_true[1] =	2673.361e3;
	Y_true[2] =	3440.304e3;
	Y_true[3] =	4.324207e3;
	Y_true[4] = -1.924299e3;
	Y_true[5] = -5.728216e3;

	printf("\nError de Estimación de Posición\n");
	printf("dX%10.1f [m]\n",Y[0]-Y_true[0]);
	printf("dY%10.1f [m]\n",Y[1]-Y_true[1]);
	printf("dZ%10.1f [m]\n",Y[2]-Y_true[2]);
	printf("\nError de Estimación de Velocidad\n");
	printf("dVx%8.1f [m/s]\n",Y[3]-Y_true[3]);
	printf("dVy%8.1f [m/s]\n",Y[4]-Y_true[4]);
	printf("dVz%8.1f [m/s]\n",Y[5]-Y_true[5]);

	system("PAUSE");

	return EXIT_SUCCESS;
}
