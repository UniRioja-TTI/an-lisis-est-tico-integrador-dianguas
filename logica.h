/*
 * logica.h
 *
 *  Created on: 1 abr 2021
 *      Author: Administrador
 */

#ifndef LOGICA_H_
#define LOGICA_H_

double** LTC(double lon, double lat); // Devuelve una matriz de rotación para pasar del sistema de coordenadas de longitud y latitud al sistema (Este, Norte, Cénit) tangente local.
double** PoleMatrix(double xp, double yp); // Devuelve una matriz para transformar de coordenadas pseudofijas terrestres a coordenadas fijas terrestres.
double** PrecMatrix(double Mjd_1, double Mjd_2); // Transformación de precesión de coordenadas ecuatoriales.
double* Position(double lon, double lat, double h); // Devuelve el vector posición a partir de la longitud, latitud y altitud.
void Geodetic(double* r,/* Salida */ double* lon, double* lat, double* h); // Devuelve las coordenadas geodéticas (longitud, latitud, altura) a partir del vector posición r.
void MeasUpdate(double** x, double z, double g, double s, double* G, double*** P, int n);

// En las cuatro siguientes funciones, Mjd_TT representa la Fecha Juliana Modificada en Tiempo Terrestre
void NutAngles(double Mjd_TT, /* Salida */ double* dpsi, double* deps); // Devuelve la nutación en longitud y oblicuidad.
double MeanObliquity(double Mjd_TT); // Devuelve la oblicuidad media de la eclíptica
double** NutMatrix(double Mjd_TT); // Devuelve la matriz que transforma el ecuador y equinocio medios a reales.
double EqnEquinox(double Mjd_TT); // Devuelve el resultado de la computación de la ecuación de los equinocios.
void AzElPa(double* s, /* Salidas */ double* Az, double* El, double** dAds, double** dEds); // Computa azimut, elevación y derivadas parciales de las coordenadas tangentes locales.

void IERS(double** eop,    double Mjd_UTC,  char interp,
		  /* Salida */
		  double *x_pole,  double *y_pole,  double *UT1_UTC,
		  double *LOD,     double *dpsi,    double *deps,
		  double *dx_pole, double *dy_pole, double *TAI_UTC);  // Función de lectura de hora IERS y datos de movimientos polares.
void timediff(double UT1_UTC, double TAI_UTC,
		/* Salida */ double* UT1_TAI, double* UTC_GPS, double* UT1_GPS, double* TT_UTC, double* GPS_UTC); // Diferencias temporales. [s]

void sacarCX(int inicio, int salto, int fin, int itBucle, double* PCtemp, /* Salidas */ double** Cx_Ret, double** Cy_Ret, double** Cz_Ret); // Automatización de trámites repetitivos
void JPL_Eph_DE430(double Mjd_TDB, /* Salidas */ double** r_Mercury, double** r_Venus, double** r_Earth, double** r_Mars, double** r_Jupiter, double** r_Saturn, double** r_Uranus,
		double** r_Neptune, double** r_Pluto, double** r_Moon, double** r_Sun); // Computa la posición ecuatorial del Sol, de la Luna y de los nueve cuerpos principales utilizando las efemérides de JPL.

double EccAnom(double M, double e); // Devuelve la anomalía eccéntrica en órbitas elípticas a partir de la anomalía media M en radianes y la eccentricidad e.
double** GHAMatrix(double Mjd_UT1); // Devuelve la matriz que transforma del verdadero ecuador y equinocio a ecuador terrestre y sistema meridiano de Greenwich.

double* AccelHarmonic(double* r, double** E, int n_max, int m_max); // Computa la aceleración debida al campo gravitatorio armónico del cuerpo central.
double* AccelPointMass(double* r, double* s, double GM); // Computa la aceleración perturbacional debida a un punto de masa.

double** G_AccelHarmonic(double* r, double** U, double n_max, double m_max); // Computa el gradiente del campo gravitatorio armónico de La Tierra.
double* Accel(double x, double* Y); // Computa la aceleración de un satélite orbitando la Tierra en función del campo gravitatorio armónico de La Tierra,
									// las perturbaciones gravitacionales del Sol y la Luna, la presión por radiación solar y la fuerza de arrastre de la atmósfera.
void AccelWrapper(double x, double* Y, /* Salidas */ double** Yp);
#endif /* LOGICA_H_ */
