/*
 * fechas.h
 *
 *  Created on: 1 abr 2021
 *      Author: Administrador
 */

extern double MJD_J2000;

#ifndef FECHAS_H_
#define FECHAS_H_

double gmst(double Mjd_UT1); // Devuelve el Greenwich Mean Sidereal Time a partir de la fecha Juliana UT1 modificada
double gast(double Mjd_UT1); // Devuelve el Greenwich Apparent Sidereal Time a partir de la fecha Juliana UT1 modificada
double Mjday(int yr, int mon, int day, int hr, int min, int sec); // Devuelve la fecha Juliana modificada.
//double Mjday(int yr, int mon, int day);
double Mjday_TDB(double Mjd_TT); // Devuelve la fecha Juliana modificada para la hora baricéntrica dinámica.
void TimeUpdate(double*** P, double** Phi, int n);

#endif /* FECHAS_H_ */
