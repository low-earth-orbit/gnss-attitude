#ifndef MATHUTIL_H
#define MATHUTIL_H
#include "struct.h"
/*
	Function declarations
*/
double deg2rad(double deg);
double rad2deg(double rad);
void normalizeXyz(double *xyz);
void normalize(Sol *sol);
void ae2xyz(double az_deg, double el_deg, double *xyz);
void ae2xyzSol(double az_deg, double el_deg, Sol *sol);
void xyz2ae(double x, double y, double z, double *azel);
void xyz2aeSol(double x, double y, double z, Sol *sol);
double meanAz(double array[], int n);
double spDist(double a1, double a2, double a3, double b1, double b2, double b3);
#endif
