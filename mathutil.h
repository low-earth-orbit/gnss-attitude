// mathutil.h
#ifndef MATHUTIL_H
#define MATHUTIL_H

/*
	Function declarations
*/
double deg2rad (double deg);
double rad2deg (double rad);
void normalizeXyz (double *xyz);
void ae2xyz (double az_deg, double el_deg, double *xyz);
void xyz2ae (double x, double y, double z, double* azel);
double meanAz (double array[], int n);
double spDist(double a1, double a2, double a3, double b1, double b2, double b3);
#endif
