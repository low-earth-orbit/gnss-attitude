// convert.h
#ifndef CONVERT_H
#define CONVERT_H

/*
	Function declarations
*/
double deg2rad (double deg);
double rad2deg (double rad);
void ae2xyz (double az_deg, double el_deg, double *xyz);
void xyz2ae (double x, double y, double z, double* azel);

#endif