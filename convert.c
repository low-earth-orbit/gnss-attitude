/*
gcc convert.c -lm -o convert
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "convert.h"

double deg2rad (double deg){
    return deg * M_PI / 180.0;
}

double rad2deg (double rad) {
    return rad / M_PI * 180.0;
}

/*
	Convert (x, y, z) vector to a unit vector
*/
void normalizeXyz (double *xyz) {
	double len = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[1]*xyz[1]);
	if (len != 0) {
		xyz[0] /= len;
		xyz[1] /= len;
		xyz[2] /= len;
	}
}

/*
	Input azimuth and elevation in degrees
	Get array of the resulting UNIT vector (x, y, z)
*/
void ae2xyz (double az_deg, double el_deg, double *xyz){
	double az = deg2rad(az_deg);
	double el = deg2rad(el_deg);
	double x = cos(el)*sin(az);
	double y = cos(el)*cos(az);
	double z = sin(el);
	xyz[0] = x;
	xyz[1] = y;
	xyz[2] = z;
}

/*
	Input x, y, z local coordinates known as E(x) N(y) U(x)
		The vector can be of unit or non-unit length
	Output azimuth and elevation angle in degrees
*/
void xyz2ae (double x, double y, double z, double* azel){
	double r = sqrt(x*x + y*y + z*z);
	double az = atan(x/y);
	double el = asin(z/r);
	/* adjust for quadrant */
	if (y < 0.0)	az += M_PI;
	if (az < 0.0)	az += 2*M_PI;
	double az_deg = rad2deg(az);
	double el_deg = rad2deg(el);
	azel[0] = az_deg;
	azel[1] = el_deg;
	
}

/*
	Calculate circular mean (azimuth)
*/
double meanAz (double array[], int n){
	if (n == 1)	return array[0];
	double mean, x, y, r;
	x = 0.0;
	y = 0.0;
	r = 0.0;
	for (int i = 0; i < n; i++){
		y += cos(deg2rad(array[i]));
		x += sin(deg2rad(array[i]));
	}
	y /= (double)n;
	x /= (double)n;
	r = sqrt(x*x + y*y); // keep this in mind. can be used for uniformity test.
	mean = atan(x/y);
	/* adjust for quadrant*/
	if (y < 0.0)	mean += M_PI;
	if (mean < 0.0)	mean += 2*M_PI;
	return rad2deg(mean);
}

/*
	Calculate spherical distance (great-circle distance) between two points on UNIT sphere.
	The output is the distance, also the angular distance in rad if the two points are on UNIT sphere.
*/
double spDist(double a1, double a2, double a3, double b1, double b2, double b3){
	//printf("a1 a2 a3 = %lf %lf %lf \t b1 b2 b3 = %lf %lf %lf \n", a1, a2, a3, b1, b2, b3);
	return acos(a1*b1 + a2*b2 + a3*b3);
}

/*
int main (void) {
	
	double x = 1;
	double y = 2;
	double z = 3;
	double len = sqrt(x*x+y*y+z*z);
	printf("x = %lf\ty = %lf\tz = %lf\t\n", x/len, y/len, z/len);
	double array[2];
	xyz2ae(x, y, z, array);
	//printf("az = %lf\tel = %lf\n", array[0], array[1]);
	
	double az = array[0];
	double el = array[1];
	double vector[3];
	ae2xyz(az, el, vector);
	printf("x = %lf\ty = %lf\tz = %lf\t\n", vector[0], vector[1], vector[2]);
	
	exit(0);
}
*/