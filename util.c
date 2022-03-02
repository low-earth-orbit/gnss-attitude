/*
gcc util.c struct.c -lm -o util
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "util.h"
#include "struct.h"

/* comparison function for qsort satArray */
int cmpSatArray(const void *a, const void *b)
{
	const Sat *a1 = *(const Sat **)a;
	const Sat *b1 = *(const Sat **)b;
	return strcmp(a1->time, b1->time);
}

bool isStrInArray(char *str, char **array, long int index)
{
	for (long int i = 0; i < index; i++)
	{
		if (strcmp(str, array[i]) == 0)
			return true;
	}
	return false;
}

char *concat(const char *str1, const char *str2)
{
	char *str;
	str = (char *)malloc(strlen(str1) + strlen(str2) + 2);
	// +1 for the null-terminator; +1 for the seperator
	if (str == NULL)
	{
		fprintf(stderr, "malloc() failed in concat()\n");
	}
	else
	{
		strcpy(str, str1);
		strcat(str, " ");
		strcat(str, str2);
	}
	return str;
}

double deg2rad(double deg)
{
	return deg * M_PI / 180.0;
}

double rad2deg(double rad)
{
	return rad / M_PI * 180.0;
}

/*
	Convert (x, y, z) vector to a unit vector
*/
void normalizeXyz(double *xyz)
{
	double len = sqrt(pow(xyz[0], 2) + pow(xyz[1], 2) + pow(xyz[2], 2));
	if (len != 0)
	{
		xyz[0] /= len;
		xyz[1] /= len;
		xyz[2] /= len;
	}
}

/*
	normalize xyz in Sol
*/
void normalize(Sol *sol)
{
	double len = sqrt(pow(*sol->x, 2) + pow(*sol->y, 2) + pow(*sol->z, 2));
	// printf("%lf,%lf,%lf\n", *sol->x, *sol->y, *sol->z);
	// printf("length = %lf\n", len);
	if (len != 0)
	{
		*sol->x /= len;
		*sol->y /= len;
		*sol->z /= len;
	}
}

/*
	Input: azimuth and elevation in degrees
	Output: array of the resulting UNIT vector (x, y, z) i.e. (e, n, u)
*/
void ae2xyz(double az_deg, double el_deg, double *xyz)
{
	double az = deg2rad(az_deg);
	double el = deg2rad(el_deg);
	double x = cos(el) * sin(az);
	double y = cos(el) * cos(az);
	double z = sin(el);
	xyz[0] = x;
	xyz[1] = y;
	xyz[2] = z;
}

void ae2xyzSol(double az_deg, double el_deg, Sol *sol)
{
	double az = deg2rad(az_deg);
	double el = deg2rad(el_deg);
	double x = cos(el) * sin(az);
	double y = cos(el) * cos(az);
	double z = sin(el);
	*sol->x = x;
	*sol->y = y;
	*sol->z = z;
}

/*
	Input x, y, z local coordinates, E(x) N(y) U(x)
		Don't mistake this x, y, z as ECEF coordinates
	Output azimuth and elevation angle in degrees
*/
void xyz2ae(double x, double y, double z, double *azel)
{
	double az, el;
	double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)); // the vector can be of unit or non-unit length

	if (y == 0 && x > 0)
	{
		azel[0] = 90.0;
		azel[1] = 0.0;
	}
	else if (y == 0 && x < 0)
	{
		azel[0] = 270.0;
		azel[1] = 0.0;
	}
	else if (y == 0 && x == 0)
	{
		azel[0] = -999; // undefined
		azel[1] = 0.0;
	}
	else
	{
		az = atan2(x, y);
		/* adjust for quadrant */
		if (az < 0)
		{
			az = az + 2.0 * M_PI;
		}
		az = rad2deg(az);

		el = asin(z / r); // r should not be 0
		el = rad2deg(el);
		azel[0] = az;
		azel[1] = el;
	}
}

/*
	Input local coordinates, E(x) N(y) U(x)
		Don't mistake this x, y, z as ECEF coordinates
	Sets the values (azimuth and elevation angle in degrees) in Sol object
*/
void xyz2aeSol(double x, double y, double z, Sol *sol)
{
	double az, el;
	double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	if (y == 0 && x > 0)
	{
		*(sol->az) = 90.0;
		*(sol->el) = 0.0;
	}
	else if (y == 0 && x < 0)
	{
		*(sol->az) = 270.0;
		*(sol->el) = 0.0;
	}
	else if (y == 0 && x == 0)
	{
		*(sol->az) = -999; // undefined
		*(sol->el) = 0.0;
	}
	else
	{
		az = atan2(x, y);
		/* adjust for quadrant */
		if (az < 0)
		{
			az = az + 2.0 * M_PI;
		}
		az = rad2deg(az);

		el = asin(z / r); // r should not be 0
		el = rad2deg(el);
		*(sol->az) = az;
		*(sol->el) = el;
	}
}

/*
	Calculate circular mean of a set of angles (in degree)
*/
double cirMean(double array[], int n)
{
	if (n == 1)
		return array[0];
	double mean, x, y;
	x = 0.0;
	y = 0.0;
	for (int i = 0; i < n; i++)
	{
		y += cos(deg2rad(array[i]));
		x += sin(deg2rad(array[i]));
	}
	y /= (double)n;
	x /= (double)n;
	if (x / y == M_PI / 2)
	{
		return -999;
		printf("Found x/y = pi/2 in meanAz(). Returned -999 Undefined.\n");
	}
	mean = atan2(x, y);
	/* adjust for quadrant */
	if (mean < 0)
	{
		mean += 2.0 * M_PI;
	}
	return rad2deg(mean);
}

/*
	Calculate circular variance of a set of angles (in degree)
*/
double cirVar(double array[], int n)
{
	if (n == 1)
		return array[0];
	double x, y;
	x = 0.0;
	y = 0.0;
	for (int i = 0; i < n; i++)
	{
		y += cos(deg2rad(array[i]));
		x += sin(deg2rad(array[i]));
	}
	y /= (double)n;
	x /= (double)n;
	double r = sqrt(x * x + y * y) / (double)n;
	double var = 1.0 - r;
	return var;
}

/*
	Calculate circular variance of azimuth
	*epoch: pointer to Epoch object
	output: circular variance of azimuths
*/
double cirStdAzEpoch(Epoch *epoch)
{
	int n = *(epoch->numSat);
	double x, y;
	x = 0.0;
	y = 0.0;
	for (int i = 0; i < n; i++)
	{
		y += cos(deg2rad(*epoch->epochSatArray[i]->az));
		x += sin(deg2rad(*epoch->epochSatArray[i]->az));
	}
	y /= (double)n;
	x /= (double)n;
	double r = sqrt(x * x + y * y) / (double)n;
	double std = sqrt(-2.0 * log(r));
	return std;
}

/*
	Calculate spherical standard deviation of a set of cartesian coordinates (x, y, z) i.e. (e, n, u) in the epoch
	*epoch: pointer to Epoch object
	output: spherical standard deviation
*/
double spStdEpoch(Epoch *epoch)
{
	int n = *(epoch->numSat);
	double x, y, z;
	x = 0.0;
	y = 0.0;
	z = 0.0;
	for (int i = 0; i < n; i++)
	{
		double xyz[3];
		ae2xyz(*epoch->epochSatArray[i]->az, *epoch->epochSatArray[i]->el, xyz);
		x += xyz[0];
		y += xyz[1];
		z += xyz[2];
	}
	x /= (double)n;
	y /= (double)n;
	z /= (double)n;
	double rsq = (x * x + y * y + z * z);
	double std = sqrt(1 - rsq);
	return std;
}

/*
	Calculate spherical distance (great-circle distance) between two points on UNIT sphere.
	The output is the angular distance in rad if the two points are on UNIT sphere.
*/
double spDist(double a1, double a2, double a3, double b1, double b2, double b3)
{
	// printf("a1 a2 a3 = %lf %lf %lf \t b1 b2 b3 = %lf %lf %lf \n", a1, a2, a3, b1, b2, b3);
	return acos(a1 * b1 + a2 * b2 + a3 * b3);
}

/*
int main (void) {

	double x = 1;
	double y = 2;
	double z = 3;
	double len = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
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
