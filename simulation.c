#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "convert.h" // Coordinate transformation methods
/*
This program simulates the satellite-antenna geometry.
	a. Only the satellites above Earth's horizon can be possibly seen regardless of the direction at which the antenna is pointing. This is true when the antenna is located on the Earth's surface.
	b. Antenna can only receive signals from satellites above antenna's horizon. This is typical for GNSS antennas receiving signals. Of course, the user may set a cutoff angle for post-processing.
Therefore, there are two half spheres. When the antenna is pointing up at 90 deg elevation angle, the two half spheres coincide. However, when the antenna is not at 90 deg elevation angle only the area these two half sphere overlap have visible GNSS satellites.
The underlying assumption of this randomness simulation is that at given momoent, the satellites in sky are randomly (i.e. uniformly) distributed.

To compile: gcc simulation.c convert.c -lm -o simulation
*/

/*
	Generate satSphere; that is the half spheral above local Earth's horizon.
*/
void rpSat (double points[][5], int n)
{
	for (int i = 1; i < n; i ++){
		double u, v, theta, phi, x, y, z;
		do {
			u = ((double) rand() / RAND_MAX);
			v = ((double) rand() / RAND_MAX);
			//printf("u = %lf, v = %lf\n", u, v);
			theta = 2 * M_PI * u;
			phi = acos (2*v - 1); // lat and long
			//printf("theta = %lf, phi = %lf\n", theta, phi);
			x = sin(theta)*cos(phi);
			y = sin(theta)*sin(phi);
			z = cos(theta); // x, y, z ECEF
							// It's ok we can treat this as local coordinates xyz*/
		} while (z <= 0); // cut off to make a semisphere
		double azel[2];
		xyz2ae (x, y, z, azel);
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;
		points[i][3] = azel[0];
		points[i][4] = azel[1];
	}
}

/*
	Generate visible satellite sphere; that is antenna boresight angle depedent.
	Input boresightAngle in degrees; if antenna is pointing straight up, that is 90 deg.
	This function assumes antenna azimuth = 180 deg. That is pointing to negative y (northing).
*/
void rpVisSat (double points[][5], int n, double boresightAngle)
{
	for (int i = 1; i < n; i ++){
		double cutoff = 90.0 - boresightAngle;
		double u, v, theta, phi, x, y, z;
		double azel[2];
		do {
			u = ((double) rand() / RAND_MAX);
			v = ((double) rand() / RAND_MAX);
			//printf("u = %lf, v = %lf\n", u, v);
			theta = 2 * M_PI * u;
			phi = acos (2*v - 1); // lat and long
			//printf("theta = %lf, phi = %lf\n", theta, phi);
			x = sin(theta)*cos(phi);
			y = sin(theta)*sin(phi);
			z = cos(theta); // x, y, z ECEF
							// It's ok we can treat this as local coordinates xyz
			xyz2ae (x, y, z, azel);
		} while ((z <= 0) || (y > 0 && azel[1] <= cutoff)); // cut off to make a semisphere
		points[i][0] = x;
		points[i][1] = y;
		points[i][2] = z;
		points[i][3] = azel[0];
		points[i][4] = azel[1];
	}
}

int main (void) {
	int numPt = 1000;
	double points[numPt][5];
	rpSat(points, numPt);
	for (int i = 1; i < numPt; i ++){
		printf("xyz = %lf, %lf, %lf || azel = %lf, %lf\n", points[i][0], points[i][1], points[i][2], points[i][3], points[i][4]);
	}
	
	int numVisPt = 1000;
	double visPoints[numVisPt][5];
	rpSat(visPoints, numVisPt);
	for (int i = 1; i < numVisPt; i ++){
		printf("xyz = %lf, %lf, %lf || azel = %lf, %lf\n", visPoints[i][0], visPoints[i][1], visPoints[i][2], visPoints[i][3], visPoints[i][4]);
	}
	
	
	

	
	/*
		exit
	*/
	exit(0);
} // end of main()