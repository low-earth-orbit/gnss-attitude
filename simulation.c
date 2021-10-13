#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "convert.h" // Coordinate transformation methods
/*
This program simulates the satellite-antenna geometry.
	a. Only the satellites above Earth's horizon can be possibly seen regardless of the direction in which the antenna is pointing. This is true when the antenna is located on the Earth's surface.
	b. Antenna can only receive signals from satellites above the antenna's horizon. This is typical for GNSS antennas receiving signals. Of course, the user may set a cutoff angle for post-processing.
Therefore, there are two half spheres. When the antenna is pointing up at a 90 deg elevation angle, the two half spheres coincide. However, when the antenna is not at a 90 deg elevation angle only the area these two half-sphere overlap has visible GNSS satellites.
The underlying assumption of this randomness simulation is that at a given moment, the satellites in the sky are randomly (i.e. uniformly) distributed.

To compile: gcc simulation.c convert.c -lm -o simulation
*/

/*
	Generate points visible by the antenna, which are antenna boresight angle dependent. This function assumes antenna azimuth = 180 deg. That is pointing to negative y (northing).
	
		points[][5]: points array
		
		antEl: antenna boresight elevation angle in degrees; if the antenna is pointing straight up, that is 90 deg.
		
		n: number of points on the COMPLETE sphere
		
		returns: number of visible points
	
*/
/*
	Currently this way is not really ramdonly uniform distributed. Greater point density around pole location
*/
int rpVisSat (double points[][5], int n, double antEl)
{
	int numVisPt = 0;
	double theta, phi, x, y, z, c;
	double azel[2];
	for (int i = 1; i < n; i ++){
		theta = 2 * M_PI * (double)rand() / (double)RAND_MAX;
		phi = M_PI *((double)rand() / (double)RAND_MAX);
		x = sin(phi)*cos(theta);
		y = sin(phi)*sin(theta);
		z = cos(phi); // Treat this as local coordinates xyz
		xyz2ae (x, y, z, azel);
		c = tan(deg2rad(antEl));
		// if in visible area, record the point
		if ((z > 0) && (-y + c*z > 0)) {
			points[numVisPt][0] = x;
			points[numVisPt][1] = y;
			points[numVisPt][2] = z;
			points[numVisPt][3] = azel[0];
			points[numVisPt][4] = azel[1];
			numVisPt++;
		}
	}
	return numVisPt;
}

int main (void) {
	int numPt = 10000;
	double visPoints[numPt][5];
	int numVisPt = rpVisSat(visPoints, numPt, 45);
	for (int i = 1; i < numVisPt; i ++){
		//printf("xyz = %lf, %lf, %lf || azel = %lf, %lf\n", visPoints[i][0], visPoints[i][1], visPoints[i][2], visPoints[i][3], visPoints[i][4]);
		printf("%lf %lf %lf\n", visPoints[i][0], visPoints[i][1], visPoints[i][2]);
	}
	
	/*exit*/
	exit(0);
} // end of main()