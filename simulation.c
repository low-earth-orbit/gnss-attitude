#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "convert.h" // Coordinate transformation methods
/*
This program simulates the satellite-antenna geometry.

a. Only the satellites above Earth's horizon can be possibly seen regardless of the direction in which the antenna is pointing. This is true when the antenna is located on or close to the Earth's surface. The user may set a cutoff angle for post-processing. Assume such angle is 0 here.

b. Antenna can only receive signals from satellites above the antenna's horizon. This is typical for GNSS antennas receiving signals.

Therefore, there are two hemispheres. When the antenna is pointing up at a 90 deg elevation angle, the two coincide. However, when the antenna is not at a 90 deg elevation angle only the area these hemispheres overlap has visible GNSS satellites.

The underlying assumption of this simulation is that GNSS satellites are randomly distributed on the celestial sphere.

gcc simulation.c convert.c -lm -o simulation
./simulation > input.txt
*/

/*
	Generate points visible by the antenna, which are boresight angle dependent. This function assumes antenna azimuth = 180 deg. That is pointing to negative y (northing). Reference: https://mathworld.wolfram.com/SpherePointPicking.html
	
		points[][5]: points array
		
		antEl: antenna's elevation angle (a.k.a. boresight angle) in degrees; if the antenna is pointing straight up, that is 90 deg.
		
		n: number of points on the ENTIRE sphere
		
		returns: number of satellites visible to the antenna with specific elevation angle
	
*/

typedef struct {
	double x; 
	double y;
	double z;// x (E), y(N), z(U) in local coordinates
	double az;
	double el;
	double snr;
} SimuSat; // simulated satellite


int rpVisSat (SimuSat* sat, int n, double antEl)
{
	int numVisPt = 0;
	double theta, phi, x, y, z, c;
	double azel[2];
	for (int i = 1; i < n; i ++){
		theta = 2.0 * M_PI * (double)rand() / (double)RAND_MAX;
		phi = acos(2.0 * (double)rand() / (double)RAND_MAX -1);
		x = sin(phi)*cos(theta);
		y = sin(phi)*sin(theta);
		z = cos(phi); // Treat this as local coordinates xyz
		xyz2ae (x, y, z, azel);
		c = tan(deg2rad(antEl));
		// if in visible area, record the point
		if ((z > 0) && (-y + c*z > 0)) {
			sat[numVisPt].x = x;
			sat[numVisPt].y = y;
			sat[numVisPt].z = z;
			sat[numVisPt].az = azel[0];
			sat[numVisPt].el= azel[1];
			sat[numVisPt].snr = 0;// snr will be updated in main()
			numVisPt++;
		}
	}
	return numVisPt;
}

int main (void) {
	int antEl = 60; // Antenna boresight elevation angle
	// While elevation angle is adjustable, antenna azimuth is simulated at 180 deg by rpVisSat()
	// Boresight vector is (0, -cos(antEl), sin(antEl))
	int numEpoch = 1000; // number of simulated epoch
	int numSat = 100; // number of GNSS satellites globally available
	const double MAX_SNR = 50;
	const double MIN_SNR = 20;// Set max and min snr values for snr computation. Assume linear relationship between SNR and (off-)boresight angle
	
	double spd, snr;
	int numVisPt;
	printf("SIMULATED INPUT FILE || \"SIMUEPOCH#\" \"TIME\"Epoch# \"SAT\" AZ EL SNR SAT#\n");// print header 	
	for (int i=0; i<numEpoch; i++){ // one simulation per loop 
		
		SimuSat* visSat = (SimuSat*)malloc(numSat*sizeof(SimuSat));
		if (visSat == NULL) {
			fprintf(stderr, "malloc() failed for creating visSat\n");
			exit(-1);
		}
		srand(time(0)+i);// set seed for rand()
		numVisPt = rpVisSat(visSat, numSat, antEl);
		if (numVisPt !=0){
			for (int j = 0; j < numVisPt; j ++){
				/* compute SNR*/
				spd = spDist(visSat[j].x, visSat[j].y, visSat[j].z, 0, -cos(deg2rad(antEl)), sin(deg2rad(antEl)));
				snr = MAX_SNR - ( (M_PI*0.5 - spd) / (0.5*M_PI))*(MAX_SNR - MIN_SNR);
				visSat[j].snr = snr;
				/*
					print
				*/
				//printf("xyz = %lf, %lf, %lf || azel = %lf, %lf || snr = %lf\n", visSat[j].x, visSat[j].y, visSat[j].z, visSat[j].az, visSat[j].el, visSat[j].snr); // for check
				//printf("%lf %lf %lf\n", visSat[j].x, visSat[j].y, visSat[j].z);// for plot
				printf("SIMUEPOCH# TIME%06d SAT %lf %lf %lf %i\n", i, visSat[j].az, visSat[j].el, visSat[j].snr, j); // for output as input file
			}
		}
		
		/* free after use*/
		free(visSat);
		
	}
	/*exit*/
	exit(0);
} // end of main()