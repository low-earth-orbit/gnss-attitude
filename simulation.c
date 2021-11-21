#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "util.h"	// math utilities
#include "config.h" // configuration
#include "struct.h" // structures
#include "snr.h"	// snr adjustment & mapping

/*
This program simulates the satellite-antenna geometry.

a. Only the satellites above Earth's horizon can be possibly seen regardless of the direction in which the antenna is pointing. This is true when the antenna is located on or close to the Earth's surface. The user may set a cutoff angle for post-processing. Assume such angle is 0 here.

b. Antenna can only receive signals from satellites above the antenna's horizon. This is typical for GNSS antennas receiving signals.

Therefore, there are two hemispheres. When the antenna is pointing up at a 90 deg elevation angle, the two coincide. However, when the antenna is not at a 90 deg elevation angle only the area these hemispheres overlap has visible GNSS satellites.

The assumption of this simulation is that GNSS satellites are randomly distributed on the celestial sphere.

gcc -Wall simulation.c util.c struct.c -o simulation  -lgsl -lgslcblas -lm
./simulation > input.txt
*/

typedef struct SimuSat
{
	double x;
	double y;
	double z; // x (E), y(N), z(U) in local coordinates
	double az;
	double el;
	double snr;
} SimuSat; // simulated satellite

/*
	Generate points visible by the antenna, which are boresight angle dependent. This function assumes antenna azimuth = 180 deg. That is pointing to negative y (northing). Reference: https://mathworld.wolfram.com/SpherePointPicking.html
	
		sat: pointer to SimuSat object
		
		antEl: antenna's elevation angle (a.k.a. boresight angle) in degrees; if the antenna is pointing straight up, that is 90 deg.
		
		n: number of points on the ENTIRE sphere
		
		returns: number of satellites visible to the antenna with specific elevation angle
*/
int randSat(SimuSat *sat, int n, double antEl)
{
	int numVisPt = 0;
	double theta, phi, x, y, z, c;
	double azel[2];
	for (int i = 1; i < n; i++)
	{
		theta = 2.0 * M_PI * (double)rand() / (double)RAND_MAX;
		phi = acos(2.0 * (double)rand() / (double)RAND_MAX - 1);
		x = sin(phi) * cos(theta);
		y = sin(phi) * sin(theta);
		z = cos(phi); // Treat this as local coordinates xyz
		xyz2ae(x, y, z, azel);

		// if in visible area, record the point
		if (z > 0)
		{
			if (antEl == 90)
			{
				sat[numVisPt].x = x;
				sat[numVisPt].y = y;
				sat[numVisPt].z = z;
				sat[numVisPt].az = azel[0];
				sat[numVisPt].el = azel[1];
				sat[numVisPt].snr = 0; // snr will be updated in main()
				numVisPt++;
			} // catch tan(pi/2) situation
			else
			{
				c = tan(deg2rad(antEl));
				if (-y + c * z > 0)
				{
					sat[numVisPt].x = x;
					sat[numVisPt].y = y;
					sat[numVisPt].z = z;
					sat[numVisPt].az = azel[0];
					sat[numVisPt].el = azel[1];
					sat[numVisPt].snr = 0; // snr will be updated in main()
					numVisPt++;
				}
			}
		}
	}
	return numVisPt;
}

int main(void)
{
	// While elevation angle is adjustable in <config.h>, antenna azimuth is simulated at 180 deg by randSat()
	// Boresight vector is (0, -cos(TRUE_EL), sin(TRUE_EL))

	srand(time(NULL)); // set seed for rand()
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(T);

	double spd, snr;
	int numVisPt;
	double snrAdd;

	printf("SIMULATED INPUT FILE || \"SIMUEPOCH#\" \"TIME\"Epoch# \"SAT\" AZ EL SNR SAT#\n"); // print header

	for (int i = 0; i < NUM_EPOCH; i++)
	{ // one simulation per loop
		SimuSat *visSat = (SimuSat *)malloc(NUM_SAT_SPHERE * sizeof(SimuSat));
		if (visSat == NULL)
		{
			fprintf(stderr, "malloc() failed for creating visSat\n");
			exit(-1);
		}
		numVisPt = randSat(visSat, NUM_SAT_SPHERE, TRUE_EL);
		if (numVisPt != 0)
		{
			for (int j = 0; j < numVisPt; j++)
			{
				/* compute SNR by assumed relationship between SNR and off-boresight angle */
				spd = spDist(visSat[j].x, visSat[j].y, visSat[j].z, 0, -cos(deg2rad(TRUE_EL)), sin(deg2rad(TRUE_EL))); // spd: off-boresight angle in rad
				snr = SNR_A * cos(spd) + SNR_C;																		   // cosine relationship is default (preferred)
				visSat[j].snr = snr;

				/* apply uniform SNR variation */
				snrAdd = (double)SNR_STD * gsl_ran_gaussian_ziggurat(r, 1);

				visSat[j].snr += snrAdd;

				/*
					print
				*/
				printf("SIMUEPOCH# TIME%06d SAT %lf %lf %lf %i\n", i, visSat[j].az, visSat[j].el, visSat[j].snr, j); // for output as input file
			}
		}

		/* free after use*/
		free(visSat);
	}
	gsl_rng_free(r);

	exit(0);
} // end of main()
