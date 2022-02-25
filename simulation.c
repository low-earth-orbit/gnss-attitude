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

a. Only the satellites above Earth's horizon are visible regardless receiving antenna's direction. This is true when the antenna is located on or close to the Earth's surface. The user may set a cutoff angle for post-processing. Assume such angle is 0 here.

b. Antenna can only receive signals from satellites above the antenna's horizon plane, typical for GNSS receiving antennas. This plane is orthogonal to antenna axis.

c. It is assumed that GNSS satellites are randomly distributed on the celestial sphere.

gcc -Wall simulation.c util.c struct.c -o simulation  -lgsl -lgslcblas -lm
./simulation > input.txt
*/

typedef struct SimuSat
{
	double x;
	double y;
	double z; // x(E), y(N), z(U) in local coordinates
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
	double theta, phi;
	double azel[2];
	double xyz[3];

	/* r_e = 6370	r_s = 6370 + 21450 = 27820 */
	double h = 6370 / 27820;
	double c = tan(M_PI_2 - deg2rad(antEl));

	for (int i = 1; i < n; i++)
	{
		theta = 2.0 * M_PI * (double)rand() / (double)RAND_MAX;
		phi = acos(2.0 * (double)rand() / (double)RAND_MAX - 1);
		xyz[0] = sin(phi) * cos(theta);
		xyz[1] = sin(phi) * sin(theta);
		xyz[2] = cos(phi) / 2;

		// if in visible area, record the point
		if (xyz[2] > h)
		{
			if (antEl == 90 || (antEl > 0 && xyz[2] > c * xyz[1] + h) || (antEl < 0 && xyz[2] < c * xyz[1] + h) || (antEl == 0 && xyz[1] < 0))
			{
				xyz[2] = cos(phi) - h; // adjust by observer location
				normalizeXyz(xyz);
				xyz2ae(xyz[0], xyz[1], xyz[2], azel);

				sat[numVisPt].x = xyz[0];
				sat[numVisPt].y = xyz[1];
				sat[numVisPt].z = xyz[2];
				sat[numVisPt].az = azel[0];
				sat[numVisPt].el = azel[1];
				sat[numVisPt].snr = 0; // snr will be updated in main()
				numVisPt++;
			}
		}
	}
	return numVisPt;
}

/*
	generate n visible satellite
*/
int randSat2(SimuSat *sat, int n, double antEl)
{
	int numVisPt = 0;
	double theta, phi;
	double azel[2];
	double xyz[3];

	/* r_e = 6370	r_s = 6370 + 21450 = 27820 */
	double h = 6370 / 27820;
	double c = tan(M_PI_2 - deg2rad(antEl));

	while (numVisPt < n)
	{
		theta = 2.0 * M_PI * (double)rand() / (double)RAND_MAX;
		phi = acos(2.0 * (double)rand() / (double)RAND_MAX - 1);
		xyz[0] = sin(phi) * cos(theta);
		xyz[1] = sin(phi) * sin(theta);
		xyz[2] = cos(phi) / 2;

		// if in visible area, record the point
		if (xyz[2] > h)
		{
			if (antEl == 90 || (antEl > 0 && xyz[2] > c * xyz[1] + h) || (antEl < 0 && xyz[2] < c * xyz[1] + h) || (antEl == 0 && xyz[1] < 0))
			{
				xyz[2] = cos(phi) - h; // adjust by observer location
				normalizeXyz(xyz);
				xyz2ae(xyz[0], xyz[1], xyz[2], azel);

				sat[numVisPt].x = xyz[0];
				sat[numVisPt].y = xyz[1];
				sat[numVisPt].z = xyz[2];
				sat[numVisPt].az = azel[0];
				sat[numVisPt].el = azel[1];
				sat[numVisPt].snr = 0; // snr will be updated in main()
				numVisPt++;
			}
		}
	}
	return numVisPt;
}

int testSat(SimuSat *sat, int n, double antEl)
{
	int numVisPt = 0;
	double theta, phi;
	double xyz[3];

	/* r_e = 6370	r_s = 6370 + 21450 = 27820 */
	double h = 6370.0 / 27820.0;
	double c = tan(M_PI_2 - deg2rad(antEl));

	for (int i = 1; i < n; i++)
	{
		theta = 2.0 * M_PI * (double)rand() / (double)RAND_MAX;
		phi = acos(2.0 * (double)rand() / (double)RAND_MAX - 1);
		xyz[0] = sin(phi) * cos(theta);
		xyz[1] = sin(phi) * sin(theta);
		xyz[2] = cos(phi) / 2;

		// if in visible area, record the point
		if (xyz[2] > h)
		{
			if (antEl == 90 || (antEl > 0 && xyz[2] > c * xyz[1] + h) || (antEl < 0 && xyz[2] < c * xyz[1] + h) || (antEl == 0 && xyz[1] < 0))
			{
				// xyz[2] = cos(phi) - h; // adjust by observer location
				// normalizeXyz(xyz);
				// xyz2ae(xyz[0], xyz[1], xyz[2], azel);

				sat[numVisPt].x = xyz[0];
				sat[numVisPt].y = xyz[1];
				sat[numVisPt].z = xyz[2];
				sat[numVisPt].az = 0;
				sat[numVisPt].el = 0;
				sat[numVisPt].snr = 0;
				numVisPt++;
			}
		}
	}
	return numVisPt;
}

int main(void)
{
	// Elevation angle is adjustable in <config.h>, antenna azimuth is simulated at 180 deg by randSat()

	FILE *fpw = NULL;
	fpw = fopen("input.txt", "w");

	srand(time(NULL)); // set seed for rand()
	gsl_rng_env_setup();
	const gsl_rng_type *T = gsl_rng_default;
	gsl_rng *r = gsl_rng_alloc(T);
	// gsl_rng *r2 = gsl_rng_alloc(T);

	double spd, snr;
	int numVisPt;
	double snrAdd;
	fprintf(fpw, "SIMULATED INPUT FILE\n"); // print header

	for (int i = 0; i < NUM_EPOCH; i++)
	{ // one simulation per loop
		SimuSat *visSat;
		if (SIMULATE_SAT_VISIBLE)
		{
			visSat = (SimuSat *)malloc(NUM_SAT_VISIBLE * sizeof(SimuSat));

			if (visSat == NULL)
			{
				fprintf(stderr, "malloc() failed for creating visSat\n");
				exit(-1);
			}

			numVisPt = randSat2(visSat, NUM_SAT_VISIBLE, TRUE_EL);
		}
		else
		{
			visSat = (SimuSat *)malloc(NUM_SAT_SPHERE * sizeof(SimuSat));

			if (visSat == NULL)
			{
				fprintf(stderr, "malloc() failed for creating visSat\n");
				exit(-1);
			}
			numVisPt = randSat(visSat, NUM_SAT_SPHERE, TRUE_EL);
		}

		double snrSigma;
		if (numVisPt != 0)
		{
			for (int j = 0; j < numVisPt; j++)
			{
				/* compute SNR by assumed relationship between SNR and off-boresight angle */
				spd = spDist(visSat[j].x, visSat[j].y, visSat[j].z, 0, -cos(deg2rad(TRUE_EL)), sin(deg2rad(TRUE_EL))); // spd: off-boresight angle in rad
				// snr = SNR_A * cos(spd) + SNR_C;																		   // cosine relationship
				snr = -(SNR_A / 8100.0) * pow(rad2deg(spd), 2) + SNR_C; // quadratic
				visSat[j].snr = snr;

				/* calculate satellite's elevation angle */
				double alphaS = rad2deg(asin(visSat[j].z));

				snrSigma = SNR_STD_MAX + ((SNR_STD_MIN - SNR_STD_MAX) / 90.0) * alphaS;

				snrAdd = gsl_ran_gaussian_ziggurat(r, snrSigma);
				visSat[j].snr += snrAdd;

				/*
					print
				*/
				fprintf(fpw, "SIMUEPOCH# TIME%06d SAT %lf %lf %lf %i\n", i, visSat[j].az, visSat[j].el, visSat[j].snr, j); // for output as input file
			}
		}
		free(visSat);
	}
	fclose(fpw);

	// x, y, z for visualization, 1 epoch
	FILE *fpwt = NULL;
	fpwt = fopen("data.txt", "w");
	int numTestSat = 10000;
	SimuSat *testVisSat = (SimuSat *)malloc(numTestSat * sizeof(SimuSat));
	numVisPt = testSat(testVisSat, numTestSat, TRUE_EL);
	if (numVisPt != 0)
	{
		for (int j = 0; j < numVisPt; j++)
		{
			fprintf(fpwt, "%lf %lf %lf\n", testVisSat[j].x, testVisSat[j].y, testVisSat[j].z);
		}
	}
	free(testVisSat);
	fclose(fpwt);
	// end of test

	gsl_rng_free(r);

	exit(0);
} // end of main()
