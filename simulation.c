#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "mathutil.h"
#include "truth.h"

/*
This program simulates the satellite-antenna geometry.

a. Only the satellites above Earth's horizon can be possibly seen regardless of the direction in which the antenna is pointing. This is true when the antenna is located on or close to the Earth's surface. The user may set a cutoff angle for post-processing. Assume such angle is 0 here.

b. Antenna can only receive signals from satellites above the antenna's horizon. This is typical for GNSS antennas receiving signals.

Therefore, there are two hemispheres. When the antenna is pointing up at a 90 deg elevation angle, the two coincide. However, when the antenna is not at a 90 deg elevation angle only the area these hemispheres overlap has visible GNSS satellites.

The assumption of this simulation is that GNSS satellites are randomly distributed on the celestial sphere.

gcc -Wall simulation.c mathutil.c -lm -o simulation
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
	
		points[][5]: points array
		
		antEl: antenna's elevation angle (a.k.a. boresight angle) in degrees; if the antenna is pointing straight up, that is 90 deg.
		
		n: number of points on the ENTIRE sphere
		
		returns: number of satellites visible to the antenna with specific elevation angle
	
*/
int rpVisSat(SimuSat *sat, int n, double antEl)
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

/*
	Return a normally distributed random number 
*/
double randNormal()
{
	double x, y, rsq, f;
	do
	{
		x = 2.0 * (double)rand() / (double)RAND_MAX - 1.0;
		y = 2.0 * (double)rand() / (double)RAND_MAX - 1.0;
		rsq = x * x + y * y;
	} while (rsq >= 1.0 || rsq == 0.0);
	f = sqrt(-2.0 * log(rsq) / rsq);
	return x * f;
}

int main(void)
{
	int antEl = TRUE_EL; // Antenna boresight elevation angle
	// While elevation angle is adjustable, antenna azimuth is simulated at 180 deg by rpVisSat()
	// Boresight vector is (0, -cos(antEl), sin(antEl))
	int numEpoch = NUM_EPOCH; // number of simulated epoch
	int numSat = 100;	  // number of GNSS satellites globally available
	const double MAX_SNR = 50;
	const double MIN_SNR = 35; // Set max and min snr values for snr computation
	const double MAX_SNR_STD = 3;
	const double MIN_SNR_STD = 0.5; // set max and min snr value variation (standard deviation)

	srand(time(NULL)); // set seed for rand()

	double spd, snr;
	int numVisPt;
	double snrAdd;
	printf("SIMULATED INPUT FILE || \"SIMUEPOCH#\" \"TIME\"Epoch# \"SAT\" AZ EL SNR SAT#\n"); // print header
	for (int i = 0; i < numEpoch; i++)
	{ // one simulation per loop
		SimuSat *visSat = (SimuSat *)malloc(numSat * sizeof(SimuSat));
		if (visSat == NULL)
		{
			fprintf(stderr, "malloc() failed for creating visSat\n");
			exit(-1);
		}
		numVisPt = rpVisSat(visSat, numSat, antEl);
		if (numVisPt != 0)
		{
			for (int j = 0; j < numVisPt; j++)
			{
				/* compute SNR by assumed relationship between SNR and off-boresight angle */
				spd = spDist(visSat[j].x, visSat[j].y, visSat[j].z, 0, -cos(deg2rad(antEl)), sin(deg2rad(antEl))); // this spd is off-boresight angle in rad
				snr = (MAX_SNR - MIN_SNR) * cos(spd) + MIN_SNR;													   // cosine relationship is default (preferred), other options below
				//snr = ((MIN_SNR - MAX_SNR)/8100)*pow(rad2deg(spd),2) + MAX_SNR; // quadratic
				//snr = MAX_SNR - ( (M_PI*0.5 - spd) / (0.5*M_PI))*(MAX_SNR - MIN_SNR); // linear
				//printf("%lf %lf\n", rad2deg(spd), snr);
				visSat[j].snr = snr;

				/* apply SNR variation */
				//printf("%lf\n", randNormal());
				snrAdd = randNormal() * (MIN_SNR_STD + (spd / (0.5 * M_PI)) * (MAX_SNR_STD - MIN_SNR_STD)); // assume linear relationship between variation of SNR and off-boresight angle
				visSat[j].snr += snrAdd;
				//printf("%lf, %lf\n", spd, visSat[j].snr);
				/*
					print
				*/
				//printf("xyz = %lf, %lf, %lf || azel = %lf, %lf || snr = %lf\n", visSat[j].x, visSat[j].y, visSat[j].z, visSat[j].az, visSat[j].el, visSat[j].snr); // for check
				//printf("%lf %lf %lf\n", visSat[j].x, visSat[j].y, visSat[j].z);// for plot
				//printf("%lf %lf\n", visSat[j].az, visSat[j].el);// az, el for plot
				printf("SIMUEPOCH# TIME%06d SAT %lf %lf %lf %i\n", i, visSat[j].az, visSat[j].el, visSat[j].snr, j); // for output as input file
			}
		}

		/* free after use*/
		free(visSat);
	}
	/*exit*/
	exit(0);
} // end of main()