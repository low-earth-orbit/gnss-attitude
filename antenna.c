#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include "util.h"	// math utilities
#include "config.h" // configuration
#include "struct.h" // structures
#include "snr.h"	// snr adjustment & mapping

/*

Information for users of this program:
1. Modify the below configuration if necessary
2. Input file should not contain satellites without azimuth, elevation or SNR information.

sudo apt-get install libgsl-dev
gcc -Wall antenna.c util.c struct.c snr.c -o antenna -lgsl -lgslcblas -lm
./antenna > output.txt
valgrind --tool=memcheck --leak-check=yes --leak-check=full -s --track-origins=yes --show-leak-kinds=all ./antenna
*/

int main(void)
{
	/*
		open file for reading
	*/
	FILE *fp = NULL;
	fp = fopen(INPUT_FILE_PATH, "r");
	if (fp == NULL)
	{
		perror("Error opening file\n");
		return (-1);
	}

	/*
		record satArray
	*/
	Sat **satArray; // A sat array storing all satellite signals
	satArray = (Sat **)malloc(MAX_NUM_SIGNAL * sizeof(Sat *));
	long int satArrayIndex = 0;

	char *line = malloc(sizeof(char) * (MAX_NUM_CHAR_LINE));
	fgets(line, MAX_NUM_CHAR_LINE, fp); // Skip header
	while (fgets(line, MAX_NUM_CHAR_LINE, fp) != NULL)
	{
		char *time1 = (char *)malloc(sizeof(char) * (NUM_CHAR_DATE + 1));
		char *time2 = (char *)malloc(sizeof(char) * (NUM_CHAR_TIME + 1));
		char *prn = (char *)malloc(sizeof(char) * (NUM_CHAR_SAT + 1));
		double *az = (double *)malloc(sizeof(double));
		double *el = (double *)malloc(sizeof(double));
		double *snr = (double *)malloc(sizeof(double));

		sscanf(line, "%s %s %s %lf %lf %lf", time1, time2, prn, az, el, snr);
		char *time = concat(time1, time2);
		/* adjust SNR here */
		if (!SIMULATION)
		{
			printf("SNR (meas) = %lf\t", *snr);
			adjSnr(prn, el, snr);
			printf("SNR (adj) = %lf\n", *snr);
		}
		satArray[satArrayIndex] = createSat(time, prn, az, el, snr); // Add each sat to sat array
		satArrayIndex++;

		free(time1);
		free(time2); // free memory because time1 and time2 are concatenated to a new char* time
	}

	/*
		close file pointer
	*/
	free(line);
	fclose(fp);

	/*
		catch empty data error
	*/
	if (satArrayIndex == 0)
	{
		perror("No received signals found in file\n");
		return (-1);
	}

	if (satArrayIndex < 3)
	{
		perror("Less than 3 signals found in file\n");
		return (-1);
	}

	/*
		sort satArray by time
	*/
	qsort(satArray, satArrayIndex, sizeof(Sat *), cmpSatArray);
	/* 	for (int j = 0; j < satArrayIndex; j++)
	{
		printf("After sort: %s\t%s\n", satArray[j]->time, satArray[j]->prn);
	} */

	/* 
		form epochArray
	*/
	Epoch **epochArray; // A epoch array storing all epochs
	epochArray = (Epoch **)malloc(MAX_NUM_EPOCH * sizeof(Epoch *));

	long int *epochArrayIndex = (long int *)malloc(sizeof(long int));
	*epochArrayIndex = 0;

	long int *i = malloc(sizeof(long int)); // satArray counter
	*i = 0;

	while (true)
	{
		Sat **epochSatArray = (Sat **)malloc(MAX_NUM_SAT_EPOCH * sizeof(Sat *));
		int *epochSatArrayIndex = (int *)malloc(sizeof(int));
		*epochSatArrayIndex = 0;

		while (*i < satArrayIndex)
		{

			if (*epochSatArrayIndex == 0 && *i != satArrayIndex - 1) // first sat in epoch && not the last sat
			{
				epochSatArray[*epochSatArrayIndex] = satArray[*i];
				*epochSatArrayIndex += 1;
			}

			else if (*epochSatArrayIndex == 0 && *i == satArrayIndex - 1) // first sat in epoch && the last sat
			{
				free(epochSatArrayIndex);
				free(epochSatArray);
			}
			else if (strcmp(satArray[*i]->time, satArray[*i - 1]->time) == 0) // current one belongs to the same epoch
			{
				epochSatArray[*epochSatArrayIndex] = satArray[*i];
				*epochSatArrayIndex += 1;

				if (*i == satArrayIndex - 1)
				{
					if (*epochSatArrayIndex >= 3)
					{
						epochArray[*epochArrayIndex] = createEpoch(satArray[*i - 1]->time, epochSatArray, epochSatArrayIndex);
						*epochArrayIndex += 1;
					}
					else // do not record to array
					{
						free(epochSatArrayIndex);
						free(epochSatArray);
					}
				}
			}
			else if (strcmp(satArray[*i]->time, satArray[*i - 1]->time) != 0) // current one belongs to a new epoch
			{
				if (*epochSatArrayIndex >= 3)
				{
					epochArray[*epochArrayIndex] = createEpoch(satArray[*i - 1]->time, epochSatArray, epochSatArrayIndex);
					*epochArrayIndex += 1;
				}
				else // do not record to array
				{
					free(epochSatArrayIndex);
					free(epochSatArray);
				}
				if (*i != satArrayIndex - 1)
				{
					break; // if not the last sat: break the inner while loop; skip increment *i
				}
			}
			*i += 1;
		}

		if (*i == satArrayIndex)
		{
			break; // break the outer while loop
		}
	}
	free(i);

	/*
		print epoch array to check file input read
	*/
	printEpochArray(epochArray, *epochArrayIndex);

	/*
		Duncan's method (Duncan & Dunn, 1998) -- Vector sum of signal-to-noise ratio (SNR) weighted line-of-sight (LOS) vectors. Duncan's method is biased toward the spherical area where the satellite signals come from. If the antenna's elevation angle is negative, Duncan's method yields a positive elevation angle estimate.
	*/
	Sol **dunSolArray = malloc(sizeof(Sol *) * *epochArrayIndex);
	// print header of the output
	printf("================== Duncan's method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		int n = *(epochArray[i]->numSat);
		Sol *dunSol = malloc(sizeof(Sol));
		dunSol->x = calloc(1, sizeof(double));
		dunSol->y = calloc(1, sizeof(double));
		dunSol->z = calloc(1, sizeof(double));
		dunSol->az = calloc(1, sizeof(double));
		dunSol->el = calloc(1, sizeof(double));

		for (int j = 0; j < n; j++)
		{
			/* calculate LOS vector from azimuth and elevation*/
			double xyz[3];
			ae2xyz(*(*epochArray[i]).epochSatArray[j]->az, *(*epochArray[i]).epochSatArray[j]->el, xyz);

			/* weight LOS vector by SNR */
			xyz[0] *= *(*epochArray[i]).epochSatArray[j]->snr;
			xyz[1] *= *(*epochArray[i]).epochSatArray[j]->snr;
			xyz[2] *= *(*epochArray[i]).epochSatArray[j]->snr;

			/* add to vector sum */
			*(dunSol->x) += xyz[0];
			*(dunSol->y) += xyz[1];
			*(dunSol->z) += xyz[2];
		}

		/* convert the sum to unit vector to get xyz solution*/
		normalize(dunSol);

		/* from xyz solution derive azimuth-elevation solution */
		xyz2aeSol(*(dunSol->x), *(dunSol->y), *(dunSol->z), dunSol);

		/* print result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", (*epochArray[i]).time, *(*epochArray[i]).numSat, *(dunSol->x), *(dunSol->y), *(dunSol->z), *(dunSol->az), *(dunSol->el));

		dunSolArray[i] = dunSol;
	}

	/*
		LOS with geometric correction made to elevation angle -- vector sum of non-weighted line-of-sight (LOS) vectors with geometry adjustment for elevation angle. This method is designed by the author of the program to address the geometry issue existing in Duncan's method.
	*/
	Sol **geoSolArray = malloc(sizeof(Sol *) * *epochArrayIndex);
	// print header of the output
	printf("================== LOS (Geometry) method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		int n = *(epochArray[i]->numSat);
		Sol *geoSol = malloc(sizeof(Sol));
		geoSol->x = calloc(1, sizeof(double));
		geoSol->y = calloc(1, sizeof(double));
		geoSol->z = calloc(1, sizeof(double));
		geoSol->az = calloc(1, sizeof(double));
		geoSol->el = calloc(1, sizeof(double));

		for (int j = 0; j < n; j++)
		{
			/* calculate LOS vector from azimuth and elevation*/
			double xyz[3];
			ae2xyz(*(*epochArray[i]).epochSatArray[j]->az, *(*epochArray[i]).epochSatArray[j]->el, xyz);

			/* add to vector sum */
			*(geoSol->x) += xyz[0];
			*(geoSol->y) += xyz[1];
			*(geoSol->z) += xyz[2];
		}

		/* convert the sum to unit vector to get xyz solution*/
		normalize(geoSol);

		/* from xyz solution derive azimuth-elevation solution */
		xyz2aeSol(*(geoSol->x), *(geoSol->y), *(geoSol->z), geoSol);

		/* adjust the elevation angle el = (2 * el - 90) */
		*(geoSol->el) = *(geoSol->el) * 2.0 - 90.0;

		/* recompute xyz solution using adjusted elevation angle */
		ae2xyzSol(*(geoSol->az), *(geoSol->el), geoSol);

		/* print result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", (*epochArray[i]).time, *(*epochArray[i]).numSat, *(geoSol->x), *(geoSol->y), *(geoSol->z), *(geoSol->az), *(geoSol->el));

		/* save to array */
		geoSolArray[i] = geoSol;
	}

	/*
		LOS with elevation angle determined by azimuth dispersion -- vector sum of non-weighted line-of-sight (LOS) vectors for azimuth. elevaltion angle is from the correlation between {circular standard deviation of azimuth} and {elevation angle}. The relationship should be built through simulation or calibration data set. This method is designed by the author of the program to address the geometry issue existing in Duncan's method.
	*/
	Sol **statSolArray = malloc(sizeof(Sol *) * *epochArrayIndex);
	// print header of the output
	printf("================== LOS (Statistics) method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		int n = *(epochArray[i]->numSat);
		Sol *statSol = malloc(sizeof(Sol));
		statSol->x = calloc(1, sizeof(double));
		statSol->y = calloc(1, sizeof(double));
		statSol->z = calloc(1, sizeof(double));
		statSol->az = calloc(1, sizeof(double));
		statSol->el = calloc(1, sizeof(double));

		for (int j = 0; j < n; j++)
		{
			/* calculate LOS vector from azimuth and elevation*/
			double xyz[3];
			ae2xyz(*(*epochArray[i]).epochSatArray[j]->az, *(*epochArray[i]).epochSatArray[j]->el, xyz);

			/* add to vector sum */
			*(statSol->x) += xyz[0];
			*(statSol->y) += xyz[1];
			*(statSol->z) += xyz[2];
		}

		/* convert the sum to unit vector to get xyz solution*/
		normalize(statSol);

		/* from xyz solution derive azimuth-elevation solution */
		xyz2aeSol(*(statSol->x), *(statSol->y), *(statSol->z), statSol);

		/* elevation angle */
		double std = cirStdAzEpoch(epochArray[i]);
		double elFunc = std * 98.323 - 262.77; // <- this relationship is built through simulation or calibration data set
		if (elFunc > 90)
		{ // catch overflow
			elFunc = 90;
		}
		else if (elFunc <= -90)
		{
			elFunc = -90;
		}
		*(statSol->el) = elFunc;
		//double std = spStdEpoch(epochArray[i]);
		//*(statSol->el) = 383.71 * std * std - 26.155 * std - 167.54;

		/* recompute xyz solution using new elevation angle */
		ae2xyzSol(*(statSol->az), *(statSol->el), statSol);

		/* print result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", (*epochArray[i]).time, *(*epochArray[i]).numSat, *(statSol->x), *(statSol->y), *(statSol->z), *(statSol->az), *(statSol->el));

		/* save to array */
		statSolArray[i] = statSol;
	}

	/*
		Axelrad's method (Axelrad & Behre, 1999) -- Compared to Duncan's method, this is the proper use of SNR in determining antenna boresight vector. It requires antenna gain mapping (the relationship between off-boresight angle and SNR for the antenna) and adjustment to measured SNR.
	*/
	Sol **axelSolArray = malloc(sizeof(Sol *) * *epochArrayIndex);

	printf("================== Axelrad's method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");
	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		/*	
			Observation equation X*b = cos(a) corresponds to X*c = y below
			See GNU Scientific Library Reference Manual for more: https://www.gnu.org/software/gsl/doc/html/lls.html
		*/
		int n = *(epochArray[i]->numSat); // number of observations in the epoch
		Sol *axelSol = malloc(sizeof(Sol));
		axelSol->x = calloc(1, sizeof(double));
		axelSol->y = calloc(1, sizeof(double));
		axelSol->z = calloc(1, sizeof(double));
		axelSol->az = calloc(1, sizeof(double));
		axelSol->el = calloc(1, sizeof(double));

		/* set up size of matrices for least squares (LS) regression */
		double chisq;
		gsl_matrix *X, *cov;
		gsl_vector *y, *w, *c;
		X = gsl_matrix_alloc(n, 3);
		y = gsl_vector_alloc(n);	  // n*1 matrix
		w = gsl_vector_alloc(n);	  // n diagonal elements of n*n weight matrix
		c = gsl_vector_alloc(3);	  // coefficients (x, y, z) -> the boresight vector
		cov = gsl_matrix_alloc(3, 3); // cov = inverse(transpose(X) W X)

		for (int j = 0; j < n; j++)
		{
			/* calculate LOS vector from azimuth and elevation*/
			double xyz[3];
			ae2xyz(*(*epochArray[i]).epochSatArray[j]->az, *(*epochArray[i]).epochSatArray[j]->el, xyz);

			double cosA;
			if (SIMULATION)
			{
				cosA = (*(*epochArray[i]).epochSatArray[j]->snr - SNR_C) / SNR_A; // find cosA from the mapping function snr = (MAX_SNR-MIN_SNR)*cos(A)+MIN_SNR;
			}
			else // real data, not simulation
			{
				double alpha = getAlpha((epochArray[i])->epochSatArray[j]->prn, (*epochArray[i]).epochSatArray[j]->snr);
				/*
				double alphaSq = (*(*epochArray[i]).epochSatArray[j]->snr - MAP_B) / MAP_A;
				if (alphaSq < 0)
				{
					alphaSq = 0;
				}
				*/
				cosA = cos(deg2rad(alpha));
			}

			if (cosA > 1)
			{ // catch the case that cosA > 1
				cosA = 1;
			}
			else if (cosA < 0)
			{ // catch the case that cosA < 0
				cosA = 0;
			}

			double sigmaSnr = ((SNR_STD_MAX - SNR_STD_MIN) / 8100.0) * pow(rad2deg(acos(cosA)), 2) + SNR_STD_MIN;
			//printf("alpha = %lf snr_std = %lf\n", rad2deg(acos(cosA)), sigmaSnr);
			double sigma = sigmaSnr / SNR_A;

			// Set each observation equation
			gsl_matrix_set(X, j, 0, xyz[0]); // coefficient c0 = x
			gsl_matrix_set(X, j, 1, xyz[1]); // coefficient c1 = y
			gsl_matrix_set(X, j, 2, xyz[2]); // coefficient c2 = z
			gsl_vector_set(y, j, cosA);
			gsl_vector_set(w, j, 1.0 / (sigma * sigma));
		}

		/* run multi-parameter regression */
		gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, 3);
		gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
		gsl_multifit_linear_free(work);
		/* clang-format off */

		#define C(i) (gsl_vector_get(c, (i)))
		#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))
		/* clang-format on */

		/* save best fit */
		*(axelSol->x) = C(0);
		*(axelSol->y) = C(1);
		*(axelSol->z) = C(2);

		/* free matrices for LS */
		gsl_matrix_free(X);
		gsl_vector_free(y);
		gsl_vector_free(w);
		gsl_vector_free(c);
		gsl_matrix_free(cov);

		/* normalize the resulting vector to get xyz solution*/
		normalize(axelSol);

		/* from xyz solution derive azimuth-elevation solution */
		xyz2aeSol(*(axelSol->x), *(axelSol->y), *(axelSol->z), axelSol);

		/* print result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", (*epochArray[i]).time, *(*epochArray[i]).numSat, *(axelSol->x), *(axelSol->y), *(axelSol->z), *(axelSol->az), *(axelSol->el));

		/* save to array */
		axelSolArray[i] = axelSol;
	}

	/*
		Statistics if true antenna boresight (E, N, U) is known 
	*/
	if (TRUE_EL >= -90 && TRUE_EL <= 90 && TRUE_AZ <= 360 && TRUE_AZ >= 0)
	{
		double rmsDun;
		double sumDun = 0;

		double rmsGeo;
		double sumGeo = 0;

		double rmsStat;
		double sumStat = 0;

		double rmsAxel;
		double sumAxel = 0;

		for (long int i = 0; i < *epochArrayIndex; i++)
		{
			double trueAntennaXyz[3];
			ae2xyz(TRUE_AZ, TRUE_EL, trueAntennaXyz);
			sumDun += pow(spDist(*dunSolArray[i]->x, *dunSolArray[i]->y, *dunSolArray[i]->z, trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]), 2);
			sumGeo += pow(spDist(*geoSolArray[i]->x, *geoSolArray[i]->y, *geoSolArray[i]->z, trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]), 2);
			sumStat += pow(spDist(*statSolArray[i]->x, *statSolArray[i]->y, *statSolArray[i]->z, trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]), 2);
			sumAxel += pow(spDist(*axelSolArray[i]->x, *axelSolArray[i]->y, *axelSolArray[i]->z, trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]), 2);
		}

		rmsDun = sqrt(sumDun / *epochArrayIndex);
		rmsDun = rad2deg(rmsDun);

		rmsGeo = sqrt(sumGeo / *epochArrayIndex);
		rmsGeo = rad2deg(rmsGeo);

		rmsStat = sqrt(sumStat / *epochArrayIndex);
		rmsStat = rad2deg(rmsStat);

		rmsAxel = sqrt(sumAxel / *epochArrayIndex);
		rmsAxel = rad2deg(rmsAxel);

		printf("================== Statistics ==================\n%li epochs, antenna @ %i deg\nRMS Duncan's = %lf deg\nRMS LOS (Geometry) = %lf deg\nRMS LOS (Statistics) = %lf deg\nRMS Axelrad's = %lf deg\n", *epochArrayIndex, (int)TRUE_EL, rmsDun, rmsGeo, rmsStat, rmsAxel);
	}

	/*
		free() file input 
	*/
	for (long int i = 0; i < satArrayIndex; i++)
	{
		free(satArray[i]->time);
		free(satArray[i]->prn);
		free(satArray[i]->az);
		free(satArray[i]->el);
		free(satArray[i]->snr); // free attributes
		free(satArray[i]);		// free Sat
	}
	free(satArray); // free satArray

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		// time freed in satArray[i]
		free(epochArray[i]->epochSatArray); // epochSatArray[i] freed in satArray[i]
		free(epochArray[i]->numSat);
		free(epochArray[i]); // free Epoch
	}
	free(epochArray);

	/*
		free() solutions
	*/
	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		free(dunSolArray[i]->x);
		free(dunSolArray[i]->y);
		free(dunSolArray[i]->z);
		free(dunSolArray[i]->az);
		free(dunSolArray[i]->el);
		free(dunSolArray[i]);
	}
	free(dunSolArray);

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		free(geoSolArray[i]->x);
		free(geoSolArray[i]->y);
		free(geoSolArray[i]->z);
		free(geoSolArray[i]->az);
		free(geoSolArray[i]->el);
		free(geoSolArray[i]);
	}
	free(geoSolArray);

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		free(statSolArray[i]->x);
		free(statSolArray[i]->y);
		free(statSolArray[i]->z);
		free(statSolArray[i]->az);
		free(statSolArray[i]->el);
		free(statSolArray[i]);
	}
	free(statSolArray);

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		free(axelSolArray[i]->x);
		free(axelSolArray[i]->y);
		free(axelSolArray[i]->z);
		free(axelSolArray[i]->az);
		free(axelSolArray[i]->el);
		free(axelSolArray[i]);
	}
	free(axelSolArray);
	free(epochArrayIndex);
	/*
		exit
	*/
	exit(0);

} // end of main()