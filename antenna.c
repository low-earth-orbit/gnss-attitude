#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include "util.h"
#include "truth.h"
#include "struct.h"

/*

Information for users of this program:
1. Modify the below configuration if necessary
2. Input file should not contain satellites without azimuth, elevation or SNR information.

gcc -Wall antenna.c mathutil.c struct.c -o antenna -lgsl -lgslcblas -lm
./antenna > output.txt
valgrind --leak-check=full -s ./antenna
*/

/* Configuration */
#define INPUT_FILE_PATH "input.txt"
#define MAX_NUM_EPOCH NUM_EPOCH

/* Usually no need to change*/
#define MAX_NUM_SAT_EPOCH 100 // maximum number of satellites visible in an epoch
#define MAX_NUM_SIGNAL (MAX_NUM_EPOCH * MAX_NUM_SAT_EPOCH)
#define MAX_NUM_CHAR_LINE 100 // num of char in each record or line
#define NUM_CHAR_DATE 10
#define NUM_CHAR_TIME 10
#define NUM_CHAR_SAT 3
#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))

int main(void)
{
	/*
		create arrays
	*/
	char **timeArray; // A time array storing all unique times
	timeArray = (char **)malloc(MAX_NUM_EPOCH * sizeof(char *));

	Sat **satArray; // A sat array storing all satellite signals
	satArray = (Sat **)malloc(MAX_NUM_SIGNAL * sizeof(Sat *));

	Epoch **epochArray; // A epoch array storing all epochs
	epochArray = (Epoch **)malloc(MAX_NUM_EPOCH * sizeof(Epoch *));

	long int satArrayIndex = 0;
	long int timeArrayIndex = 0;

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
		File processing to get timeArray and satArray
	*/
	char *line = malloc(sizeof(char) * (MAX_NUM_CHAR_LINE));
	fgets(line, MAX_NUM_CHAR_LINE, fp); // Skip header
	while (fgets(line, MAX_NUM_CHAR_LINE, fp) != NULL)
	{
		char *time1 = (char *)malloc(sizeof(char) * (NUM_CHAR_DATE + 1));
		char *time2 = (char *)malloc(sizeof(char) * (NUM_CHAR_TIME + 1));
		char *satName = (char *)malloc(sizeof(char) * (NUM_CHAR_SAT + 1));
		double *az = (double *)malloc(sizeof(double));
		double *el = (double *)malloc(sizeof(double));
		double *snr = (double *)malloc(sizeof(double));

		sscanf(line, "%s %s %s %lf %lf %lf", time1, time2, satName, az, el, snr);
		char *time = concat(time1, time2);

		if (!isStrInArray(time, timeArray, timeArrayIndex))
		{ // Check if the current time is unique
			timeArray[timeArrayIndex] = time;
			timeArrayIndex++;
		}

		satArray[satArrayIndex] = createSat(time, satName, az, el, snr); // Add each sat to sat array
		satArrayIndex++;

		free(time1);
		free(time2); // free memory because time1 and time2 are concatenated to a new char* time
	}
	fclose(fp);
	free(line);

	/* 
		epochArray
		for each unique time with num of sat >= 3
			1) create an Epoch object;
			2) load it to epochArray
	*/
	long int epochArrayIndex = 0; // num of epoch can be less than number of unique time
	for (long int i = 0; i < timeArrayIndex; i++)
	{
		Sat **epochSatArray = (Sat **)malloc(MAX_NUM_SAT_EPOCH * sizeof(Sat *));

		/*
			Count number of sat in the array
		*/
		int *epochSatArrayIndex = (int *)malloc(sizeof(int));
		*epochSatArrayIndex = 0;
		for (long int j = 0; j < satArrayIndex; j++)
		{
			if (strcmp(timeArray[i], satArray[j]->time) == 0)
			{
				epochSatArray[*epochSatArrayIndex] = satArray[j];
				*epochSatArrayIndex = *epochSatArrayIndex + 1;
			}
		}

		if (*epochSatArrayIndex >= 3)
		{
			epochArray[epochArrayIndex] = createEpoch(timeArray[i], epochSatArray, epochSatArrayIndex);
			epochArrayIndex++;
		}
		else
		{ // if an epoch has less than 3 sat, do not record the epoch
			free(epochSatArray);
			free(epochSatArrayIndex);
		}
	}

	/*
		print epoch array to check file input read
	*/
	//printEpochArray(epochArray, epochArrayIndex);

	/*
		Duncan's method (Duncan & Dunn, 1998) -- Vector sum of signal-to-noise ratio (SNR) weighted line-of-sight (LOS) vectors
	*/
	Sol **dunSolArray = malloc(sizeof(Sol *) * epochArrayIndex);
	// print header of the output
	printf("================== Duncan's method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (long int i = 0; i < epochArrayIndex; i++)
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
			//printf("%lf,%lf,%lf\n", xyz[0], xyz[1], xyz[2]);

			/* weight LOS vector by SNR */
			xyz[0] *= *(*epochArray[i]).epochSatArray[j]->snr;
			xyz[1] *= *(*epochArray[i]).epochSatArray[j]->snr;
			xyz[2] *= *(*epochArray[i]).epochSatArray[j]->snr;
			//printf("snr = %lf, %lf,%lf,%lf\n", *(*epochArray[i]).epochSatArray[j]->snr, xyz[0], xyz[1], xyz[2]);

			/* add to vector sum */
			*(dunSol->x) += xyz[0];
			*(dunSol->y) += xyz[1];
			*(dunSol->z) += xyz[2];
			//printf("%lf,%lf,%lf\n", *(dunSol->x), *(dunSol->y), *(dunSol->z));
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
		Geometry method -- vector sum of non-weighted line-of-sight (LOS) vectors with geometry adjustment for elevation angle. Duncan's method is biased toward the spherical area where the satellite signals come from. If the antenna's elevation angle is negative, Duncan's method yields a positive elevation angle estimate. This method is designed by the author of the program to address the geometry issue existing in Duncan's method.
	*/
	Sol **geoSolArray = malloc(sizeof(Sol *) * epochArrayIndex);
	// print header of the output
	printf("================== Geometry method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (long int i = 0; i < epochArrayIndex; i++)
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
			//printf("%lf,%lf,%lf\n", xyz[0], xyz[1], xyz[2]);

			/* add to vector sum */
			*(geoSol->x) += xyz[0];
			*(geoSol->y) += xyz[1];
			*(geoSol->z) += xyz[2];
			//printf("%lf,%lf,%lf\n", *(geoSol->x), *(geoSol->y), *(geoSol->z));
		}

		/* convert the sum to unit vector to get xyz solution*/
		normalize(geoSol);

		/* from xyz solution derive azimuth-elevation solution */
		xyz2aeSol(*(geoSol->x), *(geoSol->y), *(geoSol->z), geoSol);

		/* adjust the elevation angle el = (2 * el - 90) */
		//printf("Before: az = %lf, el = %lf\n", *(geoSol->az), *(geoSol->el));
		*(geoSol->el) = *(geoSol->el) * 2.0 - 90.0;
		//printf("After: az = %lf, el = %lf\n", *(geoSol->az), *(geoSol->el));

		/* recompute xyz solution using adjusted elevation angle */
		//printf("Before: x = %lf, y = %lf, z = %lf\n", *(geoSol->x), *(geoSol->y), *(geoSol->z));
		ae2xyzSol(*(geoSol->az), *(geoSol->el), geoSol);
		//printf("After: x = %lf, y = %lf, z = %lf\n", *(geoSol->x), *(geoSol->y), *(geoSol->z));

		/* print result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", (*epochArray[i]).time, *(*epochArray[i]).numSat, *(geoSol->x), *(geoSol->y), *(geoSol->z), *(geoSol->az), *(geoSol->el));

		/* save to array */
		geoSolArray[i] = geoSol;
	}
	//////////////// above OK
	/*
		Axelrad's method (Axelrad & Behre, 1999) -- Compared to Duncan's method, this is the proper use of SNR in determining antenna boresight vector. It requires antenna gain mapping (the relationship between off-boresight angle and SNR for the antenna) and adjustment to measured SNR.
	*/
	Sol **axelSolArray = malloc(sizeof(Sol *) * epochArrayIndex);

	//double xyz[epochArrayIndex][MAX_NUM_SAT_EPOCH][3]; // line of sight vectors
	//double xyzSol[epochArrayIndex][3];
	//double aeSol[epochArrayIndex][2];

	printf("================== Axelrad's method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");
	for (long int i = 0; i < epochArrayIndex; i++)
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

			/* get sigma of cos(a) to be used in weight matrix W */
			double sigmaSnr = 0.5 + (3 - 0.5) * (*(*epochArray[i]).epochSatArray[j]->snr - 35) / (50 - 35);
			if (sigmaSnr < 0.5)														  // catch the case that sigma <= 0
				sigmaSnr = 0.5;														  // set to minumum sigmaSnr
			double sigma = sigmaSnr / 15.0;											  // uncertainty is a function of SNR. divide sigmaSnr by the coefficient of sigma_cosA
			double cosA = (*(*epochArray[i]).epochSatArray[j]->snr - 35) / (50 - 35); // find cosA from the mapping function snr = (MAX_SNR-MIN_SNR)*cos(A)+MIN_SNR;

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

		double rmsAxel;
		double sumAxel = 0;

		for (long int i = 0; i < epochArrayIndex; i++)
		{
			double trueAntennaXyz[3];
			ae2xyz(TRUE_AZ, TRUE_EL, trueAntennaXyz);
			sumDun += pow(spDist(*dunSolArray[i]->x, *dunSolArray[i]->y, *dunSolArray[i]->z, trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]), 2);
			sumGeo += pow(spDist(*geoSolArray[i]->x, *geoSolArray[i]->y, *geoSolArray[i]->z, trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]), 2);
			sumAxel += pow(spDist(*axelSolArray[i]->x, *axelSolArray[i]->y, *axelSolArray[i]->z, trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]), 2);
		}

		rmsDun = sqrt(sumDun / epochArrayIndex);
		rmsDun = rad2deg(rmsDun);

		rmsGeo = sqrt(sumGeo / epochArrayIndex);
		rmsGeo = rad2deg(rmsGeo);

		rmsAxel = sqrt(sumAxel / epochArrayIndex);
		rmsAxel = rad2deg(rmsAxel);

		printf("================== Statistics ==================\n%li epochs, antenna @ %i deg\nRMS Duncan's method = %lf deg\nRMS Geometry method = %lf deg\nRMS Axelrad's method = %lf deg\n", epochArrayIndex, TRUE_EL, rmsDun, rmsGeo, rmsAxel);
	}

	/*
		free() file input 
	*/
	free(timeArray); // element of this element, time, will be freed in satArray[i]->time

	for (long int i = 0; i < satArrayIndex; i++)
	{
		free(satArray[i]->time);
		free(satArray[i]->satName);
		free(satArray[i]->az);
		free(satArray[i]->el);
		free(satArray[i]->snr); // free attributes
		free(satArray[i]);		// free Sat
	}
	free(satArray); // free satArray

	for (long int i = 0; i < epochArrayIndex; i++)
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
	for (long int i = 0; i < epochArrayIndex; i++)
	{
		free(dunSolArray[i]->x);
		free(dunSolArray[i]->y);
		free(dunSolArray[i]->z);
		free(dunSolArray[i]->az);
		free(dunSolArray[i]->el);
		free(dunSolArray[i]);
	}
	free(dunSolArray);

	for (long int i = 0; i < epochArrayIndex; i++)
	{
		free(geoSolArray[i]->x);
		free(geoSolArray[i]->y);
		free(geoSolArray[i]->z);
		free(geoSolArray[i]->az);
		free(geoSolArray[i]->el);
		free(geoSolArray[i]);
	}
	free(geoSolArray);

	for (long int i = 0; i < epochArrayIndex; i++)
	{
		free(axelSolArray[i]->x);
		free(axelSolArray[i]->y);
		free(axelSolArray[i]->z);
		free(axelSolArray[i]->az);
		free(axelSolArray[i]->el);
		free(axelSolArray[i]);
	}
	free(axelSolArray);
	/*
		exit
	*/
	exit(0);

} // end of main()