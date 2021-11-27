#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_multifit.h>
#include "util.h"	// math utilities
#include "config.h" // configuration
#include "struct.h" // structures
#include "snr.h"	// snr adjustment & mapping

/*
sudo apt-get install libgsl-dev
gcc -Wall antenna.c util.c struct.c snr.c -o antenna -lgsl -lgslcblas -lm
valgrind --tool=memcheck --leak-check=yes --leak-check=full -s --track-origins=yes --show-leak-kinds=all ./antenna
*/
#define C(i) (gsl_vector_get(c, (i)))
#define COV(i, j) (gsl_matrix_get(cov, (i), (j)))

int main(int argc, char **argv)
{
	clock_t startTime = clock();

	/*
		File I/O
	*/
	FILE *fp = NULL;
	FILE *fp2 = NULL;
	FILE *fpw = NULL;

	if (argc == 1)
	{
		fprintf(stderr, "No input file specified. Using default input file \"%s\"\n", INPUT_FILE_PATH_DEFAULT);
		fp = fopen(INPUT_FILE_PATH_DEFAULT, "r");
		if (fp == NULL)
		{
			fprintf(stderr, "Error opening input file \"%s\"\n", argv[1]);
			return (-1);
		}
	}
	else if (argc == 2 || argc == 3)
	{
		fp = fopen(argv[1], "r");
		if (fp == NULL)
		{
			fprintf(stderr, "Error opening input file \"%s\".\n", argv[1]);
			return (-1);
		}

		if (argc == 3)
		{
			fp2 = fopen(argv[2], "r");
			if (fp2 == NULL)
			{
				fprintf(stderr, "Error opening input file \"%s\".\n", argv[2]);
				return (-1);
			}
		}
	}
	else
	{
		fprintf(stderr, "Too many arguments. Usage:\n./antenna input.txt (input2.txt)\n");
		return (-1);
	}

	fpw = fopen(OUTPUT_FILE_PATH_DEFAULT, "w");
	if (fpw == NULL)
	{
		fprintf(stderr, "Error writing output file \"output.txt\".\n");
		return (-1);
	}

	/*
		record satArray
	*/
	Sat **satArray; // A sat array storing all satellite signals
	satArray = (Sat **)malloc(MAX_NUM_SIGNAL * sizeof(Sat *));
	long int satArrayIndex = 0;
	char *line = malloc(sizeof(char) * (MAX_NUM_CHAR_LINE));

	/* L1 input file */
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

		if (!SIMULATION)
		{
			adjSnr(prn, el, snr);
		}
		double *snr2 = (double *)malloc(sizeof(double));
		*snr2 = -1.0;

		satArray[satArrayIndex] = createSat(time, prn, az, el, snr, snr2); // Add each sat to sat array
		satArrayIndex++;

		free(time1);
		free(time2); // free memory because time1 and time2 are concatenated to a new char* time
	}

	/* L2 input file */
	if (argc == 3)
	{
		fgets(line, MAX_NUM_CHAR_LINE, fp2); // Skip header
		while (fgets(line, MAX_NUM_CHAR_LINE, fp2) != NULL)
		{
			char *time1 = (char *)malloc(sizeof(char) * (NUM_CHAR_DATE + 1));
			char *time2 = (char *)malloc(sizeof(char) * (NUM_CHAR_TIME + 1));
			char *prn = (char *)malloc(sizeof(char) * (NUM_CHAR_SAT + 1));
			double *az = (double *)malloc(sizeof(double));
			double *el = (double *)malloc(sizeof(double));

			double *snr2 = (double *)malloc(sizeof(double));

			sscanf(line, "%s %s %s %lf %lf %lf", time1, time2, prn, az, el, snr2);
			char *time = concat(time1, time2);
			/* adjust SNR here */
			if (!SIMULATION)
			{
				//printf("SNR (meas) = %.2f\t", *snr2);
				adjSnr2(prn, el, snr2);
				//printf("SNR (adj) = %.2f\n", *snr2);
			}
			double *snr = (double *)malloc(sizeof(double));
			*snr = -1.0;
			//printf("snr = %.2f\n", *snr);
			satArray[satArrayIndex] = createSat(time, prn, az, el, snr, snr2); // Add each sat to sat array
			satArrayIndex++;

			free(time1);
			free(time2); // free memory because time1 and time2 are concatenated to a new char* time
		}
	}

	/*
		close input file pointer
	*/
	free(line);
	fclose(fp);
	if (argc == 3)
		fclose(fp2);

	/*
		catch empty data error
	*/
	if (satArrayIndex == 0)
	{
		fprintf(stderr, "No received signals found in file\n");
		return (-1);
	}

	if (satArrayIndex < 3)
	{
		fprintf(stderr, "Less than 3 signals found in file\n");
		return (-1);
	}

	/*
		sort satArray by time
	*/
	qsort(satArray, satArrayIndex, sizeof(Sat *), cmpSatArray);

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
		Sat **epochSatArray = (Sat **)malloc(MAX_NUM_SIGNAL_EPOCH * sizeof(Sat *));
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
	//printEpochArray(epochArray, *epochArrayIndex);

	/*
		Axelrad's method (Axelrad & Behre, 1999) -- Compared to Duncan's method, this is the proper use of SNR in determining antenna boresight vector. It requires antenna gain mapping (the relationship between off-boresight angle and SNR for the antenna) and adjustment to measured SNR.
	*/
	Sol **axelSolArray = malloc(sizeof(Sol *) * *epochArrayIndex);

	fprintf(fpw, "Epoch (GPST), # of Signal, X(E), Y(N), Z(U), Az(deg), El(deg)\n");

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
				if (cosA > 1)
				{ // catch the case that cosA > 1
					cosA = 1;
				}
				else if (cosA < 0)
				{ // catch the case that cosA < 0
					cosA = 0;
				}
			}
			else // real data, not simulation
			{
				if (*(*epochArray[i]).epochSatArray[j]->snr < 0)
				{
					cosA = getCosA2((epochArray[i])->epochSatArray[j]->prn, (*epochArray[i]).epochSatArray[j]->snr2);
				}
				else
				{
					cosA = getCosA((epochArray[i])->epochSatArray[j]->prn, (*epochArray[i]).epochSatArray[j]->snr);
				}
			}

			//double sigma = (3.0 / 81000.0) * pow(*(*epochArray[i]).epochSatArray[j]->el, 2) + 0.7;
			double sigma = 1;
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

		/* least squares stats */
		//printf("# covariance matrix:\n");
		//printf("[ %+.5e, %+.5e, %+.5e  \n", COV(0, 0), COV(0, 1), COV(0, 2));
		//printf("  %+.5e, %+.5e, %+.5e  \n", COV(1, 0), COV(1, 1), COV(1, 2));
		//printf("  %+.5e, %+.5e, %+.5e ]\n", COV(2, 0), COV(2, 1), COV(2, 2));
		//printf("# chisq = %g\n", chisq);
		/*
		double redChiSq = chisq / (n - 3); // reduced chisq = chisq / (# of signals - 3)
		double lsStdX = rad2deg(asin(sqrt(redChiSq * COV(0, 0))));
		double lsStdY = rad2deg(asin(sqrt(redChiSq * COV(1, 1))));
		double lsStdZ = rad2deg(asin(sqrt(redChiSq * COV(2, 2))));
		double lsStdA = rad2deg(asin(sqrt(redChiSq * COV(0, 0) + redChiSq * COV(1, 1) + redChiSq * COV(2, 2))));
		printf("standard deviation of x, y, z from least squares\nx = %.2f°, y = %.2f°, z = %.2f°, 3D = %.2f°\n", lsStdX, lsStdY, lsStdZ, lsStdA);
		*/

		/* free matrices for LS */
		gsl_matrix_free(X);
		gsl_vector_free(y);
		gsl_vector_free(w);
		gsl_vector_free(c);
		gsl_matrix_free(cov);

		/* normalize the resulting vector to get xyz solution*/
		normalize(axelSol);

		/* from xyz solution derive azimuth-elevation solution*/
		xyz2aeSol(*(axelSol->x), *(axelSol->y), *(axelSol->z), axelSol);

		/* apply convergence correction built by simulation */
		if (CONVERGENCE_CORRECTION)
		{
			if (*(axelSol->el) > 40.181)
			{
				*(axelSol->el) += 0.0000022678 * pow(*(axelSol->el), 3) - 0.0005537871 * pow(*(axelSol->el), 2) + 0.0455341172 * *(axelSol->el) - 1.2636551736;
			}
			else if (*(axelSol->el) > 0.562466)
			{
				*(axelSol->el) += 0.0000029756 * pow(*(axelSol->el), 3) - 0.0002119836 * pow(*(axelSol->el), 2) + 0.0133919561 * *(axelSol->el) - 0.5684236546;
			}
			else if (*(axelSol->el) > -33.9139)
			{
				*(axelSol->el) += -0.0000374714 * pow(*(axelSol->el), 3) - 0.0058097797 * pow(*(axelSol->el), 2) + 0.0088882136 * *(axelSol->el) - 0.5690671978;
			}
			else if (*(axelSol->el) > -61.2864)
			{
				*(axelSol->el) += 0.0000074641 * pow(*(axelSol->el), 3) + 0.0074059911 * pow(*(axelSol->el), 2) + 0.7488207967 * *(axelSol->el) + 11.0845716711;
			}
			else if (*(axelSol->el) > -73.799)
			{
				*(axelSol->el) += 0.0029469680 * pow(*(axelSol->el), 3) + 0.5621635563 * pow(*(axelSol->el), 2) + 35.6911758009 * *(axelSol->el) + 745.5426340882;
			}
			else if (*(axelSol->el) <= -73.799) // at this elevation angle, the result is not reliable anyway
			{
				*(axelSol->el) += -0.0418102295 * pow(*(axelSol->el), 2) - 5.4402816535 * *(axelSol->el) - 185.1646192938;
			}

			if (*(axelSol->el) > 90)
			{
				*(axelSol->el) = 90;
			}
			else if (*(axelSol->el) < -90)
			{
				*(axelSol->el) = -90;
			}
			// recompute xyz
			ae2xyzSol(*(axelSol->az), *(axelSol->el), axelSol);
		}

		/* save result */
		fprintf(fpw, "%s,%i,%lf,%lf,%lf,%lf,%lf\n", (*epochArray[i]).time, *(*epochArray[i]).numSat, *(axelSol->x), *(axelSol->y), *(axelSol->z), *(axelSol->az), *(axelSol->el));

		/* save to array */
		axelSolArray[i] = axelSol;
	}

	/*
		Statistics
	*/
	/* convergence */
	double mAxelX = 0, mAxelY = 0, mAxelZ = 0;

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		mAxelX += *axelSolArray[i]->x;
		mAxelY += *axelSolArray[i]->y;
		mAxelZ += *axelSolArray[i]->z;
	}

	double mXyzAxel[3] = {mAxelX, mAxelY, mAxelZ};
	normalizeXyz(mXyzAxel);
	double mAeAxel[2] = {0.0, 0.0};
	xyz2ae(mXyzAxel[0], mXyzAxel[1], mXyzAxel[2], mAeAxel);

	/* RMSE and standard deviation */
	double trueAntennaXyz[3];
	ae2xyz(TRUE_AZ, TRUE_EL, trueAntennaXyz);

	double rmsAxelX, rmsAxelY, rmsAxelZ;
	double sumAxelX = 0;
	double sumAxelY = 0;
	double sumAxelZ = 0;

	double stdAxelX, stdAxelY, stdAxelZ;
	double sumAxelX2 = 0;
	double sumAxelY2 = 0;
	double sumAxelZ2 = 0;

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		sumAxelX += pow((*axelSolArray[i]->x - trueAntennaXyz[0]), 2);
		sumAxelY += pow((*axelSolArray[i]->y - trueAntennaXyz[1]), 2);
		sumAxelZ += pow((*axelSolArray[i]->z - trueAntennaXyz[2]), 2);

		sumAxelX2 += pow((*axelSolArray[i]->x - mXyzAxel[0]), 2);
		sumAxelY2 += pow((*axelSolArray[i]->y - mXyzAxel[1]), 2);
		sumAxelZ2 += pow((*axelSolArray[i]->z - mXyzAxel[2]), 2);
	}

	/* rmse by component */
	rmsAxelX = sqrt(sumAxelX / *epochArrayIndex);
	rmsAxelX = rad2deg(asin(rmsAxelX));
	rmsAxelY = sqrt(sumAxelY / *epochArrayIndex);
	rmsAxelY = rad2deg(asin(rmsAxelY));
	rmsAxelZ = sqrt(sumAxelZ / *epochArrayIndex);
	rmsAxelZ = rad2deg(asin(rmsAxelZ));
	double rmsAxelA = sqrt(pow(rmsAxelX, 2) + pow(rmsAxelY, 2) + pow(rmsAxelZ, 2));
	double rmsAxelAz = sqrt(pow(rmsAxelA, 2) - pow(rmsAxelZ, 2));

	/* std by component */
	stdAxelX = sqrt(sumAxelX2 / *epochArrayIndex);
	stdAxelX = rad2deg(asin(stdAxelX));
	stdAxelY = sqrt(sumAxelY2 / *epochArrayIndex);
	stdAxelY = rad2deg(asin(stdAxelY));
	stdAxelZ = sqrt(sumAxelZ2 / *epochArrayIndex);
	stdAxelZ = rad2deg(asin(stdAxelZ));
	double stdAxelA = sqrt(pow(stdAxelX, 2) + pow(stdAxelY, 2) + pow(stdAxelZ, 2));
	double stdAxelAz = sqrt(pow(stdAxelA, 2) - pow(stdAxelZ, 2));

	printf("----------\nStatistics\n----------\nNumber of epochs\n%li\n", *epochArrayIndex);

	if (TRUE_EL >= -90 && TRUE_EL <= 90 && TRUE_AZ <= 360 && TRUE_AZ >= 0) // if antenna truth is provided by the user)
	{
		printf("\nAntenna truth by user input (E, N, U, Az, El)\n%.2f, %.2f, %.2f, %.2f°, %.2f°\n", trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2], (double)TRUE_AZ, (double)TRUE_EL);
	}

	printf("\nConvergence (E, N, U, Az, El)\n");
	printf("%.2f, %.2f, %.2f, %.2f°, %.2f°\n", mXyzAxel[0], mXyzAxel[1], mXyzAxel[2], mAeAxel[0], mAeAxel[1]);

	printf("\nStandard deviation\nE = %.2f°\nN = %.2f°\nU = %.2f°\nAz = %.2f°\nEl = %.2f°\n3D = %.2f°\n", stdAxelX, stdAxelY, stdAxelZ, stdAxelAz, stdAxelZ, stdAxelA);

	if (TRUE_EL >= -90 && TRUE_EL <= 90 && TRUE_AZ <= 360 && TRUE_AZ >= 0) // if antenna truth is provided by the user)
	{
		printf("\nRMSE\nE = %.2f°\nN = %.2f°\nU = %.2f°\nAz = %.2f°\nEl = %.2f°\n3D = %.2f°\n", rmsAxelX, rmsAxelY, rmsAxelZ, rmsAxelAz, rmsAxelZ, rmsAxelA);
	}

	/* close output file */
	fclose(fpw);

	/*
		free() file input 
	*/
	for (long int i = 0; i < satArrayIndex; i++)
	{
		free(satArray[i]->time);
		free(satArray[i]->prn);
		free(satArray[i]->az);
		free(satArray[i]->el);
		free(satArray[i]->snr);
		free(satArray[i]->snr2); // free attributes
		free(satArray[i]);		 // free Sat
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
		free(axelSolArray[i]->x);
		free(axelSolArray[i]->y);
		free(axelSolArray[i]->z);
		free(axelSolArray[i]->az);
		free(axelSolArray[i]->el);
		free(axelSolArray[i]);
	}
	free(axelSolArray);
	free(epochArrayIndex);

	clock_t endTime = clock();
	double timeSpent = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	printf("\nProgram execution time: %.2f seconds\n", timeSpent);

	/*
		exit
	*/
	return 0;

} // end of main()
