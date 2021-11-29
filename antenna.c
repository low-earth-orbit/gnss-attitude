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
		fprintf(stderr, "Too many arguments. Usage:\n./antenna [input_1.txt] [input_2.txt]\n");
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

	if (argc == 1) // default input file
	{
		argc++;
	}
	for (int i = 1; i < argc; i++)
	{
		FILE *fpThis;
		if (i == 1) // first frequency
		{
			fpThis = fp;
		}
		else if (i == 2) // second frequency
		{
			fpThis = fp2;
		}
		fgets(line, MAX_NUM_CHAR_LINE, fpThis); // Skip header
		while (fgets(line, MAX_NUM_CHAR_LINE, fpThis) != NULL)
		{
			char *time1 = (char *)malloc(sizeof(char) * (NUM_CHAR_DATE + 1));
			char *time2 = (char *)malloc(sizeof(char) * (NUM_CHAR_TIME + 1));
			char *prn = (char *)malloc(sizeof(char) * (NUM_CHAR_SAT + 1));
			double *az = (double *)malloc(sizeof(double));
			double *el = (double *)malloc(sizeof(double));
			double *snrThis = (double *)malloc(sizeof(double));
			sscanf(line, "%s %s %s %lf %lf %lf", time1, time2, prn, az, el, snrThis);
			if (*az < 0 || *az > 360 || *el > 90 || *snrThis <= 0) // sanity check
			{
				printf("Skipped invalid data found in Az, El or SNR.\n");
				free(prn);
				free(az);
				free(el);
				free(snrThis);
			}
			else
			{
				char *time = concat(time1, time2);
				if (!SIMULATION)
				{
					if (i == 1)
					{
						adjSnr(prn, el, snrThis);
					}
					else if (i == 2)
					{
						adjSnr2(prn, el, snrThis);
					}
				}
				double *snrOther = (double *)malloc(sizeof(double));
				*snrOther = -1.0;
				if (i == 1)
				{
					satArray[satArrayIndex] = createSat(time, prn, az, el, snrThis, snrOther);
				}
				else if (i == 2)
				{
					satArray[satArrayIndex] = createSat(time, prn, az, el, snrOther, snrThis);
				}
				satArrayIndex++;
			}
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
	Sol **solArray = malloc(sizeof(Sol *) * *epochArrayIndex);

	fprintf(fpw, "Epoch (GPS Time), # of Signals, E, N, U, Az (deg), El (deg)\n");

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		/*	
			Observation equation X*b = cos(a) corresponds to X*c = y below
			See GNU Scientific Library Reference Manual for more: https://www.gnu.org/software/gsl/doc/html/lls.html
		*/
		int n = *(epochArray[i]->numSat); // number of observations in the epoch
		Sol *sol = malloc(sizeof(Sol));
		sol->x = calloc(1, sizeof(double));
		sol->y = calloc(1, sizeof(double));
		sol->z = calloc(1, sizeof(double));
		sol->az = calloc(1, sizeof(double));
		sol->el = calloc(1, sizeof(double));

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
			double xyz[3]; // x is E, y is N, z is U
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
				else if (*(*epochArray[i]).epochSatArray[j]->snr2 < 0)
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
		*(sol->x) = C(0);
		*(sol->y) = C(1);
		*(sol->z) = C(2);

		/* least squares stats */
		//printf("# covariance matrix:\n");
		//printf("[ %+.5e, %+.5e, %+.5e  \n", COV(0, 0), COV(0, 1), COV(0, 2));
		//printf("  %+.5e, %+.5e, %+.5e  \n", COV(1, 0), COV(1, 1), COV(1, 2));
		//printf("  %+.5e, %+.5e, %+.5e ]\n", COV(2, 0), COV(2, 1), COV(2, 2));
		//printf("# chisq = %g\n", chisq);
		/*
		double redChiSq = chisq / (n - 3); // reduced chisq = chisq / (# of signals - 3)
		double lsStdX = sqrt(redChiSq * COV(0, 0))));
		double lsStdY = sqrt(redChiSq * COV(1, 1))));
		double lsStdZ = sqrt(redChiSq * COV(2, 2))));
		*/

		/* free matrices for LS */
		gsl_matrix_free(X);
		gsl_vector_free(y);
		gsl_vector_free(w);
		gsl_vector_free(c);
		gsl_matrix_free(cov);

		/* normalize the resulting vector to get xyz solution*/
		normalize(sol);

		/* from xyz solution derive azimuth-elevation solution*/
		xyz2aeSol(*(sol->x), *(sol->y), *(sol->z), sol);

		/* apply convergence correction built by simulation */
		if (CONVERGENCE_CORRECTION)
		{
			if (*(sol->el) > 40.181)
			{
				*(sol->el) += 0.0000022678 * pow(*(sol->el), 3) - 0.0005537871 * pow(*(sol->el), 2) + 0.0455341172 * *(sol->el) - 1.2636551736;
			}
			else if (*(sol->el) > 0.562466)
			{
				*(sol->el) += 0.0000029756 * pow(*(sol->el), 3) - 0.0002119836 * pow(*(sol->el), 2) + 0.0133919561 * *(sol->el) - 0.5684236546;
			}
			else if (*(sol->el) > -33.9139)
			{
				*(sol->el) += -0.0000374714 * pow(*(sol->el), 3) - 0.0058097797 * pow(*(sol->el), 2) + 0.0088882136 * *(sol->el) - 0.5690671978;
			}
			else if (*(sol->el) > -61.2864)
			{
				*(sol->el) += 0.0000074641 * pow(*(sol->el), 3) + 0.0074059911 * pow(*(sol->el), 2) + 0.7488207967 * *(sol->el) + 11.0845716711;
			}
			else if (*(sol->el) > -73.799)
			{
				*(sol->el) += 0.0029469680 * pow(*(sol->el), 3) + 0.5621635563 * pow(*(sol->el), 2) + 35.6911758009 * *(sol->el) + 745.5426340882;
			}
			else if (*(sol->el) <= -73.799) // at this elevation angle, the result is not reliable anyway
			{
				*(sol->el) += -0.0418102295 * pow(*(sol->el), 2) - 5.4402816535 * *(sol->el) - 185.1646192938;
			}

			if (*(sol->el) > 90)
			{
				*(sol->el) = 90;
			}
			else if (*(sol->el) < -90)
			{
				*(sol->el) = -90;
			}
			// recompute xyz
			ae2xyzSol(*(sol->az), *(sol->el), sol);
		}

		/* save result */
		fprintf(fpw, "%s,%i,%lf,%lf,%lf,%lf,%lf\n", (*epochArray[i]).time, *(*epochArray[i]).numSat, *(sol->x), *(sol->y), *(sol->z), *(sol->az), *(sol->el));

		/* save to array */
		solArray[i] = sol;
	}

	/*
		Statistics
	*/
	/* convergence */
	double mX = 0.0, mY = 0.0, mZ = 0.0;

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		mX += *solArray[i]->x;
		mY += *solArray[i]->y;
		mZ += *solArray[i]->z;
	}
	mX /= *epochArrayIndex;
	mY /= *epochArrayIndex;
	mZ /= *epochArrayIndex;
	double mXyz[3] = {mX, mY, mZ};
	normalizeXyz(mXyz);
	double mAe[2] = {-1.0, -1.0};
	xyz2ae(mXyz[0], mXyz[1], mXyz[2], mAe);

	/* RMSE and standard deviation */
	double trueAntennaXyz[3];
	ae2xyz(TRUE_AZ, TRUE_EL, trueAntennaXyz);

	double rmsX, rmsY, rmsZ;
	double sumX = 0;
	double sumY = 0;
	double sumZ = 0;

	double stdX, stdY, stdZ;
	double sumX2 = 0;
	double sumY2 = 0;
	double sumZ2 = 0;

	for (long int i = 0; i < *epochArrayIndex; i++)
	{
		sumX += pow((*solArray[i]->x - trueAntennaXyz[0]), 2);
		sumY += pow((*solArray[i]->y - trueAntennaXyz[1]), 2);
		sumZ += pow((*solArray[i]->z - trueAntennaXyz[2]), 2);

		sumX2 += pow((*solArray[i]->x - mXyz[0]), 2);
		sumY2 += pow((*solArray[i]->y - mXyz[1]), 2);
		sumZ2 += pow((*solArray[i]->z - mXyz[2]), 2);
	}

	/* RMSE by component */
	rmsX = sqrt(sumX / *epochArrayIndex);
	rmsY = sqrt(sumY / *epochArrayIndex);
	rmsZ = sqrt(sumZ / *epochArrayIndex);
	double rmsA = rad2deg(asin(sqrt(pow(rmsX, 2) + pow(rmsY, 2) + pow(rmsZ, 2))));

	/* STD by component */
	stdX = sqrt(sumX2 / (*epochArrayIndex - 1)); // minus 1 since using sample mean
	stdY = sqrt(sumY2 / (*epochArrayIndex - 1));
	stdZ = sqrt(sumZ2 / (*epochArrayIndex - 1));
	double stdA = rad2deg(asin(sqrt(pow(stdX, 2) + pow(stdY, 2) + pow(stdZ, 2))));

	printf("----------\nStatistics\n----------\nNumber of epochs\n%li\n", *epochArrayIndex);

	if (TRUE_EL >= -90 && TRUE_EL <= 90 && TRUE_AZ <= 360 && TRUE_AZ >= 0) // if antenna truth is provided by the user)
	{
		printf("\nAntenna truth by user input (E, N, U, Az, El)\n%.2f, %.2f, %.2f, %.2f°, %.2f°\n", trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2], (double)TRUE_AZ, (double)TRUE_EL);
	}

	printf("\nConvergence (E, N, U, Az, El)\n");
	printf("%.2f, %.2f, %.2f, %.2f°, %.2f°\n", mXyz[0], mXyz[1], mXyz[2], mAe[0], mAe[1]);

	printf("\nStandard deviation\nE = %.4f\nN = %.4f\nU = %.4f\n3D = %.2f°\n", stdX, stdY, stdZ, stdA);

	if (TRUE_EL >= -90 && TRUE_EL <= 90 && TRUE_AZ <= 360 && TRUE_AZ >= 0) // if antenna truth is provided by the user)
	{
		printf("\nRoot-mean-square deviation\nE = %.4f\nN = %.4f\nU = %.4f\n3D = %.2f°\n", rmsX, rmsY, rmsZ, rmsA);
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
		free(solArray[i]->x);
		free(solArray[i]->y);
		free(solArray[i]->z);
		free(solArray[i]->az);
		free(solArray[i]->el);
		free(solArray[i]);
	}
	free(solArray);
	free(epochArrayIndex);

	clock_t endTime = clock();
	double timeSpent = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	printf("\nProgram execution time: %.2f seconds\n", timeSpent);

	/*
		exit
	*/
	return 0;

} // end of main()
