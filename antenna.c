#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include "mathutil.h"
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

bool isStrInArray(char *str, char **array, long int index)
{
	for (long int i = 0; i < index; i++)
	{
		if (strcmp(str, array[i]) == 0)
			return true;
	}
	return false;
}

char *concat(const char *str1, const char *str2)
{
	char *str;
	str = (char *)malloc(strlen(str1) + strlen(str2) + 2);
	// +1 for the null-terminator; +1 for the seperator
	if (str == NULL)
	{
		fprintf(stderr, "malloc() failed in concat()\n");
	}
	else
	{
		strcpy(str, str1);
		strcat(str, " ");
		strcat(str, str2);
	}
	return str;
}

void printEpochArray(Epoch **epochArray, long int numEpoch)
{
	for (long int i = 0; i < numEpoch; i++)
	{
		int n = *(epochArray[i]->numSat);
		printf("======== Epoch %s contains %i satellite signals ========\n", (*epochArray[i]).time, n + 1);
		for (int j = 0; j < n; j++)
			printf("%s\t%s\t%lf\t%lf\t%lf\n", (*epochArray[i]).epochSatArray[j]->time, (*epochArray[i]).epochSatArray[j]->satName, *(*epochArray[i]).epochSatArray[j]->az, *(*epochArray[i]).epochSatArray[j]->el, *(*epochArray[i]).epochSatArray[j]->snr);
	}
}

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
	/*
	double xyzDun[epochArrayIndex][MAX_NUM_SAT_EPOCH][3];
	double xyzDunSol[epochArrayIndex][3]; // E, N, U solution array
	double aeDunSol[epochArrayIndex][2];  // azimuth elevation solution array
	*/
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
			double xyz[3];
			/* calculate LOS vector from azimuth and elevation*/
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
	/////////////// OK above

	/////////////////////////////////////////////////////////////////////// deleted
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
		free() Duncan's solution
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

	/*
		exit
	*/
	exit(0);

} // end of main()