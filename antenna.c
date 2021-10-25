#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include "mathutil.h"
#include "truth.h"

/*

Information for users of this program:
1. Modify the below configuration if necessary
2. Input file should not contain satellites without azimuth, elevation or SNR information.

gcc -Wall antenna.c mathutil.c -o antenna -lgsl -lgslcblas -lm
./antenna > output.txt
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

typedef struct Sat
{
	char *time; // GPS time
	char *satName;
	double *az;
	double *el;
	double *snr;
} Sat;

typedef struct Epoch
{
	char *time;
	Sat **satArrayInEpoch;
	int numSat;
	// Sol* solutionArray;
} Epoch;

Sat *createSat(char *time, char *satName, double *az, double *el, double *snr)
{
	Sat *satObj = malloc(sizeof(Sat));
	(*satObj).time = time;
	(*satObj).satName = satName;
	(*satObj).az = az;
	(*satObj).el = el;
	(*satObj).snr = snr;
	return satObj;
}

Epoch *createEpoch(char *time, Sat **satArrayInEpoch, int numSat)
{
	Epoch *epochObj = malloc(sizeof(Epoch));
	(*epochObj).time = time;
	(*epochObj).satArrayInEpoch = satArrayInEpoch;
	(*epochObj).numSat = numSat;
	return epochObj;
}

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
		printf("======== Epoch %s contains %i satellites/signals ========\n", (*epochArray[i]).time, (*epochArray[i]).numSat);
		double n = (*epochArray[i]).numSat;
		for (int j = 0; j < n; j++)
			printf("%s\t%s\t%lf\t%lf\t%lf\n", (*epochArray[i]).satArrayInEpoch[j]->time, (*epochArray[i]).satArrayInEpoch[j]->satName, *(*epochArray[i]).satArrayInEpoch[j]->az, *(*epochArray[i]).satArrayInEpoch[j]->el, *(*epochArray[i]).satArrayInEpoch[j]->snr);
	}
}

int main(void)
{
	/*
		create arrays
	*/
	char **timeArray; // A time array storing all unique times
	timeArray = (char **)malloc(MAX_NUM_EPOCH * sizeof(char *));
	/*
	for (long int i = 0; i < MAX_NUM_EPOCH; i++)
	{
		timeArray[i] = (char *)malloc(NUM_CHAR_DATE + NUM_CHAR_TIME + 2);
		timeArray[i][NUM_CHAR_DATE + NUM_CHAR_TIME + 1] = '\0';
	}
	*/
	Sat **satArray; // A sat array storing all satellite signals
	satArray = (Sat **)malloc(MAX_NUM_SIGNAL * sizeof(Sat *));
	/*
	for (long int i = 0; i < MAX_NUM_SIGNAL; i++)
		satArray[i] = (Sat *)malloc(sizeof(Sat));
	*/

	Epoch **epochArray; // A epoch array storing all epochs
	epochArray = (Epoch **)malloc(MAX_NUM_EPOCH * sizeof(Epoch *));
	/*
	for (long int i = 0; i < MAX_NUM_EPOCH; i++)
		epochArray[i] = (Epoch *)malloc(sizeof(Epoch));
	*/
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
			//printf("timeArray[timeArrayIndex] = %s\n", timeArray[timeArrayIndex]);
			timeArrayIndex++;
		}

		Sat *satObj = createSat(time, satName, az, el, snr);
		satArray[satArrayIndex] = satObj; // Add each sat to sat array

		satArrayIndex++;

		free(time1);
		free(time2);
		/*
		time1 = NULL;
		time2 = NULL;
		time = NULL;
		satName = NULL;
		az = NULL;
		el = NULL;
		snr = NULL;
		satObj = NULL;
		*/
		//free(satName);
		//free(az);
		//free(el);
		//free(snr);
		//free(satObj);
	}
	fclose(fp);
	free(line);
	printf("timeArrayIndex = %li satArrayIndex = %li\n", timeArrayIndex, satArrayIndex);

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
		for (int i = 0; i < MAX_NUM_SAT_EPOCH; i++)
			epochSatArray[i] = (Sat *)malloc(sizeof(Sat));
		*/

		/*
			Count number of sat in the array
		*/
		int epochSatArrayIndex = 0;
		for (long int j = 0; j < satArrayIndex; j++)
		{
			if (strcmp(timeArray[i], satArray[j]->time) == 0)
			{
				//printf("i = %li j = %li\n", i, j);
				//printf("epoch = %s\n", satArray[j]->time);
				epochSatArray[epochSatArrayIndex] = satArray[j];
				//printf("epoch = %lf\n", *epochSatArray[epochSatArrayIndex]->snr);

				epochSatArrayIndex++;
			}
		}

		if (epochSatArrayIndex >= 3)
		{ // if an epoch has fewer than 3 sat, do not record the epoch
			epochArray[epochArrayIndex] = createEpoch(timeArray[i], epochSatArray, epochSatArrayIndex);
			printf("epochArrayIndex = %i\n", epochSatArrayIndex);
			//printf("%s\t%s\t%lf\t%lf\t%lf\n", (*epochArray[i]).satArrayInEpoch[0].time, (*epochArray[i]).satArrayInEpoch[0].satName, *(*epochArray[i]).satArrayInEpoch[0].az, *(*epochArray[i]).satArrayInEpoch[0].el, *(*epochArray[i]).satArrayInEpoch[0].snr);
			epochArrayIndex++;
		}
	}

	//for (long int i = 0; i < epochArrayIndex; i++)
	//{
	//printf("======== Epoch %s contains %i satellites/signals ========\n", (*epochArray[i]).time, (*epochArray[i]).numSat);
	//double n = (*epochArray[i]).numSat;

	/******************** seg fault here *******************************/
	//printf("%lf\n", *(*epochArray[0]).satArrayInEpoch[1]->snr);
	//for (int j = 0; j < n; j++)
	//{
	//printf("%lf\n", *(*epochArray[i]).satArrayInEpoch[j].snr);
	//printf("%s\t%s\t%lf\t%lf\t%lf\n", (*epochArray[i]).satArrayInEpoch[j].time, (*epochArray[i]).satArrayInEpoch[j].satName, *(*epochArray[i]).satArrayInEpoch[j].az, *(*epochArray[i]).satArrayInEpoch[j].el, *(*epochArray[i]).satArrayInEpoch[j].snr);
	//}
	//}

	printEpochArray(epochArray, epochArrayIndex);

	//	printf("%s\t%s\t%lf\t%lf\t%lf\n", epochArray[i]->satArrayInEpoch[j].time, epochArray[i]->satArrayInEpoch[j].satName, *(epochArray[i]->satArrayInEpoch[j].az), *(epochArray[i]->satArrayInEpoch[j].el), *(epochArray[i]->satArrayInEpoch[j].snr));
	///////////////////////////////////////////////////////////////////////deleted
	/*
		free()
			This notice from valgrind is normal: Conditional jump or move depends on uninitialized value(s)
			This is because the arrays are not fully propagated; e.g. MAX_NUM_SIGNAL < actual num of signals
			Could use calloc() instead of malloc()
	
	for (long int i = 0; i < MAX_NUM_EPOCH; i++)
		free(timeArray[i]);
	free(timeArray);

	for (long int i = 0; i < MAX_NUM_SIGNAL; i++)
	{
		free(satArray[i]->satName);
		free(satArray[i]->az);
		free(satArray[i]->el);
		free(satArray[i]->snr);
		free(satArray[i]);
	}

	free(satArray);

	for (long int i = 0; i < MAX_NUM_EPOCH; i++)
		free(epochArray[i]->satArrayInEpoch);
	free(epochArray);
*/
	/*
		exit
	*/
	exit(0);

} // end of main()