#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "util.h"
#include "struct.h"

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

Epoch *createEpoch(char *time, Sat **epochSatArray, int *numSat)
{
	Epoch *epochObj = malloc(sizeof(Epoch));
	(*epochObj).time = time;
	(*epochObj).epochSatArray = epochSatArray;
	(*epochObj).numSat = numSat;
	return epochObj;
}

Sol *createSol(double *x, double *y, double *z, double *az, double *el)
{
	Sol *solObj = malloc(sizeof(Sol));
	(*solObj).x = x;
	(*solObj).y = y;
	(*solObj).z = z;
	(*solObj).az = az;
	(*solObj).el = el;
	return solObj;
}

void printEpochArray(Epoch **epochArray, long int numEpoch)
{
	for (long int i = 0; i < numEpoch; i++)
	{
		int n = *(epochArray[i]->numSat);
		printf("======== Epoch %s contains %i satellite signals ========\n", (*epochArray[i]).time, n);
		for (int j = 0; j < n; j++)
			printf("%s\t%s\t%lf\t%lf\t%lf\n", (*epochArray[i]).epochSatArray[j]->time, (*epochArray[i]).epochSatArray[j]->satName, *(*epochArray[i]).epochSatArray[j]->az, *(*epochArray[i]).epochSatArray[j]->el, *(*epochArray[i]).epochSatArray[j]->snr);
	}
}