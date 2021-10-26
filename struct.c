#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "mathutil.h"
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
