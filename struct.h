#ifndef STRUCT_H
#define STRUCT_H

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
	Sat **epochSatArray;
	int *numSat;
} Epoch;

typedef struct Sol
{
	double *x, *y, *z, *az, *el;
} Sol;

Sat *createSat(char *time, char *satName, double *az, double *el, double *snr);
Epoch *createEpoch(char *time, Sat **epochSatArray, int *numSat);
Sol *createSol(double *x, double *y, double *z, double *az, double *el);

#endif
