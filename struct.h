#ifndef STRUCT_H
#define STRUCT_H

typedef struct Sat
{
	char *time; // GPS time
	char *prn;
	double *az;
	double *el;
	double *snr;
	double *snr2;
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

Sat *createSat(char *time, char *prn, double *az, double *el, double *snr, double *snr2);
Epoch *createEpoch(char *time, Sat **epochSatArray, int *numSat);
Sol *createSol(double *x, double *y, double *z, double *az, double *el);
void printEpochArray(Epoch **epochArray, long int numEpoch);

#endif
