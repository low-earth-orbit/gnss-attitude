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
#define MAX_NUM_EPOCHES 1000

/* Usually no need to change*/
#define MAX_NUM_SAT_EPOCH 100 // maximum number of satellites visible in an epoch
#define MAX_NUM_SIGNALS MAX_NUM_EPOCHES*MAX_NUM_SAT_EPOCH
#define MAX_NUM_CHAR_LINE 100 // num of char in each record or line
#define NUM_CHAR_DATE 10
#define NUM_CHAR_TIME 10
#define NUM_CHAR_SAT 3

typedef struct {
	char* time; // GPS time
	char* satName;
	double az;
	double el;
	double snr;
} Sat;

typedef struct {
	char* time;
	Sat* satArrayInEpoch;
	int numSat;
	// Sol* solutionArray;
} Epoch;

Sat createSat(char* time, char* satName, double az, double el, double snr) {
	Sat satObj;
	satObj.time = time;
	satObj.satName = satName;
	satObj.az = az;
	satObj.el = el;
	satObj.snr = snr;
	return satObj;
}

Epoch createEpoch(char* time, Sat* satArrayInEpoch, int numSat) {
	Epoch epochObj;
	epochObj.time = time;
	epochObj.satArrayInEpoch = satArrayInEpoch;
	epochObj.numSat = numSat;
	return epochObj;
}

bool isStrInArray(char* str, char** array, int index) {
	for (int i = 0; i < index; i++) {
		if(strcmp(str, array[i]) == 0)
			return true;
	}
	return false;
}

char* concat(const char* str1, const char* str2)
{
	char* str;
	str = (char*)malloc(strlen(str1) + strlen(str2) +2); 
	// +1 for the null-terminator; +1 for the seperator
	if (str == NULL) {
		fprintf(stderr, "malloc() failed in concat()\n");
	}
	else {
		strcpy(str, str1);
		strcat(str, " ");
		strcat(str, str2);
	}
	return str;
}

/*
void printEpochArray(Epoch* epochArray, int numEpoch) {
 	for (int i = 0; i < numEpoch; i++) {
		printf("======== Epoch %s contains %i satellites/signals ========\n", epochArray[i].time, epochArray[i].numSat);
		for (int j = 0; j < epochArray[i].numSat; j++){
			printf("%s\t%s\t%lf\t%lf\t%lf\n", epochArray[i].satArrayInEpoch[j].time, epochArray[i].satArrayInEpoch[j].satName, epochArray[i].satArrayInEpoch[j].az, epochArray[i].satArrayInEpoch[j].el, epochArray[i].satArrayInEpoch[j].snr);
		}
	}
}
*/

int main (void) {
	/*
		create arrays
	*/
	char** timeArray; // A time array storing all unique times
	Sat* satArray; // A sat array storing all satellite signals
	Epoch* epochArray; // A epoch array storing all epochs

	timeArray = (char**)malloc(MAX_NUM_EPOCHES*sizeof(char*));
	if (timeArray == NULL) {
		fprintf(stderr, "malloc() failed for creating timeArray\n");
		exit(-1);
	}
	for (int i = 0; i < MAX_NUM_EPOCHES; i++) {
		timeArray[i] = (char*)malloc(NUM_CHAR_DATE + NUM_CHAR_TIME + 2);
		if (timeArray[i] == NULL) {
			fprintf(stderr, "malloc() failed for creating timeArray element\n");
			exit(-1);
		}
	}
	
	satArray = (Sat*)malloc(MAX_NUM_SIGNALS*sizeof(Sat));
	if (satArray == NULL) {
		fprintf(stderr, "malloc() failed for creating satArray\n");
		exit(-1);
	}
	
	epochArray = (Epoch*)malloc(MAX_NUM_EPOCHES*sizeof(Epoch));
	if (epochArray == NULL) {
		fprintf(stderr, "malloc() failed for creating epochArray\n");
		exit(-1);
	}
	
	int satArrayIndex = 0;
	int timeArrayIndex = 0;
	
	/*
		open file for reading
	*/
	FILE *fp = NULL;
	fp = fopen(INPUT_FILE_PATH, "r");
	if(fp == NULL) {
		perror("Error opening file\n");
		return(-1);
	}

	/*
		File processing to get timeArray and satArray
	*/
	char line[MAX_NUM_CHAR_LINE+1];
	fgets(line, sizeof(line), fp); // Skip header
	while(fgets(line, sizeof(line), fp) != NULL) {
		char* time1 = (char*)malloc(sizeof(char)*(NUM_CHAR_DATE +1));
		if (time1 == NULL) {
			fprintf(stderr, "malloc() failed for creating time1\n");
			exit(-1);
		}
		
		char* time2 = (char*)malloc(sizeof(char)*(NUM_CHAR_TIME +1));
		if (time2 == NULL) {
			fprintf(stderr, "malloc() failed for creating time2\n");
			exit(-1);
		}
		
		char* satName = (char*)malloc(sizeof(char)*(NUM_CHAR_SAT +1));
		if (satName == NULL) {
			fprintf(stderr, "malloc() failed for creating satName\n");
			exit(-1);
		}
		
		double az;
		double el;
		double snr;

		sscanf(line, "%s %s %s %lf %lf %lf", time1, time2, satName, &az, &el, &snr);
		
		char* time = concat(time1, time2);
		if (time == NULL)
			exit(-1);//failed to allocate memory in concat() above
		
		if(!isStrInArray(time, timeArray, timeArrayIndex)) { // Check if the current time is unique
			strcpy (timeArray[timeArrayIndex], time);
			timeArrayIndex++;
		}
		
		Sat satObj = createSat(time, satName, az, el, snr);
		satArray[satArrayIndex] = satObj; // Add each sat to sat array

		satArrayIndex++;

		free(time1);
		free(time2);
	}
	fclose(fp);
	
	/* 
		epoch array
		for each unique time:
			1) create an Epoch object;
			2) load it to epochArray
	*/
	int epochArrayIndex = 0; // number epoch can be less than number of unique time 
 	for (int i = 0; i < timeArrayIndex; i++) {

		Sat* epochSatArray = (Sat*)malloc(MAX_NUM_SAT_EPOCH*sizeof(Sat));
		if (epochSatArray == NULL) {
			fprintf(stderr, "malloc() failed for creating epochSatArray\n");
			exit(-1);
		}
		/*
			Count number of sat in the array
		*/
		int epochSatArrayIndex = 0;
		for (int j = 0; j < satArrayIndex; j++){
			if (strcmp(timeArray[i], satArray[j].time)==0) {
				epochSatArray[epochSatArrayIndex] = satArray[j];
				epochSatArrayIndex++;
			}
		}
		
		if (epochSatArrayIndex >= 3) {// if an epoch has fewer than 3 sat, do not record the epoch
			Epoch epochObj = createEpoch(timeArray[i], epochSatArray, epochSatArrayIndex); 
			epochArray[epochArrayIndex] = epochObj;
			epochArrayIndex++;
		}
		else {
			free(epochSatArray);
		}

	}
	//printEpochArray(epochArray, epochArrayIndex);
	
	/*
		Duncan's method (Duncan & Dunn, 1998) -- Vector sum of signal-to-noise ratio (SNR) weighted line-of-sight (LOS) vectors
	*/
	double xyzDun[epochArrayIndex][MAX_NUM_SAT_EPOCH][3];
	double xyzDunSol[epochArrayIndex][3]; // xyz solution array
	double aeDunSol[epochArrayIndex][2]; // ae solution array

	// print header of the output
	printf("================== Duncan's method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (int i = 0; i < epochArrayIndex; i++) {
		for (int j = 0; j < epochArray[i].numSat; j++){
			/* calculate LOS vector from azimuth and elevation*/
			ae2xyz(epochArray[i].satArrayInEpoch[j].az, epochArray[i].satArrayInEpoch[j].el, xyzDun[i][j]); 
			
			/* weight LOS vector by SNR */
			xyzDun[i][j][0] *= epochArray[i].satArrayInEpoch[j].snr;
			xyzDun[i][j][1] *= epochArray[i].satArrayInEpoch[j].snr; 
			xyzDun[i][j][2] *= epochArray[i].satArrayInEpoch[j].snr;

			/* add to vector sum */
			xyzDunSol[i][0] += xyzDun[i][j][0];
			xyzDunSol[i][1] += xyzDun[i][j][1]; 
			xyzDunSol[i][2] += xyzDun[i][j][2];
		}

		/* convert the sum to unit vector to get xyz solution*/
		normalizeXyz(xyzDunSol[i]);

		/* from xyz solution derive azimuth-elevation solution */
		xyz2ae(xyzDunSol[i][0], xyzDunSol[i][1], xyzDunSol[i][2], aeDunSol[i]);

		/* print Duncan's method result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", epochArray[i].time, epochArray[i].numSat, xyzDunSol[i][0], xyzDunSol[i][1], xyzDunSol[i][2], aeDunSol[i][0], aeDunSol[i][1]);
		
	}
	
	/*
		Geometry method -- Vector sum of non-weighted line-of-sight (LOS) vectors with geometry adjustment for elevation angle. Duncan's method is biased toward the spherical area where the satellite signals come from. If the antenna's elevation angle is negative, Duncan's method yields a positive elevation angle estimate. This method is designed by the author of the program to address the geometry issue existing in Duncan's method.
	*/
	double xyzGeo[epochArrayIndex][MAX_NUM_SAT_EPOCH][3];
	double xyzGeoSol[epochArrayIndex][3]; // xyz solution array
	double aeGeoSol[epochArrayIndex][2]; // ae solution array

	// print header of the output
	printf("================== Geometry method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (int i = 0; i < epochArrayIndex; i++) {
		for (int j = 0; j < epochArray[i].numSat; j++){
			/* calculate LOS vector from azimuth and elevation*/
			ae2xyz(epochArray[i].satArrayInEpoch[j].az, epochArray[i].satArrayInEpoch[j].el, xyzGeo[i][j]); 

			/* add to vector sum */
			xyzGeoSol[i][0] += xyzGeo[i][j][0];
			xyzGeoSol[i][1] += xyzGeo[i][j][1]; 
			xyzGeoSol[i][2] += xyzGeo[i][j][2];
		}

		/* convert the sum to unit vector to get xyz solution*/
		normalizeXyz(xyzGeoSol[i]);
		
		/* from xyz solution derive azimuth-elevation solution */
		xyz2ae(xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2], aeGeoSol[i]);

		/* adjust the elevation angle by (2 * el - 90) */
		aeGeoSol[i][1] = aeGeoSol[i][1]*2.0 - 90.0;
		
		/* recompute xyz solution using adjusted elevation angle */
		ae2xyz(aeGeoSol[i][0], aeGeoSol[i][1], xyzGeoSol[i]);

		/* print Geometry method result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", epochArray[i].time, epochArray[i].numSat, xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2], aeGeoSol[i][0], aeGeoSol[i][1]);
		
	}
	
	/*
		Axelrad's method (1999) -- Compared to Duncan's method, this is the proper use of SNR in determining antenna boresight vector. It requires antenna gain mapping (the relationship between off-boresight angle and SNR for the antenna) and adjustment to measured SNR.
	*/
	double xyz[epochArrayIndex][MAX_NUM_SAT_EPOCH][3];// line of sight vectors 
	double xyzSol[epochArrayIndex][3]; // xyz solution array
	double aeSol[epochArrayIndex][2]; // ae solution array
	
	// print header of the output
	printf("================== Axelrad's method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");
	for (int i = 0; i < epochArrayIndex; i++) {
		
		/*	
			Observation equation: transpose(s) b = cos(alpha) 
										corresponds to X c = y below
			See GNU Scientific Library Reference Manual for more: https://www.gnu.org/software/gsl/doc/html/lls.html#examples
		*/
		int n = epochArray[i].numSat; // number of observations
		double chisq;
		gsl_matrix *X, *cov;
		gsl_vector *y, *w, *c;
		X = gsl_matrix_alloc (n, 3);
		y = gsl_vector_alloc (n); // n*1 matrix
		w = gsl_vector_alloc (n); // n*n matrix
		c = gsl_vector_alloc (3); // (x, y, z) this is the boresight vector
		cov = gsl_matrix_alloc (3, 3); // cov = inverse(transpose(X) W X)
		
		for (int j = 0; j < epochArray[i].numSat; j++){
			/* calculate LOS vector from azimuth and elevation*/
			ae2xyz(epochArray[i].satArrayInEpoch[j].az, epochArray[i].satArrayInEpoch[j].el, xyz[i][j]);
			
			double sigmaSnr = 0.5+(3-0.5)*(epochArray[i].satArrayInEpoch[j].snr - 35)/(50-35);
			if (sigmaSnr < 0.5)
				sigmaSnr = 0.5; // catch the case that sigma <= 0
			double sigma = sigmaSnr/15.0; // uncertainty is a function of SNR
			double cosAlpha = (epochArray[i].satArrayInEpoch[j].snr - 35)/(50-35); // the mapping function is snr = (MAX_SNR-MIN_SNR)*cos(spd)+MIN_SNR; 
			
			// Set each observation equation 
			gsl_matrix_set (X, j, 0, xyz[i][j][0]);
			gsl_matrix_set (X, j, 1, xyz[i][j][1]);
			gsl_matrix_set (X, j, 2, xyz[i][j][2]);
			gsl_vector_set (y, j, cosAlpha);
			gsl_vector_set (w, j, 1.0/(sigma*sigma));
		}

		gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 3);
		gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
		gsl_multifit_linear_free (work);
		
		#define C(i) (gsl_vector_get(c,(i)))
		#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

		/* save best fit */
		xyzSol[i][0] = C(0);
		xyzSol[i][1] = C(1);
		xyzSol[i][2] = C(2);
		
		gsl_matrix_free (X);
		gsl_vector_free (y);
		gsl_vector_free (w);
		gsl_vector_free (c);
		gsl_matrix_free (cov);
		
		/* normalize the resulting vector to get xyz solution*/
		normalizeXyz(xyzSol[i]);
		
		/* from xyz solution derive azimuth-elevation solution */
		xyz2ae(xyzSol[i][0], xyzSol[i][1], xyzSol[i][2], aeSol[i]);
		
		/* print Axelrad's method result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", epochArray[i].time, epochArray[i].numSat, xyzSol[i][0], xyzSol[i][1], xyzSol[i][2], aeSol[i][0], aeSol[i][1]);
	}
	
	
	/*
		Statistics if true antenna boresight (E, N, U) is known 
	*/
	if (TRUE_EL >= -90 && TRUE_EL <= 90 && TRUE_AZ <= 360 && TRUE_AZ>=0){
		double rmsDun;
		double sumDun = 0;

		double rmsGeo;
		double sumGeo = 0;
		
		double rms;
		double sum = 0;
		
		for (int i = 0; i < epochArrayIndex; i++) {
			double trueAntennaXyz[3];
			ae2xyz(TRUE_AZ,TRUE_EL,trueAntennaXyz);
			sumDun += pow(spDist(xyzDunSol[i][0], xyzDunSol[i][1], xyzDunSol[i][2], trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]),2);
			sumGeo += pow(spDist(xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2], trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]),2);
			sum += pow(spDist(xyzSol[i][0], xyzSol[i][1], xyzSol[i][2], trueAntennaXyz[0], trueAntennaXyz[1], trueAntennaXyz[2]),2);
			
		}
		
		rmsDun = sqrt(sumDun/epochArrayIndex);
		rmsDun = rad2deg(rmsDun);
		
		rmsGeo = sqrt(sumGeo/epochArrayIndex);
		rmsGeo = rad2deg(rmsGeo);
		
		rms = sqrt(sum/epochArrayIndex);
		rms = rad2deg(rms);
		
		printf ("================== Statistics ==================\n%i epochs, antenna @ %i deg\nRMS Duncan's method = %lf deg\nRMS Geometry method = %lf deg\nRMS Axelrad's method = %lf deg\n", epochArrayIndex, TRUE_EL, rmsDun, rmsGeo, rms);
	}


	/*
		free()
			This notice from valgrind is normal:
				Conditional jump or move depends on uninitialised value(s)
			This is because the arrays are not fully propagated; e.g. MAX_NUM_SIGNALS < actual num of signals
			Could use calloc() instead of malloc()
	*/
	for (int i = 0; i < MAX_NUM_EPOCHES; i++)
		free(timeArray[i]);
	free(timeArray);
	
 	for (int i = 0; i < MAX_NUM_SIGNALS; i++) {
		free(satArray[i].time);
		free(satArray[i].satName);
	}
	free(satArray);

	for (int i = 0; i < MAX_NUM_EPOCHES; i++)
		free(epochArray[i].satArrayInEpoch);
	free(epochArray);
	
	/*
		exit
	*/
	exit(0);
} // end of main()