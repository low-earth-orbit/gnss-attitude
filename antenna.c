#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "convert.h" // Coordinate transformation methods
/*

Information for users of this program:
1. Modify the below configuration if necessary
2. Do a sanity check for the input file. Input file should closely resemble "input_example.txt".
	Input file should not contain satellites without azimuth, elevation or SNR information.

gcc antenna.c convert.c -lm -o antenna
./antenna > output.txt
*/

/* Configuration */
#define INPUT_FILE_PATH "input.txt"
#define MAX_NUM_EPOCHES 86400
#define TRUE_EL -80 // if true elevation angle is know, input here for statistics (in degrees) 

/* Usually no need to change*/
#define EXCLUDE_GEO_SAT true // true: exclude geostationary satellites
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
} Sat; // satellite signal

typedef struct {
	char* time;
	Sat* satArrayInEpoch; // array of satellites with the same time
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
		if(strcmp(str, array[i]) == 0) {
			// If duplicate str entry found
			return true;
		}
	}
	return false;
}

/*
bool isTimeInSatArray(char* time, Sat* satArray, int index) {
	for (int i = 0; i < index; i++) {
		if(strcmp(time, satArray[i]) == 0) {
			// If duplicate time entry found
			return true;
		}
	}
	return false;
}
*/

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

void printEpochArray(Epoch* epochArray, int numEpoch) {
 	for (int i = 0; i < numEpoch; i++) {
		printf("======== Epoch %s contains %i satellites/signals ========\n", epochArray[i].time, epochArray[i].numSat);
		for (int j = 0; j < epochArray[i].numSat; j++){
			printf("%s\t%s\t%lf\t%lf\t%lf\n", epochArray[i].satArrayInEpoch[j].time, epochArray[i].satArrayInEpoch[j].satName, epochArray[i].satArrayInEpoch[j].az, epochArray[i].satArrayInEpoch[j].el, epochArray[i].satArrayInEpoch[j].snr);
		}
	}
}

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
	
	/*
		array indexes
	*/
	int satArrayIndex = 0;
	int timeArrayIndex = 0; // epoch array index is the same, or number of epochs
	
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
	char line[MAX_NUM_CHAR_LINE+1];//+1 in case the max provided does not include the null-terminator
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
			/*
			printf("timeArray[timeArrayIndex]: %s\n", timeArray[timeArrayIndex]);
			printf("timeArrayIndex: %d\n", timeArrayIndex);
			*/
			timeArrayIndex++;
		}
		
		Sat satObj = createSat(time, satName, az, el, snr);
		satArray[satArrayIndex] = satObj; // Add each sat to sat array
		
		/*
		printf("satArray[satArrayIndex].time: %s\t", satArray[satArrayIndex].time);
		printf("satArray[satArrayIndex].satName: %s\t", satArray[satArrayIndex].satName);
		printf("satArrayIndex: %d\n", satArrayIndex);
		*/
		satArrayIndex++;

		free(time1);
		free(time2);// time1 and time2 used temporarily
	}
	
	/*
		close file
	*/
	fclose(fp);
	
	/* 
		epoch array
		for each unique time:
			1) create an Epoch object;
			2) load it to epochArray
	*/
 	for (int i = 0; i < timeArrayIndex; i++) {
		Sat* epochSatArray = (Sat*)malloc(MAX_NUM_SAT_EPOCH*sizeof(Sat));
		if (epochSatArray == NULL) {
			fprintf(stderr, "malloc() failed for creating epochSatArray\n");
			exit(-1);
		}
		
		int epochSatArrayIndex = 0; //num of sat in an epoch initially 0
		for (int j = 0; j < satArrayIndex; j++){
			//printf("%i %s %s\n", satArrayIndex, timeArray[i], satArray[j].time);
			if (strcmp(timeArray[i], satArray[j].time)==0) {
				epochSatArray[epochSatArrayIndex] = satArray[j];
				epochSatArrayIndex++;
			}
		}
		Epoch epochObj = createEpoch(timeArray[i], epochSatArray, epochSatArrayIndex); 
		epochArray[i] = epochObj;
	}
	//print lines read in by the program
	//printEpochArray(epochArray, timeArrayIndex);
	
	/*
		Duncan method (Duncan & Dunn, 1998) -- Vector sum of signal-to-noise ratio (SNR) weighted line-of-sight (LOS) vectors
	*/
	double xyzDun[timeArrayIndex][MAX_NUM_SAT_EPOCH][3];
	double xyzDunSol[timeArrayIndex][3]; // xyz solution array
	double aeDuncanSol[timeArrayIndex][2]; // ae solution array

	// print header of the output
	printf("================== Duncan method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (int i = 0; i < timeArrayIndex; i++) {
		for (int j = 0; j < epochArray[i].numSat; j++){
			/* calculate LOS vector from azimuth and elevation*/
			ae2xyz(epochArray[i].satArrayInEpoch[j].az, epochArray[i].satArrayInEpoch[j].el, xyzDun[i][j]); 
			//printf("epoch index = %i || sat index = %i || LOS vector (%lf, %lf %lf)\n", i, j, xyzDun[i][j][0], xyzDun[i][j][1], xyzDun[i][j][2]);
			
			/* weight LOS vector by SNR */
			xyzDun[i][j][0] *= epochArray[i].satArrayInEpoch[j].snr;
			xyzDun[i][j][1] *= epochArray[i].satArrayInEpoch[j].snr; 
			xyzDun[i][j][2] *= epochArray[i].satArrayInEpoch[j].snr;
			//printf("epoch index = %i || sat index = %i || weighted LOS vector (%lf, %lf %lf)\n", i, j, xyzDun[i][j][0], xyzDun[i][j][1], xyzDun[i][j][2]);

			/* add to vector sum */
			xyzDunSol[i][0] += xyzDun[i][j][0];
			xyzDunSol[i][1] += xyzDun[i][j][1]; 
			xyzDunSol[i][2] += xyzDun[i][j][2];
			//printf("epoch index = %i || sat index = %i || Accumulated sum of weighted LOS vector (%lf, %lf %lf)\n", i, j, xyzDunSol[i][0], xyzDunSol[i][1], xyzDunSol[i][2]);
		}

		/* convert the sum to unit vector to get xyz solution*/
		//printf("Before normalize %s,%i,%lf,%lf,%lf\n", epochArray[i].time, epochArray[i].numSat, xyzDunSol[i][0], xyzDunSol[i][1], xyzDunSol[i][2]);
		normalizeXyz(xyzDunSol[i]);

		/* from xyz solution derive azimuth-elevation solution */
		xyz2ae(xyzDunSol[i][0], xyzDunSol[i][1], xyzDunSol[i][2], aeDuncanSol[i]);

		/* print Duncan method result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", epochArray[i].time, epochArray[i].numSat, xyzDunSol[i][0], xyzDunSol[i][1], xyzDunSol[i][2], aeDuncanSol[i][0], aeDuncanSol[i][1]);
		
	}//End of Duncan method
	
	/*
		Geometry method -- Vector sum of non-weighted line-of-sight (LOS) vectors with geometry adjustment for elevation angle
	*/
	double xyzGeo[timeArrayIndex][MAX_NUM_SAT_EPOCH][3];
	double xyzGeoSol[timeArrayIndex][3]; // xyz solution array
	double aeGeoSol[timeArrayIndex][2]; // ae solution array

	// print header of the output
	printf("================== Geometry method ==================\nEpoch(GPST),#Sat,X(E),Y(N),Z(U),Az(deg),El(deg)\n");

	for (int i = 0; i < timeArrayIndex; i++) {
		for (int j = 0; j < epochArray[i].numSat; j++){
			/* calculate LOS vector from azimuth and elevation*/
			ae2xyz(epochArray[i].satArrayInEpoch[j].az, epochArray[i].satArrayInEpoch[j].el, xyzGeo[i][j]); 
			//printf("epoch index = %i || sat index = %i || LOS vector (%lf, %lf %lf)\n", i, j, xyzGeo[i][j][0], xyzGeo[i][j][1], xyzGeo[i][j][2]);

			/* add to vector sum */
			xyzGeoSol[i][0] += xyzGeo[i][j][0];
			xyzGeoSol[i][1] += xyzGeo[i][j][1]; 
			xyzGeoSol[i][2] += xyzGeo[i][j][2];
			//printf("epoch index = %i || sat index = %i || Accumulated sum of LOS vector (%lf, %lf %lf)\n", i, j, xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2]);
		}

		/* convert the sum to unit vector to get xyz solution*/
		//printf("Before normalize %s,%i,%lf,%lf,%lf\n", epochArray[i].time, epochArray[i].numSat, xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2]);
		normalizeXyz(xyzGeoSol[i]);
		//printf("After normalize %s,%i,%lf,%lf,%lf\n", epochArray[i].time, epochArray[i].numSat, xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2]);
		
		/* from xyz solution derive azimuth-elevation solution */
		xyz2ae(xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2], aeGeoSol[i]);

		/* adjust the elevation angle by (2 * el - 90) */
		//printf("Before adjustment el = %lf\n", aeGeoSol[i][1]);
		aeGeoSol[i][1] = aeGeoSol[i][1]*2.0 - 90.0;
		//printf("After adjustment el = %lf\n", aeGeoSol[i][1]);
		
		/* recompute xyz solution using adjusted elevation angle */
		ae2xyz(aeGeoSol[i][0], aeGeoSol[i][1], xyzGeoSol[i]);

		/* print Geometry method result */
		printf("%s,%i,%lf,%lf,%lf,%lf,%lf\n", epochArray[i].time, epochArray[i].numSat, xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2], aeGeoSol[i][0], aeGeoSol[i][1]);
		
	}//End of Geometry method
	
	/*
		Statistics if true antenna boresight (E, N, U) is known 
	*/
	if (TRUE_EL >= -90 && TRUE_EL <= 90){
		double rmsDun;
		double sumDun = 0;

		double rmsGeo;
		double sumGeo = 0;
		
		for (int i = 0; i < timeArrayIndex; i++) {
			sumDun += pow(spDist(xyzDunSol[i][0], xyzDunSol[i][1], xyzDunSol[i][2], 0, -cos(deg2rad(TRUE_EL)), sin(deg2rad(TRUE_EL))),2);
			sumGeo += pow(spDist(xyzGeoSol[i][0], xyzGeoSol[i][1], xyzGeoSol[i][2], 0, -cos(deg2rad(TRUE_EL)), sin(deg2rad(TRUE_EL))),2);
		}
		
		rmsDun = sqrt(sumDun/timeArrayIndex);
		rmsDun = rad2deg(rmsDun);
		
		rmsGeo = sqrt(sumGeo/timeArrayIndex);
		rmsGeo = rad2deg(rmsGeo);
		
		printf ("================== Statistics ==================\n%i epochs, antenna @ %i deg\nRMS Duncan method = %lf deg\nRMS Geometry method = %lf deg\n", timeArrayIndex, TRUE_EL, rmsDun, rmsGeo);
	}


	/*
		free()
			This notice from valgrind is normal:
				Conditional jump or move depends on uninitialised value(s)
			This is because the arrays are not fully propagated; e.g. MAX_NUM_SIGNALS < actual num of signals
			Could use calloc() instead of malloc()
	*/
	for (int i = 0; i < MAX_NUM_EPOCHES; i++) {
		free(timeArray[i]);
	}
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