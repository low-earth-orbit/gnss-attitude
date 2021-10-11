/*
天地玄黃宇宙洪荒日月盈昃辰宿列張
*/

/* To compile
gcc convert.c antenna.c -lm -o antenna
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "convert.h"

/*
	Input file information
*/
#define INPUT_FILE_PATH "input.txt"
#define MAX_NUM_SIGNALS 10000
#define MAX_NUM_EPOCHES 1000
#define MAX_NUM_CHAR_LINE 100 // num of char in each record or line
#define MAX_NUM_SAT_EPOCH 100 // maximum number of satellites visible in an epoch
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
	Sat* satArray; // array of satellites with the same time
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

Epoch createEpoch(char* time, Sat* satArray, int numSat) {
	Epoch epochObj;
	epochObj.time = time;
	epochObj.satArray = satArray;
	epochObj.numSat = numSat;
	return epochObj;
}

bool isTimeInTimeArray(char* time, char** timeArray, int index) {
	for (int i = 0; i < index; i++) {
		if(strcmp(time, timeArray[i]) == 0) {
			// If duplicate time entry found
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
			printf("%s\t%s\t%lf\t%lf\t%lf\n", epochArray[i].satArray[j].time, epochArray[i].satArray[j].satName, epochArray[i].satArray[j].az, epochArray[i].satArray[j].el, epochArray[i].satArray[j].snr);
		}
	}
}
	
int main (void) {
	/*
		create arrays
	*/
	char** timeArray; // A time array storing all unique times
	Sat* satArray; // A sat array storing all satellite signals
	Epoch* epochArray; // A epoch array storing all epoches
	
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
	int timeArrayIndex = 0; // epoch array index is the same, or number of epoches
	
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
	char line[MAX_NUM_CHAR_LINE	+1];//+1 in case the max provided does not include the null-terminator
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

		free(time1);
		free(time2);// time1 and time2 used temporarily
		
		if(!isTimeInTimeArray(time, timeArray, timeArrayIndex)) { // Check if the current time is unique
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
		
		int epochSatArrayIndex = 0;
		for (int j = 0; j < satArrayIndex; j++){
			//printf("%i %s %s\n", satArrayIndex, timeArray[i], satArray[j].time);
			if (strcmp(timeArray[i], satArray[j].time)==0) {
				epochSatArray[epochSatArrayIndex] = satArray[j];
				epochSatArrayIndex++;
			}
		}
		Epoch epochObj = createEpoch(timeArray[i], epochSatArray, epochSatArrayIndex); //numSat initially 0
		epochArray[i] = epochObj;
	}
	//printEpochArray(epochArray, timeArrayIndex);
	
	/*
		SNR method calculations
	*/
	double xyz[timeArrayIndex][MAX_NUM_SAT_EPOCH][3];
	double xyzSol[timeArrayIndex][3]; // xyz solution array
	for (int i = 0; i < timeArrayIndex; i++) {
		for (int j = 0; j < epochArray[i].numSat; j++){
			/* get LOS vector */
			ae2xyz(epochArray[i].satArray[j].az, epochArray[i].satArray[j].el, xyz[i][j]); 
			//printf("epoch index = %i || sat index = %i || LOS vector (%lf, %lf %lf)\n", i, j, xyz[i][j][0], xyz[i][j][1], xyz[i][j][2]);
			
			/* weight by snr */ 
			xyz[i][j][0] *= epochArray[i].satArray[j].snr;
			xyz[i][j][1] *= epochArray[i].satArray[j].snr; 
			xyz[i][j][2] *= epochArray[i].satArray[j].snr;			
			//printf("epoch index = %i || sat index = %i || weighted LOS vector (%lf, %lf %lf)\n", i, j, xyz[i][j][0], xyz[i][j][1], xyz[i][j][2]);
			
			/* add to sum*/
			xyzSol[i][0] += xyz[i][j][0];
			xyzSol[i][1] += xyz[i][j][1]; 
			xyzSol[i][2] += xyz[i][j][2];
			//printf("epoch index = %i || sat index = %i || Accumulated sum of weighted LOS vector (%lf, %lf %lf)\n", i, j, xyzSol[i][0], xyzSol[i][1], xyzSol[i][2]);
		}
	}
	
	/*
		convert the sum to unit vector, finish xyz solution
	*/
	for (int i = 0; i < timeArrayIndex; i++) {
		normalizeXyz(xyzSol[i]);
		printf("epoch = %s || antenna boresight unit vector (x, y, z) = (%lf, %lf, %lf)\n", epochArray[i].time, xyzSol[i][0], xyzSol[i][1], xyzSol[i][2]);
	}
	
	/*
		from xyz solution derive azimuth-elevation solution
	*/
	double aeSol[timeArrayIndex][2]; // ae solution array
	for (int i = 0; i < timeArrayIndex; i++) {
		xyz2ae(xyzSol[i][0], xyzSol[i][1], xyzSol[i][2], aeSol[i]);
		printf("epoch = %s || antenna boresight azimuth and elevation angle in degrees (az, el) = (%lf, %lf)\n", epochArray[i].time, aeSol[i][0], aeSol[i][1]);
	}
	
	
	/*
		free()
			This notice from valgrind is normal:
				Conditional jump or move depends on uninitialised value(s)
			This is because the arrays are not fully allocated; e.g. MAX_NUM_SIGNALS < actual num of signals
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
		free(epochArray[i].satArray);
	free(epochArray);
	
	/*
		exit
	*/
	exit(0);
} // end of main()