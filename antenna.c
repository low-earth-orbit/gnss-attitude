#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*
	Input file information
*/
#define INPUT_FILE_PATH "new.txt"
#define MAX_NUM_SIGNALS 10000
#define MAX_NUM_EPOCHES 1000
#define MAX_NUM_CHAR_EACH_LINE 100 // num of char in each record or line
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

// Epoch createEpoch(Sat* satArray) {
//	 Epoch epochObj;
//	 for
//	 epochObj.satArray = satArray;
//	 epochObj.numberOfSatellie = numberOfSatellie;

//	 return epochObj;
// }

bool isTimeInArray(char* time, char** timeArray, int timeArrayIndex) {
	//printf("time: %s	size: %i\n", time, strlen(time));
	for (int i = 0; i < timeArrayIndex; i++) {
		if(strcmp(time, timeArray[i]) == 0) {
			// If duplicate time entry found
			
			return true;
		}
	}
	//printf("time: %s\n", time);
	//return true;
	return false;
}

char* concat(const char* str1, const char* str2)
{
	char* str;
	str = (char*)malloc(strlen(str1) + strlen(str2) +2); 
	// +1 for the null-terminator; +1 for the seperator
	if (str == NULL) {
		fprintf(stderr, "malloc() failed in concat()\n");
	} // check for errors in malloc here
	else {
		strcpy(str, str1);
		strcat(str, " ");
		strcat(str, str2);
	}
	return str;
}

int main (void) {
	char** timeArray; // A time array storing all unique times
	Sat* satArray; // A sat array storing all satellite signals
	//Epoch* epochArray;
	
	timeArray = malloc(MAX_NUM_EPOCHES*sizeof(char*));
	if (timeArray == NULL) {
		fprintf(stderr, "malloc() failed for creating timeArray\n");
	}
	for (int i = 0; i < MAX_NUM_EPOCHES; i++) {
		timeArray[i] = (char*)malloc(NUM_CHAR_DATE + NUM_CHAR_TIME + 2);
		if (timeArray[i] == NULL) {
			fprintf(stderr, "malloc() failed for creating timeArray element\n");
		}
	}
	
	satArray = (Sat*)malloc(MAX_NUM_SIGNALS*sizeof(Sat));
	if (satArray == NULL) {
		fprintf(stderr, "malloc() failed for creating satArray\n");
	}
	
	//epochArray = (Sat*)malloc(MAX_NUM_EPOCHES*sizeof(Sat));
	
	int satArrayIndex = 0;
	int timeArrayIndex = 0;
	
	FILE *fp = NULL;
	char line[MAX_NUM_CHAR_EACH_LINE];
	
	/* opening file for reading */
	fp = fopen(INPUT_FILE_PATH, "r");
	if(fp == NULL) {
		perror("Error opening file\n");
		return(-1);
	}
	
	fgets(line, sizeof(line), fp); // Skip header
	while(fgets(line, sizeof(line), fp) != NULL) {


		char* time1 = (char*)malloc(sizeof(char)*(NUM_CHAR_DATE +1));
		if (time1 == NULL) {
			fprintf(stderr, "malloc() failed for creating time1\n");
		}
		
		char* time2 = (char*)malloc(sizeof(char)*(NUM_CHAR_TIME +1));
		if (time2 == NULL) {
			fprintf(stderr, "malloc() failed for creating time2\n");
		}
		
		char* satName = (char*)malloc(sizeof(char)*(NUM_CHAR_SAT +1));
		if (satName == NULL) {
			fprintf(stderr, "malloc() failed for creating satName\n");
		}
		
		double az;
		double el;
		double snr;

		sscanf(line, "%s %s %s %lf %lf %lf", time1, time2, satName, &az, &el, &snr);
		char* time = concat(time1, time2);
		
		if(!isTimeInArray(time, timeArray, timeArrayIndex)) { // Check if the current time is unique
			strcpy (timeArray[timeArrayIndex], time);
			//
			printf("timeArray[timeArrayIndex]: %s\n", timeArray[timeArrayIndex]);
			printf("timeArrayIndex: %d\n", timeArrayIndex);
			//
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
		
		/*
			free()
		*/
		free(time1);
		free(time2);
		free(satName);
		free(time);
	}

	/*
		free()
	*/
 	for (int i = 0; i < MAX_NUM_EPOCHES; i++) {
		free(timeArray[i]);
	}
	free(timeArray);
	free(satArray);

	fclose(fp); 
	
	exit(0);
}
