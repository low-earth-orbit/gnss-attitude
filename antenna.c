#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#define MAX_NUM_SIGNALS 10000
#define MAX_NUM_EPOCHES 1000
#define LEN_EACH_RECORD 1000

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
//     Epoch epochObj;
//     for
//     epochObj.satArray = satArray;
//     epochObj.numberOfSatellie = numberOfSatellie;

//     return epochObj;
// }

/*
bool isTimeInSatArray(char* time, Sat* satArray, int satArrayIndex) {
    for (int i = 0; i < satArrayIndex; i++) {
        if(strcmp(satArray[i].time, time) == 0) {
            // If duplicate time entry found
            //printf("%s\n", satArray[i].sat);
            return true;
        }
    //printf("time: %s\n", time);
    //return false;
	}
}
*/

bool isTimeInArray(char* time, char** uniqueTimeArray, int timeArrayIndex) {
    for (int i = 0; i < timeArrayIndex; i++) {
        if(strcmp(time, uniqueTimeArray[i]) == 0) {
            // If duplicate time entry found
            return true;
        }
    }
    //printf("time: %s\n", time);
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
	char** uniqueTimeArray; // A time array storing all unique times
	Sat* satArray; // A sat array storing all satelite signals
	//Epoch* epochArray;
	
	uniqueTimeArray = (char**)malloc(MAX_NUM_EPOCHES*sizeof(char*));
	satArray = (Sat*)malloc(MAX_NUM_SIGNALS*sizeof(Sat));
	//epochArray = (Sat*)malloc(MAX_NUM_EPOCHES*sizeof(Sat));
	
	int satArrayIndex = 0;
	int timeArrayIndex = 0;
	
	FILE *fp = NULL;
	char line[LEN_EACH_RECORD];
	
	/* opening file for reading */
	fp = fopen("new.txt" , "r");
	if(fp == NULL) {
		perror("Error opening file\n");
		return(-1);
	}
	
	fgets(line, sizeof(line), fp); // Skip header
	while(fgets(line, sizeof(line), fp) != NULL) {

		char* time1 = (char*)malloc(sizeof(char)*11);
		char* time2 = (char*)malloc(sizeof(char)*11);
		char* satName = (char*)malloc(sizeof(char)*4);
		double az;
		double el;
		double snr;

		sscanf(line, "%s %s %s %lf %lf %lf", time1, time2, satName, &az, &el, &snr);
		char* time = concat(time1, time2);
		//printf("time: %s\n", time);
		
		/*
 		if(!isTimeInSatArray(time, satArray, satArrayIndex)) { // Check if the current time is unique
			uniqueTimeArray[timeArrayIndex] = time;
			printf("uniqueTimeArray[timeArrayIndex]: %s\n", uniqueTimeArray[timeArrayIndex]);
			printf("timeArrayIndex: %d\n", timeArrayIndex);
			timeArrayIndex++;
        }
		*/
  		if(!isTimeInArray(time, uniqueTimeArray, timeArrayIndex)) { // Check if the current time is unique
			uniqueTimeArray[timeArrayIndex] = time;
			printf("uniqueTimeArray[timeArrayIndex]: %s\n", uniqueTimeArray[timeArrayIndex]);
			printf("timeArrayIndex: %d\n", timeArrayIndex);
			timeArrayIndex++;
        }
		
		Sat satObj = createSat(time, satName, az, el, snr);
		satArray[satArrayIndex] = satObj; // Add each sat to sat array
		//printf("satArrayIndex: %d\n", satArrayIndex);
		//printf("satArray[satArrayIndex].time: %s\n", satArray[satArrayIndex].time);
        satArrayIndex++;
        
        free(time1);
        free(time2);
        free(satName);
		free(time);
    }
    
	free(uniqueTimeArray);
	free(satArray);
	fclose(fp); 
	exit(0);
}
