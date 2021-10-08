#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#define ARRAY_SIZE 8192

typedef struct {
    char* time;
    char* sat;
    double az;
    double el;
    double snr;
    double l1_mp;
}Sat;

char* concat(const char *time1, const char *time2)
{
    char *time = malloc(strlen(time1) + strlen(time2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(time, time1);
    strcat(time, "\t");
    strcat(time, time2);
    return time;
}

Sat createSat(char* time1, char* time2, char* sat, double az, double el, double snr, double l1_mp) {
    Sat satObj;
    satObj.time = concat(time1, time2);
    satObj.sat = sat;
    satObj.az = az;
    satObj.el = el;
    satObj.snr = snr;
    satObj.l1_mp = l1_mp;

    return satObj;
}

typedef struct {
    Sat* satArray; // Time are the same
    int numberOfSatellie;
    // Sol* solutionArray;
}Epoch;


// Epoch createEpoch(Sat* satArray) {
//     Epoch epochObj;
//     for
//     epochObj.satArray = satArray;
//     epochObj.numberOfSatellie = numberOfSatellie;

//     return epochObj;
// }

//36x36 2d array

char* time1;
char* time2;
char* sat;
double az;
double el;
double snr;
double l1_mp;

char** uniqueTime; // A time array storing all unique times
Sat* SatArray;
Epoch* epochArray;

bool checkDuplicateTime(char* time, Sat* satArray, int satArrayIndex) {
    for (int i = 0; i < satArrayIndex; i++) {
        if(strcmp(satArray[i].time, time) == 0) {
            // If duplicate time entry found
            //printf("%s\n", satArray[i].sat);
            return true;
        }
    }
    //printf("time: %s\n", time);
    return false;
}

int main (void) {
   FILE *fp = NULL;
   char line[1024];
   
   uniqueTime = (char**)malloc(sizeof(char*));
   SatArray = (Sat*)malloc(5000*sizeof(Sat));
   //epochArray = (Sat*)malloc(5000*sizeof(Sat));
   int SatArrayIndex = 0;
   int TimeArrayIndex = 0;

   /* opening file for reading */
   fp = fopen("sat_pos_az_el.txt" , "r");
   if(fp == NULL) {
      perror("Error opening file");
      return(-1);
   }

   fgets(line, sizeof(line), fp); // Skip header
   while(fgets(line, sizeof(line), fp) != NULL) {
        
        time1 = (char*)malloc(sizeof(char*));
        time2 = (char*)malloc(sizeof(char*));
        sat = (char*)malloc(sizeof(char)*20);

        sscanf(line, "%s %s %s %lf %lf %lf %lf", time1, time2, sat, &az, &el, &snr, &l1_mp);

        Sat satObj = createSat(time1, time2, sat, az, el, snr, l1_mp);
        SatArray[SatArrayIndex] = satObj; // Add each sat to sat array
        
        
        if(!checkDuplicateTime(concat(time1, time2), SatArray, SatArrayIndex)) { // Check if the current time is unique
            uniqueTime[TimeArrayIndex] = concat(time1, time2);
            printf("uniqueTime[TimeArrayIndex]: %s\n", uniqueTime[TimeArrayIndex]);
            printf("SatArray[SatArrayIndex].sat: %s\n", SatArray[SatArrayIndex].sat);
            printf("TimeArrayIndex: %d\n", TimeArrayIndex);
            TimeArrayIndex++;
        }  
        SatArrayIndex++;
    }

   fclose(fp);
   
   exit(0);
}
