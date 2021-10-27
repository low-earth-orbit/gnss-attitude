#ifndef CONFIG_H
#define CONFIG_H

/* Configuration */
#define TRUE_EL 90
#define TRUE_AZ 180	  //If true azimuth and elevation angle of the antenna are known, input here for statistics (in degrees)
#define NUM_EPOCH 100 // number of satellites to simulate
#define INPUT_FILE_PATH "input.txt"

/* Usually no need to change*/
#define MAX_NUM_EPOCH 86400	  // up to 24h data @ 1hz rate
#define MAX_NUM_SAT_EPOCH 100 // maximum number of satellites visible in an epoch
#define MAX_NUM_SIGNAL (MAX_NUM_EPOCH * MAX_NUM_SAT_EPOCH)
#define MAX_NUM_CHAR_LINE 100 // num of char in each record or line
#define NUM_CHAR_DATE 10
#define NUM_CHAR_TIME 10
#define NUM_CHAR_SAT 3

#endif
