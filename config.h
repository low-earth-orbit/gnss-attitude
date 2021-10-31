#ifndef CONFIG_H
#define CONFIG_H

/* antenna truth */
#define TRUE_EL 90
#define TRUE_AZ 180 //If true azimuth and elevation angle of the antenna are known, input here for statistics (in degrees)

/* file I/O */
#define INPUT_FILE_PATH "input.txt"

/* mapping function */
#define SNR_A 13.36 // amplitude A in SNR mapping function
#define SNR_C 36.98 // constant c in SNR mapping function
#define SNR_STD 2.0 // SNR variation (standard deviation)

/* for simulation only */
#define NUM_EPOCH 86400 // number of epochs to simulate

/* usually no need to change*/
#define MAX_NUM_EPOCH 86400	  // up to 24h data @ 1hz rate
#define MAX_NUM_SAT_EPOCH 100 // maximum number of satellites visible in an epoch
#define MAX_NUM_SIGNAL (MAX_NUM_EPOCH * MAX_NUM_SAT_EPOCH)
#define MAX_NUM_CHAR_LINE 100 // num of char in each record or line
#define NUM_CHAR_DATE 10
#define NUM_CHAR_TIME 10
#define NUM_CHAR_SAT 3

#endif
