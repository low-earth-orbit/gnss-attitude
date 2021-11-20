#ifndef CONFIG_H
#define CONFIG_H

/* file IO */
#define INPUT_FILE_PATH_DEFAULT "input.txt"
#define OUTPUT_FILE_PATH_DEFAULT "output.txt"

/* if true azimuth and elevation angle of the antenna are known, input here for statistics (in degrees) */
#define TRUE_EL 90
#define TRUE_AZ 180

/* used only for simulation: default off */
#define SIMULATION false   // simulation on or off
#define NUM_EPOCH 10000	   // number of epochs to simulate
#define NUM_SAT_SPHERE 100 // number of GNSS satellites globally available
#define SNR_A 14.73		   // mapping function SNR = A cos(a) + c
#define SNR_C 36.33
#define SNR_STD 2
//#define SNR_STD_MAX 1.0 // SNR variation (standard deviation)
//#define SNR_STD_MIN 1.0 // default value is 1, meaning non-weighted least squares

/* normally no need to change */
#define MAX_NUM_EPOCH 86400	  // up to 24h data @ 1hz rate, change if your file has number of epochs more than this
#define MAX_NUM_SIGNAL_EPOCH 200 // maximum number of satellite signals visible in an epoch
#define MAX_NUM_SIGNAL (MAX_NUM_EPOCH * MAX_NUM_SIGNAL_EPOCH)
#define MAX_NUM_CHAR_LINE 100 // num of char in each record or line
#define NUM_CHAR_DATE 10
#define NUM_CHAR_TIME 10
#define NUM_CHAR_SAT 3

#endif
