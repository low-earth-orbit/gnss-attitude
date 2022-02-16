#ifndef CONFIG_H
#define CONFIG_H

/* file IO */
#define INPUT_FILE_PATH_DEFAULT "input.txt"
#define OUTPUT_FILE_PATH_DEFAULT "output.txt"

/* if true azimuth and elevation angle of the antenna are known, input here for statistics (in degrees) */
#define TRUE_EL 30
#define TRUE_AZ 180

/* processing options */
#define CONVERGENCE_CORRECTION false // elevation corrector to determined elevation angle

/* simulation mode: false (turned off) by default */
#define SIMULATION true	  // simulation switch
#define NUM_EPOCH 10000	  // number of epochs to simulate
#define NUM_SAT_SPHERE 48 // number of GNSS satellites globally available
#define SNR_A 10		  // quadratic mapping function SNR = -( SNR_A /8100) (off-boresight angle in degrees)^2 + SNR_C
#define SNR_C 50		  // Max SNR value
#define SNR_STD_MIN 1	  // minimum snr standard deviation i.e. at 0 deg off-boresight angle
#define SNR_STD_MAX 4.5	  // max snr standard deviation i.e. at 90 deg off-boresight angle
#define SKEWNESS false	  // simulate multipath effect dragging down SNR values at lower elev sat

/* file and array sizes */
#define MAX_NUM_EPOCH 86400		 // supports up to 24h data @ 1hz rate, change if your file has number of epochs more than this
#define MAX_NUM_SIGNAL_EPOCH 200 // maximum number of satellite signals visible in an epoch
#define MAX_NUM_SIGNAL (MAX_NUM_EPOCH * MAX_NUM_SIGNAL_EPOCH)
#define MAX_NUM_CHAR_LINE 100 // num of char in each record or line
#define NUM_CHAR_DATE 10
#define NUM_CHAR_TIME 10
#define NUM_CHAR_SAT 3

#endif
