#ifndef CONFIG_H
#define CONFIG_H

/* Antenna parameters */
#define SNR_CAP 55
#define SNR_FLOOR 30
#define ANT_SNR_ADJ_MIN 42 // average min adjusted SNR
#define ANT_SNR_ADJ_MAX 50 // 50 // average max adjusted SNR
#define SAT_CUTOFF 0	   // satellite cutoff elevation angle

/* antenna truth for statistics or simulation */
#define TRUE_AZ 35
#define TRUE_EL 30 // in degrees

/* processing options */
#define CONVERGENCE_CORRECTION false // elevation corrector to determined elevation angle
#define WEIGHT_METHOD 0				 // 0: off (uniform weight for all observations)
									 // 1: satellite elevation angle
									 // 2: multipath linear combination

/* simulation mode */
#define SIMULATION false  // simulation switch
#define NUM_EPOCH 10000	  // number of epochs to simulate
#define NUM_SAT_SPHERE 24 // number of GNSS satellites globally available
#define SIMULATE_SAT_VISIBLE false
#define NUM_SAT_VISIBLE 0 // number of satellite in the epoch to simulate
#define SNR_A 10		  // quadratic mapping function SNR = -( SNR_A /8100) (off-boresight angle in degrees)^2 + SNR_C
#define SNR_C 100		  // max SNR value
#define SNR_STD_MIN 1	  // minimum snr standard deviation i.e. at 0 deg off-boresight angle
#define SNR_STD_MAX 1	  // max snr standard deviation i.e. at 90 deg off-boresight angle

/* file IO */
#define DEBUG false // full terminal output for debug
#define RMS true	// output rms (requires antenna truth)
#define INPUT_FILE_PATH_DEFAULT "input.txt"
#define OUTPUT_FILE_PATH_DEFAULT "output.txt"
#define MAX_NUM_EPOCH 86400		 // supports up to 24h data @ 1hz rate, change if your file has number of epochs more than this
#define MAX_NUM_SIGNAL_EPOCH 200 // maximum number of satellite signals visible in an epoch
#define MAX_NUM_SIGNAL (MAX_NUM_EPOCH * MAX_NUM_SIGNAL_EPOCH)
#define MAX_NUM_CHAR_LINE 100 // num of char in each record or line
#define NUM_CHAR_DATE 10
#define NUM_CHAR_TIME 10
#define NUM_CHAR_SAT 3

#endif
