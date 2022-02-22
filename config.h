#ifndef CONFIG_H
#define CONFIG_H

/* Antenna parameters */
#define ANT_SNR_STD_MIN 1.27 // 1.02 //(1.02 + 1.52) / 2
#define ANT_SNR_STD_MAX 4.14 // 4.53 //(4.53 + 3.75) / 2
#define ANT_SNR_RAW_MIN 10	 // absolute min measured SNR
#define ANT_SNR_RAW_MAX 60	 // absolute maximum measured SNR
#define ANT_SNR_ADJ_MIN 42	 // average min adjusted SNR
#define ANT_SNR_ADJ_MAX 50	 // 50 // average max adjusted SNR
#define SAT_CUTOFF 0		 // satellite cutoff elevation angle

/* antenna truth for statistics or simulation */
#define TRUE_EL -60 // in degrees
#define TRUE_AZ -180

/* processing options */
#define CONVERGENCE_CORRECTION true // elevation corrector to determined elevation angle
#define OUTLIER_FACTOR 3			// sigma factor SNR outlier factor

/* simulation mode */
#define SIMULATION false // simulation switch
#define NUM_EPOCH 10000	 // number of epochs to simulate
#define SIMULATE_SAT_VISIBLE true
#define NUM_SAT_SPHERE 100 // number of GNSS satellites globally available
#define NUM_SAT_VISIBLE 30 // number of satellite in the epoch to simulate
#define SNR_A 10		   // quadratic mapping function SNR = -( SNR_A /8100) (off-boresight angle in degrees)^2 + SNR_C
#define SNR_C 100		   // max SNR value
#define SNR_STD_MIN 0.5	   // minimum snr standard deviation i.e. at 0 deg off-boresight angle
#define SNR_STD_MAX 0.5	   // max snr standard deviation i.e. at 90 deg off-boresight angle
#define SKEW false		   // simulate multipath

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
