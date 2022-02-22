#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "util.h"
#include "struct.h"
#include "config.h" // configuration

/* THIS FILE IS FOR Ublox antenna
	Selected band 1
		GPS				SIC
		GLONASS			SIC

	No band 2
	No band 3
*/

/* GPS
char *IIRLegacy[5] = {"G13", "G16", "G20", "G21", "G28"};
char *IIRImproved[3] = {"G22", "G19", "G02"};
char *IIRM[7] = {"G05", "G07", "G12", "G15", "G17", "G29", "G31"};
char *III[5] = {"G04", "G11", "G14", "G18", "G23"};
*/

/* GPS SIC */
char *gps1[4] = {"G13", "G16", "G20", "G21"};
char *gps2[2] = {"G10", "G32"};

/* GPS S5X
char *l5gps1[7] = {"G05", "G07", "G12", "G15", "G17", "G29", "G31"};
*/

/* GLONASS SIC */
char *glo1[5] = {"R19", "R22", "R06", "R13", "R20"};
char *glo2[1] = {"R16"};
char *glo3[1] = {"R01"};
char *glo4[3] = {"R18", "R10", "R08"};

/* GLONASS S2P
char *l2glo1[2] = {"R01", "R13"};
char *l2glo2[2] = {"R16", "R22"};
char *l2glo3[1] = {"R19"};
*/

/* BeiDou S1X & S5X
char *bdsIGSO[3] = {"C38", "C39", "C40"};
*/

/* Galileo S1X
char *gal1[1] = {"E18"};
char *gal2[3] = {"E11", "E12", "E19"};
*/

/* Galileo S5X
char *l5gal1[2] = {"E11", "E19"};
char *l5gal2[2] = {"E14", "E18"};
char *l5gal3[1] = {"E12"};
*/

/*
	Selected band 1 signals
	Makes adjustment to measured SNR
*/
void adjSnr(char *prn, double *el, double *snr)
{
	if (prn[0] == 'G') // if GPS
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((20200 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		if (isStrInArray(prn, gps1, 4))
		{
			/* power adjustment */
			*snr -= 0.24852350057561;
			/* off-nadir adjustment */
			*snr += (1 / 641.212390401401) * pow((*el - 35.658512773265), 2);
		}
		else if (isStrInArray(prn, gps2, 2))
		{
			*snr -= 2.79209176750643;
			*snr += (1 / 21773.5082364001) * pow((*el - 194.183097860413), 2);
		}
		else
		{
			*snr += (1 / 2580.72026229872) * pow((*el - 34.5483935101112), 2);
		}
		*snr -= 137.57710845457 - ANT_SNR_ADJ_MAX; // normalize so max is ANT_SNR_ADJ_MAX
	}
	else if (prn[0] == 'R') // if GLO
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((19140 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		/* off nadir adjustment */
		*snr += (1 / 568.800755288789) * pow((*el - 33.1297078304606), 2);

		if (isStrInArray(prn, glo1, 5))
		{
			*snr -= -10.91272021662;
		}
		else if (isStrInArray(prn, glo2, 1))
		{
			*snr -= -4.99255336287507;
		}
		else if (isStrInArray(prn, glo3, 1))
		{
			*snr -= -3.26493391767487;
		}
		else if (isStrInArray(prn, glo4, 3))
		{
			*snr -= -1.56287012522724;
		}
		// no power adjustment for 03,07,02,17,14,15,21,11,09,24,04,12,05,23

		*snr -= 139.571690147127 - ANT_SNR_ADJ_MAX; // normalize so max is ANT_SNR_ADJ_MAX
	}
}

/*
	Selected band 1 signals
	based on PRN and SNR, look up off-boresight angle from the corresponding mapping function
	Ouput: cosine alpha
*/
double getCosA(char *prn, double *snr)
{
	double a, c;
	c = ANT_SNR_ADJ_MAX;
	if (prn[0] == 'G') // if GPS
	{
		a = -(ANT_SNR_ADJ_MAX - ANT_SNR_ADJ_MIN) / 8100.0;
		// a = -0.00110280101958772; // coefficient A in SNR mapping function SNR = A a^2 + c
	}
	else if (prn[0] == 'R') // if GLO
	{
		a = -(ANT_SNR_ADJ_MAX - ANT_SNR_ADJ_MIN) / 8100.0;
		// a = -0.000870181419750473; // coefficient A in SNR mapping function SNR = A a^2 + c
	}
	else
	{
		fprintf(stderr, "PRN not recognized; check your input file.\n");
	}

	// if (*snr > ANT_SNR_ADJ_MAX)
	// {
	// 	*snr = ANT_SNR_ADJ_MAX - (ANT_SNR_STD_MIN * OUTLIER_FACTOR - (*snr - c)); // exceeds less, substract more
	// }
	// else if (*snr < ANT_SNR_ADJ_MIN)
	// {
	// 	*snr = ANT_SNR_ADJ_MIN + (ANT_SNR_STD_MAX * OUTLIER_FACTOR - (*snr - c)); // falls behind less, add more
	// }

	double alphaSq = (*snr - c) / a;
	double alpha;
	double cosA;

	if (alphaSq < 0)
	{
		alpha = 0;
	}
	else
	{
		alpha = sqrt(alphaSq);
	}

	if (alpha > 90)
	{
		alpha = 90;
	}

	cosA = cos(deg2rad(alpha));
	return cosA;
}

/*
	Selected band 2 signals
*/
void adjSnr2(char *prn, double *el, double *snr)
{
	printf("No band 3 signal configured in SNR module.\n");
}

/*
	Selected band 2 signals
*/
double getCosA2(char *prn, double *snr)
{
	printf("No band 3 signal configured in SNR module.\n");
	return -1.0;
}

/*
	Selected band 3 signals
*/
void adjSnr3(char *prn, double *el, double *snr)
{
	printf("No band 3 signal configured in SNR module.\n");
}

/*
	Selected band 3 signals
*/
double getCosA3(char *prn, double *snr)
{
	printf("No band 3 signal configured in SNR module.\n");
	return -1.0;
}
