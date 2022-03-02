#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "util.h"
#include "struct.h"

/* THIS FILE IS FOR UNB3 STATION 
	Selected band 1 
		GPS				SIC
		GLONASS			SIC
		BeiDou			S1X
		Galileo			S1X
	Selected band 2
		GPS				S5X
		GLONASS			S2P
		BeiDou			S5X
		Galileo			S5X
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

/* GPS S5X */
char *l5gps1[7] = {"G05", "G07", "G12", "G15", "G17", "G29", "G31"};

/* GLONASS SIC */
char *glo1[5] = {"R19", "R22", "R06", "R13", "R20"};
char *glo2[1] = {"R16"};
char *glo3[1] = {"R01"};
char *glo4[3] = {"R18", "R10", "R08"};

/* GLONASS S2P */
char *l2glo1[2] = {"R01", "R13"};
char *l2glo2[2] = {"R16", "R22"};
char *l2glo3[1] = {"R19"};

/* BeiDou S1X & S5X*/
char *bdsIGSO[3] = {"C38", "C39", "C40"};

/* Galileo S1X */
char *gal1[1] = {"E18"};
char *gal2[3] = {"E11", "E12", "E19"};

/* Galileo S5X */
char *l5gal1[2] = {"E11", "E19"};
char *l5gal2[2] = {"E14", "E18"};
char *l5gal3[1] = {"E12"};

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
	}
	else if (prn[0] == 'C') // if BDS
	{
		*snr += (1 / 755.252712796527) * pow((*el - 33.5734523918311), 2);

		if (isStrInArray(prn, bdsIGSO, 3))
		{
			*snr += 20.0 * log10(sqrt(pow((35786 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
			*snr -= 1.52720412762558;
		}
		else
		{
			*snr += 20.0 * log10(sqrt(pow((21500 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		}
	}
	else if (prn[0] == 'E') // if GAL
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((23222 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		if (isStrInArray(prn, gal1, 1))
		{
			*snr -= 2.71604252526557;
			*snr += (1 / 557.479567031803) * pow((*el - 49.6690514642272), 2);
		}
		else if (isStrInArray(prn, gal2, 3))
		{
			*snr -= -2.94713629817769;
			*snr += (1 / 475.512462961052) * pow((*el - 36.3152172021346), 2);
		}
		else
		{
			*snr += (1 / 799.106246264129) * pow((*el - 34.6062359758989), 2);
		}
	}
}

/*
	Selected band 1 signals
	based on PRN and SNR, look up off-boresight angle from the corresponding mapping function
	Ouput: cosine alpha
*/
double getCosA(char *prn, double *snr)
{
	double a, b, c;
	if (prn[0] == 'G') // if GPS
	{
		a = -0.0000582954475421757; // coefficient A in SNR mapping function SNR = A a^2 + b
		b = 2.73344499843439;
		c = 135.992224960431; // constant b in SNR mapping function
	}
	else if (prn[0] == 'R') // if GLO
	{
		a = -0.00152704703004618;
		b = 2;
		c = 139.547595133865;
	}
	else if (prn[0] == 'C') // if BDS
	{
		a = -0.00158688882127597;
		b = 2;
		c = 138.366163043981;
	}
	else if (prn[0] == 'E') // if GAL
	{
		a = -0.00159161009195648;
		b = 2;
		c = 140.498788381318;
	}
	else
	{
		fprintf(stderr, "PRN not recognized; check your input file.\n");
	}

	double alphaSq = (*snr - c) / a;
	double alpha;
	double cosA;

	if (alphaSq < 0)
	{
		alpha = 0;
	}
	else
	{
		alpha = pow(alphaSq, (1 / b));
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
	if (prn[0] == 'G') // if GPS
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((20200 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		if (isStrInArray(prn, l5gps1, 7))
		{
			/* power adjustment */
			*snr -= -4.13204675579881;
			/* off-nadir adjustment */
			*snr += (1 / 2601.38891417737) * pow((*el - 119.455650845675), 2);
		}
		else
		{
			*snr += (1 / 1356.10630933808) * pow((*el - 63.7616736614594), 2);
		}
	}
	else if (prn[0] == 'R') // if GLO
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((19140 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		/* off nadir adjustment */
		*snr += (1 / 829.990053190839) * pow((*el - 42.9605167180683), 2);

		if (isStrInArray(prn, l2glo1, 2))
		{
			*snr -= -5.11534971992524;
		}
		else if (isStrInArray(prn, l2glo2, 2))
		{
			*snr -= -12.084126134283;
		}
		else if (isStrInArray(prn, l2glo3, 1))
		{
			*snr -= -2.56182116293503;
		}
	}
	else if (prn[0] == 'C') // if BDS
	{

		if (isStrInArray(prn, bdsIGSO, 3))
		{
			*snr += 20.0 * log10(sqrt(pow((35786 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
			*snr -= 3.37609418826693;
		}
		else
		{
			*snr += 20.0 * log10(sqrt(pow((21500 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		}

		*snr += (1 / 1148.03420582938) * pow((*el - 43.2370679814993), 2);
	}
	else if (prn[0] == 'E') // if Galileo
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((23222 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		if (isStrInArray(prn, l5gal1, 2))
		{
			*snr -= -4.04061445720716;
			*snr += (1 / 1031.68361033166) * pow((*el - 29.7836305800164), 2);
		}
		else if (isStrInArray(prn, l5gal2, 2))
		{
			*snr -= 1.28670728269114;
			*snr += (1 / 1031.68361033166) * pow((*el - 29.7836305800164), 2);
		}
		else if (isStrInArray(prn, l5gal3, 1))
		{
			*snr -= -3.12308535024438;
			*snr += (1 / 1031.68361033166) * pow((*el - 29.7836305800164), 2);
		}
		else
		{
			*snr += (1 / 1358.93545020396) * pow((*el - 27.6007889933465), 2);
		}
	}
}

/*
	Selected band 2 signals
*/
double getCosA2(char *prn, double *snr)
{
	double a, b;
	if (prn[0] == 'G') // if GPS
	{
		a = -0.00131698575710085; // coefficient A in snr mapping function snr2 = A a^2 + b
		b = 141.915876909848;	  // constant b in snr mapping function
	}
	else if (prn[0] == 'R') // if GLO
	{
		a = -0.00154145910190414;
		b = 137.75957995882;
	}
	else if (prn[0] == 'C') // if BDS
	{
		a = -0.00158898922644623;
		b = 141.54284589507;
	}
	else if (prn[0] == 'E') // if GAL
	{
		a = -0.00188873929797853;
		b = 142.564060886881;
	}
	else
	{
		fprintf(stderr, "PRN not recognized; check your input file.\n");
	}

	double alphaSq = (*snr - b) / a;
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
