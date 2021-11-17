#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "util.h"
#include "struct.h"

/* GPS groups */
char *IIRLegacy[5] = {"G13", "G20", "G28", "G16", "G21"};
char *IIRImproved[3] = {"G22", "G19", "G02"};
char *IIRM[7] = {"G17", "G31", "G12", "G15", "G29", "G07", "G05"};
char *III[5] = {"G04", "G11", "G14", "G18", "G23"};

/* GLONASS groups */
char *glo1[5] = {"R19", "R22", "R06", "R13", "R20"};
char *glo2[1] = {"R16"};
char *glo3[1] = {"R01"};
char *glo4[3] = {"R18", "R10", "R08"};

/* BeiDou */
char *bdsIGSO[3] = {"C38", "C39", "C40"};

/* Galileo */
char *gal1[1] = {"E18"};
char *gal2[3] = {"E11", "E12", "E19"};

/*
	L1 signals
*/
void adjSnr(char *prn, double *el, double *snr)
{
	if (prn[0] == 'G') // if GPS
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((20200 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		if (isStrInArray(prn, IIRLegacy, 5))
		{
			/* power adjustment */
			*snr += 0.217500466451852;
			/* off-nadir adjustment */
			*snr += (1 / 363.421120729045) * pow((*el - 42.5153301150902), 2);
		}
		else if (isStrInArray(prn, IIRImproved, 3))
		{
			*snr += 0.217500466451852;
			*snr += (1 / 621.389380226501) * pow((*el - 43.7356043887304), 2);
		}
		else if (isStrInArray(prn, IIRM, 7))
		{
			*snr += 0.124778481882819;
			*snr += (1 / 621.389380226501) * pow((*el - 43.7356043887304), 2); // IIM are all improved antenna panel
		}
		else if (isStrInArray(prn, III, 5))
		{
			*snr += 0.787531267514597;
			*snr += (1 / 642.782267958422) * pow((*el - 42.6823843084975), 2);
		}
		else // the rest are IIF: no adjustment for IIF EIRP; it's used as the reference
		{
			*snr += (1 / 875.058989953713) * pow((*el - 60.4752792046367), 2);
		}
	}
	else if (prn[0] == 'R') // if GLO
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((19140 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		/* off nadir adjustment */
		*snr += (1 / 521.007362094454) * pow((*el - 36.550985427478), 2);

		if (isStrInArray(prn, glo1, 5))
		{
			*snr += 10.9064243930786;
		}
		else if (isStrInArray(prn, glo2, 1))
		{
			*snr += 4.96681858814593;
		}
		else if (isStrInArray(prn, glo3, 1))
		{
			*snr += 3.25537068647712;
		}
		else if (isStrInArray(prn, glo4, 3))
		{
			*snr += 1.55168919985381;
		}
		// no power adjustment for 03,07,02,17,14,15,21,11,09,24,04,12,05,23
	}
	else if (prn[0] == 'C') // if BDS
	{
		*snr += (1 / 691.303223036127) * pow((*el - 42.072884954551), 2);

		if (isStrInArray(prn, bdsIGSO, 3))
		{
			*snr += 20.0 * log10(sqrt(pow((35786 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		}
		else
		{
			*snr += 20.0 * log10(sqrt(pow((21500 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		}
	}
	else if (prn[0] == 'E') // if Galileo
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((23222 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		if (isStrInArray(prn, gal1, 1))
		{
			*snr -= 2.71468820075693;
			*snr += (1 / 557.983981371012) * pow((*el - 49.6691259258094), 2);
		}
		else if (isStrInArray(prn, gal2, 3))
		{
			*snr -= -2.9480183968319;
			*snr += (1 / 475.995769065121) * pow((*el - 36.3059407271915), 2);
		}
		else
		{
			*snr += (1 / 818.653254377219) * pow((*el - 32.7874643347818), 2);
		}
	}
}

/*
	based on PRN and SNR, look up off-boresight angle from the corresponding mapping function
	Ouput: cosine alpha
*/
double getCosA(char *prn, double *snr)
{
	double a, b;
	if (prn[0] == 'G') // if GPS
	{
		a = -0.00126170488051369; // coefficient A in SNR mapping function SNR = A a^2 + b
		b = 137.103669190613;	  // constant b in SNR mapping function
	}
	else if (prn[0] == 'R') // if GLO
	{
		a = -0.0014250672233639;
		b = 139.250622304925;
	}
	else if (prn[0] == 'C') // if BDS
	{
		a = -0.00133450520582277;
		b = 137.639888386765;
	}
	else if (prn[0] == 'E') // if Galileo
	{
		a = -0.00159161009195648;
		b = 140.498788381318;
	}
	else
	{
		fprintf(stderr, "PRN not recognized; check your input file.\n");
	}

	double alphaSq = (*snr - b) / a;
	if (alphaSq < 0)
	{
		alphaSq = 0;
	}
	double alpha = sqrt(alphaSq);
	if (alpha > 90)
	{
		alpha = 90;
	}
	else if (alpha < 0)
	{
		alpha = 0;
	}
	double cosA = cos(deg2rad(alpha));
	return cosA;
}

/*
	L2 signals
*/
void adjSnr2(char *prn, double *el, double *snr2)
{
	if (prn[0] == 'G') // if GPS
	{
		/* path loss adjustment */
		*snr2 += 20.0 * log10(sqrt(pow((20200 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		if (isStrInArray(prn, IIRLegacy, 5))
		{
			/* power adjustment */
			*snr2 += 0.217500466451852;
			/* off-nadir adjustment */
			*snr2 += (1 / 363.421120729045) * pow((*el - 42.5153301150902), 2);
		}
		else if (isStrInArray(prn, IIRImproved, 3))
		{
			*snr2 += 0.217500466451852;
			*snr2 += (1 / 621.389380226501) * pow((*el - 43.7356043887304), 2);
		}
		else if (isStrInArray(prn, IIRM, 7))
		{
			*snr2 += 0.124778481882819;
			*snr2 += (1 / 621.389380226501) * pow((*el - 43.7356043887304), 2); // IIM are all improved antenna panel
		}
		else if (isStrInArray(prn, III, 5))
		{
			*snr2 += 0.787531267514597;
			*snr2 += (1 / 642.782267958422) * pow((*el - 42.6823843084975), 2);
		}
		else // the rest are IIF: no adjustment for IIF EIRP; it's used as the reference
		{
			*snr2 += (1 / 875.058989953713) * pow((*el - 60.4752792046367), 2);
		}
	}
	else if (prn[0] == 'R') // if GLO
	{
		/* path loss adjustment */
		*snr2 += 20.0 * log10(sqrt(pow((19140 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		/* off nadir adjustment */
		*snr2 += (1 / 521.007362094454) * pow((*el - 36.550985427478), 2);

		if (isStrInArray(prn, glo1, 5))
		{
			*snr2 += 10.9064243930786;
		}
		else if (isStrInArray(prn, glo2, 1))
		{
			*snr2 += 4.96681858814593;
		}
		else if (isStrInArray(prn, glo3, 1))
		{
			*snr2 += 3.25537068647712;
		}
		else if (isStrInArray(prn, glo4, 3))
		{
			*snr2 += 1.55168919985381;
		}
		// no power adjustment for 03,07,02,17,14,15,21,11,09,24,04,12,05,23
	}
	else if (prn[0] == 'C') // if BDS
	{
		*snr2 += (1 / 691.303223036127) * pow((*el - 42.072884954551), 2);

		if (isStrInArray(prn, bdsIGSO, 3))
		{
			*snr2 += 20.0 * log10(sqrt(pow((35786 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		}
		else
		{
			*snr2 += 20.0 * log10(sqrt(pow((21500 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		}
	}
	else if (prn[0] == 'E') // if Galileo
	{
		/* path loss adjustment */
		*snr2 += 20.0 * log10(sqrt(pow((23222 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		if (isStrInArray(prn, gal1, 1))
		{
			*snr2 -= 2.71468820075693;
			*snr2 += (1 / 557.983981371012) * pow((*el - 49.6691259258094), 2);
		}
		else if (isStrInArray(prn, gal2, 3))
		{
			*snr2 -= -2.9480183968319;
			*snr2 += (1 / 475.995769065121) * pow((*el - 36.3059407271915), 2);
		}
		else
		{
			*snr2 += (1 / 818.653254377219) * pow((*el - 32.7874643347818), 2);
		}
	}
}

double getCosA2(char *prn, double *snr2)
{
	double a, b;
	if (prn[0] == 'G') // if GPS
	{
		a = -0.00126170488051369; // coefficient A in snr2 mapping function snr2 = A a^2 + b
		b = 137.103669190613;	  // constant b in snr2 mapping function
	}
	else if (prn[0] == 'R') // if GLO
	{
		a = -0.0014250672233639;
		b = 139.250622304925;
	}
	else if (prn[0] == 'C') // if BDS
	{
		a = -0.00133450520582277;
		b = 137.639888386765;
	}
	else if (prn[0] == 'E') // if Galileo
	{
		a = -0.00159161009195648;
		b = 140.498788381318;
	}
	else
	{
		fprintf(stderr, "PRN not recognized; check your input file.\n");
	}

	double alphaSq = (*snr2 - b) / a;
	if (alphaSq < 0)
	{
		alphaSq = 0;
	}
	double alpha = sqrt(alphaSq);
	if (alpha > 90)
	{
		alpha = 90;
	}
	else if (alpha < 0)
	{
		alpha = 0;
	}
	double cosA = cos(deg2rad(alpha));
	return cosA;
}
