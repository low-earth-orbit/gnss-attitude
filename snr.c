#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "util.h"
#include "struct.h"

/* GPS groups */
/*
char *IIRLegacy[5] = {"G13", "G20", "G28", "G16", "G21"};
char *IIRImproved[3] = {"G22", "G19", "G02"};
char *IIRM[7] = {"G17", "G31", "G12", "G15", "G29", "G07", "G05"};
char *III[5] = {"G04", "G11", "G14", "G18", "G23"};
*/
char *gps1[4] = {"G13", "G20", "G16", "G21"};
char *gps2[2] = {"G10", "G32"};

/* GLONASS groups */
char *glo1[5] = {"R19", "R22", "R06", "R13", "R20"};
char *glo2[1] = {"R16"};
char *glo3[1] = {"R01"};
char *glo4[3] = {"R18", "R10", "R08"};
char *l2glo1[2] = {"R01", "R13"}; // L2
char *l2glo2[2] = {"R16", "R22"};
char *l2glo3[1] = {"R19"};

/* BeiDou */
char *bdsIGSO[3] = {"C38", "C39", "C40"};
char *l2bds1[8] = {"C06", "C08", "C09", "C11", "C12", "C13", "C14", "C16"}; // L2

/* Galileo */
char *gal1[1] = {"E18"};
char *gal2[3] = {"E11", "E12", "E19"};

void adjSnr(char *prn, double *el, double *snr)
{
	if (prn[0] == 'G') // if GPS
	{
		/* path loss adjustment */
		*snr += 20.0 * log10(sqrt(pow((20200 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		if (isStrInArray(prn, gps1, 4))
		{
			/* power adjustment */
			*snr -= 0.197857694272294;
			/* off-nadir adjustment */
			*snr += (1 / 413.926109730221) * pow((*el - 38.1146746000328), 2);
		}
		else if (isStrInArray(prn, gps2, 2))
		{
			*snr -= 1.78065087076735;
			*snr += (1 / 1080.8849341692) * pow((*el - 50.5096233934499), 2);
		}
		else
		{
			*snr += (1 / 766.360450955615) * pow((*el - 40.0844218330287), 2);
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
		*snr += (1 / 691.189635777277) * pow((*el - 42.0691729658449), 2);

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
			*snr += (1 / 817.933673844063) * pow((*el - 32.8019067709832), 2);
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
		a = -0.0014824959162717; // coefficient A in SNR mapping function SNR = A a^2 + b
		b = 137.134400410517;	 // constant b in SNR mapping function
	}
	else if (prn[0] == 'R') // if GLO
	{
		a = -0.00152704703004618;
		b = 139.547595133865;
	}
	else if (prn[0] == 'C') // if BDS
	{
		a = -0.0013345837922313;
		b = 137.640254265316;
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
	Due to a bug in RTKLIB no GAL L2 signals in the input file
*/
void adjSnr2(char *prn, double *el, double *snr2)
{
	if (prn[0] == 'G') // if GPS
	{

		/* path loss adjustment */
		*snr2 += 20.0 * log10(sqrt(pow((20200 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		if (isStrInArray(prn, gps1, 4))
		{
			/* power adjustment */
			*snr2 -= -0.217500466451852;
			/* off-nadir adjustment */
			*snr2 += (1 / 363.421120729045) * pow((*el - 42.5153301150902), 2);
		}
		else if (isStrInArray(prn, gps2, 2))
		{
			*snr2 -= -0.217500466451852;
			*snr2 += (1 / 621.389380226501) * pow((*el - 43.7356043887304), 2);
		}
		else
		{
			*snr2 += (1 / 875.058989953713) * pow((*el - 60.4752792046367), 2);
		}
	}
	else if (prn[0] == 'R') // if GLO
	{
		/* path loss adjustment */
		*snr2 += 20.0 * log10(sqrt(pow((19140 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));

		/* off nadir adjustment */
		*snr2 += (1 / 829.990053190839) * pow((*el - 42.9605167180683), 2);

		if (isStrInArray(prn, l2glo1, 2))
		{
			*snr2 -= -5.11534971992524;
		}
		else if (isStrInArray(prn, l2glo2, 2))
		{
			*snr2 -= -12.084126134283;
		}
		else if (isStrInArray(prn, l2glo3, 1))
		{
			*snr2 -= -2.56182116293503;
		}
	}
	else if (prn[0] == 'C') // if BDS
	{

		if (isStrInArray(prn, bdsIGSO, 3))
		{
			*snr2 += 20.0 * log10(sqrt(pow((35786 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		}
		else
		{
			*snr2 += 20.0 * log10(sqrt(pow((21500 + 6370), 2) - pow(6370, 2) * cos(pow(deg2rad(*el), 2)) - 6370 * sin(deg2rad(*el))));
		}

		if (isStrInArray(prn, l2bds1, 8))
		{
			*snr2 -= -1.75287484083274;
			*snr2 += (1 / 296.58285396992) * pow((*el - 36.4468079717406), 2);
		}
		else
		{
			*snr2 += (1 / 824.031034556182) * pow((*el - 35.0932731272401), 2);
		}
	}
	//Due to a bug in RTKLIB no GAL L2 signals in the input file
}

/*
	L2 Signals
*/
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
		a = -0.00154145910190414;
		b = 137.75957995882;
	}
	else if (prn[0] == 'C') // if BDS L2
	{
		a = -0.0015152861764974;
		b = 139.39109924668;
	}
	//Due to a bug in RTKLIB no GAL L2 signals in the input file
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
