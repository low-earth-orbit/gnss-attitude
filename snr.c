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

void adjSnr(char *prn, double *el, double *snr)
{
	if (prn[0] == 'G') // if GPS
	{
		/* no path loss adjustment for GPS sats */
		//*snr += 20.0 * log10(sqrt(pow(26558.1, 2) - pow(6378.1, 2) * cos(pow(deg2rad(*el), 2)) - 6378.1 * sin(deg2rad(*el))));
		if (isStrInArray(prn, IIRLegacy, 5))
		{
			*snr += 0.201399479;
			*snr += (1 / 445.2543247) * pow((*el - 41.89688572), 2);
		}
		else if (isStrInArray(prn, IIRImproved, 3))
		{
			*snr += 0.201399479;
			*snr += (1 / 754.2717845) * pow((*el - 42.324671908748), 2);
		}
		else if (isStrInArray(prn, IIRM, 7))
		{
			*snr += 0.049106072;
			*snr += (1 / 754.2717845) * pow((*el - 42.324671908748), 2); // IIM are all improved antenna panel
		}
		else if (isStrInArray(prn, III, 5))
		{
			*snr += 0.903867316;
			*snr += (1 / 943.164578563126) * pow((*el - 41.1209273878527), 2);
		}
		else // the rest are IIF: no adjustment for IIF EIRP; it's used as the reference
		{
			*snr += (1 / 1058.16544398738) * pow((*el - 56.5593739795749), 2);
		}
	}
}
