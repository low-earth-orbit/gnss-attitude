#ifndef SNR_H
#define SNR_H

void adjSnr(char *prn, double *el, double *snr);
double getCosA(char *prn, double *snr);
void adjSnr2(char *prn, double *el, double *snr2);
double getCosA2(char *prn, double *snr2);

#endif
