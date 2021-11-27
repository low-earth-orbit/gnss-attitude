# gnss-ads
gnss-ads is a GNSS-based single antenna attitude determination software implementing the signal-to-noise ratio (SNR) based algorithm originally proposed by Axelrad and Behre (1999) and later improved by Wang et al. (2009) and acting as an extension of RTKLIB, an open-source code library for GNSS processing. This software takes processing results from RTKLIB and outputs the determined antenna boresight vectors in local coordinates (ENU).

The AD algorithm can be summarized as follows:
1. Before determining the attitude using observation data, a calibration data set is collected. For each constellation, using the calibration data set, a multiparameter nonlinear regression is used for obtaining SNR adjustment parameters for each satellite group and the SNR mapping function.
2. For a given epoch, the line-of-sight (LOS) vectors from satellite to receiver are calculated from satellite and receiver locations. The satellite and receiver locations are derived from the GNSS navigation and observation files. 
3. SNR values from the observation file are adjusted according to the adjustment terms developed in Step 1.
4. Using adjusted SNR values from Step 3, the off-boresight angles can be found by the SNR mapping function developed in Step 1.
5. The LOS vectors from Step 2 and off-boresight angles from Step 4 are put into a multiple linear regression to determine the antenna boresight vector.
