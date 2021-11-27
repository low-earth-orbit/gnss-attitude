# gnss-ads
gnss-ads is a GNSS-based single antenna attitude determination software implementing the SNR-based algorithm originally proposed by Axelrad and Behre (1999)[^1] and later improved by Wang et al. (2005)[^2]. Acting as an extension for RTKLIB[^3], this software takes processing results from RTKLIB and outputs the determined antenna boresight vectors in local coordinates (ENU). 

gnss-ads is being developed for future use in VIOLET, a nano-satellite by CubeSat NB[^4]. The author introduced 2 innovations:
1. The adjustment terms and the SNR mapping function are determined simultaneously, using a statistical model.
2. All 4 globally operating navigation satellite systems (GPS, BeiDoue, GLONASS and Galileo) are incoporated as well as support for multiple frequencies.

Here is the summary of designed AD algorithm:
1. Before determining the attitude using observation data, a calibration data set is collected. For each constellation, using the calibration data set, a multiparameter nonlinear regression is used for obtaining SNR adjustment parameters for each satellite group and the SNR mapping function.
3. For a given epoch, the line-of-sight (LOS) vectors from satellite to receiver are calculated from satellite and receiver locations. The satellite and receiver locations are derived from the GNSS navigation and observation files. 
4. SNR values from the observation file are adjusted according to the adjustment terms developed in Step 1.
5. Using adjusted SNR values from Step 3, the off-boresight angles can be found by the SNR mapping function developed in Step 1.
6. The LOS vectors from Step 2 and off-boresight angles from Step 4 are put into a multiple linear regression to determine the antenna boresight vector.

[^1]: https://ieeexplore.ieee.org/abstract/document/736346
[^2]: https://eprints.qut.edu.au/2938/
[^3]: http://www.rtklib.com/
[^4]: https://www.unb.ca/initiatives/cubesat/
