# gnss-ads
gnss-ads is a GNSS-based single antenna attitude determination software implementing the SNR-based algorithm proposed initially by Axelrad and Behre (1999)[^1] and later improved by Wang et al. (2005)[^2]. Acting as an extension for RTKLIB[^3], this software takes processing results from RTKLIB and outputs the determined antenna boresight vectors in local coordinates (ENU). 

gnss-ads is being developed for future use in VIOLET, a nano-satellite by CubeSat NB[^4]. The author introduced a few innovations:
1. The adjustment terms and the SNR mapping function are determined simultaneously, using a statistical model.
2. All four globally operating navigation satellite systems (GPS, BeiDou, GLONASS and Galileo) are incorporated.
3. Support for multiple frequencies.

Here is the summary of the designed AD algorithm:
1. Before determining the attitude using observation data, a calibration data set is collected. Using the calibration data set for each constellation, a multiparameter nonlinear regression is used to obtain SNR adjustment parameters for each satellite group and the SNR mapping function.
3. For a given epoch, the line-of-sight (LOS) vectors from satellite to receiver are calculated from satellite and receiver locations. The satellite and receiver locations are derived from the GNSS navigation and observation files. 
4. SNR values from the observation file are adjusted according to the adjustment terms developed in Step 1.
5. The off-boresight angles can be found by the SNR mapping function developed in Step 1 plus the adjusted SNR values from Step 3.
6. The LOS vectors from Step 2 and off-boresight angles from Step 4 are put into a multiple linear regression to determine the antenna boresight vector.

[^1]: https://ieeexplore.ieee.org/abstract/document/736346
[^2]: https://eprints.qut.edu.au/2938/
[^3]: http://www.rtklib.com/
[^4]: https://www.unb.ca/initiatives/cubesat/
