# gnss-attitude
## About
gnss-attitude is a GNSS-based single antenna attitude determination software implementing the SNR-based algorithm proposed initially by [Axelrad and Behre (1999)](https://ieeexplore.ieee.org/abstract/document/736346). Acting as an extension for [RTKLIB](http://www.rtklib.com/), this software takes processing results from RTKLIB. It outputs the determined antenna boresight vectors in local coordinates (ENU).

gnss-attitude is being developed for future use in VIOLET, a nano-satellite by [CubeSat NB](https://www.unb.ca/initiatives/cubesat/). The author introduced a few innovations:
1. The adjustment terms and the SNR mapping function are determined simultaneously, using a statistical model.
2. All four globally operating navigation satellite systems (GPS, BeiDou, GLONASS and Galileo) are incorporated.
3. Support for multiple frequencies.

## Under Development
The software is under development. No releases have been published.

## Designed AD algorithm
1. Before determining the attitude using observation data, a calibration data set is collected. A multiparameter nonlinear regression obtains SNR adjustment parameters for each satellite group and the SNR mapping function.
2. For a given epoch, the line-of-sight (LOS) vectors from satellite to receiver are calculated from satellite and receiver locations. The satellite and receiver locations are derived from the GNSS navigation and observation files.
3. SNR values from the observation file are adjusted according to the adjustment terms developed in Step 1.
4. The off-boresight angles can be found by the SNR mapping function developed in Step 1 and the adjusted SNR values from Step 3.
5. The LOS vectors from Step 2 and off-boresight angles from Step 4 are put into a multiple linear regression to determine the antenna boresight vector.

## Performance
* Using a general SNR mapping function, the system delivers an accuracy of 5° - 20° (RMS).
* Better performance can be achieved if calibration is performed.

## Prerequisites
1. Linux environment
2. [RTKLIB](http://www.rtklib.com/)
3. GNU Scientific Library

        sudo apt-get install libgsl-dev

4. GNU C Compiler

        sudo apt install build-essential

## License
Distributed under the GNU General Public License v3.0. See `LICENSE.txt` for more information.
