# gnss-attitude
## Introduction
gnss-attitude is a GNSS-based single antenna attitude determination software implementing the SNR-based algorithm first proposed by [Axelrad and Behre (1999)](https://ieeexplore.ieee.org/abstract/document/736346). As part of my undergraduate capstone design project, gnss-attitude is being developed for future use in VIOLET, a nanosatellite by [CubeSat NB](https://www.unb.ca/initiatives/cubesat/), to support its camera function. Acting as an extension for [RTKLIB](http://www.rtklib.com/), this software accepts processing results from RTKLIB and outputs the determined antenna boresight vectors in local coordinates (ENU). It can be used for both spacecrafts within GNSS space service range and ground vehicles. **This software is under initial development. Anything is subject to substantial change and considered unstable.**

The missing parts are the calibration for the specific antenna-receiver pair you're using and some treatments according to the orbital altitude; for these, you have to do them yourself. The complete methodologies will be made publicly available in a technical report, expected April 2022.

The author introduced a few innovations:
1. Statistical model for one-step determination of the adjustment terms and the SNR mapping function
2. Use of all 4 globally operating navigation satellite systems: GPS, BeiDou, GLONASS and Galileo
3. Support for multiple frequencies

## Algorithm
1. Before determining the attitude using observation data, a calibration data set is collected. A multiparameter nonlinear regression obtains SNR adjustment parameters for each satellite group and the SNR mapping function.
2. For a given epoch, the line-of-sight (LOS) vectors from satellite to receiver are calculated from satellite and receiver locations. The satellite and receiver locations are derived from the GNSS navigation and observation files.
3. SNR values from the observation file are adjusted according to the adjustment terms developed in Step 1.
4. The off-boresight angles can be found by the SNR mapping function developed in Step 1 and the adjusted SNR values from Step 3.
5. The LOS vectors from Step 2 and off-boresight angles from Step 4 are put into a multiple linear regression to determine the antenna boresight vector.

## Performance
Using a general SNR mapping function, the system delivers an accuracy of about 5° - 20° (RMS). Better performance can be achieved if calibration is performed (under evaluation).

## Limitations
1. Accuracy is subject to the number of signals received and satellite geometry, particularly when the antenna points down, in the woods, outside of GNSS service volume, etc.

2. The single antenna algorithm as implemented here can only determine the boresight vector. Rotation around the boresight axis is undetectable.

3. Would not work with antennas not following the typical gain pattern of a GNSS antenna, such as chip antenna.

## Prerequisites
1. Linux machine (update package list)

        sudo apt update

2. [RTKLIB](http://www.rtklib.com/)
3. GNU Scientific Library

        sudo apt-get install libgsl-dev

4. Build-essential (gcc & make)

        sudo apt install build-essential

## Usage
1. Save AZ/EL/SNR/MP from RTKLIB as `input.txt`

2. Edit `config.c`

3. Compile and run

   Option A: custom input file & multiple frequencies

        make antenna
        ./antenna [input_1] [input_2] [input_3]

   Option B: default input file path & single frequency
 
        make run

## License
gnss-attitude is licensed under the GNU General Public License v3.0. See `LICENSE.txt` for more information.
