# gnss-attitude
## Introduction
gnss-attitude is a GNSS-based single antenna attitude determination processor implementing the SNR-based algorithm first proposed by [Axelrad and Behre (1999)](https://doi.org/10.1109/5.736346). As part of my undergraduate capstone design project, gnss-attitude is being developed for future use in VIOLET, a nanosatellite by [CubeSat NB](https://www.unb.ca/initiatives/cubesat/), to support its camera function. Acting as an extension of [RTKLIB](http://www.rtklib.com/), this application accepts processing results from RTKLIB and outputs the determined antenna boresight vectors in local coordinates (ENU).

This program enables absolute boresight determination by a single GNSS antenna directly, without the use of an IMU. With several new techniques employed, preliminary results show that the accuracy is a few degrees RMS, in contrast to 15° RMS reported by [Wang et al. (2005)](https://doi.org/10.2514/6.2005-5993) and within 10° degrees using a combination of two opposite-pointing antennas reported by [Eagleson et al. (2018)](https://digitalcommons.usu.edu/smallsat/2018/all2018/424/).

For potential testers: 
1. This software is under development for evaluating the performance of such a system. Code is subject to change and considered unstable.
2. The missing part is the calibration for the specific antenna-receiver pair you're using. For this, you have to do it yourself. 
3. The complete methodologies will be made publicly available in a technical report, expected April 2022.

## Algorithm
1. A calibration data set is collected. For each constellation at each signal channel, a multiparameter nonlinear regression obtains SNR adjustment parameters for each satellite and the SNR mapping function.
2. The line-of-sight (LOS) vectors from satellite to receiver are calculated based on satellite and receiver locations for a given epoch.
3. SNR values from the observation file are adjusted according to the adjustment terms developed in [1].
4. The off-boresight angles can be found by the SNR mapping function developed in [1] and the adjusted SNR values from [3].
5. Boresight can be determined by least squares using LOS vectors [2] and off-boresight angles [4].

## Performance
If multipath is not a concern, the algorithm can achieve degree-level RMS accuracy. Bias exists for most ground environments.

## Limitations
1. Highly susceptible to multipath effects.

2. Rotation around the boresight axis is undetectable.

3. Does not work with GNSS antennas not following the ideal gain pattern.

4. Requires one-time calibration.

## Prerequisites
1. Linux machine (update package list)

        sudo apt update

2. [RTKLIB](http://www.rtklib.com/)
3. GNU Scientific Library

        sudo apt-get install libgsl-dev

4. Build-essential (gcc & make)

        sudo apt install build-essential

## Usage
1. In RTKLIB, save AZ/EL/SNR/MP of the selected frequency as a text file such as `input.txt`

2. Update `snr.c` with parameters obtained from calibration

3. Edit `config.c` if necessary

4. Compile and run

   Option A: custom input file & multiple frequencies

        make antenna
        ./antenna [input_1] [input_2] [input_3]

   Option B: default input file path `input.txt` & single frequency
 
        make run

## Road map
2022/03/01 - initial development release

## Planned


## License
gnss-attitude is licensed under the GNU General Public License v3.0. See `LICENSE.txt` for more information.

