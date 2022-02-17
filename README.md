# gnss-attitude
## Introduction
gnss-attitude is a GNSS-based single antenna attitude determination processor implementing the SNR-based algorithm first proposed by [Axelrad and Behre (1999)](https://doi.org/10.1109/5.736346). As part of my undergraduate capstone design project, gnss-attitude is being developed for future use in VIOLET, a nanosatellite by [CubeSat NB](https://www.unb.ca/initiatives/cubesat/), to support its camera function. Acting as an extension of [RTKLIB](http://www.rtklib.com/), this application accepts processing results from RTKLIB and outputs the determined antenna boresight vectors in local coordinates (ENU).

This program enables absolute boresight determination by a single GNSS antenna directly, without the use of an IMU. With several new techniques employed, preliminary results show that the accuracy is a few degrees RMS, in contrast to 15° RMS reported by [Wang et al. (2005)](https://doi.org/10.2514/6.2005-5993) and within 10° degrees using a combination of two opposite-pointing antennas reported by [Eagleson et al. (2018)](https://digitalcommons.usu.edu/smallsat/2018/all2018/424/).

For potential testers: 
1. This software is under initial development. Anything is subject to substantial change and considered unstable. Initial development release is expected Feburary 2022.
2. The missing part is the calibration for the specific antenna-receiver pair you're using. For this, you have to do it yourself. 
3. The complete methodologies will be made publicly available in a technical report, expected April 2022.

## Algorithm
1. Before determining the attitude using observation data, a calibration data set is collected. A multiparameter nonlinear regression obtains SNR adjustment parameters for each satellite group and the SNR mapping function.
2. For a given epoch, the line-of-sight (LOS) vectors from satellite to receiver are calculated from satellite and receiver locations. The satellite and receiver locations are derived from the GNSS navigation and observation files.
3. SNR values are adjusted according to the adjustment terms developed in Step 1.
4. The off-boresight angles can be found using the SNR mapping function developed in Step 1 and the adjusted SNR values from Step 3.
5. The LOS vectors from Step 2 and off-boresight angles from Step 4 are put into an ordinary least squares for the antenna boresight vector.

## Performance
Using a general SNR mapping function, the system delivers an accuracy of about 5° - 20° (RMS). Better performance can be achieved if calibration is performed (under evaluation).

## Limitations
1. Accuracy depends on the noise in SNR and a number of other factors.

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
Under development

## Planned
2022/02/  - initial development release

## License
gnss-attitude is licensed under the GNU General Public License v3.0. See `LICENSE.txt` for more information.

