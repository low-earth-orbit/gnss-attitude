#!/bin/bash
echo "GNSS Single Antenna Attitude Determination"
echo "Simulation Running"
echo "Data Generation ..."
SECONDS=0
gcc -Wall simulation.c util.c struct.c -o simulation  -lgsl -lgslcblas -lm
./simulation
echo "Completed in $SECONDS seconds"
echo "Processing ..."
SECONDS=0
gcc -Wall antenna.c util.c struct.c snr.c -o antenna -lgsl -lgslcblas -lm
./antenna
echo "Completed in $SECONDS seconds"
echo "Exiting ..."