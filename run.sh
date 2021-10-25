#!/bin/bash
echo "GNSS Single Antenna Attitude Determination"
echo "Data Generation ..."
SECONDS=0
gcc -Wall simulation.c mathutil.c -lm -o simulation
./simulation > input.txt
echo "Completed in $SECONDS seconds"
echo "Processing ..."
SECONDS=0
gcc -Wall antenna.c mathutil.c -o antenna -lgsl -lgslcblas -lm
./antenna > output.txt
echo "Completed in $SECONDS seconds"
echo "Exiting ..."