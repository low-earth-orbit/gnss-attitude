#!/bin/bash
echo "GNSS Single Antenna Attitude Determination"
echo "============== compile simulation =============="
gcc -Wall simulation.c mathutil.c -lm -o simulation
sleep 1s
echo "============== DONE! =============="
echo "============== run simulation =============="
./simulation > input.txt
sleep 2s
echo "============== DONE! =============="
echo "============== compile antenna =============="
gcc -Wall antenna.c mathutil.c -o antenna -lgsl -lgslcblas -lm
sleep 1s
echo "============== DONE! =============="
echo "============== run antenna =============="
./antenna > output.txt
sleep 1s
echo "============== DONE! =============="
echo "============== EXITING =============="
sleep 1s