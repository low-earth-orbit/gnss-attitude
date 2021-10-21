#!/bin/bash
echo "GNSS Single Antenna Attitude Determination"
echo "============== compile simulation =============="
gcc simulation.c convert.c -lm -o simulation
sleep 1s
echo "============== DONE! =============="
echo "============== run simulation =============="
./simulation > input.txt
sleep 2s
echo "============== DONE! =============="
echo "============== compile antenna =============="
gcc antenna.c convert.c -lm -o antenna
sleep 1s
echo "============== DONE! =============="
echo "============== run antenna =============="
./antenna > output.txt
sleep 1s
echo "============== DONE! =============="
echo "============== EXITING =============="