GCC = gcc
CFLAGS = -g -Wall 
EXTRAFLAGS =  -lgsl -lgslcblas -lm

antenna: clean util.o struct.o snr.o antenna.o
			$(GCC) $(CFLAGS) util.o struct.o snr.o antenna.o -o antenna $(EXTRAFLAGS)

antenna.o: antenna.c config.h
			$(GCC) $(CFLAGS) -c antenna.c $(EXTRAFLAGS)

snr.o: snr.c
			$(GCC) $(CFLAGS) -c snr.c $(EXTRAFLAGS)

util.o: util.c
			$(GCC) $(CFLAGS) -c util.c $(EXTRAFLAGS)

struct.o: struct.c
			$(GCC) $(CFLAGS) -c struct.c $(EXTRAFLAGS)

simulation: clean util.o struct.o simulation.o
			$(GCC) $(CFLAGS) util.o struct.o simulation.o -o simulation $(EXTRAFLAGS)
			
simulation.o: simulation.c config.h
			$(GCC) $(CFLAGS) -c simulation.c $(EXTRAFLAGS)

run: clean antenna
			./antenna

run_simulation: clean antenna simulation
			./simulation
			./antenna

clean:
			rm -f *.o antenna simulation
