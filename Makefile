GCC = gcc
CFLAGS = -g -Wall 
EXTRAFLAGS =  -lgsl -lgslcblas -lm

antenna: clean util.o struct.o snr.o antenna.o config.h
			$(GCC) $(CFLAGS) util.o struct.o snr.o antenna.o -o antenna $(EXTRAFLAGS)

simulation: clean util.o struct.o simulation.o
			$(GCC) $(CFLAGS) util.o struct.o simulation.o -o simulation $(EXTRAFLAGS)

simulation.o: simulation.c
			$(GCC) $(CFLAGS) -c simulation.c $(EXTRAFLAGS)

antenna.o: antenna.c
			$(GCC) $(CFLAGS) -c antenna.c $(EXTRAFLAGS)

snr.o: snr.c
			$(GCC) $(CFLAGS) -c snr.c $(EXTRAFLAGS)

util.o: util.c
			$(GCC) $(CFLAGS) -c util.c $(EXTRAFLAGS)

struct.o: struct.c
			$(GCC) $(CFLAGS) -c struct.c $(EXTRAFLAGS)

run_antenna: clean antenna
			./antenna

run_simulation: clean simulation antenna
			./simulation
			./antenna

clean:
			/bin/rm -f *.o antenna simulation
