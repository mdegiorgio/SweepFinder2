#COPTS=-g -pg
COPTS = -O3 -Wall

SweepFinder2: SweepFinder2.o freq.o factorials.o bfgs.o sort.o my_rand.o
	gcc -o SweepFinder2 SweepFinder2.o $(COPTS) freq.o factorials.o bfgs.o sort.o my_rand.o -lm

SweepFinder2.o: SweepFinder2.c SweepFinder2.h
	gcc -c SweepFinder2.c $(COPTS)

freq.o: freq.c freq.h
	gcc -c freq.c $(COPTS)

factorials.o: factorials.c factorials.h
	gcc -c factorials.c $(COPTS)

bfgs.o: bfgs.c bfgs.h
	gcc -c bfgs.c $(COPTS)

sort.o: sort.c sort.h
	gcc -c sort.c $(COPTS)

my_rand.o: my_rand.c my_rand.h
	gcc -c my_rand.c $(COPTS)
