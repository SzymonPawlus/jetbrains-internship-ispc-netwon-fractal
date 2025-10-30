all: main

newton.o: newton.ispc
	ispc newton.ispc -o newton.o -h newton.h

main: main.cpp newton.o
	g++ -o main newton_serial.cpp main.cpp newton.o -lm
