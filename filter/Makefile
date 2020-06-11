CC=g++
CFLAGS=-g -Wall
INCLUDE=-I./bamtools_CPP/include
### Lab
LIBPATH=-L./bamtools_CPP/lib
lib=-lbamtools

VAR=ONLY_MUTATION_DIR_READ

filter:filter.cpp filter.h
	#$(CC) -D$(VAR) $(CFLAGS) $(INCLUDE) $(LIBPATH) $(lib) -o filter filter.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBPATH) $(lib) -o filter filter.cpp
