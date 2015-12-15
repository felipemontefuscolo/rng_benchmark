CXXFLAGS=-g -std=c++11 -fopenmp -O3 -march=native -mtune=native
#LDFLAGS= -static -ltrng4 -dynamic
LDFLAGS=/usr/local/lib/libtrng4.a -lboost_program_options
EFLAGS=-I/home/fuscolo/src/Random123-Boost
do:
	g++ ${CXXFLAGS} ${EFLAGS} bench.cpp -o bench ${LDFLAGS}
