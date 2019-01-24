# This is a makefile.
# This makes a parallel simulation for same dielectric problem
# Use option -p in CC for profiling with gprof

PROG = confined_ions

SRC = main.cpp functions.cpp interface.cpp  parallel_precal.cpp cpmd.cpp realforces.cpp fakeforces.cpp fmd.cpp energies.cpp

OBJ = ${SRC:.cpp=.o}

.SUFFIXES: .cpp .o

#CC = g++ -O3 -g -Wall
# do not use -fopenmp unless parallelizing properly. 
CC = g++ -O3 -g -Wall -fopenmp
# profiling. testing speed
#CC = g++ -static -O3 -g -Wall -fopenmp -p

LFLAG = -lgsl -lgslcblas -lboost_program_options -lm

#LIBS = -I/shared/gsl-1.16/include/ -L/shared/gsl-1.16/lib/ -I/usr/include/ -L/usr/lib64/
#LIBS = -I/shared/gsl-1.16/include/ -L/shared/gsl-1.16/lib/ 
#LIBS = -I/software/gsl/1.14/include/ -L/software/gsl/1.14/lib -I/software/boost/1.41.0/include/ -L/software/boost/1.41.0/lib/ 
#LIBS = -I/usr/include/ -L/usr/lib
LIBS = -I/shared/gsl-1.16/include/ -L/shared/gsl-1.16/lib/ -I/shared/boost/include/ -L/shared/boost/lib/

CFLAG = -c

OFLAG = -o

$(PROG) : $(OBJ)
#	$(CC) $(OFLAG) $@ $^ $(LFLAG)
	$(CC) $(OFLAG) $@ $^ $(LIBS) $(LFLAG)
#	$(CC) $(LIBS) $(LFLAG) $(OFLAG) $@ $^

.cpp.o :
#	$(CC) $(CFLAG) $(OFLAG) $@ $< 
	$(CC) $(CFLAG) $(OFLAG) $@ $< $(LIBS) $(LFLAG)
#	$(CC) $(LIBS) $(CFLAG) $(OFLAG) $@ $<

clean:
	rm -f *.o

dataclean:
	rm -f outfiles/*.dat outfiles/*.xyz  outfiles/*.lammpstrj  datafiles/*.dat outfiles/verifiles/*.dat
