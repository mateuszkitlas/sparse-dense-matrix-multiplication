CXX     = mpic++ --std=c++11
CC     = mpicc
CXXFLAGS= -Wall -c -O3 ${ARGS}
LDFLAGS = -Wall -O3 ${ARGS}
ALL     = matrixmul

HEADERS = metgen.h sparse.h fullsparse.h common.h dense.h compute.h densematgen.h

all: $(ALL)

$(ALL): %: %.o densematgen.o sparse.o sparse_mpi.o fullsparse.o common.o dense.o compute.o
	$(CXX) $(LDFLAGS) $^ -o $@

%.o: %.cpp %.c $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $@ $<

clean:
	rm -f *.o *core *~ *.out *.err $(ALL)
