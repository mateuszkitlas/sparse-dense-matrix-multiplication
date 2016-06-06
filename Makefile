CXX     = mpic++
CXXFLAGS= -Wall -c -O3 ${ARGS}
LDFLAGS = -Wall -O3 ${ARGS}
ALL     = matrixmul

HEADERS = metgen.h sparse.h fullsparse.h

all: $(ALL)

$(ALL): %: %.o densematgen.o sparse.o sparse_mpi.o fullsparse.o
	$(CXX) $(LDFLAGS) $^ -o $@

%.o: %.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $@ $<


clean:
	rm -f *.o *core *~ *.out *.err $(ALL)
