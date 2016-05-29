CXX     = mpic++
CC      = mpic++
CXXFLAGS= -Wall -c -O3
CFLAGS  = -Wall -c --std=c99 -O3
LDFLAGS = -Wall -O3 --std=c99
ALL     = matrixmul
MATGENFILE = densematgen.o

HEADERS = metgen.h

all: $(ALL)

$(ALL): %: %.o $(MATGENFILE)
	$(CXX) $(LDFLAGS) $^ -o $@

%.o: %.c $(HEADERS) Makefile
	$(CC) $(CFLAGS) $@ $<

%.o: %.cpp $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $@ $<


clean:
	rm -f *.o *core *~ *.out *.err $(ALL)
