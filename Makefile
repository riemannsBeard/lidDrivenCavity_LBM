CC = g++
CFLAGS = -Wall -g -O1 -std=c++11
LFLAGS = 
LIBS = -larmadillo
SRCS = main.C LBM.C
OBJECTS = $(SRCS:.C=.o)
PROGRAM = lbm

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
		$(CC) $(CFLAGS) $(LIBS) -o $(PROGRAM) $(OBJECTS) $(LIBS)

.C.o:
	$(CC) $(CFLAGS) $(LIBS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf *.exe *.o *.out 
