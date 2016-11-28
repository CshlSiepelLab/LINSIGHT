# compiler
CPP_C = g++

# the compiler flag
CFLAG = -g -O2 -Wall -fno-math-errno -std=c++11

# the linking flag
# The static link is bad for Valgrind
# LIFLAG = $(CFLAG) -static
LIFLAG = $(CFLAG)

# end flag
# LEFLAG = -lgsl -lgslcblas -lm
LEFLAG =

# final target
all: LINSIGHT-fit  LINSIGHT-score LINSIGHT-prep

LINSIGHT-fit: L-INSIGHT-fit.o
	$(CPP_C) $(LIFLAG) -o $@ $^ $(LEFLAG)

LINSIGHT-score: L-INSIGHT-pred.o
	$(CPP_C) $(LIFLAG) -o $@ $^ $(LEFLAG)

LINSIGHT-prep: L-INSIGHT-prep.o
	$(CPP_C) $(LIFLAG) -o $@ $^ $(LEFLAG)

L-INSIGHT-fit.o L-INSIGHT-pred.o L-INSIGHT-prep.o: %.o: %.cpp
	$(CPP_C) $(CFLAG) -o $@ -c $<
