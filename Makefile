include /usr/local/include/cantera/Cantera.mak

CC=gcc
CXX=g++
RM=rm -f
CCFLAGS=-g -std=c++11
CFLAGS=-O3 -c
CPPFLAGS=$(CANTERA_INCLUDES) -g -std=c++0x
#CPPFLAGS=-I/usr/local/include -g
LDFLAGS=-g 
LDLIBS=$(CANTERA_LIBS)
#LDLIBS=-L/usr/local/lib -L/home/wwang/Soft/sundials/lib

LIBS=$(wildcard ./src/*.h)
SRCS=$(wildcard ./src/*.cpp)
DIR=$(notdir $(SRCS))
OBJS=$(addprefix ./obj/,$(patsubst %.cpp,%.o,$(DIR)))
#SRCS=Reaction.cpp ReimannSolver.cpp exchanging.cpp
#OBJS=$(subst .cpp,.o,$(SRCS))

all: Reaction

Reaction: $(OBJS) 
	$(CXX) $(LDFLAGS) $(CCFLAGS) -Xpreprocessor -fopenmp -o Reaction $(OBJS) $(LDLIBS) 
#	$(CXX) $(LDFLAGS) $(CCFLAGS) -o Reaction $(OBJS) $(LDLIBS) -lsundials_cvodes -lsundials_nvecserial -lpthread -lcantera -lcantera_fortran

$(OBJS): ./obj/%.o: ./src/%.cpp $(LIBS)
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) *~ 
