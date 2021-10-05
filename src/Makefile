CC = g++
CFLAGS = -Wall -g -O3 -std=c++17 -pg
RM = rm -rf
AR = ar rc
TARGET = thesis
all: $(TARGET)

OBJS = globals.o main.o cplex.o sparse_vector.o subgradient.o sparse_matrix.o  problem.o step.o
HEADERS =


CPLEX_HOME = /opt/ibm/ILOG/CPLEX_Studio201
LIBS = -L${CPLEX_HOME}/cplex/lib/x86-64_linux/static_pic -L${CPLEX_HOME}/concert/lib/x86-64_linux/static_pic -L. -lilocplex -lconcert -lcplex -lm -lpthread -ldl

EIGEN_HOME = /usr/include/eigen3
INC = -I../include -I${CPLEX_HOME}/cplex/include -I${CPLEX_HOME}/concert/include -I${EIGEN_HOME}


.SUFFIXES:
.SUFFIXES: .o .c .cpp .cc
.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
.cpp.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
.cc.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(TARGET): $(OBJS) $(LIBUTILS)
	$(CC) -o $(TARGET) $(CFLAGS) $(OBJS) $(LIBS)

$(OBJS) : $(HEADERS)

$(LIBUTILS): $(OBJS_LIBUTILS)
	$(AR) $(LIBUTILS) $(OBJS_LIBUTILS)

$(LIBUTILS) : $(HEADERS_LIBUTILS)

clean:
	rm -rf *.o
	rm thesis

again:
	make clean
	make

wow:
	@echo "                                      W O W W W W WWWWWWWWWWWWWWWWWWW"

who:
	@echo "you are user $(USER) with uname `uname` (OS = $(OS)) and you working with compiler setting $(setting)"


#CC = g++
#IDIR = ../include
#CFLAGS = -std=c++0x -Wall -O0 -g -I../include -I$(DIR)

#ODIR = obj

#LIBS = -lm
#_DEPS = parsers.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

#_OBJ = main.o parsers.o
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

#$(ODIR)/%.o: %.c $(DEPS)
	#$(CC) -c -o $@ $< $(CFLAGS)

#main: $(OBJ)
	#$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

#.PHONY: clean

#clean:
	#rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~


#CPLEX_HOME = /opt/ibm/ILOG/CPLEX_Studio201/cplex
#AR = ar rc
#LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl
#INCS = -I../include -I${CPLEX_HOME}/include/ilcplex


#RM = rm -rf

#$(TARGET): $(OBJS) $(LIBUTILS)
	#$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

#main.o: ../src/main.cc $(INCLUDE)
	#$(CXX) $(CXXFLAGS) -c ../src/main.cc

#$(OBJS) : $(HEADERS)

#$(LIBUTILS): $(OBJS_LIBUTILS)
	#$(AR) $(LIBUTILS) $(OBJS_LIBUTILS)

#$(LIBUTILS) : $(HEADERS_LIBUTILS)
