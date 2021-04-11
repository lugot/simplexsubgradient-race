#CC = g++
#CFLAGS = -Wall -O0 -g
RM = rm -rf
AR = ar rc

SRC_DIR := ../src
OBJ_DIR := obj
BIN_DIR := .

EXE := $(BIN_DIR)/main
SRC := $(wildcard $(SRC_DIR)/*.cc)
OBJ := $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

CPLEX_HOME := /opt/ibm/ILOG/CPLEX_Studio201
EIGEN_HOME := /usr/include/eigen3

CPPFLAGS := -I../include -I${CPLEX_HOME}/cplex/include -I${CPLEX_HOME}/concert/include -I${EIGEN_HOME} -MMD -MP
CXXFLAGS := -Wall -g
LDFLAGS  := -L${CPLEX_HOME}/cplex/lib/x86-64_linux/static_pic -L${CPLEX_HOME}/concert/lib/x86-64_linux/static_pic -L.
LDLIBS   := -lm -lilocplex -lconcert -lcplex -lm -lpthread -ldl
CXX = g++

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJS) $(LIBUTILS) | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@


clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)

#-include $(OBJ:.o=.d)
