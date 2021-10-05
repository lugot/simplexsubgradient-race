#ifndef INCLUDE_CPLEX_H_
#define INCLUDE_CPLEX_H_

#include <string>

#include "ilconcert/iloenv.h"
#include "ilcplex/ilocplexi.h"

class Cplex {
   public:
    Cplex(){};
    explicit Cplex(const std::string& model_name) { importModel(model_name); }
    Cplex(const Cplex& other);
    ~Cplex();

    void importModel(const std::string& filename);
    bool solve();

    IloEnv env;
    IloModel model;
    IloObjective obj;
    IloNumVarArray var;
    IloRangeArray rng;
    IloCplex cplex;

    std::string model_name;
    double solvetime;

    // wrapper of IloAlgorithm::Status
    enum class SolutionStatus {
        Unknown,
        Feasible,
        Optimal,
        Infeasible,
        Unbounded,
        InfeasibleOrUnbounded,
        Error,
        EmptyInstance,
        NeverSolved,
        ModelChanged,
    };

    std::ofstream log;

   protected:
    SolutionStatus status;
};

#endif  // INCLUDE_CPLEX_H_
