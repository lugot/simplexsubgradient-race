#ifndef INCLUDE_CPLEX_INSTANCE_H_
#define INCLUDE_CPLEX_INSTANCE_H_

#include <string>

#include "../include/instance.h"
#include "ilconcert/iloenv.h"
#include "ilcplex/ilocplexi.h"

class CplexInstance : public Instance {
   public:
    CplexInstance();
    explicit CplexInstance(const std::string& model_name);
    CplexInstance(const CplexInstance& other);
    ~CplexInstance();

    void importModel(const std::string& filename);
    bool solve();

    IloEnv env;
    IloModel model;
    IloObjective obj;
    IloNumVarArray var;
    IloRangeArray rng;
    IloCplex cplex;

    std::string model_name;
};

#endif  // INCLUDE_CPLEX_INSTANCE_H_
