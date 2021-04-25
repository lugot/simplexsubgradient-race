#ifndef INCLUDE_CPLEX_INSTANCE_H_
#define INCLUDE_CPLEX_INSTANCE_H_

#include <string>

#include "../include/instance.h"
#include "Eigen/Sparse"
#include "ilconcert/iloenv.h"
#include "ilcplex/ilocplexi.h"
//#include "ilconcert/iloenv.h"
//#include "ilcplex/ilocplexi.h"

typedef std::tuple<Eigen::SparseMatrix<double>, Eigen::SparseVector<double>,
                   Eigen::SparseVector<double>, Eigen::SparseVector<double>,
                   Eigen::SparseVector<double>>
    CanonicalForm;

/*
 * Just a wrapper.
 */
class CplexInstance : public Instance {
   public:
    CplexInstance();
    explicit CplexInstance(const std::string& model_name);
    CplexInstance(const CplexInstance& other);

    // TODO(lugot): LEARN const
    void importModel(const std::string& model_name);
    bool solve();
    void printStatus();

    CanonicalForm getCanonicalForm();

   private:
    // TODO(lugot): LEARN better understand handler thinghy
    IloEnv env;
    IloModel model;
    IloObjective obj;
    IloNumVarArray var;
    IloRangeArray rng;
    IloCplex cplex;
};

#endif  // INCLUDE_CPLEX_INSTANCE_H_
