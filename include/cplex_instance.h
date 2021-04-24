#ifndef INCLUDE_CPLEX_INSTANCE_H_
#define INCLUDE_CPLEX_INSTANCE_H_

#include <string>

#include "Eigen/Sparse"
#include "ilconcert/iloenv.h"
#include "ilcplex/ilocplexi.h"

typedef std::tuple<Eigen::SparseMatrix<double>*, Eigen::SparseVector<double>*,
                   Eigen::SparseVector<double>*, Eigen::SparseVector<double>*,
                   Eigen::SparseVector<double>*>
    CanonicalForm;

/*
 * Just a wrapper.
 */
class CplexInstance {
   public:
    CplexInstance();
    // TODO(lugot): LEARN explicit constructors
    explicit CplexInstance(const std::string& model_name);
    explicit CplexInstance(const CplexInstance& other);

    // wrapper of IloAlgorithm::Status
    enum class Status {
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

    // TODO(lugot): LEARN const
    void importModel(const std::string& model_name);
    bool solve();
    bool isSolved();
    Status getStatus();
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

    Status status;
};

// TODO(lugot): LEARN friend
std::ostream& operator<<(std::ostream& os, const CplexInstance::Status& s);

#endif  // INCLUDE_CPLEX_INSTANCE_H_
