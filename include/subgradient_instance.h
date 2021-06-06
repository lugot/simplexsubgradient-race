#ifndef INCLUDE_SUBGRADIENT_INSTANCE_H_
#define INCLUDE_SUBGRADIENT_INSTANCE_H_

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "../include/cplex_instance.h"
#include "../include/sparse_matrix.h"
#include "../include/sparse_vector.h"

class SubgradientInstance {
   public:
    enum class Methods {
        Pure,
        Deflected,
        Conditional,
        Hybrid,
    };

    enum class Status {
        Error,
        NeverSolved,
        LowLambdaSoft,
        LowLambdaHard,
        ReachedMaxIterations,
        ReachedTimelimit,
        Optimal,
        EpsilonOptimal
    };

    explicit SubgradientInstance(const CplexInstance& cinst);
    Status solve(Methods method);

   private:
    // model
    int m, n;
    SparseVector c, b;
    SparseMatrix A;
    std::vector<double> lb, ub;
    std::vector<bool> equal;  // TODO(lugot): CHANGE in classic sparse vector

    std::string model_name;

    // IterationStatus iter_status;
    Status status;
    SparseVector xstar, ustar;
    double fstar, phistar, gap;

    // correctness
    bool bestsol_avaiable;
    SparseVector best_ustar, best_xstar;

    // method variances
    Status solvePure();
    Status solveDeflected();
    Status solveConditional();
    Status solveHybrid();

    void solveSP(SparseVector* x, const SparseVector& u);
    bool directionInfeasible(const SparseVector& u, const SparseVector& s);

    // save solution
    void saveSolutions();

    // instance converter methods
    enum Sense { GREATER, EQUAL, LESS };
    std::map<std::string, int> varindexer;
    void extractConstraint(const IloRange& r, int row,
                           std::vector<std::tuple<int, int, double>>* a,
                           double* b);
    Sense constraintSense(const IloRange& r);
};

/* std::ostream& operator<<(std::ostream& os,
                         const SubgradientInstance::IterationStatus& s); */

/* std::ostream& operator<<(std::ostream& os,
                         const SubgradientInstance::Iteration& i); */

#endif  // INCLUDE_SUBGRADIENT_INSTANCE_H_
