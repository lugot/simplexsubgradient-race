#ifndef INCLUDE_SUBGRADIENT_INSTANCE_H_
#define INCLUDE_SUBGRADIENT_INSTANCE_H_

#include <Eigen/Sparse>

#include "../include/cplex_instance.h"

class SubgradientInstance : public Instance {
   public:
    SubgradientInstance(CplexInstance& cinst);
    void solve();

    enum class IterationStatus {
        Unknown,
        ConstraintsInfeasible,
        VariableBOundsInfeasible,
        Feasible
    };

   private:
    int m, n;
    Eigen::SparseMatrix<double> A;
    Eigen::SparseVector<double> b;
    Eigen::SparseVector<double> c;
    Eigen::SparseVector<double> lb;
    Eigen::SparseVector<double> ub;

    CplexInstance qp;

    IterationStatus iter_status;

    Eigen::SparseVector<double> xstar;
    bool isFeasible(Eigen::SparseVector<double>& x);
};

#endif  // INCLUDE_SUBGRADIENT_INSTANCE_H_
