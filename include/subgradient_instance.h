#ifndef INCLUDE_SUBGRADIENT_INSTANCE_H_
#define INCLUDE_SUBGRADIENT_INSTANCE_H_

#include <Eigen/Sparse>

#include "../include/cplex_instance.h"

class SubgradientInstance : public Instance {
   public:
    explicit SubgradientInstance(CplexInstance& cinst);
    void solve();

    enum class IterationStatus {
        Unknown,
        ConstraintViolated,
        LowerBoundViolated,
        UpperBoundViolated,
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
    int violated_idx;

    Eigen::SparseVector<double> xstar;
    bool isFeasible(Eigen::SparseVector<double>& x);
};

std::ostream& operator<<(std::ostream& os, const SubgradientInstance::IterationStatus & s);

#endif  // INCLUDE_SUBGRADIENT_INSTANCE_H_
