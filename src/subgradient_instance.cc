#include "../include/subgradient_instance.h"

#include "../include/globals.h"

SubgradientInstance::SubgradientInstance(CplexInstance& cinst) : Instance() {
    std::tie(A, b, c, lb, ub) = cinst.getCanonicalForm();
    m = A.rows();
    n = A.cols();
    assert(b.rows() == m);
    assert(c.rows() == n);
    assert(lb.rows() == n);
    assert(ub.rows() == n);

    qp = CplexInstance(cinst);
    status = SolutionStatus::NeverSolved;
    iter_status = IterationStatus::Unknown;
}

void SubgradientInstance::solve() {
    Eigen::SparseVector<double> x(n);
    bool feasible = isFeasible(x);

    std::cout << feasible << std::endl;
}

bool SubgradientInstance::isFeasible(Eigen::SparseVector<double>& x) {
    assert(x.rows() == n);

    // check for constraints: merge-like code structure for constraint
    // O( nnz(lhs) + nnz(b) ) = O ( 2m )
    // anyway, nnz(b) should ~ m
    Eigen::SparseVector<double> lhs = A * x;
    std::cout << x << " " << lhs << " " << b << std::endl;

    Eigen::SparseVector<double>::InnerIterator itl(lhs);
    Eigen::SparseVector<double>::InnerIterator itb(b);
    while (itl && itb) {
        if (itl.index() < itb.index()) {
            // since itb.index is zero
            // TODO(lugot): CHECK floating point safety
            if (itl.value() < -EPS) return false;
            ++itl;
        } else if (itl.index() > itb.index()) {
            // since itl.index is zero
            if (itb.value() > EPS) return false;
            ++itb;
        } else {  // itl.index() == itb.index()
            // itl.value() > itb.value() -> infeasible
            if (fabs(itb.value() - itl.value()) < EPS) return false;
            ++itl;
            ++itb;
        }
    }
    while (itl) {
        // since itb.index is zero
        // TODO(lugot): CHECK floating point safety
        if (itl.value() < -EPS) return false;
        ++itl;
    }
    while (itb) {
        // since itl.index is zero
        if (itb.value() > EPS) return false;
        ++itb;
    }

    return true;
}
