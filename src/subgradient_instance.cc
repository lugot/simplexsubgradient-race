#include "../include/subgradient_instance.h"

#include "../include/globals.h"
#include "../include/utils.h"

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
    for (int i = 0; i < n; i++) x.coeffRef(i) = 1.0;
    bool feasible = isFeasible(x);

    std::cout << feasible << " " << iter_status << " " << violated_idx
              << std::endl;

    std::cout << "HERE" << std::endl;
    qp.updateObjective(x);
    qp.solve();
}

bool SubgradientInstance::isFeasible(Eigen::SparseVector<double>& x) {
    assert(x.rows() == n);
    // TODO(lugot): TEST
    // TODO(lugot): ORDER of checks for performance improvement
    // TODO(lugot): ADD subfunctions to avoid projections
    // .. maybe do the check is slower than perform projection on variables?

    // check for constraints
    Eigen::SparseVector<double> lhs = A * x;
    lhs.prune(EPS);  // TODO(lugot): CHECK if is worth it
    
    std::cout << lhs << "\n" << b << std::endl;

    if ((violated_idx = sparseVectorCompareZero(lhs, b)) != m) {
        // TODO(lugot): UNK track index?
        iter_status = IterationStatus::ConstraintViolated;
        return false;
    }
    lhs.finalize();

    // check for variables lower bounds
    if ((violated_idx = sparseVectorCompareZero(lb, x)) != m) {
        // TODO(lugot): UNK track index?
        iter_status = IterationStatus::LowerBoundViolated;
        return false;
    }

    // check for variables upper ounds
    if ((violated_idx = sparseVectorCompareZero(lb, x)) != m) {
        // TODO(lugot): UNK track index?
        iter_status = IterationStatus::LowerBoundViolated;
        return false;
    }

    iter_status = IterationStatus::Feasible;
    return true;
}

std::ostream& operator<<(std::ostream& os,
                         const SubgradientInstance::IterationStatus& s) {
    switch (s) {
        case SubgradientInstance::IterationStatus::Unknown:
            os << "Unknownk";
            break;
        case SubgradientInstance::IterationStatus::ConstraintViolated:
            os << "ConstraintViolated";
            break;
        case SubgradientInstance::IterationStatus::LowerBoundViolated:
            os << "LowerBoundViolated";
            break;
        case SubgradientInstance::IterationStatus::UpperBoundViolated:
            os << "UpperBoundViolated";
            break;
        case SubgradientInstance::IterationStatus::Feasible:
            os << "Feasible";
            break;
    }

    return os;
}
