#include "../include/subgradient_instance.h"

#include "../include/globals.h"
#include "../include/utils.h"

using std::endl;

SubgradientInstance::SubgradientInstance(CplexInstance& cinst) : Instance() {
    std::tie(A, b, c, lb, ub) = cinst.getCanonicalForm();
    m = A.rows();
    n = A.cols();
    assert(b.rows() == m);
    assert(c.rows() == n);
    assert(lb.rows() == n);
    assert(ub.rows() == n);

    status = SolutionStatus::NeverSolved;
    iter_status = IterationStatus::Unknown;
}

void SubgradientInstance::solve() {
    Eigen::SparseVector<double> x(n);

    std::cout << "A:" << A;
    std::cout << "b: " << b;
    std::cout << lb << std::endl << ub << std::endl;

    Eigen::SparseVector<double> u(m);

    std::vector<double> phis(P);
    double mu;
    mu = 0.1;
    int lowlambda_lifespan = 0;

    for (int k = 0; k < 100; k++) {
        Eigen::SparseVector<double> temp(n);
        // Eigen::SparseVector<double> temp = c;

        for (int i = 0; i < n; i++) {
            for (Eigen::SparseVector<double>::InnerIterator it(u); it; ++it) {
                temp.coeffRef(i) -= it.value() * A.coeff(it.row(), i);
            }
        }
        temp = temp + c;
        temp.prune(XSMALL);
        std::cout << "temp: " << temp;

        for (int i = 0; i < n; i++) {
            if (temp.coeff(i) < 0) {
                x.coeffRef(i) = ub.coeff(i);
            } else if (temp.coeff(i) > 0) {
                x.coeffRef(i) = lb.coeff(i);
            }
        }
        x.prune(XSMALL);

        Eigen::SparseVector<double> s = b - A * x;
        s.prune(XSMALL);

        double phi_k = phi(u, x);
        phis[k % P] = phi_k;

        double phimax, phimin;  // TODO(lugot): change
        if (k > P) {
            phimin = *std::min(phis.begin(), phis.end());
            phimax = *std::max(phis.begin(), phis.end());

            if (phimax - phimin > 0.01 * phimin) {
                mu = 0.5 * mu;
            } else if (phimax - phimin < 0.001 * phimin) {
                mu = std::min(2.0, 1.5 * mu);
            }
            /* else mu remains mu */
        }
        /* else mu remains mu */

        double lambda = mu * (449.606 - phi_k) / (s.norm() * s.norm());
        // double lambda = 1.0 / (k+1);
        //
        if (lambda < LOWLAMBDA_TRESHOLD) {
            lowlambda_lifespan++;
            if (lowlambda_lifespan >= LOWLAMBDA_MAXLIFE) {
                std::cout << "Exit: lambda less than " << LOWLAMBDA_TRESHOLD
                          << "for more than " << LOWLAMBDA_MAXLIFE
                          << "iterations." << endl;

                return;
            }
        } else {
            lowlambda_lifespan = 0;
        }

        u = u + lambda * s;
        for (Eigen::SparseVector<double>::InnerIterator it(u); it; ++it) {
            it.valueRef() = std::max(0.0, it.valueRef());
        }
        u.prune(XSMALL);

        iterations.push_back(new Iteration(k, u, x, s, phi_k, lambda, mu));
        std::cout << *iterations[iterations.size() - 1] << std::endl;
    }
}

double SubgradientInstance::phi(const Eigen::SparseVector<double>& u,
                                const Eigen::SparseVector<double>& x) {
    // TODO(lugot): PERFORMANCE
    return c.dot(x) + u.dot(b - A * x);
}
/* double step_size(const std::vector<SubgradientInstance::Iteration>&
   iterations, int k) {} */

double SubgradientInstance::constraintsUpperBound() {
    Eigen::SparseVector<double> x(n);
    for (int i = 0; i < n; ++i) x.coeffRef(i) = 1.0;

    // TODO(lugot): BISECTION
    // TODO(lugot): IMPLEMENT equalities
    while (!isConstraintFeasible(x)) x *= 2.0;

    return x.coeff(0);
}

bool SubgradientInstance::isConstraintFeasible(
    const Eigen::SparseVector<double>& x) {
    Eigen::SparseVector<double> Ax = A * x;
    Eigen::SparseVector<double>::InnerIterator itAx(Ax), itb(b);

    while (itAx && itAx) {
        if (itAx.col() == itb.col()) {
            // (Ax)i < (b)i -> false
            if (itAx.value() - itb.value() < EPSILON) return false;

            ++itAx;
            ++itb;
        } else if (itAx.col() < itb.col()) {
            if (itAx.value() < EPSILON) return false;

            ++itAx;
        } else {
            if (-itb.value() < EPSILON) return false;

            ++itb;
        }
    }
    while (itAx) {
        if (itAx.value() < EPSILON) return false;

        ++itAx;
    }
    while (itb) {
        if (-itb.value() < EPSILON) return false;

        ++itb;
    }

    return true;
}

SubgradientInstance::Iteration::Iteration(int k, Eigen::SparseVector<double> u,
                                          Eigen::SparseVector<double> x,
                                          Eigen::SparseVector<double> s,
                                          double dual_objective,
                                          double step_size,
                                          double step_size_param) {
    this->k = k;
    this->u = u;
    this->x = x;
    this->s = s;
    this->dual_objective = dual_objective;
    this->step_size = step_size;
    this->step_size_param = step_size_param;
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

std::ostream& operator<<(std::ostream& os,
                         const SubgradientInstance::Iteration& i) {
    os << "Iteration " << i.k << ":\n";
    if (VERBOSE) {
        os << "x: " << i.x << "\n";
    } else {
        os << "x: [...]"
           << "\n";
    }
    if (VERBOSE) {
        os << "u: " << i.u << " " << i.u.norm() << "\n";
    } else {
        os << "u: [...] " << i.u.norm() << "\n";
    }
    if (VERBOSE) {
        os << "s: " << i.s << " " << i.s.norm() << "\n";
    } else {
        os << "s: [...] " << i.s.norm() << "\n";
    }
    os << "phi(u): " << i.dual_objective << "\n";
    os << "step size (lambda): " << i.step_size << "\n";
    os << "step param (mu): " << i.step_size_param << "\n";

    return os;
}
