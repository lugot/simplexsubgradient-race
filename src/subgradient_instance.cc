#include "../include/subgradient_instance.h"

#include <iomanip>

#include "../include/globals.h"
#include "../include/utils.h"

SubgradientInstance::SubgradientInstance(const CplexInstance& cinst) {
    // model sizes
    n = cinst.var.getSize();
    m = cinst.rng.getSize();

    // check if some constraits are equalities, add rows in case
    for (int i = 0; i < cinst.rng.getSize(); ++i) {
        if (constraintSense(cinst.rng[i]) == EQUAL) m++;
    }
    ustar = SparseVector(m);
    xstar = SparseVector(n);

    // vectors & friends
    c = SparseVector(n);
    b = SparseVector(m);
    A = SparseMatrix(m, n);
    lb = std::vector<double>(n, LOWERBOUND);
    ub = std::vector<double>(n, UPPERBOUND);

    // model name for file retreiving
    model_name = cinst.model_name;

    // track the names of the variables for consistency
    varindexer = std::map<std::string, int>();

    // fill the objective sparse vector
    IloExpr::LinearIterator ito;
    for (ito = cinst.obj.getLinearIterator(); ito.ok(); ++ito) {
        std::string varname = ito.getVar().getName();
        if (varindexer.find(varname) == varindexer.end()) {
            varindexer[varname] = varindexer.size();
        }

        c.push_back({varindexer[varname], ito.getCoef()});
    }

    // fill A and b (formulation Ax >= b)
    std::vector<std::tuple<int, int, double>> lhs;
    double rhs;
    int row = 0;
    for (int i = 0; i < cinst.rng.getSize(); ++i) {
        IloRange r = cinst.rng[i];

        // extract lhs and rhs from IloRange object ..
        extractConstraint(r, row, &lhs, &rhs);
        // .. and add them to the sparse vector and matrix
        A.addRow(lhs);
        b.push_back({row, rhs});
        row++;

        if (constraintSense(r) == EQUAL) {
            // swap the coefficents in sign, both lhs and rhs
            std::vector<std::tuple<int, int, double>>::iterator itl;
            for (itl = lhs.begin(); itl != lhs.end(); ++itl) {
                std::get<0>(*itl) = row;  // next row
                // std::get<1>(*it) untouched
                std::get<2>(*itl) = -std::get<2>(*itl);  // swap sign
            }
            rhs = -rhs;

            A.addRow(lhs);
            b.push_back({row, rhs});
            row++;
        }
    }

    // change the variable bounds, if present
    for (int i = 0; i < cinst.var.getSize(); ++i) {
        lb[i] = std::max(lb[i], cinst.var[i].getLB());
        ub[i] = std::min(ub[i], cinst.var[i].getUB());
    }

    // iniitialize the other members
    status = Status::NeverSolved;

    // initialize to suppress warning
    fstar = phistar = gap = 0.0;

    // retreive best optimal solution if avaiable
    std::ifstream optvectors("../subgradient_solutions/" + model_name);
    if (optvectors.is_open()) {
        double nnz;

        optvectors >> nnz;
        best_ustar = SparseVector(m);
        while (nnz--) {
            int pos;
            double val;
            optvectors >> pos >> val;

            best_ustar.push_back({pos, val});
        }

        optvectors >> nnz;
        best_xstar = SparseVector(n);
        while (nnz--) {
            int pos;
            double val;
            optvectors >> pos >> val;

            best_xstar.push_back({pos, val});
        }

        bestsol_avaiable = true;
    } else {
        bestsol_avaiable = false;
        best_xstar = best_ustar = SparseVector(0);
    }
}

SubgradientInstance::Status SubgradientInstance::solve(Methods method) {
    if (EXTRA) {
        std::cout << "\n\nSolving on model (Ax >= b formulation)" << std::endl;
        std::cout << "A " << A;
        std::cout << "b " << b;
        std::cout << "c " << c;
        std::cout << "\n\n";
    }

    switch (method) {
        case Methods::Pure:
            status = solvePure();
            break;
        case Methods::Deflected:
            status = solveDeflected();
            break;
        case Methods::Conditional:
            status = solveConditional();
            break;
        case Methods::Hybrid:
            status = solveHybrid();
            break;
    }

    switch (status) {
        case Status::NeverSolved:
        case Status::Error:
            assert(status != Status::NeverSolved);
            assert(status != Status::Error);
            break;

        case Status::LowLambdaSoft:
            std::cout << "Exit: lambda less than " << SOFT_LAMBDA_TRESHOLD
                      << " for more than " << SOFT_LAMBDA_MAXLIFE
                      << " iterations" << std::endl;
            break;

        case Status::LowLambdaHard:
            std::cout << "Exit: lambda less than " << HARD_LAMBDA_TRESHOLD
                      << std::endl;
            break;

        case Status::ReachedMaxIterations:
            std::cout << "Exit: reached max iterations (" << MAX_ITERATIONS
                      << ")" << std::endl;
            break;

        case Status::EpsilonOptimal:
            std::cout << "Exit: reached eps-optimality (eps = abs(ustar * s) = "
                      << gap << ")" << std::endl;
            saveSolutions();
            break;

        case Status::Optimal:
            std::cout << "Exit: reached optimality (d < 0)" << std::endl;
            saveSolutions();
            break;
    }

    return status;
}

// method variances
SubgradientInstance::Status SubgradientInstance::solvePure() {
    SparseVector u(m), s(m), x(n);
    const SparseVector zero(m);

    // step size params
    double lambda = -1.0;  // first iteration only
    double mu = 0.1;
    std::vector<double> dual_objectives(P);
    int soft_lambda_lifespan = 0;

    status = Status::NeverSolved;
    int k;
    for (k = 0; k < MAX_ITERATIONS; k++) {
        // compute n vector c - uT A
        SparseVector temp = c;
        temp -= u * A;

        // check i-th component: if > 0 set LB else UB
        x.data.clear();  // TODO(lugot): performance
        int i = 0;
        std::vector<std::pair<int, double>>::iterator it;
        for (it = temp.data.begin(); it != temp.data.end(); ++it) {
            while (i < it->first) {
                x.push_back({i, lb[i]});
                ++i;
            }

            if (it->second < EPSILON) {
                x.push_back({i, ub[i]});
            } else {
                // set to LB also if it->second is 0
                x.push_back({i, lb[i]});
            }
            ++i;
        }
        while (i < n) {
            x.push_back({i, lb[i]});
            ++i;
        }

        // compute the subgradient s as b - A x
        s = b - (A * x);

        // compute objs
        double dual_objective, primal_objective;
        dual_objectives[k % P] = dual_objective = c * x + u * s;
        primal_objective = c * x;

        // track the iteration status
        if (TRACKING) {
            std::cout << 'P' << " " << k << " " << dual_objective << " "
                      << primal_objective;
            if (bestsol_avaiable) {
                std::cout << " " << SparseVector::dist(u, ustar);
            }
            std::cout << " " << lambda << std::endl;
        }
        if (EXTRA) {
            std::cout << "u: " << u << std::endl;
            std::cout << "s: " << s << std::endl;
            std::cout << "x: " << x << std::endl;
        }

        // optimality check
        if (s < zero) {
            xstar = x;
            ustar = u;
            fstar = primal_objective;
            phistar = dual_objective;

            if ((gap = fabs(ustar * s)) < EPSILON) {
                return Status::Optimal;
            } else {
                return Status::EpsilonOptimal;
            }
            break;
        }

        // compute the step size lambda
        if (k > P) {
            double phimin = *std::min_element(dual_objectives.begin(),
                                              dual_objectives.end());
            double phimax = *std::max_element(dual_objectives.begin(),
                                              dual_objectives.end());

            if (phimax - 1.01 * phimin > -EPSILON) {
                // phimax - phimin > 0.01 * phimin
                mu = 0.5 * mu;
            } else if (phimax - 1.001 * phimin < EPSILON) {
                // phimax - phimin < 0.001 * phimin
                mu = std::min(2.0, 1.5 * mu);
            }
        }
        lambda = mu * (449.606 - dual_objective) / SparseVector::squaredNorm(s);

        // check if lambda is too small
        if (lambda < HARD_LAMBDA_TRESHOLD) {
            return Status::LowLambdaHard;
        } else if (lambda < SOFT_LAMBDA_TRESHOLD) {
            soft_lambda_lifespan++;
            if (soft_lambda_lifespan >= SOFT_LAMBDA_MAXLIFE) {
                return Status::LowLambdaSoft;
            } else {
                soft_lambda_lifespan = 0;
            }
        }

        // update the u vector by subgradient
        u += s * lambda;
        // project into R^m+, by removing those elements from the sparse vector
        u.data.erase(std::remove_if(u.data.begin(), u.data.end(),
                                    [](const std::pair<int, double>& p) {
                                        return p.second < EPSILON;
                                    }),
                     u.data.end());
    }
    if (k == MAX_ITERATIONS) {
        return Status::ReachedMaxIterations;
    }

    return Status::Error;
}

SubgradientInstance::Status SubgradientInstance::solveDeflected() {
    SparseVector u(m), s(m), d(m), x(n);
    const SparseVector zero(m);

    // step size params
    double lambda = -1.0;
    double mu = 0.1;
    std::vector<double> dual_objectives(P);
    int soft_lambda_lifespan = 0;

    status = Status::NeverSolved;
    int k;
    for (k = 0; k < MAX_ITERATIONS; k++) {
        // compute n vector c - uT A
        SparseVector temp = c;
        temp -= u * A;

        // check i-th component: if > 0 set LB else UB
        x.data.clear();  // TODO(lugot): performance
        int i = 0;
        std::vector<std::pair<int, double>>::iterator it;
        for (it = temp.data.begin(); it != temp.data.end(); ++it) {
            while (i < it->first) {
                x.push_back({i, lb[i]});
                ++i;
            }

            if (it->second < EPSILON) {
                x.push_back({i, ub[i]});
            } else {
                // set to LB also if it->second is 0
                x.push_back({i, lb[i]});
            }
            ++i;
        }
        while (i < n) {
            x.push_back({i, lb[i]});
            ++i;
        }

        // compute the subgradient s as b - A x as well as d
        s = b - (A * x);

        double delta_selector = s * d;
        double delta;
        if (delta_selector < EPSILON) {
            delta = -TAU * delta_selector;

            if (d.data.size() != 0) delta /= SparseVector::squaredNorm(d);

            // d = s + d * delta;
            d *= delta;
            d += s;
        } else {
            delta = 0;
            d = s;
        }

        // compute objs
        double dual_objective, primal_objective;
        dual_objectives[k % P] = dual_objective = c * x + u * s;
        primal_objective = c * x;

        // track the iteration status
        if (TRACKING) {
            std::cout << 'P' << " " << k << " " << dual_objective << " "
                      << primal_objective;
            if (bestsol_avaiable) {
                std::cout << " " << SparseVector::dist(u, ustar);
            }
            std::cout << " " << lambda << std::endl;
        }
        if (EXTRA) {
            std::cout << "u: " << u << std::endl;
            std::cout << "s: " << s << std::endl;
            std::cout << "d: " << d << std::endl;
            std::cout << "x: " << x << std::endl;
        }

        // optimality check
        if (s < zero) {
            xstar = x;
            ustar = u;
            fstar = primal_objective;
            phistar = dual_objective;

            if ((gap = fabs(ustar * s)) < EPSILON) {
                return Status::Optimal;
            } else {
                return Status::EpsilonOptimal;
            }
            break;
        }

        // compute the step size lambda
        if (k > P) {
            double phimin = *std::min_element(dual_objectives.begin(),
                                              dual_objectives.end());
            double phimax = *std::max_element(dual_objectives.begin(),
                                              dual_objectives.end());

            if (phimax - 1.01 * phimin > -EPSILON) {
                // phimax - phimin > 0.01 * phimin
                mu = 0.5 * mu;
            } else if (phimax - 1.001 * phimin < EPSILON) {
                // phimax - phimin < 0.001 * phimin
                mu = std::min(2.0, 1.5 * mu);
            }
        }
        lambda = mu * (500 - dual_objective) / SparseVector::squaredNorm(d);

        // check if lambda is too small
        if (lambda < HARD_LAMBDA_TRESHOLD) {
            return Status::LowLambdaHard;
        } else if (lambda < SOFT_LAMBDA_TRESHOLD) {
            soft_lambda_lifespan++;
            if (soft_lambda_lifespan >= SOFT_LAMBDA_MAXLIFE) {
                return Status::LowLambdaSoft;
            } else {
                soft_lambda_lifespan = 0;
            }
        }

        // update the u vector by subgradient
        u += d * lambda;
        // project into R^m+, by removing those elements from the sparse vector
        u.data.erase(std::remove_if(u.data.begin(), u.data.end(),
                                    [](const std::pair<int, double>& p) {
                                        return p.second < EPSILON;
                                    }),
                     u.data.end());
    }
    if (k == MAX_ITERATIONS) {
        return Status::ReachedMaxIterations;
    }

    return Status::Error;
}

SubgradientInstance::Status SubgradientInstance::solveConditional() {
    SparseVector u(m), s(m), cond_s(m), x(n);
    const SparseVector zero(m);

    // step size params
    double lambda = -1.0;
    double mu = 0.1;
    std::vector<double> dual_objectives(P);
    int soft_lambda_lifespan = 0;

    status = Status::NeverSolved;
    int k;
    for (k = 0; k < MAX_ITERATIONS; k++) {
        // compute n vector c - uT A
        SparseVector temp = c;
        temp -= u * A;

        // check i-th component: if > 0 set LB else UB
        x.data.clear();  // TODO(lugot): performance
        int i = 0;
        std::vector<std::pair<int, double>>::iterator it;
        for (it = temp.data.begin(); it != temp.data.end(); ++it) {
            while (i < it->first) {
                x.push_back({i, lb[i]});
                ++i;
            }

            if (it->second < EPSILON) {
                x.push_back({i, ub[i]});
            } else {
                // set to LB also if it->second is 0
                x.push_back({i, lb[i]});
            }
            ++i;
        }
        while (i < n) {
            x.push_back({i, lb[i]});
            ++i;
        }

        // compute the i-th comp of conditional subgradient s as
        // 0 if bi - A^i x <= 0 and  ui = 0, bi - A^i x otherwise
        // compute s as usual and then remove if satisfy condition
        s = b - (A * x);
        // TODO(lugot): SPEEDUP this is O(nlogn) -> O(n), n = sum(nonzeroes)
        s.data.erase(std::remove_if(s.data.begin(), s.data.end(),
                                    [&u](const std::pair<int, double>& p) {
                                        return p.second < EPSILON &&
                                               u.isZero(p.first);
                                    }),
                     s.data.end());

        // compute objs
        double dual_objective, primal_objective;
        dual_objectives[k % P] = dual_objective = c * x + u * s;
        primal_objective = c * x;

        // track the iteration status
        if (TRACKING) {
            std::cout << 'P' << " " << k << " " << dual_objective << " "
                      << primal_objective;
            if (bestsol_avaiable) {
                std::cout << " " << SparseVector::dist(u, ustar);
            }
            std::cout << " " << lambda << std::endl;
        }
        if (EXTRA) {
            std::cout << "u: " << u << std::endl;
            std::cout << "s: " << s << std::endl;
            std::cout << "x: " << x << std::endl;
        }

        // optimality check
        if (s < zero) {
            xstar = x;
            ustar = u;
            fstar = primal_objective;
            phistar = dual_objective;

            if ((gap = fabs(ustar * s)) < EPSILON) {
                return Status::Optimal;
            } else {
                return Status::EpsilonOptimal;
            }
            break;
        }

        // compute the step size lambda
        if (k > P) {
            double phimin = *std::min_element(dual_objectives.begin(),
                                              dual_objectives.end());
            double phimax = *std::max_element(dual_objectives.begin(),
                                              dual_objectives.end());

            if (phimax - 1.01 * phimin > -EPSILON) {
                // phimax - phimin > 0.01 * phimin
                mu = 0.5 * mu;
            } else if (phimax - 1.001 * phimin < EPSILON) {
                // phimax - phimin < 0.001 * phimin
                mu = std::min(2.0, 1.5 * mu);
            }
        }
        lambda = mu * (449.606 - dual_objective) / SparseVector::squaredNorm(s);

        // check if lambda is too small
        if (lambda < HARD_LAMBDA_TRESHOLD) {
            return Status::LowLambdaHard;
        } else if (lambda < SOFT_LAMBDA_TRESHOLD) {
            soft_lambda_lifespan++;
            if (soft_lambda_lifespan >= SOFT_LAMBDA_MAXLIFE) {
                return Status::LowLambdaSoft;
            } else {
                soft_lambda_lifespan = 0;
            }
        }

        // update the u vector by subgradient
        u += s * lambda;
        // project into R^m+, by removing those elements from the sparse vector
        size_t unnz = static_cast<int>(u.data.size());
        u.data.erase(std::remove_if(u.data.begin(), u.data.end(),
                                    [](const std::pair<int, double>& p) {
                                        return p.second < EPSILON;
                                    }),
                     u.data.end());
        if ((u.data.size() != unnz) && VERBOSE) {
            std::cout << "Zigzagging King II eliminated" << std::endl;
            // TODO(lugot): MODIFY
        }
    }
    if (k == MAX_ITERATIONS) {
        return Status::ReachedMaxIterations;
    }

    return Status::Error;
}

SubgradientInstance::Status SubgradientInstance::solveHybrid() {
    return status;
}

void SubgradientInstance::saveSolutions() {
    std::ofstream outfile;
    outfile.open("../subgradient_solutions/" + model_name);
    outfile << ustar << "\n" << xstar << "\n";
    outfile.close();
}

// helpers
/* double SubgradientInstance::phi(const SparseVector& u, const SparseVector& x)
{ return c * x + (b - (A * x)) * u;
} */
// double SubgradientInstance::obj(const SparseVector& x) { return c * x; }

// instance converter methods
void SubgradientInstance::extractConstraint(
    const IloRange& r, int row, std::vector<std::tuple<int, int, double>>* a,
    double* b) {
    // get number of non zero coefficents
    int nnz = 0;
    IloExpr::LinearIterator it;
    for (it = r.getLinearIterator(); it.ok(); ++it) nnz++;
    // TODO(lugot): seriously cplex?
    a->resize(nnz);  // TODO(lugot): reserve better

    // understand the sense of constraint, set a mutiplicative factor, we need
    // Ax >= b format
    Sense s = constraintSense(r);
    int mult = s == LESS ? -1 : 1;

    int i = 0;
    for (it = r.getLinearIterator(); it.ok(); ++it) {
        std::string varname = it.getVar().getName();
        if (varindexer.find(varname) == varindexer.end()) {
            varindexer[varname] = varindexer.size();
        }

        a->at(i) = {row, varindexer[varname], mult * it.getCoef()};

        i++;
    }

    if (s == LESS) {
        *b = -r.getUB();
    } else {
        *b = r.getLB();
    }
}

SubgradientInstance::Sense SubgradientInstance::constraintSense(
    const IloRange& r) {
    // aT x = b -> ub == lb
    if (r.getUB() == r.getLB()) return EQUAL;
    // aT x >= b -> ub == inf
    if (r.getUB() == IloInfinity) return GREATER;
    // aT x <= b -> lb == inf
    return LESS;
}

// Iteration tracker methods
/* SubgradientInstance::Iteration::Iteration(int k, const SparseVector& u,
                                          const SparseVector& x,
                                          const SparseVector& s,
                                          const SparseVector& d, double phi,
                                          double obj, double lambda, double mu,
                                          double delta) {
    this->k = k;
    this->u = u;
    this->x = x;
    this->s = s;
    this->d = d;
    this->phi = phi;
    this->obj = obj;
    this->lambda = lambda;
    this->mu = mu;
    this->delta = delta;
} */

/* std::ostream& operator<<(std::ostream& os,
                         const SubgradientInstance::Iteration& i) {
    os << "phi: ";
    if (!std::signbit(i.phi)) os << "+";
    os << std::defaultfloat << i.phi << " ";

    os << "obj: ";
    if (!std::signbit(i.obj)) os << "+";
    os << std::defaultfloat << i.obj << " ";

    os << "  [iteration " << i.k << "]";

    if (EXTRA) {
        os << "\n";
        os << "lambda: " << i.lambda << " mu: " << i.mu << " delta: " << i.delta
           << "\n";
        os << "x: " << i.x;
        os << "u: " << i.u;
        os << "s: " << i.s;
        os << "d: " << i.d;
    }

    return os;
} */
