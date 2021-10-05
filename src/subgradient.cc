#include "../include/subgradient.h"

#include <Eigen/Sparse>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <random>

#include "../include/globals.h"
#include "../include/logline.h"
#include "../include/step.h"

bool timelimitkill(std::chrono::steady_clock::time_point begin,
                   double timelimit) {
    std::chrono::steady_clock::time_point end;
    end = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
            .count();

    if (!IGNORETIMELIMIT && elapsed > timelimit) return true;
    return false;
}

Subgradient::Subgradient() {
    lambda = 0.0;
    mu = 0.2;
    f = INF;
    phi = -INF;
    phistar_UB = phibar = -INF;

    fstar = INF;
    phistar = -INF;
    gap = INF;

    iszig1 = false;
    iszig2 = false;

    alpha = 1;
    delta = 0;
}

void Subgradient::solve(const Problem& lp, double UB, double timelimit,
                        double mu0) {
    // std::cout << "solving problem on dim " << lp.n << " x " << lp.m
    //           << std::endl;

    this->mu = mu0;
    // std::cout << mu << std::endl;
    std::chrono::steady_clock::time_point begin, end;
    begin = std::chrono::steady_clock::now();

    if (LOGGING) {
        log.open("../results/" + lp.instance_name + ".subgr_log");

        LogLine line(',', false);
        std::vector<std::string> logline = {
            "method", "k",  "conditional", "deflection",  "dual", "primal",
            "lambda", "mu", "delta",       "lambda_life", "time"};

        line >> logline;
        // std::clog << line;
        log << line;
    }

    phistar_UB = UB;

    // mu evaluation
    // double up = mueval(lp, 1.0, begin, timelimit);
    // if (timelimitkill(begin, timelimit)) return;
    // if (VERBOSE) std::cout << "init mu=" << 1.0 << " -> " << up << std::endl;
    //
    // double middle = mueval(lp, 0.01, begin, timelimit);
    // if (timelimitkill(begin, timelimit)) return;
    // if (VERBOSE) std::cout << "init mu=" << 0.01 << " -> " << middle <<
    // std::endl;
    //
    // double down = mueval(lp, 0.0001, begin, timelimit);
    // if (timelimitkill(begin, timelimit)) return;
    // if (VERBOSE) std::cout << "init mu=" << 0.0001 << " -> " << down <<
    // std::endl;

    initialize(lp, begin, timelimit);

    end = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
            .count();
    LogLine line;

    // std::cout << "INITIAL PHI " << phi << std::endl;
    line >> 3 >> 0 >> iszig1 >> iszig2 >> phi >> f >> lambda >> mu >> delta >>
        0 >> elapsed;
    log << line;

    if (!IGNORETIMELIMIT && elapsed > timelimit) return;
    if (optimal) return;
    if (EXTRA) {
        std::cout << "vectors after initialization:" << std::endl;
        std::cout << "x: " << x << std::endl;
        std::cout << "s: " << s << std::endl;
        std::cout << "d: " << d << std::endl;
        std::cout << "u: " << u << std::endl;
        std::cout << std::endl;
    }

    int k = 1;
    while (iterate(lp, k, begin, timelimit)) {
        end = std::chrono::steady_clock::now();
        double elapsed =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
                .count();

        // if (k == 5) mu = 0.3;

        // track the iteration status
        LogLine line;
        switch (mode) {
            case Step::Mode::None:
                line >> 0;
                break;
            case Step::Mode::Deflect:
                line >> 1;
                break;
            case Step::Mode::Condition:
                line >> 2;
                break;
            case Step::Mode::Hybrid:
                line >> 3;
                break;
        }
        line >> k >> iszig2 >> iszig1 >> phi >> f >> lambda >> mu >> delta >>
            0 >> elapsed;
        log << line;

        if (!IGNORETIMELIMIT && elapsed > timelimit) {
            break;
        }

        if (LOGGING) {
            line.setDelim(',');
            line.setFancy(false);
            // std::clog << line;
        }
        if (!(k % ITERINFO)) {
            line.setFancy(true);
            std::cout << line;
        }
        if (EXTRA) {
            std::cout << "x: " << x << std::endl;
            std::cout << "s: " << s << std::endl;
            std::cout << "d: " << d << std::endl;
            std::cout << "u: " << u << std::endl;
        }

        k++;

        if (k == MAXITER) break;
        if (optimal) break;
    }

    // std::cout << "phistar: " << phistar << std::endl;
    // std::cout << "done " << k << " iterations" << std::endl;

    if (LOGGING) {
        log.close();
    }
}

double Subgradient::mueval(const Problem& lp, double mu,
                           std::chrono::steady_clock::time_point begin,
                           double timelimit) {
    // retreive problem
    SparseVector b = lp.b, c = lp.c;
    SparseMatrix A = lp.A;

    int n = lp.n;
    int m = lp.m;

    SparseVector u = SparseVector(m);
    SparseVector s = SparseVector(m);
    SparseVector x = SparseVector(n);

    solveSP(lp, &x, u, 0);

    std::chrono::steady_clock::time_point end;
    end = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
            .count();
    if (!IGNORETIMELIMIT && elapsed > timelimit) return false;

    // compute the subgradient s as b - A x
    s = b - (A * x);
    // compute objs
    double phi = c * x + u * s;

    SparseVector d = s;

    std::cout << phistar_UB << " " << phi << " " << SparseVector::squaredNorm(d)
              << std::endl;
    lambda = mu * (phistar_UB - phi) / SparseVector::squaredNorm(s);
    std::cout << lambda << std::endl;
    Step step = Step(lp, mode);
    std::tie(u, d) = step.performStep(s, d, u, lambda);

    project(lp, &u);

    x = Subgradient::solveSP(lp, u);

    std::cout << x << std::endl;
    return c * x + u * (b - A * x);
}

bool Subgradient::initialize(const Problem& lp,
                             std::chrono::steady_clock::time_point begin,
                             double timelimit) {
    // retreive problem
    SparseVector b = lp.b, c = lp.c;
    SparseMatrix A = lp.A;
    std::vector<bool> equal = lp.equal;

    int n = lp.n;
    int m = lp.m;

    u = SparseVector(m);
    s = SparseVector(m);
    d = SparseVector(m);
    x = SparseVector(n);

    fstar = INF;
    phistar = -INF;

    solveSP(lp, &x, u, 0);

    std::chrono::steady_clock::time_point end;
    end = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
            .count();
    if (!IGNORETIMELIMIT && elapsed > timelimit) return false;

    // compute the subgradient s as b - A x
    s = b - (A * x);
    // compute objs
    f = c * x;
    phis[0] = phi = f + u * s;

    d = s;

    iszig1 = Step::zigzagg1(s, d);
    iszig2 = Step::zigzagg2(s, u, lp);

    // optimality check
    if (s < zero) {
        optimal = true;
        std::cout << "SUBGR  OPTIMAL" << std::endl;
        return false;
    }

    // set phistar as this one
    fstar = f;
    phistar = phi;

    xstar = x;
    ustar = u;
    dstar = d;

    iszig1 = Step::zigzagg1(s, d);
    iszig2 = Step::zigzagg2(s, u, lp);

    lambda = mu * (phistar_UB - phi) / SparseVector::squaredNorm(d);
    Step step = Step(lp, mode);
    std::tie(u, d) = step.performStep(s, d, u, lambda);

    project(lp, &u);

    return true;
}

bool Subgradient::iterate(const Problem& lp, int k,
                          std::chrono::steady_clock::time_point begin,
                          double timelimit) {
    // retreive problem
    SparseVector b = lp.b, c = lp.c;
    SparseMatrix A = lp.A;
    std::vector<bool> equal = lp.equal;

    // solve the subproblem by putting variables to their upper/lower bounds
    solveSP(lp, &x, u, k);

    std::chrono::steady_clock::time_point end;
    end = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
            .count();
    if (!IGNORETIMELIMIT && elapsed > timelimit) return false;

    // compute the subgradient s as b - A x
    s = b - (A * x);
    // compute objs
    f = c * x;
    phis[k % P] = phi = f + u * s;

    // monotonicity check
    if (phi - phistar > -EPS) {
        // is monotone, update vectors

        fstar = f;
        phistar = phi;

        xstar = x;
        ustar = u;
        dstar = d;  // note it's not s! it's previous stepdir
    }

    // optimality check
    if (s == zero) {
        std::cout << "SUBGR  OPTIMAL" << std::endl;
        gap = 0.0;
        optimal = true;
        return false;
    }
    if (s < zero) {
        std::cout << "SUBGR  OPTIMAL" << std::endl;
        fstar = f;
        phistar = phi;

        xstar = x;
        ustar = u;
        dstar = d;

        gap = fabs(ustar * s);
        optimal = true;
        return false;
    }

    iszig1 = Step::zigzagg1(s, d);
    iszig2 = Step::zigzagg2(s, u, lp);

    // if (k >= P) {
    if (k >= P && k % P == 0) {
        mu = muKarp(phi, phis, mu);
        // mu = muCaprara(phis, mu);
    }

    alpha = 1;
    phibar = alpha * phistar_UB + (1 - alpha) * phi;
    lambda = mu * (phibar - phi) / SparseVector::squaredNorm(d);
    // if (lambda < 1e-15) return false;

    SparseVector oldu = u;

    Step step = Step(lp, mode);
    std::tie(u, d) = step.performStep(s, d, u, lambda);

    project(lp, &u);

    return true;
}

SparseVector* Subgradient::solveSP(const Problem& lp, SparseVector* x,
                                   const SparseVector& u, int k) {
    std::chrono::steady_clock::time_point begin, end;

    // compute n vector c - uT A
    begin = std::chrono::steady_clock::now();
    SparseVector temp = lp.c;
    temp -= const_cast<SparseVector&>(u) * lp.A;
    temp.prune(HIGHEPS);

    // check i-th component: if > 0 set LB else UB
    x->data.clear();
    int i = 0;
    std::vector<std::pair<int, double>>::iterator it;
    for (it = temp.data.begin(); it != temp.data.end(); ++it) {
        while (i < it->first) {
            double topush = lp.lb[i];
            // topush += (lp.ub[i] - lp.lb[i]) / 4;

            x->push_back({i, topush});
            ++i;
        }

        if (it->second < EPS) {
            x->push_back({i, lp.ub[i]});
        } else {
            // set to lp.getLowerBounds() also if it->second is 0
            x->push_back({i, lp.lb[i]});
        }
        ++i;
    }
    while (i < lp.n) {
        double topush = lp.lb[i];
        // topush += (lp.ub[i] - lp.lb[i]) / 4;
        x->push_back({i, topush});

        ++i;
    }

    end = std::chrono::steady_clock::now();
    if (EXTRA) {
        double elapsed =
            std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
                .count();
        // std::cout << "SolveSP(u) time: " << elapsed << "[µs]" << std::endl;
    }

    return x;
}

SparseVector Subgradient::solveSP(const Problem& lp, const SparseVector& u) {
    SparseVector x = SparseVector(lp.n);
    solveSP(lp, &x, u, 1000);

    return x;
}

Subgradient::DirectionStatus Subgradient::directionInfeasible(
    const Problem& lp, const SparseVector& u, const SparseVector& s) {
    std::vector<std::pair<int, double>>::const_iterator itu, its;
    itu = u.data.begin();
    its = s.data.begin();
    while (itu != u.data.end() && its != s.data.end()) {
        if (itu->first == its->first) {
            ++itu;
            ++its;
        } else if (itu->first < its->first) {
            ++itu;
        } else {
            if (its->second < EPS && !lp.equal[its->first]) {
                // here if ui is zero and si < 0, need to condtion
                std::cout << itu->first << " " << its->first << " "
                          << its->second << std::endl;
                IMHERE
                return Subgradient::DirectionStatus::Conditional;
            }
            ++its;
        }
    }

    while (its != s.data.end()) {
        if (its->second < EPS && !lp.equal[its->first]) {
            // here if ui is zero and si < 0, need to condtion
            return Subgradient::DirectionStatus::Conditional;
        }

        ++its;
    }
    return Subgradient::DirectionStatus::Deflected;
}

void Subgradient::project(const Problem& lp, SparseVector* u) {
    // project into R^m+, by removing those elements from the sparse vector
    u->data.erase(std::remove_if(u->data.begin(), u->data.end(),
                                 [&lp](const std::pair<int, double>& p) {
                                     return p.second < EPS &&
                                            !lp.equal[p.first];
                                 }),
                  u->data.end());
}


double Subgradient::muKarp(double phi, std::vector<double> phis, double mu) {
    std::cout << phi << std::endl;
    if (phi < *std::max_element(phis.begin(), phis.end())) {
        return mu / 2;
    }
    return mu;
}

double Subgradient::muCaprara(std::vector<double> phis, double mu) {
    double phimin = *std::min_element(phis.begin(), phis.end());
    double phimax = *std::max_element(phis.begin(), phis.end());

    // for (double phi: phis) std::cout << phi << " ";
    // std::cout << std::endl;
    // std::cout << "max, min: " << phimax << " " << phimin << std::endl;
    // std::cout << "mu: " << mu << std::endl;

    if (phimax - 1.01 * phimin > -EPS) {
        // phimax - phimin > 0.01 * phimin
        mu = std::max(MU_MIN, 0.5 * mu);
    } else if (phimax - 1.001 * phimin < EPS) {
        // phimax - phimin < 0.001 * phimin
        mu = std::min(MU_MAX, 1.5 * mu);
    }

    // std::cout << "mu: " << mu << std::endl;

    return mu;
}

// CPLEX SOLVE SP
// SparseVector* Subgradient::solveSP(const Problem& lp, SparseVector* x,
//                                    const SparseVector& u) {
//     std::chrono::steady_clock::time_point begin, end;
//     begin = std::chrono::steady_clock::now();
//
//     Eigen::SparseVector<double> eb, ec, ex, eu;
//     Eigen::SparseMatrix<double> eA;
//
//     eb = lp.b.toEigen();
//     ec = lp.c.toEigen();
//     eA = lp.A.toEigen();
//     eu = u.toEigen();
//     ex = x->toEigen();
//
//     Eigen::SparseVector<double> etemp = ec.transpose() - eu.transpose() * eA;
//     SparseVector temp2 = SparseVector(etemp);
//     temp2.prune(HIGHEPS);
//
//     end = std::chrono::steady_clock::now();
//     std::cout << "eigen_time: "
//               << std::chrono::duration_cast<std::chrono::microseconds>(end -
//                                                                        begin)
//                      .count()
//               << "[µs]" << std::endl;
//
//     // compute n vector c - uT A
//     begin = std::chrono::steady_clock::now();
//     SparseVector temp = lp.c;
//     temp -= const_cast<SparseVector&>(u) * lp.A;
//     temp.prune(HIGHEPS);
//     end = std::chrono::steady_clock::now();
//     std::cout << "mine_time: "
//               << std::chrono::duration_cast<std::chrono::microseconds>(end -
//                                                                        begin)
//                      .count()
//               << "[µs]" << std::endl;
//
//     assert(temp == temp2);
//
//     // check i-th component: if > 0 set LB else UB
//     begin = std::chrono::steady_clock::now();
//     x->data.clear();
//     int i = 0;
//
//     std::vector<std::pair<int, double>>::iterator it;
//     for (it = temp.data.begin(); it != temp.data.end(); ++it) {
//         while (i < it->first) {
//             x->push_back({i, lp.lb[i]});
//             ++i;
//         }
//
//         if (it->second < EPS) {
//             x->push_back({i, lp.ub[i]});
//         } else {
//             // set to lp.getLowerBounds() also if it->second is 0
//             x->push_back({i, lp.lb[i]});
//         }
//         ++i;
//     }
//     while (i < lp.n) {
//         x->push_back({i, lp.lb[i]});
//         ++i;
//     }
//     end = std::chrono::steady_clock::now();
//     std::cout << "mine_time: "
//               << std::chrono::duration_cast<std::chrono::microseconds>(end -
//                                                                        begin)
//                      .count()
//               << "[µs]" << std::endl;
//
//     // std::cout << "x: " << *x << std::endl;
//
//     begin = std::chrono::steady_clock::now();
//     IloEnv env;
//     IloModel model(env);
//     IloNumVarArray var(env);
//     IloRangeArray con(env);
//
//
//     for (int i = 0; i < lp.n; ++i) var.add(IloNumVar(env, lp.lb[i],
//     lp.ub[i]));
//
//     IloObjective obj = IloMinimize(env);
//     std::vector<std::pair<int, double>>::iterator it2;
//     for (it2 = temp.data.begin(); it2 != temp.data.end(); ++it2) {
//         obj.setLinearCoef(var[it2->first], it2->second);
//     }
//
//     for (int i = 0; i < lp.n; ++i) {
//         con.add(var[i] <= lp.ub[i]);
//     }
//     model.add(obj);
//     model.add(con);
//
//     IloCplex cplex(model);
//     cplex.solve();
//
//     SparseVector cpxx = SparseVector(lp.n);
//     for (int i = 0; i < lp.n; ++i) {
//         double xi = cplex.getValue(var[i]);
//         cpxx.push_back({i, xi});
//     }
//     end = std::chrono::steady_clock::now();
//     std::cout << "cplex_time: "
//               << std::chrono::duration_cast<std::chrono::microseconds>(end -
//                                                                        begin)
//                      .count()
//               << "[µs]" << std::endl;
//
//     // std::cout << "cplex: " << cpxx << std::endl;
//     // std::cout << "mine: " << *x << std::endl;
//
//     assert(cpxx == *x);
//
//
//
//     return x;
// }
