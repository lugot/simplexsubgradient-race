#include "../include/step.h"

#include <../include/subgradient.h>

#include <algorithm>
#include <random>

#include "../include/globals.h"

const double Step::TAU = 1.5;

Step::Step(const Problem& lp, Mode mode) {
    this->lp = lp;
    this->mode = mode;
}

std::pair<SparseVector, SparseVector> Step::performStep(const SparseVector& s,
                                                        const SparseVector& d,
                                                        const SparseVector& u,
                                                        double lambda) {
    if (mode == Step::Mode::Hybrid && !zigzagg1(s, d)) {
        mode = Step::Mode::Condition;
    }
    if (mode == Step::Mode::Hybrid && !zigzagg2(s, u, lp)) {
        mode = Step::Mode::Deflect;
    }

    switch (mode) {
        case Mode::None: {
            return {u + s * lambda, s};
        }
        case Mode::Condition: {
            if (zigzagg2(s, u, lp)) {
                SparseVector d = s;
                feasibilize(&d, u);

                return {u + d * lambda, d};
            }
        }
        case Mode::Deflect: {
            double delta = deltaBelgachem(d, s);
            SparseVector dir = s + d * delta;
            return {u + dir * lambda, dir};
        }
        case Mode::Hybrid: {
            SparseVector bars = s;
            feasibilize(&bars, u);

            SparseVector bard = d;
            feasibilize(&bard, u);

            double delta = deltaBelgachem(bard, bars);
            SparseVector dir = s + d * delta;
            feasibilize(&dir, u);
            return {u + dir * lambda, dir};
        }
    }

    return {u, s};
}

double Step::deltaCamerini(const SparseVector& d, const SparseVector& s) {
    double delta;

    double deflection_selector = d * s;
    std::cout << "selector: " << deflection_selector << std::endl;
    if (deflection_selector < EPS) {
        delta = -Step::TAU * deflection_selector;

        if (d.data.size() != 0) {
            delta /= SparseVector::squaredNorm(d);
        }
    } else {
        delta = 0.0;
    }

    return delta;
}
double Step::deltaSherali(const SparseVector& d, const SparseVector& s) {
    return SparseVector::norm(s) / SparseVector::norm(d);
}
double Step::deltaBelgachem(const SparseVector& d, const SparseVector& s) {
    double delta;
    double deflection_selector = d * s;
    if (deflection_selector < EPS) {
        double alpha = -std::cos(deflection_selector);
        double tau = std::min(Step::TAU, 1 / (2 - alpha) - 1e-4);  // eps: 1e-4

        delta = -tau * (1 - alpha) * (deflection_selector);
        delta += alpha * SparseVector::norm(d) * SparseVector::norm(s);

        if (d.data.size() != 0) {
            delta /= SparseVector::squaredNorm(d);
        }
    } else {
        delta = 0.0;
    }

    return delta;
}
// double Step::alphaRandomCombination(SparseVector* d, const SparseVector& s) {
//     double alpha;
//
//     double deflection_selector = (*d) * s;
//     if (deflection_selector < EPS) {
//         double stepdir_sqnorm = SparseVector::squaredNorm(*d);
//
//         double alpha_UB =
//             std::min(1.0, std::sqrt(stepdir_sqnorm) /
//                               (stepdir_sqnorm - deflection_selector));
//
//         std::random_device rd;
//         std::default_random_engine eng(rd());
//         std::uniform_real_distribution<double> distr(0.0, alpha_UB);
//
//         alpha = distr(eng);
//     } else {
//         alpha = 1.0;
//     }
//
//     return alpha;
// }

void Step::feasibilize(SparseVector* g, const SparseVector& u) {
    g->data.erase(
        std::remove_if(g->data.begin(), g->data.end(),
                       [this, u](const std::pair<int, double>& p) {
                           return !this->lp.equal[p.first] &&
                                  (p.second < EPS && u.isZero(p.first));
                       }),
        g->data.end());
}

bool Step::zigzagg1(const SparseVector& s, const SparseVector& d) {
    return (d * s < EPS);
}

bool Step::zigzagg2(const SparseVector& s, const SparseVector& u,
                    const Problem& lp) {
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
                return true;
            }
            ++its;
        }
    }
    while (its != s.data.end()) {
        if (its->second < EPS && !lp.equal[its->first]) {
            // here if ui is zero and si < 0, need to condtion
            return true;
        }

        ++its;
    }
    return false;
}

// OLD HYBRID
// SparseVector ans, dirans;
// double bestphi = -INF;
// int bestscheme = 0;
//
// std::vector<double> objs(8);
//
// for (int scheme = 0; scheme < 8; scheme++) {
//     if (EXTRA) {
//         std::cout << std::endl << "scheme" << scheme << std::endl;
//     }
//     SparseVector bars = s, bard = d;
//
//     if (scheme & 1) {
//         if (EXTRA) std::cout << "cleannup on s: ";
//         feasibilize(&bars, u);
//     }
//     if (scheme & 2) {
//         if (EXTRA) std::cout << "cleanup on d -1: ";
//         feasibilize(&bard, u);
//     }
//
//     double delta = deltaBelgachem(bard, bars);
//     SparseVector dir = bars + bard * delta;
//
//     if (scheme & 4) {
//         if (EXTRA) std::cout << "cleanup on d: ";
//         feasibilize(&dir, u);
//     }
//
//     if (scheme == 0) {
//         SparseVector newu = u + dir * lambda;
//         Subgradient::project(lp, &newu);
//         SparseVector x = Subgradient::solveSP(lp, newu);
//         double phi = lp.c * x + newu * (lp.b - lp.A * x);
//         objs[scheme] = phi;
//
//         ans = u + dir * lambda;
//         dirans = dir;
//     } else {
//         SparseVector newu = u + dir * lambda;
//         Subgradient::project(lp, &newu);
//         SparseVector x = Subgradient::solveSP(lp, newu);
//         double phi = lp.c * x + newu * (lp.b - lp.A * x);
//         objs[scheme] = phi;
//         if (EXTRA) std::cout << "phi: " << phi << std::endl;
//
//         if (phi > bestphi) {
//             ans = newu;
//             dirans = dir;
//             bestphi = phi;
//             bestscheme = scheme;
//         }
//     }
// }
//
// std::ofstream schemes;
// schemes.open("schemes", std::ios_base::app);
// // schemes << lp.instance_name << ",";
// for (double phi : objs) schemes << phi << ",";
// schemes << std::endl;
// schemes.close();
//
// if (EXTRA) {
//     std::cout << "bestphi: " << bestphi
//               << " bestscheme: " << bestscheme << std::endl;
// }
//
// return {ans, dirans};
