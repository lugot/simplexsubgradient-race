#ifndef INCLUDE_SUBGRADIENT_H_
#define INCLUDE_SUBGRADIENT_H_

#include <chrono>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "../include/globals.h"
#include "../include/problem.h"
#include "../include/step.h"

class Subgradient {
   public:
    friend class Step;

    enum class ObjectiveStatus { MonotoneIncreasing, Increasing, Decreasing };

    enum class DirectionStatus { Deflected, Conditional };

    Subgradient();
    void solve(const Problem& lp, double UB = INF, double timelimit = 1000, double mu = 0.3);

   private:
    bool optimal = false;
    double mueval(const Problem& lp, double mu,
                    std::chrono::steady_clock::time_point begin,
                    double timelimit);
    bool initialize(const Problem& lp,
                    std::chrono::steady_clock::time_point begin,
                    double timelimit);
    bool iterate(const Problem& lp, int k,
                 std::chrono::steady_clock::time_point begin, double timelimit);

    // iteration information
    // u -> lagrangian multipliers
    // s -> subgradient of L(u)
    // d -> stepdirection
    // x -> primal variables
    SparseVector u, s, d, x;
    // lambda -> stepsize
    // mu -> stepsize control parameter
    // delta -> deflection parameter
    double lambda, mu, delta;
    // f -> primal objective
    // phi -> dual objective
    // phistar_UB -> upper bound of phi star
    // phibar -> variable target for phistar UB
    double f, phi, phistar_UB, phibar;
    // z -> stucked counter
    // v -> stucked counter
    // alpha -> fixed upperbound vs actual phistar tradeoff parameter
    double alpha;

    std::vector<double> phis = std::vector<double>(P);
    std::vector<double> mean_phis = std::vector<double>(P);

    // constants
    const SparseVector zero = SparseVector(0);

    // best one
    SparseVector ustar, dstar, xstar;
    double phistar, gap, fstar;
    bool iszig1, iszig2;

    Step::Mode mode = Step::Mode::Hybrid;

    static SparseVector* solveSP(const Problem& lp, SparseVector* x,
                                 const SparseVector& u, int k);
    static SparseVector solveSP(const Problem& lp, const SparseVector& u);
    static DirectionStatus directionInfeasible(const Problem& lp,
                                               const SparseVector& u,
                                               const SparseVector& s);
    static void project(const Problem& lp, SparseVector* u);

    static double muKarp(double phi, std::vector<double> phis, double mu);
    static double muCaprara(std::vector<double> phis, double mu);

    std::ofstream log;
};

#endif  // INCLUDE_SUBGRADIENT_H_
