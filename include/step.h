#ifndef INCLUDE_STEP_H_
#define INCLUDE_STEP_H_

#include "../include/globals.h"
#include "../include/problem.h"
#include "../include/sparse_vector.h"

class Step {
   public:
    enum class Mode { None, Deflect, Condition, Hybrid };
    explicit Step(const Problem& lp, Mode mode = Step::Mode::Hybrid);
    std::pair<SparseVector, SparseVector> performStep(const SparseVector& s, const SparseVector& d,
                           const SparseVector& u, double lambda);

    double getDelta() { return delta; }
    int getStatus() { return status; }

    static bool zigzagg1(const SparseVector& s, const SparseVector& u);
    static bool zigzagg2(const SparseVector& s, const SparseVector& u, const Problem& lp);

   private:
    static double deltaCamerini(const SparseVector& d, const SparseVector& s);
    static double deltaSherali(const SparseVector& d, const SparseVector& s);
    static double deltaBelgachem(
        const SparseVector& d,
        const SparseVector&
            s);  
    // static double alphaRandomCombination(SparseVector* d,
    //                                      const SparseVector& s);

    void feasibilize(SparseVector* g, const SparseVector& u);

    Problem lp;
    double delta;
    int status;
    Mode mode;
    static const double TAU;
};

#endif  // INCLUDE_STEP_H_
