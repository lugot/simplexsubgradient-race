#ifndef INCLUDE_PROBLEM_H_
#define INCLUDE_PROBLEM_H_

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "../include/cplex.h"
#include "../include/sparse_matrix.h"
#include "../include/sparse_vector.h"

class Problem {
   public:
    friend class Subgradient;
    friend class Step;
    Problem() = default;
    explicit Problem(const Cplex& cplex_instance);

    std::string instance_name;

    int nnz() { return b.nnz() + A.nnz(); }

   private:
    int m, n;
    SparseVector c, b;
    SparseMatrix A;
    std::vector<double> lb, ub;
    std::vector<bool> equal;

    // instance converter methods
    enum Sense { GREATER, EQUAL, LESS };
    std::map<std::string, int> varindexer;
    void extractConstraint(const IloRange& r, int row,
                           std::vector<std::tuple<int, int, double>>* a,
                           double* b);
    Sense constraintSense(const IloRange& r);
};

#endif  // INCLUDE_PROBLEM_H_
