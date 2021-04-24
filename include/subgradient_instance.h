#ifndef INCLUDE_SUBGRADIENT_INSTANCE_H_
#define INCLUDE_SUBGRADIENT_INSTANCE_H_

#include <Eigen/Sparse>

#include "../include/cplex_instance.h"

class SubgradientInstance {
   public:
    SubgradientInstance(CplexInstance& cinst);

   private:
    Eigen::SparseMatrix<double> A;
    Eigen::SparseVector<double> b;
    Eigen::SparseVector<double> c;
};

#endif  // INCLUDE_SUBGRADIENT_INSTANCE_H_
