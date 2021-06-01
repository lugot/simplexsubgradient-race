#ifndef INCLUDE_SPARSE_MATRIX_H_
#define INCLUDE_SPARSE_MATRIX_H_

#include <ostream>
#include <tuple>
#include <vector>

#include "../include/sparse_vector.h"

class SparseMatrix {
   public:
    SparseMatrix() = default;
    explicit SparseMatrix(int nrows, int ncols);

    SparseVector operator*(const SparseVector& x);

    void addRow(const std::vector<std::tuple<int, int, double>>& row);

    int nrows, ncols;
    std::vector<SparseVector> rows, cols;

    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix& A);
};


#endif  // INCLUDE_SPARSE_MATRIX_H_
