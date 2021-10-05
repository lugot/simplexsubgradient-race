#ifndef INCLUDE_SPARSE_MATRIX_H_
#define INCLUDE_SPARSE_MATRIX_H_

#include <ostream>
#include <tuple>
#include <vector>
#include <Eigen/Sparse>

#include "../include/sparse_vector.h"

class SparseMatrix {
   public:
    SparseMatrix() = default;
    explicit SparseMatrix(int nrows, int ncols);
    Eigen::SparseMatrix<double> toEigen() const;

    SparseVector operator*(const SparseVector& x);

    void addRow(const std::vector<std::tuple<int, int, double>>& row);

    int nrows, ncols;
    std::vector<SparseVector> rows, cols;

    int nnz() {
        int ans = 0;
        for (SparseVector& row: rows) ans += row.nnz();

        return ans;
    }

    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix& A);
};


#endif  // INCLUDE_SPARSE_MATRIX_H_
