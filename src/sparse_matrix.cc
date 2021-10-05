#include "../include/sparse_matrix.h"

#include <cassert>
#include <iostream>

#include "../include/globals.h"

SparseMatrix::SparseMatrix(int nrows, int ncols) {
    this->nrows = nrows;
    this->ncols = ncols;

    this->rows = std::vector<SparseVector>(nrows, SparseVector(ncols));
    // for (int i = 0; i < nrows; ++i) this->rows[i] = SparseVector(ncols);
    this->cols = std::vector<SparseVector>(ncols, SparseVector(nrows));
    // for (int j = 0; j < ncols; ++j) this->cols[j] = SparseVector(nrows);
}

Eigen::SparseMatrix<double> SparseMatrix::toEigen() const {
    Eigen::SparseMatrix<double> ret = Eigen::SparseMatrix<double>(nrows, ncols);

    std::cout << "dimensions " << nrows << " " << ncols << std::endl;

    std::vector<Eigen::Triplet<double>> triplet_list;
    for (int i = 0; i < nrows; ++i) {
        std::vector<std::pair<int, double>>::const_iterator it;
        for (it = rows[i].data.begin(); it != rows[i].data.end(); ++it) {
            triplet_list.push_back(
                Eigen::Triplet<double>(i, it->first, it->second));
        }
    }

    ret.setFromTriplets(triplet_list.begin(), triplet_list.end());

    return ret;
}

SparseVector SparseMatrix::operator*(const SparseVector& x) {
    assert(ncols == x.size);

    SparseVector ans(nrows);
    for (int i = 0; i < nrows; ++i) ans.push_back({i, rows[i] * x});

    return ans;
}

void SparseMatrix::addRow(
    const std::vector<std::tuple<int, int, double>>& row) {
    std::vector<std::tuple<int, int, double>>::const_iterator it;

    for (it = row.begin(); it != row.end(); ++it) {
        int i = std::get<0>(*it);
        int j = std::get<1>(*it);
        double data = std::get<2>(*it);

        rows[i].push_back({j, data});
        cols[j].push_back({i, data});
    }
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix& m) {
    os << "SparseMatrix dim " << m.nrows << " x " << m.ncols << "\n";
    if (VERBOSE) {
        for (int i = 0; i < m.nrows; i++) {
            os << "Row" << i + 1 << ": " << m.rows[i] << std::endl;
        }
        for (int j = 0; j < m.ncols; j++) {
            os << "Col" << j + 1 << ": " << m.cols[j] << std::endl;
        }
    }

    return os;
}
