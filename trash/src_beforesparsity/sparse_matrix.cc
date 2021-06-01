#include "../include/sparse_matrix.h"

#include "../include/globals.h"

SparseMatrix::SparseMatrix(int nrows, int ncols) {
    this->nrows = nrows;
    this->ncols = ncols;

    this->rows = std::vector<SparseVector>(nrows, SparseVector(ncols));
    // for (int i = 0; i < nrows; ++i) this->rows[i] = SparseVector(ncols);
    this->cols = std::vector<SparseVector>(ncols, SparseVector(nrows));
    // for (int j = 0; j < ncols; ++j) this->cols[j] = SparseVector(nrows);
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
            os << "Row" << i + 1 << ": " << m.rows[i];
        }
        for (int j = 0; j < m.ncols; j++) {
            os << "Col" << j + 1 << ": " << m.cols[j];
        }
    }

    return os;
}
