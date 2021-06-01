#include "../include/sparse_matrix.h"

#include "../include/globals.h"

Node::Node(int row, int col, double data) {
    this->row = row;
    this->col = col;
    this->data = data;

    this->next[0] = nullptr;
    this->next[1] = nullptr;
}

SparseMatrix::SparseMatrix(int nrows, int ncols) {
    this->nrows = nrows;
    this->ncols = ncols;

    this->rows = std::vector<Node*>(nrows);
    for (int i = 0; i < nrows; ++i) this->rows[i] = new Node(i, -1, 0.0);
    this->cols = std::vector<Node*>(ncols);
    for (int j = 0; j < ncols; ++j) this->cols[j] = new Node(-1, j, 0.0);
}

void SparseMatrix::addRow(
    const std::vector<std::tuple<int, int, double>>& row) {
    std::vector<std::tuple<int, int, double>>::const_iterator it;

    for (it = row.begin(); it != row.end(); ++it) {
        int i = std::get<0>(*it);
        int j = std::get<1>(*it);
        double data = std::get<2>(*it);

        Node* n = new Node(i, j, data);

        Node* row_end = advanceRow(i);  // TODO(lugot): CAN BE AVOIDED
        Node* col_end = descendColumn(j);

        row_end->next[0] = n;
        col_end->next[1] = n;
    }
}
Node* SparseMatrix::advanceRow(int i) {
    Node* n = rows[i];
    while (n->next[0] != nullptr) n = n->next[0];

    return n;
}

Node* SparseMatrix::descendColumn(int j) {
    Node* n = cols[j];
    while (n->next[1] != nullptr) n = n->next[1];

    return n;
}

std::ostream& operator<<(std::ostream& os, const SparseMatrix& m) {
    os << "SparseMatrix dim " << m.nrows << " x " << m.ncols << "\n";
    if (VERBOSE) {
        for (int i = 0; i < m.nrows; i++) {
            os << "Row" << i + 1 << ": ";
            Node* n = m.rows[i]->next[0];

            while (n != nullptr) {
                os << n->data << " ";
                n = n->next[0];
            }
            os << "\n";
        }
        os << "\n";
        for (int i = 0; i < m.ncols; i++) {
            os << "Col" << i + 1 << ": ";
            Node* n = m.cols[i]->next[1];

            while (n != nullptr) {
                os << n->data << " ";
                n = n->next[1];
            }
            os << "\n";
        }
    }

    return os;
}
