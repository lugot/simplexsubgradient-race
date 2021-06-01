#ifndef INCLUDE_SPARSE_MATRIX_H_
#define INCLUDE_SPARSE_MATRIX_H_

#include <ostream>
#include <tuple>
#include <vector>

class Node {
   public:
    explicit Node(int row, int col, double data);

    int row, col;
    double data;
    Node* next[2];
};

class SparseMatrix {
   public:
    explicit SparseMatrix(int nrows, int ncols);
    void addRow(const std::vector<std::tuple<int, int, double>>& row);

    int nrows, ncols;
    std::vector<Node*> rows, cols;

    friend std::ostream& operator<<(std::ostream& os, const SparseMatrix& A);

   private:
    Node* advanceRow(int i);
    Node* descendColumn(int j);
};


#endif  // INCLUDE_SPARSE_MATRIX_H_
