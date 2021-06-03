#ifndef INCLUDE_SPARSE_VECTOR_H_
#define INCLUDE_SPARSE_VECTOR_H_

#include <string>
#include <vector>

// forward-declaration TODO(lugot): LEARN
class SparseMatrix;

class SparseVector {
   public:
    SparseVector() = default;
    explicit SparseVector(int size);

    void push_back(std::pair<int, double> p);
    void push_back(std::pair<int, double> p, double tollerance);

    int getSize();

    /* std::vector<std::pair<int, double>>::iterator begin();
    std::vector<std::pair<int, double>>::iterator end(); */

    SparseVector operator+(const SparseVector& x);
    SparseVector operator+=(const SparseVector& x);
    SparseVector operator-(const SparseVector& x);
    SparseVector operator-=(const SparseVector& x);
    double operator*(const SparseVector& x);
    SparseVector operator*(const SparseMatrix& A);
    SparseVector operator*(double scalar);
    SparseVector operator*=(double scalar);
    bool operator<(const SparseVector& x);
    bool operator==(const SparseVector& x);

    void prune(double tollerance);
    bool isZero(int idx);

    static double squaredNorm(const SparseVector& x);
    static double dist(const SparseVector& x, const SparseVector& y);

    friend std::ostream& operator<<(std::ostream& os, const SparseVector& x);

    // actual data
    int size;
    std::vector<std::pair<int, double>> data;
    // bool positive;

    /* private:
     SparseVector sum(const SparseVector& x);
     SparseVector sumAndProject(const SparseVector& x);
     SparseVector accumulate(const SparseVector& x);
     SparseVector accumulateAndProject(const SparseVector& x); */
};

#endif  // INCLUDE_SPARSE_VECTOR_H_
