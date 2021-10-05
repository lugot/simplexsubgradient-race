#ifndef INCLUDE_SPARSE_VECTOR_H_
#define INCLUDE_SPARSE_VECTOR_H_

#include <math.h>

#include <Eigen/Sparse>
#include <string>
#include <utility>
#include <vector>

class SparseMatrix;

class SparseVector {
   public:
    SparseVector() = default;
    SparseVector(const SparseVector& x) {
        this->size = x.size;
        this->data = x.data;
    }
    explicit SparseVector(int size) { this->size = size; }
    explicit SparseVector(Eigen::SparseVector<double> x);
    Eigen::SparseVector<double> toEigen() const;
    SparseVector operator=(const SparseVector& x) {
        this->size = x.size;
        this->data = x.data;
        return *this;
    }
    ~SparseVector() {}

    // void push_back(std::pair<int, double> p);
    void push_back(std::pair<int, double> p, double tollerance = 0.0) {
        if (tollerance == 0.0 && p.second == 0.0) return;
        if (tollerance > 0.0 && fabs(p.second) < tollerance) return;

        data.push_back(p);
    }
    // void push_back(std::pair<int, double> p, double tollerance);
    // void push_back(std::pair<int, double> p);

    int getSize() { return this->size; }
    int nnz() { return this->data.size(); }

    SparseVector operator+(const SparseVector& x) const;
    SparseVector operator+=(const SparseVector& x);
    SparseVector operator-(const SparseVector& x);
    SparseVector operator-=(const SparseVector& x);
    double operator*(const SparseVector& x) const;
    SparseVector operator*(const SparseMatrix& A);
    SparseVector operator*(double scalar) const;
    SparseVector operator*=(double scalar);
    bool operator<(const SparseVector& x);
    bool operator>(const SparseVector& x);
    bool operator==(const SparseVector& x) const;

    void prune(double tollerance);
    bool isZero(int idx) const;
    int nonZeroes();

    static double squaredNorm(const SparseVector& x);
    static double norm(const SparseVector& x);
    static double dist(const SparseVector& x, const SparseVector& y);
    static SparseVector merge(const SparseVector& x, const SparseVector& y);

    friend std::ostream& operator<<(std::ostream& os, const SparseVector& x);

    // actual data
    int size;
    std::vector<std::pair<int, double>> data;
};

#endif  // INCLUDE_SPARSE_VECTOR_H_
