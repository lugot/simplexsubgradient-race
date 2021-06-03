#include "../include/sparse_vector.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

#include "../include/globals.h"
#include "../include/sparse_matrix.h"

SparseVector::SparseVector(int size) {
    this->size = size;
    // this->positive = false;
}
/* SparseVector::SparseVector(int size, bool positive) {
    this->size = size;
    this->positive = positive;
} */

void SparseVector::push_back(std::pair<int, double> p) {
    if (p.second == 0.0) return;

    data.push_back(p);
}
void SparseVector::push_back(std::pair<int, double> p, double tollerance) {
    if (tollerance == 0.0 && p.second == 0.0) {
        return;
    } else if (!std::signbit(tollerance) && fabs(p.second) < tollerance) {
        return;
    }

    data.push_back(p);
}

int SparseVector::getSize() { return this->size; }

SparseVector SparseVector::operator+(const SparseVector& x) {
    assert(size == x.size);

    SparseVector ans = SparseVector(size);

    std::vector<std::pair<int, double>>::const_iterator it, itx;
    it = data.begin();
    itx = x.data.begin();

    while (it != data.end() && itx != x.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (it->first == itx->first) {
            ans.push_back({it->first, it->second + itx->second});

            ++it;
            ++itx;
        } else if (it->first < itx->first) {
            ans.push_back({it->first, it->second});

            ++it;
        } else {
            ans.push_back({itx->first, itx->second});

            ++itx;
        }
    }
    while (it != data.end()) {
        ans.push_back({it->first, it->second});

        ++it;
    }
    while (itx != x.data.end()) {
        ans.push_back({itx->first, itx->second});

        ++itx;
    }

    return ans;
}

SparseVector SparseVector::operator+=(const SparseVector& x) {
    assert(size == x.size);

    std::vector<std::pair<int, double>> buffer;

    std::vector<std::pair<int, double>>::iterator it;
    std::vector<std::pair<int, double>>::const_iterator itx;
    it = data.begin();
    itx = x.data.begin();

    while (it != data.end() && itx != x.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (it->first == itx->first) {
            it->second += itx->second;

            ++it;
            ++itx;
        } else if (it->first < itx->first) {
            ++it;
        } else {
            if (itx->second != 0.0) buffer.push_back({itx->first, itx->second});

            ++itx;
        }
    }
    while (itx != x.data.end()) {
        if (itx->second != 0.0) buffer.push_back({itx->first, itx->second});

        ++itx;
    }

    // add the missing value to the actual sparse vector
    // (used a buffer to not invalidate the iterator)
    data.insert(data.end(), buffer.begin(), buffer.end());
    std::sort(data.begin(), data.end());  // TODO(lugot): PERFORMANCE

    return *this;
}

SparseVector SparseVector::operator-(const SparseVector& x) {
    assert(size == x.size);

    SparseVector ans = SparseVector(size);

    std::vector<std::pair<int, double>>::const_iterator it, itx;
    it = data.begin();
    itx = x.data.begin();

    while (it != data.end() && itx != x.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (it->first == itx->first) {
            ans.push_back({it->first, it->second - itx->second});

            ++it;
            ++itx;
        } else if (it->first < itx->first) {
            ans.push_back({it->first, it->second});

            ++it;
        } else {
            ans.push_back({itx->first, -itx->second});

            ++itx;
        }
    }
    while (it != data.end()) {
        ans.push_back({it->first, it->second});

        ++it;
    }
    while (itx != x.data.end()) {
        ans.push_back({itx->first, -itx->second});

        ++itx;
    }

    return ans;
}

SparseVector SparseVector::operator-=(const SparseVector& x) {
    assert(size == x.size);

    std::vector<std::pair<int, double>> buffer;

    std::vector<std::pair<int, double>>::iterator it;
    std::vector<std::pair<int, double>>::const_iterator itx;
    it = data.begin();
    itx = x.data.begin();

    while (it != data.end() && itx != x.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (it->first == itx->first) {
            it->second -= itx->second;

            ++it;
            ++itx;
        } else if (it->first < itx->first) {
            ++it;
        } else {
            if (itx->second != 0.0) {
                buffer.push_back({itx->first, -itx->second});
            }

            ++itx;
        }
    }
    while (itx != x.data.end()) {
        if (itx->second != 0.0) buffer.push_back({itx->first, -itx->second});

        ++itx;
    }

    // add the missing value to the actual sparse vector
    // (used a buffer to not invalidate the iterator)
    data.insert(data.end(), buffer.begin(), buffer.end());
    std::sort(data.begin(), data.end());  // TODO(lugot): PERFORMANCE

    return *this;
}

double SparseVector::operator*(const SparseVector& x) {
    assert(size == x.size);

    double ans = 0.0f;

    std::vector<std::pair<int, double>>::const_iterator it, itx;
    it = data.begin();
    itx = x.data.begin();

    while (it != data.end() && itx != x.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (it->first == itx->first) {
            ans += it->second * itx->second;

            ++it;
            ++itx;
        } else if (it->first < itx->first) {
            ++it;
        } else {
            ++itx;
        }
    }

    return ans;
}

SparseVector SparseVector::operator*(const SparseMatrix& A) {
    assert(size == A.nrows);

    SparseVector ans(A.ncols);
    for (int i = 0; i < A.ncols; ++i) ans.push_back({i, (*this) * A.cols[i]});

    return ans;
}

SparseVector SparseVector::operator*(double scalar) {
    SparseVector ans(size);

    std::vector<std::pair<int, double>>::iterator it;
    for (it = data.begin(); it != data.end(); ++it) {
        ans.push_back({it->first, scalar * it->second});
    }

    return ans;
}

SparseVector SparseVector::operator*=(double scalar) {
    std::vector<std::pair<int, double>>::iterator it;
    for (it = data.begin(); it != data.end(); ++it) {
        it->second *= scalar;
    }
    return *this;
}

bool SparseVector::operator<(const SparseVector& x) {
    assert(size == x.size);

    std::vector<std::pair<int, double>>::const_iterator it, itx;
    it = data.begin();
    itx = x.data.begin();

    while (it != data.end() && itx != x.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (it->first == itx->first) {
            if (itx->second - it->second < EPS) return false;

            ++it;
            ++itx;
        } else if (it->first < itx->first) {
            if (-it->second < EPS) return false;

            ++it;
        } else {
            if (itx->second < EPS) return false;

            ++itx;
        }
    }
    while (it != data.end()) {
        if (-it->second < EPS) return false;

        ++it;
    }
    while (itx != x.data.end()) {
        if (itx->second < EPS) return false;

        ++itx;
    }

    return true;
}
bool SparseVector::operator==(const SparseVector& x) {
    assert(size == x.size);

    std::vector<std::pair<int, double>>::const_iterator it, itx;
    it = data.begin();
    itx = x.data.begin();

    while (it != data.end() && itx != x.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (it->first == itx->first) {
            if (fabs(it->second - itx->second) > EPS) return false;

            ++it;
            ++itx;
        } else if (it->first < itx->first) {
            if (fabs(it->second) > EPS) return false;

            ++it;
        } else {
            if (fabs(itx->second) > EPS) return false;

            ++itx;
        }
    }
    while (it != data.end()) {
        if (fabs(it->second) > EPS) return false;

        ++it;
    }
    while (itx != x.data.end()) {
        if (fabs(itx->second) > EPS) return false;

        ++itx;
    }

    return true;
}

double SparseVector::squaredNorm(const SparseVector& x) {
    double ans = 0.0f;

    std::vector<std::pair<int, double>>::const_iterator it;
    for (it = x.data.begin(); it != x.data.end(); ++it) {
        ans += (it->second) * (it->second);
    }

    return ans;
}

double SparseVector::dist(const SparseVector& x, const SparseVector& y) {
    assert(x.size == y.size);

    double ans = 0.0f;

    std::vector<std::pair<int, double>>::const_iterator itx, ity;
    itx = x.data.begin();
    ity = y.data.begin();

    while (itx != x.data.end() && ity != y.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (itx->first == ity->first) {
            ans += (itx->second - ity->second) * (itx->second - ity->second);

            ++itx;
            ++ity;
        } else if (itx->first < ity->first) {
            ans += (itx->second) * (itx->second);

            ++itx;
        } else {
            ans += (ity->second) * (ity->second);

            ++ity;
        }
    }
    while (itx != x.data.end()) {
        ans += (itx->second) * (itx->second);

        ++itx;
    }
    while (ity != y.data.end()) {
        ans += (ity->second) * (ity->second);

        ++ity;
    }

    return ans;
}

bool SparseVector::isZero(int idx) {
    std::vector<std::pair<int, double>>::iterator it;
    it = std::lower_bound(data.begin(), data.end(), idx,
                          [](const std::pair<int, double>& p, int target) {
                              return p.first < target;
                          });

    if (it == data.end() || it->first != idx || it->second == 0.0) return true;
    return false;
}

void SparseVector::prune(double tollerance) {
    data.erase(std::remove_if(data.begin(), data.end(),
                              [tollerance](const std::pair<int, double>& p) {
                                  return fabs(p.second) < tollerance;
                              }),
               data.end());
}

std::ostream& operator<<(std::ostream& os, const SparseVector& x) {
    os << x.size << "(" << x.data.size() << ") ";
    if (VERBOSE) {
        std::vector<std::pair<int, double>>::const_iterator it;
        for (it = x.data.begin(); it != x.data.end(); ++it) {
            os << "(" << it->first << "," << it->second << ") ";
        }
    }

    return os;
}
