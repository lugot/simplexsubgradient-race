#include "../include/sparse_vector.h"

#include <cassert>
#include <cmath>

#include "../include/globals.h"

SparseVector::SparseVector(int size) { this->size = size; }

void SparseVector::push_back(std::pair<int, double> p) {
    data.push_back(p);
}

int SparseVector::getSize() { return this->size; }

bool SparseVector::negative() {
    std::vector<std::pair<int, double>>::iterator it;
    for (it = data.begin(); it != data.end(); ++it) {
        if (!std::signbit(it->second)) return false;
    }

    return true;
}

void accumulate(double scalar);

double SparseVector::squaredNorm() {
    double ans = 0.0f;

    std::vector<std::pair<int, double>>::iterator it;
    for (it = data.begin(); it != data.end(); ++it) {
        ans += (it->second) * (it->second);
    }

    return ans;
}

SparseVector SparseVector::sum(const SparseVector& x, const SparseVector& y) {
    assert(x.size == y.size);

    SparseVector ans = SparseVector(x.size);

    std::vector<std::pair<int, double>>::const_iterator itx, ity;
    itx = x.data.begin();
    ity = y.data.begin();

    while (itx != x.data.end() && ity != y.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (itx->first == ity->first) {
            ans.data.push_back({itx->first, itx->second + ity->second});

            ++itx;
            ++ity;
        } else if (itx->first < ity->first) {
            ans.data.push_back({itx->first, itx->second});

            ++itx;
        } else {
            ans.data.push_back({ity->first, ity->second});

            ++ity;
        }
    }
    while (itx != x.data.end()) {
        ans.data.push_back({itx->first, itx->second});

        ++itx;
    }
    while (ity != y.data.end()) {
        ans.data.push_back({ity->first, ity->second});

        ++ity;
    }

    return ans;
}

SparseVector SparseVector::dot(const SparseVector& x, const SparseVector& y) {
    assert(x.size == y.size);

    SparseVector ans = SparseVector(x.size);

    std::vector<std::pair<int, double>>::const_iterator itx, ity;
    itx = x.data.begin();
    ity = y.data.begin();

    while (itx != x.data.end() && ity != y.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (itx->first == ity->first) {
            ans.data.push_back({itx->first, itx->second * ity->second});

            ++itx;
            ++ity;
        } else if (itx->first < ity->first) {
            ++itx;
        } else {
            ++ity;
        }
    }

    return ans;
}

bool SparseVector::negativeDot(const SparseVector& x, const SparseVector& y) {
    assert(x.size == y.size);

    std::vector<std::pair<int, double>>::const_iterator itx, ity;
    itx = x.data.begin();
    ity = y.data.begin();

    while (itx != x.data.end() && ity != y.data.end()) {
        // TODO(lugot): PERFORMANCE
        if (itx->first == ity->first) {
            if (std::signbit(itx->second) == std::signbit(ity->second)) {
                // signbit returns if is negative
                // ++ and -- are positive
                return false;
            }

            ++itx;
            ++ity;
        } else if (itx->first < ity->first) {
            ++itx;
        } else {
            ++ity;
        }
    }

    return true;
}

std::ostream& operator<<(std::ostream& os, const SparseVector& x) {
    os << "dim " << x.size << ": ";
    if (VERBOSE) {
        std::vector<std::pair<int, double>>::const_iterator it;
        for (it = x.data.begin(); it != x.data.end(); ++it) {
            os << "(" << it->first << ", " << it->second << ") ";
        }
        os << "\n";
    }

    return os;
}
