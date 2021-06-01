#include "../include/utils.h"

#include "../include/globals.h"

int sparseVectorCompareZero(const Eigen::SparseVector<double>& l,
                            const Eigen::SparseVector<double>& r) {
    assert(l.rows() == r.rows());

    // merge-like code structure
    // O( nnz(lhs) + nnz(b) ) = O ( 2m ), m = l.rows() = r.rows()
    Eigen::SparseVector<double>::InnerIterator itl(l);
    Eigen::SparseVector<double>::InnerIterator itr(r);
    while (itl && itr) {
        if (itl.index() < itr.index()) {
            // since itr.index is zero
            // TODO(lugot): CHECK floating point safety
            if (-itl.value() < EPSILON) return itl.index();
            ++itl;
        } else if (itl.index() > itr.index()) {
            // since itl.index is zero
            if (itr.value() < EPSILON) return itr.index();
            ++itr;
        } else {  // itl.index() == itb.index()
            // itl.value() > itr.value() -> infeasible
            if (itr.value() - itl.value() < EPSILON) return itl.index();
            ++itl;
            ++itr;
        }
    }
    while (itl) {
        // since itb.index is zero
        // TODO(lugot): CHECK floating point safety
        if (-itl.value() < EPSILON) return itl.index();
        ++itl;
    }
    while (itr) {
        // since itl.index is zero
        if (itr.value() < EPSILON) return itr.index();
        ++itr;
    }

    return l.rows();
}

int sparseVectorCompareInf(const Eigen::SparseVector<double>& l,
                           const Eigen::SparseVector<double>& r) {
    assert(l.rows() == r.rows());

    // merge-like code structure
    // O( nnz(lhs) + nnz(b) ) = O ( 2m ), m = l.rows() = r.rows()
    Eigen::SparseVector<double>::InnerIterator itl(l);
    Eigen::SparseVector<double>::InnerIterator itr(r);
    while (itl && itr) {
        if (itl.index() < itr.index()) {
            // no need to check since itr.value() is inf
            ++itl;
        } else if (itl.index() > itr.index()) {
            // since itl.index is zero
            if (itr.value() < EPSILON) return itr.index();
            ++itr;
        } else {  // itl.index() == itb.index()
            // itl.value() > itr.value() -> infeasible
            if (itr.value() - itl.value() < EPSILON) return itl.index();
            ++itl;
            ++itr;
        }
    }
    // skipped
    // while (itl) {
    //++itl;
    // }
    while (itr) {
        // since itl.index is zero
        if (itr.value() < EPSILON) return itr.index();
        ++itr;
    }

    return l.rows();
}
