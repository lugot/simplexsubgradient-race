#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include <Eigen/Sparse>

/**
 * Perform l < r pointwise. Consider nonzeros as zeroes.
 *
 * @param l left sparse vector
 * @param r right sparse vector
 * @return first variable that not satisfy l < r, vector lenght otherwise.
 */
int sparseVectorCompareZero(const Eigen::SparseVector<double>& l,
                            const Eigen::SparseVector<double>& r);

/**
 * Perform l < r pointwise. Consider nonzeros of r  as IloInfinite (or double
 * inf).
 *
 * @param l left sparse vector
 * @param r right sparse vector
 * @return first variable that not satisfy l < r, vector lenght otherwise.
 */
int sparseVectorCompareInf(const Eigen::SparseVector<double>& l,
                           const Eigen::SparseVector<double>& r);

#endif /* INCLUDE_UTILS_H_ */
