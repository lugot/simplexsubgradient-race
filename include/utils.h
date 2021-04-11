#ifndef _UTILS_H_
#define _UTILS_H_

#include <Eigen/Sparse>

template<class T>
const std::ostream& operator<<(std::ostream& os, const Eigen::Triplet<T>& t);
std::vector<std::string> split(std::string line);

#endif /* _UTILS_H_ */
