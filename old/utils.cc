#include "../include/utils.h"

template <class T>
const std::ostream& operator<<(std::ostream& os, const Eigen::Triplet<T>& t) {
    os << "pos: (" << t.row() << ". " << t.col() << ") value: " << t.value()
       << std::endl;
    return os;
}

std::vector<std::string> split(std::string line) {
    std::stringstream ss(line);
    std::vector<std::string> tokenized_line;
    std::string token;

    while (ss >> token) tokenized_line.push_back(token);

    return tokenized_line;
}
