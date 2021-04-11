#include "utils.h"

using namespace std;
using namespace Eigen;

template<class T>
const ostream& operator<<(ostream& os, const Triplet<T>& t) {
    os << "pos: (" << t.row() << ". " << t.col() << ") value: " << t.value() << endl;
    return os;
}

vector<string> split(string line) {

      stringstream ss(line);
      vector<string> tokenized_line;
      string token;

      while (ss >> token) tokenized_line.push_back(token);

      return tokenized_line;
}
