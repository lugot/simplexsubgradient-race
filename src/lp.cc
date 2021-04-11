#include "lp.h"

using namespace std;
using namespace Eigen;

//lp::lp(string model_name,
        //SparseVector<double> c, SparseMatrix<double> A, SparseVector<double> b,
        //map<string, int> name_hasher, map<int, string> var_hasher) {

    //this->model_name_ = model_name;

    //this->c_ = c;
    //this->A_ = A;
    //this->b_ = b;

    //this->name_hasher_ = name_hasher;
    //this->var_hasher_ = var_hasher;
//}
//lp::lp(string model_name) {

    //this->model_name_ = model_name;
//}

lp::lp() {
}

string& lp::getModelName() { return this->model_name_; }

SparseVector<double>& lp::getc() { return this->c_; }
SparseMatrix<double>& lp::getA() { return this->A_; }
SparseVector<double>& lp::getb() { return this->b_; }
