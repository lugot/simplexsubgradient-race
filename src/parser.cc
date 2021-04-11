#include "lp.h"
#include "globals.h"
#include <string>
#include <Eigen/Sparse>
#include <ilcplex/ilocplex.h>

using namespace std;
using namespace Eigen;

void lp::importModel(string model_name) {

    int buf[10];
    for (int i=0; i<100; i++) buf[i] = 0;

    env_ = IloEnv();

    try {
        model_ = IloModel(env_);
        cplex_ = IloCplex(env_);
        cplex_.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);

        obj_ = IloObjective(env_);
        var_ = IloNumVarArray(env_);
        rng_ = IloRangeArray(env_);
        string filename = "../data/" + model_name;
        cplex_.importModel(model_, filename.c_str(), obj_, var_, rng_);

        // parse the model file
        cplex_.extract(model_);

        // retreive the size of the model
        int n = cplex_.getNcols();
        int m = cplex_.getNrows();

        // equality to inequality transformation
        for (int i=0; i<rng_.getSize(); ++i) if (rng_[i].getLB() == rng_[i].getUB()) m++;
        // variable domain to constraints transformation
        for (int i=0; i<var_.getSize(); i++) if (var_[i].getUB() != numeric_limits<double>::infinity()) m++;
        m += n;

        // populate objective
        c_.resize(n);
        for (IloExpr::LinearIterator it = obj_.getLinearIterator(); it.ok(); ++it) {
            c_.coeffRef(it.getVar().getId()-1) = it.getCoef();
        }

        // populate constraints
        A_.resize(m, n);
        b_.resize(m);

        cout << "n: " << n << ", m: " << m << endl;

        //alternative if too slow:
        //1: SparseMatrix<double> mat(rows,cols);         // default is column major
        //2: mat.reserve(VectorXi::Constant(cols,6));
        //3: for each i,j such that v_ij != 0
        //4:   mat.insert(i,j) = v_ij;                    // alternative: mat.coeffRef(i,j) += v_ij;
        //5: mat.makeCompressed();                        // optional

        vector<Triplet<double>> triplets;
        int row_idx = 0;

        // iterate over constraints
        vector<Triplet<double>> triplets_row;
        for (int i=0; i<rng_.getSize(); ++i) {
            IloRange& r = rng_[i]; //TODO: cant belive

            if (r.getLB() != -IloInfinity) {
                //throw IloException("Cannot handle LB on constraints");
            }

            for (IloExpr::LinearIterator it = r.getLinearIterator(); it.ok(); ++it) {
                Triplet<double> t = Triplet<double>(row_idx, it.getVar().getId()-1, it.getCoef());
                triplets_row.push_back(t);
            }

            triplets.insert(triplets.end(), triplets_row.begin(), triplets_row.end());
            if (r.getUB() != 0.0) b_.coeffRef(row_idx) = r.getUB();

            if (rng_[i].getLB() == rng_[i].getUB()) {

                // duplicate constraint with swapped signs
                for (Triplet<double>& t: triplets_row) t = Triplet<double>(t.row()+1, t.col(), -t.value());

                triplets.insert(triplets.end(), triplets_row.begin(), triplets_row.end());
                if (r.getUB() != 0.0) b_.coeffRef(row_idx+1) = -r.getUB();

                row_idx++;
            }
            row_idx++;
        }

        cout << triplets.size() << endl;

        int ma, mb;
        ma = mb = 0;

        // iterate over variables
        for (int i=0; i<var_.getSize(); i++) {
            IloNumVar& v = var_[i];

            if (v.getUB() != numeric_limits<double>::infinity()) {
                triplets.push_back(Triplet<double>(row_idx, v.getId()-1, 1));
                if (v.getUB() != 0.0) b_.coeffRef(row_idx) = v.getUB();
                row_idx++;
            }

            triplets.push_back(Triplet<double>(row_idx, v.getId()-1, -1));
            if (v.getLB() != 0.0) b_.coeffRef(row_idx) = -v.getLB();
            row_idx++;

            ma = max(ma, row_idx);
            mb = max(mb, (int)v.getId()-1);
        }

        cout << row_idx << " " << ma << " " << mb << endl;

        // fill A
        A_.setFromTriplets(triplets.begin(), triplets.end());
   }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }

    //env_.end();
}

//pair<string, double> next_variable(vector<string>::iterator& it, vector<string>::iterator enditer);
//bool rhs(string varname);

// how much stupid can i be?
//lp parse_lpmanual(string model_name) {

    //// build model
    //lp model = lp(model_name);

    //SparseVector<double>& c = model.getc();
    //SparseMatrix<double>& A = model.getA();
    //SparseVector<double>& b = model.getb();

    //map<string, int>& name_hasher = model.get_name_hasher();
    //map<int, string>& var_hasher = model.get_var_hasher();


    //// number of variables and constraints
    //int n, m;
    //n = m = 0;


    //// first passage: variables name parser,
    //// faster and cleaner wrt single passage build model
    //ifstream lpfile;
    //lpfile.open("../data/" +  model_name, ios::in);

    //string line;
    //while (getline(lpfile, line)) {
        //m++;

        //vector<string> tokenized_line = split(line);

        //vector<string>::iterator it = tokenized_line.begin()+1;  // avoid header
        //while (it != tokenized_line.end()) {
            //string varname;
            //int coeff;

            //tie(varname, coeff) = next_variable(it, tokenized_line.end());

            //if (varname == "=") m++;
            //if (rhs(varname)) continue;
            //if (!name_hasher.count(varname)) {
                //name_hasher[varname] = n;
                //var_hasher[n++] = varname;
            //}
        //}
    //}
    //lpfile.close();


    ////second passage: store in sparse matrix
    //lpfile.open("../data/" +  model_name, ios::in);

    //// resize and prepare in oneshot the spare matrix and vectors
    //c.resize(n);
    //A.resize(m, m);
    //b.resize(m);

    //// sore the triplets for fast sparse matrix filling
    //vector<Triplet<double>> triplets_row;
    //vector<Triplet<double>> triplets;
    //triplets.reserve(m*n/4);  // totally guess on num of nonzero elements

    //int row = 0;
    //while (getline(lpfile, line)) {

        //vector<string> tokenized_line = split(line);

        //triplets_row.clear();
        //triplets_row.reserve(m);

        //vector<string>::iterator it = tokenized_line.begin()+1;  // avoid header
        //while (it != tokenized_line.end()) {
            //string varname;
            //double coeff;

            //tie(varname, coeff) = next_variable(it, tokenized_line.end());

            //if (!rhs(varname)) {
                //Triplet<double> t = Triplet<double>(row, name_hasher[varname], coeff);
                //triplets_row.push_back(t);
            //}

            //if (rhs(varname)) {
                //// parsing constraint: fill b and copy triplets row in general one
                //triplets.insert(triplets.end(), triplets_row.begin(), triplets_row.end());
                //if (coeff != 0.0) b.coeffRef(row) = coeff;

                //if (varname == "=") {
                    //// duplicate constraint with swapped signs
                    //for (Triplet<double>& t: triplets_row) t = Triplet<double>(t.row()+1, t.col(), -t.value());

                    //triplets.insert(triplets.end(), triplets_row.begin(), triplets_row.end());
                    //if (coeff != 0.0) b.coeffRef(row+1) = -coeff;

                    //row++;
                //}
                //row++;
            //}
            //else if (it == tokenized_line.end()) {
                //// parsing objective: fill the c cost vector with triplets
                //for (Triplet<double> t: triplets_row) c.coeffRef(t.col()) = t.value();
            //}
        //}
    //}
    //// finally fill A
    //A.setFromTriplets(triplets.begin(), triplets.end());
    //lpfile.close();

    //return model;
//}


//pair<string, double> next_variable(vector<string>::iterator& it, vector<string>::iterator enditer) {
    //if (it == enditer) return {"", 0.0};

    //string varname;
    //double coeff;

   //string token = *it;

    //if (next(it) == enditer) {
        //varname = token.substr(1);
        //coeff = token[0] == '+' ? +1 : -1;

        //// end of objective function line
        //if (varname[varname.length()-1] == ';') varname.pop_back();

    //} else {
        //string next_token = *next(it);

        //// end of constraint line
        //if (next_token[next_token.length()-1] ==  ';') next_token.pop_back();

        //if (token == "<" or token == "=") {
            //// parsing <eq/diseq> <rhs>[;]
            //varname = token;
            //coeff = stod(next_token);
            //it++;
        //}
        //else if (next_token[0] == '+' or next_token[0] == '-' or
                //next_token[0] == '=' or next_token[0] == '<') {
            //// parsing <sign><varname> ..
            //varname = token.substr(1);
            //coeff = token[0] == '+' ? +1 : -1;
        //}
        //else {
            //// parsing <coeff> <varname> ..
            //varname = next_token;
            //coeff = stod(token);
            //it++;
        //}
    //}

    //it++;
    //if (VERBOSE) cout << "[DEBUG] parsed " << varname << " " << coeff << endl;

    //return {varname, coeff};
//}


//bool rhs(string varname) {
    //return (varname == "<" or varname == "=");
//}
