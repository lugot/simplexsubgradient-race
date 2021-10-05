#include "../include/problem.h"

#include "../include/globals.h"

Problem::Problem(const Cplex& cinst) {
    // model sizes
    n = cinst.var.getSize();
    m = cinst.rng.getSize();

    // vectors & friends
    c = SparseVector(n);
    b = SparseVector(m);
    A = SparseMatrix(m, n);
    lb = std::vector<double>(n, LOWERBOUND);
    ub = std::vector<double>(n, UPPERBOUND);
    equal = std::vector<bool>(m);

    // model name for file retreiving
    instance_name = cinst.model_name;

    // track the names of the variables for consistency
    varindexer = std::map<std::string, int>();

    // fill the objective sparse vector
    IloExpr::LinearIterator ito;
    for (ito = cinst.obj.getLinearIterator(); ito.ok(); ++ito) {
        std::string varname = ito.getVar().getName();
        if (varindexer.find(varname) == varindexer.end()) {
            varindexer[varname] = varindexer.size();
        }

        c.push_back({varindexer[varname], ito.getCoef()});
    }

    // fill A and b (formulation Ax >= b)
    std::vector<std::tuple<int, int, double>> lhs;
    double rhs;
    int row = 0;
    for (int i = 0; i < cinst.rng.getSize(); ++i) {
        IloRange r = cinst.rng[i];

        // extract lhs and rhs from IloRange object ..
        extractConstraint(r, row, &lhs, &rhs);
        // .. and add them to the sparse vector and matrix
        A.addRow(lhs);
        b.push_back({row, rhs});

        if (constraintSense(r) == EQUAL) {
            equal[row] = true;
        }

        row++;
    }

    // change the variable bounds, if present
    for (int i = 0; i < cinst.var.getSize(); ++i) {
        lb[i] = std::max(lb[i], cinst.var[i].getLB());
        ub[i] = std::min(ub[i], cinst.var[i].getUB());
    }
}

void Problem::extractConstraint(const IloRange& r, int row,
                                std::vector<std::tuple<int, int, double>>* a,
                                double* b) {
    // get number of non zero coefficents
    int nnz = 0;
    IloExpr::LinearIterator it;
    for (it = r.getLinearIterator(); it.ok(); ++it) nnz++;
    a->resize(nnz);  // TODO(lugot): reserve better

    // understand the sense of constraint, set a mutiplicative factor, we
    // need Ax >= b format
    Sense s = constraintSense(r);
    int mult = s == LESS ? -1 : 1;

    int i = 0;
    for (it = r.getLinearIterator(); it.ok(); ++it) {
        std::string varname = it.getVar().getName();
        if (varindexer.find(varname) == varindexer.end()) {
            varindexer[varname] = varindexer.size();
        }

        a->at(i) = {row, varindexer[varname], mult * it.getCoef()};

        i++;
    }

    if (s == LESS) {
        *b = -r.getUB();
    } else {
        *b = r.getLB();
    }
}

Problem::Sense Problem::constraintSense(const IloRange& r) {
    // aT x = b -> ub == lb
    if (r.getUB() == r.getLB()) return EQUAL;
    // aT x >= b -> ub == inf
    if (r.getUB() == IloInfinity) return GREATER;
    // aT x <= b -> lb == inf
    return LESS;
}
