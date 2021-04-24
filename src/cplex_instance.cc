#include "../include/cplex_instance.h"

#include "../include/globals.h"
#include "ilcplex/ilocplex.h"

CplexInstance::CplexInstance() { status = CplexInstance::Status::NeverSolved; }

CplexInstance::CplexInstance(const std::string& model_name) : CplexInstance() {
    // do nothing
    importModel(model_name);
}

CplexInstance::CplexInstance(const CplexInstance& other) {
    status = other.status;

    if (status != CplexInstance::Status::EmptyInstance) {
        // clone all the CPLEX things, everything can be extracteb my the model
        // I guess but let's hard clone everything
        // still dont know if exists a better way to do this

        env = IloEnv();
        model = IloModel(env);

        var = IloNumVarArray(env);
        for (int i = 0; i < other.var.getSize(); ++i) {
            IloNum lb = other.var[i].getLB();
            IloNum ub = other.var[i].getLB();
            std::string varname = other.var[i].getName();

            var.add(IloNumVar(env, lb, ub, varname.c_str()));
        }

        rng = IloRangeArray(env);
        for (int i = 0; i < other.rng.getSize(); ++i) {
            rng.add(IloRange(env, other.rng[i].getExpr()));
        }

        obj = IloObjective(env, other.obj.getExpr());

        cplex = IloCplex(model);
    }
}

void CplexInstance::importModel(const std::string& model_name) {
    env = IloEnv();

    try {
        model = IloModel(env);
        cplex = IloCplex(env);

        obj = IloObjective(env);
        var = IloNumVarArray(env);
        rng = IloRangeArray(env);

        // file not found managed by exception
        std::string filename = "../data/" + model_name;
        cplex.importModel(model, filename.c_str(), obj, var, rng);

        // parse the model file
        cplex.extract(model);
    } catch (IloException& e) {
        std::cerr << "Concert exception caught: " << e << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception caught" << std::endl;
    }
}

bool CplexInstance::solve() {
    // TODO(lugot): check for status
    bool run_status;

    try {
        run_status = static_cast<bool>(cplex.solve());
    } catch (IloException& e) {
        std::cerr << "Concert exception caught: " << e << std::endl;
    }

    // set internal status (cool conversion due to order)
    status = static_cast<Status>(cplex.getStatus());

    return run_status;
}

CplexInstance::Status CplexInstance::getStatus() { return status; }

void CplexInstance::printStatus() {
    IloNumArray vals(env);

    env.out() << "Solution status = " << cplex.getStatus() << std::endl;
    env.out() << "Solution value  = " << cplex.getObjValue() << std::endl;

    // TODO(lugot): do it better
    cplex.getValues(vals, var);
    env.out() << "Values        = " << vals << std::endl;
    cplex.getSlacks(vals, rng);
    env.out() << "Slacks        = " << vals << std::endl;
    cplex.getDuals(vals, rng);
    env.out() << "Duals         = " << vals << std::endl;
    cplex.getReducedCosts(vals, var);
    env.out() << "Reduced Costs = " << vals << std::endl;
}

bool CplexInstance::isSolved() {
    return status == Status::Feasible || status == Status::Optimal ||
           status == Status::Infeasible || status == Status::Unbounded ||
           status == Status::InfeasibleOrUnbounded;
}

CanonicalForm CplexInstance::getCanonicalForm() {
    int n = cplex.getNcols();
    int m = cplex.getNrows();

    // add one rows for each equality constraint
    for (int i = 0; i < rng.getSize(); ++i) {
        if (rng[i].getLB() == rng[i].getUB()) m++;
    }

    Eigen::SparseMatrix<double>* A = new Eigen::SparseMatrix<double>(m, n);
    Eigen::SparseVector<double>* b = new Eigen::SparseVector<double>(m);
    Eigen::SparseVector<double>* c = new Eigen::SparseVector<double>(n);
    Eigen::SparseVector<double>* lb = new Eigen::SparseVector<double>(n);
    Eigen::SparseVector<double>* ub = new Eigen::SparseVector<double>(n);

    // fill the objective. The names of variables are tracked down
    std::map<std::string, int> var_indexer;  // Variable To Column
    // TODO(lugot): abort if the problem is not LP
    for (IloExpr::LinearIterator it = obj.getLinearIterator(); it.ok(); ++it) {
        std::string varname = it.getVar().getName();
        if (var_indexer[varname] == 0)
            var_indexer[varname] = var_indexer.size();

        c->coeffRef(var_indexer[varname]) = it.getCoef();
    }

    // fill the constraints using triplets for speedup TODO(lugot) VERIFY
    std::vector<Eigen::Triplet<double>> triplets;
    int row;
    for (int i = 0; i < rng.getSize(); i++) {
        IloRange& r = rng[i];

        // get the lhs coefficents
        std::vector<Eigen::Triplet<double>> constraint;
        for (IloExpr::LinearIterator it = r.getLinearIterator(); it.ok();
             ++it) {
            std::string varname = it.getVar().getName();
            if (var_indexer[varname] == 0)
                var_indexer[varname] = var_indexer.size();

            constraint.push_back(Eigen::Triplet<double>(
                row, var_indexer[varname], it.getCoef()));
        }

        // add the constaint to the general triplets
        triplets.insert(triplets.end(), constraint.begin(), constraint.end());
        if (r.getUB() != 0.0) b->coeffRef(row) = r.getUB();
        row++;

        // check for equality constraint
        if (r.getLB() == r.getUB()) {
            // swap the coefficent sign and
            for (Eigen::Triplet<double>& t : constraint) {
                t = Eigen::Triplet<double>(t.row() + 1, t.col(), -t.value());
            }

            // add the new row
            triplets.insert(triplets.end(), constraint.begin(),
                            constraint.end());
            if (r.getUB() != 0.0) b->coeffRef(row) = r.getUB();
            row++;
        }
    }

    // fill the LBs and UBs
    for (int i = 0; i < var.getSize(); ++i) {
        if (var[i].getLB() != 0.0) lb->coeffRef(i) = var[i].getLB();
        if (var[i].getUB() != 0.0) ub->coeffRef(i) = var[i].getUB();
    }

    return {A, b, c, lb, ub};
}

std::ostream& operator<<(std::ostream& os, const CplexInstance::Status& s) {
    switch (s) {
        case CplexInstance::Status::EmptyInstance:
            os << "EmptyInstance";
            break;
        case CplexInstance::Status::NeverSolved:
            os << "NeverSolved";
            break;
        case CplexInstance::Status::ModelChanged:
            os << "ModelChanged";
            break;
        default:
            // get back to IloAlgorithm::Status << overloaded operator
            os << static_cast<IloAlgorithm::Status>(s);
            break;
    }

    return os;
}
