#include "../include/cplex_instance.h"

#include "../include/globals.h"
#include "ilcplex/ilocplex.h"

CplexInstance::CplexInstance() : Instance() {}

CplexInstance::CplexInstance(const std::string& model_name) : Instance() {
    // do nothing
    importModel(model_name);
}

CplexInstance::CplexInstance(const CplexInstance& other) {
    status = other.status;

    if (status != SolutionStatus::EmptyInstance) {
        // clone all the CPLEX things, everything can be extracteb my the model
        // I guess but let's hard clone everything
        // still dont know if exists a better way to do this

        env = IloEnv();
        model = IloModel(env);

        var = IloNumVarArray(env);
        for (int i = 0; i < other.var.getSize(); ++i) {
            IloNum lb = other.var[i].getLB();
            IloNum ub = other.var[i].getUB();
            std::string varname = other.var[i].getName();

            var.add(IloNumVar(env, lb, ub, varname.c_str()));
        }

        rng = IloRangeArray(env);
        for (int i = 0; i < other.rng.getSize(); ++i) {
            IloNum lb = other.rng[i].getLB();
            IloNum ub = other.rng[i].getUB();
            IloExpr expr = other.rng[i].getExpr();
            std::string conname = other.var[i].getName();

            rng.add(IloRange(env, lb, expr, ub, conname.c_str()));
        }

        obj = IloObjective(env, other.obj.getExpr(), other.obj.getSense(),
                           other.obj.getName());

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
    // TODO(lugot): CHECK for status
    bool run_status;

    try {
        std::cout << obj.getExpr() << std::endl;

        run_status = static_cast<bool>(cplex.solve());

        // set internal status (cool conversion due to order)
        status = static_cast<SolutionStatus>(cplex.getStatus());

        // TODO(lugot): EXTRACT in a function
        IloNumArray vals(env);
        env.out() << "Solution status = " << cplex.getStatus() << std::endl;
        env.out() << "Solution value  = " << cplex.getObjValue() << std::endl;
        cplex.getValues(vals, var);
        env.out() << "Values        = " << vals << std::endl;
        cplex.getSlacks(vals, rng);
        env.out() << "Slacks        = " << vals << std::endl;
        cplex.getDuals(vals, rng);
        env.out() << "Duals         = " << vals << std::endl;
        cplex.getReducedCosts(vals, var);
        env.out() << "Reduced Costs = " << vals << std::endl;
    } catch (IloException& e) {
        std::cerr << "Concert exception caught: " << e << std::endl;
    }

    return run_status;
}

void CplexInstance::updateObjective(const Eigen::SparseVector<double>& x) {
    assert(x.rows() == var.getSize());

    //model.remove(obj);
    obj.end();


    IloExpr qexpr = IloExpr(env);

    // TODO(lugot): PERFORMANCE check other qp entering
    for (int i = 0; i < var.getSize(); ++i) {
        qexpr += 0.5 * var[i] * var[i];
    }
    for (Eigen::SparseVector<double>::InnerIterator it(x); it; ++it) {
        qexpr += -it.value() * var[it.index()];
        qexpr += 0.5 * it.value() * it.value();
    }

    obj = IloObjective(env, qexpr);
    // obj.setExpr(qexpr);
    model.add(obj);
    model.add(var);
    model.add(rng);
    std::cout << "model: " << model << std::endl;
}

CanonicalForm CplexInstance::getCanonicalForm() {
    // TODO(lugot): UNDERSTAND why cplex obj here has n=0, m=0
    int n = var.getSize();
    int m = rng.getSize();

    // add one rows for each equality constraint
    for (int i = 0; i < rng.getSize(); ++i) {
        if (rng[i].getLB() == rng[i].getUB()) m++;
    }

    Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>(m, n);
    Eigen::SparseVector<double> b = Eigen::SparseVector<double>(m);
    Eigen::SparseVector<double> c = Eigen::SparseVector<double>(n);
    Eigen::SparseVector<double> lb = Eigen::SparseVector<double>(n);
    Eigen::SparseVector<double> ub = Eigen::SparseVector<double>(n);

    // fill the objective. The names of variables are tracked down
    std::map<std::string, int> var_indexer;  // Variable To Column
    // TODO(lugot): abort if the problem is not LP
    for (IloExpr::LinearIterator it = obj.getLinearIterator(); it.ok(); ++it) {
        std::string varname = it.getVar().getName();
        if (var_indexer.find(varname) == var_indexer.end())
            var_indexer[varname] = var_indexer.size();

        c.coeffRef(var_indexer[varname]) = it.getCoef();
    }

    // fill the constraints using triplets for speedup TODO(lugot) VERIFY
    std::vector<Eigen::Triplet<double>> triplets;
    int row = 0;
    for (int i = 0; i < rng.getSize(); i++) {
        IloRange& r = rng[i];

        // get the lhs coefficents
        std::vector<Eigen::Triplet<double>> constraint;
        for (IloExpr::LinearIterator it = r.getLinearIterator(); it.ok();
             ++it) {
            std::string varname = it.getVar().getName();
            if (var_indexer.find(varname) == var_indexer.end())
                var_indexer[varname] = var_indexer.size();

            constraint.push_back(Eigen::Triplet<double>(
                row, var_indexer[varname], it.getCoef()));
        }

        // add the constaint to the general triplets
        triplets.insert(triplets.end(), constraint.begin(), constraint.end());
        if (r.getUB() != 0.0) b.coeffRef(row) = r.getUB();
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
            if (r.getUB() != 0.0) b.coeffRef(row) = r.getUB();
            row++;
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());

    // fill the LBs and UBs
    for (int i = 0; i < var.getSize(); ++i) {
        // lb: nonzero element are actyally zero
        if (var[i].getLB() != 0.0) lb.coeffRef(i) = var[i].getLB();
        // ub: nozero element are inf: skip those checks
        if (var[i].getUB() != IloInfinity) ub.coeffRef(i) = var[i].getUB();
    }

    // TODO(lugot): UNDERSTAND copy elision and return value optimization
    return {A, b, c, lb, ub};
}
