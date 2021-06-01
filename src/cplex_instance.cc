#include "../include/cplex_instance.h"

#include "../include/globals.h"
#include "ilcplex/ilocplex.h"

CplexInstance::CplexInstance() : Instance() {}

CplexInstance::CplexInstance(const std::string& filename) : Instance() {
    // do nothing
    importModel(filename);

    int i = 0;
    while (filename[i] != '.') i++;
    this->model_name = filename.substr(0, i);
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

CplexInstance::~CplexInstance() {
    cplex.end();
    rng.end();
    var.end();
    obj.end();
    model.end();
    env.end();
}

void CplexInstance::importModel(const std::string& model_name) {
    env = IloEnv();

    try {
        IloModel modeltmp = IloModel(env);
        cplex = IloCplex(env);

        obj = IloObjective(env);
        var = IloNumVarArray(env);
        rng = IloRangeArray(env);

        // file not found managed by exception
        std::string filename = "../data/" + model_name;
        cplex.importModel(modeltmp, filename.c_str(), obj, var, rng);

        // parse the model file
        cplex.extract(modeltmp);

        model = IloModel(env);
        model.add(modeltmp);

        model.add(IloConversion(env, var, IloNumVar::Type::Float));
        cplex.extract(model);

        // TODO(lugot): FIX MEMLEAK

        if (EXTRA) {
            std::cout << "[VERBOSE] Variable bounds" << std::endl;
            for (int i = 0; i < var.getSize(); ++i) {
                std::cout << var[i].getName() << ": " << var[i].getLB() << " "
                          << var[i].getUB() << " " << var[i].getType()
                          << std::endl;
                /* assert(var[i].getLB() != -IloInfinity);
                assert(var[i].getUB() != IloInfinity); */
            }
        }

        // modeltmp.end();
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

        vals.end();
    } catch (IloException& e) {
        std::cerr << "Concert exception caught: " << e << std::endl;
    }

    return run_status;
}
