#include "../include/cplex.h"

#include <chrono>
#include <iostream>
#include <ostream>

#include "../include/globals.h"
#include "ilcplex/ilocplex.h"

Cplex::Cplex(const Cplex& other) {
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

Cplex::~Cplex() {
    log.close();

    cplex.end();
    rng.end();
    var.end();
    obj.end();
    model.end();
    env.end();
}

void Cplex::importModel(const std::string& model_name) {
    // also fix model name
    int i = 0;
    while (model_name[i] != '.') i++;
    this->model_name = model_name.substr(0, i);

    env = IloEnv();

    std::string path = "../results/" + this->model_name + ".cplex_log";
    // std::cout << "path: " << path << std::endl;
    log.open(path.c_str());

    try {
        IloModel modeltmp = IloModel(env);
        cplex = IloCplex(env);
        cplex.setOut(log);

        obj = IloObjective(env);
        var = IloNumVarArray(env);
        rng = IloRangeArray(env);

        // file not found managed by exception
        std::string filename = "../testbed/" + model_name;
        // std::cout << "filename: " << filename << std::endl;
        cplex.importModel(modeltmp, filename.c_str(), obj, var, rng);

        // parse the model file
        cplex.extract(modeltmp);

        model = IloModel(env);
        model.add(modeltmp);

        model.add(IloConversion(env, var, IloNumVar::Type::Float));
        cplex.extract(model);

        // TODO(lugot): FIX MEMLEAK

        if (EXTRA) {
            // std::cout << "[VERBOSE] Variable bounds" << std::endl;
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
        exit(EXIT_FAILURE);
    } catch (...) {
        std::cerr << "Unknown exception caught" << std::endl;
    }

    cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
    cplex.setParam(IloCplex::Param::Threads, 1);
    // cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
    cplex.setParam(IloCplex::Param::Output::WriteLevel, 1);
    cplex.setParam(IloCplex::Param::Simplex::Display, 2);
    // cplex.setParam(IloCplex::Param::Sifting::Simplex, 0);
    // cplex.setParam(IloCplex::Param::Sifting::Display, 2);
    // cplex.setParam(IloCplex::Param::Sifting::Iterations, 0);
    // cplex.setParam(IloCplex::Param::Read::Scale, -1);
}

bool Cplex::solve() {
    bool run_status;
    std::chrono::steady_clock::time_point begin, end;
    begin = std::chrono::steady_clock::now();

    cplex.setOut(log);
    cplex.setWarning(log);
    cplex.setError(log);

    try {
        // std::cout << obj.getExpr() << std::endl;

        run_status = static_cast<bool>(cplex.solve());

        // set internal status (cool conversion due to order)
        status = static_cast<SolutionStatus>(cplex.getStatus());

        // TODO(lugot): EXTRACT in a function
        IloNumArray vals(env);
        // env.out() << "Solution status = " << cplex.getStatus() << std::endl;
        // env.out() << "Solution value  = " << cplex.getObjValue() <<
        // std::endl;

        vals.end();
    } catch (IloException& e) {
        std::cerr << "Concert exception caught: " << e << std::endl;
    }

    end = std::chrono::steady_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
            .count();
    solvetime = elapsed;

    log << "solvetime: " << elapsed;
    if (cplex.getStatus() != IloAlgorithm::Status::Optimal)
        log << " notopt ";
    else
        log << " opt ";
    log << cplex.getObjValue() << std::endl;

    log << "dimension(n-m): " << var.getSize() << " " << rng.getSize();

    return run_status;
}
