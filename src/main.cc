#include <stdlib.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "../include/cplex.h"
#include "../include/globals.h"
#include "../include/subgradient.h"

void print_usage() {
    printf("Usage: ./<name_executable> [options]\n");
    printf("  -v verbose\n");
    printf("  -e extraverbose\n");
    printf("  -l activate logging\n");
    printf("  -d <delim> logging delimitator\n");
    printf("  -i <number> print more info every <number> iterations\n");
    printf("  -n <instance_name>\n");
    printf("  -h --help\n");

    exit(EXIT_SUCCESS);
}

int main(int argc, char* argv[]) {
    std::string single_instance = "mock.lp";

    std::map<std::string, double> mus, phis;

    std::string line;
    std::ifstream optmu("./mus");
    std::cout << optmu.is_open() << std::endl;
    while (getline(optmu, line)) {
        int i = 0;
        while (line[i] != ',') ++i;

        std::string instance_name = line.substr(0, i);

        ++i;
        while (line[i] != ',') ++i;

        double mu = atof(line.substr(i + 1).c_str());

        // std::cout << instance_name << " " << mu << std::endl;
        mus[instance_name + ".mps.gz"] = mu;
    }
    optmu.close();

    int ninst = 0;

    double mu = 0.3;
    double timelimit = -1;
    double solution = 0;

    int opt;
    while ((opt = getopt(argc, argv, "veld:i:n:T:M:N:Iq:w:e:")) != -1) {
        switch (opt) {
            case 'v':
                VERBOSE = true;
                break;
            case 'e':
                VERBOSE = EXTRA = true;
                break;
            case 'l':
                LOGGING = true;
                break;
            case 'I':
                IGNORETIMELIMIT = 1;
                break;
            case 'd':
                DELIM = optarg[0];
                break;
            case 'i':
                ITERINFO = atoi(optarg);
                break;
            case 'n':
                single_instance = std::string(optarg);
                break;
            case 'T':
                TIMILIMIT = atoi(optarg);
                break;
            case 'M':
                MAXITER = atoi(optarg);
                break;
            case 'N':
                ninst = atoi(optarg);
                break;
            case 'q':
                mu = atof(optarg);
                break;
            case 'w':
                timelimit = atof(optarg);
                break;
            case 'h':
            case '?':
                print_usage();
        }
    }

    if (ninst == 0) {
        // std::cout << "single execution on " << single_instance << std::endl;

        std::chrono::steady_clock::time_point begin, end;
        begin = std::chrono::steady_clock::now();

        Cplex cpx(single_instance);
        cpx.cplex.setParam(IloCplex::Param::TimeLimit, 600);
        cpx.solve();
        // std::cout << "cplex optimal solution: " << cpx.cplex.getObjValue()
        //           << std::endl;
        // std::cout << "cplex solvetime: " << cpx.solvetime << std::endl
        //           << std::endl;

        IloNumArray vals(cpx.env);
        cpx.cplex.getValues(vals, cpx.var);

        double ub = 0.0;
        for (int i = 0; i < vals.getSize(); ++i) {
            ub = std::max(ub, vals[i]);
        }
        // std::cout << "variable upper bound: " << ub << std::endl;
        double optimal = cpx.cplex.getObjValue();
        UPPERBOUND = ub;

        end = std::chrono::steady_clock::now();
        double elapsed =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
                .count();
        std::cout << "CPLEX time" << elapsed << " " << cpx.solvetime
                  << std::endl;
        begin = std::chrono::steady_clock::now();

        Problem lp(cpx);
        Subgradient slp;

        double bestmu = -1;
        double bestphi = 0;

        // optimal *= 1.05;
        slp.solve(lp, optimal, cpx.solvetime, mu);

        end = std::chrono::steady_clock::now();
        elapsed =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
                .count();
        std::cout << "subgr time" << elapsed << std::endl;
    } else {
        std::vector<std::string> instances(ninst);

        std::ofstream notopt;
        notopt.open("../notopt");

        std::fstream list;
        list.open("../testbed/list");
        for (std::string& instance_name : instances) list >> instance_name;

        for (std::string instance_name : instances) {
            std::cout << std::endl << std::endl;
            std::cout << "--- multiple execution on " << instance_name << " ---"
                      << std::endl;

            Cplex cpx(instance_name);
            cpx.cplex.setParam(IloCplex::Param::TimeLimit, 3600);
            cpx.solve();
            std::cout << "cplex optimal solution: " << cpx.cplex.getObjValue()
                      << std::endl;
            std::cout << "cplex solvetime: " << cpx.solvetime << std::endl
                      << std::endl;

            bool opt = true;
            if (cpx.cplex.getStatus() != IloAlgorithm::Status::Optimal) {
                std::cout << "NOT_OPTIMAL: " << instance_name << std::endl;
                opt = false;
                notopt << instance_name << std::endl;
            }

            IloNumArray vals(cpx.env);
            cpx.cplex.getValues(vals, cpx.var);

            double ub = 0.0;
            for (int i = 0; i < vals.getSize(); ++i) {
                ub = std::max(ub, vals[i]);
            }
            // if (EXTRA) std::cout << "variable upper bound: " << ub << std::endl;
            // UPPERBOUND = ub;
            double optimal = cpx.cplex.getObjValue();

            Problem lp(cpx);
            cpx.log << " " << lp.nnz() << std::endl;

            Subgradient slp;

            // optimal *= 1.05;
            std::cout << "using mu" << mus[instance_name] << std::endl;
            if (opt) {
                slp.solve(lp, optimal, cpx.solvetime, mus[instance_name]);
            } else {
                slp.solve(lp, optimal, cpx.solvetime, mus[instance_name]);
            }
        }

        list.close();
        notopt.close();
    }

    return EXIT_SUCCESS;
}
