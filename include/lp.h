#ifndef _LP_H_
#define _LP_H_

#include "ilconcert/ilolinear.h"
#include "ilcplex/ilocplexi.h"
#include <Eigen/Sparse>

class lp {

    public:
        lp();
        //lp(std::string model_name);
        //lp(std::string model_name,
                //Eigen::SparseVector<double> c,
                //Eigen::SparseMatrix<double> A,
                //Eigen::SparseVector<double> b,
                //std::map<std::string, int> name_hasher,
                //std::map<int, std::string> var_hasher);

        std::string& getModelName();

        Eigen::SparseVector<double>& getc();
        Eigen::SparseMatrix<double>& getA();
        Eigen::SparseVector<double>& getb();

        IloCplex& getCplex();

        void importModel(std::string model_name);
        void solveSimplex();
        void solveSubgradient();

        enum projection_algorithms {
            CPLEX,
            ITERATIVE
        };

    private:
        std::string model_name_;

        int n, m;

        Eigen::SparseMatrix<double> A_;
        Eigen::SparseVector<double> b_, c_;

        bool feasible(Eigen::SparseVector<double>& u);
        Eigen::SparseVector<double> project(Eigen::SparseVector<double>& u, projection_algorithms proj_alg );


        IloEnv env_;
        IloModel model_;
        IloObjective obj_;
        IloNumVarArray var_;
        IloRangeArray rng_;
        IloCplex cplex_;

        std::map<IloNumVar, int> var_indexer;
};

// TODO: implement
std::ostream& operator<<(std::ostream& os, const lp& model);

#endif /* _LP_H_ */
