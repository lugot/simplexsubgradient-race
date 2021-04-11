#include "ilconcert/ilolinear.h"
#include "lp.h"
#include "globals.h"
#include <ilcplex/ilocplex.h>
#include <stdexcept>

using namespace std;
using namespace Eigen;


bool feasible(SparseMatrix<double>& A, SparseVector<double>& b, SparseVector<double>& x);
static void printObjective (IloObjective obj);

void lp::solveSimplex() {

    //https://www.ibm.com/docs/en/icos/20.1.0?topic=application-constructing-environment-iloenv

    try {
        //https://www.ibm.com/docs/en/icos/20.1.0?topic=application-creating-model-ilomodel

        IloEnv env = model_.getEnv();

        if (cplex_.solve() == IloFalse) {
           env_.error() << "Failed to optimize LP" << endl;
           throw(-1);
        }
        IloNumArray vals(env_);
        env.out() << "Solution status = " << cplex_.getStatus() << endl;
        env.out() << "Solution value  = " << cplex_.getObjValue() << endl;
        cplex_.getValues(vals, var_);
        env.out() << "Values        = " << vals << endl;
        cplex_.getSlacks(vals, rng_);
        env.out() << "Slacks        = " << vals << endl;
        cplex_.getDuals(vals, rng_);
        env.out() << "Duals         = " << vals << endl;
        cplex_.getReducedCosts(vals, var_);
        env.out() << "Reduced Costs = " << vals << endl;

         env.end();

     }
     catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
     }
     catch (...) {
        cerr << "Unknown exception caught" << endl;
     }


}

void lp::solveSubgradient() {
    SparseMatrix<double>& A = A_;
    SparseVector<double>& b = b_;
    SparseVector<double>& c = c_;

    int n = A.cols();
    int m = A.rows();

    SparseVector<double> u;
    u.resize(n);
    u.coeffRef(0) = 0.0;



    cout << "feasible: " << feasible(u) << endl;
    u = u-c;
    RowVectorXd dx(u);
    cout << "x: " << dx << endl;
    cout << "feasible: " << feasible(u) << endl;

    SparseVector<double> uplus = project(u, CPLEX);
    RowVectorXd du(uplus);
    cout << "u: " << du << endl;


    //while (feasible(A, b, x)) {
        //x = x-c;
        //VectorXd dx(x);
        //cout << "x: " << dx << endl;
    //}



}

bool lp::feasible(SparseVector<double>& u) {
    // TODO: maybe this vector can be avoided
    SparseVector<double> res = A_*u;

    SparseVector<double>::InnerIterator itb(b_);
    for (SparseVector<double>::InnerIterator it(res); it; ++it) {
        while (itb && itb.index() < it) ++itb;

        double rhs = (it.index() == itb.index()) ? itb.value() : 0;
        if (rhs - it.value() > EPS) return false;
    }

    return true;
}

SparseVector<double> lp::project(Eigen::SparseVector<double>& u, lp::projection_algorithms proj_alg ) {

    if (proj_alg != lp::CPLEX) throw invalid_argument("invalid projection algoritmh");


    cout << endl << "solving qp" << endl;

    int n = u.rows();
    SparseVector<double> uplus(n);

    //IloEnv env = var_.getEnv();
    //obj_.end();
    //obj_ = IloObjective(env_);
    //
    printObjective(obj_);
    obj_.setQuadCoef(var_[0], var_[0], 0);
    printObjective(obj_);
    try {
    cout << n << " as " << var_.getSize() << endl;
    obj_.setLinearCoef(var_[0], 0);
    cout << n << " as " << var_.getSize() << endl;
        //for (int i=0; i<var_.getSize(); i++) {
            //obj_.setQuadCoef(var_[i], var_[i], 0.5);
            //obj_.setLinearCoef(var_[i], -u.coeff(i));
        //}

        //model_.add(obj);
        // Optimize the problem and obtain solution.
        if (!cplex_.solve()) {
         //env.error() << "Failed to optimize LP" << endl;
         throw(-1);
        }

        //IloNumArray vals(env);
        //env.out() << "Solution status = " << cplex.getStatus() << endl;
        //env.out() << "Solution value  = " << cplex.getObjValue() << endl;
        //cplex.getValues(vals, var);
        //env.out() << "Values        = " << vals << endl;
        //cplex.getSlacks(vals, con);
        //env.out() << "Slacks        = " << vals << endl;
        //cplex.getDuals(vals, con);
        //env.out() << "Duals         = " << vals << endl;
        //cplex.getReducedCosts(vals, var);
        //env.out() << "Reduced Costs = " << vals << endl;

   }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }


    return uplus;
}

static void
printObjective (IloObjective obj)
{
   IloEnv env = obj.getEnv();

   env.out() << "obj: " << obj << endl;

   // Count the number of linear terms
   // in the objective function.
   IloInt nlinterms = 0;
   for (IloExpr::LinearIterator lit = obj.getLinearIterator(); lit.ok(); ++lit) {
      ++nlinterms;
   }

   // Count the number of quadratic terms
   // in the objective function.
   IloInt nquadterms = 0;
   IloInt nquaddiag  = 0;
   for (IloExpr::QuadIterator qit = obj.getQuadIterator(); qit.ok(); ++qit) {
      ++nquadterms;
      if ( qit.getVar1().getImpl() == qit.getVar2().getImpl() )
         ++nquaddiag;
   }

   env.out() << "number of linear terms in the objective             : " << nlinterms  << endl;
   env.out() << "number of quadratic terms in the objective          : " << nquadterms << endl;
   env.out() << "number of diagonal quadratic terms in the objective : " << nquaddiag  << endl;
   env.out() << endl;
} // END printObjective

//TODO: REMOVE
//TODO: refactoring
//void cplex_opt(lp sparse_model) {

    //int m = sparse_model.getA().rows();
    //int n = sparse_model.getA().cols();

    ////https://www.ibm.com/docs/en/icos/20.1.0?topic=application-constructing-environment-iloenv
    //IloEnv env;

    //try {
        ////https://www.ibm.com/docs/en/icos/20.1.0?topic=application-creating-model-ilomodel
        //IloModel model(env);
        //IloNumVarArray var(env);
        //IloRangeArray con(env, m, -IloInfinity, 0);

        ////TODO: sicuro esiste un metodo che fa tutto one shot
        //// default constructor of IloNumVar: lb=0, ub=IloInfinity
        //for (int i=0; i<n; i++) var.add(IloNumVar(env));

        //// add ub of nonzero constraints
        //for (SparseVector<double>::InnerIterator it(sparse_model.getb()); it; ++it) {
            ////TODO: strange incremental iterator warning
            //con[it.index()].setUB(it.value());
        //}
        //// add nonzero coefficent constraints
        //for (int k=0; k<sparse_model.getA().outerSize(); ++k) {
            //for (SparseMatrix<double>::InnerIterator it(sparse_model.getA(), k); it; ++it) {
                //con[it.row()].setLinearCoef(var[it.col()], it.value());
            //}
        //}

        //// add objective
        //IloObjective obj = IloMinimize(env);
        //for (SparseVector<double>::InnerIterator it(sparse_model.getc()); it; ++it) {
            ////TODO: strange incremental iterator warning
            //obj.setLinearCoef(var[it.index()], it.value());
        //}

        //model.add(obj);
        //model.add(con);

        ////https://www.ibm.com/docs/en/icos/20.1.0?topic=application-solving-model-ilocplex
        //IloCplex cplex(model);

        //if (cplex.solve() == IloFalse) {
           //env.error() << "Failed to optimize LP" << endl;
           //throw(-1);
        //}

        //IloNumArray vals(env);
        //env.out() << "Solution status = " << cplex.getStatus() << endl;
        //env.out() << "Solution value  = " << cplex.getObjValue() << endl;
        //cplex.getValues(vals, var);
        //env.out() << "Values        = " << vals << endl;
        //cplex.getSlacks(vals, con);
        //env.out() << "Slacks        = " << vals << endl;
        //cplex.getDuals(vals, con);
        //env.out() << "Duals         = " << vals << endl;
        //cplex.getReducedCosts(vals, var);
        //env.out() << "Reduced Costs = " << vals << endl;


     //}
     //catch (IloException& e) {
        //cerr << "Concert exception caught: " << e << endl;
     //}
     //catch (...) {
        //cerr << "Unknown exception caught" << endl;
     //}



     //env.end();


    //return ;
//}
