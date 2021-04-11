#include "globals.h"
#include "parsers.h"
#include "solvers.h"

#include <cstdlib>
#include <iostream>

using namespace std;

int main() {

    lp model;
    model.importModel("prod.lp");
    //cplex_opt(model);

    int n = model.getA().cols();
    int m = model.getA().rows();

    cout << model.getA() << endl << n << " " << m << endl;

    model.solveSimplex();
    model.solveSubgradient();

    return EXIT_SUCCESS;
}
