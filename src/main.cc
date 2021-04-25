#include <cstdlib>
#include <iostream>
#include <vector>

#include "../include/cplex_instance.h"
#include "../include/globals.h"
#include "subgradient_instance.h"

int main() {
    // TODO(lugot): IMPLEMENT better C++ getopt
    //
    // CplexInstance lp;
    // std::cout << lp.isSolved() << std::endl;

    std::string model_name = "prod.lp";
    // lp.importModel(model_name);
    CplexInstance lp(model_name);
    CplexInstance qp(lp);
    std::cout << "lp instance status: " << lp.getStatus() << std::endl;
    lp.solve();
    lp.getStatus();

    std::cout << "qp instance status: " << qp.getStatus() << std::endl;
    qp.solve();
    qp.getStatus();

    SubgradientInstance slp(lp);
    slp.solve();

    return EXIT_SUCCESS;
}
