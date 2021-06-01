#include <cstdlib>
#include <iostream>
#include <vector>

#include "../include/cplex_instance.h"
#include "../include/globals.h"
#include "../include/subgradient_instance.h"

int main() {
    // TODO(lugot): IMPLEMENT better C++ getopt
    //
    // CplexInstance lp;
    // std::cout << lp.isSolved() << std::endl;

    std::string model_name = "mock.lp";
    // std::string model_name = "location.lp";
    // lp.importModel(model_name);
    CplexInstance lp(model_name);
    lp.solve();
    std::cout << "lp instance status: " << lp.getStatus() << std::endl;

    SubgradientInstance slp(lp);
    slp.solve(SubgradientInstance::Methods::Pure);

    return EXIT_SUCCESS;
}
