#include <cstdlib>
#include <iostream>
#include <vector>

#include "../include/cplex_instance.h"
#include "../include/globals.h"

int main() {
    // TODO(lugot): IMPLEMENT better C++ getopt
    //
    // CplexInstance lp;
    // std::cout << lp.isSolved() << std::endl;

    std::string model_name = "prod.lp";
    // lp.importModel(model_name);
    //
    CplexInstance lp(model_name);
    lp.solve();
    lp.printStatus();

    CplexInstance::Status s = lp.getStatus();

    CplexInstance other(lp);
    other.solve();
    std::cout << other.getStatus() << std::endl;

    std::cout << s << std::endl;

    return EXIT_SUCCESS;
}
