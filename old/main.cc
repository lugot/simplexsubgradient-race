#include <cstdlib>
#include <iostream>
#include <vector>

#include "../include/globals.h"
#include "cplex_instance.h"

int main() {
    CplexInstance lp;
    std::cout << lp.isSolved() << std::endl;

    std::string model_name = "arod.lp";
    lp.importModel(model_name);
    lp.solve();





    return EXIT_SUCCESS;
}
