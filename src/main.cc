#include <cstdlib>
#include <iostream>
#include <vector>

#include "../include/cplex_instance.h"
#include "../include/globals.h"
#include "../include/subgradient_instance.h"

int main(int argc, char *argv[]) {
    std::string model_name = "mock.lp";

    int opt;
    while ((opt = getopt(argc, argv, "veld:i:m:T:M:")) != -1) {
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
            case 'd':
                DELIM = optarg[0];
                break;
            case 'i':
                ITERINFO = atoi(optarg);
                break;
            case 'm':
                model_name = std::string(optarg);
                break;
            case 'T':
                TIMILIMIT = atoi(optarg);
                break;
            case 'M':
                MAXITER = atoi(optarg);
                break;
        }
    }

    CplexInstance lp(model_name);
    if (VERBOSE) {
        lp.solve();
        // std::cout << "lp instance status: " << lp.getStatus() << std::endl;
    }

    SubgradientInstance slp(lp);
    slp.solve(SubgradientInstance::Methods::Hybrid);

    return EXIT_SUCCESS;
}
