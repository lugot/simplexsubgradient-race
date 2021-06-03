#ifndef INCLUDE_GLOBALS_H_
#define INCLUDE_GLOBALS_H_

#define VERBOSE 1
#define EXTRA 0

#define TRACKING 1
#define LOG 1
#define DELIM ","
// #define DELIM " "

#define EPS 1e-9
#define INF 1e20  // TODO(lugot): TUNE to avoid double overflows

// #define UPPERBOUND std::numeric_limits<double>::max()
// #define LOWERBOUND std::numeric_limits<double>::min()
#define UPPERBOUND 1000.0
#define LOWERBOUND 0.0

#define MAX_ITERATIONS 5*1e5

#define P 20
#define TAU 1.5

#define SOFT_LAMBDA_TRESHOLD 1e-7
#define SOFT_LAMBDA_MAXLIFE 50
#define HARD_LAMBDA_TRESHOLD 1e-9

#define MU_MIN 0.001
#define MU_MAX 2.0

#define IMHERE std::cout << "--- im here ---" << std::endl;

#endif  // INCLUDE_GLOBALS_H_
