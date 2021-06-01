#ifndef INCLUDE_GLOBALS_H_
#define INCLUDE_GLOBALS_H_

#define VERBOSE 1
#define EXTRA 1
#define TRACKING 1
#define EPSILON 1e-9
#define XSMALL 1e-5

// #define UPPERBOUND std::numeric_limits<double>::max()
// #define LOWERBOUND std::numeric_limits<double>::min()
#define UPPERBOUND 500.0
#define LOWERBOUND 0.0

#define MAX_ITERATIONS 1e5

#define P 20
#define TAU 1.5

#define SOFT_LAMBDA_TRESHOLD 1e-7
#define SOFT_LAMBDA_MAXLIFE 50
#define HARD_LAMBDA_TRESHOLD 1e-9

#define IMHERE std::cout << "--- im here ---" << std::endl;

#endif  // INCLUDE_GLOBALS_H_
