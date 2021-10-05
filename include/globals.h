#ifndef INCLUDE_GLOBALS_H_
#define INCLUDE_GLOBALS_H_

// command line specified globals
extern bool VERBOSE;
extern bool EXTRA;
extern bool LOGGING;
extern char DELIM;
extern int ITERINFO;

extern int TIMILIMIT;
extern int MAXITER;
extern int UPPERBOUND;
extern int IGNORETIMELIMIT;

// globals
#define EPS 1e-9
#define HIGHEPS 1e-8
#define INF 1e20

// #define UPPERBOUND std::numeric_limits<double>::max()
// #define LOWERBOUND std::numeric_limits<double>::min()
#define LOWERBOUND 0.0

#define P 20
// #define TAU 1.5
#define WINSIZE 7

#define SOFT_LAMBDA_TRESHOLD 1e-7
#define SOFT_LAMBDA_MAXLIFE 50
#define HARD_LAMBDA_TRESHOLD 1e-9

#define MU_MIN 0.00001
#define MU_MAX 2.0

#define SLOPE_MIN 0.75

#define STUCKTHRESHOLD 1e-6
#define STUCKTIME 100

// #define R1 100

#define MAX_NITERSTUCKED1 100
#define MAX_NITERSTUCKED2 100

#define EPS_MONOTONICITY 1e-3

#define NOTMONOTIME 20

#define VTV_LIM1 30
#define VTV_LIM2 VTV_LIM1
#define VTV_MAX 50
#define VTV_EPSILON 1e-4
#define VTV_VBAR 20

#define IMHERE std::cout << "--- im here ---" << std::endl;
#define IMHERE2 std::cout << "--- im here2 ---" << std::endl;

#endif  // INCLUDE_GLOBALS_H_
