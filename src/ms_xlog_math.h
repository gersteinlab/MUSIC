#ifndef _XLOG_MATH_
#define _XLOG_MATH_

#include "math.h"

//#define EPS (0.00001f)

// Extended logarithmic/exponential functions from Mann 2006's work, numerically stable hidden markov model.
// Includes a log value for log of 0, implements PPF_SUMmation, product for implementation of forward-backward algorithm.

#define LOG_OF_ZERO (-1.0 * pow(2.0, 30)) // Log of zero: -2^(30), can be represented exactly.

// Following is the epsilon value that defines the interval to assume two xlog values are equal.
#define XLOG_EPSILON (pow(0.1, 10)) // 10^(-8), assuming exp(-8) is very low, this epsilon interval makes sure a reliable comparison.
//#define XLOG_EPSILON (-10.0) // 10^(-8), assuming exp(-8) is very low, this epsilon interval makes sure a reliable comparison.

#define MAX(x, y) ((x)>(y)?(x):(y))
#define MIN(x, y) ((x)>(y)?(y):(x))

// Convert probabilities into log space, with defaults base e.
double xlog(double prob);

double xlog_increment(double log1);

double xlog_sum(double log1, double log2);

double xlog_mul(double log1, double log2);

double xlog_div(double log1, double log2);

double xlog_sub(double log1, double log2);

double xlog_pow(double log_value, double power);

double xlog_pow(double log_value, int power);

double xlog_max(double log_val1, double log_val2);

double xexp(double log_value);

// returns log1 > log2
bool is_bigger(double log1, double log2);

// returns log1 > log2
bool is_smaller(double log1, double log2);

// Compare xlogs using intervals.
bool xlog_comp(double log1, double log2);

bool xlog_geq(double log1, double log2);
bool xlog_gt(double log1, double log2);

double get_xlog_comp_prec();

#endif // _XLOG_MATH_


