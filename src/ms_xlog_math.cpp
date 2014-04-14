#include <math.h>
#include "ms_xlog_math.h"
#include <stdio.h>
#include <stdlib.h>

//extern double n_sums;
//extern double n_muls;
//extern double n_divs;
//extern double n_subs;

// Convert probabilities into log space, with defaults base e.
double xlog(double value)
{
	if(value == 0.0)
	{
		return(LOG_OF_ZERO);
	}
	else if (value > 0.0)
	{
		return(log(value));
	}
	else
	{
		printf("log of a negative number @ %s(%d): %.6f", __FILE__, __LINE__, value);
		int* p = NULL;
		*p = 0;
		exit(0);
	}
}

double xexp(double log_value)
{
	if(log_value == xlog(0.0))
	{
		return(0);
	}
	else
	{
		return(exp(log_value));
	}
}

double xlog_increment(double log1)
{
	return(xlog_sum(log1, 0.0));
}

double xlog_sum(double log1, double log2)
{
	//n_sums++;
	//return(log1 + log2);

	if(log1 == xlog(0.0))
	{
		return(log2);
	}
	else if(log2 == xlog(0.0))
	{
		return(log1);
	}
	else
	{
		if(log1 > log2)
		{
			return( log1 + xlog(1 + xexp(log2-log1)) );
		}
		else
		{
			return( log2 + xlog(1 + xexp(log1-log2)) );
		}

		//return(xlog((log1) + (log2)));
	}
}

// Subtract two logs, return xlog((log1) - (log2)), exit if result is negative, that is, log1 < log2.
double xlog_sub(double log1, double log2)
{
	//n_subs++;
	//return(log1 - log2);
	if(log1 < log2)
	{
		printf("%.5f < %.5f in PPF_SUB\n", log1, log2);
		printf("Cannot compute logarithm of a negative number in PPF_SUB @ %s(%d)\n", __FILE__, __LINE__);
		double* p = NULL;
		*p = 0;
		exit(0);
	}
	else if(log1 == log2)
	{
		return(xlog(0.0));
	}
	else // log1 > log2.
	{
		return( log1 + xlog(1 - xexp(log2-log1)) );
	}

	//if(xlog_comp(log1, log2))
	if(log1 == log2)
	{
		return(xlog(0.0));
	}

	if(log1 < log2)
	{
		printf("%.5f < %.5f in PPF_SUB\n", log1, log2);
		printf("Cannot compute logarithm of a negative number in PPF_SUB @ %s(%d)\n", __FILE__, __LINE__);
		double* p = NULL;
		*p = 0;
		exit(0);
	}
	else if(log1 == log2)
	{
		return(xlog(0.0));
	}
	else // log1 > log2.
	{
		return( log1 + xlog(1 - xexp(log2-log1)) );
	}
}

double xlog_mul(double log1, double log2)
{
	//n_muls++;
	//return(log1 * log2);

	if(log1 == xlog(0.0) || log2 == xlog(0.0))
	{
		return(xlog(0.0));
	}
	else
	{
		return(log1 + log2);
	}
}

// Returns 0 if log1 is 0 no matter what log2 is.
double xlog_div(double log1, double log2)
{
	//n_divs++;
	//return(log1 / log2);

	if(log1 == xlog(0.0))
	{
		return(xlog(0.0));
	}

	if(log2 == xlog(0.0))
	{
		double* p = NULL;
		*p = 0;
		printf("Division by zero error at %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
	else
	{
		//return(xlog_mul(log1, -1.0 * log2));
		return(log1 - log2);
	}
}

double xlog_pow(double log_value, double pow)
{
	if(log_value == xlog(0.0))
	{
		return(xlog(0.0));
	}
	else // Can use multiplication trick.
	{
		return(log_value * pow);
	}
}

double xlog_pow(double log_value, int pow)
{
	if(log_value == xlog(0.0))
	{
		return(xlog(0.0));
	}
	else // Can use multiplication trick.
	{
		return(log_value * (double)pow);
	}
}

bool xlog_geq(double log1, double log2)
{
	if(log1 > log2)
	{
		return(true);
	}

	if(xlog_comp(log1, log2))
	{
		return(true);
	}

	return(false);
}

double xlog_max(double log_val1, double log_val2)
{
	if(xlog_comp(log_val1, log_val2))
	{
		return(log_val1);
	}
	else if(log_val1 > log_val2)
	{
		return(log_val1);
	}
	else
	{
		return(log_val2);
	}
}

bool xlog_gt(double log1, double log2)
{
	// Do not do logarithmic check if any of the values are 0.
	if(log1 == xlog(0.0) || log2 == xlog(0.0))
	{
		if(log1 == log2)
		{
			return(false);
		}

		if(log1 == xlog(0.0))
		{
			return(false);
		}
		else
		{
			return(true);
		}
	}

	// If xlogs can be equated directly, return true.
	if(log1 == log2)
	{
		return(false);
	}

	if(log1 > log2 + XLOG_EPSILON)
	{
		return(true);
	}
	else
	{
		return(false);
	}
}

// Compare two xlog values checking for intervals based on epsilon value.
// The interpretation of epsilon based comparison is: (x/y) < ~1 & (y/x) < ~1
bool xlog_comp(double log1, double log2)
{
	// If xlogs can be equated directly, return true.
	if(log1 == log2)
	{
		return(true);
	}

	if(log1 <= log2 + XLOG_EPSILON && log1 >= log2 - XLOG_EPSILON)
	{
		return(true);
	}

	//printf("Values in PPF_COMP:\n%.25f\n%.25f\nXLOG_EPSILON=%.25f\n", log1, log2, XLOG_EPSILON);

	return(false);
}

//// compare two xlog values checking for intervals based on epsilon value.
//bool xlog_comp2(double log1, double log2)
//{
//	// if xlogs can be equated directly, return true.
//	//printf("comparing %.25f, %.25f\n", log1, log2);
//	if(log1 == log2)
//	{
//		return(true);
//	}
//
//	if(log1 > log2)
//	{
//		if((log1 != xlog(0.0f) && 
//			(xlog_div(xlog_sub(log1, log2), log1)) <= -3.0f &&
//			(log2 != xlog(0.0f) && 
//			(xlog_div(xlog_sub(log1, log2), log2)) <= -3.0f))
//		{
//			return(true);
//		}
//		else
//		{
//			return(false);
//		}
//	}
//	else if(log2 > log1)
//	{
//		if((log1 != xlog(0.0f) && 
//			(xlog_div(xlog_sub(log2, log1), log1)) <= -3.0f) &&
//			(log2 != xlog(0.0f) && 
//			(xlog_div(xlog_sub(log2, log1), log2)) <= -3.0f))
//		{
//			return(true);
//		}
//		else
//		{
//			return(false);
//		}
//	}
//
//	printf("encountered the thing that should not be @ %s(%d) for comparison of %.25f %.25f\n", __FILE__, __LINE__, log1, log2);
//	exit(0);
//
//	//if(log1 <= log2 + xlog_epsilon && log1 >= log2 - xlog_epsilon)
//	//{
//	//	return(true);
//	//}
//
//	//printf("values in ppf_comp:\n%.25f\n%.25f\nxlog_epsilon=%.25f\n", log1, log2, xlog_epsilon);
//
//	return(false);
//}

// return log1 > log2.
// This function is not very dependable, that is, 
// it is not very logical to think of comparison this way,
// so it should be used cautiously.
bool is_bigger(double log1, double log2)
{
	return(log1 > (log2 + XLOG_EPSILON));
}

// Get the precision of comparison function.
double get_xlog_comp_prec()
{
	double log_prec = xlog(1.0f);
	double log_one = xlog(1.0f);
	while(1)
	{
		if(xlog_comp(log_one, xlog_sum(log_one, log_prec)))
		{
			printf("%lf = %lf + %G\n", log_one, log_one, log_prec);
			break;
		}

		log_prec = xlog_div(log_prec, xlog(2.0f));
	}

	return(log_prec);
}



