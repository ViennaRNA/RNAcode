#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "extreme_fit.h"
#include "lm.h"

/* Code from HMMER package copyright S.R. Eddy  */

void Lawless416(double *x, int *y, int n, double lambda, double *ret_f, double *ret_df);
void Lawless422(double *x, int *y, int n, int z, double c, double lambda, double *ret_f, double *ret_df);
int EVDMaxLikelyFit(double *x, int *c, int n, double *ret_mu, double *ret_lambda);


/* Function: Lawless416()
 * Date:     SRE, Thu Nov 13 11:48:50 1997 [St. Louis]
 * 
 * Purpose:  Equation 4.1.6 from [Lawless82], pg. 143, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to EVD lambda parameter.
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *           
 *           Can either deal with a histogram or an array.
 *           
 *           Warning: beware overflow/underflow issues! not bulletproof.
 *           
 * Args:     x      - array of sample values (or x-axis of a histogram)
 *           y      - NULL (or y-axis of a histogram)
 *           n      - number of samples (or number of histogram bins)
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.1.6 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.1.6 evaluated at lambda
 *           
 * Return:   (void)
 */ 
void
Lawless416(double *x, int *y, int n, double lambda, double *ret_f, double *ret_df)
{

  double esum;                  /* \sum e^(-lambda xi)      */
  double xesum;                 /* \sum xi e^(-lambda xi)   */
  double xxesum;                /* \sum xi^2 e^(-lambda xi) */
  double xsum;                  /* \sum xi                  */
  double mult;                  /* histogram count multiplier */
  double total;                 /* total samples            */
  int i;


  esum = xesum = xsum  = xxesum = total = 0.;
  for (i = 0; i < n; i++)
    {
      mult = (y == NULL) ? 1. : (double) y[i];
      xsum   += mult * x[i];
      xesum  += mult * x[i] * exp(-1. * lambda * x[i]);
      xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
      esum   += mult * exp(-1. * lambda * x[i]);
      total  += mult;
    }
  *ret_f  = 1./lambda - xsum / total + xesum / esum;
  *ret_df = ((xesum / esum) * (xesum / esum))
    - (xxesum / esum)
    - (1. / (lambda * lambda));

  return;
}



/* Function: Lawless422()
 * Date:     SRE, Mon Nov 17 09:42:48 1997 [St. Louis]
 * 
 * Purpose:  Equation 4.2.2 from [Lawless82], pg. 169, and
 *           its first derivative with respect to lambda,
 *           for finding the ML fit to EVD lambda parameter
 *           for Type I censored data. 
 *           This equation gives a result of zero for the maximum
 *           likelihood lambda.
 *           
 *           Can either deal with a histogram or an array.
 *           
 *           Warning: beware overflow/underflow issues! not bulletproof.
 *           
 * Args:     x      - array of sample values (or x-axis of a histogram)
 *           y      - NULL (or y-axis of a histogram)
 *           n      - number of observed samples (or number of histogram bins)
 *           z      - number of censored samples 
 *           c      - censoring value; all observed x_i >= c         
 *           lambda - a lambda to test
 *           ret_f  - RETURN: 4.2.2 evaluated at lambda
 *           ret_df - RETURN: first derivative of 4.2.2 evaluated at lambda
 *           
 * Return:   (void)
 */ 
void
Lawless422(double *x, int *y, int n, int z, double c,
           double lambda, double *ret_f, double *ret_df)
{
  double esum;                  /* \sum e^(-lambda xi)      + z term    */
  double xesum;                 /* \sum xi e^(-lambda xi)   + z term    */
  double xxesum;                /* \sum xi^2 e^(-lambda xi) + z term    */
  double xsum;                  /* \sum xi                  (no z term) */
  double mult;                  /* histogram count multiplier */
  double total;                 /* total samples            */
  int i;

  esum = xesum = xsum  = xxesum = total = 0.;
  for (i = 0; i < n; i++)
    {
      mult = (y == NULL) ? 1. : (double) y[i];
      xsum   += mult * x[i];
      esum   += mult *               exp(-1. * lambda * x[i]);
      xesum  += mult * x[i] *        exp(-1. * lambda * x[i]);
      xxesum += mult * x[i] * x[i] * exp(-1. * lambda * x[i]);
      total  += mult;
    }

  /* Add z terms for censored data
   */
  esum   += (double) z *         exp(-1. * lambda * c);
  xesum  += (double) z * c *     exp(-1. * lambda * c);
  xxesum += (double) z * c * c * exp(-1. * lambda * c);

  *ret_f  = 1./lambda - xsum / total + xesum / esum;
  *ret_df = ((xesum / esum) * (xesum / esum))
    - (xxesum / esum)
    - (1. / (lambda * lambda));

  return;
}



/* Function: EVDMaxLikelyFit()
 * Date:     SRE, Fri Nov 14 07:56:29 1997 [St. Louis]
 * 
 * Purpose:  Given a list or a histogram of EVD-distributed samples,
 *           find maximum likelihood parameters lambda and
 *           mu. 
 *           
 * Algorithm: Uses approach described in [Lawless82]. Solves
 *           for lambda using Newton/Raphson iterations;
 *           then substitutes lambda into Lawless' equation 4.1.5
 *           to get mu. 
 *           
 *           Newton/Raphson algorithm developed from description in
 *           Numerical Recipes in C [Press88]. 
 *           
 * Args:     x          - list of EVD distributed samples or x-axis of histogram
 *           c          - NULL, or y-axis of histogram
 *           n          - number of samples, or number of histogram bins 
 *           ret_mu     : RETURN: ML estimate of mu
 *           ret_lambda : RETURN: ML estimate of lambda
 *           
 * Return:   1 on success; 0 on any failure
 */

int
EVDMaxLikelyFit(double *x, int *c, int n, double *ret_mu, double *ret_lambda)
{
  double  lambda, mu;
  double  fx;                    /* f(x)  */
  double  dfx;                   /* f'(x) */
  double esum;                  /* \sum e^(-lambda xi) */ 
  double mult;
  double total;
  double  tol = 1e-5;
  int    i;

  /* 1. Find an initial guess at lambda: linear regression here?
   */
  lambda = 0.2;

  /* 2. Use Newton/Raphson to solve Lawless 4.1.6 and find ML lambda
   */
  for (i = 0; i < 100; i++)
    {
      Lawless416(x, c, n, lambda, &fx, &dfx);
      if (fabs(fx) < tol) break;             /* success */
      lambda = lambda - fx / dfx;            /* Newton/Raphson is simple */
      if (lambda <= 0.) lambda = 0.001;      /* but be a little careful  */
    }

  /* 2.5: If we did 100 iterations but didn't converge, Newton/Raphson failed.
   *      Resort to a bisection search. Worse convergence speed
   *      but guaranteed to converge (unlike Newton/Raphson).
   *      We assume (!?) that fx is a monotonically decreasing function of x;
   *      i.e. fx > 0 if we are left of the root, fx < 0 if we
   *      are right of the root.
   */ 
  if (i == 100)
    {
      double left, right, mid;
      //SQD_DPRINTF2(("EVDMaxLikelyFit(): Newton/Raphson failed; switchover to bisection"));

                                /* First we need to bracket the root */
      lambda = right = left = 0.2;
      Lawless416(x, c, n, lambda, &fx, &dfx);
      if (fx < 0.) 
        {                       /* fix right; search left. */
          do {
            left -= 0.1;
            if (left < 0.) { 
              printf("EVDMaxLikelyFit(): failed to bracket root"); 
              return 0; 
            }
            Lawless416(x, c, n, left, &fx, &dfx);
          } while (fx < 0.);
        }
      else
        {                       /* fix left; search right. */
          do {
            right += 0.1;
            Lawless416(x, c, n, right, &fx, &dfx);
            if (right > 100.) {
              printf("EVDMaxLikelyFit(): failed to bracket root"); 
              return 0; 
            }
          } while (fx > 0.);
        }
                        /* now we bisection search in left/right interval */
      for (i = 0; i < 100; i++)
        {
          mid = (left + right) / 2.; 
          Lawless416(x, c, n, mid, &fx, &dfx);
          if (fabs(fx) < tol) break;             /* success */
          if (fx > 0.)  left = mid;
          else          right = mid;
        }
      if (i == 100) { 
        printf("EVDMaxLikelyFit(): even the bisection search failed"); 
        return 0; 
      }
      lambda = mid;
    }

  /* 3. Substitute into Lawless 4.1.5 to find mu
   */
  esum = 0.;
  total = 0.;
  for (i = 0; i < n; i++)
    {
      mult   = (c == NULL) ? 1. : (double) c[i];
      esum  += mult * exp(-1 * lambda * x[i]);
      total += mult;
    }
  mu = -1. * log(esum / total) / lambda;

  *ret_lambda = lambda;
  *ret_mu     = mu;   
  return 1;
}




