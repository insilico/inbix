//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <errno.h>
#include <stdio.h>
#include <time.h>

#include "stats.h"

// #ifdef WITH_LAPACK
// #include "lapackf.h"
// #endif

#include "helper.h"
#include "crandom.h"
#include "options.h"
#include "plink.h"
#include "perm.h"
#include "dcdflib.h"
#include "ipmpar.h"

#include "model.h"
#include "linear.h"
#include "logistic.h"

#define FPMIN 1.0e-30

extern ofstream LOG;
extern Plink * PP;

bool realnum(double d) {
  double zero = 0;
  if(d != d || d == 1 / zero || d == -1 / zero)
    return false;
  else
    return true;
}

long double factorial(int x) {
  int i;
  long double result = 1;
  for(i = 2; i <= x; i++)
    result *= i;
  return result;
}

double normdist(double z) {
  double sqrt2pi = 2.50662827463;
  double t0, z1, p0;
  t0 = 1 / (1 + 0.2316419 * fabs(z));
  z1 = exp(-0.5 * z * z) / sqrt2pi;
  p0 = z1 * t0
          * (0.31938153 +
          t0 * (-0.356563782 +
          t0 * (1.781477937 +
          t0 * (-1.821255978 +
          1.330274429 * t0))));
  return z >= 0 ? 1 - p0 : p0;
}

double chiprobP(double x, double df) {

  if(!realnum(x)) return -9;

  double p, q;
  int st = 0; // error variable
  int w = 1; // function variable
  double bnd = 1; // boundary function

  // NCP is set to 0
  cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

  // Check status
  if(st != 0) return -9;

  // Return p-value
  return q;
}

double inverse_chiprob(double q, double df) {

  if(!realnum(q)) return -9;
  else if(q >= 1) return 0;

  double x;
  double p = 1 - q;
  int st = 0; // error variable
  int w = 2; // function variable
  double bnd = 1; // boundary function

  // NCP is set to 0
  cdfchi(&w, &p, &q, &x, &df, &st, &bnd);

  // Check status
  if(st != 0) return -9;

  // Return p-value
  return x;
}

double gammln(double xx) {
  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146, -86.50532032941677,
    24.01409824083091, -1.231739572450155,
    0.1208650973866179e-2, -0.5395239384953e-5};
  int j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for(j = 0; j <= 5; j++) ser += cof[j] / ++y;
  return -tmp + log(2.5066282746310005 * ser / x);
}

// Inverse normal distribution

/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 */

/* Coefficients in rational approximations. */

static const double a[] = {
  -3.969683028665376e+01,
  2.209460984245205e+02,
  -2.759285104469687e+02,
  1.383577518672690e+02,
  -3.066479806614716e+01,
  2.506628277459239e+00
};

static const double b[] = {
  -5.447609879822406e+01,
  1.615858368580409e+02,
  -1.556989798598866e+02,
  6.680131188771972e+01,
  -1.328068155288572e+01
};

static const double c[] = {
  -7.784894002430293e-03,
  -3.223964580411365e-01,
  -2.400758277161838e+00,
  -2.549732539343734e+00,
  4.374664141464968e+00,
  2.938163982698783e+00
};

static const double d[] = {
  7.784695709041462e-03,
  3.224671290700398e-01,
  2.445134137142996e+00,
  3.754408661907416e+00
};

#define LOW 0.02425
#define HIGH 0.97575

double ltqnorm(double p) {
  double q, r;

  errno = 0;

  if(p < 0 || p > 1) {
    return 0.0;
  } else if(p == 0) {
    return -HUGE_VAL /* minus "infinity" */;
  } else if(p == 1) {
    return HUGE_VAL /* "infinity" */;
  } else if(p < LOW) {
    /* Rational approximation for lower region */
    q = sqrt(-2 * log(p));
    return(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if(p > HIGH) {
    /* Rational approximation for upper region */
    q = sqrt(-2 * log(1 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else {
    /* Rational approximation for central region */
    q = p - 0.5;
    r = q*q;
    return(((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5])*q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }
}

double pT(double T, double df) {

  if(!realnum(T))
    return -9;

  T = abs(T);

  double p, q;
  int st = 0; // error variable
  int w = 1; // function variable
  double bnd = 1; // boundary function

  // NCP is set to 0
  cdft(&w, &p, &q, &T, &df, &st, &bnd);

  // Check status
  if(st != 0) return -9;

  // Return two-sided p-value
  return 2 * q;
}

double pF(const double F, const int df1, const int df2) {
  return betai(0.5 * df2, 0.5 * df1, (double) df2 / (double) (df2 + df1 * F));
}

double betai(const double a, const double b, const double x) {
  double bt;

  if(x < 0.0 || x > 1.0) error("Internal error: bad x in routine betai");
  if(x == 0.0 || x == 1.0) bt = 0.0;
  else
    bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
  if(x < (a + 1.0) / (a + b + 2.0))
    return bt * betacf(a, b, x) / a;
  else
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

double betacf(const double a, const double b, const double x) {
  const int MAXIT = 100;
  const double EPS = 3e-7;

  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;

  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  c = 1.0;
  d = 1.0 - qab * x / qap;
  if(fabs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  h = d;
  for(m = 1; m <= MAXIT; m++) {
    m2 = 2 * m;
    aa = m * (b - m) * x / ((qam + m2)*(a + m2));
    d = 1.0 + aa*d;
    if(fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if(fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    h *= d*c;
    aa = -(a + m)*(qab + m) * x / ((a + m2)*(qap + m2));
    d = 1.0 + aa*d;
    if(fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if(fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    del = d*c;
    h *= del;
    if(fabs(del - 1.0) <= EPS) break;
  }
  if(m > MAXIT) error("Internal error in betacf() function (please report)");
  return h;
}

vector< vector<double> > inverse(vector< vector<double> > & m) {
  double d;
  int i, j;

  if(m.size() == 0) error("Internal error: matrix with no rows (inverse function)");
  if(m.size() != m[0].size()) error("Internal error: cannot invert non-square matrix");
  int n = m.size();

  // indx is an integer array
  vector<int> indx(n);

  vector<double> col(n);
  vector<vector<double> > y(n);
  for(int i = 0; i < n; i++) y[i].resize(n);
  vector<vector<double> > tm;
  tm = m;

  ludcmp(tm, indx, d);

  for(j = 0; j < n; j++) {
    for(i = 0; i < n; i++) col[i] = 0;
    col[j] = 1;
    lubksb(tm, indx, col);
    for(i = 0; i < n; i++) y[i][j] = col[i];
  }

  return y;
}

vector<double> eigenvalues(vector<vector<double> > & a) {

  // 'a' should be a square, symmetric matrix
  int n = a.size();
  vector<double> e(n);
  vector<double> d(n);
  tred2(a, d, e);
  vector<vector<double> > z; // dummy
  tqli(d, e, z);
  return d;
}

// Householder method to reduce real, symmetric matrix
// to tridiagonal form
// Modified to return only eigenvalues.

void tred2(vector<vector<double> > & a,
        vector<double> & d,
        vector<double> &e) {
  int l, k, j, i;
  double scale, hh, h, g, f;

  int n = d.size();
  for(i = n - 1; i > 0; i--) {
    l = i - 1;
    h = scale = 0.0;
    if(l > 0) {
      for(k = 0; k < l + 1; k++)
        scale += fabs(a[i][k]);
      if(scale == 0.0)
        e[i] = a[i][l];
      else {
        for(k = 0; k < l + 1; k++) {
          a[i][k] /= scale;
          h += a[i][k] * a[i][k];
        }
        f = a[i][l];
        g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale*g;
        h -= f*g;
        a[i][l] = f - g;
        f = 0.0;
        for(j = 0; j < l + 1; j++) {
          // Next statement can be omitted if eigenvectors not wanted
          // 	  a[j][i]=a[i][j]/h;
          g = 0.0;
          for(k = 0; k < j + 1; k++)
            g += a[j][k] * a[i][k];
          for(k = j + 1; k < l + 1; k++)
            g += a[k][j] * a[i][k];
          e[j] = g / h;
          f += e[j] * a[i][j];
        }
        hh = f / (h + h);
        for(j = 0; j < l + 1; j++) {
          f = a[i][j];
          e[j] = g = e[j] - hh*f;
          for(k = 0; k < j + 1; k++)
            a[j][k] -= (f * e[k] + g * a[i][k]);
        }
      }
    } else
      e[i] = a[i][l];
    d[i] = h;
  }
  // Next statement can be omitted if eigenvectors not wanted
  //   d[0]=0.0;
  e[0] = 0.0;
  // Contents of this loop can be omitted if eigenvectors not
  //	wanted except for statement d[i]=a[i][i];
  for(i = 0; i < n; i++) {
    //     l=i;
    //     if (d[i] != 0.0) {
    //       for (j=0;j<l;j++) {
    // 	g=0.0;
    // 	for (k=0;k<l;k++)
    // 	  g += a[i][k]*a[k][j];
    // 	for (k=0;k<l;k++)
    // 	  a[k][j] -= g*a[k][i];
    //       }
    //     }
    d[i] = a[i][i];
    //     a[i][i]=1.0;
    //     for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
  }
}

// Modified to return only eigenvalues.
void tqli(vector<double> &d, vector<double>&e,
        vector<vector<double> > &z) {
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  double volatile temp;
  int n = d.size();
  for(i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0.0;
  for(l = 0; l < n; l++) {
    iter = 0;
    do {
      for(m = l; m < n - 1; m++) {
        dd = fabs(d[m]) + fabs(d[m + 1]);
        temp = fabs(e[m]) + dd;
        if(temp == dd) break;
      }
      if(m != l) {
        if(iter++ == 30) error("Internal problem in tqli routine");
        g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
        s = c = 1.0;
        p = 0.0;
        for(i = m - 1; i >= l; i--) {
          f = s * e[i];
          b = c * e[i];
          e[i + 1] = (r = pythag(f, g));
          if(r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c*b;
          d[i + 1] = g + (p = s * r);
          g = c * r - b;
          // Next loop can be omitted if eigenvectors not wanted
          /* for (k=0;k<n;k++) {
             f=z[k][i+1];
             z[k][i+1]=s*z[k][i]+c*f;
             z[k][i]=c*z[k][i]-s*f;
             } */
        }
        if(r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while(m != l);
  }
}

//////////////////////////////////////////////////
// As above, but with eigenvectors returned also
Eigen eigenvectors(vector<vector<double> > & a) {
  // 'a' should be a square, symmetric matrix
  int n = a.size();

  Eigen E;
  E.set(n);

  vector<double> e(n, 0);
  EV_tred2(a, E.d, e);
  EV_tqli(E.d, e, a);
  E.z = a;
  return E;
}

// Householder method to reduce real, symmetric matrix
// to tridiagonal form
// Modified to return both eigenvalues and eigenvectors
void EV_tred2(vector<vector<double> > & a,
        vector<double> & d,
        vector<double> &e) {
  int l, k, j, i;
  double scale, hh, h, g, f;

  int n = d.size();
  for(i = n - 1; i > 0; i--) {
    l = i - 1;
    h = scale = 0.0;
    if(l > 0) {
      for(k = 0; k < l + 1; k++)
        scale += fabs(a[i][k]);
      if(scale == 0.0)
        e[i] = a[i][l];
      else {
        for(k = 0; k < l + 1; k++) {
          a[i][k] /= scale;
          h += a[i][k] * a[i][k];
        }
        f = a[i][l];
        g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale*g;
        h -= f*g;
        a[i][l] = f - g;
        f = 0.0;
        for(j = 0; j < l + 1; j++) {
          a[j][i] = a[i][j] / h;
          g = 0.0;
          for(k = 0; k < j + 1; k++)
            g += a[j][k] * a[i][k];
          for(k = j + 1; k < l + 1; k++)
            g += a[k][j] * a[i][k];
          e[j] = g / h;
          f += e[j] * a[i][j];
        }
        hh = f / (h + h);
        for(j = 0; j < l + 1; j++) {
          f = a[i][j];
          e[j] = g = e[j] - hh*f;
          for(k = 0; k < j + 1; k++)
            a[j][k] -= (f * e[k] + g * a[i][k]);
        }
      }
    } else
      e[i] = a[i][l];
    d[i] = h;
  }

  d[0] = 0.0;
  e[0] = 0.0;

  for(i = 0; i < n; i++) {
    l = i;
    if(d[i] != 0.0) {
      for(j = 0; j < l; j++) {
        g = 0.0;
        for(k = 0; k < l; k++)
          g += a[i][k] * a[k][j];
        for(k = 0; k < l; k++)
          a[k][j] -= g * a[k][i];
      }
    }
    d[i] = a[i][i];
    a[i][i] = 1.0;
    for(j = 0; j < l; j++) a[j][i] = a[i][j] = 0.0;
  }
}

// Modified to return eigenvalues and eigenvectors
void EV_tqli(vector<double> &d, vector<double>&e,
        vector<vector<double> > &z) {
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  int n = d.size();
  for(i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0.0;
  for(l = 0; l < n; l++) {
    iter = 0;
    do {
      for(m = l; m < n - 1; m++) {
        dd = fabs(d[m]) + fabs(d[m + 1]);
        if(fabs(e[m]) + dd == dd) break;
      }
      if(m != l) {
        if(iter++ == 30) error("Internal problem in tqli routine");
        g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
        s = c = 1.0;
        p = 0.0;
        for(i = m - 1; i >= l; i--) {
          f = s * e[i];
          b = c * e[i];
          e[i + 1] = (r = pythag(f, g));
          if(r == 0.0) {
            d[i + 1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c*b;
          d[i + 1] = g + (p = s * r);
          g = c * r - b;

          for(k = 0; k < n; k++) {
            f = z[k][i + 1];
            z[k][i + 1] = s * z[k][i] + c*f;
            z[k][i] = c * z[k][i] - s*f;
          }
        }
        if(r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while(m != l);
  }
}

/////////////////////////
// Romberg integration
double qromb(double func(const double), double a, double b) {
  const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
  const double EPS = 1.0e-10;
  double ss, dss;
  vector_t s(JMAX), h(JMAXP), s_t(K), h_t(K);
  int i, j;

  h[0] = 1.0;
  for(j = 1; j <= JMAX; j++) {
    s[j - 1] = trapzd(func, a, b, j);
    if(j >= K) {
      for(i = 0; i < K; i++) {
        h_t[i] = h[j - K + i];
        s_t[i] = s[j - K + i];
      }
      polint(h_t, s_t, 0.0, ss, dss);
      if(fabs(dss) <= EPS * fabs(ss)) return ss;
    }
    h[j] = 0.25 * h[j - 1];
  }
  error("Internal error: too many steps in routine qromb");
  return 0.0;
}

void polint(vector_t &xa, vector_t &ya, const double x, double &y, double &dy) {
  int i, m, ns = 0;
  double den, dif, dift, ho, hp, w;

  int n = xa.size();
  vector_t c(n), d(n);
  dif = fabs(x - xa[0]);
  for(i = 0; i < n; i++) {
    if((dift = fabs(x - xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  y = ya[ns--];
  for(m = 1; m < n; m++) {
    for(i = 0; i < n - m; i++) {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      if((den = ho - hp) == 0.0) error("Error in routine polint");
      den = w / den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    y += (dy = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
  }

}

double trapzd(double func(const double), const double a, const double b, const int n) {
  double x, tnm, sum, del;
  static double s;
  int it, j;

  if(n == 1) {
    return(s = 0.5 * (b - a)*(func(a) + func(b)));
  } else {
    for(it = 1, j = 1; j < n - 1; j++) it <<= 1;
    tnm = it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    for(sum = 0.0, j = 0; j < it; j++, x += del) sum += func(x);
    s = 0.5 * (s + (b - a) * sum / tnm);
    return s;
  }
}

/////////////////////////
void svdvar(vector<vector<double> > & v,
        vector<double> & w,
        vector<vector<double> > & cvm) {
  int i, j, k;
  double sum;

  int ma = w.size();
  vector<double> wti(ma);
  for(i = 0; i < ma; i++) {
    wti[i] = 0.0;
    if(w[i] != 0.0) wti[i] = 1.0 / (w[i] * w[i]);
  }
  for(i = 0; i < ma; i++) {
    for(j = 0; j < i + 1; j++) {
      sum = 0.0;
      for(k = 0; k < ma; k++)
        sum += v[i][k] * v[j][k] * wti[k];
      cvm[j][i] = cvm[i][j] = sum;
    }
  }
}

void svbksb(vector<vector<double> > &u,
        vector<double> &w,
        vector<vector<double> > &v,
        vector<double> &b,
        vector<double> &x) {
  int jj, j, i;
  double s;

  //   int us = u.size()>0 ? u[0].size() : 0;
  //   int vs = v.size()>0 ? v[0].size() : 0;
  //   cout << "U = " << u.size() << " " << us<< "\n";
  //   cout << "V = " << v.size() << " " << vs << "\n";
  //   cout << "w = " << w.size() << "\n";
  //   cout << "b = " << b.size() << "\n";
  //   cout << "x = " << x.size() << "\n";

  int m = u.size();
  int n = u[0].size();
  vector<double> tmp(n);
  for(j = 0; j < n; j++) {
    s = 0.0;
    if(w[j] != 0.0) {
      for(i = 0; i < m; i++) s += u[i][j] * b[i];
      s /= w[j];
    }
    tmp[j] = s;
  }

  for(j = 0; j < n; j++) {
    s = 0.0;
    for(jj = 0; jj < n; jj++) s += v[j][jj] * tmp[jj];
    x[j] = s;
  }

}

bool svd(matrix_t & u, vector_t &w, matrix_t &v) {

  // #ifdef WITH_LAPACK
  //   matrix_t u2;
  //   svd_lapack(u,w,u2,v);
  //   u = u2;
  //#else

  const double eps = 1e-12;

  if(u.size() == 0)
    error("Internal problem: matrix with no rows in svd()");

  int r = u.size();
  int c = u[0].size();
  w.resize(c);
  sizeMatrix(v, c, c);

  bool flag = svdcmp(u, w, v);

  return flag;

  // Look for singular values
  //   double wmax = 0;
  //   for (int i=0; i<n; i++)
  //     wmax = w[i] > wmax ? w[i] : wmax;
  //   double wmin = wmax * eps;
  //   for (int i=0; i<n; i++)
  //     {
  //       w[i] = w[i] < wmin ? 0 : 1/w[i];
  //     }  

  // #endif
}

vector< vector<double> > svd_inverse(vector< vector<double> > & u, bool & flag) {

  const double eps = 1e-24;

  if(u.size() == 0)
    error("Internal problem: matrix with no rows (inverse function)");
  if(u.size() != u[0].size())
    error("Internal problem: Cannot invert non-square matrix");
  int n = u.size();

  vector<double> w(n, 0);

  vector<vector<double> > v(n);
  for(int i = 0; i < n; i++)
    v[i].resize(n, 0);

  flag = svdcmp(u, w, v);

  // Look for singular values
  double wmax = 0;
  for(int i = 0; i < n; i++)
    wmax = w[i] > wmax ? w[i] : wmax;
  double wmin = wmax * eps;
  for(int i = 0; i < n; i++) {
    w[i] = w[i] < wmin ? 0 : 1 / w[i];
  }

  // u w t(v)

  // row U * 1/w

  // results matrix
  vector<vector<double> > r(n);
  for(int i = 0; i < n; i++) {
    r[i].resize(n, 0);
    for(int j = 0; j < n; j++)
      u[i][j] = u[i][j] * w[j];
  }

  // [nxn].[t(v)] 
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      for(int k = 0; k < n; k++)
        r[i][j] += u[i][k] * v[j][k];

  return r;
}

// Matrix square root function
vector<vector<double> > msqrt(vector<vector<double> > & u) {

  // Using SVD, square root is U . sqrt(D) . V_T

  //  msqrt <- function(m) {
  //  m <- svd(m)
  //  m$u %*% sqrt(diag(m$d)) %*% t(m$v) }

  const double eps = 1e-12;

  int n = u.size();
  vector<double> d(n, 0);
  vector<vector<double> > v(n);
  for(int i = 0; i < n; i++)
    v[i].resize(n, 0);

  svdcmp(u, d, v);

  // Take square root of diagonal values
  for(int i = 0; i < n; i++)
    d[i] = sqrt(d[i]);

  // Multiple to reconstruct original

  vector<vector<double> > r(n);
  for(int i = 0; i < n; i++)
    r[i].resize(n, 0);

  vector<vector<double> > r2 = r;

  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      r[i][j] = u[i][j] * d[j];

  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      for(int k = 0; k < n; k++)
        r2[i][j] += r[i][k] * v[j][k];

  return r2;
}

void ludcmp(vector<vector<double> > &a, vector<int> &indx, double &d) {
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  int n = a.size();
  vector<double> vv(n);
  d = 1;

  for(i = 0; i < n; i++) {
    big = 0;
    for(j = 0; j < n; j++)
      if((temp = fabs(a[i][j])) > big) big = temp;
    if(big == 0) error("singular matrix in ludcmp");
    vv[i] = 1 / big;
  }

  for(j = 0; j < n; j++) {
    for(i = 0; i < j; i++) {
      sum = a[i][j];
      for(k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0;
    for(i = j; i < n; i++) {
      sum = a[i][j];
      for(k = 0; k < j; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if((dum = vv[i] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if(j != imax) {
      for(k = 0; k < n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a[j][j] == 0) a[j][j] = 1.0e-20;

    if(j != n - 1) {
      dum = 1 / (a[j][j]);
      for(i = j + 1; i < n; i++) a[i][j] *= dum;
    }
  }
}

void lubksb(vector<vector<double> > &a, vector<int> &indx, vector<double> &b) {

  int i, ii = 0, ip, j;
  double sum;

  int n = a.size();

  for(i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if(ii != 0)
      for(j = ii - 1; j < i; j++) sum -= a[i][j] * b[j];
    else if(sum != 0.0) ii = i + 1;
    b[i] = sum;
  }
  for(i = n - 1; i >= 0; i--) {
    sum = b[i];
    for(j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

bool svdcmp(vector<vector<double> > & a,
        vector<double> & w,
        vector<vector<double> > &v) {
  bool flag;
  int i, its, j, jj, k, l, nm;
  double anorm, c, f, g, h, s, scale, x, y, z;
  double volatile temp;

  int m = a.size();
  if(m == 0) error("Internal problem in SVD function (no observations left?)");
  int n = a[0].size();

  vector<double> rv1(n);
  g = scale = anorm = 0.0;
  for(i = 0; i < n; i++) {
    l = i + 2;
    rv1[i] = scale*g;
    g = s = scale = 0.0;
    if(i < m) {
      for(k = i; k < m; k++) scale += fabs(a[k][i]);
      if(scale != 0.0) {
        for(k = i; k < m; k++) {
          a[k][i] /= scale;
          s += a[k][i] * a[k][i];
        }
        f = a[i][i];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        a[i][i] = f - g;
        for(j = l - 1; j < n; j++) {
          for(s = 0.0, k = i; k < m; k++) s += a[k][i] * a[k][j];
          f = s / h;
          for(k = i; k < m; k++) a[k][j] += f * a[k][i];
        }
        for(k = i; k < m; k++) a[k][i] *= scale;
      }
    }
    w[i] = scale *g;
    g = s = scale = 0.0;
    if(i + 1 <= m && i + 1 != n) {
      for(k = l - 1; k < n; k++) scale += fabs(a[i][k]);
      if(scale != 0.0) {
        for(k = l - 1; k < n; k++) {
          a[i][k] /= scale;
          s += a[i][k] * a[i][k];
        }
        f = a[i][l - 1];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        a[i][l - 1] = f - g;
        for(k = l - 1; k < n; k++) rv1[k] = a[i][k] / h;
        for(j = l - 1; j < m; j++) {
          for(s = 0.0, k = l - 1; k < n; k++) s += a[j][k] * a[i][k];
          for(k = l - 1; k < n; k++) a[j][k] += s * rv1[k];
        }
        for(k = l - 1; k < n; k++) a[i][k] *= scale;
      }
    }
    anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }
  for(i = n - 1; i >= 0; i--) {
    if(i < n - 1) {
      if(g != 0.0) {
        for(j = l; j < n; j++)
          v[j][i] = (a[i][j] / a[i][l]) / g;
        for(j = l; j < n; j++) {
          for(s = 0.0, k = l; k < n; k++) s += a[i][k] * v[k][j];
          for(k = l; k < n; k++) v[k][j] += s * v[k][i];
        }
      }
      for(j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  for(i = MIN(m, n) - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];
    for(j = l; j < n; j++) a[i][j] = 0.0;
    if(g != 0.0) {
      g = 1.0 / g;
      for(j = l; j < n; j++) {
        for(s = 0.0, k = l; k < m; k++) s += a[k][i] * a[k][j];
        f = (s / a[i][i]) * g;
        for(k = i; k < m; k++) a[k][j] += f * a[k][i];
      }
      for(j = i; j < m; j++) a[j][i] *= g;
    } else for(j = i; j < m; j++) a[j][i] = 0.0;
    ++a[i][i];
  }
  for(k = n - 1; k >= 0; k--) {
    for(its = 0; its < 30; its++) {
      flag = true;
      for(l = k; l >= 0; l--) {
        nm = l - 1;
        temp = fabs(rv1[l]) + anorm;
        if(temp == anorm) {
          flag = false;
          break;
        }
        temp = fabs(w[nm]) + anorm;
        if(temp == anorm) break;
      }
      if(flag) {
        c = 0.0;
        s = 1.0;
        for(i = l; i < k + 1; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          temp = fabs(f) + anorm;
          if(temp == anorm) break;
          g = w[i];
          h = pythag(f, g);
          w[i] = h;
          h = 1.0 / h;
          c = g*h;
          s = -f*h;
          for(j = 0; j < m; j++) {
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y * c + z*s;
            a[j][i] = z * c - y*s;
          }
        }
      }
      z = w[k];
      if(l == k) {
        if(z < 0.0) {
          w[k] = -z;
          for(j = 0; j < n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if(its == 29)
        return false; // cannot converge: multi-collinearity?
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z)*(y + z)+(g - h)*(g + h)) / (2.0 * h * y);
      g = pythag(f, 1.0);
      f = ((x - z)*(x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
      c = s = 1.0;
      for(j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s*g;
        g = c*g;
        z = pythag(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g*s;
        g = g * c - x*s;
        h = y*s;
        y *= c;
        for(jj = 0; jj < n; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z*s;
          v[jj][i] = z * c - x*s;
        }
        z = pythag(f, h);
        w[j] = z;
        if(z) {
          z = 1.0 / z;
          c = f*z;
          s = h*z;
        }
        f = c * g + s*y;
        x = c * y - s*g;
        for(jj = 0; jj < m; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z*s;
          a[jj][i] = z * c - y*s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  return true;
}

double pythag(const double a, const double b) {
  double absa, absb;

  absa = fabs(a);
  absb = fabs(b);
  if(absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
  else return(absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

double SQR(double a) {
  return a*a;
}

void multMatrix(vector<vector<double> > & a,
        vector<vector<double> > & b,
        vector<vector<double> > & c) {

  int ar = a.size();
  int br = b.size();
  if(ar == 0 || br == 0)
    error("Internal error: multiplying 0-sized matrices");

  int ac = a[0].size();
  int bc = b[0].size();
  if(ac != br)
    error("Internal error: non-conformable matrices in multMatrix()");

  int cr = ar;
  int cc = bc;

  c.clear();
  sizeMatrix(c, cr, cc);

  for(int i = 0; i < ar; i++)
    for(int j = 0; j < bc; j++)
      for(int k = 0; k < ac; k++)
        c[i][j] += a[i][k] * b[k][j];
}

double symTable(table_t t) {
  // Test for symmetry in a n x n table
  int a = t.size();
  if(a == 0)
    return -1;
  int b = t[0].size();
  if(a != b)
    return -1;

  int df = a * (a - 1) / 2;
  double x = 0;
  for(int i = 0; i < a; i++)
    for(int j = 0; j < i; j++) {
      double tmp = t[i][j] - t[j][i];
      tmp *= tmp;
      tmp /= t[i][j] + t[j][i];
      x += tmp;
    }
  return chiprobP(x, df);
}

double chiTable(table_t t) {
  int a = t.size();
  if(a == 0)
    return -1;
  int b = t[0].size();

  vector_t rows(a);
  vector_t cols(b);
  double sum = 0;

  for(int i = 0; i < a; i++)
    for(int j = 0; j < b; j++) {
      rows[i] += t[i][j];
      cols[j] += t[i][j];
      sum += t[i][j];
    }

  // Sum (O-E)^2 / E 
  matrix_t exp;
  sizeMatrix(exp, a, b);
  double chisq = 0;
  for(int i = 0; i < a; i++)
    for(int j = 0; j < b; j++) {
      exp[i][j] += (rows[i] * cols[j]) / sum;
      double tmp = t[i][j] - exp[i][j];
      tmp *= tmp;
      tmp /= exp[i][j];
      chisq += tmp;
    }
  return chisq;
}

double chi2x2(table_t t) {
  return chi2x2(t[0][0], t[0][1], t[1][0], t[1][1]);
}

double chi2x2(matrix_t t) {
  return chi2x2(t[0][0], t[0][1], t[1][0], t[1][1]);
}

double chi2x2(double a, double b, double c, double d) {
  double row1 = a + b;
  double row2 = c + d;
  double col1 = a + c;
  double col2 = b + d;

  if(row1 * row2 * col1 * col2 == 0)
    return 0;

  double total = col1 + col2;

  double E_a = (row1 * col1) / total;
  double E_b = (row1 * col2) / total;
  double E_c = (row2 * col1) / total;
  double E_d = (row2 * col2) / total;

  return((a - E_a)*(a - E_a)) / E_a +
          ((b - E_b)*(b - E_b)) / E_b +
          ((c - E_c)*(c - E_c)) / E_c +
          ((d - E_d)*(d - E_d)) / E_d;
}

int pca(matrix_t & x,
        boolmatrix_t & mask,
        vector_t & p,
        matrix_t & s,
        matrix_t & v,
        bool mean_centre = true) {

  // Missing values in bmask (F)

  // Populate s with scores ( U.D ) for significant PCs

  // Center each column
  // SVD -> U.D.Vt
  // Return UD, and # of PCs we should look at (0 means error)
  // Handle missing data by mean imputation (i.e. set to 0 after centering)

  // Edit g to equal U.W.V'  (after any editing)
  int nrow = x.size();
  if(nrow == 0)
    return 0;
  int ncol = x[0].size();

  vector_t means(ncol);
  vector<int> cnt(nrow, 0);

  for(int r = 0; r < nrow; r++) {
    for(int c = 0; c < ncol; c++) {
      if(!mask[r][c]) {
        means[c] += x[r][c];
        ++cnt[c];
      }
    }
  }

  for(int c = 0; c < ncol; c++)
    means[c] /= (double) cnt[c];

  // Center on column means
  if(mean_centre) {
    for(int r = 0; r < nrow; r++)
      for(int c = 0; c < ncol; c++) {
        if(mask[r][c])
          x[r][c] = 0;
        else
          x[r][c] -= means[c];
      }
  } else {
    // ALTERNATE: no mean centering
    for(int r = 0; r < nrow; r++)
      for(int c = 0; c < ncol; c++) {
        if(mask[r][c])
          x[r][c] = means[c];
      }
  }

  // Perform SVD on X
  vector_t p2;
  svd(x, p2, v);

  // Figure out how many component to return
  // Use Dunn & Everitt (2001) rule of s^2 above 0.7/n
  double thresh = 0.7 / (double) ncol;
  double totvar = 0;

  vector<int> keep;
  for(int i = 0; i < ncol; i++)
    totvar += p2[i] * p2[i];

  for(int i = 0; i < ncol; i++) {
    if((p2[i] * p2[i]) / totvar >= thresh) {
      p2[i] = 1;
      keep.push_back(i);
    } else
      p2[i] = 0;
  }

  // What to return? If in 2sided mode, just return all 
  // PC scores that meet criterion

  // Return PC scores in s, PCs in v
  matrix_t w = vec2diag(p2);
  matrix_t z, z2;

  // S = X %*% P
  multMatrix(x, w, z);

  if(par::elf_pcmode_2sided) {
    // Calculate scores, then return
    sizeMatrix(s, nrow, keep.size());

    for(int r = 0; r < nrow; r++)
      for(int c = 0; c < keep.size(); c++) {
        s[r][c] = z[r][keep[c]];
      }

    p.resize(keep.size());
    for(int c = 0; c < keep.size(); c++)
      p[c] = (p2[keep[c]] * p2[keep[c]]) / totvar;

    //      cout << "sizes = " << keep.size() << " " << s[0].size() << " " << ncol << "\n";
    return keep.size();

  }

  // Otherwise, reconstruct X as U.W.V'
  if(!par::elf_pcmode_2sided) {

    multMatrix(z, v, x);

    // For now return everything --- add back in the generic PCA 
    // later -- but for now, we will use this special version just
    // for ELF calculations

    return ncol;
  }

  // Otherwise, returned pruned x
  multMatrix(z, v, z2);

  sizeMatrix(x, nrow, keep.size());

  for(int r = 0; r < nrow; r++)
    for(int c = 0; c < keep.size(); c++) {
      x[r][c] = z2[r][keep[c]];
    }

  p.resize(keep.size());
  for(int c = 0; c < keep.size(); c++)
    p[c] = (p2[keep[c]] * p2[keep[c]]) / totvar;

  return keep.size();
}

matrix_t vec2diag(vector_t & v) {
  matrix_t d;
  sizeMatrix(d, v.size(), v.size());
  for(int i = 0; i < v.size(); i++)
    d[i][i] = v[i];
  return d;
}

double rnorm() {
  double u1 = CRandom::rand();
  double u2 = CRandom::rand();
  return sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
  // z2 = sqrt(-2*log(u1)) * sin(2*M_PI*u2);
}

// ----------------------------------------------------------------------------
// some handy matrix functions
// bcw - 5/10/13
// ----------------------------------------------------------------------------

bool matrixComputeCovariance(matrix_t X, matrix_t& covMatrix, 
        matrix_t& corMatrix) {

  int n = X.size();

	cout << "matrixComputeCovariance, n = " << n << endl;
	
  // compute covariances
	PP->printLOG("Computing covariance\n");
  matrix_t P;
  sizeMatrix(P, n, n);
  matrixFill(P, 1);
  matrixDivideScalar(P, n);
  vector_t one(n, 1);
  matrix_t diag = vec2diag(one);
  matrix_t Q;
  sizeMatrix(Q, n, n);
  matrixSubtract(diag, P, Q);
  matrix_t xStar;
  multMatrix(Q, X, xStar);
  matrix_t xStarT;
  matrixTranspose(xStar, xStarT);
  multMatrix(xStarT, xStar, covMatrix);
  matrixDivideScalar(covMatrix, n - 1);

  // compute correlations from covariances
	PP->printLOG("Computing correlation\n");
  vector_t covDiag;
  matrixGetDiag(covMatrix, covDiag);
  // standard deviations from variances
  for(int i=0; i < covDiag.size(); ++i) {
    covDiag[i] = 1.0 / sqrt(covDiag[i]);
  }
  matrix_t D = vec2diag(covDiag);
  matrix_t temp;
  multMatrix(D, covMatrix, temp);
  multMatrix(temp, D, corMatrix);
  
  return true;
}

bool matrixFill(matrix_t& m, double f) {
  int i, j;
  for(i=0; i < m.size(); ++i) {
    for(j=0; j < m[0].size(); ++j) {
      m[i][j] = f;
    }
  }
  return true;
}

bool matrixMultiplyScalar(matrix_t& m, double s) {
  int i, j;
  for(i=0; i < m.size(); ++i) {
    for(j=0; j < m[0].size(); ++j) {
      m[i][j] *= s;
    }
  }
  return true;
}

bool matrixDivideScalar(matrix_t& m, double s) {
  int i, j;
  for(i=0; i < m.size(); ++i) {
    for(j=0; j < m[0].size(); ++j) {
      m[i][j] /= s;
    }
  }
  return true;
}

bool matrixAdd(matrix_t m1, matrix_t m2, matrix_t& result) {
  int m1Rows = m1.size();
  int m2Rows = m2.size();
  int m1Cols = m1[0].size();
  int m2Cols = m2[0].size();

  if(m1Rows != m2Rows) {
    error("Error in matrixSubtract, rows sizes are not compatible");
  }
  
  if(m1Cols != m2Cols) {
    error("Error in matrixSubtract, rows sizes are not compatible");
  }
  
  int i, j;
  for(i=0; i < m1Rows; ++i) {
    for(j=0; j < m1Cols; ++j) {
      result[i][j] = m1[i][j] + m2[i][j];
    }
  }
  return true;
}

bool matrixSubtract(matrix_t m1, matrix_t m2, matrix_t& result) {
  int m1Rows = m1.size();
  int m2Rows = m2.size();
  int m1Cols = m1[0].size();
  int m2Cols = m2[0].size();

  if(m1Rows != m2Rows) {
    error("Error in matrixSubtract, rows sizes are not compatible");
  }
  
  if(m1Cols != m2Cols) {
    error("Error in matrixSubtract, rows sizes are not compatible");
  }

  int i, j;
  for(i=0; i < m1Rows; ++i) {
    for(j=0; j < m1Cols; ++j) {
      result[i][j] = m1[i][j] - m2[i][j];
    }
  }
  return true;
}

bool matrixTranspose(matrix_t in, matrix_t& out) {
  int inRows = in.size();
  int inCols = in[0].size();
  sizeMatrix(out, inCols, inRows);
  
  int i, j;
  for(i=0; i < inRows; ++i) {
    for(j=0; j < inCols; ++j) {
      out[j][i] = in[i][j];
    }
  }
  return true;
}

bool matrixGetDiag(matrix_t m, vector_t& d) {
  d.clear();
  for(int i=0; i < m.size(); ++i) {
    d.push_back(m[i][i]);
  }
  return true;
}

bool matrixSetDiag(matrix_t& m, vector_t d) {
  for(int i=0; i < d.size(); ++i) {
    m[i][i] = d[i];
  }
  return true;
}

bool matrixRead(string mFilename, matrix_t& m, vector<string>& variableNames) {
    // open the numeric attributes file if possible
  checkFileExists(mFilename);
  ifstream matrixFile(mFilename.c_str(), ios::in);
  if(matrixFile.fail()) {
    return false;
  }

  bool readHeader = false;
  int rows = 0;
  int cols = 0;
  while(!matrixFile.eof()) {

    char nline[par::MAX_LINE_LENGTH];
    matrixFile.getline(nline, par::MAX_LINE_LENGTH, '\n');

    // convert to string
    string sline = nline;
    if(sline == "") continue;

    // read line from text file into a vector of tokens
    string buf;
    stringstream ss(sline);
    vector<string> tokens;
    while(ss >> buf) {
      tokens.push_back(buf);
    }

    // parse header if not parsed already
    if(!readHeader) {
      // save numeric attribute names = tokens minus FID and IID
      cols = tokens.size() - 2;
      variableNames.resize(cols);
      copy(tokens.begin() + 2, tokens.end(), variableNames.begin());
      readHeader = true;
      continue;
    } else {
      if(tokens.size() != (cols + 2)) {
        matrixFile.close();
        cerr << "Line:\n" + sline + "\n";
        return false;
      }
    }
    
    ++rows;
    
    // Add numeric attribute values to data matrix
    vector_t dataValues;
    bool okay = true;
    dataValues.clear();
    for(int c = 2; c < cols + 2; c++) {
      double t = 0;
      if(!from_string<double>(t, tokens[c], std::dec))
        okay = false;
      dataValues.push_back(t);
    }
    if(okay) {
			m.push_back(dataValues);
    }
    else {
      cerr << "Error reading data values from line " << rows << endl;
      return false;
    }

  }
  matrixFile.close();

  PP->printLOG("Read matrix from [" + mFilename + "]: " + 
  int2str(rows) + " rows x " + int2str(cols) + " columns\n");
    
  return true;
}

bool matrixWrite(matrix_t m, string mFilename, 
        vector<string> variableNames) {
  PP->printLOG("Writing matrix [ " + mFilename + " ]\n");
  ofstream outFile(mFilename.c_str());
  if(outFile.fail()) {
    return false;
  }
  outFile.precision(6);
  outFile.fixed;

  // write header
  int hIdx = 0;
  for(vector<string>::const_iterator hIt = variableNames.begin();
          hIt != variableNames.end(); ++hIt, ++hIdx) {
    if(hIdx) {
      outFile << "\t" << *hIt;
    }
    else {
      outFile << *hIt;
    }
  }
  outFile << endl;

  // write the covariance matrix to a file
  for(int i=0; i < m.size(); ++i) {
    for(int j=0; j < m[i].size(); ++j) {
      if(j) {
        outFile << "\t" << m[i][j];
      }
      else {
        outFile << m[i][j];
      }
    }
    outFile << endl;
  }
  
  outFile.close();
  
  return true;
}

bool matrixSums(matrix_t m, vector_t& sums, int dim) {
  sums.clear();
  if(dim == 0) {
    sums.resize(m.size(), 0);
  } else {
    sums.resize(m[0].size(), 0);
  }
  for(int i=0; i < m.size(); ++i) {
    for(int j=0; j < m[0].size(); ++j) {
      if(dim == 0) {
        sums[j] += m[i][j];
      } else {
        sums[i] += m[i][j];
      }
    }
  }

  return true;
}

bool matrixComputeNodeDegrees(matrix_t a, vector_t& ad) {
  if(a.size() != a[0].size()) {
    error("Cannot compute degrees of non-square network matrix");
  }
  ad.resize(a.size(), 0);
  for(int i=0; i < a.size(); ++i) {
    for(int j=0; j < a[0].size(); ++j) {
      if(a[i][j] != 0) {
        ad[i] += 1;
      }
    }
  }
  
  return true;
}

bool matrixExtractRowColIdx(matrix_t m, intvec_t rowIdx, intvec_t colIdx, 
        matrix_t& nm) {
  int nr = rowIdx.size();
  int nc = colIdx.size();
  
  sizeMatrix(nm, nr, nc);
  
  for(int i=0; i < nr; ++i) {
    for(int j=0; j < nc; ++j) {
      nm[i][j] = m[rowIdx[i]][colIdx[j]];
    }
  }
  
  return true;
}

bool matrixGetTrace(matrix_t m, double& t) {
  t = 0;
  for(int i=0; i < m.size(); ++i) {
    t += m[i][i];
  }
  return true;
}

bool matrixElementWiseMultiply(matrix_t m, matrix_t n, matrix_t& m_out) {
  int mRow = m.size();
  int mCol = m[0].size();
  int nRow = n.size();
  int nCol = n[0].size();
  if((mRow == nRow) && (mCol == nCol)) {
    sizeMatrix(m_out, nRow, nCol);

    int i, j;
    for(i=0; i < nRow; ++i) {
      for(j=0; j < nCol; ++j) {
        m_out[i][j] = m[i][j] * n[i][j];
      }
    }
  } else {
    error("matrixElementWiseMultiply has incompatible matrices");
  }
  
  return true;
}

bool matrixConnectivityThreshold(matrix_t& m, double t, bool binary) {
	// threshold the adjacency matrix
  for(unsigned int i=0; i < m.size(); ++i) {
    for(unsigned int j=0; j < m[0].size(); ++j) {
      if(m[i][j] <= t) {
        m[i][j] = 0.0;
      }
      else {
        if(binary) {
          m[i][j] = 1.0;
        }
      }
    }
  }
  
  return true;
}

bool reportNumericSummaryStats() {
  int numInd = PP->sample.size();
  int numVar = PP->nlistname.size();
  
  string outFilename = par::output_file_name + ".numeric.summary.ind.tab";
  PP->printLOG("Writing individuals summary [ " + outFilename + " ]\n");
  ofstream outFileI(outFilename.c_str());
  if(outFileI.fail()) {
    return false;
  }
  outFileI << "FID\tIID\tMean\tVar\tSd" << endl;
  for(int i=0; i < numInd; i++) {
    vector_t values;
    for(int j=0; j < numVar; ++j) {
      values.push_back(PP->sample[i]->nlist[j]);
    }
    vector_t summary;
    vectorSummary(values, summary);
    double coefOfVar = summary[2] / summary[0];
    outFileI 
            << PP->sample[i]->fid << "\t" 
            << PP->sample[i]->iid << "\t"
            << summary[0] << "\t" 
            << summary[1] << "\t" 
            << summary[2] 
            << endl;
  }
  outFileI.close();
  
  outFilename = par::output_file_name + ".numeric.summary.var.tab";
  PP->printLOG("Writing variables summary [ " + outFilename + " ]\n");
  ofstream outFileV(outFilename.c_str());
  if(outFileV.fail()) {
    return false;
  }
  outFileV << "Variable\tMean\tVar\tSd\tCV" << endl;
  for(int i=0; i < numVar; i++) {
    vector_t values;
    for(int j=0; j < numInd; ++j) {
      values.push_back(PP->sample[j]->nlist[i]);
    }
    vector_t summary;
    vectorSummary(values, summary);
    double coefOfVar = summary[2] / summary[0];
    outFileV
            << PP->nlistname[i] << "\t"
            << summary[0] << "\t" 
            << summary[1] << "\t" 
            << summary[2] << "\t"
            << coefOfVar
            << endl;
  }
  outFileV.close();
  
  return true;
}

bool vectorSummary(vector_t values, vector_t& summary) {
  double n = values.size();
  if(!n) {
    return false;
  }
  double sum = 0;
  for(int i=0; i < values.size(); ++i) {
    sum += values[i];
  }
  double mean = sum / n;
  double SSE = 0;
  for(int i=0; i < values.size(); ++i) {
    double deviation = values[i] - mean;
    SSE += deviation * deviation;
  }
  double var = SSE / (values.size() - 1);
  double sd = sqrt(var);

  summary.resize(3);
  summary[0] = mean;
  summary[1] = var;
  summary[2] = sd;
}

bool rankByRegression(RegressionRankType rankType, rankedlist_t& ranks,
				RegressionRankResults& results) {

  // should the sorted results be reversed, i.e., by p-value
  bool sortDescending = false;
  double mainEffectValueP = 1.0;
  double mainEffectValueB = 0.0;
  double mainEffectValueS = 0.0;
  double rankValue = 0;
  
  // model SNPs
  for(int i=0; i < PP->nl_all; ++i) {
    string coefLabel = PP->locus[i]->name;
    Model *mainEffectModel;
    if(par::bt) {
      mainEffectModel = new LogisticModel(PP);
    } else {
      mainEffectModel = new LinearModel(PP);
    }
  	mainEffectModel->setMissing();
    pair<double, double> snpResult;
    mainEffectModel->addAdditiveSNP(i);
    mainEffectModel->label.push_back(coefLabel);

    // add covariates if specified
    if(par::covar_file) {
      for(int i = 0; i < par::clist_number; i++) {
        // add covariate to the model
        mainEffectModel->addCovariate(i);
        mainEffectModel->label.push_back(PP->clistname[i]);
      }
    }
    
    // Build design matrix
    mainEffectModel->buildDesignMatrix();

    // Fit linear model
    int tp = 1; 
    mainEffectModel->testParameter = tp; // single variable main effect
    mainEffectModel->fitLM();
    
    mainEffectValueP = 1.0;
    mainEffectValueB = 0.0;
    mainEffectValueS = 0.0;
    rankValue = 0;
    if(!mainEffectModel->isValid()) {
      string failMsg = "WARNING: Invalid main effect regression fit for variable [" +
        coefLabel + "]";
      PP->printLOG(failMsg + "\n");
    } else {
      // Obtain estimates and statistics
      vector_t betaMainEffectCoefs = mainEffectModel->getCoefs();
      // p-values don't include intercept term
      vector_t betaMainEffectCoefPvals = mainEffectModel->getPVals();
      mainEffectValueP = betaMainEffectCoefPvals[tp - 1];
      vector_t mainEffectModelSE = mainEffectModel->getSE();
      // always use first coefficient after intercept as main effect term
      mainEffectValueB = betaMainEffectCoefs[tp];
      // t- or z-statistic
      mainEffectValueS = betaMainEffectCoefs[tp] / mainEffectModelSE[tp];
    }

    switch(rankType) {
      case REGRESSION_RANK_STAT:
        rankValue = mainEffectValueS;
        sortDescending = true;
        break;
      case REGRESSION_RANK_BETA:
        rankValue = mainEffectValueB;
        sortDescending = true;
        break;
      case REGRESSION_RANK_PVAL:
        // rankValue = -log10(mainEffectValueP);
        rankValue = mainEffectValueP;
        break;
    }

    ranks.push_back(make_pair(rankValue, coefLabel));

    results.vars.push_back(coefLabel);
    results.coefs.push_back(mainEffectValueB);
    results.pvals.push_back(mainEffectValueP);
    results.stats.push_back(mainEffectValueS);

    delete mainEffectModel;
  }

  // model numerics
  for(int i=0; i < PP->nlistname.size(); ++i) {
    string coefLabel = PP->nlistname[i];
    Model *mainEffectModel;
    if(par::bt) {
      mainEffectModel = new LogisticModel(PP);
    } else {
      mainEffectModel = new LinearModel(PP);
    }
  	mainEffectModel->setMissing();
    pair<double, double> numResult;
    mainEffectModel->addNumeric(i);
    mainEffectModel->label.push_back(coefLabel);

    // add covariates if specified
    if(par::covar_file) {
      for(int i = 0; i < par::clist_number; i++) {
        // add covariate to the model
        mainEffectModel->addCovariate(i);
        mainEffectModel->label.push_back(PP->clistname[i]);
      }
    }

    // Build design matrix
    mainEffectModel->buildDesignMatrix();

    // Fit linear model
    int tp = 1; 
    mainEffectModel->testParameter = tp; // single variable main effect
    mainEffectModel->fitLM();

    mainEffectValueP = 1.0;
    mainEffectValueB = 0.0;
    mainEffectValueS = 0.0;
    rankValue = 0;
    if(!mainEffectModel->isValid()) {
      string failMsg = "WARNING: Invalid main effect regression fit for variable [" +
        coefLabel + "]";
      PP->printLOG(failMsg + "\n");
    } else {
      // Obtain estimates and statistics
      vector_t betaMainEffectCoefs = mainEffectModel->getCoefs();
      // p-values don't include intercept term
      vector_t betaMainEffectCoefPvals = mainEffectModel->getPVals();
      mainEffectValueP = betaMainEffectCoefPvals[tp - 1];
      vector_t mainEffectModelSE = mainEffectModel->getSE();

      // always use first coefficient after intercept as main effect term
      mainEffectValueB = betaMainEffectCoefs[tp];
  		// t- or z-statistic
      mainEffectValueS = betaMainEffectCoefs[tp] / mainEffectModelSE[tp];
    }

    switch(rankType) {
      case REGRESSION_RANK_STAT:
        rankValue = mainEffectValueS;
        sortDescending = true;
        break;
      case REGRESSION_RANK_BETA:
        rankValue = mainEffectValueB;
        sortDescending = true;
        break;
      case REGRESSION_RANK_PVAL:
        // rankValue = -log10(mainEffectValueP);
        rankValue = mainEffectValueP;
        break;
    }

    ranks.push_back(make_pair(rankValue, coefLabel));

		results.vars.push_back(coefLabel);
		results.coefs.push_back(mainEffectValueB);
		results.pvals.push_back(mainEffectValueP);
		results.stats.push_back(mainEffectValueS);
		
    delete mainEffectModel;
  }

  // sort results in reverse order
	sort(ranks.begin(), ranks.end());
  if(sortDescending) {
    reverse(ranks.begin(), ranks.end());
  }
  
  return true;
}

pair<double, double> fitModel(Model* mainEffectModel) {
  pair<double, double> result;

  // Build design matrix
  mainEffectModel->buildDesignMatrix();

	// Fit linear model
	int tp = 1; 
	mainEffectModel->testParameter = tp; // single variable main effect
	mainEffectModel->fitLM();

	// Obtain estimates and statistics
	vector_t betaMainEffectCoefs = mainEffectModel->getCoefs();
	// p-values don't include intercept term
	vector_t betaMainEffectCoefPvals = mainEffectModel->getPVals();
  double mainEffectPval = betaMainEffectCoefPvals[tp - 1];
	vector_t mainEffectModelSE = mainEffectModel->getSE();

  // always use first coefficient after intercept as main effect term
	// double mainEffectValue = betaMainEffectCoefs[tp];
	double mainEffectValue = betaMainEffectCoefs[tp] /
						mainEffectModelSE[tp];
  
  result.first = mainEffectValue;
  result.second = mainEffectPval;
  
  return result;
}

bool numericLowValueFilter(double percentile, boolvec_t& varFlags) {
  // get all values into a vector
  vector_t allValues;
  for(int i=0; i < PP->nlistname.size(); ++i) {
    for(int j=0; j < PP->sample.size(); ++j) {
      allValues.push_back(PP->sample[j]->nlist[i]);
    }
  }

  // find percentile cutoff threshold
  sort(allValues.begin(), allValues.end());
  int thresholdIndex = floor(percentile * (double) allValues.size());
  double threshold = allValues[thresholdIndex];
  PP->printLOG("Detecting values below: " + dbl2str(threshold) + "\n");
  
  // for each numeric attribute, if all values below threshold, remove it
  int numSamples = PP->sample.size();
  int numPass = 0;
  int numFail = 0;
  for(int i=0; i < PP->nlistname.size(); ++i) {
    int lowCount = 0;
    for(int j=0; j < numSamples; ++j) {
      if(PP->sample[j]->nlist[i] < threshold) {
        ++lowCount;
      }
    }
    if(lowCount == numSamples) {
//      cout << "numericLowValueFilter: " << PP->nlistname[i] 
//        << " failed low value test" << endl;
      ++numFail;
      varFlags[i] = false;
    }
    else {
      ++numPass;
    }
  }

  cout << "numericLowValueFilter: " 
    << numPass << " passed, " << numFail << " failed" << endl;
  
  return true;
}

bool numericVarianceFilter(double percentile, boolvec_t& varFlags) {
  vector_t variances;
  for(int i=0; i < varFlags.size(); ++i) {
    if(varFlags[i]) {
      vector_t varValues;
      for(int j=0; j < PP->sample.size(); ++j) {
        varValues.push_back(PP->sample[j]->nlist[i]);
      }
      vector_t summary;
      vectorSummary(varValues, summary);
      variances.push_back(summary[1]);
    }
  }
  
  // find percentile cutoff threshold
  vector_t variancesCopy(variances.size());
  copy(variances.begin(), variances.end(), variancesCopy.begin());
  sort(variancesCopy.begin(), variancesCopy.end());
  int thresholdIndex = floor(percentile * (double) variancesCopy.size());
  double threshold = variancesCopy[thresholdIndex];
  PP->printLOG("Detecting variance below: " + dbl2str(threshold) + "\n");

  // for each numeric attribute, if all values below threshold, remove it
  int numFail = 0;
  int numPass = 0;
  for(int i=0; i < variances.size(); ++i) {
    if(variances[i] < threshold) {
//      cout << "numericLowVarianceFilter: " << PP->nlistname[i] 
//        << " failed low variance test with variance: " 
//        << variances[i] << endl;
      ++numFail;
      varFlags[i] = false;
    } else {
      ++numPass;
    }
  }
  
  cout << "numericLowVarianceFilter: " 
    << numPass << " passed, " << numFail << " failed" << endl;

  return true;
}

bool quantile(vector_t values, double percentile, double& percentileValue) {

  if(values.size() < 1) {
    cerr << "quantile error: no values" << endl;
    return false;
  }

  sort(values.begin(), values.end());
  int percentileIndex = floor(percentile * (double) values.size());
  percentileValue = values[percentileIndex];
  
  return true;
}

// data set transforms prior to analysis - bcw - 10/30/13
bool numericMeanCenter() {
  for(int i=0; i < PP->nlistname.size(); ++i) {
    double varSum = 0;
    for(int j=0; j < PP->sample.size(); ++j) {
      varSum += PP->sample[j]->nlist[i];
    }
    double varAvg = varSum / (double) PP->sample.size();
    for(int j=0; j < PP->sample.size(); ++j) {
      PP->sample[j]->nlist[i] -= varAvg;
    }
  }  
  return true;
}

bool numericStandardize() {
  for(int i=0; i < PP->nlistname.size(); ++i) {
    vector_t varVals;
    for(int j=0; j < PP->sample.size(); ++j) {
      varVals.push_back(PP->sample[j]->nlist[i]);
    }
    vector_t summary;
    vectorSummary(varVals, summary);
    if(summary[2]) {
      for(int j=0; j < PP->sample.size(); ++j) {
        PP->sample[j]->nlist[i] /= summary[2];
      }
    }
    else {
      cout << "WARNING: standard deviation 0 in numericStandardize" 
        << " for variable index: " << i << endl;
    }
  }  
  return true;
}

// get cases and controls into vectors for a numeric variable index
bool getNumericCaseControl(int varIndex, vector_t& cases, vector_t& controls) {
	// determine the number of affected and unaffected individuals
	int nAff = 0;
	int nUnaff = 0;
	for(int i=0; i < PP->sample.size(); i++) {
		if(PP->sample[i]->aff) {
			++nAff;
		}
		else {
 			if(PP->sample[i]->missing) {
  			// cerr << "Missing phenotype(s) for numeric data" << endl;
        return false;
  		} else {
        ++nUnaff;
      }
		}
	}

	// size return vectors
	cases.resize(nAff);
	controls.resize(nUnaff);
	
	// load numerics into passed matrices
	int aIdx = 0;
	int uIdx = 0;
	for(int i=0; i < PP->sample.size(); i++) {
			if(PP->sample[i]->aff) {
				cases[aIdx] = PP->sample[i]->nlist[varIndex];
			} else {
			  // handle missing phenos - bcw - 2/5/15
  			if(!PP->sample[i]->missing) {
  				controls[uIdx] = PP->sample[i]->nlist[varIndex];
  			}
			}
		if(PP->sample[i]->aff) {
			++aIdx;
		}
		else {
			if(!PP->sample[i]->missing) {
  			++uIdx;
  		}
		}
	}
  
  return true;
}

// t-test of two groups - bcw - 10/30/13
bool tTest(int varIndex, double& t) {
  // get the group data values into vectors
  vector_t g1_data;
  vector_t g2_data;
  getNumericCaseControl(varIndex, g1_data, g2_data);
  if(g1_data.size() == 0) {
    cerr << "Group 1 size = 0 in tTest" << endl;
    return false;
  }
  double g1_n = (double) g1_data.size();
  if(g2_data.size() == 0) {
    cerr << "Group 2 size = 0 in tTest" << endl;
    return false;
  }
  double g2_n = (double) g2_data.size();

  // gets summary information for each group (case - control)
  vector_t g1_summary;
  vectorSummary(g1_data, g1_summary);
  double g1_mean = g1_summary[0];
  double g1_var = g1_summary[1];
  vector_t g2_summary;
  vectorSummary(g2_data, g2_summary);
  double g2_mean = g2_summary[0];
  double g2_var = g2_summary[1];
  t = 
    (g2_mean - g1_mean) / 
    sqrt((g1_var / g1_n) + (g2_var / g2_n));

  return true;
}

// z-test of two groups - bcw - 7/31/14
bool zTest(int varIndex, double& z) {
  // get the group data values into vectors
  vector_t g1_data;
  vector_t g2_data;
  getNumericCaseControl(varIndex, g1_data, g2_data);
  if(g1_data.size() == 0) {
    cerr << "Group 1 size = 0 in zTest" << endl;
    return false;
  }
  double g1_n = (double) g1_data.size();
  if(g2_data.size() == 0) {
    cerr << "Group 2 size = 0 in zTest" << endl;
    return false;
  }
  double g2_n = (double) g2_data.size();

  // gets summary information for each group (case - control)
  vector_t g1_summary;
  vectorSummary(g1_data, g1_summary);
  double g1_mean = g1_summary[0];
  double g1_var = g1_summary[1];
  double z_i_1 = 0.5 * log((abs((1 + g1_mean) / (1 - g1_mean))));

  vector_t g2_summary;
  vectorSummary(g2_data, g2_summary);
  double g2_mean = g2_summary[0];
  double g2_var = g2_summary[1];
  double z_i_2 = 0.5 * log((abs((1 + g2_mean) / (1 - g2_mean))));

  z = abs(z_i_1 - z_i_2) / sqrt((1/(g1_n - 3) + 1 / (g2_n - 3)));

  if(std::isinf(z) || std::isnan(z)) {
    error("Infinite or NaN in Z-test zTest() calculation");
  }

  return true;
}
