//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#ifndef __STATS_H__
#define __STATS_H__

#include <string>
#include <vector>
#include <cstdio>

#include "plink.h"
#include "model.h"

using namespace std;

void sizeMatrix(matrix_t &, int, int);
void sizeTable(table_t &, int, int);
void multMatrix(matrix_t & a,
        matrix_t & b,
        matrix_t & c);
matrix_t vec2diag(vector_t &);

// ----------------------------------------------------------------------------
// some handy matrix functions
// bcw - 5/10/13
// divide a matrix by a scalar
bool matrixDivideScalar(matrix_t& m, double s);
// fill a matrix with a value
bool matrixFill(matrix_t& m, double f);
// subtract one matrix from another
bool matrixSubtract(matrix_t m1, matrix_t m2, matrix_t& result);
// transpose a matrix
bool matrixTranspose(matrix_t in, matrix_t& out);
// extract the matrix diagonal
bool matrixGetDiag(matrix_t m, vector_t& d);
// extract the matrix diagonal
bool matrixSetDiag(matrix_t& m, vector_t d);
// compute a covariance and correlation matrices from a file of numeric data
bool matrixComputeCovariance(matrix_t X, matrix_t& covMatrix, 
        matrix_t& corMatrix);
// read/write a matrix from/to a file with variable names header
bool matrixRead(string mFilename, matrix_t& m, vector<string>& variableNames);
bool matrixWrite(matrix_t m, string mFilename, vector<string> variableNames);

// added for modularity - bcw - 5/13/13
// get the sums of columns or rows, indicated by dim (0=columns, 1=rows), 
// into a vector
bool matrixSums(matrix_t m, vector_t& sums, int dim);
// extract a new matrix from a list of row and column indices
bool matrixExtractRowColIdx(matrix_t m, intvec_t rowIdx, intvec_t colIdx, 
        matrix_t& nm);

// added for SNPrank/CentralityRanker - bcw - 5/13/13
// get the trace of a matrix, sum of the diagonal
bool matrixGetTrace(matrix_t m, double& t);
// multiply two matrices element-by-element
bool matrixElementWiseMultiply(matrix_t m, matrix_t n, matrix_t& out);
// divide a matrix by a scalar
bool matrixMultiplyScalar(matrix_t& m, double s);
// add one matrix to another
bool matrixAdd(matrix_t m1, matrix_t m2, matrix_t& result);
// compute the degrees of the nodes in an network adjacency matrix
bool matrixComputeNodeDegrees(matrix_t a, vector_t& ad);
// threshold a matrix into a connectivity (0/1) matrix)
bool matrixConnectivityThreshold(matrix_t& m, double t, bool binary);

// added for numeric file summary stats - bcw - 5/23/13
bool reportNumericSummaryStats();
// return summary stats for a vector of values: mean, var, sd
bool vectorSummary(vector_t values, vector_t& summary);

// added for regression ranker - bcw - 5/27/13
bool rankByRegression(RegressionRankType rankType, rankedlist_t& ranks);
pair<double, double> fitModel(Model* m);

// my own matrix multiply using openmp
bool matrixMultiply(matrix_t m1, matrix_t m2, matrix_t& result);

// ----------------------------------------------------------------------------

class Eigen {
public:

  void set(int n) {
    d.resize(n, 0);
    sizeMatrix(z, n, n);
  }

  vector_t d; // eigenvalues
  matrix_t z; // eigenvectors
};

bool realnum(double);

long double factorial(int);
double normdist(double);
double ltqnorm(double);
double chi2x2(double, double, double, double);
double chi2x2(table_t);
double chi2x2(matrix_t);
double chiTable(table_t);
double chiprobP(double, double);
double symTable(table_t);
double inverse_chiprob(double, double);
double gammp(double a, double x);
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double gammln(double xx);

double rnorm();

void lubksb(vector<vector<double> > &a, vector<int> &indx, vector<double> &b);
void ludcmp(vector<vector<double> > &a, vector<int> &indx, double &d);
vector< vector<double> > inverse(vector< vector<double> > & m);

vector<double> eigenvalues(vector<vector<double> > & a);
void tred2(vector<vector<double> >&, vector<double> &, vector<double> &);
void tqli(vector<double> &d, vector<double>&e, vector<vector<double> > &z);

Eigen eigenvectors(vector<vector<double> > & a);
void EV_tred2(vector<vector<double> >&, vector<double> &, vector<double> &);
void EV_tqli(vector<double> &d, vector<double>&e, vector<vector<double> > &z);

vector< vector<double> > svd_inverse(vector< vector<double> > &, bool &);
bool svd(matrix_t &, vector_t &, matrix_t &);
bool svdcmp(vector<vector<double> > &,
        vector<double> &,
        vector<vector<double> > &);
void svbksb(vector<vector<double> > &u,
        vector<double> &w,
        vector<vector<double> > &v,
        vector<double> &b,
        vector<double> &x);
vector<vector<double> > msqrt(vector<vector<double> > & u);

double qromb(double func(const double), double a, double b);
void polint(vector_t &xa, vector_t &ya, const double x, double &y, double &dy);
double trapzd(double func(const double), const double a, const double b, const int n);

void svdvar(vector<vector<double> > & v,
        vector<double> & w,
        vector<vector<double> > & cvm);

int pca(matrix_t & x, boolmatrix_t & mask, vector_t & p, matrix_t & s, matrix_t & v, bool);

double pythag(const double a, const double b);

double betacf(const double a, const double b, const double x);
double betai(const double a, const double b, const double x);
double pF(const double F, const int df1, const int df2);
double pT(const double T, const double df);

#endif
