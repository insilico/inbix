/*
 * CentralityRanker.cpp - Bill White - 5/15/13
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "plink.h"
#include "stats.h"
#include "helper.h"

#include "StringUtils.h"
#include "CentralityRanker.h"

using namespace std;
using namespace Insilico;

// Plink object
extern Plink* PP;

CentralityRanker::CentralityRanker(string gainFileParam, bool isUpperTriangular)
{

#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: CentralityRanker constructor, GAIN file: "
			<< gainFileParam << endl;
#endif

	// read the file and save the header (markers) and GAIN matrix
	if(!ReadGainFile(gainFileParam, isUpperTriangular)) {
		error("Reading GAIN file: " + gainFileParam + "\n");
	}
	gainFile = gainFileParam;
	gamma = 0.0;
	n = G[0].size();
  sizeMatrix(Gdiag, n, n);
}

CentralityRanker::CentralityRanker(double** variablesMatrix, unsigned int dim,
		vector<string>& variableNames)
{
	// setup G matrix for  algorithm
	sizeMatrix(G, dim, dim);
	// copy variable matrix values into G
	for(unsigned int i=0; i < dim; ++i) {
		for(unsigned int j=0; j < dim; ++j) {
			G[i][j] = G[j][i] = variablesMatrix[i][j];
		}
	}
	for(unsigned int i=0; i < dim; ++i) {
		variableNames.push_back(variableNames[i]);
	}
	// set default values
	gainFile = "";
	gamma = 0.0;
	n = dim;
  sizeMatrix(Gdiag, n, n);
}

CentralityRanker::CentralityRanker(matrix_t& A, vector<string>& variableNames)
{
	unsigned int dim = A[0].size();
	// setup G matrix for  algorithm
	sizeMatrix(G, dim, dim);
	// copy variable matrix values into G
	for(unsigned int i=0; i < dim; ++i) {
		for(unsigned int j=0; j < dim; ++j) {
			G[i][j] = G[j][i] = A[i][j];
		}
	}
	for(unsigned int i=0; i < dim; ++i) {
		variableNames.push_back(variableNames[i]);
	}
	// set default values
	gainFile = "";
	gamma = 0.0;
	n = G[0].size();
  sizeMatrix(Gdiag, n, n);
}

CentralityRanker::~CentralityRanker()
{}

bool CentralityRanker::CalculateCentrality(SolverMethod method)
{
	// NOTE: This code attempts to match the SnpRankNextGen code, which is
	// based on the original Matlab .m code.

	// G diagonal is main effects
  vector_t diag;
  matrixGetDiag(G, diag);
  matrixSetDiag(Gdiag, diag);
#ifdef DEBUG_CENTRALITY
	display(diag);
#endif

	// sum of diagonal elements is Gtrace
	matrixGetTrace(G, Gtrace);
#ifdef DEBUG_CENTRALITY
	cout << "Gtrace: " << Gtrace << endl;;
#endif

#ifdef DEBUG_CENTRALITY
	cout << "n: " << n << endl;
#endif

	// colsum = degree of each variable(in undirected graphs, out-degree = in-degree)
	// 1 x n (row) vector of column sums (d_j in Eqs. 4, 5 from SNPRank paper)
	matrixSums(G, colsum, 0);
#ifdef DEBUG_CENTRALITY
	display(colsum);
#endif

	if(method == POWER_METHOD) {
		PowerMethodSolver();
	}
	else {
		GaussEliminationSolver();
	}

	return true;
}

void CentralityRanker::SetGlobalGamma(double gammaParam)
{
	gamma = gammaParam;
}

bool CentralityRanker::SetGammaVector(vector_t& gammaVectorValues)
{
	if(gammaVectorValues.size() != n) {
		cerr << "CentralityRanker::SetGammaVector vector size mismatch: "
				<< gammaVectorValues.size() << " vs. " << n << endl;
		return false;
	}
	gammaVector.clear();
	for(vector_t::const_iterator it=gammaVectorValues.begin();
			it != gammaVectorValues.end(); ++it) {
		gammaVector.push_back(*it);
	}
	return true;
}

void CentralityRanker::WriteToFile(string outFile)
{
  PP->printLOG("Writing centrality scores to [" + outFile + "]\n");
	// r indices sorted in descending order
	sort(ranks.begin(), ranks.end());
  reverse(ranks.begin(), ranks.end());
	// output r (centrality rankings) to file, truncating to 6 decimal places
	ofstream outputFileHandle(outFile.c_str());
	outputFileHandle << "SNP\tCentrality\tdiag\tdegree" << endl;
	for(int i = 0; i < r.size(); i++) {
    double score = ranks[i].first;
    int index = ranks[i].second;
    outputFileHandle << variableNames[index] << "\t"
         << score << "\t"
         << Gdiag[index][index] << "\t"
         << colsum[index];
    outputFileHandle << endl;
	}
	outputFileHandle.close();
}

void CentralityRanker::WriteToConsole()
{
	// r indices sorted in descending order
	sort(ranks.begin(), ranks.end());
  reverse(ranks.begin(), ranks.end());
	// output r (centrality rankings) to file, truncating to 6 decimal places
	cout << "SNP\tCentrality\tdiag\tdegree" << endl;
	for(int i = 0; i < r.size(); i++) {
    double score = ranks[i].first;
    int index = ranks[i].second;
    cout << variableNames[index] << "\t"
         << score << "\t"
         << Gdiag[index][index] << "\t"
         << colsum[index];
    cout << endl;
	}
}

// ------------------ P R I V A T E   M E T H O D S --------------------------

bool CentralityRanker::ReadGainFile(string gainFilename, bool isUpperTriangular)
{
	ifstream gainFileHandle(gainFilename.c_str());
	if(!gainFileHandle.is_open()) {
		cerr << "ERROR: Could not open (re)GAIN file: " << gainFilename << endl;
		return false;
	}

	string line;

	// read first line (header)
	getline(gainFileHandle, line);
	string trimmedLine = trim(line);

	// tokenize header
	split(variableNames, trimmedLine);

	// initialize dimensions of data matrix
	size_t numVariables = variableNames.size();
	if(!numVariables) {
		cerr << "ERROR: Could not parse SNP names from (re)GAIN file header" << endl;
		return false;
	}
	sizeMatrix(G, numVariables, numVariables);

	// read numeric data into G
	size_t row = 0;
	vector<string> lineTokens;
	unsigned int tokensExpected = numVariables;
	while(getline(gainFileHandle, line)) {
		trim(line);
		lineTokens.clear();
		split(lineTokens, line);
		// cout << "Read " << lineTokens.size() << " expected: " << tokensExpected << endl;
		if(lineTokens.size() != tokensExpected) {
			cerr << "ERROR line:" << endl << endl << line << endl << endl;
			cerr << "ERROR: Could not parse (re)GAIN file row: " << (row+2) << endl
					<< "Expecting " << tokensExpected << " values, got "
					<< lineTokens.size() << endl;
			return false;
		}
		size_t startIndex = numVariables - tokensExpected;
		for (size_t tokenIndex=0, col=startIndex; col < numVariables; ++tokenIndex, ++col) {
			string token = lineTokens[tokenIndex];
      trim(token);
			if(token == "") {
				cerr << "ERROR: parsing line " << row << " col " << col
						<< " of (re)Gain file, token [" << token << "]" << endl;
				return false;
			}
			// cout << "row: " << row << ", col: " << col << endl;
			double t;
			if(!from_string<double>(t, token, std::dec)) {
				error("Parsing REGAIN line " + line);
			}
			G[row][col] = t;
			if((row != col) && isUpperTriangular) {
				if(!from_string<double>(t, token, std::dec)) {
					error("Parsing REGAIN line " + line);
				}
				G[col][row] = t;
			}
		}
		row++;
		if(isUpperTriangular) {
			--tokensExpected;
		}
	}

	gainFileHandle.close();

	return true;
}

bool CentralityRanker::GaussEliminationSolver()
{
#ifdef DEBUG_CENTRALITY
	cout << "GAUSS METHOD" << endl;
#endif

	// zero the G diagonal
  vector_t diagZ(n, 0);
	matrixSetDiag(G, diagZ);

#ifdef DEBUG_CENTRALITY
	display(G);
#endif

	// by bam in Matlba/Octave - added here by bcw - 10/26/12
	// including the factor of n (below) reduces the effect of main effects.
	// if you remove the n, the centrality is highly correlated with main effect
	// too correlated I think. I think the n creates a good balance between
	// correlation of centrality with main effect and degree.
	Gtrace *= n;

	// create gamma vector/matrix
	vector_t colsum_G;
  matrixSums(G, colsum_G, 0);
#ifdef DEBUG_CENTRALITY
	display(colsum_G);
#endif
	vector_t rowsum_G;
  matrixSums(G, rowsum_G, 1);
#ifdef DEBUG_CENTRALITY
	display(rowsum_G);
#endif

	vector_t rowsum_denom(n);
	for(unsigned int i=0; i < n; ++i) {
		double localSum = 0;
		for(unsigned int j=0; j < n; ++j) {
			// added this check per snprank_nextgen_octave_3b.m
			// from bam 12/13/12
			double factor = G[i][j]? 1: 0;
			 localSum += (factor * colsum_G[j]);
		}
		rowsum_denom[i] = localSum;
	}
#ifdef DEBUG_CENTRALITY
	display(rowsum_denom);
#endif

	vector_t gamma_vec(n);
	if(gammaVector.size()) {
		for(unsigned int i=0; i < n; ++i) {
			gamma_vec[i] = gammaVector[i];
		}
	}
	else {
		for(unsigned int i=0; i < n; ++i) {
			if(gamma != 0) {
				gamma_vec[i] = gamma;
			}
			else {
				if(rowsum_denom[i] != 0) {
					gamma_vec[i] = rowsum_G[i] / rowsum_denom[i];
				}
				else {
					gamma_vec[i] = 0;
				}
			}
		}
	}
	if(gamma == 0) {
		double averageGamma = 0;
    for(int i=0; i < gamma_vec.size(); ++i) {
      averageGamma += gamma_vec[i];
    }
    averageGamma /= (double) n;
		cout << "Average gamma = " << averageGamma << endl;
	}
  
  // repmat(gamma_vec, n, 1);
	matrix_t gamma_matrix;
  sizeMatrix(gamma_matrix, n, n);
  for(int i=0; i < n; ++i) {
    for(int j=0; j < n; ++j) {
      gamma_matrix[i][j] = gamma_vec[j];
    }
  }
  
#ifdef DEBUG_CENTRALITY
	display(gamma_vec);
	display(gamma_matrix);
#endif

	// b is the "charity vector"
	vector_t b(n);
	for(unsigned int i=0; i < n; ++i) {
		if(Gtrace != 0) {
			b[i] = ((1.0 - gamma_vec[i]) / n) + (Gdiag[i][i] / Gtrace);
		}
		else {
			b[i] = (1.0 - gamma_vec[i]) / n;
		}
	}
#ifdef DEBUG_CENTRALITY
	display(b);
#endif

	// D = n x n sparse matrix with 1/colsum for nonzero elements of colsum
  matrix_t D;
  sizeMatrix(D, n, n);
	for (size_t i = 0; i < n; i++){
		if(colsum_G[i] != 0) {
			D[i][i] = 1 / colsum_G[i];
		}
		else {
			D[i][i] = 0;
		}
	}
#ifdef DEBUG_CENTRALITY
	display(D);
#endif

	// SOLVE: Ax = b system of equations using gaussian elimination
	// A = (I - gamma_matrix % G * D)
	// b = b, x = r = centrality scores
	// r = solve((I - gamma_matrix % G * D), b);

  // -------------- replaces Armadillo one-liner above ! ----------------------
	matrix_t I;
  sizeMatrix(I, n, n); 
  vector_t eye(n, 1);
  matrixSetDiag(I, eye);
  matrix_t tmpMatrix;
  matrixElementWiseMultiply(gamma_matrix, G, tmpMatrix);
  matrix_t tmpMatrix2;
  multMatrix(tmpMatrix, D, tmpMatrix2);
  matrix_t u;
  sizeMatrix(u, n, n);
  matrixSubtract(I, tmpMatrix2, u);
  // solve Ar=b
  vector_t w(n);
  matrix_t v;
	sizeMatrix(v, n, n);
	const double TOL = 1.0e-13;
#ifdef DEBUG_CENTRALITY_SOLVER
  display(u);
  display(v);
  display(w);
#endif
	if(!svdcmp(u, w, v)) {
    display(u);
    display(w);
    display(v);
    error("SVD solver failed");
	}
#ifdef DEBUG_CENTRALITY_SOLVER
  display(u);
  display(v);
  display(w);
  cout << "Performing back substitution" << endl;
#endif

	double wmax = 0.0;
	for(int j = 0; j < n; j++) {
		if(w[j] > wmax) {
			wmax = w[j];
		}
	}
	double thresh = TOL * wmax;
	for(int j = 0; j < n; j++) {
		if(w[j] < thresh) {
			w[j] = 0.0;
		}
	}
  r.resize(n);
	svbksb(u, w, v, b, r);
#ifdef DEBUG_CENTRALITY_SOLVER
	display(r);
#endif
  // --------------------------------------------------------------------------

  //	r = r / sum(r);
  double rSum = 0;
  for(int i=0; i < r.size(); ++i) {
    rSum += r[i];
  }
  for(int i=0; i < r.size(); ++i) {
    r[i] /= rSum;
    ranks.push_back(make_pair(r[i], i));
  }
#ifdef DEBUG_CENTRALITY_SOLVER
	display(r);
#endif

	return true;
}

bool CentralityRanker::PowerMethodSolver()
{
	// compute initial T, Markov chain transition matrix
	// first term (gamma * G * D) is zero when d_j = 0
	// for right hand term of T, Gdiag is nx1 and T_nz is 1xn
	// initialize second term multiplier as a row vector of 1s
	vector_t diag_T_nz(n, 1);
  for(int i=0; i < n; ++i) {
    diag_T_nz[i] /= n;
  }

	// (A) if row degree is 0, this will give 1/n gives more charity
  // if no incoming connections
	// (B) if row degree is non-0, gives (1-gamma)/n
	// indices of colsum vector that are non-zero
	intvec_t colsum_nzidx;
  for(int i=0; i < colsum.size(); ++i) {
    double thisSum = colsum[i];
    if(thisSum) {
      colsum_nzidx.push_back(i);
    }
  }

	// D is an n x n sparse matrix with values 1/colsum for nonzero elements
	// of colsum indices given by colsum_nzidx.  Other elements are zero.
  sizeMatrix(D, n, n);
  matrixFill(D, 0);
	// D = zeros<matrix_t>(n, n);
	// matrix_t fillD = ones(colsum_nzidx.n_elem, 1) / colsum.elem(colsum_nzidx);
  //	for (int i = 0; i < (int) fillD.size(); i++){
  //		D(colsum_nzidx[i], colsum_nzidx[i]) = fillD[i];
  //	}
  for(int i=0; i < colsum_nzidx.size(); ++i) {
    int colIndex = colsum_nzidx[i];
    D[i][i] = 1 / colsum[i];
  }
#ifdef DEBUG_CENTRALITY
	display(D);
#endif

	for (int i = 0; i < colsum_nzidx.size(); i++){
		diag_T_nz[colsum_nzidx[i]] = (1 - gamma) / n;
	}

	// find the diagonal elements that are non-zero
	intvec_t diag_nzidx;
  for(int i=0; i < n; ++i) {
    if(Gdiag[i][i]) {
      diag_nzidx.push_back((1 - gamma) * Gdiag[i][i] / n);
    }
  }
	// (C) for non-zero diag elements, gives (1-gamma)gii/n,
	// which supercedes (A) and (B)

	// this !replaces! the constant charity term if there are non-zero main effects
  //	for (int i = 0; i < (int) diag_nzidx.n_elem; i++){
  //		diag_T_nz(diag_nzidx[i]) = (1 - gamma) * Gdiag(diag_nzidx[i]) / n;
  //	}

	vector_t e(n, 1);
	matrix_t T;
  sizeMatrix(T, n, n);
  
  matrix_t GD;
  multMatrix(G, D, GD);
  matrixMultiplyScalar(GD, gamma);
  matrix_t ed(n);
  for(int i=0; i < n; ++i) {
    for(int j=0; j < n; ++j) {
      ed[i][j] = e[i] * diag_T_nz[j];
    }
  }
	if(Gtrace == 0.0) {
    matrixAdd(GD, ed, T);
	}
	else {
    matrixAdd(GD, ed, T);
		matrixDivideScalar(T, Gtrace);
	}

	// initialize size of vector r to store centrality scores
  //	r.set_size(G.size());
  //	r.fill(1.0 / G.size());
  r.resize(n, 1 / n);
	double threshold = 1.0E-4;
	double lambda = 0.0;
	bool converged = false;
	vector_t r_old = r;

	// if the absolute value of the difference between old and current r
	// vector is < threshold, we have converged
	while(!converged) {
		r_old = r;
		// r = T * r;
    for(int i=0; i < n; ++i) {
      double tempSum = 0;
      for(int j=0; j < n; ++j) {
        tempSum += T[i][j] * r[j];
      }
      r[i] = tempSum;
    }

		// sum of r elements
		lambda = 0;
    for(int i=0; i < r.size(); ++i) {
      lambda += r[i];
    }

		// normalize eigenvector r so sum(r) == 1
		// r = r / lambda;
    for(int i=0; i < r.size(); ++i) {
      r[i] /= lambda;;
    }

		// check convergence, ensure all elements of r - r_old < threshold
    //		if(min((intvec_t) (abs(r - r_old) < threshold)) == 1) {
    //			converged = true;
    //		}
    //		else {
    //			converged = false;
    //		}
    for(int i=0; i < r.size(); ++i) {
      if(fabs(r[i] - r_old[i]) > threshold) {
        converged = false;
        continue;
      }
    }

	}

	return converged;
}
