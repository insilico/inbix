/*
 * CentralityRanker.cpp - Bill White - 5/15/13
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <ios>
#include <float.h>

#include <armadillo>

#include "plink.h"
#include "stats.h"
#include "helper.h"

#include "StringUtils.h"
#include "CentralityRanker.h"

using namespace std;
using namespace Insilico;
using namespace arma;

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
	n = G.n_cols;
  Gdiag.resize(n, n);
}

CentralityRanker::CentralityRanker(double** variablesMatrix, unsigned int dim,
		vector<string>& variableNames)
{
	// setup G matrix for  algorithm
	G.resize(dim, dim);
	// copy variable matrix values into G
	for(unsigned int i=0; i < dim; ++i) {
		for(unsigned int j=0; j < dim; ++j) {
			G(i, j) = G(j, i) = variablesMatrix[i][j];
		}
	}
	for(unsigned int i=0; i < dim; ++i) {
		variableNames.push_back(variableNames[i]);
	}
	// set default values
	gainFile = "";
	gamma = 0.0;
	n = dim;
  Gdiag.resize(n, n);
}

CentralityRanker::CentralityRanker(mat& A, vector<string>& variableNames)
{
	unsigned int dim = A.n_cols;
	// setup G matrix for SNPrank algorithm
	G.resize(dim, dim);
	// copy variable matrix values into G
	for(unsigned int i=0; i < dim; ++i) {
		for(unsigned int j=0; j < dim; ++j) {
			G(i, j) = G(j, i) = A(i, j);
		}
	}
	for(unsigned int i=0; i < dim; ++i) {
		variableNames.push_back(variableNames[i]);
	}
	// set default values
	gainFile = "";
	gamma = 0.0;
	n = G.n_cols;
  Gdiag.resize(n, n);
}

CentralityRanker::~CentralityRanker()
{}

bool CentralityRanker::CalculateCentrality(SolverMethod method)
{
	// NOTE: This code attempts to match the original Matlab .m code.

	// G diagonal is main effects
	Gdiag = G.diag();
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: Gdiag: " << endl << Gdiag << endl;
#endif

	//sum of diagonal elements is Gtrace
	Gtrace = trace(G);
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: Gtrace: "	<< Gtrace << endl;
#endif

	//dimension of data matrix
	n = G.n_rows;
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: n: " << n << endl;
#endif

	// colsum = degree of each SNP(in undirected graphs, out-degree = in-degree)
	// 1 x n (row) vector of column sums (d_j in Eqs. 4, 5 from SNPRank paper)
	colsum = sum(G);
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: colsum: " << endl << colsum << endl;
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

void CentralityRanker::WriteToFile(string outFile, int topN)
{
  PP->printLOG("Writing centrality scores to [" + outFile + "]\n");
	// r indices sorted in descending order
	uvec r_indices = sort_index(r, 1);
	int index = 0;

	// output r (SNPrank rankings) to file, truncating to 6 decimal places
	ofstream outputFileHandle(outFile.c_str());
	streamsize savedPrecision = cout.precision();
	cout.precision(numeric_limits<double>::digits10);
	outputFileHandle << "SNP\tSNPrank\tdiag\tdegree" << endl;
	for(int i = 0; i < (int)r.n_elem; i++) {
			index = r_indices[i];
			outputFileHandle << variableNames[index] << "\t"
					 << r[index] << "\t"
					 << Gdiag(index) << "\t"
					 << colsum(index);
			outputFileHandle << endl;
	}
	outputFileHandle.close();
	cout.precision(savedPrecision);
}

void CentralityRanker::WriteToConsole(int topN)
{
	// r indices sorted in descending order
	uvec r_indices = sort_index(r, 1);
	int index = 0;
	// output r (SNPrank rankings) to file with full precision
	streamsize savedPrecision = cout.precision();
	cout.precision(numeric_limits<double>::digits10);
	cout << "SNP\tSNPrank\tdiag\tdegree" << endl;
	for(int i = 0; i < (int) r.n_elem; i++) {
			index = r_indices[i];
			cout << variableNames[index] << "\t"
					 << r[index] << "\t"
					 << Gdiag(index) << "\t"
					 << colsum(index);
			cout << endl;
	}
	cout.precision(savedPrecision);
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
	G.set_size(numVariables, numVariables);

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
			G(row, col) = t;
			if((row != col) && isUpperTriangular) {
				if(!from_string<double>(t, token, std::dec)) {
					error("Parsing REGAIN line " + line);
				}
				G(col, row) = t;
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
	G.diag() = zeros<vec>(n);
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: new G = G - Gdiag: " << endl << G << endl;
#endif

	// by bam in Matlba/Octave - added here by bcw - 10/26/12
	// including the factor of n (below) reduces the effect of main effects.
	// if you remove the n, the snprank is highly correlated with main effect
	// too correlated I think. I think the n creates a good balance between
	// correlation of snprank with main effect and degree.
	Gtrace *= n;

	// create gamma vector/matrix
	rowvec colsum_G = sum(G);
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: colsum: " << endl << colsum_G << endl;
#endif
	rowvec rowsum_G = sum(trans(G));
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: rowsum: " << endl << rowsum_G << endl;
#endif

	colvec rowsum_denom(n);
	for(unsigned int i=0; i < n; ++i) {
		double localSum = 0;
		for(unsigned int j=0; j < n; ++j) {
			// added this check per snprank_nextgen_octave_3b.m
			// from bam 12/13/12
			double factor = G(i, j)? 1: 0;
			 localSum += (factor * colsum_G(j));
		}
		rowsum_denom(i) = localSum;
	}
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: rowsum_denom: " << endl << rowsum_denom << endl;
#endif

	rowvec gamma_vec(n);
	if(gammaVector.size()) {
		for(unsigned int i=0; i < n; ++i) {
			gamma_vec(i) = gammaVector[i];
		}
	}
	else {
		for(unsigned int i=0; i < n; ++i) {
			if(gamma != 0) {
				gamma_vec(i) = gamma;
			}
			else {
				if(rowsum_denom(i) != 0) {
					gamma_vec(i) = rowsum_G(i) / rowsum_denom(i);
				}
				else {
					gamma_vec(i) = 0;
				}
			}
		}
	}
	if(gamma == 0) {
		double averageGamma = sum(gamma_vec) / (double) n;
		cout << "Average gamma = " << averageGamma << endl;
	}
	mat gamma_matrix = repmat(gamma_vec, n, 1);
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: gamma_vec: " << endl << gamma_vec << endl;
	cout << "DEBUG: gamma_matrix: " << endl << gamma_matrix << endl;
#endif

	// b is the "charity vector"
	colvec b(n);
	for(unsigned int i=0; i < n; ++i) {
		if(Gtrace != 0) {
			b(i) = ((1.0 - gamma_vec(i)) / n) + (Gdiag(i) / Gtrace);
		}
		else {
			b(i) = (1.0 - gamma_vec(i)) / n;
		}
	}
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: b: " << endl << b << endl;
#endif

	// D = n x n sparse matrix with 1/colsum for nonzero elements of colsum
	D = zeros<mat>(n, n);
	for (size_t i = 0; i < n; i++){
		if(colsum_G(i) != 0) {
			D(i, i) = 1 / colsum_G(i);
		}
		else {
			D(i, i) = 0;
		}
	}
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: D: " << endl << D << endl;
#endif

	// SOLVE: Ax = b system of equations using gaussian elimination
	// A = (I - gamma_matrix % G * D)
	// b = b, x = r = snprank scores
	mat I(n, n); I.eye();
  mat temp = I - gamma_matrix % G * D;
  // cout << temp << endl;
	r = solve(temp, b);
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: r: " << endl << r << endl;
#endif
	r = r / sum(r);
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: r (normalized): " << endl << r << endl;
#endif

	return true;
}

bool CentralityRanker::PowerMethodSolver()
{
	// compute initial T, Markov chain transition matrix
	// first term (gamma * G * D) is zero when d_j = 0
	// for right hand term of T, Gdiag is nx1 and T_nz is 1xn
	// initialize second term multiplier as a row vector of 1s

	rowvec diag_T_nz;
	diag_T_nz = ones(1, n) / n;
	// (A) if row degree is 0, this will give 1/n gives more charity
  // if no incoming connections

	// (B) if row degree is non-0, gives (1-gamma)/n
	// indices of colsum vector that are non-zero
	uvec colsum_nzidx = find(colsum);
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: colsum_nzidx: " << endl << colsum_nzidx << endl;
#endif

	// D is an n x n sparse matrix with values 1/colsum for nonzero elements
	// of colsum indices given by colsum_nzidx.  Other elements are zero.
	D = zeros<mat>(n, n);
	mat fillD = ones(colsum_nzidx.n_elem, 1) / colsum.elem(colsum_nzidx);
	for (int i = 0; i < (int) fillD.size(); i++){
		D(colsum_nzidx[i], colsum_nzidx[i]) = fillD[i];
	}
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: D: " << endl << D << endl;
  cout << "DEBUG: Gdiag: " << endl << Gdiag << endl;
	cout << "DEBUG: G: " << endl << G << endl;
#endif

	for (int i = 0; i < (int) colsum_nzidx.n_elem; i++){
		diag_T_nz(colsum_nzidx[i]) = (1 - gamma) / n;
	}

	// find the diagonal elements that are zero
	uvec diag_nzidx = find(Gdiag);
	// (C) for non-zero diag elements, gives (1-gamma)gii/n,
	// which supercedes (A) and (B)

	// this !replaces! the constant charity term if there are non-zero main effects
	for (int i = 0; i < (int) diag_nzidx.n_elem; i++){
		diag_T_nz(diag_nzidx[i]) = (1 - gamma) * Gdiag(diag_nzidx[i]) / n;
	}
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: diag_T_nz: " << endl << diag_T_nz << endl;
	cout << "DEBUG: diag_nzidx: " << endl << diag_nzidx << endl;
#endif

	arma::colvec e = ones(n, 1);
	mat T = zeros<mat>(n, n);
  mat GD = G * D;
  mat ed = (e * diag_T_nz);
  mat gammaGD = gamma * GD;
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: GD: " << endl << GD << endl;
	cout << "DEBUG: gammaGD: " << endl << gammaGD << endl;
	cout << "DEBUG: ed: " << endl << ed << endl;
#endif
  
	if(Gtrace == 0.0) {
		T = (gamma * G * D) + (e * diag_T_nz);
	}
	else {
		T = (gamma * G * D) + (e * diag_T_nz) / Gtrace;
	}
#ifdef DEBUG_CENTRALITY
	cout << "DEBUG: T: " << endl << T << endl;
#endif

	// initialize size of vector r to store snprank scores
	r.set_size(n);
	r.fill(1.0 / n);

	double threshold = 1.0E-4;
	double lambda = 0.0;
	bool converged = false;
	vec r_old = r;

	// if the absolute value of the difference between old and current r
	// vector is < threshold, we have converged
  int iterations = 0;
  cout << "Entering convergence loop" << endl;
	while(!converged) {
    ++iterations;
    
		r_old = r;
		r = T * r;
    
		// sum of r elements
		lambda = sum(r);

		// normalize eigenvector r so sum(r) == 1
		r = r / lambda;

		// check convergence, ensure all elements of r - r_old < threshold
		if(min((uvec) (abs(r - r_old) < threshold)) == 1) {
			converged = true;
		}
		else {
			converged = false;
		}
	}

	return converged;
}
