/*
 * CentralityRanker.h - Bill White - 5/15/13
 */

#ifndef CENTRALITY_RANKER_H
#define CENTRALITY_RANKER_H

#include <string>
#include <vector>
#include <cstring>

#include <armadillo>

#include "plink.h"

enum SolverMethod { POWER_METHOD, GAUSS_ELIMINATION };

class CentralityRanker {
public:
	// construct using a file representing the variable interactions matrix
	CentralityRanker(std::string gainFileParam, bool isUpperTriangular=false);
	// matrix constructor for calling as a library method
	CentralityRanker(double** variablesMatrix, unsigned int dim,
			std::vector<std::string>& variableNames);
	CentralityRanker(arma::mat& A, std::vector<std::string> varNames);
	CentralityRanker(arma::sp_mat& A, std::vector<std::string> varNames);
	virtual ~CentralityRanker();

	bool CalculateCentrality(SolverMethod method);
	void SetGlobalGamma(double gammaParam);
	bool SetGammaVector(vector_t& gammaVectorValues);
	bool Permute(int numPerms);
	// write results
	void WriteToFile(std::string outfile, int topN = -1);
	void WriteToConsole(int topN = -1);
	arma::vec GetResultsByVariable() { return r; }
private:
	// reads the GAIN input file
	bool ReadGainFile(std::string gainFilename, bool isUpperTriangular=false);

	// use system of equations solver (Gaussian elimination)
	bool GaussEliminationSolver();
	// use power method solver
	bool PowerMethodSolver();

	// data file
	std::string gainFile;

	// gamma parameter
	double gamma;
	std::vector<double> gammaVector;

	// vector that stores variable names
	std::vector<std::string> variableNames;
	// GAIN matrix
	arma::mat G;
	// intermediate results of centrality calculations
	arma::rowvec colsum;
	size_t n;
	arma::mat Gdiag;
	double Gtrace;
	arma::mat D;

	// centrality rank scores
	arma::vec r;
  vector<pair<double, int> > ranks;
  
  int topN;
};

#endif /* CENTRALITY_RANKER_H */
