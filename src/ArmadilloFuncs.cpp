/* 
 * File:   ArmadilloFuncs.cpp
 * Author: bwhite
 *
 * Created on June 12, 2013, 9:37 AM
 * 
 * Functions that use Armadillo Linear Algebra.
 */

#include <string>
#include <vector>

#include <armadillo>

#include "plink.h"
#include "helper.h"

#include "ArmadilloFuncs.h"

using namespace arma;
using namespace std;

extern Plink* PP;

bool armaComputeCovariance(mat X, mat& covMatrix, mat& corMatrix) {

//	N <- nrow(X)
//	p <- ncol(X)
//	# n x 1 summing vector used to create deviation score form
//	one <- matrix(rep(1, N), nrow=N)
//	P <- one %*% t(one) / N
//	Q <- diag(N) - P
//	Xstar <- Q %*% X
//	# variance-covariance
//	V <- (t(Xstar) %*% Xstar) / (N-1)
//	# compute the correlation from the covariance
//	D <- sqrt(diag(V))
//	R <- D %*% V %*% D
	
  int n = X.n_rows;

  // compute covariances
	PP->printLOG("Computing covariance matrix\n");
	vec one = ones<vec>(n);
  mat P = one * one.t() / n;
	mat diag1(n, n);
	diag1.eye();
	mat Q = diag1 - P;
  mat xStar = Q * X;
	covMatrix = xStar.t() * xStar / (n - 1);

  // compute correlations from covariances
	PP->printLOG("Computing correlation matrix\n");
	mat D = zeros<mat>(covMatrix.n_cols, covMatrix.n_cols);
	for(int i=0; i < covMatrix.n_cols; ++i) {
		D(i, i) = 1.0 / sqrt(covMatrix(i, i));
	}
	corMatrix = D * covMatrix * D;
  
  return true;
}

bool armaReadMatrix(string mFilename, mat& m, vector<string>& variableNames) {
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
      m.resize(rows, variableNames.size());
			for(int i=0; i < dataValues.size(); ++i) {
				m(rows-1, i) = dataValues[i];
			}
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

bool armaWriteMatrix(mat& m, string mFilename, vector<string> variableNames) {
  PP->printLOG("Writing matrix [ " + mFilename + " ]\n");
  ofstream outFile(mFilename.c_str());
  if(outFile.fail()) {
    return false;
  }
  outFile.precision(6);
  outFile.fixed;

  // write the variables header
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

  // write the matrix
  for(int i=0; i < m.n_rows; ++i) {
    for(int j=0; j < m.n_cols; ++j) {
      if(j) {
        outFile << "\t" << m(i, j);
      }
      else {
        outFile << m(i, j);
      }
    }
    outFile << endl;
  }
  
  outFile.close();

	return true;
}

bool armaGetPlinkNumericToMatrixAll(mat& X) {
	int numNumerics = PP->nlistname.size();
	X.resize(numNumerics, numNumerics);
	
	// load numerics into passed matrix
	for(int i=0; i < PP->sample.size(); i++) {
		for(int j=0; j < numNumerics; ++j) {
			X(i, j) = PP->sample[i]->nlist[j];
		}
	}
	
	return true;
}

bool armaGetPlinkNumericToMatrixCaseControl(mat& X, mat& Y) {
	
	// determine the number of affected and unaffected individuals
	int nAff = 0;
	int nUnaff = 0;
	for(int i=0; i < PP->sample.size(); i++) {
		if(PP->sample[i]->aff) {
			++nAff;
		}
		else {
      if(!PP->sample[i]->missing) {
		    ++nUnaff;
		  }
		}
	}
	PP->printLOG("Detected " + int2str(nAff) + " affected and " + 
					int2str(nUnaff) + " unaffected individuals\n");
	// size matrices
	int numNumerics = PP->nlistname.size();
	X.resize(nAff, numNumerics);
	Y.resize(nUnaff, numNumerics);
	
	// load numerics into passed matrices
	PP->printLOG("Loading case and control matrices\n");
	int aIdx = 0;
	int uIdx = 0;
	for(int i=0; i < PP->sample.size(); i++) {
		for(int j=0; j < numNumerics; ++j) {
			if(PP->sample[i]->aff) {
				X(aIdx, j) = PP->sample[i]->nlist[j];
			} else {
				Y(uIdx, j) = PP->sample[i]->nlist[j];
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
