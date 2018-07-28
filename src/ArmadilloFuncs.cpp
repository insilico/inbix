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
#include <cmath>

#include <armadillo>

#include "plink.h"
#include "helper.h"
#include "stats.h"

#include "ArmadilloFuncs.h"
#include "Insilico.h"

using namespace arma;
using namespace std;

// differential coexpression
bool armaDcgain(sp_mat& zvals, mat& pvals, bool computeDiagonal) {
  // only works on 'numeric' data types added to PLINK by inbix
  if(PP->nl_all) {
    PP->printLOG("ERROR: armaDcgain requires 'numeric' attributes only\n");
    return false;
  }
  // phenotypes
  int nAff = 0;
  int nUnaff = 0;
  for(int i=0; i < PP->sample.size(); i++) {
    if(PP->sample[i]->aff) {
      ++nAff;
    }
    else {
      ++nUnaff;
    }
  }
  if((nAff == 0) || (nUnaff == 0)) {
    PP->printLOG("ERROR: Single phenotype detected");
    return false;
  }
  if((nAff < 4) || (nUnaff < 4)) {
    PP->printLOG("ERROR: zTest requires at least 4 individuals in a each phenotype group");
    return false;
  }
  double df = nAff + nUnaff - 2;
  PP->printLOG("Performing z-tests with " + dbl2str(df) + " degrees of freedom\n");
  PP->printLOG("WARNING: all main effect (matrix diagonal) p-values are set to 1.\n");
  uint numVars = PP->nlistname.size();
  bool readyToRun = true;
  uint infinityCount = 0;
  uint nanCount = 0;
  uint totalTests = 0;
  if(computeDiagonal) {
    PP->printLOG("Performing Z-tests for zVals and pVals matrix diagonals\n");
  } else {
    PP->printLOG("Setting matrix diagonals zVals to 0.0 and pVals to 1.0\n");
  }
  // main effects loop
  #pragma omp parallel for
  for(uint i=0; i < numVars; ++i) {
    // t-test for diagonals?
    // double t;
    // tTest(i, t);
    // double p = pT(t, df);
    // zvals(i, i) = t;
    double z = 0.0;
    double p = 1.0;
    if(computeDiagonal) {
      if(!zTest(i, z)) {
        error("Z-test failed for variable index [ " + int2str(i) + " ]");
      }
    }
    #pragma omp critical
    {
      ++totalTests;
      if(std::isnan(z)) {
        ++nanCount;
      }
      if(std::isinf(z)) {
        ++infinityCount;
      }
      zvals(i, i) = z;
      pvals(i, i) = p;
    }
  }
  // z-test for off-diagonal elements
  PP->printLOG("Computing coexpression for CASES and CONTROLS.\n");
  mat X;
  mat Y;
  if(!armaGetPlinkNumericToMatrixCaseControl(X, Y)) {
    PP->printLOG("WARNING: Cannot read numeric data into case-control matrices");
    readyToRun = false;
  }
  if(!X.is_finite()) {
    PP->printLOG("WARNING: armaGetPlinkNumericToMatrixCaseControl(X, Y) matrix X is not finite\n");
    readyToRun = false;
  }
  if(!Y.is_finite()) {
    PP->printLOG("WARNING: armaGetPlinkNumericToMatrixCaseControl(X, Y) matrix Y is not finite");
    readyToRun = false;
  }
  // cout << "X: " << X.n_rows << " x " << X.n_cols << endl;
  // cout << "Y: " << Y.n_rows << " x " << Y.n_cols << endl;
  // cout << "X" << endl << X.submat(0,0,4,4) << endl;
  // cout << "Y" << endl << Y.submat(0,0,4,4) << endl;
  // compute covariances/correlations
  mat covMatrixX;
  mat corMatrixX;
  if(readyToRun && !armaComputeCovariance(X, covMatrixX, corMatrixX)) {
    PP->printLOG("WARNING: Could not compute coexpression matrix for cases\n");
    readyToRun = false;
  }
  if(readyToRun && !covMatrixX.is_finite()) {
    PP->printLOG("WARNING: armaComputeCovariance(X, covMatrixX, corMatrixX) covMatrixX matrix is not finite");
    readyToRun = false;
  }
  if(readyToRun && !corMatrixX.is_finite()) {
    PP->printLOG("WARNING: armaComputeCovariance(X, corMatrixX, corMatrixX) corMatrixX matrix is not finite");
    readyToRun = false;
  }
  mat covMatrixY;
  mat corMatrixY;
  if(readyToRun && !armaComputeCovariance(Y, covMatrixY, corMatrixY)) {
    PP->printLOG("WARNING: Could not compute coexpression matrix for controls");
    readyToRun = false;
  }
  if(readyToRun && !covMatrixY.is_finite()) {
    PP->printLOG("WARNING: armaComputeCovariance(Y, covMatrixY, corMatrixY) covMatrixX matrix is not finite");
    readyToRun = false;
  }
  if(readyToRun && !corMatrixY.is_finite()) {
    PP->printLOG("WARNING: armaComputeCovariance(Y, covMatrixY, corMatrixY) corMatrixX matrix is not finite");
    readyToRun = false;
  }
  // are all the prerequisites satisfied?
  if(!readyToRun) {
    return readyToRun;
  }
  // DEBUG
  // cout << corMatrixX.n_rows << " x " << corMatrixX.n_cols << endl;
  // cout << "cor(X)" << endl << corMatrixX.submat(0,0,4,4) << endl;
  // cout << "cor(Y)" << endl << corMatrixY.submat(0,0,4,4) << endl;

  // algorithm ported from the R script z_test.R
  PP->printLOG("Performing Z-tests for interactions\n");
  uint i, j;
#pragma omp parallel for private(i, j) collapse(2)
  for(i=0; i < numVars; ++i) {
    for(j=0; j < numVars; ++j) {
      if(i <= j) {
        continue;
      }
      double r_ij_1 = corMatrixX(i, j);
      double r_ij_2 = corMatrixY(i, j);
      if((r_ij_1 == 1) || (r_ij_2 == 1)) {
        continue;
      }
      double z_ij_1 = 0.5 * log((abs((1 + r_ij_1) / (1 - r_ij_1))));
      double z_ij_2 = 0.5 * log((abs((1 + r_ij_2) / (1 - r_ij_2))));
      double Z_ij = abs(z_ij_1 - z_ij_2) / sqrt((1.0 / (nAff - 3.0) + 1.0 / (nUnaff - 3.0)));
      double p = 2 * normdist(-abs(Z_ij)); 
      #pragma omp critical
      {
        // happy lights
        if(i && ((i % 1000) == 0) && (j == (i + 1))) {
          PP->printLOG(Timestamp() + " - " + int2str(i) + " of " + int2str(numVars) + "\n");
        }
        // update shared memory variables
        ++totalTests;
        if(std::isinf(Z_ij)) {
          ++infinityCount;
          zvals(i, j) = 0.0;
          zvals(j, i) = 0.0;
          pvals(i, j) = 1.0;
          pvals(j, i) = 1.0;
        } else {
          if(std::isnan(Z_ij)) {
            ++nanCount;
            zvals(i, j) = 0.0;
            zvals(j, i) = 0.0;
            pvals(i, j) = 1.0;
            pvals(j, i) = 1.0;
          } else {
            zvals(i, j) = Z_ij;
            zvals(j, i) = Z_ij;
            pvals(i, j) = p;
            pvals(j, i) = p;
          }
        }
      } // OpenMP critical section
    }
  } // end OpenMP for loop
  PP->printLOG("" + int2str(totalTests) + " total number of tests\n");
  bool returnValue = true;
  if(returnValue && infinityCount) {
    PP->printLOG("ERROR(S): " + int2str(infinityCount) + " infinite Z values found\n");
    returnValue = false;
  }
  if(returnValue && nanCount) {
    PP->printLOG("ERROR(S): " + int2str(nanCount) + " nan Z values found\n");
    returnValue = false;
  }
  if(returnValue && !zvals.is_finite()) {
    PP->printLOG("ERROR(S): armaComputeCovariance(Y, covMatrixY, corMatrixY) covMatrixX matrix is not finite\n");
    returnValue = false;
  }
  if(returnValue && !pvals.is_finite()) {
    PP->printLOG("ERROR(S): armaComputeCovariance(Y, covMatrixY, corMatrixY) corMatrixX matrix is not finite\n");
    returnValue = false;
  }

  return returnValue;
}

bool armaComputeCovariance(mat X, mat& covMatrix, mat& corMatrix) {

//  bcw R implementation:
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

bool armaComputeSparseCovariance(mat X, sp_mat& covMatrix, sp_mat& corMatrix) {

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
        cerr << "Row " << (rows + 1) << ":\n" + sline + "\n";
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
  ofstream outFile(mFilename);
  if(outFile.fail()) {
    return false;
  }
  // outFile.precision(6);
  // outFile.fixed;

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

bool armaWriteSparseMatrix(sp_mat& m, string mFilename, vector<string> variableNames) {
  PP->printLOG("Writing matrix [ " + mFilename + " ]\n");
  ofstream outFile(mFilename);
  if(outFile.fail()) {
    return false;
  }
  // outFile.precision(6);
  // outFile.fixed;

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
      if(PP->sample[i]->missing) {
        error("PLINK SNP file has missing phenotype(s)");
		  } else {
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
    //printf("sample: %3d  fid: %s iid: %s\n", i, PP->sample[i]->fid.c_str(), PP->sample[i]->iid.c_str());
	}

	return true;
}
