/* 
 * File:   ArmadilloFuncs.h
 * Author: bwhite
 *
 * Created on June 12, 2013, 9:37 AM
 */

#ifndef ARMADILLOFUNCS_H
#define	ARMADILLOFUNCS_H

#include <string>
#include <vector>
#include <armadillo>

// compute a covariance and correlation matrices from a file of numeric data
bool armaComputeCovariance(arma::mat X, arma::mat& covMatrix, arma::mat& corMatrix);
bool armaComputeSparseCovariance(arma::mat X, arma::sp_mat& covMatrix, arma::sp_mat& corMatrix);

// read an Armadillo matrix from a tab-delimited text file
bool armaReadMatrix(std::string mFilename, arma::mat& m, 
				std::vector<std::string>& variableNames);

// write an Armadillo matrix to a tab-delimited text file
bool armaWriteMatrix(arma::mat& m, std::string mFilename, 
				std::vector<std::string> variableNames);
bool armaWriteSparseMatrix(arma::sp_mat& m, std::string mFilename, 
				std::vector<std::string> variableNames);

// get PLINK numeric data to an Armadillo matrix for all individuals
bool armaGetPlinkNumericToMatrixAll(arma::mat& X);
bool armaGetPlinkNumericToMatrixCaseControl(arma::mat& X, arma::mat& Y);

bool armaDcgain(arma::sp_mat& results, arma::mat& pvals, bool computeDiagonal=false);

#endif	/* ARMADILLOFUNCS_H */
