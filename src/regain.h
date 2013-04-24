/*==============================================================================
 *
 * Filename:  regain.h - Bill White - 4/23/13
 *
 * Description:  Regression GAIN calculations.  Uses linear or logistic 
 * regression in calculating main effects and interactions of SNPs and 
 * numeric attributes.
 * 
 * Originaly Created:  02/02/2012
 * Author:  Nick Davis, nick-davis@utulsa.edu
 * =============================================================================
 */

#ifndef __REGAIN_H__
#define __REGAIN_H__

#include <iostream>
#include <vector>

#include "model.h"
#include "zfstream.h"

using namespace std;

// type for storing p-value and matrix position (row,col) of
// reGAIN interaction terms
typedef pair< double, pair<int, int> > mat_el;

class Regain {
public:
	Regain(bool compr, double sifthr, bool integrative, bool compo,
					bool fdrpr = false);
	~Regain();
	// iterate over all SNPs and numeric attributes (if present), calculating
	// regression for main effect and interaction terms
	void run();
	// calculate main effect regression coefficients for diagonal terms	
	void mainEffect(int e1, bool numeric);
	// add covariate terms to model, handling multiple covariates
	void addCovariates(Model &m);
	// calculate epistatic interaction between two SNPs or numeric attributes,
	// or a SNP and a numeric attribute
	void interactionEffect(int e1, bool numeric1, int e2, bool numeric2);
	// write contents of reGAIN or reGAIN p-value matrix to file.  Integrative
	// reGAIN files have a '.int.' in the filename, p-values have a '.pvals.'
	// If fdr is true, the filename is .pruned.regain, and the log output
	// reflects this.
	void writeRegain(bool pvalues, bool fdrprune = false);
	// Benjamini Hochberg FDR pruning - removes interaction 
	// terms from reGAIN matrix based on BH FDR threshold
	// code based on method described in All of Statistics p. 167
	void fdrPrune(double fdr);
	// write R commands to plot FDR graph
	void writeRcomm(double T, double fdr);
	// comparison fnc for mat_el types
	static bool mecomp(const mat_el &l, const mat_el &r);
private:
	// integrative regain mode
	bool intregain;
	// use zlib compression?
	bool compressed;
	// apply FDR pruning matrix?
	bool fdrprune;
	// write out component matrices
	bool component;
	// num attributes (SNPs + numeric for integrative, SNPs for normal regain)
	int numattr;
	// SIF interaction threshold
	double sif_thresh;
	// Output matrix files (used for writing regain and p-values files)
	ZOutput REGAIN_MATRIX;
	ZOutput SNP_MATRIX;
	ZOutput NUM_MATRIX;
	ZOutput INT_MATRIX;
	// additional output files
	ofstream MEBETAS;
	ofstream BETAS;
	ofstream SIF;
	ofstream SNP_SIF;
	ofstream NUM_SIF;
	ofstream INT_SIF;
	// in memory arrays
	double** regainMatrix;
	double** regainPMatrix;
	// collection of all interaction terms as mat_el types
	vector<mat_el> gainPint;
};
#endif
