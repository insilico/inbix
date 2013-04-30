//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>

#include "plink.h"

using namespace std;

class Model {
public:
	Model();

	virtual ~Model() {
	};

	virtual void setDependent() = 0;
	virtual void pruneY() = 0;
	virtual void fitLM() = 0;
	virtual vector_t getCoefs() = 0;
	virtual vector_t getVar() = 0;
	virtual vector_t getSE() = 0;
	virtual vector_t getPVals() = 0;
	virtual void displayResults(ofstream &, Locus *) = 0;
	virtual void fitUnivariateLM() = 0;

	void setMissing();
	vector<bool> getMissing();
	void setMissing(vector<bool>&);
	void yokeMissing(Model *);
	void setHaploid();
	void setX();
	void setDominant();
	void setRecessive();
	void hasSNPs(bool);
	void addAdditiveSNP(int);
	void addDominanceSNP(int);
	void addHaplotypeDosage(set<int>&);
	void addSexEffect();
	bool isSexInModel();
	void addCovariate(int);
	void addNumeric(int);
	void addInteraction(int, int);
	void buildDesignMatrix();
	bool checkVIF();
	vector<bool> validParameters();

	bool isValid() {
		return all_valid;
	}
	double getStatistic();
	//  double getPValue();
	double linearHypothesis(matrix_t &, vector_t &);

	int Ysize() {
		return nind;
	}

	int getNP() {
		return np;
	}

	void setValid() {
		all_valid = true;
	}

	void noCluster();
	void setCluster();
	virtual void HuberWhite() = 0;

	// get model fitting information - bcw - 4/24/13
	bool fitConverged();
	int fitNumIterations();

	/////////////////////////////////////////////////////////////////////////
	// WHY ARE THESE PUBLIC MEMBERS? bcw - 4/29/13
	/////////////////////////////////////////////////////////////////////////

	// Independent variables (can be directly manipulated...)
	vector<vector<double> > X;

	// publicly accessible parameter meta data
	vector<string> label;
	vector<int> order;
	vector<int> type;

	// index of the parameter to perform statistical tests
	int testParameter;
protected:
	Plink * P;

	// Missing flag
	vector<bool> miss;

	// number of individuals
	int nind;
	// number of parameters/coefficients: intercept + main effects + interaction
	int np; 

	bool has_snps;

	vector<bool> xchr;
	vector<bool> haploid;

	bool sex_effect;

	vector<bool> valid;
	bool all_valid;

	// beta coefficients for each term/parameter in the model
	vector_t coef; 
	// Sigma? TODO: define this! bcw - 4/29/13
	matrix_t S; 

	// Term types
	enum terms {
		INTERCEPT,
		ADDITIVE,
		DOMDEV,
		HAPLOTYPE,
		SEX,
		COVARIATE,
		INTERACTION,
		QFAM,
		NUMERIC
	};

	double buildIntercept();
	double buildAdditive(Individual *, int);
	double buildDominance(Individual *, int);
	double buildHaplotype(int, int);
	double buildSex(Individual *);
	double buildCovariate(Individual *, int);
	double buildNumeric(Individual *, int);
	double buildInteraction(Individual *, int, vector_t &);
	double buildQFAM(Individual *);

	bool skip;

	// List of additive SNP effects
	// assuming SNP major mode
	vector<int> additive;

	int mAA;
	int mAB;
	int mBB;

	double mA, mB;

	// List of dominance deviation SNP effects
	vector<int> dominance;

	// List of covariates (clist)
	vector<int> covariate;

	// List of numeric attributes - added for inbix - bcw - 4/20/13
	vector<int> numeric;

	// List of pairwise interactions
	// ( indexing previously specified components, 1,2,..)
	vector<int2> interaction;

	// List of sets of haplotypes
	vector<set<int> > haplotype;

	// Clustering information
	bool cluster;
	vector<int> clst;
	int nc;

	// new variables for model fitting information - bcw - 4/24/13
	bool converged;
	int numIterations;
};

#endif
