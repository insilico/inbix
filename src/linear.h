//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#ifndef __LINEAR_H__
#define __LINEAR_H__

#include <vector>

#include "plink.h"
#include "model.h"

using namespace std;

class LinearModel : public Model {
public:
	LinearModel(Plink *);
	~LinearModel() {
	};

	void setDependent();
	void fitLM();
	void fitUnivariateLM();

	void reset();
	void pruneY();
	vector_t getCoefs();
	vector_t getVar();
	vector_t getSE();
	void displayResults(ofstream &, Locus *);
	vector_t getPVals();
	double getPValue();
	void HuberWhite();

	void standardise();
	double calculateRSS();
	double calculateRSquared();
	double calculateAdjustedRSquared();
	double calculateMallowC(LinearModel *);
	double calculateFTest(LinearModel *);
private:
	vector_t Y;
	vector<int> C; // THIS IS NEVER USED! bcw 4/29/13
	vector<double> se; // THIS IS NEVER USED! bcw 4/29/13
	double chisq;

	vector<double> sig;
	vector<double> w;
	vector<vector<double> > u;
	vector<vector<double> > v;

	double varY;
	double meanY;

	double RSS;

	void function(const int i, vector<double> & p); // <- NEVER USED bcw 4/29/13
	void setVariance();
};


#endif
