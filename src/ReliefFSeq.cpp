/* 
 * File:   ReliefFSeq.cpp
 * Author: Bill White
 * 
 * Created on: 7/21/12
 */

#include <iostream>
#include <iomanip>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <gsl/gsl_cdf.h>

#include "plink.h"
#include "helper.h"

#include "ReliefF.h"
#include "ReliefFSeq.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "DistanceMetrics.h"
#include "Insilico.h"
#include "Statistics.h"
#include "options.h"

using namespace std;
using namespace boost;

ReliefFSeq::ReliefFSeq(Dataset* ds, Plink* plinkPtr) :
		ReliefF::ReliefF(ds, plinkPtr, RNASEQ_ANALYSIS) {
	string configValue;
	/// set the various algorithm modes to one of four combinations:
	/// snr-snr, snr-relieff, tstat-pval, tstat-abst
	mode = par::algorithmSeqMode;
	snrMode = par::algorithmSnrMode;
	tstatMode = par::algorithmTstatMode;
  if((mode == "snr") || (mode == "tstat")) {
    cout << Timestamp() << "ReliefFSeq mode set to: " << mode << endl;
    if(mode == "snr") {
      if((snrMode != "snr") && (snrMode != "relieff")) {
        error("ERROR: Unrecognized ReliefFSeq SNR mode: " + snrMode);
      }
      cout << Timestamp() << "ReliefFSeq SNR mode set to: " << snrMode << endl;
    }
    if(mode == "tstat") {
      if((tstatMode != "pval") && ((tstatMode != "abst") || (tstatMode != "rawt"))) {
        error("ERROR: Unrecognized ReliefFSeq t-statistic mode: " + tstatMode);
      }
      cout << Timestamp() << "ReliefFSeq t-statistic mode set to: " << tstatMode << endl;
    }
  }
  else {
    error("ERROR: Unrecognized ReliefFSeq mode: " + mode);
    exit(EXIT_FAILURE);
  }
	
	/// set the s0 value
  s0 = lexical_cast<double>(par::algorithmSeqS0);
  if((s0 >= 0) || (s0 <= 1.0)) {
    cout << Timestamp() << "ReliefFSeq s0 set to: " << s0 << endl;
  }
  else {
    error("ERROR: ReliefFSeq s0 out of range (0, 1): " + dbl2str(par::algorithmSeqS0));
  }

	cout << Timestamp() << "ReliefFSeq initializing with config map" << endl;
}

ReliefFSeq::~ReliefFSeq() {
}

bool ReliefFSeq::ComputeAttributeScores() {
	// preconditions:
	// 1. case-control data
	// 2. all numeric variables

	W.resize(dataset->NumNumerics(), 0.0);

	// pre-compute all instance-to-instance distances and get nearest neighbors
	PreComputeDistances();

	// using pseudo-code notation from white board discussion - 7/21/12
	// changed to use Brett's email (7/21/12) equations - 7/23/12
	cout << Timestamp() << "Running ReliefFSeq algorithm" << endl;
	vector<string> numNames;
	numNames = dataset->GetNumericsNames();
	vector<unsigned int> numericIndices =
			dataset->MaskGetAttributeIndices(NUMERIC_TYPE);
	// DEBUG
	//	string rawScoresFileName = dataset->GetNumericsFilename() + "_rawscores.tab";
	//	ofstream outFile(rawScoresFileName.c_str());
	//	outFile << "gene\tmuMiss\tmuHit\tsigmaMiss\tsigmaHit\tnum\tden\tdms0\tsnr" << endl;

	/// run this loop on as many cores as possible through OpenMP
#pragma omp parallel for
	for (unsigned int numIdx = 0; numIdx < numericIndices.size();
			++numIdx) {
		unsigned int alpha = numericIndices[numIdx];
		pair<double, double> muDeltaAlphas = MuDeltaAlphas(alpha);
		double muDeltaHitAlpha = muDeltaAlphas.first;
		double muDeltaMissAlpha = muDeltaAlphas.second;

		pair<double, double> sigmaDeltaAlphas = SigmaDeltaAlphas(alpha,
				muDeltaHitAlpha, muDeltaMissAlpha);
		double sigmaDeltaHitAlpha = sigmaDeltaAlphas.first;
		double sigmaDeltaMissAlpha = sigmaDeltaAlphas.second;
		double snrNum = 0.0, snrDen = 0.0;
		double tstatNum = 0.0, tstatDen = 0.0;
		double alphaWeight = 0.0;
		if(mode == "snr") {
			// mode: snr (signal to noise ratio)
			snrNum = fabs(muDeltaMissAlpha - muDeltaHitAlpha);
			snrDen = sigmaDeltaMissAlpha + sigmaDeltaHitAlpha;
//			outFile << numNames[numIdx]
//					<< "\t" << muDeltaMissAlpha << "\t" << muDeltaHitAlpha
//					<< "\t" << sigmaDeltaMissAlpha << "\t" << sigmaDeltaHitAlpha
//					<< "\t" << num << "\t" << den << "\t" << (den + s0);
			if(snrMode == "snr") {
				alphaWeight = snrNum / (snrDen + s0);
			}
			else {
				alphaWeight = snrNum;
			}
		}
		else {
			// mode: tstat (t-statistic)
			// from Brett's email - 8/15/12
			// Also we could change the score to a real t-statistic:
			// (xbar1 – xbar2)/(Sp*sqrt(1/n1 + 1/n2)), 
			// where Sp = pooled standard deviation=
			// sqrt(((n1-1)*variance1 + (n2-1)*variance2)/(n1+n2-2)).
			double n1, n2;
			n1 = n2 = m * k;
			double variance1 = sigmaDeltaHitAlpha;
			double variance2 = sigmaDeltaMissAlpha;
			double pooledStdDev =
					sqrt(((n1 - 1) * variance1 + (n2 - 1) * variance2) / (n1 + n2 - 2));
			tstatNum = muDeltaMissAlpha - muDeltaHitAlpha;
			tstatDen = pooledStdDev * sqrt((1.0 / n1) + (1.0 / n2));
			// make into a t-statistic and use for pvalue
			double t = tstatNum / (tstatDen + s0);
//#pragma omp critical
//			cout << numNames[numIdx]
//					<< "\t" << muDeltaMissAlpha << "\t" << muDeltaHitAlpha
//					<< "\t" << sigmaDeltaMissAlpha << "\t" << sigmaDeltaHitAlpha
//					<< "\t" << tstatNum << "\t" << tstatDen << "\t" << t << endl;
			double df =  n1 + n2 - 2;
			double gslPval = 1.0;
			if(t < 0) {
				gslPval = gsl_cdf_tdist_P(-t, df);
			}
			else {
				gslPval = gsl_cdf_tdist_P(t, df);
			}
			// assign the variable a weight for ReliefF
			if(tstatMode == "pval") {
				// use 1-pvalue as the attribute scrore
				alphaWeight = 1.0 - (2.0 * (1.0 - gslPval));
			}
			else {
				if(tstatMode == "abst") {
					// use absolute value of the t statistic as the weight
					alphaWeight = fabs(t);
				}
				else {
					// use raw value of the t statistic as the weight
					alphaWeight = t;
				}
			}
		}
	
		/// assign a weight to this variable index
		W[numIdx] = alphaWeight;
		// DEBUG
		// outFile << "\t" << W[numIdx] << endl;
	} // for all gene alpha

	// DEBUG
	// outFile.close();

	return true;
}

AttributeScores ReliefFSeq::GetScores() {
	AttributeScores returnScores;
	vector<string> numNames;
	numNames = dataset->GetNumericsNames();
	for(unsigned int alpha = 0; alpha < dataset->NumNumerics(); ++alpha) {
		returnScores.push_back(make_pair(W[alpha], numNames[alpha]));
	}
	return returnScores;
}

pair<double, double> ReliefFSeq::MuDeltaAlphas(unsigned int alpha) {

	// for all instances
	double missSum = 0.0;
	double hitSum = 0.0;
	for(unsigned int i = 0; i < m; ++i) {
		DatasetInstance* S_i = dataset->GetInstance(i);

		// get hits and misses for this instance
		vector<unsigned int> hits(k);
		vector<unsigned int> misses(k);
		map<ClassLevel, vector<unsigned int> > allMisses;
		S_i->GetNNearestInstances(k, hits, allMisses);
		misses = allMisses.begin()->second;

		// sum over all nearest hits and misses neighbors
		for(unsigned int j = 0; j < k; ++j) {
			hitSum += diffManhattan(alpha, S_i, dataset->GetInstance(hits[j]));
			missSum += diffManhattan(alpha, S_i, dataset->GetInstance(misses[j]));
		}
	}

	// return the average of the hit and miss diffs
	pair<double, double> returnValues;
	double avgFactor = 1.0 / ((double) m * (double) k);
	returnValues.first = hitSum * avgFactor;
	returnValues.second = missSum * avgFactor;

	return returnValues;
}

pair<double, double> ReliefFSeq::SigmaDeltaAlphas(unsigned int alpha,
		double muDeltaHit, double muDeltaMiss) {

	double missSum = 0.0;
	double hitSum = 0.0;

	// for all instances
	for(unsigned int i = 0; i < m; ++i) {
		DatasetInstance* S_i = dataset->GetInstance(i);
		/// get hits and misses for this instance
		vector<unsigned int> hits(k);
		map<ClassLevel, vector<unsigned int> > allMisses;
		S_i->GetNNearestInstances(k, hits, allMisses);
		// for all nearest neighbor hits
		for(unsigned int j = 0; j < k; ++j) {
			hitSum +=
					pow(
							(diffManhattan(alpha, S_i, dataset->GetInstance(hits[j]))
									- muDeltaHit), 2);
		}
		// for all nearest neighbor misses
		// assume only one other miss class
		vector<unsigned int> misses(k);
		misses = allMisses.begin()->second;
		for(unsigned int j = 0; j < k; ++j) {
			missSum += pow(
					(diffManhattan(alpha, S_i, dataset->GetInstance(misses[j]))
							- muDeltaMiss), 2);
		}
	}

	// return the standard deviation of the hit and miss diffs
	pair<double, double> returnValues;
	double avgFactor = 1.0 / ((double) m * (double) k);
	returnValues.first = sqrt(hitSum * avgFactor);
	returnValues.second = sqrt(missSum * avgFactor);

	return returnValues;
}