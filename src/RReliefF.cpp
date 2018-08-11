/* 
 * File:   RReliefF.cpp
 * Author: billwhite
 * 
 * Created on September 27, 2011, 9:21 PM
 */

#include <iostream>
#include <vector>

#include "plink.h"
#include "helper.h"

#include "ReliefF.h"
#include "RReliefF.h"
#include "Dataset.h"
#include "DistanceMetrics.h"
#include "Insilico.h"

using namespace std;

RReliefF::RReliefF(Dataset* ds, Plink* plinkPtr) :
		ReliefF::ReliefF(ds, plinkPtr, REGRESSION_ANALYSIS) {
	PP->printLOG(Timestamp() + "RReliefF initialization\n");
	if (!ds->HasContinuousPhenotypes()) {
		error("ERROR: Attempting to construct RReliefF object without a continuous phenotype data set");
	}
}

RReliefF::~RReliefF() {
}

bool RReliefF::ComputeAttributeScores() {
  PP->printLOG(Timestamp() + "Regression Relief-F ComputeAttributeScores() START\n");
	// pre-compute all instance-to-instance distances and get nearest neighbors
	PreComputeDistances();

	// results are stored in scores, raw weights in W
	W.resize(dataset->NumVariables(), 0.0);
  scores.clear();

	double dblM = static_cast<double>(m);
  vector<uint> attributeIndices = dataset->MaskGetAttributeIndices(DISCRETE_TYPE);
  vector<string> attributeNames = dataset->MaskGetAttributeNames(DISCRETE_TYPE);
  vector<uint> numericIndices = dataset->MaskGetAttributeIndices(NUMERIC_TYPE);
  vector<string> numericNames = dataset->MaskGetAttributeNames(NUMERIC_TYPE);
  uint numAttributes = attributeIndices.size();
  uint numNumerics = numericIndices.size();
  
  /* using pseudocode notation from paper
	 *
	 * Used to hold the probability of a different class value given nearest
	 * instances (numeric class)
	 */
	double ndc = 0.0;
	/**
	 * Used to hold the probability of different value of an attribute given
	 * nearest instances (numeric class case)
	 */
	vector<double> nda;
	nda.resize(dataset->NumVariables(), 0.0);
	/**
	 * Used to hold the probability of a different class value and different 
   * attribute value given nearest instances (numeric class case)
	 */
	vector<double> ndcda;
	ndcda.resize(dataset->NumVariables(), 0.0);

	// pointer to the (possibly random) instance being sampled
	PP->printLOG(Timestamp() + "Running RRelief-F algorithm:\n");
	vector<string> instanceIds = dataset->GetInstanceIds();
  #pragma omp parallel for
	for (uint i = 0; i < (uint) m; i++) {

    DatasetInstance* R_i = NULL;
		if (randomlySelect) {
			// randomly sample an instance (without replacement?)
			R_i = dataset->GetRandomInstance();
		} else {
			// deterministic/indexed instance sampling, ie, every instance against
			// every other instance
			uint instanceIndex;
			dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
			R_i = dataset->GetInstance(instanceIndex);
		}
		if (!R_i) {
			error("ERROR: Random or indexed instance count not be found for index: ["
					+ int2str(i) + "]\n");
		}

		// K NEAREST NEIGHBORS
		// find k nearest neighbors
		vector<uint> nNearestNeighbors;
		bool canGetNeighbors = R_i->GetNNearestInstances(k, nNearestNeighbors);
		if (!canGetNeighbors) {
			error("ERROR: Cannot get " + int2str(k) + " nearest neighbors\n");
		}
		// check algorithm preconditions
		if (nNearestNeighbors.size() < 1) {
			error("ERROR: No nearest hits found\n");
		}
		if (nNearestNeighbors.size() < k) {
			error("ERROR: Could not find enough neighbors\n");
		}

		// update: using pseudocode notation
		for (uint j = 0; j < k; ++j) {
			// get the jth nearest neighbor
			DatasetInstance* I_j = dataset->GetInstance(nNearestNeighbors[j]);
			double diffPredicted = diffPredictedValueTau(R_i, I_j);
			double d_ij = R_i->GetInfluenceFactorD(j);
#pragma omp critical
      {
        ndc += (diffPredicted * d_ij);
      }
			uint scoresIndex = 0;
      // ---------------------
			// attributes
			for (uint attrIdx = 0; attrIdx < numAttributes;
					++attrIdx) {
				uint A = attributeIndices[attrIdx];
				double attrScore = snpDiffFuncPtr(A, R_i, I_j) * d_ij;
#pragma omp critical
        {
          nda[scoresIndex] += attrScore;
          ndcda[scoresIndex] += (diffPredicted * attrScore);
        }
        if(par::algorithm_verbose) {
          cout << "(i, j) = (" << i << "," << j << ") =>"
                  << " diff predicted: " << diffPredicted
                  << ", d_ij: " << d_ij
                  << ", ndc: " << ndc
                  << ", A: " << A
                  << ", snpDiff: " << snpDiffFuncPtr(A, R_i, I_j)
                  << ", nda[A}: " << nda[scoresIndex]
                  << " ndcda[A]: " << ndcda[scoresIndex]
                  << endl;
        }
				++scoresIndex;
			}
      // ---------------------
			// numerics
			for (uint numIdx = 0; numIdx < numericIndices.size(); ++numIdx) {
				uint N = numericIndices[numIdx];
				double numScore = numDiffFuncPtr(N, R_i, I_j) * d_ij;
#pragma omp critical
        {
          nda[scoresIndex] += numScore;
          ndcda[scoresIndex] += (diffPredicted * numScore);
        }
        if(par::algorithm_verbose) {
          cout << "(i, j) = (" << i << "," << j << ") =>"
                  << " diff predicted: " << diffPredicted
                  << ", d_ij: " << d_ij
                  << ", N: " << N
                  << ", snpDiff: " << numDiffFuncPtr(N, R_i, I_j)
                  << ", nda[N}: " << nda[scoresIndex]
                  << " ndcda[N]: " << ndcda[scoresIndex]
                  << endl;
        }
				++scoresIndex;
			}
		}
		// happy lights
		if (i && ((i % 100) == 0)) {
			PP->printLOG(Timestamp() + int2str(i) + "/" + int2str(m)  + "\n");
		}
	}
	PP->printLOG(Timestamp() + int2str(m) + "/" + int2str(m) + " done\n");

	PP->printLOG(Timestamp() + "Computing final scores\n");
  uint scoresIdx = 0;
  #pragma omp parralel for
  for (uint attIdx = 0; attIdx < numAttributes; ++attIdx) {
    uint A = attributeIndices[attIdx];
    string attributeName = attributeNames[A];
#pragma omp critical
    {
      double tempW = (ndcda[scoresIdx] / ndc) - 
        ((nda[scoresIdx] - ndcda[scoresIdx]) / (dblM - ndc));
      if(std::isnan(tempW)) {
        error("WARNING: detected [NaN] in weight calculation, using zero instead\n");
        W[scoresIdx] = 0.0;
      } else {
        W[scoresIdx] = tempW;
      }
      scores.push_back(make_pair(W[scoresIdx], attributeName));
      ++scoresIdx;
      // happy lights
      if (attIdx && ((attIdx % 100) == 0)) {
        PP->printLOG(Timestamp() + int2str(attIdx) + "/" + int2str(numAttributes)  + "\n");
      }
    }
  }
  PP->printLOG(Timestamp() + int2str(numAttributes) + "/" + int2str(numAttributes)  + "\n");
  #pragma omp parralel for
  for (uint numIdx = 0; numIdx < numNumerics; ++numIdx) {
    uint N = numericIndices[numIdx];
    string numericName = numericNames[N];
#pragma omp critical
    {
      double tempW = 
        (ndcda[scoresIdx] / ndc) - ((nda[scoresIdx] - ndcda[scoresIdx]) / 
          (dblM - ndc));
      if(std::isnan(tempW)) {
        error("WARNING: detected [NaN] in weight calculation, using zero instead\n");
        W[scoresIdx] = 0.0;
      } else {
        W[scoresIdx] = tempW;
      }
      scores.push_back(make_pair(W[scoresIdx], numericName));
      ++scoresIdx;
      // happy lights
      if (numIdx && ((numIdx % 100) == 0)) {
        PP->printLOG(Timestamp() + int2str(numIdx) + "/" + int2str(numNumerics)  + "\n");
      }
    }
  }
  PP->printLOG(Timestamp() + int2str(numNumerics) + "/" + int2str(numNumerics)  + "\n");

  // normalize scores if flag set
  if(normalizeScores) {
    PP->printLOG(Timestamp() + "Normalizing scores\n");
    NormalizeScores();
  }

  PP->printLOG(Timestamp() + "Relief-F ComputeAttributeScores() END\n");
  
	return true;
}
