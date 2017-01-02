/*
 * ReliefF.cpp - Bill White - 7/16/05
 * 
 * ReliefF algorithm implementation.
 *
 * Totally redone for the McKinney In Silico Lab in 2011.
 * Using OpenMP for multi-core parallelization - April 2011.
 * Integration into inbix Fall 2016.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <sstream>

#include <omp.h>

#include "ReliefF.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "StringUtils.h"
#include "DistanceMetrics.h"
#include "Insilico.h"

#include "plink.h"
#include "options.h"
#include "helper.h"

using namespace std;
using namespace insilico;

/// scores map: score->attribute index
typedef vector<pair<double, unsigned int> > ScoresMap;
/// scores map iterator
typedef vector<pair<double, unsigned int> >::iterator ScoresMapIt;
/// attribute index map: attribute index->score
typedef vector<pair<unsigned int, double> > AttributeIndex;
/// attribute index map iterator
typedef vector<pair<unsigned int, double> >::const_iterator AttributeIndexIt;

/// attribute score sorting functor
bool scoreSort(const ScoreVarPair& p1, const ScoreVarPair& p2) {
  return p1.first < p2.first;
}

/// attribute index sorting functor
bool attributeSort(const pair<unsigned int, double>& p1,
        const pair<unsigned int, double>& p2) {
  return p1.first < p2.first;
}

/// functor for T comparison - instance-to-instance distance
typedef pair<unsigned int, DatasetInstance*> T;
class deref_less : public std::binary_function<T, T, bool> {
public:
  bool operator()(const T a, const T b) const {
    return(a.first < b.first);
  }
};

ReliefF::ReliefF(Dataset* ds, Plink* plinkPtr, AnalysisType anaType):
AttributeRanker::AttributeRanker(ds) {
  PP->printLOG(Timestamp() + "ReliefF initialization from Plink parameters and Dataset pointer\n");
  PP = plinkPtr;
  if(ds) {
    dataset = ds;
  } else {
    error("ReliefF constructor dataset pointer is not initialized\n");
  }
  analysisType = anaType;

  m = dataset->NumInstances();
  randomlySelect = true;
  PP->printLOG(Timestamp() + "Number of samples: m = " + int2str(m) + "\n");
  if(m == 0 || m == ds->NumInstances()) {
    // sample deterministically unless a sample size has been set
    PP->printLOG(Timestamp() + "Sampling all instances deterministically\n");
    randomlySelect = false;
    m = ds->NumInstances();
  } else {
    PP->printLOG(Timestamp() + "Sampling instances randomly\n");
    randomlySelect = true;
  }

  weightByDistanceMethod = par::weightByDistanceMethod;

  // default k, in options.h/.cpp
  k = par::k;
  if(k) {
    PP->printLOG(Timestamp() + "Number of nearest neighbors: k = " + int2str(k) + "\n");
    // k nearest neighbors and m randomly selected instances
    // spread differences and thus weight updates
    // over (m x k) iterations
    // double one_over_m_times_k = 1.0 / (((double) m) * ((double) k));
    //                               m       *           k
  } else {
    PP->printLOG(Timestamp() + "k nearest neighbors will be optimized\n");
  }

  numTarget = ds->NumVariables();
  PP->printLOG(Timestamp() + "Number of attributes: " + int2str(numTarget) + "\n");

  if(par::do_iterative_removal) {
    numTarget = par::relieffNumTarget;
    unsigned int numPredictors = dataset->NumVariables();
    if((numTarget < 1) || (numTarget > numPredictors)) {
      if(numTarget == 0) {
        numTarget = numPredictors;
      } else {
            error("Target number of variables out of range: " + int2str(numTarget) + "\n");
      }
    }
  
    if(par::relieffIterNumToRemove) {
      removePerIteration = par::relieffIterNumToRemove;
      PP->printLOG(Timestamp() + "Iteratively removing " + int2str(removePerIteration) + "\n");
      doRemovePercent = false;
    } else {
      removePercentage = ((double) par::relieffIterPercentToRemove) / 100.0;
      removePerIteration =
               (unsigned int) ((double) dataset->NumAttributes()
               * removePercentage + 0.5);
      PP->printLOG(Timestamp() + "Iteratively removing " + 
              int2str(removePercentage * 100) + "% = " + 
              int2str(removePerIteration) + "\n");
      doRemovePercent = true;
    }
    if((removePerIteration < 1)
            || (removePerIteration >= numPredictors)) {
      error("Number to remove per iteration [" + int2str(removePerIteration) +
            "] not in valid range 1 < n < " + int2str(numPredictors));
    }
  }

  /// set the SNP metric function pointer based on command line params or defaults
  snpMetric = par::snpMetric;
  numMetric = par::numMetric;
  bool snpMetricFunctionUnset = true;
  if(snpMetricFunctionUnset && to_upper(snpMetric) == "GM") {
    snpDiff = diffGMM;
    snpMetricFunctionUnset = false;
  }
  if(snpMetricFunctionUnset && to_upper(snpMetric) == "AM") {
    snpDiff = diffAMM;
    snpMetricFunctionUnset = false;
  }
  if(snpMetricFunctionUnset && to_upper(snpMetric) == "NCA") {
    snpDiff = diffNCA;
    snpMetricFunctionUnset = false;
  }
  if(snpMetricFunctionUnset && to_upper(snpMetric) == "NCA6") {
    snpDiff = diffNCA6;
    snpMetricFunctionUnset = false;
  }
  if(snpMetricFunctionUnset && to_upper(snpMetric) == "GRM") {
    // no need to set a function pointer here for GRM
    error("GCTA GRM metric is not allowed in weight update metric, only nearest neighbors\n");
  }
  if(snpMetricFunctionUnset && to_upper(snpMetric) == "KM") {
    error("ERROR: KM is not supported as a ReliefF metric\n");
    // snpDiff = diffKM;
    // snpMetricFunctionUnset = false;
  }
  if(snpMetricFunctionUnset) {
    error("Cannot set SNP metric to [" + snpMetric + "]\n");
  }
  if(to_upper(numMetric) == "MANHATTAN") {
    numDiff = diffManhattan;
  } else {
    if(to_upper(numMetric) == "EUCLIDEAN") {
      numDiff = diffEuclidean;
    } else {
      error("[" + numMetric + "] is not a valid numeric metric type\n");
    }
  }

  PP->printLOG(Timestamp() + "ReliefF SNP distance metric: " + snpMetric + "\n");
  PP->printLOG(Timestamp() + "ReliefF continuous distance metric: " + numMetric + "\n");
  
  weightByDistanceMethod = par::weightByDistanceMethod;
  if((weightByDistanceMethod != "exponential")
          && (weightByDistanceMethod != "equal")) {
    error("ERROR: Invalid --weight-by-distance-method: " + weightByDistanceMethod);
  }
  weightByDistanceSigma = static_cast<double>(par::weightByDistanceSigma);
  PP->printLOG(Timestamp() + "Weight by distance method: " + weightByDistanceMethod);
  if(weightByDistanceMethod == "exponential") {
    PP->printLOG(", using sigma = " + dbl2str(weightByDistanceSigma) + "\n");
  } else {
    PP->printLOG("\n");
  }

	PP->printLOG(Timestamp() + "ReliefF has " + int2str(omp_get_num_procs()) + " threads\n");

  vector<string> atrNames = dataset->GetAttributeNames();
  vector<string> numNames = dataset->GetNumericsNames();
  scoreNames.resize(atrNames.size() + numNames.size());
  copy(atrNames.begin(), atrNames.end(), scoreNames.begin());
  copy(numNames.begin(), numNames.end(), scoreNames.begin() + atrNames.size());
}

ReliefF::~ReliefF() {
}

bool ReliefF::ComputeAttributeScores() {
  PP->printLOG(Timestamp() + "Relief-F ComputeAttributeScores() START\n");
  
  // changed from matrix to map for ID matching - November 2011
  PreComputeDistances();

  /// algorithm line 1
	W.clear();
  W.resize(dataset->NumVariables(), 0.0);

  // pointer to the instance being sampled
  DatasetInstance* R_i;
  int i = 0;
  
  double one_over_m_times_k = 1.0 / (((double) m) * ((double) k));
  PP->printLOG(Timestamp() + "Averaging factor 1/(m*k): "  + 
    dbl2str(one_over_m_times_k) + "\n");

  vector<string> instanceIds = dataset->GetInstanceIds();
//  for(i = 0; i < (int) m; i++) {
//    cout << instanceIds[i] + "\n");
//  }
  /// algorithm line 2
  for(i = 0; i < (int) m; i++) {
    // algorithm line 3
    if(randomlySelect) {
      // randomly sample an instance (without replacement?)
      R_i = dataset->GetRandomInstance();
    } else {
      // deterministic/indexed instance sampling, ie, every instance against
      // every other instance
      unsigned int instanceIndex;
      dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
      R_i = dataset->GetInstance(instanceIndex);
    }
    if(!R_i) {
      error("Random or indexed instance count not be found for index: [" + int2str(i) + "]\n");
    }
    ClassLevel class_R_i = R_i->GetClass();

    /// algorithm lines 4, 5 and 6
    // find k nearest hits and nearest misses
    vector<unsigned int> hits;
    map<ClassLevel, vector<unsigned int> > misses;
    bool canGetNeighbors = false;
    canGetNeighbors = R_i->GetNNearestInstances(k, hits, misses);
    //    cout << "Instance class: " << R_i->GetClass() << ", hits: ";
    //    for(unsigned int ii = 0; ii < hits.size(); ++ii) {
    //      cout << hits[ii] << " ";
    //    }
    //    cout + "\n") << "Misses:\n");
    //    map<ClassLevel, vector<unsigned int> >::const_iterator iit;
    //    for(iit = misses.begin(); iit != misses.end(); ++iit) {
    //      cout << "Class: " << iit->first << ", misses: ";
    //      vector<unsigned int> ids = iit->second;
    //      for(unsigned int jj = 0; jj < ids.size(); ++jj) {
    //        cout << ids[jj] << " ";
    //      }
    //      cout + "\n");
    //    }

    if(!canGetNeighbors) {
      error("relieff cannot get " + int2str(k) + " nearest neighbors\n");
    }

    // check algorithm preconditions
    if(hits.size() < 1) {
      error("No nearest hits found\n");
    }
    if(hits.size() < k) {
      error("Could not find enough neighbors that are hits\n");
    }
    map<ClassLevel, vector<unsigned int> >::const_iterator it;
    for(it = misses.begin(); it != misses.end(); ++it) {
      vector<unsigned int> missIds = it->second;
      if(missIds.size() < 1) {
        error("No nearest misses found\n");
      }
      if(missIds.size() < k) {
        error("Could not find enough neighbors that are misses\n");
      }
      if(missIds.size() != hits.size()) {
        error("Could not find equal number of neighbors for hits and misses:" +
                int2str(hits.size()) + " vs. " + int2str(misses.size()) + "\n");
      }
    }

    // UPDATE WEIGHTS FOR ATTRIBUTE 'A' BASED ON THIS AND NEIGHBORING INSTANCES
    // update weights/relevance scores for each attribute averaged
    // across k nearest neighbors and m (possibly randomly) selected instances
    unsigned int A = 0;
    unsigned int scoresIdx = 0;
    if(dataset->HasGenotypes()) {
      vector<unsigned int> attributeIndicies =
              dataset->MaskGetAttributeIndices(DISCRETE_TYPE);
      /// algorithm line 7
      for(unsigned int attrIdx = 0; attrIdx < attributeIndicies.size();
              ++attrIdx) {
        A = attributeIndicies[attrIdx];
        double hitSum = 0.0, missSum = 0.0;
        /// algorithm line 8
        for(unsigned int j = 0; j < k; j++) {
          DatasetInstance* H_j = dataset->GetInstance(hits[j]);
          double rawDistance = snpDiff(A, R_i, H_j);
          hitSum += (rawDistance * one_over_m_times_k);
        }
        /// algorithm line 9
        map<ClassLevel, vector<unsigned int> >::const_iterator mit;
        for(mit = misses.begin(); mit != misses.end(); ++mit) {
          ClassLevel C = mit->first;
          vector<unsigned int> missIds = mit->second;
          double P_C = dataset->GetClassProbability(C);
          double P_C_R = dataset->GetClassProbability(class_R_i);
          double adjustmentFactor = P_C / (1.0 - P_C_R);
          double tempSum = 0.0;
          for(unsigned int j = 0; j < k; j++) {
            DatasetInstance* M_j = dataset->GetInstance(missIds[j]);
            double rawDistance = snpDiff(A, R_i, M_j);
            tempSum += (rawDistance * one_over_m_times_k);
          } // nearest neighbors
          missSum += (adjustmentFactor * tempSum);
        }

        W[scoresIdx] = W[scoresIdx] - hitSum + missSum;
        ++scoresIdx;
      } // all attributes
    } // has genotypes

    // loop here for numeric attributes if they exist - 6/19/11
    if(dataset->HasNumerics()) {
      vector<unsigned int> numericIndices =
              dataset->MaskGetAttributeIndices(NUMERIC_TYPE);
      for(unsigned int numIdx = 0; numIdx < numericIndices.size();
              ++numIdx) {
        A = numericIndices[numIdx];
        double hitSum = 0.0, missSum = 0.0;
        for(unsigned int j = 0; j < k; j++) {
          DatasetInstance* H_j = dataset->GetInstance(hits[j]);
          hitSum += (numDiff(A, R_i, H_j) * one_over_m_times_k);
        }

        map<ClassLevel, vector<unsigned int> >::const_iterator mit;
        for(mit = misses.begin(); mit != misses.end(); ++mit) {
          ClassLevel C = mit->first;
          vector<unsigned int> missIds = mit->second;
          double P_C = dataset->GetClassProbability(C);
          double P_C_R = dataset->GetClassProbability(class_R_i);
          double adjustmentFactor = P_C / (1.0 - P_C_R);
          double tempSum = 0.0;
          for(unsigned int j = 0; j < k; j++) {
            DatasetInstance* M_j = dataset->GetInstance(missIds[j]);
            tempSum += (numDiff(A, R_i, M_j) * one_over_m_times_k);
          } // nearest neighbors
          missSum += (adjustmentFactor * tempSum);
        }
        W[scoresIdx] = W[scoresIdx] - hitSum + missSum;
        ++scoresIdx;
      }
    } // has numerics

    // happy lights
    if(i && ((i % 100) == 0)) {
      PP->printLOG(Timestamp() + int2str(i) + "/" + int2str(m) + "\n");
    }

  } // number to randomly select -or all instances
  PP->printLOG(Timestamp() + int2str(i) + "/" + int2str(m) + " done\n");

  // superclass variable scores is what gets returned to callers of 
  // AttributeRanker and subclasses
  scores.clear();
  for(uint i=0; i < W.size(); ++i) {
    scores.push_back(make_pair(W[i], scoreNames[i]));
  }

  // normalize scores is flag set
  if(normalizeScores) {
    PP->printLOG(Timestamp() + "Normalizing scores\n");
    NormalizeScores();
  }

  PP->printLOG(Timestamp() + "Relief-F ComputeAttributeScores() END\n");

  return true;
}

bool ReliefF::ComputeAttributeScoresIteratively() {
  // final scores after all iterations
  std::map<std::string, double> finalScores;

  // save the current dataset mask
  dataset->MaskPushAll();

  // IterativeReliefF or TuRF (Tuned Relief-F)
  unsigned int iterations = 1;
  while(dataset->NumVariables() > 0) {

    PP->printLOG(Timestamp() +
            "------------------------------------------------------------" +
            "-----------------------------------------\n");
    PP->printLOG(Timestamp() + "[" + int2str(iterations) + "] Working attributes: " +
            int2str(dataset->NumVariables()) + "\n");

    ComputeAttributeScores();
    vector<ScoreVarPair > attributeScores = GetScores();

    // save worst attributes and remove from consideration on next iteration
    sort(attributeScores.begin(), attributeScores.end(), scoreSort);
    unsigned int removeThisIteration = 0;
    if(dataset->NumVariables() < removePerIteration) {
      removeThisIteration = dataset->NumVariables();
    } else {
      if(doRemovePercent) {
        removeThisIteration =
                (unsigned int) ((double) dataset->NumAttributes()
                * removePercentage + 0.5);
      } else {
        removeThisIteration = removePerIteration;
      }
    }
    for(unsigned int i = 0; i < removeThisIteration; ++i) {
      string attributeToDelete = attributeScores[i].second;
      // cout << "\t\t\tremoving attribute: " << attributeToDelete + "\n");

      if(!dataset->MaskRemoveVariable(attributeToDelete)) {
        error(string("ERROR: ReliefF::ComputeAttributeScoresIteratively: ") +
              string("could not find attribute name in data set: ") + 
              attributeToDelete + "\n");
      }
      finalScores[attributeToDelete] = attributeScores[i].first;
    }

    ++iterations;
    //ResetForNextIteration();
  } // iterate

  // populate finalScores with remaining scores
  vector<double>::const_iterator scoresIt;
  vector<string> varNames = dataset->GetVariableNames();
  for(unsigned int i = 0; i < varNames.size(); ++i) {
    if(par::verbose) PP->printLOG(varNames[i] + " => " + dbl2str(scoresIt[i]) + "\n");
    finalScores[varNames[i]] = W[i];
  }

  W.resize(scoreNames.size());
  for(unsigned int i = 0; i < scoreNames.size(); ++i) {
    if(finalScores.find(scoreNames[i]) == finalScores.end()) {
      error("Variable name " + scoreNames[i] + " could not be looked up");
    }
    W[i] = finalScores[scoreNames[i]];
  }

  // restore the dataset attribute mask
  dataset->MaskPopAll();

  scores.clear();
  for(auto i=0; i < W.size(); ++i) {
    scores.push_back(make_pair(W[i], scoreNames[i]));
  }

  return true;
}

bool ReliefF::ComputeAttributeScoresKopt() {
  PP->printLOG(Timestamp() + "Running Relief-F with kopt to determine best k\n");
  // set the optimization parameters from the command line parameters
  if(!SetKoptParameters()) {
    return false;
  }

  // iterate over all k's
  //vector<map<string, double> > allScores;
  vector<unsigned int> koptValues;
	bool hasNames = false;
	vector<vector<double> > allScores;
	vector<string> scoreNames;
  for(unsigned int thisK = koptBegin; thisK <= koptEnd; thisK += koptStep) {
    // run ReliefF on this k
    PP->printLOG(Timestamp() + "--------------------------\n");
    PP->printLOG(Timestamp() + "Running ReliefSeq for k=" + int2str(thisK) + "\n");
    k =thisK;
    koptValues.push_back(thisK);
    scores.clear();
		dataset->ResetNearestNeighbors();
    scores = ComputeScores();
	  sort(scores.begin(), scores.end(), scoresSortAscByName);
		vector<double> thisScores;
		AttributeScoresCIt scoresIt = scores.begin();
		for(; scoresIt != scores.end(); ++scoresIt) {
			if(!hasNames) {
				scoreNames.push_back(scoresIt->second);
			}
			thisScores.push_back(scoresIt->first);
		}
		allScores.push_back(thisScores);
	
		// I/O
    if(par::do_write_each_k_scores) {
      stringstream filePrefix;
      filePrefix << par::output_file_name << "." << thisK;
      WriteAttributeScores(filePrefix.str());
    }
    // PrintScores();
    hasNames = true;
  }

  // print allScores
//	for(unsigned int i=0; i < koptValues.size(); ++i) {
//		for(unsigned int j=0; j < scoreNames.size(); ++j) {
//      cout << allScores[i][j] << " ";
//    }
//    cout + "\n");
//  }
  
  // pick best scores and k's for each attribute
  scores.clear();
	for(unsigned int i=0; i < scoreNames.size(); ++i) {
    string thisVar = scoreNames[i];
		unsigned int bestK = koptValues[0];
    // cout << thisVar;
		double bestScore = -1.0;
		for(unsigned int j=0; j < koptValues.size(); ++j) {
			unsigned int thisK = koptValues[j];
			double thisScore = allScores[j][i];
      if(thisScore > bestScore) {
        bestScore = thisScore;
        bestK = thisK;
      }
    }
    // cout << "\t" << bestScore << " (" << bestK << ")\n");
    scores.push_back(make_pair(bestScore, thisVar));
    bestKs[thisVar] = bestK;
  }

  // removed sort bcw 12/21/16
  //sort(scores.begin(), scores.end(), scoresSortDesc);
  
  if(par::do_write_best_k) {
    WriteBestKs(par::output_file_name);
  }
  
  return true;
}

bool ReliefF::ResetForNextIteration() {
  PP->printLOG(Timestamp() + "***** ResetForNextIteration *****\n");
  PreComputeDistances();
  return true;
}

void ReliefF::PrintAttributeScores(ofstream& outFile) {
  unsigned int nameIdx = 0;
  AttributeScores scoresMap;
  vector<double>::const_iterator scoresIt = W.begin();
  for(; scoresIt != W.end(); ++scoresIt) {
    scoresMap.push_back(make_pair(*scoresIt, scoreNames[nameIdx]));
    ++nameIdx;
  }
  sort(scoresMap.begin(), scoresMap.end(), scoresSortDesc);
  AttributeScoresCIt smIt=scoresMap.begin();
  for(; smIt != scoresMap.end(); ++smIt) {
    outFile << smIt->first << "\t" << smIt->second << endl;
  }
}

void ReliefF::WriteAttributeScores(string baseFilename) {
  string resultsFilename = baseFilename;
  if(dataset->HasContinuousPhenotypes()) {
    resultsFilename += ".rrelieff.tab";
  } else {
    resultsFilename += ".relieff.tab";
  }
	PP->printLOG(Timestamp() + "Writing Relief-F results to: " + resultsFilename + "\n");

  ofstream outFile;
  outFile.open(resultsFilename.c_str());
  if(outFile.bad()) {
    error("ERROR: Could not open scores file " + resultsFilename + " for writing\n");
  }
  //sort(W.begin(), W.end(), std::greater<double>());
  PrintAttributeScores(outFile);
  outFile.close();
}

bool ReliefF::PreComputeDistances() {
  PP->printLOG(Timestamp() + "Precomputing instance distances\n");
  map<string, unsigned int> instanceMask = dataset->MaskGetInstanceMask();
  vector<string> instanceIds = dataset->MaskGetInstanceIds();
  int numInstances = instanceIds.size();

  // create a distance matrix
  PP->printLOG(Timestamp() + "Allocating distance matrix...");
  double** distanceMatrix;
  distanceMatrix = new double*[numInstances];
  for(int i = 0; i < numInstances; ++i) {
    distanceMatrix[i] = new double[numInstances];
    for(int j = 0; j < numInstances; ++j) {
      distanceMatrix[i][j] = 0.0;
    }
  }
  PP->printLOG(" done\n");

  // TCGA genetic relationship matrix (GRM))
  if(par::snpMetricNN == "grm") {
    if(dataset->NumNumerics()) {
      error("GRM distance metric is not available for numeric data");
    }
    PP->printLOG(Timestamp() + "1) Computing instance-to-instance distances " +
      "with GCTA genetic relationship matrix (GRM)\n");
    vector<double> p = dataset->GetMAFs();
    unsigned int N = dataset->NumAttributes();
    double A_jk = 0;
#pragma omp parallel for schedule(dynamic,1)
    for(int j = 0; j < numInstances; ++j) {
      // NOTE index variable names chosen to match GCTA paper
      for(int k = j; k < numInstances; ++k) {
//        unsigned int dsi1Index;
//        dataset->GetInstanceIndexForID(instanceIds[j], dsi1Index);
//        unsigned int dsi2Index;
//        dataset->GetInstanceIndexForID(instanceIds[k], dsi2Index);
        double sum = 0.0;
        for(int i = 0; i < N; ++i) {
          AttributeLevel x_ij = dataset->GetInstance(j)->GetAttribute(i);
          AttributeLevel x_ik = dataset->GetInstance(k)->GetAttribute(i);
          double p_i = p[i];
          double two_p_i = 2 * p_i;
          double summation_expr = 0;
          if(j == k) {
            summation_expr = 
                    ((x_ij * x_ij - ((1 + two_p_i)) * x_ij) + (two_p_i * two_p_i)) / 
                    (two_p_i * (1 - p_i));
          } else {
            summation_expr = 
                    ((x_ij - two_p_i) * (x_ik - two_p_i)) / 
                    (two_p_i * (1 - p_i));
          }
          sum += summation_expr;
        }
        if(j == k) {
          A_jk = 1 + (sum / N);
          distanceMatrix[j][k] = (1- A_jk);
        } else {
          A_jk = sum / N;
          distanceMatrix[j][k] = distanceMatrix[k][j] = (1 - A_jk);
        }
      }
      if(j && (j % 100 == 0)) {
        PP->printLOG(Timestamp() + int2str(j) + "/" + int2str(numInstances) + "\n");
      }
    }
    PP->printLOG(Timestamp() + int2str(numInstances) + "/" + int2str(numInstances) + " done\n");
    
    // write GRM matrix to file with output prefix
    PP->printLOG(Timestamp() + "[ " + par::output_file_name + ".grm.tab ]\n");
    ofstream outFile(par::output_file_name + ".grm.tab");
    for(int i=0; i < numInstances; ++i) {
      for(int j=0; j < numInstances; ++j) {
        if(j) {
          outFile << "\t" << distanceMatrix[i][j];  
        } else {
          outFile << distanceMatrix[i][j];  
        }
      }
      outFile << endl;
    }
    outFile.close();
  } else {
    // populate the matrix - upper triangular
    // NOTE: make complete symmetric matrix for neighbor-to-neighbor sums
    PP->printLOG(Timestamp() + "1) Computing instance-to-instance distances\n");
#pragma omp parallel for schedule(dynamic,1)
    for(int i=0; i < numInstances; ++i) {
      //cout << "Computing instance to instance distances. Row: " << i + "\n");
      // #pragma omp parallel for
      for(int j = i + 1; j < numInstances; ++j) {
        unsigned int dsi1Index;
        dataset->GetInstanceIndexForID(instanceIds[i], dsi1Index);
        unsigned int dsi2Index;
        dataset->GetInstanceIndexForID(instanceIds[j], dsi2Index);
        /// be sure to call Dataset::ComputeInstanceToInstanceDistance
        distanceMatrix[i][j] = distanceMatrix[j][i] =
                dataset->ComputeInstanceToInstanceDistance(
                dataset->GetInstance(dsi1Index),
                dataset->GetInstance(dsi2Index));
        //cout << i << ", " << j << " => " << distanceMatrix[i][j] + "\n");
      }
      if(i && (i % 100 == 0)) {
        PP->printLOG(Timestamp() + int2str(i) + "/" + int2str(numInstances) + "\n");
      }
    }
    PP->printLOG(Timestamp() + int2str(numInstances) + "/" + int2str(numInstances) + " done\n");
  }
  //  DEBUG
//    ofstream outFile;
//    outFile.open("distanceMatrix.csv");
//    for(unsigned int i=0; i < dataset->NumInstances(); ++i) {
//      for(unsigned int j=0; j < dataset->NumInstances(); ++j) {
//        if(j)
//          outFile << "," << distanceMatrix[i][j];
//        else
//          outFile << distanceMatrix[i][j];
//      }
//      outFile + "\n");
//    }
//    outFile.close();
  //  DEBUG

  // for each instance: if discrete class, store the distance sums for same
  // and different classes, else store distances to all other instances
  // (regression ReliefF)
  if(dataset->HasContinuousPhenotypes()) {
    PP->printLOG(Timestamp() + "2) Calculating continuous phenotype nearest neighbors... ");
  } else {
    // multiclass - 12/1/11
    if(dataset->NumClasses() > 2) {
      PP->printLOG(Timestamp() + "2) Calculating same and different classes nearest neighbors... ");
    } else {
      PP->printLOG(Timestamp() + "2) Calculating same and different class nearest neighbors... ");
    }
  }
  PP->printLOG("\n");

  for(int i = 0; i < numInstances; ++i) {
    unsigned int thisInstanceIndex = instanceMask[instanceIds[i]];
    DatasetInstance* thisInstance = dataset->GetInstance(thisInstanceIndex);

    DistancePair nearestNeighborInfo;
    if(dataset->HasContinuousPhenotypes()) {
      DistancePairs instanceDistances;
      for(int j = 0; j < numInstances; ++j) {
        if(i == j)
          continue;
        double instanceToInstanceDistance = distanceMatrix[i][j];
        nearestNeighborInfo = make_pair(instanceToInstanceDistance,
                instanceIds[j]);
        instanceDistances.push_back(nearestNeighborInfo);
      }
      thisInstance->SetDistanceSums(k, instanceDistances);
    } else {
      ClassLevel thisClass = thisInstance->GetClass();
      DistancePairs sameSums;
      // changed to an array for multiclass - 12/1/11
      map<ClassLevel, DistancePairs> diffSums;
      for(int j = 0; j < numInstances; ++j) {
        if(i == j)
          continue;
        double instanceToInstanceDistance = distanceMatrix[i][j];
        unsigned int otherInstanceIndex = instanceMask[instanceIds[j]];
        DatasetInstance* otherInstance = dataset->GetInstance(
                otherInstanceIndex);
        nearestNeighborInfo = make_pair(instanceToInstanceDistance, instanceIds[j]);
        if(otherInstance->GetClass() == thisClass) {
          sameSums.push_back(nearestNeighborInfo);
        } else {
          ClassLevel otherClass = otherInstance->GetClass();
          diffSums[otherClass].push_back(nearestNeighborInfo);
        }
      }
      thisInstance->SetDistanceSums(k, sameSums, diffSums);
    }

    if(i && (i % 100 == 0)) {
      PP->printLOG(Timestamp() + int2str(i) + "/" + int2str(numInstances) + "\n");
    }
  }
  PP->printLOG(Timestamp() + int2str(numInstances) + "/" + int2str(numInstances) + " done\n");

  PP->printLOG(Timestamp() + "3) Calculating weight by distance factors for nearest neighbors... \n");
  ComputeWeightByDistanceFactors();

  // release the dynamically-allocated distance matrix
  PP->printLOG(Timestamp() + "Freeing distance matrix memory...");
  for(int i = 0; i < numInstances; ++i) {
    delete[] distanceMatrix[i];
  }
  delete[] distanceMatrix;
  PP->printLOG(" done\n");

  return true;
}

// this is an adapter method that calls ComputeAttributeScores
AttributeScores ReliefF::ComputeScores() {
  ComputeAttributeScores();
  return scores;
}

bool ReliefF::ComputeWeightByDistanceFactors() {
  vector<string> instanceIds = dataset->GetInstanceIds();
  for(unsigned int i = 0; i < dataset->NumInstances(); ++i) {
    // this instance
    unsigned int instanceIndex;
    dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
    DatasetInstance* dsi = dataset->GetInstance(instanceIndex);
    vector<double> d1_ij;
    double d1_ij_sum = 0.0;
    for(unsigned int rank_j = 1; rank_j <= k; ++rank_j) {
      double d1_ij_value = 0.0;
      if(weightByDistanceMethod == "exponential") {
        double exponentArg = (double) rank_j / weightByDistanceSigma;
        d1_ij_value = exp(-(exponentArg * exponentArg));
      } else {
        if(weightByDistanceMethod == "one_over_k") {
          d1_ij_value = 1.0 / (double) rank_j;
        } else {
          // equal
          d1_ij_value = 1.0 / (double) k;
        }
      }
      d1_ij.push_back(d1_ij_value);
      d1_ij_sum += d1_ij_value;
    }
    // "normalize" the factors - divide through by the total/sum
    dsi->ClearInfluenceFactors();
    for(unsigned int neighborIdx = 0; neighborIdx < k; ++neighborIdx) {
      double influenceFactorD = d1_ij[neighborIdx] / d1_ij_sum;
      //      cout << "d_ij: " << d1_ij[neighborIdx]
      //              << ", cummulative sum: " << d1_ij_sum
      //              << ", normalized value: " << influenceFactorD + "\n");
      dsi->AddInfluenceFactorD(influenceFactorD);
    }
    //    cout << "---------------------------------------------------------\n");
  } // end all instances

  return true;
}

bool ReliefF::SetKoptParameters() {
  // TODO: wrap these conversions in try/catch exception handler
  unsigned int tempKoptBegin = par::koptBegin;
  unsigned int tempKoptEnd = par::koptEnd;
  unsigned int tempKoptStep = par::koptStep;
  // changed for continuous phenos - 7/26/15
  unsigned int kmax = dataset->NumInstances();

  // error conditions
  if(tempKoptBegin > tempKoptEnd) {
    error("k optimization begin [" + int2str(tempKoptBegin)
            + "] is greater than end [" + int2str(tempKoptEnd) + "]\n");
    return false;
  }
  if(tempKoptEnd > kmax) {
    error("k optimization end [" + int2str(tempKoptEnd) +
          "] is greater than maximum k [" + int2str(kmax) + "]\n");
    return false;
  }
  if((tempKoptBegin == tempKoptEnd) == tempKoptStep) {
    error("k optimization specified but the range "
            "and step values do not specify any iterations"
            "\n");
    return false;
  }

  // passed all error checks
  koptBegin = tempKoptBegin;
  koptEnd = tempKoptEnd;
  koptStep = tempKoptStep;
  PP->printLOG(Timestamp() + "k optimization parameters: begin: " + int2str(koptBegin)
          + ", kopt end: " + int2str(koptEnd) + ", step: " + int2str(koptStep) + "\n");

  return true;
}

unsigned int ReliefF::GetKmax() {
  map<ClassLevel, vector<unsigned int> > classIdxMap =
          dataset->GetClassIndexes();
  map<ClassLevel, vector<unsigned int> >::const_iterator mit =
          classIdxMap.begin();
  unsigned int minClassCount = mit->second.size();
  ++mit;
  for(; mit != classIdxMap.end(); ++mit) {
    if(mit->second.size() < minClassCount) {
      minClassCount = mit->second.size();
    }
  }
  return minClassCount - 1;
}

void ReliefF::PrintBestKs() {
  for(map<string, unsigned int>::const_iterator kIt = bestKs.begin();
      kIt != bestKs.end(); ++kIt) {
    PP->printLOG(kIt->first + "\t" + int2str(kIt->second) + "\n");
  }
}

void ReliefF::WriteBestKs(string baseFilename) {
  string resultsFilename = baseFilename;
  ofstream outFile;
  resultsFilename = baseFilename + ".bestk";
  outFile.open(resultsFilename.c_str());
  if(outFile.bad()) {
    error("Could not open scores file " + resultsFilename + "for writing\n");
    exit(1);
  }
  PP->printLOG(Timestamp() + "Writing Relief-F best k's to [" + resultsFilename + "]\n");
  for(map<string, unsigned int>::const_iterator kIt = bestKs.begin();
    kIt != bestKs.end(); ++kIt) {
    outFile << kIt->first << "\t" << kIt->second << endl;
  }

  outFile.close();
}

bool ReliefF::RemoveWorstAttributes(unsigned int numToRemove) {
  unsigned int numToRemoveAdj = numToRemove;
  unsigned int numAttr = dataset->NumAttributes();
  if((numAttr - numToRemove) < numTarget) {
    PP->printLOG(Timestamp() + "WARNING: attempt to remove " + 
      int2str(numToRemove) + " attributes which will remove more than target " + 
      "number of attributes " + int2str(numTarget) + ". Adjusting\n");
    numToRemoveAdj = numAttr - numTarget;
  }
  PP->printLOG(Timestamp() + "Removing " + int2str(numToRemoveAdj) + " attributes\n");
  sort(scores.begin(), scores.end(), scoresSortAsc);
  for(unsigned int i = 0; i < numToRemoveAdj; ++i) {
    // worst score and attribute name
    ScoreVarPair worst = scores[i];
    if(par::verbose) {
        PP->printLOG("\t\t\t\tReliefF removing: " + worst.second + 
          " (" + dbl2str(worst.first) + ")\n");
    }
    // save worst
    removedAttributes.push_back(worst);
    // remove the attribute from those under consideration
    if(!dataset->MaskRemoveVariableType(worst.second, DISCRETE_TYPE)) {
      error("Could not remove worst attribute: " + worst.second + "\n");
    }
  }

  return true;
}
