/* 
 * File:   EvaporativeCoolingPrivacy.cpp
 * Author: bwhite
 * 
 * Created on October 18, 2016, 10:57 PM
 */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <iterator>
#include <vector>
#include <armadillo>
#include <random>
#include <cmath>
#include <ctgmath>
#include <limits>
        
#include <boost/format.hpp>

// inbix PLINK base
#include "plink.h"
#include "options.h"
#include "helper.h"

// inbix additions to PLINK
#include "Insilico.h"
#include "ArmadilloFuncs.h"
#include "Dataset.h"
#include "EvaporativeCoolingPrivacy.h"
#include "ReliefF.h"
#include "RandomForest.h"
#include "StringUtils.h"

using namespace std;
using namespace arma;

// ----------------------------------------------------------------------------
// public methods

EvaporativeCoolingPrivacy::EvaporativeCoolingPrivacy(Dataset* trainset, 
                                                     Dataset* holdoset, 
                                                     Dataset* testset, 
                                                     Plink* plinkPtr, 
                                                     bool datasetsAreSims) {
  // --------------------------------------------------------------------------
  PP = plinkPtr;
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy constructor START\n");
  // data
  train = trainset;
  holdout = holdoset;
  test = testset;
  // pointer to a PLINK environment
  dataIsSimulated = datasetsAreSims;
  numInstances = 0;
  curVarNames = trainset->GetVariableNames();
  numVariables = trainset->NumVariables();
  numSignalsInData = numVariables;
  // --------------------------------------------------------------------------
  // if the passed data sets are simulated (for paper)
  if(datasetsAreSims) {
    PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy constructor: simulated data\n");
    if(trainset->NumVariables() != holdoset->NumVariables()) {
      if(testset) {
        if(trainset->NumVariables() != testset->NumVariables() ||
            holdoset->NumVariables() != testset->NumVariables()) {
                error("Training, holdout and testing data sets must have the same number of variables\n");
        }
      }
    }
    // special variables for simulated data sets
    PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy variables: " + int2str(numVariables) + "\n");
    PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy percent signal: " + dbl2str(par::ecPrivacyPercentSignal) + "\n");
    numSignalsInData = (uint) ((double) numVariables * par::ecPrivacyPercentSignal);
    PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy signals: " + int2str(numSignalsInData) + "\n");
    signalNames.clear();
    for(uint sigNum=0; sigNum < numSignalsInData; ++sigNum) {
      signalNames.push_back(curVarNames[sigNum]);
    }
  }
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy initialize\n");
  // --------------------------------------------------------------------------
  // algorithm
  Q_EPS = 0.005;
  MAX_ITERATIONS = 1e6;
  updateInterval = par::ecPrivacyUpdateFrequency;
  iteration = 0;
  update = 0;
  deltaQ = 0;
  threshold = 0;
  tolerance = 0;
  startTemp = par::ecPrivacyStartTemp;
  currentTemp = startTemp;
  finalTemp = par::ecPrivacyFinalTemp;
  // Trang: larger tau takes longer to get to Tmin (finalTemp)
  // > tau <- d^2/2 # larger tau takes longer to get to Tmin
  // tau = (train->NumVariables() * train->NumVariables()) / 2.0; 
  tau = par::ecPrivacyTau;
  summedProbabilities = 0;
  randUniformValue = 0;
  randomForestPredictError = 0;
  trainError = 0;
  holdError = 0;
  trainError = 0;
  // --------------------------------------------------------------------------
  // end of constructor
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy END\n");
  PrintState();
}

EvaporativeCoolingPrivacy::~EvaporativeCoolingPrivacy() {
}

bool EvaporativeCoolingPrivacy::ComputeScores() {
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() START\n");
  PP->printLOG(Timestamp() + "Running temperature schedule until final temperature\n");

  // write statistics for every iteration to file
  string iterationOutputFile = par::output_file_name + ".privacy.iterations.tab";
  PP->printLOG(Timestamp() + "Writing iteration error results to [ " + 
    iterationOutputFile +  " ]\n");
 	ofstream iterationOutputStream(iterationOutputFile);
	if(!iterationOutputStream.is_open()) {
		error("Could not open iteration output file [ " + iterationOutputFile + " ]\n");
	}
  if(UsingSimData()) {
    iterationOutputStream 
            << "Iteration\tTemperature\tKeep\tTrainAcc\tHoldoutAcc\tTestAcc\tLastRemoved\tCorrect" << endl;
  } else {
    iterationOutputStream 
            << "Iteration\tTemperature\tKeep\tTrainAcc\tHoldoutAcc\tTestAcc\tLastRemoved" << endl;
  }

  // initialize all masks to contain all variables
  PP->printLOG(Timestamp() + "Initializing variable masks\n");
  train->MaskIncludeAllAttributes(NUMERIC_TYPE);
  holdout->MaskIncludeAllAttributes(NUMERIC_TYPE);

  // compute train and holdout importance with Relief-F
  PP->printLOG(Timestamp() + "Calculating initial importance scores and delta Q\n");
  this->ComputeImportance();
  PP->printLOG(Timestamp() + "\tscore count: " + int2str(trainImportance.size()) + "\n");
  PP->printLOG(Timestamp() + "\tdelta Q: " + dbl2str(deltaQ) + "\n");

  // main optimization loop
  PP->printLOG(Timestamp() + "Entering EVAPORATIVE COOLING cooling schedule loop\n");
  startTemp = par::ecPrivacyStartTemp;
  finalTemp = par::ecPrivacyFinalTemp;
  currentTemp = startTemp;
  PP->printLOG(Timestamp() + "Starting temperature: " + dbl2str(startTemp) + "\n");
  numInstances = train->NumInstances();
  // TODO: make these two parameters?
  threshold = 4.0 / sqrt(numInstances);
  tolerance = 1.0 / sqrt(numInstances);
  tau = par::ecPrivacyTau;
  removeAttrs.clear();
  keepAttrs = train->GetVariableNames();
  iteration = 1;
  update = 0;
  while((currentTemp > finalTemp) && (train->NumVariables() > 0)) {
    if(iteration > MAX_ITERATIONS) {
      error("EvaporativeCoolingPrivacy::ComputeScores() Maximum iterations reached " 
              + int2str(iteration) + "\n");
    }
    PP->printLOG(Timestamp() + "--------------------------------------------\n");
    PP->printLOG(Timestamp() + "iteration = " + int2str(iteration) + "\n");
    train->PrintMaskStats();
    PP->printLOG(Timestamp() + "Number of mask variables: " + 
      int2str(curVarNames.size()) + "\n");
    // importance
    // recompute p_t, p_h
    PP->printLOG(Timestamp() + "Computing p_t, p_h and delta_t\n");
    this->ComputeAttributeProbabilities();
    // evaporate worst
    uint attributesRemoved = this->EvaporateWorstAttributes(1);
    if(attributesRemoved) {
      PP->printLOG(Timestamp() + "Dataset now has: " + int2str(train->NumVariables()) + "\n");
    } else {
      error("EvaporativeCoolingPrivacy::ComputeScores() could not remove worst");
    }
    if((iteration % updateInterval) == 1) {
      PP->printLOG(Timestamp() + "Calculating importance scores and delta Q\n");
      this->ComputeImportance();
      PP->printLOG(Timestamp() + "\tscore count: " + int2str(trainImportance.size()) + "\n");
      PP->printLOG(Timestamp() + "\tdelta Q: " + dbl2str(deltaQ) + "\n");
      PP->printLOG(Timestamp() + "Running classifiers on train, holdout and test\n");
      this->ComputeBestAttributesErrors();
      PP->printLOG(Timestamp() + "Running temperature update\n");
      this->UpdateTemperature();
      PP->printLOG(Timestamp() + "New temperature: " + dbl2str(currentTemp) + "\n");
      if(par::verbose) this->PrintState();
      //      cout << "Accuracies:" << endl 
      //              << iteration << "\t"
      //              << keepAttrs.size() << "\t"
      //              << (1 - trainError) << "\t" 
      //              << (1 - holdError) << "\t" 
      //              << (1 - testError) << "\t"
      //              << endl;
      // ----------------------------------------------------------------------
      // write iteration update results to iterationOutputFile
      // if simulated data, report the number of correctly detected attributes
      string lastGeneRemoved = removeAttrs[removeAttrs.size()-1];
      if(UsingSimData()) {
        uint numSignalsFound = CurrentNumberCorrect(keepAttrs);
        iterationOutputStream 
                << iteration << "\t"
                << currentTemp << "\t"
                << keepAttrs.size() << "\t"
                << (1 - trainError) << "\t" 
                << (1 - holdError) << "\t" 
                << (1 - testError) << "\t"
                << lastGeneRemoved << "\t"
                << numSignalsFound
                << endl;
      } else {
        iterationOutputStream 
                << iteration << "\t"
                << currentTemp << "\t"
                << keepAttrs.size() << "\t"
                << (1 - trainError) << "\t" 
                << (1 - holdError) << "\t" 
                << (1 - testError) << "\t"
                << lastGeneRemoved
                << endl;
        // << insilico::join(keepAttrs.begin(), keepAttrs.end(), ",") 
      }
      ++update;
    }
    ++iteration;
  }
  // end temperature schedule
  iterationOutputStream.close();
  if(par::verbose) this->PrintState();
  PP->printLOG(Timestamp() + "Exiting EVAPORATIVE COOLING cooling schedule after " + int2str(iteration) + " iterations\n");
  PP->printLOG(Timestamp() + "Final temperature:              " + dbl2str(currentTemp) + "\n");
  PP->printLOG(Timestamp() + "Final kept attributes set size: " + int2str(keepAttrs.size()) + "\n");
  PP->printLOG(Timestamp() + "Final classification errors:    " + dbl2str(trainError) + 
               "\t" + dbl2str(holdError) + "\t" + dbl2str(testError) + "\n");

  // clean up temporary Ranger prediction files
  string removeFile = par::output_file_name + ".forest";
  unlink(removeFile.c_str());
  removeFile = par::output_file_name + ".prediction";
  unlink(removeFile.c_str());

  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() END\n");
  
  return true;
}

uint EvaporativeCoolingPrivacy::CurrentNumberCorrect(vector<string> testSet) {
  uint returnValue = 0;
  for(uint i=0; i < testSet.size(); ++i) {
    string candidate = testSet[i];
    bool found = false;
    for(uint j=0; (j < signalNames.size()) && !found; ++j) {
      if(candidate == signalNames[j]) {
        found = true;
        ++returnValue;
      }
    }
  }
  return returnValue;
}

void EvaporativeCoolingPrivacy::PrintState() {
  cout << "**********************************************************" << endl;
  cout << "PRIVACY EC CURRENT STATE" << endl;
  cout << "Q_EPS:              " << Q_EPS << endl;
  cout << "iteration:          " << iteration << endl;
  cout << "updates:            " << update << endl;
  cout << "max iterations:     " << MAX_ITERATIONS << endl;
  cout << "deltaQ:             " << deltaQ << endl;
  cout << "start temp:         " << startTemp << endl;
  cout << "current temp:       " << currentTemp << endl;
  cout << "final temp:         " << finalTemp << endl;
  cout << "tau:                " << tau << endl;
  cout << "update frequency:   " << updateInterval << endl;
  cout << "num attributes:     " << train->NumVariables() << " (train)" << endl;
  cout << "num instances:      " << train->NumInstances() << " (train)" << endl;
  cout << "data simulated?     " << (UsingSimData()? "true": "false") << endl;
  if(UsingSimData()) {
    cout << "num sim signals:    " << numSignalsInData << endl;
  }
  cout << "last train error:   " << trainError << endl;
  cout << "last holdout error: " << holdError << endl;
  if(test) {
    cout << "last test error:    " << testError << endl;  
  }

  cout << "**********************************************************" << endl;
}

ResultsLists EvaporativeCoolingPrivacy::GetKeptRemoved() {
  ResultsLists returnLists;
  returnLists.first = keepAttrs;
  returnLists.second = removeAttrs;
  
  return returnLists;
}

bool EvaporativeCoolingPrivacy::WriteBestAttributes(string fileSuffix) {
  if(!keepAttrs.size()) {
    PP->printLOG("WARNING EvaporativeCoolingPrivacy::WriteBestAttributes: "
            "No attributes to write\n");
    return false;
  }
  string resultsFilename = fileSuffix + ".selected.attributes.tab";
	PP->printLOG(Timestamp() + "Writing Privacy EC results to: " + resultsFilename + "\n");
  ofstream outFile;
  outFile.open(resultsFilename.c_str());
  if(outFile.bad()) {
    error("ERROR: Could not open scores file " + resultsFilename + " for writing\n");
  }
  for(uint i=0; i < keepAttrs.size(); ++i) {
    outFile << keepAttrs[i] << endl;
  }
  outFile.close();

  return true;
}

pair<uint, double> EvaporativeCoolingPrivacy::CheckDetectedAttributes() {
  pair<uint, double>  retValues;
  uint goodCount = 0;
  for(uint i=0; i < keepAttrs.size(); ++i) {
    string thisAttr = keepAttrs[i];
    if(find(signalNames.begin(), signalNames.end(), thisAttr) != signalNames.end()) {
      ++goodCount;
    }
  }
  retValues.first = goodCount;
  retValues.second = (double) goodCount / numSignalsInData;
  
  return retValues;
}

// ----------------------------------------------------------------------------
// private methods

bool EvaporativeCoolingPrivacy::ComputeImportance() {
  // --------------------------------------------------------------------------
  ReliefF* relief = new ReliefF(train, PP, NUMERIC_ONLY_ANALYSIS);
  // relief->SetNormalize(true);
  if(relief->ComputeAttributeScores()) {
    //relief->PrintScores(std::cout);
    trainImportance.clear();
    trainImportance = relief->GetScores();
  } else {
    error("ReliefF::ComputeAttributeScores() failed on training set\n");
  }
  if(trainImportance.size() != train->NumVariables()) {
    error("EvaporativeCoolingPrivacy::ComputeImportance size mismatch in importance");
  }
  sort(trainImportance.begin(), trainImportance.end(), scoresSortAscByName);
  delete relief;
  // --------------------------------------------------------------------------
  relief = new ReliefF(holdout, PP, NUMERIC_ONLY_ANALYSIS);
  // relief->SetNormalize(true);
  if(relief->ComputeAttributeScores()) {
    holdoutImportance.clear();
    holdoutImportance = relief->GetScores();
  } else {
    error("ReliefF::ComputeAttributeScores() failed on holdout set\n");
  }
  if(holdoutImportance.size() != train->NumVariables()) {
    error("EvaporativeCoolingPrivacy::ComputeImportance size mismatch in importance");
  }
  sort(holdoutImportance.begin(), holdoutImportance.end(), scoresSortAscByName);
  delete relief;
  // --------------------------------------------------------------------------
  PP->printLOG(Timestamp() + "Find delta Q = maximum absolute difference in training and " 
                "holdout importance scores\n");
  AttributeScoresCIt trainIt = trainImportance.begin();
  AttributeScoresCIt holdoIt = holdoutImportance.begin();
  diffScores.clear();
  for(; trainIt != trainImportance.end(); ++trainIt, ++holdoIt) {
    double absDiff = fabs((*trainIt).first - (*holdoIt).first);
    diffScores.push_back(absDiff);
    string varName = (*trainIt).second;
    diffImportance[varName] = absDiff;
  }
  deltaQ = *max_element(diffScores.begin(), diffScores.end());
  PP->printLOG(Timestamp() + "deltaQ: " + dbl2str(deltaQ) + "\n");

  //  cout << endl
  //          << "EvaporativeCoolingPrivacy::ComputeImportance "
  //          << "train importance size: " 
  //          << trainImportance.size() << endl << endl;
  
  return true;
}

bool EvaporativeCoolingPrivacy::ComputeAttributeProbabilities() {
    PP->printLOG(Timestamp() + "Computing probabilities of attributes P(a)\n");
    //    >   PAs <- exp(-q1.scores/(2*delta.q*k*myT))
    //    >   sumPAs <- sum(PAs)
    //    >   scaled.PAs <- PAs/sum(PAs)
    //    >   cum.scaled.PAs <- cumsum(scaled.PAs)
    PP->printLOG(Timestamp() + "\tdiff)\n");
    diff.clear();
    for(uint i=0; i < trainImportance.size(); ++i) {
      string key = trainImportance[i].second;
      if(diffImportance[key] > 0.01) {
        diff.push_back(diffImportance[key]);
      } else {
        diff.push_back(0.01);
      }
    }
    PP->printLOG(Timestamp() + "\tP(a)) and sum probs\n");
    // compute attribute probabilities and sum of probabilities
    attributeProbabilty.clear();
    summedProbabilities = 0;
    AttributeScoresCIt trainIt = trainImportance.begin();
    for(uint i=0; trainIt != trainImportance.end(); ++trainIt, ++i) {
      double q1Score = (*trainIt).first;
      double prob = exp(-q1Score / (2.0 * diff[i] * currentTemp));
      if(std::isinf(prob)) {
        PP->printLOG("WARNING: infinity detected; setting to MAX_DOUBLE\n");
        prob = numeric_limits<double>::max();
      }
      if(par::algorithm_verbose) {
        cout << prob << "\t" << q1Score << "\t" << diff[i] << "\t" 
                << currentTemp << endl;
      }
      attributeProbabilty.push_back(prob);
      summedProbabilities += prob;
    }
    PP->printLOG(Timestamp() + "\tscale probs\n");
    // normalize/scale to (0, 1)
    scaledProbabilities.clear();
    for(uint i=0; i < attributeProbabilty.size(); ++i) {
      scaledProbabilities.push_back(attributeProbabilty[i] / summedProbabilities);
    }
    PP->printLOG(Timestamp() + "\tcumsum\n");
    // cumulative probabilities
    cummulativeProbabilities.clear();
    double runningProbability = 0;
    for(uint i=0; i < scaledProbabilities.size(); ++i) {
      runningProbability += scaledProbabilities[i];
      //      cout << i << "\t"
      //              << curVarNames[i] << "\t"
      //              << scaledProbabilities[i] << "\t" 
      //              << runningProbability << endl;
      cummulativeProbabilities.push_back(runningProbability);
    }

    // DEBUG ASSERT
    if(runningProbability < 0.999999) {
      error("EvaporativeCoolingPrivacy::ComputeAttributeProbabilities "
            "Bill your stats are whack bud!\n");
    }
    
    return true;
}

bool EvaporativeCoolingPrivacy::GenerateRandomUniformProbabilities() {
  // prob.rands <- runif(1, min = 0, max = 1)
  PP->printLOG(Timestamp() + "Generating uniform probabilities in (0, 1)\n");
  uniform_real_distribution<double> runif(0, 1);
  // generate random uniform probabilities
  randUniformProbs.clear();
  for(uint prIdx=0; prIdx < trainImportance.size(); ++prIdx) {
    randUniformProbs.push_back(runif(engine));
  }
  // this sort is IMPORTANT! bcw 12/9/16
  sort(randUniformProbs.begin(), randUniformProbs.end());
  
  // TODO: bcw - 12/31/16 - use one random value not a vector?
  // prob.rands <- runif(1, min = 0, max = 1)
  randUniformValue = runif(engine);
          
  return true;
}

// remove worst attributes: has 's' in the method name
uint EvaporativeCoolingPrivacy::EvaporateWorstAttributes(uint numToRemove) {
  PP->printLOG(Timestamp() + "Evaporating the worst " + int2str(numToRemove) + " attributes\n");
  // find which attributes to remove and those to keep
  uniform_real_distribution<double> runif(0, 1);
  randUniformValue = runif(engine);
  bool found = false;
  string thisVar;
  if(par::verbose) PP->printLOG(Timestamp() + "Scanning probabilities\n");
  possiblyRemove.clear();
  for(uint curVarIdx=0; curVarIdx < curVarNames.size(); ++curVarIdx) {
    thisVar = curVarNames[curVarIdx];
    if(randUniformValue < cummulativeProbabilities[curVarIdx]) {
      possiblyRemove.push_back(thisVar);
      found = true;
    } 
  }
  if(!found) {
    PP->printLOG("WARNING: No variables found to remove\n");
    return 0;
  }
  // remove numToRemove of the toRemove variables
  uint loopToRemove = 0;
  if(found) {
    PP->printLOG(Timestamp() + "Found " + int2str(possiblyRemove.size()) + " candidates to remove\n");
    loopToRemove = numToRemove;
    if(possiblyRemove.size() < loopToRemove) {
      PP->printLOG("WARNING: Number to remove is larger than to remove candidates\n");
      loopToRemove = possiblyRemove.size();
    }
    for(uint i=0; i < loopToRemove; ++i) {
      string removeVar = possiblyRemove[i];
      PP->printLOG(Timestamp() + "Evaporating/removing variable: " + removeVar + "\n");
      train->MaskRemoveVariable(removeVar);
      holdout->MaskRemoveVariable(removeVar);
      if(test) {
        test->MaskRemoveVariable(removeVar);
      }
      removeAttrs.push_back(removeVar);
      // remove the importance scores for the variables removed
      this->RemoveImportanceScore(removeVar);
    }
    keepAttrs = train->MaskGetAllVariableNames();
  } else {
    PP->printLOG(Timestamp() + "WARNING: Could not remove any variables in EvaporateWorstAttributes\n");
    return 0;
  }
    
  // update the algorithm tracking variables
  curVarNames = train->MaskGetAllVariableNames();
  curVarMap = train->MaskGetAttributeMask(NUMERIC_TYPE);

  return loopToRemove;
}

bool EvaporativeCoolingPrivacy::RemoveImportanceScore(std::string varToRemove) {
  bool found = false;
  AttributeScoresIt citTrain = trainImportance.begin();
  AttributeScoresIt citHoldout = holdoutImportance.begin();
  while(!found && 
        (citTrain != trainImportance.end())&& 
        (citHoldout != holdoutImportance.end())) {
    if(((*citTrain).second == varToRemove) &&
       ((*citHoldout).second == varToRemove)) {
      found = true;
      trainImportance.erase(citTrain);
      holdoutImportance.erase(citHoldout);
    }
    ++citTrain;
    ++citHoldout;
  }

  return found;
}

// remove *the* worst attribute: no 's' in the method name
//bool EvaporativeCoolingPrivacy::EvaporateWorstAttribute() {
//  if(par::verbose) 
//    PP->printLOG(Timestamp() + "Evaporating single worst attribute\n");
//  sort(trainImportance.begin(), trainImportance.end(), scoresSortAsc);
//  pair<double, string> worstAttr = *(trainImportance.begin());
//  string thisVar = worstAttr.second;
//  double thisVarValue = worstAttr.first;
//  sort(trainImportance.begin(), trainImportance.end(), scoresSortAscByName);
//  if(par::verbose) {
//    PP->printLOG(Timestamp() + "Worst attribute [" + thisVar + "]\n");
//    PP->printLOG(Timestamp() + thisVar + " (" + int2str(thisVarValue) + ")\n");
//  }
//  // adjust masks
//  train->MaskRemoveVariableType(thisVar, NUMERIC_TYPE);
//  holdout->MaskRemoveVariableType(thisVar, NUMERIC_TYPE);
//  test->MaskRemoveVariableType(thisVar, NUMERIC_TYPE);
//  // keep/remove accounting
//  removeAttrs.push_back(thisVar);
//  keepAttrs = train->MaskGetAllVariableNames();
//
//  // update the algorithm tracking variables
//  curVarNames = train->MaskGetAllVariableNames();
//  curVarMap = train->MaskGetAttributeMask(NUMERIC_TYPE);
//
//  return true;
//}

bool EvaporativeCoolingPrivacy::ComputeBestAttributesErrors() {
  PP->printLOG(Timestamp() + "Compute errors for the best attributes\n");
  // get training, holdout and prediction errors with best attributes (keepAttrs)
  trainError = ClassifyAttributeSet(keepAttrs, TRAIN);
  trainErrors.push_back(trainError);
  holdError = ClassifyAttributeSet(keepAttrs, HOLDOUT);
  holdoutErrors.push_back(holdError);
  uniform_real_distribution<double> runif(0, tolerance);
  double lilBit = runif(engine);
  if(fabs(trainError - holdError) < (threshold + lilBit)) {
    holdError = trainError;
  } else {
    PP->printLOG(Timestamp() + "adjusting small difference in holdout: " + dbl2str(lilBit) + "\n");
    holdError = holdError + lilBit;
  }
  testError = 1.0;
  if(test) {
    testError = ClassifyAttributeSet(keepAttrs, TEST);
  }
  testErrors.push_back(testError);

  PP->printLOG(Timestamp() + "* Error\n");
  PP->printLOG(Timestamp() + "* train:   " + dbl2str(trainError) + "\n" +
               Timestamp() + "* holdout: " + dbl2str(holdError) + "\n");
  if(test) {
    PP->printLOG(Timestamp() + "* test:    " + dbl2str(testError) + "\n");
  }
  PP->printLOG(Timestamp() + "* Accuracy\n");
  PP->printLOG(Timestamp() + "* train:   " + dbl2str(1 - trainError) + "\n" +
               Timestamp() + "* holdout: " + dbl2str(1 - holdError) + "\n");
  if(test) {
    PP->printLOG(Timestamp() + "* test:    " + dbl2str(1 - testError) + "\n");
  }
}

double 
EvaporativeCoolingPrivacy::ClassifyAttributeSet(vector<string> attrs, 
                                                DATASET_TYPE dataType) {
  double retError = 1.0;
  RandomForest* randomForest = 0;
  string outfileTrain, outfileTest, outfileHoldout;
  switch(dataType) {
    case TRAIN:
      PP->printLOG(Timestamp() + "Classify best attributes for TRAINING data\n");
      outfileTrain = par::output_file_name + ".train.tmp";
      train->WriteNewDataset(outfileTrain, TAB_DELIMITED_DATASET);
      par::ecPrivacyTrainFile = outfileTrain;
      randomForest = new RandomForest(train, par::ecPrivacyTrainFile, attrs, false);
      randomForest->ComputeScores();
      retError = randomForest->GetClassificationError();
      randomForest->SaveForest();
      unlink(outfileTrain.c_str());
      break;
    case HOLDOUT:
      PP->printLOG(Timestamp() + "Classify best attributes for HOLDOUT data\n");
      outfileHoldout = par::output_file_name + ".holdout.tmp";
      train->WriteNewDataset(outfileHoldout, TAB_DELIMITED_DATASET);
      par::ecPrivacyTrainFile = outfileHoldout;
      randomForest = new RandomForest(holdout, par::ecPrivacyHoldoutFile, attrs, false);
      randomForest->ComputeScores();
      retError = randomForest->GetClassificationError();
      unlink(outfileHoldout.c_str());
      break;
    case TEST:
      PP->printLOG(Timestamp() + "Predicting best attributes for TESTING data\n");
      outfileTest = par::output_file_name + ".test.tmp";
      train->WriteNewDataset(outfileTest, TAB_DELIMITED_DATASET);
      par::ecPrivacyTrainFile = outfileTest;
      randomForest = new RandomForest(test, par::ecPrivacyTestFile, attrs, false);
      retError = randomForest->Predict();
      unlink(outfileTest.c_str());
      break;
    default:
      error("EvaporativeCoolingPrivacy::ClassifyAttributeSet Dataset type not recognized");
  }
  if(randomForest) {
    delete randomForest;
  }
  
  return retError;
}

bool EvaporativeCoolingPrivacy::UpdateTemperature() {
  PP->printLOG(Timestamp() + "Updating temperature\n");
  // update 11/1/16 from Trang:
  // myT <- myT*exp(-i/tau)
  currentTemp *= exp(-1 / tau);
  
  return true;
}

bool EvaporativeCoolingPrivacy::ComputeInverseImportance() {
  PP->printLOG(Timestamp() + "Computing inverse importance scores for training and holdout sets\n");

  ReliefF relief(train, PP, NUMERIC_ONLY_ANALYSIS);
  relief.ComputeAttributeScores();
  trainInvImportance = relief.GetScores();
  sort(trainInvImportance.begin(), trainInvImportance.end(), scoresSortAscByName);
  vector<double> trainRelief;
  vector<string> trainVars;
  for_each(trainInvImportance.begin(), trainInvImportance.end(),
           [&trainRelief, &trainVars](const map<double, string>::value_type& p) 
           { trainRelief.push_back(p.first); trainVars.push_back(p.second); });
  double minTrain = *(min_element(trainRelief.begin(), trainRelief.end()));
  trainInvImportance.clear();
  for(uint idx=0; idx < trainRelief.size(); ++idx) {
    double x = trainRelief[idx];
    trainInvImportance.push_back(make_pair(1 / (x - minTrain + Q_EPS), 
            trainVars[idx]));
  }

  ReliefF rf(holdout, PP, NUMERIC_ONLY_ANALYSIS);
  rf.ComputeAttributeScores();
  holdoutInvImportance = rf.GetScores();
  sort(holdoutInvImportance.begin(), holdoutInvImportance.end(), scoresSortAscByName);
  vector<double> holdoRelief;
  vector<string> holdoVars;
  for_each(holdoutInvImportance.begin(), holdoutInvImportance.end(),
           [&holdoRelief, &holdoVars](const map<double, string>::value_type& p) 
           { holdoRelief.push_back(p.first); holdoVars.push_back(p.second); });
  double min_holdo = *min_element(holdoRelief.begin(), holdoRelief.end());
  holdoutInvImportance.clear();
  for(uint idx=0; idx < holdoRelief.size(); ++idx) {
    double x = holdoRelief[idx];
    holdoutInvImportance.push_back(make_pair(1 / (x - min_holdo + Q_EPS), 
            holdoVars[idx]));
  }
           
  return true;
}
