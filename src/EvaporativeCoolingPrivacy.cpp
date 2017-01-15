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
        Dataset* holdoset, Dataset* testset, Plink* plinkPtr) {
  // pointer to a PLINK environment
  PP = plinkPtr;
  Q_EPS = 0.005;
  MAX_ITERATIONS = 1e6;
  minRemainAttributes = par::ecPrivacyMinVars;
  updateInterval = par::ecPrivacyUpdateFrequency;
  iteration = 0;
  update = 0;
  train = trainset;
  holdout = holdoset;
  test = testset;
  numInstances = trainset->NumInstances();
  numSignalsInData = static_cast<uint>(train->NumVariables() * 
          par::ecPrivacyPercentSignal);
  curVarNames = trainset->MaskGetAllVariableNames();
  curVarMap = trainset->MaskGetAttributeMask(NUMERIC_TYPE);
  for(uint i=0; i < numSignalsInData; ++i) {
    signalNames.push_back(curVarNames[i]);
  }
  deltaQ = 0;
  threshold = 4.0 / sqrt(numInstances);
  tolerance = 1.0 / sqrt(numInstances);
  startTemp = par::ecPrivacyStartTemp;
  currentTemp = startTemp;
  finalTemp = par::ecPrivacyFinalTemp;
  tau = par::ecPrivacyTau;
  
  numToRemovePerIteration = par::ecPrivacyRemovePerIteration;
  // Trang: larger tau takes longer to get to Tmin (finalTemp)
  // > tau <- d^2/2 # larger tau takes longer to get to Tmin
  // tau = (train->NumVariables() * train->NumVariables()) / 2.0; 
  summedProbabilities = 0;
  randUniformValue = 0;
  randomForestPredictError = 0;
  trainError = 0;
  holdError = 0;
  trainError = 0;

  this->PrintState();
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
  iterationOutputStream << "Iteration\tUpdate\tTemperature\tKeep\tRemove\tTrainAcc\tHoldoutAcc\tTestAcc\tCorrect\tRemoved\tRemain" << endl;

  // initialize all masks to contain all variables
  PP->printLOG(Timestamp() + "Initializing variable masks\n");
  train->MaskIncludeAllAttributes(NUMERIC_TYPE);
  holdout->MaskIncludeAllAttributes(NUMERIC_TYPE);

  // compute train and holdout importance with Relief-F
  PP->printLOG(Timestamp() + "Calculating initial importance scores and delta Q\n");
  this->ComputeImportance();
  PP->printLOG(Timestamp() + "\tscore count: " + int2str(trainImportance.size()) + "\n");
  PP->printLOG(Timestamp() + "\tdelta Q: " + dbl2str(deltaQ) + "\n");

  // initialize loop control variables
  iteration = 1;
  uint tail1 = 0;
  uint tail2 = train->NumVariables();
  vector<uint> prevNumInUpdate;
  prevNumInUpdate.push_back(tail1);
  prevNumInUpdate.push_back(tail2);
  if(par::verbose) {
    cout << endl << "i: " << iteration << " prevNumInUpdate tail1: " 
         << tail1 << "\t" << "tail2: " << tail2 << endl << endl;
  }

  // main optimization loop
  PP->printLOG(Timestamp() + "Entering EVAPORATIVE COOLING cooling schedule loop\n");
  startTemp = par::ecPrivacyStartTemp;
  finalTemp = par::ecPrivacyFinalTemp;
  currentTemp = startTemp;
  PP->printLOG(Timestamp() + "Starting temperature: " + dbl2str(startTemp) + "\n");
  threshold = 4.0 / sqrt(numInstances);
  tolerance = 1.0 / sqrt(numInstances);
  tau = par::ecPrivacyTau;
  removeAttrs.clear();
  keepAttrs = train->GetVariableNames();
  while((currentTemp > finalTemp) && (tail1 != tail2) && 
          (train->NumVariables() >= minRemainAttributes)) {
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
    uint numRemoved = this->EvaporateWorstAttributes(numToRemovePerIteration);
    PP->printLOG(Timestamp() + "Removed        : " + int2str(numRemoved) + "\n");
    PP->printLOG(Timestamp() + "Dataset now has: " + int2str(train->NumVariables()) + "\n");
    // if we need to update temperature
    if((iteration % updateInterval) == 1) {
      if(numRemoved) {
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
      }
      prevNumInUpdate.push_back(train->NumVariables());
      //      cout << "Accuracies:" << endl 
      //              << iteration << "\t"
      //              << keepAttrs.size() << "\t"
      //              << (1 - trainError) << "\t" 
      //              << (1 - holdError) << "\t" 
      //              << (1 - testError) << "\t"
      //              << endl;
      ++update;
      // ----------------------------------------------------------------------
      // write iteration update results to iterationOutputFile
      iterationOutputStream 
              << iteration << "\t"
              << update << "\t"
              << currentTemp << "\t"
              << keepAttrs.size() << "\t"
              << updateInterval << "\t"
              << (1 - trainError) << "\t" 
              << (1 - holdError) << "\t" 
              << (1 - testError) << "\t"
              << this->CurrentNumberCorrect(keepAttrs) << "\t"
              << removeAttrs[0] << "\t"
              << insilico::join(keepAttrs.begin(), keepAttrs.end(), ",") 
              << endl;
      // is this number of attributes left same as last update? 'while' above cond
      tail1 = prevNumInUpdate[prevNumInUpdate.size() - 2];
      tail2 = prevNumInUpdate[prevNumInUpdate.size() - 1];
      //    cout << endl << "Update History i: " << iteration << " prevNumInUpdate tail1: " 
      //         << tail1 << "\t" << "tail2: " << tail2 << endl << endl;
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

  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() END\n");

  return true;
}

uint EvaporativeCoolingPrivacy::CurrentNumberCorrect(vector<string> testSet) {
  uint returnValue = 0;
  for(auto i=0; i < testSet.size(); ++i) {
    string candidate = testSet[i];
    bool found = false;
    for(auto j=0; (j < signalNames.size()) && !found; ++j) {
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
  cout << "deltaQ:             " << deltaQ << endl;
  cout << "threshold:          " << threshold << endl;
  cout << "tolerance:          " << tolerance << endl;
  cout << "remove per:         " << numToRemovePerIteration << endl;
  cout << "start temp:         " << startTemp << endl;
  cout << "current temp:       " << currentTemp << endl;
  cout << "final temp:         " << finalTemp << endl;
  cout << "tau:                " << tau << endl;
  cout << "update frequency:   " << updateInterval << endl;
  cout << "num attributes:     " << train->NumVariables() << endl;
  cout << "num instances:      " << numInstances << endl;
  cout << "num sim signals:    " << numSignalsInData << endl;
  cout << "last train error:   " << trainError << endl;
  cout << "last holdout error: " << holdError << endl;
  cout << "last test error:    " << testError << endl;
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
      if(diffImportance[key] > 0.001) {
        diff.push_back(diffImportance[key]);
      } else {
        diff.push_back(0.001);
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
      error("Bill your stats are whack bud!\n");
    }
    
    return true;
}

bool EvaporativeCoolingPrivacy::GenerateRandomUniformProbabilities() {
  //    >   prob.rands <- runif(1, min = 0, max = 1)
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
  //    >   prob.rands <- runif(1, min = 0, max = 1)
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
  // removeAttrs.clear();
  // keepAttrs.clear();
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
      test->MaskRemoveVariable(removeVar);
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
  AttributeScoresCIt citTrain = trainImportance.begin();
  AttributeScoresCIt citHoldout = holdoutImportance.begin();
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
bool EvaporativeCoolingPrivacy::EvaporateWorstAttribute() {
  //    >   num.remv <- 1 # only remove 1 attribute  
  //    >   remv.atts <- kept.atts[prob.rands < cum.scaled.PAs][1] 
  if(par::verbose) 
    PP->printLOG(Timestamp() + "Evaporating single worst attribute\n");
  sort(trainImportance.begin(), trainImportance.end(), scoresSortAsc);
  pair<double, string> worstAttr = *(trainImportance.begin());
  string thisVar = worstAttr.second;
  double thisVarValue = worstAttr.first;
  sort(trainImportance.begin(), trainImportance.end(), scoresSortAscByName);
  if(par::verbose) {
    PP->printLOG(Timestamp() + "Worst attribute [" + thisVar + "]\n");
    PP->printLOG(Timestamp() + thisVar + " (" + int2str(thisVarValue) + ")\n");
  }
  // adjust masks
  train->MaskRemoveVariableType(thisVar, NUMERIC_TYPE);
  holdout->MaskRemoveVariableType(thisVar, NUMERIC_TYPE);
  test->MaskRemoveVariableType(thisVar, NUMERIC_TYPE);
  // keep/remove accounting
  removeAttrs.push_back(thisVar);
  keepAttrs = train->MaskGetAllVariableNames();

  return true;
}

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
  testError = ClassifyAttributeSet(keepAttrs, TEST);
  testErrors.push_back(testError);
  PP->printLOG(Timestamp() + "* Error\n");
  PP->printLOG(Timestamp() + "* train:   " + dbl2str(trainError) + "\n" +
               Timestamp() + "* holdout: " + dbl2str(holdError) + "\n" +
               Timestamp() + "* test:    " + dbl2str(testError) + "\n");
  PP->printLOG(Timestamp() + "* Accuracy\n");
  PP->printLOG(Timestamp() + "* train:   " + dbl2str(1-trainError) + "\n" +
               Timestamp() + "* holdout: " + dbl2str(1-holdError) + "\n" +
               Timestamp() + "* test:    " + dbl2str(1-testError) + "\n");
}

double 
EvaporativeCoolingPrivacy::ClassifyAttributeSet(vector<string> attrs, 
                                                DATASET_TYPE dataType) {
  double retError = 1.0;
  RandomForest* randomForest = 0;
  switch(dataType) {
    case TRAIN:
      PP->printLOG(Timestamp() + "Classify best attributes for TRAINING data\n");
      randomForest = new RandomForest(train, par::ecPrivacyTrainFile, attrs, false);
      randomForest->ComputeScores();
      retError = randomForest->GetClassificationError();
      randomForest->SaveForest();
      break;
    case HOLDOUT:
      PP->printLOG(Timestamp() + "Classify best attributes for HOLDOUT data\n");
      randomForest = new RandomForest(holdout, par::ecPrivacyHoldoutFile, attrs, false);
      randomForest->ComputeScores();
      retError = randomForest->GetClassificationError();
      break;
    case TEST:
      PP->printLOG(Timestamp() + "Predicting best attributes for TESTING data\n");
      randomForest = new RandomForest(test, par::ecPrivacyTestFile, attrs, true);
      retError = randomForest->Predict();
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
