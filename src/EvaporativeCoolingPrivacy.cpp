/* 
 * File:   EvaporativeCoolingPrivacy.cpp
 * Author: bwhite
 * 
 * Created on October 18, 2016, 10:57 PM
 */

#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <armadillo>
#include <random>
#include <cmath>

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

using namespace std;
using namespace arma;

// ----------------------------------------------------------------------------
// public methods

EvaporativeCoolingPrivacy::EvaporativeCoolingPrivacy(Dataset* trainset, 
        Dataset* holdoset, Dataset* testset, Plink* plinkPtr) {
  // static class variables
  Q_EPS = 0.005;
  MAX_ITERATIONS = 1e6;
  // pointer to a PLINK environment
  PP = plinkPtr;
  train = trainset;
  holdout = holdoset;
  test = testset;
  numInstances = trainset->NumInstances();
  curVarNames = trainset->MaskGetAllVariableNames();
  curVarMap = trainset->MaskGetAttributeMask(NUMERIC_TYPE);

  deltaQ = 0;
  threshold = 4.0 / sqrt(numInstances);
  tolerance = 1.0 / sqrt(numInstances);
  
  numToRemovePerIteration = par::ecPrivacyRemovePerIteration;
  startTemp = par::ecPrivacyStartTemp;
  currentTemp = startTemp;
  finalTemp = par::ecPrivacyFinalTemp;
  // Trang: larger tau takes longer to get to Tmin (finalTemp)
  // > tau <- d^2/2 # larger tau takes longer to get to Tmin
  // tau = (train->NumVariables() * train->NumVariables()) / 2.0; 
  tau = par::ecPrivacyTau;
  minRemainAttributes = par::ecPrivacyMinVars;
  updateInterval = par::ecPrivacyUpdateFrequency;
  
  summedProbabilities = 0;
  randomForestPredictError = 0;
  trainError = 0;
  holdError = 0;
  trainError = 0;
  
  iteration = 0;
  this->PrintState();
}

EvaporativeCoolingPrivacy::~EvaporativeCoolingPrivacy() {
}

bool EvaporativeCoolingPrivacy::ComputeScores() {
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() START\n");
  PP->printLOG(Timestamp() + "Running temperature schedule until final temperature\n");

  // write statistics for every iteration
  string iterationOutputFile = par::output_file_name + ".privacy.iterations.tab";
  PP->printLOG(Timestamp() + "Writing iteration error results to [ " + 
    iterationOutputFile +  " ]\n");
 	ofstream iterationOutputStream(iterationOutputFile);
	if(!iterationOutputStream.is_open()) {
		error("Could not open iteration output file [ " + iterationOutputFile + " ]\n");
	}
  iterationOutputStream << "Iteration\tTemperature\tKeep\tRemove\tTrain\tHoldout\tTest" << endl;

  // initialize
  PP->printLOG(Timestamp() + "Initializing variable masks\n");
  train->MaskIncludeAllAttributes(NUMERIC_TYPE);
  holdout->MaskIncludeAllAttributes(NUMERIC_TYPE);
  origVarNames = train->MaskGetAllVariableNames();
  origVarMap = train->MaskGetAttributeMask(NUMERIC_TYPE);

  // calculate importance by running Relief-F on train and holdout sets
  PP->printLOG(Timestamp() + "Calculating importance scores and delta Q\n");
  ComputeImportance();
    
  // main optimization loop
  PP->printLOG(Timestamp() + "Entering EVAPORATIVE COOLING cooling schedule\n");
  startTemp = par::ecPrivacyStartTemp;
  finalTemp = par::ecPrivacyFinalTemp;
  currentTemp = startTemp;
  PP->printLOG(Timestamp() + "Starting temperature: " + dbl2str(startTemp) + "\n");
  tau = par::ecPrivacyTau;
  iteration = 0;
  while((currentTemp > finalTemp) && (train->NumVariables() >= minRemainAttributes)) {
    ++iteration;
    if(iteration > MAX_ITERATIONS) {
      error("EvaporativeCoolingPrivacy::ComputeScores() Maximum ticks reached " 
              + int2str(iteration) + "\n");
    }
    PP->printLOG(Timestamp() + "--------------------------------------------\n");
    PP->printLOG(Timestamp() + "iteration = " + int2str(iteration) + "\n");
    // PrintState();
    curVarNames = train->MaskGetAllVariableNames();
    curVarMap = train->MaskGetAttributeMask(NUMERIC_TYPE);
    // recompute p_t, p_h and delta_t
    PP->printLOG(Timestamp() + "Computing p_t, p_h and delta_t\n");
    ComputeAttributeProbabilities();
    GenerateRandomUniformProbabilities();
    EvaporateWorstAttributes(numToRemovePerIteration);
    PP->printLOG(Timestamp() + "Keeping : " + int2str(keepAttrs.size()) + "\n");
    PP->printLOG(Timestamp() + "Removing: " + int2str(removeAttrs.size()) + "\n");
    PP->printLOG(Timestamp() + "Dataset : " + int2str(train->NumVariables()) + "\n");
    if((iteration % 50) == 1) {
      PP->printLOG(Timestamp() + "***** DEBUG: start mod 50 UPDATE\n");
      // calculate importance by running Relief-F on train and holdout sets
      PP->printLOG(Timestamp() + "Calculating importance scores and delta Q\n");
      ComputeImportance();
      PP->printLOG(Timestamp() + "Running classifiers on train, holdout and test\n");
      ComputeBestAttributesErrors();
      PP->printLOG(Timestamp() + "Updating temperature\n");
      UpdateTemperature();
      PP->printLOG(Timestamp() + "New temperature: " + dbl2str(currentTemp) + "\n");
      PP->printLOG(Timestamp() + "***** DEBUG: end mod 50 UPDATE\n");
    }
    // write iteration results to iterationOutputFile
    iterationOutputStream 
            << iteration << "\t"
            << currentTemp << "\t"
            << keepAttrs.size() << "\t"
            << removeAttrs.size() << "\t"
            << trainError << "\t" 
            << holdError << "\t" 
            << testError << endl;
  } // end temperature schedule
  iterationOutputStream.close();

  PP->printLOG(Timestamp() + "Exiting EVAPORATIVE COOLING cooling schedule after " + int2str(iteration) + " iterations\n");
  PP->printLOG(Timestamp() + "Final temperature:              " + dbl2str(currentTemp) + "\n");
  PP->printLOG(Timestamp() + "Final kept attributes set size: " + int2str(keepAttrs.size()) + "\n");
  PP->printLOG(Timestamp() + "Final classification errors: " + dbl2str(trainError) + 
               "\t" + dbl2str(holdError) + "\t" + dbl2str(testError) + "\n");

  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() END\n");
  
  return true;
}

void EvaporativeCoolingPrivacy::PrintState() {
  cout << "**********************************************************" << endl;
  cout << "CURRENT STATE" << endl;
  //cout << "Q_EPS:              " << Q_EPS << endl;
  cout << "iteration:          " << iteration << endl;
  cout << "deltaQ:             " << deltaQ << endl;
  cout << "threshold:          " << threshold << endl;
  cout << "tolerance:          " << tolerance << endl;
  cout << "remove per:         " << numToRemovePerIteration << endl;
  cout << "start temp:         " << startTemp << endl;
  cout << "current temp:       " << currentTemp << endl;
  cout << "final temp:         " << finalTemp << endl;
  cout << "tau:                " << tau << endl;
  cout << "num attributes:     " << train->NumVariables() << endl;
  cout << "num instances:      " << numInstances << endl;
  cout << "last train error:   " << trainError << endl;
  cout << "last holdout error: " << holdError << endl;
  cout << "last test error:    " << testError << endl;
  cout << "**********************************************************" << endl;
}

// ----------------------------------------------------------------------------
// private methods

bool EvaporativeCoolingPrivacy::ComputeImportance() {
  // --------------------------------------------------------------------------
  ReliefF* relief = new ReliefF(train, PP, NUMERIC_ONLY_ANALYSIS);
  // relief->SetNormalize(true);
  if(relief->ComputeAttributeScores()) {
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
  for(; trainIt != trainImportance.end(); ++trainIt, ++holdoIt) {
    double absDiff = fabs((*trainIt).first - (*holdoIt).first);
    diffScores.push_back(absDiff);
    string varName = (*trainIt).second;
    diffImportance[varName] = absDiff;
  }
  deltaQ = *max_element(diffScores.begin(), diffScores.end());

  // increasing importance score index to remove next
  curDeleteImportanceIndex = 0;
  
  return true;
}

bool EvaporativeCoolingPrivacy::ComputeAttributeProbabilities() {
    PP->printLOG(Timestamp() + "Computing probabilities of attributes P(a)\n");
    //    >   PAs <- exp(-q1.scores/(2*delta.q*k*myT))
    //    >   sumPAs <- sum(PAs)
    //    >   scaled.PAs <- PAs/sum(PAs)
    //    >   cum.scaled.PAs <- cumsum(scaled.PAs)
    diff.resize(trainImportance.size());
    for(uint i=0; i < trainImportance.size(); ++i) {
      string key = trainImportance[i].second;
      if(diffImportance[key] > 0.001) {
        diff[i] = diffImportance[key];
      } else {
        diff[i] = 0.001;
      }
    }
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
    // normalize/scale to (0, 1)
    scaledProbabilities.clear();
    for(uint i=0; i < attributeProbabilty.size(); ++i) {
      scaledProbabilities.push_back(attributeProbabilty[i] / summedProbabilities);
    }
    // cumulative probabilities
    cummulativeProbabilities.clear();
    double runningProbability = 0;
    for(uint i=0; i < scaledProbabilities.size(); ++i) {
      runningProbability += scaledProbabilities[i];
      cout << i << "\t" << scaledProbabilities[i] << "\t" << runningProbability << endl;
      cummulativeProbabilities.push_back(runningProbability);
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
  
  return true;
}

// remove worst attributes: has 's' in the method name
bool EvaporativeCoolingPrivacy::EvaporateWorstAttributes(uint numToRemove) {
  PP->printLOG(Timestamp() + "Evaporating the worst attributes\n");
  // find which attributes to remove and those to keep
  removeAttrs.clear();
  keepAttrs.clear();
  uint removedSoFar = 0;
  vector<string> toRemove;
  for(uint pIdx=0; pIdx < randUniformProbs.size(); ++pIdx) {
    string thisVar = curVarNames[pIdx];
    uint thisVarIdx = curVarMap[thisVar];
    if(par::verbose) {
      PP->printLOG(Timestamp() + thisVar + " (" + int2str(thisVarIdx) + ")\n");
      PP->printLOG(Timestamp() + dbl2str(randUniformProbs[pIdx]) + 
        " < "  + dbl2str(cummulativeProbabilities[pIdx]) + "?");
    }
    if((removedSoFar < numToRemove) && 
       (randUniformProbs[pIdx] < cummulativeProbabilities[pIdx])) {
      if(par::verbose) PP->printLOG(Timestamp() + " => remove\n");
      toRemove.push_back(thisVar);
      removeAttrs.push_back(thisVar);
      ++removedSoFar;
    } else {
      if(par::verbose) PP->printLOG(Timestamp() + " => keep\n");
      keepAttrs.push_back(thisVar);
    }
  }
  // remove the toRemove variables
  if(toRemove.size()) {
    for(uint i=0; i < toRemove.size(); ++i) {
      string thisVar = toRemove[i];
      train->MaskRemoveVariable(thisVar);
      holdout->MaskRemoveVariable(thisVar);
      test->MaskRemoveVariable(thisVar);
    }
  } else {
    PP->printLOG("WARNING: Could not remove any variables in EvaporateWorstAttributes\n");
    return false;
  }

  return true;
}

// remove *the* worst attribute: no 's' in the method name
bool EvaporativeCoolingPrivacy::EvaporateWorstAttribute() {
  //    >   num.remv <- 1 # only remove 1 attribute  
  //    >   remv.atts <- kept.atts[prob.rands < cum.scaled.PAs][1] 
  if(par::verbose) 
    PP->printLOG(Timestamp() + "Evaporating single worst attribute\n");
  sort(trainImportance.begin(), trainImportance.end(), scoresSortAsc);
  pair<double, string> worstAttr = *(trainImportance.begin() + curDeleteImportanceIndex);
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
  keepAttrs = train->GetVariableNames();

  // number of attributes removed after last modulus update step
  ++curDeleteImportanceIndex;
  
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
  if(fabs(trainError - holdError) < threshold + lilBit) {
    holdError = trainError;
  } else {
    PP->printLOG(Timestamp() + "adjusting small difference in holdout: " + dbl2str(lilBit) + "\n");
    holdError = holdError + lilBit;
  }
  testError = ClassifyAttributeSet(keepAttrs, TEST);
  testErrors.push_back(testError);
  PP->printLOG(Timestamp() + "* train:   " + dbl2str(trainError) + "\n" +
               Timestamp() + "* holdout: " + dbl2str(holdError) + "\n" +
               Timestamp() + "* test:    " + dbl2str(testError) + "\n");
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
      error("EvaporativeCoolingPrivacy::ClassifyAttributeSet Dataset type no recognized");
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
  currentTemp = currentTemp * exp(-1 / tau);
  
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
