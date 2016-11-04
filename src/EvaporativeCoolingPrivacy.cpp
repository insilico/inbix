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

#include "plink.h"
#include "options.h"
#include "helper.h"

#include "Insilico.h"
#include "Dataset.h"

#include "EvaporativeCoolingPrivacy.h"
#include "ReliefF.h"
#include "RandomForest.h"
#include "ArmadilloFuncs.h"

using namespace std;

EvaporativeCoolingPrivacy::EvaporativeCoolingPrivacy(Dataset* trainset, 
        Dataset* holdoset, Dataset* testset, Plink* plinkPtr) {
  // static class variables
  Q_EPS = 0.005;
  MAX_ITERATIONS = 100;
  // pointer to a PLINK environment
  PP = plinkPtr;
  // TODO: sim parameters; do i need these?
  numInstances = trainset->NumInstances();
  probBiological = 0.1;
  // algorithm parameters
  train = trainset;
  holdout = holdoset;
  test = testset;
  curVarNames = trainset->GetVariableNames();
  curVarMap = trainset->MaskGetAttributeMask(NUMERIC_TYPE);
  numAttributes = trainset->NumVariables();
  minRemainAttributes = 2;
  startTemp = par::ecStartTemp;
  finalTemp = par::ecFinalTemp;
  kConstant = 1;
  // Trang: larger tau takes longer to get to Tmin (finalTemp)
  // > tau <- d^2/2 # larger tau takes longer to get to Tmin
  tau = (numAttributes * numAttributes) / 2.0; 
  numToRemovePerIteration = 1;
  threshold = 1 / sqrt(numInstances);
  tolerance = 0.2 / sqrt(numInstances);
  iteration = 0;
  PrintState();
}

EvaporativeCoolingPrivacy::~EvaporativeCoolingPrivacy() {
}

bool EvaporativeCoolingPrivacy::ComputeScores() {
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() START\n");
  PP->printLOG(Timestamp() + "Running temperature schedule until final temperature\n");

  string iterationOutputFile = par::output_file_name + ".privacy.iterations.tab";
  PP->printLOG(Timestamp() + "Writing iteration error results to [ " + 
    iterationOutputFile +  " ]\n");
 	ofstream iterationOutputStream(iterationOutputFile);
	if(!iterationOutputStream.is_open()) {
		error("Could not open iteration output file [ " + iterationOutputFile + " ]\n");
	}
  iterationOutputStream << "Iteration\tTemperature\tKeep\tRemove\tTrain\tHoldout\tTest" << endl;

  PP->printLOG(Timestamp() + "Entering EVAPORATIVE COOLING cooling schedule\n");
  attributeProbabilty.clear();
  currentTemp = startTemp;
  PP->printLOG(Timestamp() + "Starting temperature: " + dbl2str(startTemp) + "\n");
  iteration = 0;
  while((currentTemp > finalTemp) && (train->NumVariables() >= minRemainAttributes)) {
    ++iteration;
    if(iteration > MAX_ITERATIONS) {
      error("EvaporativeCoolingPrivacy::ComputeScores() Maximum ticks reached " 
              + int2str(iteration) + "\n");
    }
    PP->printLOG(Timestamp() + "--------------------------------------------\n");
    PP->printLOG(Timestamp() + "iteration = " + int2str(iteration) + "\n");
    PrintState();
    curVarNames = train->GetVariableNames();
    curVarMap = train->MaskGetAttributeMask(NUMERIC_TYPE);
    // recompute p_t, p_h and delta_t
    PP->printLOG(Timestamp() + "Recomputing p_t, p_h and delta_t\n");
    ComputeImportance();
    //ComputeInverseImportance();
    ComputeDeltaQ();
    ComputeAttributeProbabilities();
    GenerateRandomUniformProbabilities();
    // EvaporateWorstAttribute();
    EvaporateWorstAttributes();
    PP->printLOG(Timestamp() + "Keeping : " + int2str(keepAttrs.size()) + "\n");
    PP->printLOG(Timestamp() + "Removing: " + int2str(removeAttrs.size()) + "\n");
    PP->printLOG(Timestamp() + "Dataset : " + int2str(train->NumVariables()) + "\n");
    ComputeBestAttributesErrors();
    iterationOutputStream 
            << iteration << "\t"
            << currentTemp << "\t"
            << keepAttrs.size() << "\t"
            << removeAttrs.size() << "\t"
            << trainError << "\t" 
            << holdError << "\t" 
            << testError << endl;
    // update current temperature T_t
    UpdateTemperature();
    PP->printLOG(Timestamp() + "New temperature: " + dbl2str(currentTemp) + "\n");
  } // end temperature schedule
  PP->printLOG(Timestamp() + "Exiting EVAPORATIVE COOLING cooling schedule after " + int2str(iteration) + " iterations\n");
  PP->printLOG(Timestamp() + "Final temperature:              " + dbl2str(currentTemp) + "\n");
  PP->printLOG(Timestamp() + "Final kept attributes set size: " + int2str(keepAttrs.size()) + "\n");
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() END\n");

  iterationOutputStream.close();
  
  return true;
}

void EvaporativeCoolingPrivacy::PrintState() {
  cout << "CURRENT STATE" << endl;
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
  cout << "num attributes:     " << numAttributes << endl;
  cout << "num instances:      " << numInstances << endl;
  cout << "k:                  " << kConstant << endl;
  cout << "last train error:   " << trainError << endl;
  cout << "last holdout error: " << holdError << endl;
  cout << "last test error:    " << testError << endl;
}

bool EvaporativeCoolingPrivacy::ComputeImportance() {
  ReliefF* relief = new ReliefF(train, PP, NUMERIC_ONLY_ANALYSIS);
  relief->ComputeAttributeScores();
  trainImportance = relief->GetScores();
  sort(trainImportance.begin(), trainImportance.end(), scoresSortAscByName);
  delete relief;

  ReliefF* rf = new ReliefF(holdout, PP, NUMERIC_ONLY_ANALYSIS);
  holdoutImportance = rf->ComputeScores();
  sort(holdoutImportance.begin(), holdoutImportance.end(), scoresSortAscByName);
  delete rf;
           
  return true;
}

bool EvaporativeCoolingPrivacy::ComputeDeltaQ() {
  PP->printLOG(Timestamp() + "Find delta Q = maximum absolute difference in training and " 
                "holdout inverse importance scores\n");
  AttributeScoresCIt trainIt = trainImportance.begin();
  AttributeScoresCIt holdoIt = holdoutImportance.begin();
//  AttributeScoresCIt trainIt = trainInvImportance.begin();
//  AttributeScoresCIt holdoIt = holdoutInvImportance.begin();
  vector<double> diffScores;
//  for(; trainIt != trainInvImportance.end(); ++trainIt, ++holdoIt) {
  for(; trainIt != trainImportance.end(); ++trainIt, ++holdoIt) {
    // FIXME! use key of iterator to lookup in other map
    double absDiff = fabs((*trainIt).first - (*holdoIt).first);
    diffScores.push_back(absDiff);
  }
  deltaQ = *max_element(diffScores.begin(), diffScores.end());
  
  return true;
}

bool EvaporativeCoolingPrivacy::ComputeAttributeProbabilities() {
    PP->printLOG(Timestamp() + "Computing probabilities of attributes P(a)\n");
    uint n = curVarNames.size();
    attributeProbabilty.clear();
    attributeProbabilty.resize(n);
    summedProbabilities.clear();
    summedProbabilities.resize(n);
    //    AttributeScoresCIt trainIt = trainInvImportance.begin();
    //    |   PAs <- exp(-q1.scores/(2*delta.q*k*myT))
    //    >   sumPAs <- sum(PAs)
    //    >   scaled.PAs <- PAs/sum(PAs)
    //    >   cum.scaled.PAs <- cumsum(scaled.PAs)
    double denom = 2.0 * deltaQ * kConstant * currentTemp;
    AttributeScoresCIt trainIt = trainImportance.begin();
    for(uint i=0; trainIt != trainImportance.end(); ++trainIt, ++i) {
      double q1Score = (*trainIt).first;
      double PA = -q1Score / denom;
      attributeProbabilty[i] = -PA / denom;
      summedProbabilities[i] += attributeProbabilty[i];
    }
    scaledSummedProbabilities.clear();
    scaledSummedProbabilities.resize(n);
    cummulativeProbabilities.clear();
    cummulativeProbabilities.resize(n);
    for(uint i=0; trainIt != trainImportance.end(); ++trainIt, ++i) {
      scaledSummedProbabilities[i] = attributeProbabilty[i] / summedProbabilities[i];
      cummulativeProbabilities[i] += scaledSummedProbabilities[i];
    }
    return true;
}

bool EvaporativeCoolingPrivacy::GenerateRandomUniformProbabilities() {
  //    >   prob.rands <- runif(1, min = 0, max = 1)
  PP->printLOG(Timestamp() + "Generating uniform probabilities in (0, 1)\n");
  uniform_real_distribution<double> runif(0, 1);
  // generate random uniform probabilities
  randUniformProbs.clear();
  for(uint prIdx=0; prIdx < trainInvImportance.size(); ++prIdx) {
    randUniformProbs.push_back(runif(engine));
  }
  return true;
}

bool EvaporativeCoolingPrivacy::EvaporateWorstAttributes() {
  removeAttrs.clear();
  keepAttrs.clear();
  PP->printLOG(Timestamp() + "Evaporating the worst attributes\n");
  uint removedSoFar = 0;
  for(uint pIdx=0; pIdx < randUniformProbs.size(); ++pIdx) {
    string thisVar = curVarNames[pIdx];
    uint thisVarIdx = curVarMap[thisVar];
    if(par::verbose) {
      PP->printLOG(Timestamp() + thisVar + " (" + int2str(thisVarIdx) + ")\n");
      PP->printLOG(Timestamp() + dbl2str(randUniformProbs[pIdx]) + 
        " < "  + dbl2str(cummulativeProbabilities[pIdx]) + "?\n");
    }
    if((removedSoFar < numToRemovePerIteration) && 
            randUniformProbs[pIdx] < cummulativeProbabilities[pIdx]) {
      if(par::verbose) PP->printLOG(Timestamp() + " => remove\n");
      train->MaskRemoveVariable(thisVar);
      holdout->MaskRemoveVariable(thisVar);
      test->MaskRemoveVariable(thisVar);
      removeAttrs.push_back(thisVar);
      ++removedSoFar;
    } else {
      if(par::verbose) PP->printLOG(Timestamp() + " => keep\n");
      keepAttrs.push_back(thisVar);
    }
  }
  return true;
}

bool EvaporativeCoolingPrivacy::EvaporateWorstAttribute() {
  //    >   num.remv <- 1 # only remove 1 attribute  
  //    >   remv.atts <- kept.atts[prob.rands < cum.scaled.PAs][1] 
  removeAttrs.clear();
  keepAttrs.clear();
  sort(trainImportance.begin(), trainImportance.end(), scoresSortAsc);
  pair<double, string> worstAttr = trainImportance[0];
  string thisVar = worstAttr.second;
  double thisVarValue = worstAttr.first;
  sort(trainImportance.begin(), trainImportance.end(), scoresSortAscByName);
  PP->printLOG(Timestamp() + "Evaporating the single worst attribute [" + thisVar + "]\n");
  if(par::verbose) PP->printLOG(Timestamp() + thisVar + " (" + int2str(thisVarValue) + ")\n");
  removeAttrs.push_back(thisVar);
  train->MaskRemoveVariable(thisVar);
  holdout->MaskRemoveVariable(thisVar);
  test->MaskRemoveVariable(thisVar);
  keepAttrs = train->GetVariableNames();

  return true;
}

bool EvaporativeCoolingPrivacy::ComputeBestAttributesErrors() {
  PP->printLOG(Timestamp() + "Compute errors for the best attributes\n");
  // get training, holdout and prediction errors with best attributes (keepAttrs)
  trainError = ClassifyAttributeSet(keepAttrs, TRAIN);
  trainErrors.push_back(trainError);
  holdError = ClassifyAttributeSet(keepAttrs, HOLDOUT);
  holdoutErrors.push_back(holdError);
  
//  if(fabs(trainError - holdError) < threshold + rnorm(1, 0, tolerance)) {
//    holdError <- trainError
//  } else {
//    holdError <- holdError + rnorm(1, 0, tolerance)
//  }
  
  testError = ClassifyAttributeSet(keepAttrs, TEST);
  testErrors.push_back(testError);
  PP->printLOG(Timestamp() + "* train:   " + dbl2str(trainError) + "\n" +
               Timestamp() + "* holdout: " + dbl2str(holdError) + "\n" +
               Timestamp() + "* test:    " + dbl2str(testError) + "\n");
}

double EvaporativeCoolingPrivacy::ClassifyAttributeSet(vector<string> attrs, 
        DATASET_TYPE dataType) {

  double retError = 1.0;
  RandomForest* randomForest = 0;
  switch(dataType) {
    case TRAIN:
      PP->printLOG(Timestamp() + "Classify best attributes for TRAINING data\n");
      randomForest = new RandomForest(train, par::trainFile, attrs, false);
      randomForest->ComputeScores();
      retError = randomForest->GetClassificationError();
      randomForest->SaveForest();
      break;
    case HOLDOUT:
      PP->printLOG(Timestamp() + "Classify best attributes for HOLDOUT data\n");
      randomForest = new RandomForest(holdout, par::holdoutFile, attrs, false);
      randomForest->ComputeScores();
      retError = randomForest->GetClassificationError();
      break;
    case TEST:
      PP->printLOG(Timestamp() + "Predicting best attributes for TESTING data\n");
      randomForest = new RandomForest(test, par::testFile, attrs, true);
      randomForest->Predict();
      retError = randomForest->GetClassificationError();
      break;
    default:
      error("EvaporativeCoolingPrivacy::ClassifyAttributeSet Dataset type no recognized");
  }
  delete randomForest;
  
  return retError;
}

bool EvaporativeCoolingPrivacy::UpdateTemperature() {
  PP->printLOG(Timestamp() + "Updating temperature\n");
  // myT <- mean(q1.scores[, 1] ^ 2) / (k * delta.q)
//  double T_t_sum = 0;
//  AttributeScoresCIt trainIt = trainInvImportance.begin();
//  for(uint idx=0; trainIt != trainInvImportance.end(); ++trainIt, ++idx) {
//    double PA = (*trainIt).first;
//    T_t_sum += ((PA * PA) / ( kConstant * deltaQ));
//  }
//  currentTemp = (T_t_sum / trainInvImportance.size());
  // update 11/1/16 from Trang:
  // myT <- myT*exp(-i/tau)
  currentTemp = currentTemp * exp(-iteration / tau);
  
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
