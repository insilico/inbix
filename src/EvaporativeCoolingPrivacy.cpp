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

#include "plink.h"
#include "options.h"
#include "helper.h"

#include "Insilico.h"
#include "Dataset.h"

#include "EvaporativeCoolingPrivacy.h"
#include "ReliefF.h"
#include "RandomForest.h"
#include "ArmadilloFuncs.h"

EvaporativeCoolingPrivacy::EvaporativeCoolingPrivacy(Dataset* trainset, 
        Dataset* holdoset, Dataset* testset, Plink* plinkPtr) {
  // static class variables
  Q_EPS = 0.005;
  MAX_TICKS = 100;
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
  // TODO: what is this? 'm' in Trang's code
  minRemainAttributes = 2;
  startTemp = par::ecStartTemp;
  finalTemp = par::ecFinalTemp;
  kConstant = 1;
}

EvaporativeCoolingPrivacy::~EvaporativeCoolingPrivacy() {
}

bool EvaporativeCoolingPrivacy::ComputeScores() {
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() START\n");
  PP->printLOG(Timestamp() + "Running temperature schedule until final temperature\n");
  ComputeInverseImportance();
  deltaQ = DeltaQ();
  uint tick = 1;
  selectProbabilty.clear();
  currentTemp = startTemp;
  PP->printLOG(Timestamp() + "Starting temperature: " + dbl2str(startTemp) + "\n");
  while((currentTemp > finalTemp) && (train->NumVariables() >= minRemainAttributes)) {
    PP->printLOG(Timestamp() + "--------------------------------------------\n");
    PP->printLOG(Timestamp() + "time tick = " + int2str(tick) + "\n");
    curVarNames = train->GetVariableNames();
    curVarMap = train->MaskGetAttributeMask(NUMERIC_TYPE);
    ComputeProbabilities();
    GenerateUniformRands();
    EvaporateWorst();
    PP->printLOG(Timestamp() + "Keeping : " + int2str(keepAttrs.size()) + "\n");
    PP->printLOG(Timestamp() + "Removing: " + int2str(removeAttrs.size()) + "\n");
    PP->printLOG(Timestamp() + "Dataset : " + int2str(train->NumVariables()) + "\n");
    ComputeBestAttributesErrors();
    // recompute p_t, p_h and delta_t
    PP->printLOG(Timestamp() + "Recomputing p_t, p_h and delta_t\n");
    ComputeInverseImportance();
    deltaQ = DeltaQ();
    // update current temperature T_t
    UpdateTemperature();
    PP->printLOG(Timestamp() + "New temperature: " + dbl2str(currentTemp) + "\n");
    // tick
    ++tick;
    if(tick > MAX_TICKS) {
      error("EvaporativeCoolingPrivacy::ComputeScores() Maximum ticks reached " 
              + int2str(tick) + "\n");
    }
  } // end temperature schedule
  PP->printLOG(Timestamp() + "Final temperature:              " + dbl2str(currentTemp) + "\n");
  PP->printLOG(Timestamp() + "Final kept attributes set size: " + int2str(keepAttrs.size()) + "\n");
  PP->printLOG(Timestamp() + "EvaporativeCoolingPrivacy::ComputeScores() END\n");

  return true;
}

bool EvaporativeCoolingPrivacy::ComputeInverseImportance() {
  PP->printLOG(Timestamp() + "Computing inverse importance scores for training and holdout sets\n");

  ReliefF relief(train, PP, NUMERIC_ONLY_ANALYSIS);
  relief.ComputeAttributeScores();
  trainInvImportance = relief.GetScores();
  sort(trainInvImportance.begin(), trainInvImportance.end(), scoresSortAscByName);
  vector<double> train_relief;
  vector<string> train_vars;
  for_each(trainInvImportance.begin(), trainInvImportance.end(),
           [&train_relief, &train_vars](const map<double, string>::value_type& p) 
           { train_relief.push_back(p.first); train_vars.push_back(p.second); });
  double min_train = *(min_element(train_relief.begin(), train_relief.end()));
  trainInvImportance.clear();
  for(uint idx=0; idx < train_relief.size(); ++idx) {
    double x = train_relief[idx];
    trainInvImportance.push_back(make_pair(1 / (x - min_train + Q_EPS), 
            train_vars[idx]));
  }

  ReliefF rf(holdout, PP, NUMERIC_ONLY_ANALYSIS);
  rf.ComputeAttributeScores();
  holdoutInvImportance = rf.GetScores();
  sort(holdoutInvImportance.begin(), holdoutInvImportance.end(), scoresSortAscByName);
  vector<double> holdo_relief;
  vector<string> holdo_vars;
  for_each(holdoutInvImportance.begin(), holdoutInvImportance.end(),
           [&holdo_relief, &holdo_vars](const map<double, string>::value_type& p) 
           { holdo_relief.push_back(p.first); holdo_vars.push_back(p.second); });
  double min_holdo = *min_element(holdo_relief.begin(), holdo_relief.end());
  holdoutInvImportance.clear();
  for(uint idx=0; idx < holdo_relief.size(); ++idx) {
    double x = holdo_relief[idx];
    holdoutInvImportance.push_back(make_pair(1 / (x - min_holdo + Q_EPS), 
            holdo_vars[idx]));
  }
           
  return true;
}

double EvaporativeCoolingPrivacy::DeltaQ() {
  PP->printLOG(Timestamp() + "Find delta Q = maximum absolute difference in training and " 
                "holdout inverse importance scores\n");
  AttributeScoresCIt trainIt = trainInvImportance.begin();
  AttributeScoresCIt holdoIt = holdoutInvImportance.begin();
  vector<double> diffScores;
  for(; trainIt != trainInvImportance.end(); ++trainIt, ++holdoIt) {
    // FIXME! use key of iterator to lookup in other map
    double absDiff = fabs((*trainIt).first - (*holdoIt).first);
    diffScores.push_back(absDiff);
  }
  return *max_element(diffScores.begin(), diffScores.end());
}

bool EvaporativeCoolingPrivacy::ComputeProbabilities() {
    PP->printLOG(Timestamp() + "Computing probabilities of attributes P(a)\n");
    selectProbabilty.clear();
    AttributeScoresCIt trainIt = trainInvImportance.begin();
    double denom = 2.0 * deltaQ * kConstant * currentTemp;
    for(; trainIt != trainInvImportance.end(); ++trainIt) {
      double PA = (*trainIt).first;
      selectProbabilty.push_back(exp(-(PA * PA) / denom));
    }
    return true;
}

bool EvaporativeCoolingPrivacy::GenerateUniformRands() {
  PP->printLOG(Timestamp() + "Generating uniform probabilities in (0, 1)\n");
  uniform_real_distribution<double> runif(0, 1);
  // generate random uniform probabilities
  randUniformProbs.clear();
  for(uint prIdx=0; prIdx < trainInvImportance.size(); ++prIdx) {
    randUniformProbs.push_back(runif(engine));
  }
  return true;
}

bool EvaporativeCoolingPrivacy::EvaporateWorst() {
  removeAttrs.clear();
  keepAttrs.clear();
  PP->printLOG(Timestamp() + "Evaporating the worst attributes\n");
  for(uint pIdx=0; pIdx < selectProbabilty.size(); ++pIdx) {
    string thisVar = curVarNames[pIdx];
    uint thisVarIdx = curVarMap[thisVar];
    if(par::verbose) PP->printLOG(Timestamp() + thisVar + " (" + int2str(thisVarIdx) + ")\n");
    if(selectProbabilty[pIdx] > randUniformProbs[pIdx]) {
      if(par::verbose) PP->printLOG(Timestamp() + " => remove\n");
      train->MaskRemoveVariable(thisVar);
      holdout->MaskRemoveVariable(thisVar);
      test->MaskRemoveVariable(thisVar);
      removeAttrs.push_back(thisVar);
    } else {
      if(par::verbose) PP->printLOG(Timestamp() + " => keep\n");
      keepAttrs.push_back(thisVar);
    }
  }
  return true;
}

bool EvaporativeCoolingPrivacy::UpdateTemperature() {
  PP->printLOG(Timestamp() + "Updating temperature\n");
  double T_t_sum = 0;
  AttributeScoresCIt trainIt = trainInvImportance.begin();
  for(uint idx=0; trainIt != trainInvImportance.end(); ++trainIt, ++idx) {
    double PA = (*trainIt).first;
    T_t_sum += ((PA * PA) / ( kConstant * deltaQ));
  }
  currentTemp = (T_t_sum / trainInvImportance.size());

  return true;
}

bool EvaporativeCoolingPrivacy::ComputeBestAttributesErrors() {
  PP->printLOG(Timestamp() + "Compute errors for the best attributes\n");
  // get training, holdout and prediction errors with best attributes (keepAttrs)
  trainError = ClassifyAttributeSet(keepAttrs, TRAIN);
  trainErrors.push_back(trainError);
  holdError = ClassifyAttributeSet(keepAttrs, HOLDOUT);
  holdoutErrors.push_back(holdError);
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
