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
        Dataset* holdoset, Dataset* testset, double t0, double tf, 
        Plink* plinkPtr) {
  // static class variables
  EvaporativeCoolingPrivacy::Q_EPS = 0.005;
  EvaporativeCoolingPrivacy::MAX_TICKS = 100;
  // pointer to a PLINK environment
  PP = plinkPtr;
  // algorithm parameters
  train = trainset;
  holdout = holdoset;
  test = testset;
  S_t = trainset->GetNumericMatrix();
  S_h = holdoset->GetNumericMatrix();
  S_A = testset->GetVariableNames();
  d = trainset->NumVariables();
  m = d;
  T_0 = t0;
  T_f = tf;
}

EvaporativeCoolingPrivacy::~EvaporativeCoolingPrivacy() {
}

pair<vector<double>, vector<double> > 
EvaporativeCoolingPrivacy::ComputeImportance() {
  ReliefF* relief = new ReliefF(train, PP, NUMERIC_ONLY_ANALYSIS);
  relief->ComputeAttributeScores();
  trainImportance = relief->GetScores();
  sort(trainImportance.begin(), trainImportance.end(), scoresSortAscByName);
  vector<double> train_relief;
  train_relief.reserve(trainImportance.size());
  for_each(trainImportance.begin(), trainImportance.end(),
           [&train_relief](const map<double, string>::value_type& p) 
           { train_relief.push_back(p.first); });
  double min_train = *(min_element(train_relief.begin(), train_relief.end()));
  vector<double> returnTrain;
  for(uint idx=0; idx < train_relief.size(); ++idx) {
    double x = train_relief[idx];
    returnTrain.push_back(1 / (x - min_train + Q_EPS));
  }

  ReliefF* rf = new ReliefF(holdout, PP, NUMERIC_ONLY_ANALYSIS);
  holdoutImportance = rf->ComputeScores();
  sort(holdoutImportance.begin(), holdoutImportance.end(), scoresSortAscByName);
  vector<double> holdo_relief;
  holdo_relief.reserve(holdoutImportance.size());
  for_each(holdoutImportance.begin(), holdoutImportance.end(),
           [&holdo_relief](const map<double, string>::value_type& p) 
           { holdo_relief.push_back(p.first); });
  double min_holdo = *min_element(holdo_relief.begin(), holdo_relief.end());
  vector<double> returnHoldo;
  returnHoldo.reserve(holdo_relief.size());
  for(uint idx=0; idx < holdo_relief.size(); ++idx) {
    double x = holdo_relief[idx];
    returnHoldo.push_back(1 / (x - min_holdo + Q_EPS));
  }
  
  delete relief;
  delete rf;
           
  return make_pair(returnTrain, returnHoldo);
}

bool EvaporativeCoolingPrivacy::ComputeImportanceInternal() {
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

bool EvaporativeCoolingPrivacy::ComputeScores() {
  PP->printLOG("EvaporativeCoolingPrivacy::ComputeScores() START\n");
  uniform_real_distribution<double> runif(0, 1);
  uint n = 100;
  double pb = 0.1;
  m = d;
  vector<string> A = S_A;
  PP->printLOG("Computing importance scores for training and holdout sets\n");
  ComputeImportanceInternal();
  // absolute value of the difference
  PP->printLOG("Find delta Q = maximum absolute difference in training and " 
                "holdout importance scores\n");
  AttributeScoresCIt trainIt = trainImportance.begin();
  AttributeScoresCIt holdoIt = holdoutImportance.begin();
  vector<double> diff_scores(trainImportance.size());
  uint idx = 0;
  for(idx=0; trainIt != trainImportance.end(); ++trainIt, ++holdoIt, ++idx) {
    // FIXME! use key of iterator to lookup in other map
    double abs_diff = abs((*trainIt).first - (*holdoIt).first);
    diff_scores[idx] = abs_diff;
  }
  double delta_q = *max_element(diff_scores.begin(), diff_scores.end());

  PP->printLOG("Running temperature schedule until final temperature\n");
  vector<string> removeAttrs;
  vector<string> keepAttrs;
  uint k = 1;
  double T_t = T_0;
  PP->printLOG("Starting temperature: " + dbl2str(T_0) + "\n");
  uint t = 1;
  vector<double> PAs;
  while((T_t > T_f) && train->NumAttributes()) {
    PP->printLOG("t = " + int2str(t) + "\n");
    // probabilities of attributes P(a)
    PP->printLOG("Computing probabilities of attributes P(a)\n");
    idx = 0;
    PAs.clear();
    PAs.resize((trainImportance.size()));
    AttributeScoresCIt trainIt = trainImportance.begin();
    for(idx=0; trainIt != trainImportance.end(); ++trainIt, ++idx) {
      double PA = (*trainIt).first;
      PAs[idx] = exp(-(PA * PA) / (2 * delta_q * k * T_t));
    }
    // generate random uniform probabilities
    vector<double> prob_rands(trainImportance.size());
    for(uint prIdx=0; prIdx <- trainImportance.size(); ++ prIdx) {
      prob_rands.push_back(runif(engine));
    }

    // evaporate the worst
    PP->printLOG("Evaporating the worst attributes\n");
    for(uint pIdx=0; pIdx < PAs.size(); ++pIdx) {
      string thisVar = S_A[pIdx];
      if(PAs[pIdx] > prob_rands[pIdx]) {
        train->MaskRemoveVariable(thisVar);
        test->MaskRemoveVariable(thisVar);
        removeAttrs.push_back(thisVar);
      } else {
        keepAttrs.push_back(thisVar);
      }
    }
    PP->printLOG("Keeping : " + int2str(keepAttrs.size()) + "\n");
    PP->printLOG("Removing: " + int2str(removeAttrs.size()) + "\n");
    
    // classify independent test data set with random forest
    PP->printLOG("Classifying independent test data set with random forest\n");
    RandomForest* randomForest = new RandomForest(train, PP);
    randomForestPredictError = randomForest->Predict(test);
    PP->printLOG("New prediction error: " + dbl2str(randomForestPredictError) + "\n");
  
    // recompute p_t, p_h and delta_t
    PP->printLOG("Recomputing p_t, p_h and delta_t\n");
    ComputeImportance();
    vector<double> diff_scores;
    trainIt = trainImportance.begin();
    holdoIt = holdoutImportance.begin();
    for(idx=0; trainIt != trainImportance.end(); ++trainIt, ++holdoIt, ++idx) {
      diff_scores.push_back((*trainIt).first - (*holdoIt).first);
    }
    double delta_q = *max_element(diff_scores.begin(), diff_scores.end());

    // update current temperature T_t
    PP->printLOG("Updating temperature\n");
    double T_t_sum = 0;
    trainIt = trainImportance.begin();
    for(idx=0; trainIt != trainImportance.end(); ++trainIt, ++idx) {
      double PA = (*trainIt).first;
      T_t_sum += ((PA * PA) / ( k * delta_q));
    }
    double T_t_mean = T_t_sum / trainImportance.size();
    T_t = T_t_mean;
    PP->printLOG("New temperature: " + dbl2str(T_t) + "\n");

    // tick
    ++t;
    if(t > MAX_TICKS) {
      error("EvaporativeCoolingPrivacy::ComputeScores() Maximum ticks reached " 
              + int2str(t) + "\n");
    }
  }

  PP->printLOG("Final temperature: " + dbl2str(T_t) + "\n");
  
  PP->printLOG("EvaporativeCoolingPrivacy::ComputeScores() END\n");

  return true;
}

bool EvaporativeCoolingPrivacy::ComputeScoresTrang() {
//  uniform_real_distribution<double> runif(0, 1);
//  // S_A = set of all predictors
//  uint n = 100;
//  uint d = 5000;
//  double pb = 0.1;
//
//  pair<vector<double>, vector<double> > 
//    important_scores = ComputeImportance(X_train, X_holdo);
//  vec q1_scores = important_scores.first;
//  vec q2_scores = important_scores.second;
//  
//  
//  double fholds = 0.5;
//  double ftests = 0.5;
//  double ftrains = 0.5;
//  vector<double> correct_detect_ec;
//  //double oldAccuracy = 0.5;
//  uint num_att = att_names.size(); // number of attributes
//  uint n = att_names.size();
//  double threshold = 1.0 / sqrt(n);
//  double tolerance = 0.2 / sqrt(n);
//  vector<uint> num_atts = {99999999, num_att - 1};
//  vector<string> kept_atts = att_names;
//  //var_names = list();
//  uint m = d; // max number of removed variables after each iteration
//
//  uint T0 = 10;
//  uint Tmin = 1;
//  uint tau = 100; // or 100, larger tau takes longer to get to Tmin
//  uint i = 1;
//  uint k = 1;
//  uint myT = T0;
//  vec q1_scores_plot = q1_scores;
//  while (myT > Tmin) {
//    uint lastIdx = num_atts.size() - 1;
//    uint nextLastIdx = lastIdx - 1;
//    if((num_atts[lastIdx] < 1) && (num_atts[lastIdx] == num_atts[nextLastIdx])) {
//      continue;
//    }
//    
//    //   total_probs = sum(exp(diff_scores*q1_scores/(k*myT)))
//    //   PAs = exp(diff_scores*q1_scores/(k*myT))*num_att/total_probs
//    //   print(max(diff_scores))
//    mat PAs = exp(-q1_scores ^ 2 / (2 * delta_q * k * myT));
//    //   print(head(PAs))
//    vec prob_rands;
//    prob_rands.imbue( [&]() { return runif(engine); } );
//    mat goodness_of_atts = PAs - prob_rands;
//    goodness_of_atts = goodness_of_atts[order(goodness_of_atts[,1]), ];
//    uint num_remv = min(sum(goodness_of_atts < 0), m);    // remove m worst attributes or less
//    vector<string> remv_atts = rownames(goodness_of_atts)[1:num_remv];
//    if (num_remv >= goodness_of_atts.n_rows){
//      kept_atts = NULL;
//      break;
//    }
//    else {
//      // TODO: unpack this
//      kept_atts = rownames(goodness_of_atts)[(num_remv+1):num_att];
//    }
//
//    vector<string> keep_atts = {kept_atts, "pheno"};
//    mat new_X_train = X_train(,kept_atts);
//    mat new_X_holdo = X_holdo(,kept_atts);
//    mat new_X_test = X_test(,kept_atts);
//    pair<vector<double>, vector<double> > new_scores = 
//      ComputeImportance(new_X_train, new_X_holdo);
//    q1_scores = new_scores.first;
//    q2_scores = new_scores.second;
//    q1_scores_plot(remv_atts, ) = 1 / Q_EPS;
//    diff_scores = abs(q1_scores - q2_scores);
//    delta_q = max(diff_scores);
//    att_names = kept_atts;
//    //var_names[[i]] = rownames(q1_scores_plot)[q1_scores_plot<1/Q_EPS];
//    // ? include maximization of holdo accuracy  
//
//    // Compute train and holdo accuracy for new S_A attributes: 
//    //   original privacy classifier:
//    //   weights = c(sign(q1_scores[,1]),0)
//    //   ftrain = sum(sign(data_matrix(new_X_train) %*% weights) == new_X_train$pheno)/n
//    //   fholdo = sum(sign(data_matrix(new_X_holdo) %*% weights) == new_X_holdo$pheno)/n
//    //   ftest = sum(sign(data_matrix(new_X_test) %*% weights) == new_X_test$pheno)/n
//
//    // Now, randomForest:
//    RandomForest* randomForest = new RandomForest(train, PP);
//    result_rf = randomForest(pheno ~ _, data = new_X_train, importance = T, ntree = 1000);
//  //   conf_mat = result_rf$confusion[,1:2]
//  //   ftrain = (conf_mat[1,1] + conf_mat[2,2])/sum(conf_mat)
//    ftrain = 1 - mean(result_rf$confusion[,"class_error"]);
//    holdo_pred = predict(result_rf, newdata = new_X_holdo);
//    fholdo = mean(holdo_pred == X_holdo$pheno);
//  //   print(fholdo)
//    test_pred = predict(result_rf, newdata = new_X_test);
//    ftest = mean((test_pred == X_test$pheno));
//  //   print(ftest)
//    if (abs(ftrain - fholdo) < threshold + rnorm(1, 0, tolerance)) {
//      fholdo = ftrain;
//    } else {
//      fholdo = fholdo + rnorm(1, 0, tolerance);
//    }
//
//  //   tempAccuracy = fholdo
//  //   deltaAccuracy = tempAccuracy - oldAccuracy
//  //   temp_rand = runif(1, min = 0, max = 1)
//  //   if (temp_rand < exp(deltaAccuracy/myT)){
//  //     oldAccuracy = tempAccuracy // accept S_A, else next while
//  //   }   
//    fholds = c(fholds, fholdo);
//    ftests = c(ftests, ftest);
//    ftrains = c(ftrains, ftrain);
//  //   myoldT = myT*exp(-i/tau) // cooling schedule, lower T
//    myT = mean(q1_scores[,1]^2)/(k*delta_q);
//  //   print(myT)
//  //   cat(myoldT, myT)
//  //   myT = myT*exp(-i/tau) // cooling schedule, lower T
//    num_att = att_names.size();
//    num_atts = c(num_atts, num_att);
//    correct_detect_ec = c(correct_detect_ec, sum(var_names[[i]] %in% signal_names));
//
//    i = i + 1;
//  //   print(myT)
//  }
//
//  num_atts = num_atts[-1];
//  fplots = data_frame(num_atts, ftrain = ftrains, fholdo = fholds, ftest = ftests, alg = 1);
//
//  save(fplots, melted_fs, correct_detect_ec, n, d, pb, X_train, X_holdo, type, types,
//       X_test, signal_names, myfilename, bias, shortname, threshold, tolerance, bias, biases,
//       file = paste("results/privateECresult", shortname, type, "_Rdata", sep = ""));
  
  return true;
}
