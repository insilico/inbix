#ifndef EVAPORATIVECOOLINGPRIVACY_H
#define EVAPORATIVECOOLINGPRIVACY_H

/* 
 * File:   EvaporativeCoolingPrivacy.h
 * Author: bwhite
 *
 * Created on October 18, 2016, 10:57 PM
 */

#include <string>
#include <vector>
#include <armadillo>
#include <random>

#include "plink.h"

#include "Dataset.h"

class EvaporativeCoolingPrivacy {
public:
  EvaporativeCoolingPrivacy(Dataset* trainset, Dataset* holdoset, 
                            Dataset* testset, double t0, double tf,
                            Plink* plinkPtr);
  virtual ~EvaporativeCoolingPrivacy();
  bool ComputeScores();
  bool ComputeScoresTrang();
private:
  Plink* PP;
  Dataset* train;
  Dataset* holdout;
  Dataset* test;
  std::vector<std::string> S_A;
  arma::mat S_t;
  arma::mat S_h;
  arma::mat X_test;
  std::pair<std::vector<double>, std::vector<double> > ComputeImportance();
  bool ComputeImportanceInternal();
  double T_0;
  double T_f;
  // max number of removed variables after each iteration
  uint m;
  uint d;
  
	AttributeScores trainImportance;    // p_t
	AttributeScores holdoutImportance;  // p_h
  double randomForestPredictError;

  double Q_EPS;
  uint MAX_TICKS;
  
  // random distributions
  std::mt19937_64 engine;  // Mersenne twister random number engine
};

#endif /* EVAPORATIVECOOLINGPRIVACY_H */
