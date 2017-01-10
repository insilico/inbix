#ifndef EVAPORATIVECOOLINGPRIVACY_H
#define EVAPORATIVECOOLINGPRIVACY_H

/* 
 * File:   EvaporativeCoolingPrivacy.h
 * Author: bwhite
 *
 * Created on October 18, 2016, 10:57 PM
 * THis file is a part of inbix which in turn is based on PLINK, so
 * certain assumptions may apply.
 */

#include <string>
#include <vector>
#include <armadillo>
#include <random>

#include "plink.h"

#include "Dataset.h"

enum DATASET_TYPE {
  TRAIN, 
  HOLDOUT, 
  TEST
};

typedef std::pair<std::vector<std::string>, std::vector<std::string>> ResultsLists;

class EvaporativeCoolingPrivacy {
public:
  EvaporativeCoolingPrivacy(Dataset* trainset, Dataset* holdoset, 
                            Dataset* testset, Plink* plinkPtr);
  virtual ~EvaporativeCoolingPrivacy();
  bool ComputeScores();
  void PrintState();
  ResultsLists GetKeptRemoved();
  bool WriteBestAttributes(std::string fileSuffix);
  std::pair<uint, double> CheckDetectedAttributes();
private:
  bool ComputeImportance();
  bool ComputeAttributeProbabilities();
  bool GenerateRandomUniformProbabilities();
  uint EvaporateWorstAttributes(uint numToRemove);
  bool EvaporateWorstAttribute();
  double ClassifyAttributeSet(std::vector<std::string> attrs, DATASET_TYPE);
  bool ComputeBestAttributesErrors();
  bool UpdateTemperature();
  bool ComputeInverseImportance();
  uint CurrentNumberCorrect(std::vector<std::string> testSet);
  bool RemoveImportanceScore(std::string varToRemove);
  
  // Mersenne twister random number engine - based Mersenne prime 2^19937 − 1
  std::mt19937_64 engine;  
  // PLINK environment Plink object
  Plink* PP;
  // CONSTANTS
  double Q_EPS;
  uint MAX_ITERATIONS;
  
  // algorithm
  uint minRemainAttributes;
  uint updateInterval;
  uint iteration;
  uint update;
  
  // classification data sets
  Dataset* train;
  Dataset* holdout;
  Dataset* test;
  uint numInstances;
  uint numSignalsInData;
  vector<string> signalNames;
  
  // track original and current set of variables and their indices
  // into the original data set
  std::vector<std::string> origVarNames;
  std::map<std::string, unsigned int> origVarMap;
  std::vector<std::string> curVarNames;
  std::map<std::string, unsigned int> curVarMap;
  
  // importance/quality scores
	AttributeScores trainImportance;       // q_t
	AttributeScores holdoutImportance;     // q_h
	AttributeScores trainInvImportance;    // p_t
	AttributeScores holdoutInvImportance;  // p_h

  // computing importance
  std::map<std::string, double> diffImportance;
  std::vector<double> diffScores;
  std::vector<double> diff;
  double deltaQ;
  double threshold;
  double tolerance;
    
  // temperature schedule
  double startTemp;
  double currentTemp;
  double finalTemp;
  double tau;
  
  // evaporation phase variables
  uint numToRemovePerIteration;
  std::vector<double> attributeProbabilty;
  double summedProbabilities;
  std::vector<double> scaledProbabilities;
  std::vector<double> cummulativeProbabilities;
  std::vector<double> randUniformProbs;
  double randUniformValue;
  vector<string> possiblyRemove;
  std::vector<std::string> removeAttrs;
  std::vector<std::string> keepAttrs;

  // classification/prediction errors
  double randomForestPredictError;
  double trainError;
  double holdError;
  double testError;
  std::vector<double> trainErrors;
  std::vector<double> holdoutErrors;
  std::vector<double> testErrors;
};

#endif /* EVAPORATIVECOOLINGPRIVACY_H */
