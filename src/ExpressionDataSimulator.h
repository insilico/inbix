/* 
 * File:   ExpressionDataSimulator.h
 * Author: bwhite
 *
 * Created on October 18, 2016, 8:02 PM
 */

#include <string>
#include <vector>
#include <armadillo>

#include "plink.h"

#ifndef EXPRESSIONDATASIMULATOR_H
#define EXPRESSIONDATASIMULATOR_H

/**
 * \enum SplitType.
 * Train, holdout and test splits.
 */
enum SplitType{
  TRAIN_SPLIT, /**< training data split type */
  HOLDOUT_SPLIT, /**< holdout data split type */
  TEST_SPLIT, /**< test data split type */
  NO_SPLIT /**< default no type */
};

class ExpressionDataSimulator {
public:
  ExpressionDataSimulator(Plink* plinkPtr);
  virtual ~ExpressionDataSimulator();
  bool CreateSimulation(uint n, uint d, double pb, double bias, std::string type);
  const arma::mat& GetData(SplitType splitType);
private:
  bool CreateDiffCoexpMatrixNoME(uint M, uint N, double meanExpression,
    arma::mat A, double randSdNoise, double sdNoise, 
    std::vector<uint> sampleIndicesInteraction);
  bool SimulateData(uint n_e, uint n_db, uint n_ns, 
    std::vector<std::string> sv_db, std::vector<std::string> sv_ns,
    double sd_b, double sd_gam, double sd_u, bool conf, std::string distr_db,
    double p_b, double p_gam, double p_ov);
  Plink* PP;
  utin dimM;
  uint dimN;
  uint n1;
  uint n2;
  arma::mat simulatedDataBase;
  std::vector<std::string> subIds;
  std::vector<std::string> colNames;
  std::vector<uint> phenosBase;
  arma::mat simulatedDataTrain;
  arma::mat simulatedDataHoldout;
  arma::mat simulatedDataTest;
};

#endif /* EXPRESSIONDATASIMULATOR_H */
