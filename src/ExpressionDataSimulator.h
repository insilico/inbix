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

class ExpressionDataSimulator {
public:
  ExpressionDataSimulator(Plink* plinkPtr);
  ExpressionDataSimulator(const ExpressionDataSimulator& orig);
  virtual ~ExpressionDataSimulator();
  bool Simulate(uint n, uint d, double pb, double bias, std::string type);
  const arma::mat& GetData() { return simulatedData; }
private:
  bool CreateDiffCoexpMatrixNoME(uint M, uint N, double meanExpression,
    arma::mat A, double randSdNoise, double sdNoise, 
    std::vector<uint> sampleIndicesInteraction);
  bool SimulateData(uint n_e, uint n_db, uint n_ns, 
    std::vector<std::string> sv_db, std::vector<std::string> sv_ns,
    double sd_b, double sd_gam, double sd_u, bool conf, std::string distr_db,
    double p_b, double p_gam, double p_ov);
  Plink* PP;
  arma::mat simulatedData;
};

#endif /* EXPRESSIONDATASIMULATOR_H */
