/*==============================================================================
 *
 * Filename:  RegainMinimalMinimal.h - Bill White - 10/3/18
 *
 * Description:  Regression GAIN calculation. Uses linear or logistic 
 * regression in calculating main effects and interactions of SNPs and 
 * numeric attributes.
 * =============================================================================
 */

#ifndef __REGAIN_MINIMAL_H__
#define __REGAIN_MINIMAL_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "zfstream.h"
#include "model.h"

#include "Insilico.h"

enum RegainMinimalOutputTransform {
  REGAIN_MINIMAL_OUTPUT_TRANSFORM_NONE, 
  REGAIN_MINIMAL_OUTPUT_TRANSFORM_ABS, 
  REGAIN_MINIMAL_OUTPUT_TRANSFORM_THRESH
};

class RunRecord {
public:
  RunRecord();
  ~RunRecord();
  void Print(ostream& outFile);
  std::vector<std::string> vars;
  vector_t coefs;
  vector_t pvals;
  vector_t stders;
};

class RunInfo {
public:
  RunInfo();
  ~RunInfo();
  void Print(string outFilename);
  std::vector<RunRecord> models;
};

class RegainMinimal {
public:
	RegainMinimal();
	~RegainMinimal();
  void setDefaults();
  // set the value to use when regression procedure fails
  void setFailureValue(double fValue);
  // set output threshold
  bool setOutputThreshold(double threshold);
  // set output type
  bool setOutputTransform(RegainMinimalOutputTransform transform);
  // print output options to stdout and the inbix log
  void logOutputOptions();
  bool logMatrixStats();
	// iterate over all SNPs and numeric attributes (if present), calculating
	// regression for main effect and interaction terms
	void run();
  void getRawMatrix(matrix_t& matrixRef) { matrixRef = regainMatrix; }
  // read a reGAIN file
  bool readRegainMinimalFromFile(std::string regainFilename);
  // write reGAIN matrix to a new file
  bool writeRegainMinimalToFile(std::string newRegainMinimalFilename);
  bool writeRegainMinimalPvalsToFile(std::string newRegainMinimalPvalsFilename);
  // write reGAIN matrix to a new SIF file
  bool writeRegainMinimalToSifFile(std::string newSifFilename);
  bool saveRunInfo(bool paramSaveRunInfo);
  bool writeRegainMinimalRunInfo(std::string newRunInfoFilename);
private:
  Model* createUnivariateModel(uint varIndex, bool varIsNumeric);
  Model* createInteractionModel(uint varIndex1, bool var1IsNumeric,
                                uint varIndex2, bool var2IsNumeric);
  bool fitModelParameters(Model* thisModel, uint thisCoefIdx);
  bool checkValue(std::string coefLabel, double checkVal, double checkPval,
                  double& returnVal, double& returnPVal);
	// calculate main effect regression coefficients for diagonal terms	
	void mainEffect(uint varIndex, bool varIsNumeric);
	// add covariate terms to model, handling multiple covariates
	void addCovariates(Model &m);
	// calculate epistatic interaction between two SNPs or numeric attributes,
	// or a SNP and a numeric attribute
	void interactionEffect(uint varIndex1, bool var1IsNumeric, 
                         uint varIndex2, bool var2IsNumeric);
  bool updateStats();
  void writeFailures();
  void writeWarnings();
  void printFittedModel(Model* thisModel);

  // output options - bcw - 4/30/13
  bool useOutputThreshold;
  double outputThreshold;
  RegainMinimalOutputTransform outputTransform;
  // bcw - 10/18/18
  bool saveRuninfoFlag;
  RunInfo runinfo;
	// num attributes (SNPs + numeric for integrative, SNPs for normal regain)
	uint numAttributes;
  // vector of attribute names of the regain matrix columns
  std::vector<std::string> attributeNames;
	// SIF interaction threshold
	double sifThresh;
	// in memory arrays
  matrix_t regainMatrix;
  matrix_t regainPMatrix;
  // regression warnings - bcw - 4/30/13
  std::vector<std::string> warnings;
  // regression failures - bcw - 5/29/13
  std::vector<std::string> failures;
  // value to use when regression fails
  double failureValue;
  uint nanCount;
  uint infCount;
  // some global regain calculation stats
  double minMainEffect;
  double maxMainEffect;
  double minInteraction;
  double maxInteraction;
};
#endif
