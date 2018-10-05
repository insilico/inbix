/* =============================================================================
 * Filename: RegainMinimal.cpp - Bill White - 10/3/18
 *
 * Description:  Regression GAIN calculation. Uses linear or logistic 
 * regression in calculating main effects and interactions of SNPs and 
 * numeric attributes.
 * =============================================================================
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <sstream>

#include <omp.h>

// PLINK functions and data types
#include "plink.h"
#include "options.h"
#include "model.h"
#include "logistic.h"
#include "linear.h"
#include "stats.h"
#include "helper.h"
// inbix constants and PLINK point PP
#include "RegainMinimal.h"
#include "Insilico.h"

RegainMinimal::RegainMinimal() {
  // defaults
  SetDefaults();

  numAttributes = PP->nl_all + PP->nlistname.size();
  PP->printLOG("Total number of attributes [ " + int2str(numAttributes) + " ]\n");
  
  string coefLabel = "";
  for(uint i=0; i < PP->locus.size(); ++i) {
    coefLabel = PP->locus[i]->name;
    attributeNames.push_back(coefLabel);
  }
  for(uint i=0; i < PP->nlistname.size(); ++i) {
    coefLabel = PP->nlistname[i];
    attributeNames.push_back(coefLabel);
  }
  
  sizeMatrix(regainMatrix, numAttributes, numAttributes);
  sizeMatrix(regainPMatrix, numAttributes, numAttributes);
  
  if(par::regainMatrixTransform == "abs") {
    setOutputTransform(REGAIN_MINIMAL_OUTPUT_TRANSFORM_ABS);
  } else {
    if(par::regainMatrixTransform == "thresh") {
      setOutputTransform(REGAIN_MINIMAL_OUTPUT_TRANSFORM_THRESH);
    } else {
      setOutputTransform(REGAIN_MINIMAL_OUTPUT_TRANSFORM_NONE);
    }
  }
}

RegainMinimal::~RegainMinimal() {
}

void RegainMinimal::SetDefaults() {
  useOutputThreshold = par::regainMatrixThreshold;
  outputThreshold = par::regainMatrixThresholdValue;
  failureValue = par::regainFailValue;
  nanCount = 0;
  infCount = 0;
  minMainEffect = 0;
  maxMainEffect = 0;
  minInteraction = 0;
  maxInteraction = 0;
}

void RegainMinimal::setFailureValue(double fValue) {
  failureValue = fValue;
}

bool RegainMinimal::setOutputThreshold(double threshold) {
  useOutputThreshold = true;
  outputThreshold = threshold;
  return true;
}

bool RegainMinimal::setOutputTransform(RegainMinimalOutputTransform transform) {
  outputTransform = transform;
  return true;
}

void RegainMinimal::logOutputOptions() {
  stringstream ss;
  switch(outputTransform) {
    case REGAIN_MINIMAL_OUTPUT_TRANSFORM_ABS:
      PP->printLOG("Output transform [ absolute value ]\n");
      break;
    case REGAIN_MINIMAL_OUTPUT_TRANSFORM_NONE:
      break;
    case REGAIN_MINIMAL_OUTPUT_TRANSFORM_THRESH:
      ss << "Output transform [ threshold values < " << outputThreshold
        << " => 0 ]\n";
      PP->printLOG(ss.str());
      break;
  }
  PP->printLOG("Regression failure substitution value [ " +
    dbl2str(failureValue) + " ]\n");
}

void RegainMinimal::run() {
  // reset the warnings list
  warnings.clear();
  failures.clear();
  
  // all main effects
  PP->printLOG("Run all main effects models\n");
  #pragma omp parallel for
  for(int k=0; k < numAttributes; k++) {
    mainEffect(k, k >= PP->nl_all);
  } // Next SNP/numeric attribute
  
  PP->printLOG("Run all interaction effects models\n");
  #pragma omp parallel for
  for(int k=0; k < numAttributes * (numAttributes + 1) / 2; k++) {
    uint varIndex1 = k / (numAttributes + 1);
    uint varIndex2 = k % (numAttributes + 1);
    if(varIndex2 > varIndex1) {
      varIndex1 = numAttributes - varIndex1 - 1;
      varIndex2 = numAttributes - varIndex2;
    }
    interactionEffect(varIndex1, varIndex1 >= PP->nl_all,
                      varIndex2, varIndex2 >= PP->nl_all);
  } // Next pair of SNPs/numeric attributes
  
  writeWarnings();
  writeFailures();
  if(nanCount) {
    PP->printLOG("Detected [ " + int2str(nanCount) + " ] NaN's\n");
  }
  if(infCount) {
    PP->printLOG("Detected [ " + int2str(infCount) + " ] Inf's\n");
  }
}

Model* RegainMinimal::createUnivariateModel(uint varIndex, bool varIsNumeric) {
  Model* thisModel = 0;
  // logistic regression for binary phenotypes (traits), linear otherwise
  if(par::bt) {
    LogisticModel* m = new LogisticModel(PP);
    thisModel = m;
  } else {
    LinearModel* m = new LinearModel(PP);
    thisModel = m;
  }
  // Set missing data
  thisModel->setMissing();
  // label for regression model
  string coefLabel = "";
  if(varIsNumeric) {
    coefLabel = PP->nlistname[varIndex - PP->nl_all];
  } else {
    coefLabel = PP->locus[varIndex]->name;
  }
  // Main effect of SNP/numeric attribute
  if(varIsNumeric) {
    thisModel->addNumeric(varIndex - PP->nl_all);
  } else {
    thisModel->addAdditiveSNP(varIndex);
  }
  thisModel->label.push_back(coefLabel);
  thisModel->testParameter = 1; // single variable main effect
  if(par::covar_file) {
    addCovariates(*thisModel);
  }

  return thisModel;
}

Model* RegainMinimal::createInteractionModel(uint varIndex1, bool var1IsNumeric, 
                                             uint varIndex2, bool var2IsNumeric) {
  // attempt to fit a model and retrieve the estimated parameters
  Model* thisModel = 0;
  // logistic regression for binary phenotypes (traits), linear otherwise
  if(par::bt) {
    LogisticModel* m = new LogisticModel(PP);
    thisModel = m;
  } else {
    LinearModel* m = new LinearModel(PP);
    thisModel = m;
  }
  // Set missing data
  thisModel->setMissing();
  // labels in regression model
  string coef1Label = "";
  if(var1IsNumeric) {
    coef1Label = PP->nlistname[varIndex1 - PP->nl_all];
  } else {
    coef1Label = PP->locus[varIndex1]->name;
  }
  string coef2Label = "";
  if(var2IsNumeric) {
    coef2Label = PP->nlistname[varIndex2 - PP->nl_all];
  } else {
    coef2Label = PP->locus[varIndex2]->name;
  }
  // Main effect of SNP/numeric attribute 1
  if(var1IsNumeric) {
    thisModel->addNumeric(varIndex1 - PP->nl_all);
  } else {
    thisModel->addAdditiveSNP(varIndex1);
  }
  thisModel->label.push_back(coef1Label);
  // Main effect of SNP/numeric attribute 2
  if(var2IsNumeric) {
    thisModel->addNumeric(varIndex2 - PP->nl_all);
  } else {
    thisModel->addAdditiveSNP(varIndex2);
  }
  thisModel->label.push_back(coef2Label);
  thisModel->addInteraction(1, 2);
  thisModel->label.push_back("EPI");
  thisModel->testParameter = 3;
  // add covariates if specified
  if(par::covar_file) {
    addCovariates(*thisModel);
    thisModel->testParameter += par::clist_number;
  }

  return thisModel;
}

bool RegainMinimal::fitModelParameters(Model* thisModel, uint thisCoefIdx) {
  bool success = true;
  // build design matrix and fit the model parameters
  thisModel->buildDesignMatrix();
  thisModel->fitLM();
  thisModel->testParameter <- thisCoefIdx;
  // Was the model fitting method successful?
  if(!thisModel->isValid()) {
    string failMsg = "WARNING: Invalid regression fit: ";
    RegressionInvalidType invalidReason = 
            thisModel->getRegressionFailureType();
    switch(invalidReason) {
      case REGRESSION_INVALID_NONE:
        failMsg += "Error code REGRESSION_INVALID_NONE detected";
        break;
      case REGRESSION_INVALID_SVDINV:
        failMsg += "SVD inverse failed";
        break;
      case REGRESSION_INVALID_EMPTY:
        failMsg += "Empty model, either individuals or parameters";
        break;
      case REGRESSION_INVALID_MULTICOLL:
        failMsg += "Possible multicollinearity";
        break;
      case REGRESSION_INVALID_VIF:
        failMsg += "VIF check failed";
        break;
      case REGRESSION_INVALID_LINHYPOTH:
        failMsg += "Linear model hypothesis failed";
        break;
      default:
        failMsg += "Regression invalid failure type detected: " + int2str(invalidReason);
        break;
    }
    failures.push_back(failMsg);
    success = false;
  }
  if(!thisModel->fitConverged()) {
    failures.push_back("Model failed to converge");
    success = false;
  }  
  return success;
}

bool RegainMinimal::checkValue(string coefLabel, 
                               double checkVal, double checkPval, 
                               double& returnVal, double& returnPval) {
  returnVal = checkVal;
  returnPval = checkPval;
  bool useFailureValue = false;
  // report large p-value of coefficient as a warning
  if(checkPval > par::regainLargeCoefPvalue) {
    stringstream ss;
    ss << "Large p-value [" << checkPval
      << "] on coefficient for variable [" << coefLabel << "]";
    warnings.push_back(ss.str());
  }
  if(std::isinf(checkVal)) {
    useFailureValue = true;
    ++infCount;
    stringstream ss;
    ss << "Regression test statistic is +/-infinity on coefficient "
      << "for interaction variable [ " << coefLabel << " ]\n";
    warnings.push_back(ss.str());
  } 
  if(std::isnan(checkVal)) {
    useFailureValue = true;
    ++nanCount;
    stringstream ss;
    ss << "Regression test statistic is not a number NaN on coefficient "
      << "for interaction variable [ " << coefLabel << "] \n";
    warnings.push_back(ss.str());
  }
  if(useFailureValue) {
    returnVal = failureValue;
    returnPval = 1.0;
  }
  
  return useFailureValue;  
}

void RegainMinimal::mainEffect(uint varIndex, bool varIsNumeric) {
  // setup regression model
  Model *thisModel = createUnivariateModel(varIndex, varIsNumeric);
  // attempt to fit a model and retrieve the estimated parameters
  double newVal = 0;
  double newPval = 1.0;
  // get the estimated parameters
  vector_t betaMainEffectCoefs;
  double mainEffectValue = failureValue;
  vector_t betaMainEffectCoefPvals;
  double mainEffectPval = 1.0;
  vector_t mainEffectModelSE;
  if(fitModelParameters(thisModel, 1)) {
    // Obtain estimates and statistics
    betaMainEffectCoefs = thisModel->getCoefs();
    // p-values don't include intercept term
    betaMainEffectCoefPvals = thisModel->getPVals();
    mainEffectPval = betaMainEffectCoefPvals[thisModel->testParameter - 1];
    mainEffectModelSE = thisModel->getSE();
    if(par::regainUseBetaValues) {
      mainEffectValue = betaMainEffectCoefs[thisModel->testParameter];
    } else {
      mainEffectValue = betaMainEffectCoefs[thisModel->testParameter] /
        mainEffectModelSE[thisModel->testParameter];
    }
    checkValue("", mainEffectValue, mainEffectPval, newVal, newPval);
  } 
  #pragma omp critical
  {
    // update the matrices
    if(outputTransform == REGAIN_MINIMAL_OUTPUT_TRANSFORM_ABS) {
      newVal = fabs(newVal);
    }
    regainMatrix[varIndex][varIndex] = newVal;
    regainPMatrix[varIndex][varIndex] = newPval;
    if(par::do_regain_pvalue_threshold) {
      if(newPval > par::regainPvalueThreshold) {
        regainMatrix[varIndex][varIndex] = 0;
        regainPMatrix[varIndex][varIndex] = 1;
      }
    }
  }

  // free model memory
  delete thisModel;
}

void RegainMinimal::addCovariates(Model &m) {
  for(uint i = 0; i < par::clist_number; i++) {
    // add covariate to the model
    m.addCovariate(i);
    m.label.push_back(PP->clistname[i]);
  }
}

void RegainMinimal::interactionEffect(uint varIndex1, bool var1IsNumeric,
                                      uint varIndex2, bool var2IsNumeric) {
  // labels in regression model
  string coef1Label = "";
  if(var1IsNumeric) {
    coef1Label = PP->nlistname[varIndex1 - PP->nl_all];
  } else {
    coef1Label = PP->locus[varIndex1]->name;
  }
  string coef2Label = "";
  if(var2IsNumeric) {
    coef2Label = PP->nlistname[varIndex2 - PP->nl_all];
  } else {
    coef2Label = PP->locus[varIndex2]->name;
  }
  Model *thisModel = createInteractionModel(varIndex1, var1IsNumeric,
                                            varIndex2, var2IsNumeric);
  // fit the model and get the estimated parameters
  vector_t betaInteractionCoefs;
  vector_t betaInteractionCoefPVals;
  double interactionVal = failureValue;
  double interactionPval = 1.0;
  vector_t interactionModelSE;
  vector_t::const_iterator bIt;
  vector_t::const_iterator sIt;
  double newVal = 0;
  double newPval = 1.0;
  if(fitModelParameters(thisModel, 3)) {
    // model converged, so get the estimated parameters and statistics
    betaInteractionCoefs = thisModel->getCoefs();
    interactionVal = betaInteractionCoefs[betaInteractionCoefs.size() - 1];
    betaInteractionCoefPVals = thisModel->getPVals();
    interactionPval = 
            betaInteractionCoefPVals[betaInteractionCoefPVals.size() - 1];
    interactionModelSE = thisModel->getSE();
    // calculate statistical test value from beta/SE (t-test or z-test)
    if(par::regainUseBetaValues) {
      interactionVal = betaInteractionCoefs[thisModel->testParameter];
    } else {
      interactionVal = betaInteractionCoefs[thisModel->testParameter] /
        interactionModelSE[thisModel->testParameter];
    }
    checkValue("", interactionVal, interactionPval, newVal, newPval);  
  } 

#pragma omp critical
  {
    if(outputTransform == REGAIN_MINIMAL_OUTPUT_TRANSFORM_ABS) {
      newVal = fabs(newVal);
    }
    regainMatrix[varIndex1][varIndex2] = newVal;
    regainMatrix[varIndex2][varIndex1] = newVal;
    if(par::do_regain_pvalue_threshold) {
      if(newPval > par::regainPvalueThreshold) {
        regainMatrix[varIndex1][varIndex2] = 0;
        regainMatrix[varIndex2][varIndex1] = 0;
      }
    }
    regainPMatrix[varIndex1][varIndex2] = newPval;
    regainPMatrix[varIndex2][varIndex1] = newPval;
  }
  
  // free model memory
  delete thisModel;
}

bool RegainMinimal::updateStats() {
  minMainEffect = maxMainEffect = regainMatrix[0][0];
  minInteraction = maxInteraction = regainMatrix[0][1];
  for(uint i = 0; i < numAttributes; ++i) {
    for(uint j = i; j < numAttributes; ++j) {
      if(i == j) {
        if(regainMatrix[i][j] < minMainEffect) {
          minMainEffect = regainMatrix[i][j];
        }
        if(regainMatrix[i][j] > maxMainEffect) {
          maxMainEffect = regainMatrix[i][j];
        }
      } else {
        if(regainMatrix[i][j] < minInteraction) {
          minInteraction = regainMatrix[i][j];
        }
        if(regainMatrix[i][j] > maxInteraction) {
          maxInteraction = regainMatrix[i][j];
        }
      }
    }
  }

  return true;
}

bool RegainMinimal::logMatrixStats() {
  updateStats();

  PP->printLOG("reGAIN matrix statistics:\n");
  PP->printLOG("minimum main effect [ " + dbl2str(minMainEffect) + " ]\n");
  PP->printLOG("maximum main effect [ " + dbl2str(maxMainEffect) + " ]\n");
  PP->printLOG("minimum interaction [ " + dbl2str(minInteraction) + " ]\n");
  PP->printLOG("maximum interaction [ " + dbl2str(maxInteraction) + " ]\n");

  return true;
}

void RegainMinimal::writeFailures() {
  if(failures.size()) {
    double numCombinations = static_cast<double>(numAttributes * numAttributes);
    double numFailures = (double) failures.size();
    double percentFailures = (numFailures / numCombinations) * 100.0;
    PP->printLOG(dbl2str(numFailures) + " failures in " + 
      dbl2str(numCombinations)+ " regression models "
      + dbl2str(percentFailures) + "%\n");
    string failureFilename = par::output_file_name + ".regression.failures";
    PP->printLOG("Writing failure messages to [ " + failureFilename + " ]\n");
    ofstream failureFile(failureFilename);
    for(vector<string>::const_iterator fIt = failures.begin();
      fIt != failures.end(); ++fIt) {
      failureFile << *fIt << endl;
    }
    failureFile.close();
  }
}
  
void RegainMinimal::writeWarnings() {
  // report warnings to stdout - bcw - 4/30/13
  if(warnings.size()) {
    double numCombinations = (numAttributes * (numAttributes - 1)) / 2.0;
    double numFailures = (double) warnings.size();
    double percentFailures = (numFailures / (numCombinations + numAttributes)) * 100.0;
    PP->printLOG(dbl2str(numFailures) + " warnings in " + 
      dbl2str(numCombinations + numAttributes)+ " regression models "
      + dbl2str(percentFailures) + "%\n");
    string failureFilename = par::output_file_name + ".regression.warnings";
    PP->printLOG("Writing warning messages to [ " + failureFilename + " ]\n");
    ofstream failureFile(failureFilename);
    for(vector<string>::const_iterator wIt = warnings.begin();
      wIt != warnings.end(); ++wIt) {
      failureFile << *wIt << endl;
    }
    failureFile.close();
  }
}

bool RegainMinimal::readRegainMinimalFromFile(string regainFilename) {
  checkFileExists(regainFilename);
  ifstream REGAIN_MINIMAL(regainFilename, ios::in);
  if(REGAIN_MINIMAL.fail()) {
    return false;
  }
  bool readHeader = false;
  uint matrixRow = 0;
  uint matrixCol = 0;
  double rawValue = 0;
  while(!REGAIN_MINIMAL.eof()) {
    char nline[par::MAX_LINE_LENGTH];
    REGAIN_MINIMAL.getline(nline, par::MAX_LINE_LENGTH, '\n');
    // convert to string
    string sline = nline;
    if(sline == "") continue;
    // read line from text file into a vector of tokens
    string buf;
    stringstream ss(sline);
    vector<string> tokens;
    while(ss >> buf) {
      tokens.push_back(buf);
    }
    if(!readHeader) {
      // process reGAIN file header
      numAttributes = tokens.size();
      readHeader = true;
      sizeMatrix(regainMatrix, numAttributes, numAttributes);
      for(uint i = 0; i < numAttributes; ++i) {
        attributeNames.push_back(tokens[i]);
      }
      continue;
    } else {
      matrixCol = 0;
      for(uint tokenIdx = 0; tokenIdx < tokens.size(); ++tokenIdx, ++matrixCol) {
        if(!from_string<double>(rawValue, tokens[tokenIdx], std::dec)) {
          PP->printLOG("Error parsing token:" + tokens[tokenIdx] + "\n");
          return false;
        }
        regainMatrix[matrixRow][matrixCol] = rawValue;
        regainMatrix[matrixCol][matrixRow] = rawValue;
      }
    }
    ++matrixRow;
  }
  REGAIN_MINIMAL.close();

  return true;
}

bool RegainMinimal::writeRegainMinimalToFile(string newRegainMinimalFilename) {
  PP->printLOG("Writing REGAIN_MINIMAL matrix [ " + newRegainMinimalFilename + " ]\n");
  ofstream outFile(newRegainMinimalFilename);
  if(outFile.fail()) {
    return false;
  }
  // write header
  uint hIdx = 0;
  for(vector<string>::const_iterator hIt = attributeNames.begin();
    hIt != attributeNames.end(); hIt++, hIdx++) {
    if(hIdx) {
      outFile << "\t" << *hIt;
    } else {
      outFile << *hIt;
    }
  }
  outFile << endl;
  for(uint i = 0; i < numAttributes; ++i) {
    for(uint j = 0; j < numAttributes; ++j) {
      if(j) {
        outFile << "\t" << regainMatrix[i][j];
      } else {
        outFile << regainMatrix[i][j];
      }
    }
    outFile << endl;
  }
  outFile.close();

  return true;
}

bool RegainMinimal::writeRegainMinimalPvalsToFile(string newRegainMinimalPvalsFilename) {
  PP->printLOG("Writing REGAIN_MINIMAL p-values matrix [ " + newRegainMinimalPvalsFilename + " ]\n");
  ofstream outFile(newRegainMinimalPvalsFilename);
  if(outFile.fail()) {
    return false;
  }
  // write header
  uint hIdx = 0;
  for(vector<string>::const_iterator hIt = attributeNames.begin();
    hIt != attributeNames.end(); ++hIt, ++hIdx) {
    if(hIdx) {
      outFile << "\t" << *hIt;
    } else {
      outFile << *hIt;
    }
  }
  outFile << endl;
  for(uint i = 0; i < numAttributes; ++i) {
    for(uint j = 0; j < numAttributes; ++j) {
      if(j) {
        outFile << "\t" << regainPMatrix[i][j];
      } else {
        outFile << regainPMatrix[i][j];
      }
    }
    outFile << endl;
  }
  outFile.close();

  return true;
}

bool RegainMinimal::writeRegainMinimalToSifFile(string newSifFilename) {
  PP->printLOG("Writing REGAIN_MINIMAL matrix to SIF [ " + newSifFilename + " ]\n");
  ofstream outFile(newSifFilename);
  if(outFile.fail()) {
    return false;
  }
  for(uint i = 0; i < numAttributes; ++i) {
    for(uint j = i + 1; j < numAttributes; ++j) {
      string attr1 = attributeNames[i];
      string attr2 = attributeNames[j];
      double valueToWrite = regainMatrix[i][j];
      if(valueToWrite > outputThreshold) {
        outFile << attr1 << "\t" << valueToWrite << "\t" << attr2 << endl;
      }
    }
  }
  outFile.close();

  return true;
}
