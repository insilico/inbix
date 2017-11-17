/* 
 * File:   DcVar.h
 * Author: bwhite 10/19/17
 */

#ifndef DCVAR_H
#define	DCVAR_H

#include <fstream>
#include <vector>
#include <string>
#include <map>

#include <armadillo>

#include "plink.h"
#include "Insilico.h"

// handle PLINK (Caleb, et al paper) and OMRF (Courtney) separate files
enum SNP_INPUT_TYPE {
  SNP_SRC_PLINK, SNP_SRC_FILE
};

enum CHIP_SEQ_EXTRACT_FIELD_IDX {
  CHIP_SEQ_CHROM=0, CHIP_SEQ_POS=1, CHIP_SEQ_EXPR=11, CHIP_SEQ_SNP=15
};

struct SNP_INFO {
  std::string chrom;
  uint location;
  char refAllele;
};
typedef std::map<std::string, SNP_INFO> SNP_INFO_MAP;
typedef std::map<std::string, SNP_INFO>::const_iterator SNP_INFO_MAP_IT;

// ChIP-seq expression
struct CHIP_SEQ_INFO {
  std::string chrom;
  uint position;
  uint totalRegionReads;
};
typedef std::map<std::string, CHIP_SEQ_INFO> CHIP_SEQ_INFO_MAP;

class DcVar {
public:
  DcVar(SNP_INPUT_TYPE snpInputTypeParam=SNP_SRC_PLINK,
        bool hasChipSeq=false, bool debugFlag=false);
  bool Run(bool debugFlag=false);
  virtual ~DcVar();
private:
  bool CheckInputs();
  bool SetDebugMode(bool debugFlag=true);
  void PrintState();
  bool ReadGenotypesFile();
  bool ReadSnpLocationsFile();
  bool ReadGeneExpressionFile();
  bool ReadChipSeqFile();
  bool MapPhenosToModel(std::vector<uint> phenos, std::string varModel);
  bool SplitExpressionCaseControl(arma::mat& caseMatrix, arma::mat& ctrlMatrix);
  bool ComputeDifferentialCorrelationZnaive(std::string variant, 
                                            arma::mat& X, 
                                            arma::mat& Y,
                                            arma::sp_mat& zVals);
  bool ComputeDifferentialCorrelationZ(std::string variant, 
                                       arma::mat& X, 
                                       arma::mat& Y, 
                                       double correctedP);
  bool ComputeDifferentialCorrelationZvals(std::string variant, 
                                           arma::mat& X, 
                                           arma::mat& Y, 
                                           double correctedP,
                                           arma::sp_mat& zVals,
                                           arma::sp_mat& pVals);
  bool RunPlink(bool debugFlag=false);
  bool RunOMRF(bool debugFlag=false);
  // INPUTS
  bool chipSeq;
  bool debugFlag;
  SNP_INPUT_TYPE snpInputType;
  bool chipSeqMode;
  bool debugMode;
  std::vector<std::string> snpNames;
  SNP_INFO_MAP snpLocations;
  std::vector<std::string> genotypeSubjects;
  matrix_t genotypeMatrix;
  std::vector<std::string> geneExprNames;
  std::vector<std::string> geneExprSubjects;
  matrix_t expressionMatrix;
  CHIP_SEQ_INFO_MAP chipSeqExpression;
  // ALGORITHM VARIABLES
  std::vector<uint> mappedPhenos;
  std::vector<uint> caseIdxCol;
  std::vector<uint> ctrlIdxCol;
};

#endif	/* DCVAR_H */
