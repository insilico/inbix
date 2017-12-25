/* 
 * File:   DcVar.h
 * Author: bwhite 10/19/17
 */

#ifndef DCVAR_H
#define	DCVAR_H

#include <string>
#include <map>
#include <vector>

#include <armadillo>

#include "plink.h"
// PLINK interface pointer PP defined in Insilico.h
#include "Insilico.h"

// handle both PLINK (Caleb, et al paper) BED/BIM/BAM and 
// OMRF (Courtney Montgomery) separate files
enum SNP_INPUT_TYPE {
  SNP_SRC_PLINK, SNP_SRC_FILE
};

// histone modification site records
enum CHIP_SEQ_EXTRACT_FIELD_IDX {
  CHIP_SEQ_CHROM=0, CHIP_SEQ_POS=1, CHIP_SEQ_EXPR=11, CHIP_SEQ_SNP=15
};

// SNPs/genotypes
struct SNP_INFO {
  std::string chrom;
  uint location;
  char refAllele;
};
typedef std::map<std::string, SNP_INFO> SNP_INFO_MAP;
typedef std::map<std::string, SNP_INFO>::const_iterator SNP_INFO_MAP_IT;

// ChIP-seq expression - histone modification site reads
struct CHIP_SEQ_INFO {
  std::string chrom;
  uint position;
  uint totalRegionReads;
};
typedef std::map<std::string, CHIP_SEQ_INFO> CHIP_SEQ_INFO_MAP;

const string CHECKPOINT_FILENAME = "dcvar.chk";

class DcVar {
public:
  DcVar(SNP_INPUT_TYPE snpInputTypeParam=SNP_SRC_PLINK,
        bool hasChipSeq=false, bool debugFlag=false);
  bool Run(bool debugFlag=false);
  virtual ~DcVar();
private:
  bool RunPlink(bool debugFlag=false);
  bool RunOMRF(bool debugFlag=false);
  bool CheckInputs();
  bool SetDebugMode(bool debugFlag=true);
  void PrintState();
  bool ReadGenotypesFile();
  bool ReadSnpLocationsFile();
  bool ReadGeneExpressionFile();
  bool ReadChipSeqFile();
  bool MapPhenosToModel(std::vector<uint> phenos, std::string varModel,
                        std::vector<uint>& mappedPhenos);
  bool SplitExpressionCaseControl(arma::mat& caseMatrix, arma::mat& ctrlMatrix);
  bool ComputeDifferentialCorrelationZsparse(std::string snp, 
                                            arma::mat& cases, 
                                            arma::mat& ctrls);
  bool ComputeDifferentialCorrelationZ(std::string snp, 
                                       arma::mat& cases, 
                                       arma::mat& ctrls, 
                                       double correctedP);
  bool ComputeDifferentialCorrelationZvals(std::string snp, 
                                           arma::mat& cases, 
                                           arma::mat& ctrls);
  bool FlattenPvals(vector_t& retPvals);
  bool FilterPvalues(uint& numFiltered);
  double CalculateFdrBHThreshold();
  bool WriteCheckpoint(uint snpIndex, std::string snpName);
  bool ReadCheckpoint(std::pair<uint, string>& lastSnp);
  bool WriteResults(std::string filename);
  // INPUTS
  SNP_INPUT_TYPE snpInputType;
  bool chipSeq;
  bool debugFlag;
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
  // ALGORITHM
  double numCombs;
  uint totalTests;
  std::vector<uint> caseIdxCol;
  std::vector<uint> ctrlIdxCol;
  // OUTPUTS
  arma::sp_mat zVals;
  arma::mat pVals;
};

#endif	/* DCVAR_H */
