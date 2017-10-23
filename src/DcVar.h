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

#include "Insilico.h"

// handle PLINK (Caleb, et al paper) and OMRF (Courtney) separate files
enum SNP_INPUT_TYPE {
  SNP_SRC_PLINK, SNP_SRC_FILE
};

enum CHIP_SEQ_EXTRACT_FIELD_IDX {
  CHIP_SEQ_CHROM=0, CHIP_SEQ_POS=1, CHIP_SEQ_EXPR=11, CHIP_SEQ_SNP=15
};

// SNPS/genotypes
typedef std::map<std::string, std::vector<uint>> GENOTYPE_RECORDS;
struct SNP_INFO {
  std::string chrom;
  uint location;
  char refAllele;
};
typedef std::map<std::string, SNP_INFO> SNP_INFO_MAP;

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
        bool hasChipSeq=false, bool useDebug=false);
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
  bool RunPlink(bool debugFlag=false);
  bool RunOMRF(bool debugFlag=false);
  // INPUTS
  bool chipSeq;
  bool debugFlag;
  SNP_INPUT_TYPE snpInputType;
  bool chipSeqMode;
  bool debugMode;
  std::vector<std::string> snpNames;
  std::vector<std::string> genotypeSubjects;
  std::vector<std::string> geneNames;
  std::vector<std::string> geneExprSubjects;
  GENOTYPE_RECORDS genotypes;
  SNP_INFO_MAP snpLocations;
  std::vector<std::vector<double>> exprMatrix;
  CHIP_SEQ_INFO_MAP chipSeqExpression;
  // OUTPUTS
  arma::mat resultsMatrix;
  arma::mat resultsMatrixPvals;
};

#endif	/* DCVAR_H */
