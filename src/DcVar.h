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

// handle PLINK (Caleb) and OMRF (Courtney) separate files
enum SNP_INPUT_TYPE {
  SNP_SRC_PLINK, SNP_SRC_FILE
};

// SNPS/genotypes
typedef std::map<std::string, std::vector<double>> GENOTYPE_RECORDS;
struct SNP_LOCATIONS {
  std::string chrom;
  uint location;
  char refAllele;
};
typedef std::map<std::string, SNP_LOCATIONS> SNP_INFO;

class DcVar {
public:
  DcVar(SNP_INPUT_TYPE snpInputTypeParam=SNP_SRC_PLINK, 
        bool chipSeq=false,
        bool debugFlag=false);
  virtual ~DcVar();
  bool SetDebugMode(bool debugFlag=true);
  void PrintState();
  bool ReadGenotypesFile();
  bool ReadGenotypeLocationsFile();
  bool ReadChipSeqFile();
private:
  bool CheckInputs();
  // INPUTS
  SNP_INPUT_TYPE snpInputType;
  bool chipSeqMode;
  bool debugMode;
  std::vector<std::string> snpNames;
  GENOTYPE_RECORDS genotypes;
  SNP_INFO snpLocations;
  std::vector<std::string> chipSeqNames;
  std::vector<double> chipSeqExpression;
  // OUTPUTS
  arma::mat resultsMatrix;
  arma::mat resultsMatrixPvals;
};

#endif	/* DCVAR_H */
