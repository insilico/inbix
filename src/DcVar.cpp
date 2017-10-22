/* 
 * File:   DcVar.cpp
 * Author: bwhite 10/19/17
 */

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <fstream>
// Armadillo Linear Algebra/Matrices/Vectors
#include <armadillo>
// PLINK
#include "options.h"
#include "plink.h"
#include "helper.h"
#include "zed.h"
// inbix
#include "DcVar.h"
#include "Insilico.h"
#include "StringUtils.h"

using namespace std;
using namespace insilico;

DcVar::DcVar(SNP_INPUT_TYPE snpInputTypeParam, bool hasChipSeq, bool useDebug) {
  chipSeq = hasChipSeq;
  debugMode = useDebug;
  if(snpInputTypeParam == SNP_SRC_FILE) {
    ReadGenotypesFile();
    ReadGenotypeLocationsFile();
  }
  if(!ReadGeneExpressionFile()) {
    cerr << "Cannot read gene expression. Exiting." << endl;
    exit(1);
  }
  chipSeqMode = hasChipSeq;
  if(chipSeqMode) {
    ReadChipSeqFile();
  }
  SetDebugMode(useDebug);
  if(!CheckInputs()) {
    cerr << "Inputs check failed. Exiting." << endl;
    exit(1);
  }
}

bool DcVar::Run(bool debugFlag) {
  SetDedugMode(debugFlag);
  return true;
}

DcVar::~DcVar() {
}

bool DcVar::CheckInputs() {
  // SNPs/genotype
  uint numSnp = snpLocations.size();
  // expression
  uint nunmExpr = geneNames.size();
  // chipseq
  uint numChipSeq = 0;
  if(chipSeq) {
    numChipSeq = chipSeqExpression.size();
  }
  
  return true;
}

// TODO: remove this? what was the thinking? public interface? delete?
bool DcVar::SetDebugMode(bool debugFlag) {
  debugMode = debugFlag;
  return true;
}

void DcVar::PrintState() {
  PP->printLOG("-----------------------------------------------------------\n");
  string debugFlag = debugMode? "on": "off";
  PP->printLOG("debug mode:                     " + debugFlag + "\n");
  PP->printLOG("SNPs file:                      " + par::dcvar_genotypes_file + "\n");
  PP->printLOG("SNP locations file:             " + par::dcvar_genotypes_locations_file + "\n");
  PP->printLOG("CHiP-seq expression file:       " + par::dcvar_chip_seq_file + "\n");
  PP->printLOG("p-value adjust method:          " + par::dcvar_pfilter_type + "\n");
  PP->printLOG("p-value cutoff for file output: " + dbl2str(par::dcvar_pfilter_value) + "\n");
  PP->printLOG("-----------------------------------------------------------\n");
}

bool DcVar::ReadGenotypesFile() {
  checkFileExists(par::dcvar_genotypes_file);
  PP->printLOG("Reading genotypes input from [ " + par::dcvar_genotypes_file + " ]\n");
  ZInput zin(par::dcvar_genotypes_file, compressed(par::dcvar_genotypes_file));
  uint lineCounter = 0;
  while(!zin.endOfFile()) {
    ++lineCounter;
	  vector<string> tok = zin.tokenizeLine();
    if(tok.size() < 2) {
      cerr << "Error reading line [ " << lineCounter 
              << " ] from " << par::dcvar_genotypes_file << endl;
      return(false);
    }
    vector<uint> lineGenotypes(tok.size() - 1);
	  for(int j=1; j < tok.size(); j++) {
      lineGenotypes[j - 1] = stoi(tok[j]);
	  }
    genotypes[tok[0]] = lineGenotypes;
	}
  zin.close();
  PP->printLOG("Read subject genotypes for " + int2str(lineCounter) + " SNPs\n");
  
  return true;
}

bool DcVar::ReadGenotypeLocationsFile() {
  checkFileExists(par::dcvar_genotypes_locations_file);
  PP->printLOG("Reading genotypes input from [ " + par::dcvar_genotypes_locations_file + " ]\n");
  ZInput zin(par::dcvar_genotypes_file, compressed(par::dcvar_genotypes_locations_file));
  uint lineCounter = 0;
  while(!zin.endOfFile()) {
    ++lineCounter;
	  vector<string> tok = zin.tokenizeLine();
    if(tok.size() != 5) {
      cerr << "Error reading line [ " << lineCounter 
              << " ] from " << par::dcvar_genotypes_locations_file 
              << " should have 5 columns, found " << tok.size()
              << endl;
      return(false);
    }
    SNP_INFO thisSnpInfo;
    thisSnpInfo.chrom = tok[1];
    thisSnpInfo.location = stoi(tok[2]);
    thisSnpInfo.refAllele = tok[3][0];
    snpLocations[tok[0]] = thisSnpInfo;
	}
  zin.close();
  PP->printLOG("Read subject SNP info for " + int2str(lineCounter) + " SNPs\n");

  return true;
}

bool DcVar::ReadGeneExpressionFile() {
  checkFileExists(par::dcvar_gene_expression_file);
  PP->printLOG("Reading gene expression input from [ " + par::dcvar_gene_expression_file + " ]\n");
  ifstream exprFile(par::dcvar_gene_expression_file);
  string header;
  getline(exprFile, header);
  vector<string> headerParts;
  split(headerParts, header, "\t");
  for(uint i=1; i < headerParts.size(); ++i) {
    geneExprSubjects.push_back(headerParts[i]);
  }
  string line;
  uint lineCounter = 0;
  while(getline(exprFile, line)) {
    ++lineCounter;
	  vector<string> tok;
    split(tok, line, "\t");
    if(tok.size() < 2) {
      cerr << "Error reading line [ " << lineCounter 
              << " ] from " << par::dcvar_gene_expression_file 
              << " should have more than 2 columns (subjects)"
              << endl;
      return(false);
    }
    geneNames.push_back(tok[0]);
    vector<double> newRow;
    for(uint i=1; i < tok.size(); ++i) {
      newRow.push_back(stod(tok[i]));
    }
    exprMatrix.push_back(newRow);
	}
  exprFile.close();
  PP->printLOG("Read gene expression for " + int2str(lineCounter) + " SNPs\n");
  
  return true;
}

bool DcVar::ReadChipSeqFile() {
  checkFileExists(par::dcvar_chip_seq_file);
  PP->printLOG("Reading ChIP-seq input from [ " + par::dcvar_chip_seq_file + " ]\n");
  ifstream chipSeqFile(par::dcvar_chip_seq_file);
  uint lineCounter = 0;
  string line;
  while(getline(chipSeqFile, line)) {
    ++lineCounter;
	  vector<string> tok;
    split(tok, line, "\t");
    if(tok.size() != (CHIP_SEQ_SNP + 1)) {
      cerr << "Error reading line [ " << lineCounter 
              << " ] from " << par::dcvar_chip_seq_file 
              << " should have 16 columns, found " << tok.size()
              << endl;
      return(false);
    }
    CHIP_SEQ_INFO thisChipSeqInfo;
    thisChipSeqInfo.chrom = tok[CHIP_SEQ_CHROM];
    thisChipSeqInfo.position = stoi(tok[CHIP_SEQ_POS]);
    thisChipSeqInfo.totalRegionReads = stod(tok[CHIP_SEQ_EXPR]);
    // rs28469609:38367404:C:T
    vector<string> rsnumParts;
    split(rsnumParts, tok[CHIP_SEQ_SNP], ":");
    chipSeqExpression[rsnumParts[0]] = thisChipSeqInfo;
	}
  chipSeqFile.close();
  PP->printLOG("Read ChIP-seq expression for " + int2str(lineCounter) + " SNPs\n");
  
  return true;
}
