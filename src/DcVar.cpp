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

#include <armadillo>
#include "plink.h"
#include "helper.h"

#include "options.h"
#include "DcVar.h"
#include "Insilico.h"

using namespace std;

DcVar::DcVar(SNP_INPUT_TYPE snpInputTypeParam, 
             bool chipSeq,
             bool debugFlag) {
  snpInputType = snpInputTypeParam;
  debugMode = debugFlag;
  if(snpInputType == SNP_SRC_FILE) {
    ReadGenotypesFile();
    ReadGenotypeLocationsFile();
  }
  if(chipSeq) {
    chipSeqMode = chipSeq;
    ReadChipSeqFile();
  }
  if(!CheckInputs()) {
    cerr << "Inputs check failed. Exiting." << endl;
    exit(1);
  }
}

DcVar::~DcVar() {
}

bool DcVar::CheckInputs() {
  // SNPs/genotype
  // expression
  // chipseq
  return true;
}

bool DcVar::SetDebugMode(bool debugFlag) {
  debugMode = debugFlag;
  return true;
}

void DcVar::PrintState() {
  PP->printLOG("-----------------------------------------------------------\n");
  string debugFlag = debugMode? "on": "off";
  PP->printLOG("debug mode: " + debugFlag + "\n");
  PP->printLOG("SNPs file: " + par::dcvar_genotypes_file + "\n");
  PP->printLOG("SNP locations file: " + par::dcvar_genotypes_locations_file + "\n");
  PP->printLOG("Expression file: " + par::dcvar_genotypes_locations_file + "\n");
  PP->printLOG("Expression file: " + par::dcvar_methylation_file + "\n");
  PP->printLOG("p-value adjust method: " + par::dcvar_pfilter_type + "\n");
  PP->printLOG("p-value: " + dbl2str(par::dcvar_pfilter_value) + "\n");
  PP->printLOG("p-value cutoff for file output: " + dbl2str(par::dcvar_pfilter_value) + "\n");
  PP->printLOG("-----------------------------------------------------------\n");
}

bool DcVar::ReadGenotypesFile() {
  return true;
}

bool DcVar::ReadGenotypeLocationsFile() {
  return true;
}

bool DcVar::ReadChipSeqFile() {
  return true;
}
