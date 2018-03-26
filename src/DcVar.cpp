/* 
 * File:   DcVar.cpp
 * Author: bwhite 10/19/17
 */

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>

// Armadillo Linear Algebra/Matrices/Vectors
#include <armadillo>

// Boost lexical cast
#include <boost/lexical_cast.hpp>

// PLINK
#include "options.h"
#include "plink.h"
#include "helper.h"
#include "zed.h"
#include "stats.h"

// inbix
#include "DcVar.h"
#include "Insilico.h"
#include "StringUtils.h"
#include "ArmadilloFuncs.h"

using namespace std;
using namespace insilico;
using namespace boost;
using namespace arma;

bool pvalComparatorAscending(const matrixElement& l, const matrixElement& r) {
  return l.first < r.first;
}

bool snpinfoComparatorAscending(const SNP_INFO& l, const SNP_INFO& r) {
  return l.position < r.position;
}

bool chipseqComparatorAscending(const CHIP_SEQ_INFO& l, const CHIP_SEQ_INFO& r) {
  string chromStr1(l.chrom.begin() + 3, l.chrom.end());
  uint chromNum1 = lexical_cast<uint>(chromStr1);
  string chromStr2(r.chrom.begin() + 3, r.chrom.end());
  uint chromNum2 = lexical_cast<uint>(chromStr2);
  return chromNum1 < chromNum2;
}

DcVar::DcVar(SNP_INPUT_TYPE snpInputTypeParam) {
  PP->printLOG("dcVar initializing\n");
  snpInputType = snpInputTypeParam;
  chipSeqMode = par::do_dcvar_chipseq;
  if(snpInputTypeParam == SNP_SRC_FILE) {
    if(!ReadGeneExpressionFile()) {
      error("Reading gene expression file failed. Exiting.");
    }
    if(chipSeqMode) {
      if(!ReadChipSeqFile()) {
        error("Reading ChIP-seq file failed. Exiting.");
      }
    }
    PP->printLOG("Using separate tab-delimited files for input data sets\n");
  } else {
    // assume PLINK data structures accessible through PP pointer
    PP->printLOG("Using PLINK files for input data sets\n");
  }
  if(!CheckInputs()) {
    error("Checking data sets compatability failed. Exiting.");
  }
  radius = par::dcvar_radius;
  // OpenMP parallelization
  uint numThreads = omp_get_num_threads();
  uint numProcs = omp_get_num_procs();
  PP->printLOG("OpenMP: " + int2str(numThreads) + " threads available\n");
  PP->printLOG("OpenMP: " + int2str(numProcs) + " processors available\n");
}

DcVar::~DcVar() {
}

bool DcVar::Run() {
  bool runSuccess = false;
  if(snpInputType == SNP_SRC_PLINK) {
    runSuccess = RunPlink();
  }
  if(snpInputType == SNP_SRC_FILE) {
    if(chipSeqMode) {
      runSuccess = RunOMRFChipSeq();
    } else {
      runSuccess = RunOMRF();
    }
  }
  return runSuccess;
}

pair<uint, uint> DcVar::MapSnpIndexToPlinkPhenos(uint snpIndex, string varModel) {
  pair<uint, uint> retPair;
  uint nAff = 0;
  uint nUnAff = 0;
  for(int sampleIdx=0; sampleIdx < PP->n; sampleIdx++) {
    Individual* person = PP->sample[sampleIdx];
    // cout << "variantIdx: " << variantIdx << endl;
    // cout << "sampleIdx: " << sampleIdx << endl;
    // cout << "phenotype: " << person->phenotype << endl;
    // cout << "aff: " << person->aff << endl;
    // cout << "locus size: " << PP->locus.size() << endl;
    // cout << "SNP allele1 size: " << person->one.size() << endl;
    // cout << "SNP allele2 size: " << person->two.size() << endl;
    bool i1 = person->one[snpIndex];
    bool i2 = person->two[snpIndex];
    // cout << "i1: "<< i1 << ", i2:  " << i2 << endl;
    // see Caleb's email of 2/24/15 for phenotype assignment based on var model param
    // and bit-wise genotype encoding
    double thisPheno = -9;
    bool thisAff = false;
    if(i1) {
      if(!i2) {
        // 10 het
        thisPheno = 1;
        ++nAff;
        thisAff = true;
      } else {
        // 11
        thisPheno = 1;
        ++nAff;
        thisAff = true;
      }
    } else {
      // 01 // het 
      if(i2) {
        if(varModel == "rec") {
          thisPheno = 1;
        ++nAff;
          thisAff = true;
        }
        else {
          if(varModel == "dom") {
            thisPheno = 0;
            ++nUnAff;
            thisAff = false; 
          }
         // else "hom" missing pheno = -9
        }
      }
      // 00
      else {
        thisPheno = 0; // hom
        ++nUnAff;
        thisAff = false;
      }
    }
    // cout 
    //        << "Variant index: " << variantIdx << "\t" 
    //        << "[ " << i1 << ", " << i2 << " ]" << "\t"
    //        << "Sample: " << sampleIdx << "\t" 
    //        << "Phenotype: " << thisPheno << "\t" 
    //        << "Affected: " << thisAff
    //        << endl;
    person->phenotype = thisPheno;
    person->aff = thisAff;
    person->missing = false;
  }

  retPair.first = nAff;
  retPair.second = nUnAff;
  
  return retPair;
}

bool DcVar::RunPlink() {
  if(chipSeqMode) {
    PP->printLOG("ChIP-seq not supported with PLINK files (yet)\n");
    return false;
  }
  PP->printLOG("Preparing dcVar analysis on PLINK files\n");
  // NOTE: THE SNP2Ind() CALL IS CRITICAL!!! 2/24/15
  PP->SNP2Ind();
  int numSnps = PP->nl_all;
  snpNames.resize(numSnps);
  for(uint snpNameIdx=0; snpNameIdx < numSnps; snpNameIdx++) {
    snpNames[snpNameIdx] = PP->locus[snpNameIdx]->name;
  }
  int numGenes = PP->nlistname.size();
  geneExprNames.resize(numGenes);
  copy(PP->nlistname.begin(), PP->nlistname.end(), geneExprNames.begin());
  numCombs = (numGenes * (numGenes - 1.0)) / 2.0;
  PP->printLOG("Number of genes [ " + dbl2str(numGenes) + " ]\n");
  PP->printLOG("Number of gene interactions [ " + dbl2str(numCombs) + " ]\n");
  // make sure we have variants
  if(numSnps < 1) {
    error("Variants file must specified at least one variant for this analysis!");
  }
  // make sure we have genes
  if(numGenes < 2) {
    error("Gene expression file must specified for this analysis!");
  }
  PP->printLOG(int2str(numSnps) + " variants, and " + int2str(numGenes) + " genes\n");
  // --------------------------------------------------------------------------
  // for all variants
  for(uint snpIdx=0; snpIdx < numSnps; ++snpIdx) {
    string variantName = PP->locus[snpIdx]->name;
    PP->printLOG("\n-----[ " + variantName + " ] (" + int2str(snpIdx + 1) + 
                 " of " + int2str(numSnps) + ")-----\n");
    // ------------------------------------------------------------------------
    // get variant info as case-control phenotype based on variant model
    // cout << "PP->n: " << PP->n << endl;
    // cout << "sample size: " << PP->sample.size() << endl;
    pair<uint, uint> caseControl = 
            MapSnpIndexToPlinkPhenos(snpIdx, par::dcvar_var_model);
    if((caseControl.first < 3) || (caseControl.second < 3)) {
      cerr << "Skipping variant phenotype mapping" << endl << "G1: " 
              << caseControl.first << "\tG2:" << caseControl.second << endl;
      continue;
    }
    // ------------------------------------------------------------------------
    // run dcGAIN for this variant phenotype
    zVals.zeros(numGenes, numGenes);
    pVals.zeros(numGenes, numGenes);
    if(!armaDcgain(zVals, pVals)) {
      continue;
    }
    // DEBUG
    // cout << "results" << endl << results.submat(0,0,4,4) << endl;
    // cout << "interactionPvals" << endl << interactionPvals.submat(0,0,4,4) << endl;
    // armaWriteMatrix(results, "DEBUG.dcgain", PP->nlistname);
    // armaWriteMatrix(interactionPvals, "DEBUG.interactionPvals", PP->nlistname);
    // ------------------------------------------------------------------------
    // adjust p-values
    if(par::do_dcvar_pfilter) {
      if(par::verbose) PP->printLOG("\tp-value filtering requested\n");
      uint numFiltered;
      FilterPvalues(numFiltered);
      if(par::verbose) {
        PP->printLOG("\t[ " + int2str(numFiltered) + " ] values filtered\n");
        PP->printLOG("\t[ " + int2str(zVals.n_nonzero) + " ] values pass filtering\n");
      }
    } else {
      if(par::verbose) PP->printLOG("\tNo p-value filtering requested so skipping filter\n");
    }
    // ------------------------------------------------------------------------
    // write results, if there are any to write
    if(zVals.n_nonzero) {
      PP->printLOG("\t[ " + int2str(zVals.n_nonzero) + 
                   " ] values pass filtering (if used)\n");
      string resultsFilename = 
              par::output_file_name + "." + 
              par::dcvar_pfilter_type + "." +
              variantName + 
              ".pass.tab";
      WriteResults(resultsFilename, variantName);
    } else {
      PP->printLOG("\tWARNING: nothing to write for [ " + variantName + " ]\n");
    }
    // write in case the job fails in this loop; resume with command line flag
    WriteCheckpoint(snpIdx, variantName);
  } // END all variants loop

  return true;
}

bool DcVar::RunOMRF() {
  // ---------------------------------------------------------------------------
  PP->printLOG("DcVar::RunOMRF: Performing dcVar analysis on .gz and .tab files\n");
  uint numSnps = snpNames.size();
  // expression
  uint numGenes = geneExprNames.size();
  // chipseq
  uint numChipSeq = 0;
  if(chipSeqMode) {
    numChipSeq = chipSeqExpression.size();
  }
  // make sure we have variants
  if(numSnps < 1) {
    error("SNP genotypes file must include at least one SNP for analysis!");
  }
  // make sure we have genes
  if(numGenes < MIN_NUM_GENES) {
    error("Gene expression data must include at least [ " + int2str(MIN_NUM_GENES) + " ]\n");
  }
  PP->printLOG("Read [ " + int2str(numSnps) + " ] variants, and [ " + 
               int2str(numGenes) + " ] genes\n");

  if(par::do_dcvar_pfilter) {
    PP->printLOG("Filtering p-values using [ " +  par::dcvar_pfilter_type +  " ] correction\n");
    PP->printLOG("Filtering p-values parameter [ " +  dbl2str(par::dcvar_pfilter_value) +  " ]\n");
  }
  
  // ---------------------------------------------------------------------------
  // for all genotypes/SNPs across all subjects, make genotype into binary 
  // phenotype and run differential correlation on the RNA-Seq gene pairs
  uint initSnpIdx = 0;
  if(par::dcvar_resume_snp) {
    pair <uint, string> snpInfo;
    ReadCheckpoint(snpInfo);
    initSnpIdx = snpInfo.first;
  }
  for(uint snpIdx = initSnpIdx; snpIdx < numSnps; ++snpIdx) {
    string snpName = snpNames[snpIdx];
    if(par::verbose) PP->printLOG("--------------------------------------------------------\n");
    PP->printLOG("SNP [ " + snpName + " ] " + int2str(snpIdx) + " of " + 
                 int2str(numSnps) + "\n");
    // ------------------------------------------------------------------------
    if(par::verbose) PP->printLOG("\tCreating phenotype from SNP genotypes\n");
    vector<uint> snpGenotypes;
    for(uint colIdx=0; colIdx < genotypeSubjects.size(); ++colIdx) {
      snpGenotypes.push_back(static_cast<uint>(genotypeMatrix[snpIdx][colIdx]));
    }
    // get variant genotypes for all subject and map to a genetic model
    if(par::verbose) PP->printLOG("\tGenotypes case-control status\n");    
    vector<uint> mappedPhenos;
    MapPhenosToModel(snpGenotypes, 
                     par::dcvar_var_model, 
                     mappedPhenos);
    if(par::verbose) cout << "\tCases:    " << caseIdxCol.size() << "\t";
    if(par::verbose) cout << "Controls: " << ctrlIdxCol.size() << endl;
    // ------------------------------------------------------------------------
    if(par::verbose) PP->printLOG("\tSplitting into case-control groups\n");
    uint numCases = caseIdxCol.size();
    uint numCtrls = ctrlIdxCol.size();
    if((numCases < MIN_NUM_SUBJ_PER_GROUP) || (numCtrls < MIN_NUM_SUBJ_PER_GROUP)) {
      if(par::verbose) PP->printLOG("\tWARNING: groups sizes must be greater than [ " + 
                   int2str(MIN_NUM_SUBJ_PER_GROUP - 1) + " ], skipping SNP\n");
      continue;
    }
    mat casesMatrix(numCases, numGenes);
    mat ctrlsMatrix(numCtrls, numGenes);
    // split into case-control groups for testing DC
    if(!SplitExpressionCaseControl(casesMatrix, ctrlsMatrix)) {
      error("Could not split on case control status");
    }
    // ------------------------------------------------------------------------
    if(par::verbose) PP->printLOG("\tComputeDifferentialCorrelationZals " 
                 "and first pass p-value filter [ " + 
                 dbl2str(DEFAULT_PVALUE_THRESHOLD) + " ]\n");
    // sparse matrix of significant p-values
    if(!ComputeDifferentialCorrelationZsparse(snpName, 
                                              casesMatrix, 
                                              ctrlsMatrix)) {
      error("ComputeDifferentialCorrelationZvals failed");
    }
    // ------------------------------------------------------------------------
    // adjust p-values
    if(par::do_dcvar_pfilter) {
      if(par::verbose) PP->printLOG("\tp-value filtering requested\n");
      uint numFiltered;
      FilterPvalues(numFiltered);
      if(par::verbose) {
        PP->printLOG("\t[ " + int2str(numFiltered) + " ] values filtered\n");
        PP->printLOG("\t[ " + int2str(zVals.n_nonzero) + " ] values pass filtering\n");
      }
    } else {
      if(par::verbose) PP->printLOG("\tNo p-value filtering requested so skipping filter\n");
    }
    // ------------------------------------------------------------------------
    // write results, if there are any to write
    if(zVals.n_nonzero) {
      string resultsFilename = 
              par::output_file_name + "." + 
              par::dcvar_pfilter_type + "." +
              snpName + 
              ".pass.tab";
      WriteResults(resultsFilename, snpName);
    } else {
      PP->printLOG("\tWARNING: nothing to write for [ " + snpName + " ]\n");
    }

    // write in case the job fails in this loop; resume with command line flag
    WriteCheckpoint(snpIdx, snpName);
  } // end for all SNPs
  
//  zout.close();
  
  return true;
}

bool DcVar::RunOMRFChipSeq() {
  // ---------------------------------------------------------------------------
  PP->printLOG("DcVar::RunOMRF: Performing dcVar analysis on .gz and .tab files\n");
  uint numGenes = geneExprNames.size();
  if(chipSeqMode) {
  } else {
    return false;
  }
  uint numChipseq = chipSeqExpression.size();
  // make sure we have ChIP-Seq
  if(numChipseq < 1) {
    error("ChIP-Seq expression file must include at least one for analysis!");
  }
  // make sure we have genes
  if(numGenes < MIN_NUM_GENES) {
    error("Gene expression data must include at least [ " + int2str(MIN_NUM_GENES) + " ]\n");
  }
  PP->printLOG("Read [ " + int2str(numChipseq) + " ] ChIP-Seq, and [ " + 
               int2str(numGenes) + " ] genes\n");
  // how to filter output
  if(par::do_dcvar_pfilter) {
    PP->printLOG("Filtering p-values using [ " +  par::dcvar_pfilter_type +  " ] correction\n");
    PP->printLOG("Filtering p-values parameter [ " +  dbl2str(par::dcvar_pfilter_value) +  " ]\n");
  }
  
  // ---------------------------------------------------------------------------
  // for all genotypes/SNPs across all subjects, make genotype into binary 
  // phenotype and run differential correlation on the RNA-Seq gene pairs
  uint initChipseqIdx = 0;
  if(par::dcvar_resume_snp) {
    pair <uint, string> snpInfo;
    ReadCheckpoint(snpInfo);
    initChipseqIdx = snpInfo.first;
  } else {
      initChipseqIdx = 1;
  }
  CHIP_SEQ_INFO chipseqInfo = chipSeqExpression[initChipseqIdx];
  string prevChrom = "";
  string curChrom = chipseqInfo.chrom;
  
  for(uint chipseqIdx = initChipseqIdx; chipseqIdx < numChipseq; ++chipseqIdx) {
    chipseqInfo = chipSeqExpression[chipseqIdx];
    string chipseqSnpName = chipseqInfo.rsnum;
    curChrom = chipseqInfo.chrom;
    uint curPos = chipseqInfo.position;
    if(curChrom != prevChrom) {
      PP->printLOG("\tChromosome [ " + curChrom + " ]\n");
      string chromStr(curChrom.begin() + 3, curChrom.end());
      uint chromNum = lexical_cast<uint>(chromStr);
      ReadGenotypesFile(chromNum);
      ReadSnpLocationsFile(chromNum);
      prevChrom = curChrom;
    }
    if(par::verbose) PP->printLOG("--------------------------------------------------------\n");
    PP->printLOG("\tChIP-Seq related SNP [ " + chipseqSnpName + " ] " + 
                 int2str(chipseqIdx) + " of " + int2str(numChipseq) + "\n");
    // find all SNPs within radius of the ChIP-Seq location
    if(par::verbose) 
      PP->printLOG("\tSearching for SNPs within radius [ +/- " +  int2str(radius) + " bp]\n");
    vector<uint> foundSnps;
    FindSnps(curPos, foundSnps);
    if(!foundSnps.size()) {
      if(par::verbose) PP->printLOG("\tWARNING: radius search found no SNPs\n");
      continue;
    }
    if(par::verbose) PP->printLOG("\tFound [ " + int2str(foundSnps.size()) + " ]\n");
    for(uint foundSnpIdx=0; foundSnpIdx < foundSnps.size(); ++foundSnpIdx) {
      string foundSnpName = snpNames[foundSnpIdx];
      // ------------------------------------------------------------------------
      if(par::verbose) PP->printLOG("\tCreating phenotype from SNP genotypes\n");
      vector<uint> snpGenotypes;
      for(uint colIdx=0; colIdx < genotypeSubjects.size(); ++colIdx) {
        snpGenotypes.push_back(static_cast<uint>(genotypeMatrix[foundSnpIdx][colIdx]));
      }
      // get variant genotypes for all subject and map to a genetic model
      if(par::verbose) PP->printLOG("\tGenotypes case-control status\n");    
      vector<uint> mappedPhenos;
      MapPhenosToModel(snpGenotypes, 
                       par::dcvar_var_model, 
                       mappedPhenos);
      if(par::verbose) cout << "\tCases:    " << caseIdxCol.size() << "\t";
      if(par::verbose) cout << "Controls: " << ctrlIdxCol.size() << endl;
      // ------------------------------------------------------------------------
      if(par::verbose) PP->printLOG("\tSplitting into case-control groups\n");
      uint numCases = caseIdxCol.size();
      uint numCtrls = ctrlIdxCol.size();
      if((numCases < MIN_NUM_SUBJ_PER_GROUP) || (numCtrls < MIN_NUM_SUBJ_PER_GROUP)) {
        if(par::verbose) PP->printLOG("\tWARNING: groups sizes must be greater than [ " + 
                     int2str(MIN_NUM_SUBJ_PER_GROUP - 1) + " ], skipping SNP\n");
        continue;
      }
      mat casesMatrix(numCases, numGenes);
      mat ctrlsMatrix(numCtrls, numGenes);
      // split into case-control groups for testing DC
      if(!SplitExpressionCaseControl(casesMatrix, ctrlsMatrix)) {
        error("Could not split on case control status");
      }
      // ------------------------------------------------------------------------
      if(par::verbose) PP->printLOG("\tComputeDifferentialCorrelationZals " 
                   "and first pass p-value filter [ " + 
                   dbl2str(DEFAULT_PVALUE_THRESHOLD) + " ]\n");
      // sparse matrix of significant p-values
      if(!ComputeDifferentialCorrelationZsparse(foundSnpName, 
                                                casesMatrix, 
                                                ctrlsMatrix)) {
        error("ComputeDifferentialCorrelationZvals failed");
      }
      // ------------------------------------------------------------------------
      // adjust p-values
      if(par::do_dcvar_pfilter) {
        if(par::verbose) PP->printLOG("\tp-value filtering requested\n");
        uint numFiltered;
        FilterPvalues(numFiltered);
        if(par::verbose) {
          PP->printLOG("\t[ " + int2str(numFiltered) + " ] values filtered\n");
          PP->printLOG("\t[ " + int2str(zVals.n_nonzero) + " ] values pass filtering\n");
        }
      } else {
        if(par::verbose) PP->printLOG("\tNo p-value filtering requested so skipping filter\n");
      }
      // ------------------------------------------------------------------------
      // write results, if there are any to write
      if(zVals.n_nonzero) {
        string resultsFilename = 
                par::output_file_name + "." + 
                par::dcvar_pfilter_type + "." +
                foundSnpName + 
                ".pass.tab";
        WriteResults(resultsFilename, foundSnpName);
      } else {
        PP->printLOG("\tWARNING: nothing to write for [ " + foundSnpName + " ]\n");
      }
      // write in case the job fails in this loop; resume with command line flag
      WriteCheckpoint(chipseqIdx, foundSnpName);
    } // end for all SNPs found within radius
  } // end for all ChIP-Seq SNPs
  
//  zout.close();
  
  return true;
}

bool DcVar::CheckInputs() {
  return true;
}

bool DcVar::FindSnps(uint pos, std::vector<uint>& inRadius) {
  inRadius.clear();
  bool foundAll = false;
  uint startPos = pos - radius;
  uint endPos = pos + radius;
  if(par::verbose) {
    PP->printLOG("\tChIP-Seq position: " + int2str(pos) + " => ( " + 
                 int2str(startPos) + ", " + int2str(endPos) + " )\n");
  }
  uint searched = 0;
  for(SNP_INFO_LIST_IT it=snpLocations.begin(); 
      !foundAll && it != snpLocations.end(); 
      ++it) {
    ++searched;
    uint thisPos = it->position;
    if((thisPos > startPos) && (thisPos < endPos)) {
      inRadius.push_back(thisPos);
    } else {
      if(thisPos >= endPos) {
        foundAll = true;
      }
    }
  }
  if(par::verbose) {
    PP->printLOG("\tSearched: [ " + int2str(searched) + 
                 " ] => Found: [ " + int2str(inRadius.size()) + " ]\n");
  }
  return foundAll;
}

void DcVar::PrintState() {
  PP->printLOG("-----------------------------------------------------------\n");
  PP->printLOG("p-value adjust method:          " + par::dcvar_pfilter_type + "\n");
  PP->printLOG("p-value cutoff for file output: " + dbl2str(par::dcvar_pfilter_value) + "\n");
  if(par::dcvar_chip_seq_file != "") {
    PP->printLOG("ChIP-seq expression file:       " + par::dcvar_chip_seq_file + "\n");
    PP->printLOG("ChIP-Seq search radius        : " + int2str(par::dcvar_radius) + "\n");
  }
  PP->printLOG("-----------------------------------------------------------\n");
}

bool DcVar::ReadGenotypesFile(uint chrom) {
  string genotypesFilename = "data/genotype." + int2str(chrom) + ".txt.gz";
  checkFileExists(genotypesFilename);
  PP->printLOG("Reading genotypes input from [ " + genotypesFilename + " ]\n");
  ZInput zin(genotypesFilename, compressed(genotypesFilename));
  // read header line
  PP->printLOG("Getting genotype subject names from first line header\n");
  vector<string> tok = zin.tokenizeLine();
  genotypeSubjects.clear();
	for(int i=1; i < tok.size(); i++) {
    genotypeSubjects.push_back(tok[i]);
  }
  genotypeMatrix.clear();
  snpNames.clear();
  uint lineCounter = 1;
  PP->printLOG("Reading compressed genotypes\n");
  while(!zin.endOfFile()) {
    ++lineCounter;
	  vector<string> tok = zin.tokenizeLine();
    if(tok.size() < 2) {
      cerr << "WARNING: line [ " << lineCounter 
              << " ] from [ " << genotypesFilename
              << " ] . . . skipping" << endl;
      continue;
    }
    snpNames.push_back(tok[0]);
    vector<double> lineGenotypes;
	  for(int j=1; j < tok.size(); j++) {
      lineGenotypes.push_back(lexical_cast<double>(tok[j]));
	  }
    genotypeMatrix.push_back(lineGenotypes);
	}
  zin.close();

  PP->printLOG("Read genotypes for [ " + 
  int2str(genotypeSubjects.size()) + " ] subjects and [ " + 
  int2str(snpNames.size()) + " ] SNPs\n");
  
  return true;
}

bool DcVar::ReadSnpLocationsFile(uint chrom) {
  string locationsFilename = "data/SNP_location." + int2str(chrom) + ".txt.gz";
  checkFileExists(locationsFilename);
  PP->printLOG("Reading SNP locations input from [ " + locationsFilename + " ]\n");
  ZInput zin(locationsFilename, compressed(locationsFilename));
  PP->printLOG("Reading and discarding first line header\n");
  zin.tokenizeLine();
  uint lineCounter = 1;
  snpLocations.clear();
  while(!zin.endOfFile()) {
    ++lineCounter;
	  vector<string> parsedFileLine = zin.tokenizeLine();
    if(parsedFileLine.size() != (SNP_DBSNP_ALLELE + 1)) {
      cerr << "WARNING: reading line [ " << lineCounter 
              << " ] from " << par::dcvar_snp_locations_file 
              << " should have 5 columns, found " << parsedFileLine.size()
              << ". Blank line(s)?"
              << endl;
      continue;
    }
    SNP_INFO thisSnpInfo;
    thisSnpInfo.chrom = parsedFileLine[SNP_CHROM];
    thisSnpInfo.position = lexical_cast<uint>(parsedFileLine[SNP_POS]);
    thisSnpInfo.refAllele = parsedFileLine[SNP_DBSNP_ALLELE][0];
    thisSnpInfo.rsnum = parsedFileLine[SNP_ID];
    snpLocations.push_back(thisSnpInfo);
	}
  zin.close();
  PP->printLOG("Read subject SNP location info for [ " + 
               int2str(lineCounter) + " ] SNPs\n");
  
  return true;
}

bool DcVar::ReadGeneExpressionFile() {
  checkFileExists(par::dcvar_gene_expression_file);

  PP->printLOG("Reading gene expression input from [ " + par::dcvar_gene_expression_file + " ]\n");
  ifstream exprFile(par::dcvar_gene_expression_file);
  PP->printLOG("Getting gene expression subject names from first line header\n");
  string header;
  getline(exprFile, header);
  vector<string> headerParts;
  split(headerParts, header, "\t");
  for(uint i=1; i < headerParts.size(); ++i) {
    geneExprSubjects.push_back(headerParts[i]);
  }

  PP->printLOG("Getting gene names from first column, remaining columns gene expression\n");
  string line;
  uint lineCounter = 1;
  while(getline(exprFile, line)) {
    ++lineCounter;
	  vector<string> parsedFileLine;
    split(parsedFileLine, line, "\t");
    if(parsedFileLine.size() < 2) {
      cerr << "Error reading line [ " << lineCounter 
              << " ] from " << par::dcvar_gene_expression_file 
              << " should have more than 2 columns (subjects)"
              << endl;
      continue;
    }
    geneExprNames.push_back(parsedFileLine[0]);
    vector_t thisExprRec;
    for(uint i=1; i < parsedFileLine.size(); ++i) {
      thisExprRec.push_back(lexical_cast<double>(parsedFileLine[i]));
    }
    expressionMatrix.push_back(thisExprRec);
	}
  exprFile.close();
   
  PP->printLOG("Read gene expression for [ " + 
  int2str(geneExprSubjects.size()) + " ] subjects and [ " + 
  int2str(geneExprNames.size()) + " ] genes\n");
  
  double numGenes = (double) geneExprNames.size();
  numCombs = (numGenes * (numGenes - 1.0)) / 2.0;
  PP->printLOG("Number of gene interactions [ " + dbl2str(numCombs) + " ]\n");

  return true;
}

bool DcVar::ReadChipSeqFile() {
  checkFileExists(par::dcvar_chip_seq_file);
  PP->printLOG("Reading ChIP-seq input from [ " + par::dcvar_chip_seq_file + " ]\n");
  ifstream chipSeqFile(par::dcvar_chip_seq_file);
  PP->printLOG("Reading and discarding first line header\n");
  string header;
  getline(chipSeqFile, header);
  uint lineCounter = 1;
  string line;
  while(getline(chipSeqFile, line)) {
    ++lineCounter;
	  vector<string> parsedFileLine;
    split(parsedFileLine, line, "\t");
    if(parsedFileLine.size() != (CHIP_SEQ_SNP + 1)) {
      cerr << "WARNING: reading line [ " << lineCounter 
              << " ] from " << par::dcvar_chip_seq_file 
              << " should have 16 columns, found " << parsedFileLine.size()
              << ". Blank line(s)? Attempting to continuing reading" << endl;
      continue;
    }
    CHIP_SEQ_INFO thisChipSeqInfo;
    thisChipSeqInfo.chrom = parsedFileLine[CHIP_SEQ_CHROM];
    thisChipSeqInfo.position = lexical_cast<uint>(parsedFileLine[CHIP_SEQ_POS]);
    thisChipSeqInfo.totalRegionReads = lexical_cast<uint>(parsedFileLine[CHIP_SEQ_EXPR]);
    // rs28469609:38367404:C:T
    vector<string> rsnumParts;
    split(rsnumParts, parsedFileLine[CHIP_SEQ_SNP], ":");
    thisChipSeqInfo.rsnum = rsnumParts[0];
    chipSeqExpression.push_back(thisChipSeqInfo);
	}
  chipSeqFile.close();
  PP->printLOG("Read ChIP-seq expression for " + int2str(lineCounter) + " SNPs\n");
  sort(chipSeqExpression.begin(), chipSeqExpression.end(), chipseqComparatorAscending);
  
  return true;
}

// build a new phenotype from variant genotypes
bool DcVar::MapPhenosToModel(vector<uint> phenos, string varModel,
                             vector<uint>& mappedPhenos) {
  caseIdxCol.clear();
  ctrlIdxCol.clear();
  for(uint phenoIdx=0; phenoIdx < phenos.size(); ++phenoIdx) {
    uint thisPheno = phenos[phenoIdx];
    int thisMappedPheno = -9;
    string varModel = par::dcvar_var_model;
    if(varModel == "dom") {
      thisMappedPheno = (thisPheno == 2)? 1: 0;
    } else {
      if(varModel == "rec") {
        thisMappedPheno = (thisPheno == 0)? 1: 0;
      } else {
        // hom
        thisMappedPheno = (thisPheno == 1)? 1: -9;
      }
    }
    if(thisMappedPheno) {
      caseIdxCol.push_back(phenoIdx);
    } else {
      ctrlIdxCol.push_back(phenoIdx);
    }
    mappedPhenos.push_back(thisMappedPheno);
  }
  return true;
}

bool DcVar::SplitExpressionCaseControl(mat& caseMatrix, 
                                       mat& ctrlMatrix) {
  uint nGenes = geneExprNames.size();
  uint nCases = caseIdxCol.size();
  uint nCtrls = ctrlIdxCol.size();
  for(uint col=0; col < nCases; ++col) {
    uint thisIndexCol = caseIdxCol[col];
    for(uint row=0; row < nGenes; ++row) {
      caseMatrix(col, row) = expressionMatrix[row][thisIndexCol];
    }
  }
  for(uint col=0; col < nCtrls; ++col) {
    uint thisIndexCol = ctrlIdxCol[col];
    for(uint row=0; row < nGenes; ++row) {
      ctrlMatrix(col, row) = expressionMatrix[row][thisIndexCol];
    }
  }
  
  return true;
}

bool DcVar::ComputeDifferentialCorrelationZvals(string snp, 
                                                mat& cases, 
                                                mat& ctrls) {
  if(par::verbose) PP->printLOG("\tPerforming Z-tests for all RNA-seq interactions\n");
  double n1 = static_cast<double>(caseIdxCol.size());
  double n2 = static_cast<double>(ctrlIdxCol.size());
  uint numGenes = geneExprNames.size();
  double minP = 1.0;
  double maxP = 0.0;
  uint goodPvalCount = 0;
  uint badPvalCount = 0;
  uint infCount = 0;
  double pThreshold = DEFAULT_PVALUE_THRESHOLD;
  if(par::dcvar_pfilter_type == "custom") {
    pThreshold = par::dcvar_pfilter_value;
  }
  if(par::verbose) PP->printLOG("\tFirst pass filter threshold [ " + 
     dbl2str(pThreshold) + " ]\n");
  if(par::verbose) PP->printLOG("\tEntering OpenMP parallel section for [ ");
  if(par::verbose) PP->printLOG(int2str(numCombs) + " ] dcvar combination\n");
  uint i=0, j=0;
#pragma omp parallel for schedule(dynamic, 1) private(i, j)
  for(i=0; i < numGenes; ++i) {
    for(j=i + 1; j < numGenes; ++j) {
      // correlation between this interaction pair (i, j) in cases and controls
      vec caseVarVals1 = cases.col(i);
      vec caseVarVals2 = cases.col(j);
      vec r_ij_1_v = cor(caseVarVals1, caseVarVals2);
      double r_ij_1 = (double) r_ij_1_v[0];
      vec ctrlVarVals1 = ctrls.col(i);
      vec ctrlVarVals2 = ctrls.col(j);
      vec r_ij_2_v = cor(ctrlVarVals1, ctrlVarVals2);
      double r_ij_2 = (double) r_ij_2_v[0];
      // differential correlation Z
      double z_ij_1 = 0.5 * log((abs((1 + r_ij_1) / (1 - r_ij_1))));
      double z_ij_2 = 0.5 * log((abs((1 + r_ij_2) / (1 - r_ij_2))));
      double Z_ij = abs(z_ij_1 - z_ij_2) / sqrt((1.0 / (n1 - 3.0) + 1.0 / (n2 - 3.0)));
      #pragma omp critical 
      {
        // !NOTE! critical section
        double p = DEFAULT_PVALUE;
        if(std::isinf(Z_ij)) {
          // bad Z
          ++infCount;
        } else {
          //matrixElement interactionPvalElement;
          //pair<uint, uint> indexPair = make_pair(i, j);
          p = 2 * normdist(-abs(Z_ij));
          if(p < minP) minP = p;
          if(p > maxP) maxP = p;
          if(p <= pThreshold) {
            ++goodPvalCount;
            string outString = snp +  "\t"  +
                    geneExprNames[i] + "\t"  +
                    geneExprNames[j] + "\t"  +
                    dbl2str(p);
            //zout.writeLine(outString);
          } else {
            ++badPvalCount;
          } // end else good p
        } // end else good Z
      } // end openmp critical section
    } // end for j cols
  } // end for i rows
  if(par::verbose) PP->printLOG("End OpenMP parallel section");
  if(par::verbose) PP->printLOG("\tminp [" + dbl2str(minP) + " ] "
                                "maxp [ " + dbl2str(maxP) + " ]\n");
  
  if(par::verbose) PP->printLOG("\t[ " + int2str(infCount) + " ] infinite Z values, no p-values\n");
  if(par::verbose) PP->printLOG("\t[ " + int2str(badPvalCount) + " ] p-values failed threshold test\n");
  if(par::verbose) PP->printLOG("\t[ " + int2str(goodPvalCount) + " ] p-values passed threshold test\n");
  totalTests = goodPvalCount + badPvalCount + infCount;
  if(par::verbose) PP->printLOG("\t[ " + int2str(totalTests) + " ] total tests\n");
  
  return true;
}

bool DcVar::ComputeDifferentialCorrelationZ(string snp, 
                                            mat& cases, 
                                            mat& ctrls,
                                            double correctedP) {
  PP->printLOG("\tPerforming Z-tests for all rna-seq interactions\n");
  double n1 = static_cast<double>(caseIdxCol.size());
  double n2 = static_cast<double>(ctrlIdxCol.size());
  uint numVars = geneExprNames.size();
  double minP = 1.0;
  double maxP = 0.0;
  uint goodPvalCount = 0;
#pragma omp parallel for schedule(dynamic, 1)
  for(uint i=0; i < numVars; ++i) {
    for(uint j=i + 1; j < numVars; ++j) {
      // correlation between this interaction pair (i, j) in cases and controls
      vec caseVarVals1 = cases.col(i);
      vec caseVarVals2 = cases.col(j);
      vec r_ij_1_v = cor(caseVarVals1, caseVarVals2);
      double r_ij_1 = (double) r_ij_1_v[0];
      vec ctrlVarVals1 = ctrls.col(i);
      vec ctrlVarVals2 = ctrls.col(j);
      vec r_ij_2_v = cor(ctrlVarVals1, ctrlVarVals2);
      double r_ij_2 = (double) r_ij_2_v[0];
      // differential correlation Z
      double z_ij_1 = 0.5 * log((abs((1 + r_ij_1) / (1 - r_ij_1))));
      double z_ij_2 = 0.5 * log((abs((1 + r_ij_2) / (1 - r_ij_2))));
      double Z_ij = abs(z_ij_1 - z_ij_2) / sqrt((1.0 / (n1 - 3.0) + 1.0 / (n2 - 3.0)));
      double p = 2 * normdist(-abs(Z_ij)); 
      if(std::isinf(Z_ij)) {
        cerr << "InfiniteZ" << "\t"
                << snp << "\t"
                << geneExprNames[i] << "\t" 
                << geneExprNames[j] << "\t" 
                << Z_ij << "\t"
                << p 
                << endl;
        continue;
      }
      bool writeResults = false;
      if(par::do_dcvar_pfilter) {
        if(p < correctedP) {
          writeResults = true;
        }
      } else {
        writeResults = true;
      }
      if(writeResults) {
        ++goodPvalCount;
#pragma omp critical 
{
      cout << snp << "\t"
              << geneExprNames[i] << "\t" 
              << geneExprNames[j] << "\t" 
              << Z_ij << "\t"
              << p 
              << endl;
}
      }
    }
  }

  PP->printLOG("\t[ " + int2str(goodPvalCount) + " ] p-values passed threshold test\n");
  
  return true;
}

bool DcVar::ComputeDifferentialCorrelationZsparse(string snp, 
                                                  mat& cases, 
                                                  mat& ctrls) {
  if(par::verbose) PP->printLOG("\tPerforming Z-tests for all RNA-seq interactions\n");
  double n1 = static_cast<double>(cases.n_rows);
  double n2 = static_cast<double>(ctrls.n_rows);
  uint numGenes = geneExprNames.size();
  double minP = 1.0;
  double maxP = 0.0;
  uint goodPvalCount = 0;
  uint badPvalCount = 0;
  uint infCount = 0;
  double pThreshold = DEFAULT_PVALUE_THRESHOLD;
  if(par::dcvar_pfilter_type == "custom") {
    pThreshold = par::dcvar_pfilter_value;
  }
  if(par::verbose) PP->printLOG("\tFirst pass filter threshold [ " + 
     dbl2str(pThreshold) + " ]\n");
  if(par::verbose) PP->printLOG("\tEntering OpenMP parallel section for [ "+ 
                                int2str(numCombs) + " ] dcvar combinations\n");
  zVals.set_size(numGenes, numGenes);
  pVals.ones(numGenes, numGenes);
  uint i, j;
#pragma omp parallel for private(j) collapse(2)
  for(i=0; i < numGenes; ++i) {
    for(j=0; j < numGenes; ++j) {
      // correlation between this interaction pair (i, j) in cases and controls
      vec caseVarVals1 = cases.col(i);
      vec caseVarVals2 = cases.col(j);
      vec r_ij_1_v = cor(caseVarVals1, caseVarVals2);
      double r_ij_1 = (double) r_ij_1_v[0];
      vec ctrlVarVals1 = ctrls.col(i);
      vec ctrlVarVals2 = ctrls.col(j);
      vec r_ij_2_v = cor(ctrlVarVals1, ctrlVarVals2);
      double r_ij_2 = (double) r_ij_2_v[0];
      // differential correlation Z
      double z_ij_1 = 0.5 * log((abs((1 + r_ij_1) / (1 - r_ij_1))));
      double z_ij_2 = 0.5 * log((abs((1 + r_ij_2) / (1 - r_ij_2))));
      double Z_ij = abs(z_ij_1 - z_ij_2) / sqrt((1.0 / (n1 - 3.0) + 1.0 / (n2 - 3.0)));
      #pragma omp critical 
      {
        if(par::verbose) {
          if(i && ((i % 1000) == 0) && (j == (i + 1))) {
            PP->printLOG(int2str(i) + " of " + int2str(numGenes) + "\n");
          }
        }
        // !NOTE! critical section for writing to an already opened file 'zout'!
        // and to keep track of counts, min/max p-values in public scope
        double p = DEFAULT_PVALUE;
        if(std::isinf(Z_ij)) {
          // bad Z
          ++infCount;
        } else {
          //matrixElement interactionPvalElement;
          //pair<uint, uint> indexPair = make_pair(i, j);
          p = 2 * normdist(-abs(Z_ij));
          if(p < minP) minP = p;
          if(p > maxP) maxP = p;
          if(p <= pThreshold) {
            ++goodPvalCount;
//            string outString = snp +  "\t"  +
//                    geneExprNames[i] + "\t"  +
//                    geneExprNames[j] + "\t"  +
//                    dbl2str(p);
//            zout.writeLine(outString);
            zVals(i, j) = Z_ij;
            pVals(i, j) = p;
          } else {
            ++badPvalCount;
          } // end else good p
        } // end else good Z
      } // end openmp critical section
    } // end for j cols
  } // end for i rows
  if(par::verbose) PP->printLOG("End OpenMP parallel section");
  if(par::verbose) PP->printLOG("\tminp [" + dbl2str(minP) + " ] "
                                "maxp [ " + dbl2str(maxP) + " ]\n");
  
  if(par::verbose) PP->printLOG("\t[ " + int2str(infCount) + " ] infinite Z values, no p-values\n");
  if(par::verbose) PP->printLOG("\t[ " + int2str(badPvalCount) + " ] p-values failed threshold test\n");
  if(par::verbose) PP->printLOG("\t[ " + int2str(goodPvalCount) + " ] p-values passed threshold test\n");
  totalTests = goodPvalCount + badPvalCount + infCount;
  if(par::verbose) PP->printLOG("\t[ " + int2str(totalTests) + " ] total tests\n");
  
  return true;
}

bool DcVar::FlattenPvals(vector_t& retPvals) {
  if(par::verbose) PP->printLOG("Flattening p-values list into a vector\n");
  retPvals.clear();
  for(uint i=0; i < pVals.n_rows; ++i) {
    for(uint j=i + 1; j < pVals.n_cols; ++j) {
      retPvals.push_back(pVals(i, j));
    }
  }
  
  return true;
}

bool DcVar::FilterPvalues(uint& numFiltered) {
  double filterThreshold = par::dcvar_pfilter_value;
  if(par::dcvar_pfilter_type == "fdr") {
    filterThreshold = CalculateFdrBHThreshold();
  } else {
    if(par::dcvar_pfilter_type == "bon") {
      filterThreshold = par::dcvar_pfilter_value / (numCombs * snpNames.size());
    } else {
      if(par::dcvar_pfilter_type == "custom") {
        // do nothing, but include for later mods context
        filterThreshold = par::dcvar_pfilter_value;
      } else {
        error("Unknown p-value filter type. Expects \"bon\" or \"fdr\" or \"custom\"."   
              "Got [ " + par::dcvar_pfilter_type + " ]");
      }
    }
  }

  if(par::verbose) {
    PP->printLOG("\tCustom pruning with method [ " + par::dcvar_pfilter_type + " ]\n");
    PP->printLOG("\tCustom pruning with correctedP [ " + 
                 dbl2str(filterThreshold) + " ]\n");
  }
  uint numPruned = 0;
  uint numGenes = geneExprNames.size();
  for(uint i=0; i < numGenes; ++i) {
    for(uint j=i + 1; j < numGenes; ++j) {
      if(pVals(i, j) > filterThreshold) {
        zVals(i, j) = 0;
        zVals(j, i) = 0;
        ++numPruned;
      }
    }
  }

  if(par::verbose) PP->printLOG("\tPruned [ " + int2str(numPruned) + " ] values from interaction terms\n");

  numFiltered = numPruned;
  
  return true;
}

double DcVar::CalculateFdrBHThreshold() {
  vector_t interactionPvals;
  FlattenPvals(interactionPvals);
  if(par::verbose) 
    PP->printLOG("\tCalculating FDR using Benjamini-Hochberg (BH) for pruning\n");
  uint m = interactionPvals.size() * snpNames.size();
  // sort interaction by p-value
  sort(interactionPvals.begin(), interactionPvals.end());
  // use rough FDR (RFDR) to estimate alpha based on input FDR 
  // from encore, Nick Davis code
  double alpha = 2 * m * DEFAULT_FDR / (m + 1);
  int R = -1;
  // BH method
  for(int i = 0; i < m; i++) {
    double l = (i + 1) * alpha / (double) m;
    // test whether current p-value < current l
    if(interactionPvals[i] < l) {
      R = i;
    } else {
      break;
    }
  }
  // BH threshold condition not met with any p-values, so exit
  if(R == -1) {
    if(par::verbose) PP->printLOG("\tWARNING: No p-value meets BH threshold criteria, so no pruning\n");
    return 0.0;
  }
  // BH rejection threshold
  double T = interactionPvals[R];
  if(par::verbose)  {
    PP->printLOG("\tBH rejection threshold: T = [ " + dbl2str(T) + " ], R = " + int2str(R) + "\n");
    PP->printLOG("\tPruning interactions with p-values > T [ (" + dbl2str(T) + " ])\n");
  }
  
  return T;
}

bool DcVar::WriteCheckpoint(uint snpIndex, string snpName) {
  ofstream checkpointFile(CHECKPOINT_FILENAME);
  checkpointFile << snpIndex << endl << snpName << endl;
  checkpointFile.close();
  
  return true;
}

bool DcVar::ReadCheckpoint(std::pair<uint, string>& lastSnp) {
  ifstream checkpointFile(CHECKPOINT_FILENAME);
  uint snpIndex;
  string snpName;
  string fileLine;
  getline(checkpointFile, fileLine);
  snpIndex = lexical_cast<uint>(fileLine);
  getline(checkpointFile, snpName);
  checkpointFile.close();
  lastSnp.first = snpIndex;
  lastSnp.second = snpName;
  
  return true;
}

bool DcVar::WriteResults(string filename, string curSnp) {
  if(par::verbose) 
    PP->printLOG("\tWriting interactions that passed p-value filter to [ "  + 
                 filename + " ]\n");
  // avoid writing empty matrix/zero-byte file
  // make no assumption that called has checked; display warning and return
  sp_mat::const_iterator zmatItCurrent = zVals.begin();
  sp_mat::const_iterator zmatItEnd = zVals.end();
  // no elements to write for this SNP?
  if(zmatItCurrent == zmatItEnd) {
    PP->printLOG("\tWARNING: DcVar::WriteResults method attempt to write "
                 "empty z-values sparse matrix\n");
    return false;
  }
  if((zVals.n_rows != pVals.n_rows) ||
     (zVals.n_cols != pVals.n_cols)) {
    PP->printLOG("\tWARNING: DcVar::WriteResults method attempt to write "
                 "z-values matrix dimensions not equal to the p-values matrix\n");
    PP->printLOG("\tZ: " + int2str(zVals.n_rows) + " x " + int2str(zVals.n_cols) + "\n");
    PP->printLOG("\tp: " + int2str(pVals.n_rows) + " x " + int2str(pVals.n_cols) + "\n");
    return false;
  }
  ofstream resultsFile(filename);
  resultsFile << "SNP\tGene1\tGene2\tZ\tP" << endl;
  for(; zmatItCurrent != zmatItEnd; ++zmatItCurrent) {
    double zvalue = *zmatItCurrent;
    uint row = zmatItCurrent.row();
    uint col = zmatItCurrent.col();
    double pvalue = pVals(row, col);
    resultsFile 
        << curSnp << "\t" 
        << geneExprNames[row] << "\t" 
        << geneExprNames[col] << "\t" 
        << std::scientific 
        << zvalue << "\t"
        << pvalue << "\t"
        << endl;
  }
  resultsFile.close();
  
  return true;
}

// from EpistasisEQtl - bcw - 1/3/18

bool DcVar::SetRadius(int newRadius) {
  if(newRadius < 1) {
    cerr << "Error setting cis radius to: " << newRadius << endl;
    return false;
  }
  // newRadius is in kilobases
  radius = newRadius;
  return true;
}

bool DcVar::SetLocalCis(bool localCisFlag) {
  localCis = localCisFlag;
  return true;
}

bool DcVar::SetTFRadius(int newRadius) {
  if(newRadius < 0) {
    cerr << "Error setting TF radius to: " << newRadius << endl;
    return false;
  }
  // newRadius is in kilobases, but need to store a bases
  tfRadius = newRadius * 1000;
  return true;
}

bool DcVar::SetTF(bool tfFlag) {
  tfMode = tfFlag;
  return true;
}

bool DcVar::GetSnpsForTranscript(string transcript, 
  vector<int>& snpIndices) {

  // get transcript info
  int chromosome = coordinates[transcript][COORD_CHROM];
  int bpStart = coordinates[transcript][COORD_BP_START];
  int bpEnd = coordinates[transcript][COORD_BP_END];
  int lowerThreshold = bpStart - radius;
  int upperThreshold = bpEnd + radius;
  
 // cout 
 //   << "GetSnpsForTranscript: chrom: " << chromosome  << ", radius: " 
 //   << radius  << ", " 
 //   << "(" << lowerThreshold << "), "
 //   << bpStart  << ", " 
 //   << bpEnd << ", "
 //   << "(" << upperThreshold << ")"
 //   << endl;
  
  // find SNPs matching criteria
  for(int j=0; j < PP->locus.size(); ++j) {
    Locus* thisSnp = PP->locus[j];
    if(thisSnp->chr == chromosome) {
      if(localCis) {
        // on the same chromosome and within radius of transcript
        if((thisSnp->bp >= lowerThreshold) && 
           (thisSnp->bp <= upperThreshold)) {
          snpIndices.push_back(j);
        }
      }
      else {
        // simply on the same chromosome
        snpIndices.push_back(j);
      }
    }
  }

  return true;
}

bool DcVar::GetSnpsForTFs(vector<int>& snpIndices, vector<string>& tfs) {

  int allSnps = PP->locus.size();
  PP->printLOG("Searching transcription factors in " + int2str(allSnps) + " SNPs\n");
  // for all SNPs
  for(int thisSnpIndex=0; thisSnpIndex < allSnps; ++thisSnpIndex) {
    if(thisSnpIndex && (thisSnpIndex % 100000 == 0)) {
      cout << thisSnpIndex << "/" << allSnps << endl;
    }
    Locus* thisSnp = PP->locus[thisSnpIndex];
    int chr = thisSnp->chr;
    int bp = thisSnp->bp;
    // is this bp in range of any transcription factors?
    string tf;
    if(IsSnpInTFs(chr, bp, tf)) {
      snpIndices.push_back(thisSnpIndex);
      tfs.push_back(tf);
    }
  }

  return true;
}

bool DcVar::IsSnpInTFs(int chr, int bp, string& tf) {
  // linear search through the transcription factor lookup table for SNP at
  // chr/bp, returning the transcription factor if true, else return false
  bool found = false;
  TranscriptFactorTableCIt lutIt = transcriptFactorLUT.begin();
  while((lutIt != transcriptFactorLUT.end()) && (!found)) {
    string thisTF = (*lutIt).first;
    int thisTfChr = transcriptFactorLUT[thisTF][COORD_CHROM];
    int thisTfBpBeg = transcriptFactorLUT[thisTF][COORD_BP_START];
    int thisTfBpEnd = transcriptFactorLUT[thisTF][COORD_BP_END];
    int rangeStart = thisTfBpBeg - tfRadius;
    int rangeEnd = thisTfBpEnd + tfRadius;
    if((chr == thisTfChr) && ((bp >= rangeStart) &&  (bp <= rangeEnd))) {
      tf = thisTF;
      return true;
    }
    ++lutIt;
  }
  return false;
}

bool DcVar::GetTFInfo(string tf, vector<int>& tfInfo) {
  TranscriptFactorTableCIt lutInfo = transcriptFactorLUT.find(tf);
  if(lutInfo == transcriptFactorLUT.end()) {
    cerr << "GetTFRange failed to find: " << tf << endl;
    return false;
  }
  for(int i=0; i < lutInfo->second.size(); ++i) {
    tfInfo.push_back(lutInfo->second[i]);
  }
  return true;
}
