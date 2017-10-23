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

// Boost lexical cast
#include <boost/lexical_cast.hpp>

// PLINK
#include "options.h"
#include "plink.h"
#include "helper.h"
#include "zed.h"

// inbix
#include "DcVar.h"
#include "Insilico.h"
#include "StringUtils.h"
#include "ArmadilloFuncs.h"

using namespace std;
using namespace insilico;
using namespace boost;
using namespace arma;

bool pvalComparatorAscending(const matrixElement &l, const matrixElement &r) {
  return l.first < r.first;
}

DcVar::DcVar(SNP_INPUT_TYPE snpInputTypeParam, bool hasChipSeq, bool useDebug) {
  snpInputType = snpInputTypeParam;
  chipSeq = hasChipSeq;
  debugMode = useDebug;
  if(snpInputTypeParam == SNP_SRC_FILE) {
    // OMRF data files gzipped and tab-delimited
    if(!ReadGenotypesFile()) {
      error("Reading genotypes failed. Exiting.");
    }
    if(!ReadSnpLocationsFile()) {
      error("Reading SNP location information file failed. Exiting.");
    }
    if(!ReadGeneExpressionFile()) {
      error("Reading gene expression file failed. Exiting.");
    }
    chipSeqMode = hasChipSeq;
    if(chipSeqMode) {
      if(!ReadChipSeqFile()) {
        error("Reading ChIP-seq file failed. Exiting.");
      }
    }
  } else {
    // assume PLINK data structures accessible through PP pointer
  }
  SetDebugMode(useDebug);
  if(!CheckInputs()) {
    error("Checking data sets compatability failed. Exiting.");
  }
}

bool DcVar::Run(bool debugFlag) {
  if(snpInputType == SNP_SRC_PLINK) {
    RunPlink(debugFlag);
  }
  if(snpInputType == SNP_SRC_FILE) {
    RunOMRF(debugFlag);
  }
  return true;
}

bool DcVar::RunPlink(bool debugFlag) {
  if(chipSeqMode) {
    PP->printLOG("ChIP-seq not supported with PLINK files (yet)\n");
    return false;
  }
  PP->printLOG("Preparing dcVar analysis on PLINK files\n");
  // NOTE: THE SNP2Ind() CALL IS CRITICAL!!! 2/24/15
  PP->SNP2Ind();
  int numVariants = PP->nl_all;
  int numGenes = PP->nlistname.size();
  // make sure we have variants
  if(numVariants < 1) {
    error("Variants file must specified at least one variant for this analysis!");
  }
  // make sure we have genes
  if(numGenes < 2) {
    error("Gene expression file must specified for this analysis!");
  }
  PP->printLOG(int2str(numVariants) + " variants, and " + int2str(numGenes) + " genes\n");

  // for all variants
  unsigned int variantIdx;
  for(variantIdx=0; variantIdx < numVariants; ++variantIdx) {
    string variantName = PP->locus[variantIdx]->name;
    // PP->printLOG("Variant: " + variantName + "\n");
    // get variant info as case-control phenotype based on variant model
    // cout << "PP->n: " << PP->n << endl;
    // cout << "sample size: " << PP->sample.size() << endl;
    for(int sampleIdx=0; sampleIdx < PP->n; sampleIdx++) {
      Individual* person = PP->sample[sampleIdx];
      // cout << "variantIdx: " << variantIdx << endl;
      // cout << "sampleIdx: " << sampleIdx << endl;
      // cout << "phenotype: " << person->phenotype << endl;
      // cout << "aff: " << person->aff << endl;
      // cout << "locus size: " << PP->locus.size() << endl;
      // cout << "SNP allele1 size: " << person->one.size() << endl;
      // cout << "SNP allele2 size: " << person->two.size() << endl;
      bool i1 = person->one[variantIdx];
      bool i2 = person->two[variantIdx];
      // cout << "i1: "<< i1 << ", i2:  " << i2 << endl;
      // see Caleb's email of 2/24/15 for phenotype assignment based on var model param
      // and bit-wise genotype encoding
      double thisPheno = -9;
      bool thisAff = false;
      if(i1) {
        if(!i2) {
          // 10 het
          thisPheno = 1;
          thisAff = true;
        } else {
          // 11
          thisPheno = 1;
          thisAff = true;
        }
      } else {
        // 01 // het 
        if(i2) {
          if(par::dcvar_var_model == "rec") {
            thisPheno = 1;
            thisAff = true;
          }
          else {
            if(par::dcvar_var_model == "dom") {
              thisPheno = 0;
              thisAff = false;
            }
            // else "hom" missing pheno = -9
          }
        }
        // 00
        else {
          thisPheno = 0; // hom
          thisAff = false;
        }
      }
      // cout 
      // 	<< "Variant index: " << variantIdx << "\t" 
      // 	<< "[ " << i1 << ", " << i2 << " ]" << "\t"
      // 	<< "Sample: " << sampleIdx << "\t" 
      // 	<< "Phenotype: " << thisPheno << "\t" 
      // 	<< "Affected: " << thisAff
      // 	<< endl;
      person->phenotype = thisPheno;
      person->aff = thisAff;
    }

    // run dcGAIN for this variant phenotype
    mat results(numGenes, numGenes);
    mat pvals(numGenes, numGenes);
    armaDcgain(results, pvals);
    // DEBUG
    // cout << "results" << endl << results.submat(0,0,4,4) << endl;
    // cout << "pvals" << endl << pvals.submat(0,0,4,4) << endl;
    // armaWriteMatrix(results, "DEBUG.dcgain", PP->nlistname);
    // armaWriteMatrix(pvals, "DEBUG.pvals", PP->nlistname);

    // save p-values that pass BH rejection threshold
    // setup output file
    string dcvarFilename = variantName + ".dcVarTest.txt";
    ofstream dcvarFile;
    PP->printLOG("Writing results to [ " + dcvarFilename + " ]\n");
    dcvarFile.open(dcvarFilename.c_str());
    if(dcvarFile.fail()) {
      error("Cannot open dcVar test results file for writing.");
    }
    dcvarFile.precision(6);
    dcvarFile.fixed;

    if(par::do_dcvar_pfilter) {
      double nVars = (double) numGenes;
      double nCombs = (nVars * (nVars - 1.0)) / 2.0;
      double minP = 1.0;
      double maxP = 1.0;
      int goodPvalCount = 0;
      if(par::dcvar_pfilter_type == "fdr") {
        // ------------------------------------------------------------------------
        PP->printLOG("Filtering using Benjamini-Hochberg FDR threshold\n");
        // get all p-values
        vector<matrixElement> test_pvals;
        for(int i=0; i < pvals.n_rows; ++i) {
          for(int j=0; j < pvals.n_cols; ++j) {
            if(j <= i) { continue; }
            test_pvals.push_back(make_pair(pvals(i, j), make_pair(i, j)));
          }
        }
        // sort them
        sort(test_pvals.begin(), test_pvals.end(), pvalComparatorAscending);
        // use rough FDR (RFDR) to estimate alpha based on input FDR
        int num_pvals = test_pvals.size();
        double m = (double) num_pvals * (double) numVariants;
        double alpha = 2 * m * par::dcvar_pfilter_value  / (m + 1);
        int threshold_index = -1;
        // BH method
        for(int i = 0; i < num_pvals; i++) {
          double l = (i + 1) * alpha / (double) num_pvals;
          // test whether current p-value < current l
          if(test_pvals[i].first < l) {
            threshold_index = i;
          } else {
            break;
          }
        }
        // BH threshold condition not met with any p-values, so break out of this iteration
        if(threshold_index == -1) {
          PP->printLOG("No p-value meets BH threshold criteria, so nothing saved\n");
        } else {
          // BH rejection threshold
          double T = test_pvals[threshold_index].first;
          PP->printLOG("BH rejection threshold T = " + dbl2str(T) + ", R = " +
            int2str(threshold_index) + "\n");
          goodPvalCount = threshold_index + 1;
          minP = test_pvals[0].first;
          maxP = test_pvals[test_pvals.size()-1].first;
          for(int i=0; i < threshold_index; i++) {
            double p = test_pvals[i].first;
            pair<int, int> idx = test_pvals[i].second;
            string gene1 = PP->nlistname[idx.first];
            string gene2 = PP->nlistname[idx.second];
            dcvarFile << gene1 << "\t" << gene2 << "\t" << p << endl;
          }
        }
      } else {
        PP->printLOG("Filtering using Bonferroni threshold\n");
        // insure doubles used in all intermediate calculations
        double correctedP = par::dcvar_pfilter_value / (nCombs * numVariants);
        // cout 
        // 	<< nVars << "\t" 
        // 	<< nCombs <<  "\t"
        // 	<< numVariants << "\t"
        // 	<< correctedP
        // 	<< endl;
        // printf("FDR Corrected p-value: %g\n", correctedP);
        minP = pvals(0, 0);
        maxP = pvals(0, 0);
        for(int i=0; i < pvals.n_rows; ++i) {
          for(int j=0; j < pvals.n_cols; ++j) {
            if(j <= i) { continue; }
            string gene1 = PP->nlistname[i];
            string gene2 = PP->nlistname[j];
            double p = pvals(i, j);
            // printf("p-value [%g] < [%g] ?\n", p, correctedP);
            // cout << gene1 << ", " << gene2 << " => p=" << p << " corrected=" 
            //  << correctedP << " Passed FDR test!" << endl;
            if(p < minP) minP = p;
            if(p > maxP) maxP = p;
            if(p <  correctedP) {
              ++goodPvalCount;
              // cout << gene1 << ", " << gene2 << " => p=" << p << " corrected=" 
              //  << correctedP << " Passed FDR test!" << endl;
              //printf("p-value [%g] < [%g] PASSED!\n", p, correctedP);
              dcvarFile << gene1 << "\t" << gene2 << "\t" << p << endl;
            }
          } // end pvals cols
        } // end pvals rows
      }
      PP->printLOG("Found [" + int2str(goodPvalCount) + "] tested p-values, min/max: " + 
        dbl2str(minP) + " / " + dbl2str(maxP) + "\n");
    } else {
      // no p-value filtering
      PP->printLOG("Saving ALL p-values\n");
      for(int i=0; i < pvals.n_rows; ++i) {
        for(int j=0; j < pvals.n_cols; ++j) {
          if(j <= i) { continue; }
          string gene1 = PP->nlistname[i];
          string gene2 = PP->nlistname[j];
          double p = pvals(i, j);
          dcvarFile << gene1 << "\t" << gene2 << "\t" << p << endl;
        }
      }
    }

    dcvarFile.close();

  } // END all variants loop

  return true;
}

bool DcVar::RunOMRF(bool debugFlag) {
  PP->printLOG("Performing dcVar analysis on OMRF supplied .gz and .tab files\n");
  return true;
}

DcVar::~DcVar() {
}

bool DcVar::CheckInputs() {
  // do all the input data set dimensions make sense?
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
  PP->printLOG("SNP locations file:             " + par::dcvar_snp_locations_file + "\n");
  PP->printLOG("CHiP-seq expression file:       " + par::dcvar_chip_seq_file + "\n");
  PP->printLOG("p-value adjust method:          " + par::dcvar_pfilter_type + "\n");
  PP->printLOG("p-value cutoff for file output: " + dbl2str(par::dcvar_pfilter_value) + "\n");
  PP->printLOG("-----------------------------------------------------------\n");
}

bool DcVar::ReadGenotypesFile() {
  checkFileExists(par::dcvar_genotypes_file);
  PP->printLOG("Reading genotypes input from [ " + par::dcvar_genotypes_file + " ]\n");
  ZInput zin(par::dcvar_genotypes_file, compressed(par::dcvar_genotypes_file));
  // read header line
  PP->printLOG("Getting genotype subject names from first line header\n");
  vector<string> tok = zin.tokenizeLine();
	for(int i=1; i < tok.size(); i++) {
    genotypeSubjects.push_back(tok[i]);
  }
  uint lineCounter = 1;
  while(!zin.endOfFile()) {
    ++lineCounter;
	  vector<string> tok = zin.tokenizeLine();
    if(tok.size() < 2) {
      cerr << "Error reading line [ " << lineCounter 
              << " ] from " << par::dcvar_genotypes_file << endl;
      continue;
    }
    vector<uint> lineGenotypes(tok.size() - 1);
	  for(int j=1; j < tok.size(); j++) {
      lineGenotypes[j - 1] = lexical_cast<uint>(tok[j]);
	  }
    genotypes[tok[0]] = lineGenotypes;
	}
  zin.close();
  PP->printLOG("Read subject genotypes for " + int2str(lineCounter) + " SNPs\n");
  
  return true;
}

bool DcVar::ReadSnpLocationsFile() {
  checkFileExists(par::dcvar_snp_locations_file);
  PP->printLOG("Reading SNP locations input from [ " + par::dcvar_snp_locations_file + " ]\n");
  ZInput zin(par::dcvar_snp_locations_file, compressed(par::dcvar_snp_locations_file));
  PP->printLOG("Reading and discarding first line header\n");
  zin.tokenizeLine();
  uint lineCounter = 0;
  while(!zin.endOfFile()) {
    ++lineCounter;
	  vector<string> tok = zin.tokenizeLine();
    if(tok.size() != 5) {
      cerr << "Error reading line [ " << lineCounter 
              << " ] from " << par::dcvar_snp_locations_file 
              << " should have 5 columns, found " << tok.size()
              << endl;
      continue;
    }
    SNP_INFO thisSnpInfo;
    thisSnpInfo.chrom = tok[1];
    thisSnpInfo.location = lexical_cast<uint>(tok[2]);
    thisSnpInfo.refAllele = tok[4][0];
    snpLocations[tok[0]] = thisSnpInfo;
	}
  zin.close();
  PP->printLOG("Read subject SNP location info for " + int2str(lineCounter) + " SNPs\n");

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
      newRow.push_back(lexical_cast<double>(tok[i]));
    }
    exprMatrix.push_back(newRow);
	}
  exprFile.close();
  PP->printLOG("Read gene expression for " + int2str(lineCounter) + " genes\n");
  
  return true;
}

bool DcVar::ReadChipSeqFile() {
  checkFileExists(par::dcvar_chip_seq_file);
  PP->printLOG("Reading ChIP-seq input from [ " + par::dcvar_chip_seq_file + " ]\n");
  ifstream chipSeqFile(par::dcvar_chip_seq_file);
  PP->printLOG("Reading and discarding first line header\n");
  string header;
  getline(chipSeqFile, header);
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
    thisChipSeqInfo.position = lexical_cast<uint>(tok[CHIP_SEQ_POS]);
    thisChipSeqInfo.totalRegionReads = lexical_cast<double>(tok[CHIP_SEQ_EXPR]);
    // rs28469609:38367404:C:T
    vector<string> rsnumParts;
    split(rsnumParts, tok[CHIP_SEQ_SNP], ":");
    chipSeqExpression[rsnumParts[0]] = thisChipSeqInfo;
	}
  chipSeqFile.close();
  PP->printLOG("Read ChIP-seq expression for " + int2str(lineCounter) + " SNPs\n");
  
  return true;
}
