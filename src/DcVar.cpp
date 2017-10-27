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

bool pvalComparatorAscending(const matrixElement &l, const matrixElement &r) {
  return l.first < r.first;
}

DcVar::DcVar(SNP_INPUT_TYPE snpInputTypeParam, bool hasChipSeq, bool debugFlag) {
  PP->printLOG("dcVar initializing\n");
  snpInputType = snpInputTypeParam;
  chipSeq = hasChipSeq;
  debugMode = debugFlag;
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
    PP->printLOG("Using separate tab-delimited files for input data sets\n");
  } else {
    // assume PLINK data structures accessible through PP pointer
    PP->printLOG("Using PLINK files for input data sets\n");
  }
  SetDebugMode(debugFlag);
  if(!CheckInputs()) {
    error("Checking data sets compatability failed. Exiting.");
  }
}

DcVar::~DcVar() {
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

// build a new phenotype from variant genotypes
bool DcVar::MapPhenosToModel(vector<uint> phenos, string varModel) {
  caseIdxCol.clear();
  ctrlIdxCol.clear();
  for(uint i=0; i < phenos.size(); ++i) {
    uint thisPheno = phenos[i];
    uint thisMappedPheno = 0;
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
      caseIdxCol.push_back(i);
    } else {
      ctrlIdxCol.push_back(i);
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

bool DcVar::RunOMRF(bool debugFlag) {
  PP->printLOG("Performing dcVar analysis on .gz and .tab files\n");
  uint numVariants = snpNames.size();
  // expression
  uint numGenes = geneExprNames.size();
  // chipseq
  uint numChipSeq = 0;
  if(chipSeq) {
    numChipSeq = chipSeqExpression.size();
  }
    // make sure we have variants
  if(numVariants < 1) {
    error("Variants file must specified at least one variant for this analysis!");
  }
  // make sure we have genes
  if(numGenes < 2) {
    error("Gene expression file must specified for this analysis!");
  }
  PP->printLOG("[ " + int2str(numVariants) + " variants, and [ " + 
               int2str(numGenes) + " ] genes\n");
  
  stringstream ssResultsFilename;
  stringstream ssErrorsFilename;
  ssResultsFilename << "dcvar." << par::output_file_name << ".pass.tab";
  ssErrorsFilename << "dcvar." << par::output_file_name << ".err.tab";
  string resultsFilename = ssResultsFilename.str();
  string errorsFilename = ssErrorsFilename.str();
  resultsFile.open(resultsFilename);
  errorsFile.open(errorsFilename);
  // for all variants
  for(uint snpIdx = 0; snpIdx != numVariants; ++snpIdx) {
    string snpName = snpNames[snpIdx];
    PP->printLOG("--------------------------------------------------------\n");
    PP->printLOG("SNP: " + snpName + "\n");
    // ------------------------------------------------------------------------
    PP->printLOG("\tcreating phenotype from SNP genotypes\n");
    vector<uint> snpGenotypes;
    for(uint colIdx=0; colIdx < genotypeSubjects.size(); ++colIdx) {
      snpGenotypes.push_back(static_cast<uint>(genotypeMatrix[snpIdx][colIdx]));
    }
    // get variant genotypes for all subject and map to a genetic model
    cout << "Genotypes" << endl;
    MapPhenosToModel(snpGenotypes, "dom");
    cout << "\tCases:    " << caseIdxCol.size() << "\t";
    cout << "Controls: " << ctrlIdxCol.size() << endl;
    // ------------------------------------------------------------------------
    PP->printLOG("\tsplitting into case-control groups\n");
    uint nCases = caseIdxCol.size();
    uint nCtrls = ctrlIdxCol.size();
    if((nCases == 0) || (nCtrls == 0)) {
      PP->printLOG("\tWARNING: groups size must be greater than 0\n");
      continue;
    }
    mat X(nCases, numGenes);
    mat Y(nCtrls, numGenes);
    // split into case-control groups for testing DC
    if(!SplitExpressionCaseControl(X, Y)) {
      error("Could not split on case control status");
    }
    // ------------------------------------------------------------------------
    PP->printLOG("\tComputeDifferentialCorrelationZ\n");
    if(!ComputeDifferentialCorrelationZ(snpName, X, Y)) {
      error("ComputeDifferentialCorrelationZ failed");
    }
  } // end for all variants
  resultsFile.close();
  errorsFile.close();
  
  return true;
}

bool DcVar::ComputeDifferentialCorrelationZ(string variant, 
                                            mat& cases, 
                                            mat& ctrls) {
  PP->printLOG("Performing Z-tests for interactions\n");
  double n1 = static_cast<double>(caseIdxCol.size());
  double n2 = static_cast<double>(ctrlIdxCol.size());
  uint numVars = geneExprNames.size();
  for(int i=0; i < numVars; ++i) {
    for(int j=i + 1; j < numVars; ++j) {
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
        errorsFile << "InfiniteZ" << "\t"
                << variant << "\t"
                << geneExprNames[i] << "\t" 
                << geneExprNames[j] << "\t" 
                << Z_ij << "\t"
                << p 
                << endl;
      }
      bool writeResults = false;
      if(par::do_dcvar_pfilter) {
        if(p < par::dcvar_pfilter_value) {
          writeResults = true;
        }
      } else {
        writeResults = true;
      }
      if(writeResults) {
        resultsFile 
                << variant << "\t"
                << geneExprNames[i] << "\t" 
                << geneExprNames[j] << "\t" 
                << Z_ij << "\t"
                << p 
                << endl;
      }
    }
  }

  return true;
}

bool DcVar::ComputeDifferentialCorrelationZnaive(mat& X, mat& Y) {
  // cout << "X: " << X.n_rows << " x " << X.n_cols << endl;
  // cout << "Y: " << Y.n_rows << " x " << Y.n_cols << endl;
  // cout << "X" << endl << X.submat(0,0,4,4) << endl;
  // cout << "Y" << endl << Y.submat(0,0,4,4) << endl;
  // compute covariances/correlations
  mat covMatrixX;
  mat corMatrixX;
  if(!armaComputeCovariance(X, covMatrixX, corMatrixX)) {
    error("Could not compute coexpression matrix for cases");
  }
  mat covMatrixY;
  mat corMatrixY;
  if(!armaComputeCovariance(Y, covMatrixY, corMatrixY)) {
    error("Could not compute coexpression matrix for controls");
  }

  // DEBUG
  // cout << corMatrixX.n_rows << " x " << corMatrixX.n_cols << endl;
  // cout << "cor(X)" << endl << corMatrixX.submat(0,0,4,4) << endl;
  // cout << "cor(Y)" << endl << corMatrixY.submat(0,0,4,4) << endl;

  // algorithm from R script z_test.R
  PP->printLOG("Performing Z-tests for interactions\n");
  double n1 = caseIdxCol.size();
  double n2 = ctrlIdxCol.size();
  int goodFdrCount = 0;
  double minP = 1.0;
  double maxP = 0.0;
  uint numVars = geneExprNames.size();
  for(int i=0; i < numVars; ++i) {
    for(int j=0; j < numVars; ++j) {
      if(j <= i) {
        continue;
      }
      double r_ij_1 = corMatrixX(i, j);
      double r_ij_2 = corMatrixY(i, j);
      double z_ij_1 = 0.5 * log((abs((1 + r_ij_1) / (1 - r_ij_1))));
      double z_ij_2 = 0.5 * log((abs((1 + r_ij_2) / (1 - r_ij_2))));
      double Z_ij = abs(z_ij_1 - z_ij_2) / sqrt((1.0 / (n1 - 3.0) + 1.0 / (n2 - 3.0)));
      double p = 2 * normdist(-abs(Z_ij)); 
      if(std::isinf(Z_ij)) {
        cerr << "Infinity found at (" << i << ", " << j << ")" << endl;
      }
      // if(i == 0 && j < 10) {
      //   printf("%d, %d => %10.2f %g\n", i, j, Z_ij, p);
      // }
//      resultsMatrix[i][j] = Z_ij;
//      resultsMatrix[j][i] = Z_ij;
//      if(par::do_regain_pvalue_threshold) {
//        if(p > par::regainPvalueThreshold) {
//          resultsMatrix[i][j] = 0;
//          resultsMatrix[j][i] = 0;
//        }
//      }
//      resultsMatrixPvals[i][j] = p;
//      resultsMatrixPvals[j][i] = p;
    }
  }
  
  return true;
}

bool DcVar::CheckInputs() {
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
      cerr << "WARNING: line [ " << lineCounter 
              << " ] from [ " << par::dcvar_genotypes_file << " ]" << endl;
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
      cerr << "WARNING: reading line [ " << lineCounter 
              << " ] from " << par::dcvar_snp_locations_file 
              << " should have 5 columns, found " << tok.size()
              << ". Blank line(s)?"
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
      continue;
    }
    geneExprNames.push_back(tok[0]);
    vector<double> thisExprRec;
    for(uint i=1; i < tok.size(); ++i) {
      thisExprRec.push_back(lexical_cast<double>(tok[i]));
    }
    expressionMatrix.push_back(thisExprRec);
	}
  exprFile.close();
   
  PP->printLOG("Read gene expression for [ " + 
  int2str(geneExprSubjects.size()) + " ] subjects and [ " + 
  int2str(geneExprNames.size()) + " ] genes\n");
  
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
      cerr << "WARNING: reading line [ " << lineCounter 
              << " ] from " << par::dcvar_chip_seq_file 
              << " should have 16 columns, found " << tok.size()
              << ". Blank line(s)?" << endl;
      continue;
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
