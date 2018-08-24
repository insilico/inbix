/* 
 * File:   EpistasisEQtl.cpp
 * Author: bwhite
 * 
 * Created on October 3, 2013, 11:48 AM
 */

#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <set>
#include <fstream>

#include <armadillo>

#include "plink.h"
#include "model.h"
#include "linear.h"
#include "stats.h"
#include "helper.h"

#include "EpistasisEQtl.h"
#include "Insilico.h"

using namespace std;

EpistasisEQtl::EpistasisEQtl() {
  radius = -1;
  localCis = false;
  tfMode = false;
  tfTableLoaded = false;
  tfRadius = 0;
  goodModels = 0;
  badModels = 0;
  exprFilename = par::iqtl_expression_file;
  cordFilename = par::iqtl_coord_file;
  fullInteraction = par::iqtl_interaction_full;
}

EpistasisEQtl::~EpistasisEQtl() {
}

bool EpistasisEQtl::SetDebugMode(bool debugFlag) {
  par::algorithm_verbose = debugFlag;
  return true;
}

bool EpistasisEQtl::CheckInputs() {
  // basic assumptions check
  // we have SNPs?
  uint numSnps = PP->nl_all;
  if(!numSnps) {
    error("no SNPs found \n");
  }
  // we have transcript expression levels?
  uint numTranscripts = PP->nlistname.size();
  if(!numTranscripts) {
    error("no transcript values found \n");
  }
  // we have a transcript lookup table that matches the expression data?
  uint numTranscriptInfo = coordinates.size();
  if(numTranscriptInfo != numTranscripts) {
    error("number of coordinate file entries does not match "
      "the number of transcript values found \n");
  }
  vector<string>::const_iterator cit = PP->nlistname.begin();
  for(; cit != PP->nlistname.end(); ++cit) {
    // do all the transcripts in expression values have coordinates info?
    if(coordinates.find(*cit) == coordinates.end()) {
      error("Transcript " + *cit + " not found in coordinate file\n");
    }
  }

  return true;
}

void EpistasisEQtl::PrintState() {
  PP->printLOG("-----------------------------------------------------------\n");
  PP->printLOG(Timestamp() + "EpistasisEqtl object state:\n");
  string verboseFlag = par::verbose? "on": "off";
  PP->printLOG(Timestamp() + "verbose mode: " + verboseFlag + "\n");
  string debugFlag = par::algorithm_verbose? "on": "off";
  PP->printLOG(Timestamp() + "debug mode: " + debugFlag + "\n");
  PP->printLOG(Timestamp() + "expression file: " + exprFilename + "\n");
  PP->printLOG(Timestamp() + "coordinates file: " + cordFilename + "\n");
  if(localCis) {
    PP->printLOG(Timestamp() + "local cis mode with radius: " + 
      int2str(uint(radius / 1000)) + " kilobases\n");
  }
  if(tfMode) {
    PP->printLOG(Timestamp() + "TF mode with radius: " + 
      int2str(uint(tfRadius / 1000)) + " kilobases\n");
  }
  if(fullInteraction) {
    PP->printLOG(Timestamp() + "FULL epistatic interaction mode\n");
  } else {
    if(tfMode) {
      PP->printLOG(Timestamp() + "TF/cis-trans interaction mode\n");  
    } else {
      PP->printLOG(Timestamp() + "cis-cis/cis-trans interaction mode\n");  
    }
  }
  if(tfMode) {
    PP->printLOG(Timestamp() + "transcription factor mode\n");
    PP->printLOG(Timestamp() + "search radius: " + int2str(par::iqtl_tf_radius) + "\n");
    PP->printLOG(Timestamp() + "coordinate file: " + par::iqtl_tf_coord_file + "\n");
  }
  PP->printLOG(Timestamp() + "p-value cutoff for file output: " + dbl2str(par::iqtl_pvalue) + "\n");
  PP->printLOG("-----------------------------------------------------------\n");
  
  return;
}

bool EpistasisEQtl::Run() {
  PP->printLOG(Timestamp() + "Running iQTL\n");
  if(!CheckInputs()) {
    error("iQTL::CheckInputs() failed\n");
  }
  PP->printLOG(Timestamp() + "iQTL::CheckInputs() OK\n");

  // --------------------------------------------------------------------------
  PP->printLOG(Timestamp() + "iQTL linear regression loop for all RNA-Seq transcripts\n");
  PrintState();
  
  // --------------------------------------------------------------------------
  // keep track of loop parameters and number of tests done and write to file
  string testnumbersFilename = par::output_file_name + ".testnumbers.txt";
  PP->printLOG(Timestamp() + "Writing test results to [ " + testnumbersFilename + " ]\n");
  std::ofstream TESTNUMBERS;
  TESTNUMBERS.open(testnumbersFilename, ios::out);

  string loopInfoFilename = par::output_file_name + ".loopinfo.txt";
  PP->printLOG(Timestamp() + "Writing loop information to [ " + loopInfoFilename + " ]\n");
  std::ofstream LOOPINFO;
  LOOPINFO.open(loopInfoFilename, ios::out);

  // --------------------------------------------------------------------------
  // determine SNPs in outer loop; considers transcription factor mode
  nOuterLoop = PP->nl_all;
  vector<string> thisTFSnpNames;
  if(tfMode) {
    GetSnpsForTFs(thisTFSnpIndices, thisTFSnpNames);
    nOuterLoop = thisTFSnpIndices.size();
    if(par::algorithm_verbose) {
      cout << "TF snp indices: ";
      copy(thisTFSnpIndices.begin(), thisTFSnpIndices.end(), 
              ostream_iterator<uint>(cout, "\t"));
      cout << endl;
      cout << "TFs: ";
      copy(thisTFSnpNames.begin(), 
        thisTFSnpNames.end(), 
        ostream_iterator<string>(cout, "\t"));
      cout << endl;
    }
    PP->printLOG(Timestamp() + "nOuterLoop TFs: " + int2str(nOuterLoop) + "\n");
  } else {
    PP->printLOG(Timestamp() + "nOuterLoop SNPs: " + int2str(nOuterLoop) + "\n");
  }

  // --------------------------------------------------------------------------
  // for each transcript build main effect and epistasis regression models
  string thisTranscript;
  for(uint transcriptIndex = 0; 
      transcriptIndex < PP->nlistname.size(); 
      ++transcriptIndex) {
    
    thisTranscript = PP->nlistname[transcriptIndex];
    PP->printLOG("---------------------------------------------------------\n");
    PP->printLOG(Timestamp() + "Transcript: " + thisTranscript + "\n");

    // get SNP indices for the transcript
    nInnerLoop = PP->nl_all;
    // GetSnpsForTranscript takes locl-cis into account
    if(!GetSnpsForTranscript(thisTranscript, thisTranscriptSnpIndices)) {
      error("Could not get SNPs for transcript");
    } else {
      nInnerLoop = thisTranscriptSnpIndices.size();
    }
    
    // get transcript expression vector as phenotype
    PP->setQtlPhenoFromNumericIndex(transcriptIndex);
    
    if(!(nInnerLoop * nOuterLoop)) {
      PP->printLOG("WARNING: no interactions for: " + thisTranscript + "\n");
      continue; // to next transcript
    }
    
    if(par::algorithm_verbose) {
      cout << "DEBUG outerLoopSnps set size: " << outerLoopSnps.size() << endl;
      cout << "transcript cis snp indices: ";
      copy(thisTranscriptSnpIndices.begin(), thisTranscriptSnpIndices.end(), 
              ostream_iterator<uint>(cout, "\t"));
      cout << endl;
      cout << "Loop set snp indices: ";
      copy(thisTFSnpIndices.begin(), thisTFSnpIndices.end(), 
              ostream_iterator<uint>(cout, "\t"));
      cout << endl;
    }
    
    PP->printLOG(Timestamp() + "Writing transcript loop info to: " + loopInfoFilename + "\n");
    PP->printLOG(Timestamp() + 
                 thisTranscript + "\t" + 
                 int2str(nOuterLoop) + "\t" + 
                 int2str(nInnerLoop) + "\n");
    LOOPINFO 
      << thisTranscript << "\t"
      << nInnerLoop << "\t"
      << nOuterLoop << "\t"
      << endl;

    // EQTL -------------------------------------------------------------------
    PP->printLOG(Timestamp() + "Running main effects regression models\n");
    RunEqtl(thisTranscript);
    
    // IQTL -----------------------------------------------------------------
    PP->printLOG(Timestamp() + "Running interaction effects regression models\n");
    // allocate results matrices
    resultsMatrixBetas.resize(nOuterLoop, nInnerLoop);
    resultsMatrixBetas.zeros();
    resultsMatrixPvals.resize(nOuterLoop, nInnerLoop);
    resultsMatrixPvals.ones();
    if(fullInteraction) {
      RunIqtlFull();
    } else {
      RunIqtlCisTrans();
    }

    // ------------------------------------------------------------------------
    PP->printLOG(Timestamp() + "write regression results\n");
    string iqtlFilename = par::output_file_name + "." + 
      thisTranscript + ".iqtl.txt";
    if(!WriteResults(iqtlFilename, thisTranscript, thisTFSnpNames)) {
      error("Failed to write iQTL results file: " + iqtlFilename);
    }
    
    PP->printLOG(Timestamp() + "Writing statistical test numbers to: " + testnumbersFilename + "\n");
    PP->printLOG(Timestamp() + "Total models:\t" + int2str(nOuterLoop * nInnerLoop) + "\n");
    TESTNUMBERS << thisTranscript << "\t" << (nOuterLoop * nInnerLoop) << endl;
  } // END for each transcript loop
  
  TESTNUMBERS.close();
  LOOPINFO.close();
  
  PP->printLOG(Timestamp() + "iQTL analysis finished\n");

  return true;
}

bool EpistasisEQtl::RunIqtlCisTrans() {
  PP->printLOG(Timestamp() + "iQTL Cis-Trans\n");
  if(localCis) {
    PP->printLOG(Timestamp() + "iQTL local cis SNPs: inner loop (" +  int2str(nInnerLoop) + ")\n");
  } else {
    PP->printLOG(Timestamp() + "iQTL SNPs: inner loop (" +  int2str(nInnerLoop) + ")\n");
  }
  if(tfMode) {
    PP->printLOG(Timestamp() + "iQTL TFs: outer loop (" + int2str(nOuterLoop) + ")\n");
  } else {
    PP->printLOG(Timestamp() + "iQTL SNPs: outer loop(" + int2str(nOuterLoop) + ")\n");
  }
  uint badModelsTotal = 0;
  uint goodModelsTotal = 0;
  uint threadGoodModels = 0;
  uint threadBadModels = 0;
  #pragma omp parallel for
  for(uint ii=0; ii < nOuterLoop; ++ii) {
    for(uint jj=0; jj < nInnerLoop; ++jj) {
      // cout << "Looping indices: " << "\t" << ii << "\t" << jj << endl;
      uint snpAIndex = thisTFSnpIndices[ii];
      string snpAName = PP->locus[snpAIndex]->name;
      uint snpBIndex = thisTranscriptSnpIndices[jj];
      string snpBName = PP->locus[snpBIndex]->name;
      Model* interactionModel = new LinearModel(PP);
      interactionModel->setMissing();
      interactionModel->addAdditiveSNP(snpAIndex);
      interactionModel->label.push_back(snpAName);
      interactionModel->addAdditiveSNP(snpBIndex);
      interactionModel->label.push_back(snpBName);
      if(par::covar_file) {
        for(uint kk = 0; kk < par::clist_number; kk++) {
          interactionModel->addCovariate(kk);
          interactionModel->label.push_back(PP->clistname[kk]);
        }
      }
      interactionModel->addInteraction(1, 2);
      interactionModel->label.push_back("EPI");
      interactionModel->buildDesignMatrix();
      uint modelFitParamIdx = 3;
      // add # covars to test param to get interaction param
      if(par::covar_file) {
        modelFitParamIdx += par::clist_number;
      }
      interactionModel->testParameter = modelFitParamIdx; // interaction
      interactionModel->fitLM();
      bool badModel = false;
      if(!interactionModel->isValid()) {
        if(par::verbose) {
          PP->printLOG(Timestamp() + "WARNING: outer index " + int2str(ii) + 
            " inner index " + int2str(jj) + " linear model fitLM(): invalid\n");
        }
        badModel = true;
      }
      if(!interactionModel->fitConverged()) {
        if(par::verbose) {
          PP->printLOG(Timestamp() + "WARNING: linear model fitLM(): failed to converge: SNP A: " + 
            snpAName + ", SNP B: " + snpBName + "\n");
        }
        badModel = true;
      }
      #pragma omp critical
      {
        if(!badModel) {
          ++threadGoodModels;
          ++goodModelsTotal;
          vector_t betaInteractionCoefs = interactionModel->getCoefs();
          double interactionValue = 
            betaInteractionCoefs[betaInteractionCoefs.size() - 1];
          vector_t betaInteractionCoefPVals = interactionModel->getPVals();
          double interactionPval =
            betaInteractionCoefPVals[betaInteractionCoefPVals.size() - 1];
          // double interactionPval = ((LinearModel*) interactionModel)->getPValue();
          resultsMatrixBetas(ii, jj) = interactionValue;
          resultsMatrixPvals(ii, jj) = interactionPval;
        } else {
          ++threadBadModels;
          ++badModelsTotal;
          resultsMatrixBetas(ii, jj) = 0.0;
          resultsMatrixPvals(ii, jj) = 1.0;
        }
      }
      delete interactionModel;
    }
  }
  goodModels = goodModelsTotal;
  badModels = badModelsTotal;
  PP->printLOG(Timestamp() + "iQTL Cis Trans total good models: " + int2str(goodModels) +
    ", bad models: " + int2str(badModels) + "\n");
  
  return true;
}

bool EpistasisEQtl::RunIqtlFull() {
  PP->printLOG(Timestamp() + "iQTL linear regression loop: Full SNP x SNP\n");
  uint numSnps = PP->nl_all;
#pragma omp parallel
  for(uint ii=0; ii < numSnps; ++ii) {
    if(ii && (ii % 1000 == 0)) {
      cout << ii << "/" << numSnps << endl;
    }
    for(uint jj=ii+1; jj < numSnps; ++jj) {
      uint snpAIndex = ii;
      string snpAName = PP->locus[snpAIndex]->name;
      uint snpBIndex = jj;
      string snpBName = PP->locus[snpBIndex]->name;
      Model* interactionModel = new LinearModel(PP);
      interactionModel->setMissing();
      interactionModel->addAdditiveSNP(snpAIndex);
      interactionModel->label.push_back(snpAName);
      interactionModel->addAdditiveSNP(snpBIndex);
      interactionModel->label.push_back(snpBName);
      if(par::covar_file) {
        for(uint kk = 0; kk < par::clist_number; kk++) {
          interactionModel->addCovariate(kk);
          interactionModel->label.push_back(PP->clistname[kk]);
        }
      }
      interactionModel->addInteraction(1, 2);
      interactionModel->label.push_back("EPI");
      interactionModel->buildDesignMatrix();
      uint modelFitParamIdx = 3;
      // add # covars to test param to get interaction param
      if(par::covar_file) {
        modelFitParamIdx += par::clist_number;
      }
      interactionModel->testParameter = modelFitParamIdx; // interaction
      interactionModel->fitLM();
#pragma omp critical
{          
      vector_t betaInteractionCoefs = interactionModel->getCoefs();
      double interactionValue = 
        betaInteractionCoefs[betaInteractionCoefs.size() - 1];
      vector_t betaInteractionCoefPVals = interactionModel->getPVals();
      double interactionPval =
        betaInteractionCoefPVals[betaInteractionCoefPVals.size() - 1];
      resultsMatrixBetas(ii, jj) = interactionValue;
      resultsMatrixPvals(ii, jj) = interactionPval;
}
      delete interactionModel;
    }
  }
    
  return true;
}

bool EpistasisEQtl::RunEqtl(string transcript) {
  string eqtlFilename = par::output_file_name + "." + 
    transcript + ".eqtl.txt";
  PP->printLOG(Timestamp() + "RunEqtl writing eQTL results to [ " + eqtlFilename + " ]\n");
  std::ofstream EQTL;
  EQTL.open(eqtlFilename, ios::out);
  for(uint i=0; i < thisTFSnpIndices.size(); ++i) {
    uint thisSnpIndex = thisTFSnpIndices[i];
    string thisSnpName = PP->locus[thisSnpIndex]->name;

    Model* mainEffectModel = new LinearModel(PP);
    mainEffectModel->setMissing();
    mainEffectModel->addAdditiveSNP(thisSnpIndex);
    mainEffectModel->label.push_back(thisSnpName);
    // add covariates if specified
    if(par::covar_file) {
      for(uint k=0; k < par::clist_number; k++) {
        // add covariate to the model
        mainEffectModel->addCovariate(k);
        mainEffectModel->label.push_back(PP->clistname[k]);
      }
    }
    // Build design matrix
    mainEffectModel->buildDesignMatrix();
    // Fit linear model
    uint modelFitParamIdx = 1; 
    // single variable main effect
    mainEffectModel->testParameter = modelFitParamIdx; 
    mainEffectModel->fitLM();
    // obtain estimates and statistics
    vector_t betaMainEffectCoefs = mainEffectModel->getCoefs();
    double mainEffectValue = betaMainEffectCoefs[modelFitParamIdx];
    // p-values don't include intercept term
    vector_t betaMainEffectCoefPvals = mainEffectModel->getPVals();
    double mainEffectPValue = betaMainEffectCoefPvals[modelFitParamIdx-1];
    EQTL 
      << thisSnpName << "\t"
      << transcript << "\t"
      << mainEffectValue << "\t"
      << mainEffectPValue << endl;
    delete mainEffectModel;
  }
  EQTL.close();

  return true;
}

bool EpistasisEQtl::ReadTranscriptCoordinates(string coordinatesFilename) {
  // open the numeric attributes file if possible
  checkFileExists(coordinatesFilename);
  ifstream coordinatesFile(coordinatesFilename, ios::in);
  if(coordinatesFile.fail()) {
    return false;
  }

  coordinates.clear();
  uint rows = 0;
  while(!coordinatesFile.eof()) {
    char nline[par::MAX_LINE_LENGTH];
    coordinatesFile.getline(nline, par::MAX_LINE_LENGTH, '\n');

    // convert to string
    string sline = nline;
    if(sline == "") continue;

    // read line from text file uinto a vector of tokens
    string buf;
    stringstream ss(sline);
    vector<string> tokens;
    while(ss >> buf) {
      tokens.push_back(buf);
    }
    ++rows;
    
    if(tokens.size() != 4) {
      cerr << "Error reading transcript info on line " << rows << endl;
      return false;
    }

    string gene = tokens[tokens.size()-1];

    // handle special cases of chromosome that are not integers
    uint chrom = -1;
    bool chromAssigned = false;
    if(tokens[0] == "X") { chrom = 23; chromAssigned = true; }
    if(tokens[0] == "Y") { chrom = 24; chromAssigned = true; }
    if(tokens[0] == "XY") { chrom = 25; chromAssigned = true; }
    if(tokens[0] == "MT") { chrom = 26; chromAssigned = true; }
    if(!chromAssigned) {
      uint t = 0;
      if(!from_string<uint>(t, tokens[0], std::dec)) {
        cerr << "Error parsing chromosome to integer on line " << rows << endl;
        cerr << "token: " << tokens[0] << endl;
        return false;
      }
      chrom = t;
      chromAssigned = true;
    }
    coordinates[gene].push_back(chrom);

    for(uint i=1; i < tokens.size()-1; ++i) {
      uint t = 0;
      if(!from_string<uint>(t, tokens[i], std::dec)) {
        cerr << "Error parsing transcript info to integer on line " << rows << endl;
        cerr << "token: " << tokens[i] << endl;
        return false;
      }
      coordinates[gene].push_back(t);
    }
  }
  coordinatesFile.close();

  PP->printLOG(Timestamp() + "Read " + int2str(rows) + 
    " transcript coordinates info from ["  + coordinatesFilename + "]\n");
 
  cordFilename = coordinatesFilename;
  
  return true;
}

bool EpistasisEQtl::ReadTranscriptFactorCoordinates(string coordinatesFilename) {
  // open the numeric attributes file if possible
  checkFileExists(coordinatesFilename);
  ifstream coordinatesFile(coordinatesFilename, ios::in);
  if(coordinatesFile.fail()) {
    return false;
  }

  transcriptFactorLUT.clear();
  uint rows = 0;
  while(!coordinatesFile.eof()) {
    char nline[par::MAX_LINE_LENGTH];
    coordinatesFile.getline(nline, par::MAX_LINE_LENGTH, '\n');

    // convert to string
    string sline = nline;
    if(sline == "") continue;

    // read lines from a text file uinto a vector of string tokens
    string buf;
    stringstream ss(sline);
    vector<string> tokens;
    while(ss >> buf) {
      tokens.push_back(buf);
    }
    ++rows;
    
    if(tokens.size() != 4) {
      cerr << "Error reading transcript info on line " << rows << endl;
      return false;
    }

    // gene name in the last column
    string gene = tokens[tokens.size()-1];

    // handle special cases of chromosome that are not integers
    uint chrom = -1;
    bool chromAssigned = false;
    if(tokens[0] == "X") { chrom = 23; chromAssigned = true; }
    if(tokens[0] == "Y") { chrom = 24; chromAssigned = true; }
    if(tokens[0] == "XY") { chrom = 25; chromAssigned = true; }
    if(tokens[0] == "MT") { chrom = 26; chromAssigned = true; }
    if(!chromAssigned) {
      uint t = 0;
      if(!from_string<uint>(t, tokens[0], std::dec)) {
        cerr << "Error parsing chromosome to integer on line " << rows << endl;
        cerr << "token: " << tokens[0] << endl;
        return false;
      }
      chrom = t;
      chromAssigned = true;
    }
    transcriptFactorLUT[gene].push_back(chrom);

    // start and end bp
    for(uint i=1; i < tokens.size()-1; ++i) {
      uint t = 0;
      if(!from_string<uint>(t, tokens[i], std::dec)) {
        cerr << "Error parsing transcript info to integer on line " << rows << endl;
        cerr << "token: " << tokens[i] << endl;
        return false;
      }
      transcriptFactorLUT[gene].push_back(t);
    }
  }
  coordinatesFile.close();

  PP->printLOG(Timestamp() + "Read " + int2str(rows) + 
    " transcription factor coordinates info from ["  + coordinatesFilename + "]\n");
  
  return true;
}

bool EpistasisEQtl::SetRadius(uint newRadius) {
  if(newRadius < 1) {
    cerr << "Error setting cis radius to: " << newRadius << endl;
    return false;
  }
  // newRadius is in kilobases, but need to store a bases
  radius = newRadius * 1000;
  return true;
}

bool EpistasisEQtl::SetLocalCis(bool localCisFlag) {
  localCis = localCisFlag;
  return true;
}


bool EpistasisEQtl::SetTFRadius(uint newRadius) {
  // newRadius is in kilobases, but need to store a bases
  tfRadius = newRadius * 1000;
  return true;
}

bool EpistasisEQtl::SetTF(bool tfFlag) {
  tfMode = tfFlag;
  if(tfMode && !tfTableLoaded) {
    LoadDefaultTranscriptionFactorLUT();
    tfTableLoaded = true;
  }
  return true;
}

bool EpistasisEQtl::GetSnpsForTranscript(string transcript, 
  vector<uint>& snpIndices) {
  PP->printLOG(Timestamp() + "Searching for SNPs in transcript [" + transcript + "]\n");

  // get transcript info
  uint chromosome = coordinates[transcript][COORD_CHROM];
  uint bpStart = coordinates[transcript][COORD_BP_START];
  uint bpEnd = coordinates[transcript][COORD_BP_END];
  uint lowerThreshold = bpStart - radius;
  uint upperThreshold = bpEnd + radius;
  
  if(par::algorithm_verbose) {
    cout 
      << "GetSnpsForTranscript: chrom: " << chromosome  << ", radius: " 
      << radius  << ", " 
      << "(" << lowerThreshold << "), "
      << bpStart  << ", " 
      << bpEnd << ", "
      << "(" << upperThreshold << ")"
      << endl;
  }
  
  // find SNPs matching criteria
  for(uint j=0; j < PP->locus.size(); ++j) {
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

bool EpistasisEQtl::GetSnpsForTFs(vector<uint>& snpIndices, vector<string>& tfs) {

  uint allSnps = PP->locus.size();
  PP->printLOG(Timestamp() + "Searching transcription factors in " + int2str(allSnps) + " SNPs\n");
  // for all SNPs
  for(uint thisSnpIndex=0; thisSnpIndex < allSnps; ++thisSnpIndex) {
    if(thisSnpIndex && (thisSnpIndex % 100000 == 0)) {
      cout << thisSnpIndex << "/" << allSnps << endl;
    }
    Locus* thisSnp = PP->locus[thisSnpIndex];
    uint chr = thisSnp->chr;
    uint bp = thisSnp->bp;
    // is this bp in range of any transcription factors?
    string tf;
    if(IsSnpInTFs(chr, bp, tf)) {
      snpIndices.push_back(thisSnpIndex);
      tfs.push_back(tf);
    }
  }

  return true;
}

bool EpistasisEQtl::IsSnpInTFs(uint chr, uint bp, string& tf) {
  // linear search through the transcription factor lookup table for SNP at
  // chr/bp, returning the transcription factor if true, else return false
  bool found = false;
  TranscriptFactorTableCIt lutIt = transcriptFactorLUT.begin();
  while((lutIt != transcriptFactorLUT.end()) && (!found)) {
    string thisTF = (*lutIt).first;
    uint thisTfChr = transcriptFactorLUT[thisTF][COORD_CHROM];
    uint thisTfBpBeg = transcriptFactorLUT[thisTF][COORD_BP_START];
    uint thisTfBpEnd = transcriptFactorLUT[thisTF][COORD_BP_END];
    uint rangeStart = thisTfBpBeg - tfRadius;
    uint rangeEnd = thisTfBpEnd + tfRadius;
    if((chr == thisTfChr) && ((bp >= rangeStart) &&  (bp <= rangeEnd))) {
      tf = thisTF;
      return true;
    }
    ++lutIt;
  }
  return false;
}

bool EpistasisEQtl::GetTFInfo(string tf, vector<uint>& tfInfo) {
  TranscriptFactorTableCIt lutInfo = transcriptFactorLUT.find(tf);
  if(lutInfo == transcriptFactorLUT.end()) {
    cerr << "GetTFRange failed to find: " << tf << endl;
    return false;
  }
  for(uint i=0; i < lutInfo->second.size(); ++i) {
    tfInfo.push_back(lutInfo->second[i]);
  }
  return true;
}

bool EpistasisEQtl::WriteResults(string saveFilename, 
                                 string saveTranscript,
                                 vector<string> saveTFSnpNames) {
    PP->printLOG(Timestamp() + "writing iQTL results to [ " + saveFilename + " ]\n");
    ofstream IQTL_OUT;
    IQTL_OUT.open(saveFilename, ios::out);
    if(tfMode) {
      IQTL_OUT << "SnpA\tSnpB\tTranscript\tTF\tCoef\tP" << endl;
    } else {
      IQTL_OUT << "SnpA\tSnpB\tTranscript\tCoef\tP" << endl;
    }
    for(uint kk=0; kk < nOuterLoop; ++kk) {
      for(uint ll=0; ll < nInnerLoop; ++ll) {
        double thisInteractionPval = resultsMatrixPvals(kk, ll);
        if((thisInteractionPval > 0) && 
           (thisInteractionPval < par::iqtl_pvalue)) {
          string snpTF = "NA";
          uint snpAIndex = -1;
          if(tfMode) {
            snpAIndex = thisTFSnpIndices[kk];
            snpTF = saveTFSnpNames[kk];
          } else {
            snpAIndex = thisTranscriptSnpIndices[kk];
          }
          string snpAName = PP->locus[snpAIndex]->name;
          uint snpBIndex = thisTranscriptSnpIndices[ll];
          string snpBName = PP->locus[snpBIndex]->name;
          if(tfMode) {
            IQTL_OUT
              << snpAName << "\t" 
              << snpBName << "\t"
              << saveTranscript << "\t"
              << snpTF << "\t"
              << resultsMatrixBetas(kk, ll) << "\t"
              << thisInteractionPval << endl;
          } else {
            IQTL_OUT
              << snpAName << "\t" 
              << snpBName << "\t"
              << saveTranscript << "\t"
              << resultsMatrixBetas(kk, ll) << "\t"
              << thisInteractionPval << endl;
          }
        }
      } // nInnerLoop
    } // nOuterLoop
    IQTL_OUT.close();
    
    return true;
}

bool EpistasisEQtl::LoadDefaultTranscriptionFactorLUT() {
  transcriptFactorLUT["ADNP"] = {20, 49505454, 49547527};
  transcriptFactorLUT["AFF1"] = {4, 87856153, 88062206};
  transcriptFactorLUT["AFF2"] = {23, 147582138, 148082193};
  transcriptFactorLUT["AFF3"] = {2, 100163715, 100759037};
  transcriptFactorLUT["AFF4"] = {5, 132211070, 132299354};
  transcriptFactorLUT["AHR"] = {7, 17338275, 17385775};
  transcriptFactorLUT["AHRR"] = {5, 304290, 438405};
  transcriptFactorLUT["AIRE"] = {21, 45705720, 45718102};
  transcriptFactorLUT["ALX1"] = {12, 85674035, 85695561};
  transcriptFactorLUT["ALX3"] = {1, 110602996, 110613322};
  transcriptFactorLUT["ALX4"] = {11, 44282277, 44331716};
  transcriptFactorLUT["AR"] = {23, 66763873, 66950461};
  transcriptFactorLUT["ARGFX"] = {3, 121286777, 121309469};
  transcriptFactorLUT["ARID1A"] = {1, 27022521, 27108601};
  transcriptFactorLUT["ARID1B"] = {6, 157099063, 157531913};
  transcriptFactorLUT["ARID2"] = {12, 46123619, 46301819};
  transcriptFactorLUT["ARID3A"] = {19, 926036, 972803};
  transcriptFactorLUT["ARID3B"] = {15, 74833547, 74890472};
  transcriptFactorLUT["ARID3C"] = {9, 34621454, 34628011};
  transcriptFactorLUT["ARID4A"] = {14, 58765221, 58840451};
  transcriptFactorLUT["ARID4B"] = {1, 235330209, 235491532};
  transcriptFactorLUT["ARID5A"] = {2, 97202463, 97218371};
  transcriptFactorLUT["ARID5B"] = {10, 63661012, 63856707};
  transcriptFactorLUT["ARNT"] = {1, 150782180, 150849244};
  transcriptFactorLUT["ARNT2"] = {15, 80696691, 80890277};
  transcriptFactorLUT["ARNTL"] = {11, 13299273, 13408812};
  transcriptFactorLUT["ARNTL2"] = {12, 27485786, 27578746};
  transcriptFactorLUT["ARX"] = {23, 25021812, 25034065};
  transcriptFactorLUT["ASCL1"] = {12, 103351451, 103354294};
  transcriptFactorLUT["ASCL2"] = {11, 2289727, 2292182};
  transcriptFactorLUT["ASCL3"] = {11, 8959118, 8964580};
  transcriptFactorLUT["ASCL4"] = {12, 108168161, 108170421};
  transcriptFactorLUT["ATF1"] = {12, 51157788, 51214943};
  transcriptFactorLUT["ATF2"] = {2, 175936977, 176032934};
  transcriptFactorLUT["ATF3"] = {1, 212738675, 212794119};
  transcriptFactorLUT["ATF4"] = {22, 39916568, 39918691};
  transcriptFactorLUT["ATF5"] = {19, 50431958, 50437193};
  transcriptFactorLUT["ATF6"] = {1, 161736033, 161933860};
  transcriptFactorLUT["ATF6B"] = {6, 32083044, 32096017};
  transcriptFactorLUT["ATF7"] = {12, 53905842, 54020199};
  transcriptFactorLUT["ATOH1"] = {4, 94750077, 94751142};
  transcriptFactorLUT["ATOH7"] = {10, 69990351, 69991870};
  transcriptFactorLUT["ATOH8"] = {2, 85980908, 86018506};
  transcriptFactorLUT["BACH1"] = {21, 30671115, 30718469};
  transcriptFactorLUT["BACH2"] = {6, 90636246, 91006627};
  transcriptFactorLUT["BARHL1"] = {9, 135457992, 135465640};
  transcriptFactorLUT["BARHL2"] = {1, 91177578, 91182794};
  transcriptFactorLUT["BARX1"] = {9, 96713908, 96717608};
  transcriptFactorLUT["BARX2"] = {11, 129245880, 129322174};
  transcriptFactorLUT["BATF"] = {14, 75988783, 76013334};
  transcriptFactorLUT["BATF2"] = {11, 64755416, 64757749};
  transcriptFactorLUT["BATF3"] = {1, 212859758, 212873327};
  transcriptFactorLUT["BAZ2A"] = {12, 56989379, 57030163};
  transcriptFactorLUT["BAZ2B"] = {2, 160175489, 160473112};
  transcriptFactorLUT["BBX"] = {3, 107241782, 107530176};
  transcriptFactorLUT["BCL11A"] = {2, 60684328, 60780633};
  transcriptFactorLUT["BCL11B"] = {14, 99635624, 99738050};
  transcriptFactorLUT["BCL6"] = {3, 187439164, 187454285};
  transcriptFactorLUT["BCL6B"] = {17, 6926368, 6932961};
  transcriptFactorLUT["BHLHA15"] = {7, 97841565, 97842271};
  transcriptFactorLUT["BHLHA9"] = {17, 1173857, 1174565};
  transcriptFactorLUT["BHLHE22"] = {8, 65492794, 65496191};
  transcriptFactorLUT["BHLHE23"] = {20, 61637330, 61638387};
  transcriptFactorLUT["BHLHE40"] = {3, 5021096, 5026865};
  transcriptFactorLUT["BHLHE41"] = {12, 26272958, 26278003};
  transcriptFactorLUT["BNC2"] = {9, 16409500, 16870786};
  transcriptFactorLUT["BSX"] = {11, 122848356, 122852379};
  transcriptFactorLUT["C11orf95"] = {11, 63527363, 63536113};
  transcriptFactorLUT["CAMTA1"] = {1, 6845383, 7829766};
  transcriptFactorLUT["CAMTA2"] = {17, 4871286, 4890960};
  transcriptFactorLUT["CARHSP1"] = {16, 8946798, 8962869};
  transcriptFactorLUT["CBFB"] = {16, 67063049, 67134958};
  transcriptFactorLUT["CCDC79"] = {16, 66788878, 66835523};
  transcriptFactorLUT["CDC5L"] = {6, 44355250, 44418161};
  transcriptFactorLUT["CDX1"] = {5, 149546343, 149564121};
  transcriptFactorLUT["CDX2"] = {13, 28536204, 28543505};
  transcriptFactorLUT["CDX4"] = {23, 72667089, 72674421};
  transcriptFactorLUT["CEBPA"] = {19, 33790839, 33793470};
  transcriptFactorLUT["CEBPB"] = {20, 48807119, 48809227};
  transcriptFactorLUT["CEBPD"] = {8, 48649475, 48650726};
  transcriptFactorLUT["CEBPE"] = {14, 23586514, 23588820};
  transcriptFactorLUT["CEBPG"] = {19, 33865429, 33873592};
  transcriptFactorLUT["CIC"] = {19, 42788816, 42799948};
  transcriptFactorLUT["CLOCK"] = {4, 56294067, 56413076};
  transcriptFactorLUT["CREB1"] = {2, 208394615, 208470284};
  transcriptFactorLUT["CREB3"] = {9, 35732316, 35737005};
  transcriptFactorLUT["CREB3L1"] = {11, 46299188, 46342972};
  transcriptFactorLUT["CREB3L2"] = {7, 137597556, 137686847};
  transcriptFactorLUT["CREB3L3"] = {19, 4153597, 4173051};
  transcriptFactorLUT["CREB3L4"] = {1, 153940314, 153946840};
  transcriptFactorLUT["CREB5"] = {7, 28475233, 28865511};
  transcriptFactorLUT["CREBL2"] = {12, 12764766, 12798042};
  transcriptFactorLUT["CREM"] = {10, 35456466, 35501886};
  transcriptFactorLUT["CRX"] = {19, 48325098, 48346586};
  transcriptFactorLUT["CSDC2"] = {22, 41957013, 41972670};
  transcriptFactorLUT["CSDE1"] = {1, 115259533, 115300671};
  transcriptFactorLUT["CTCF"] = {16, 67596309, 67673088};
  transcriptFactorLUT["CTCFL"] = {20, 56081797, 56100163};
  transcriptFactorLUT["CUX1"] = {7, 101459183, 101927250};
  transcriptFactorLUT["CUX2"] = {12, 111471827, 111788358};
  transcriptFactorLUT["DBP"] = {19, 49133816, 49140807};
  transcriptFactorLUT["DBX1"] = {11, 20177759, 20181870};
  transcriptFactorLUT["DBX2"] = {12, 45408538, 45444882};
  transcriptFactorLUT["DDIT3"] = {12, 57910370, 57914300};
  transcriptFactorLUT["DEAF1"] = {11, 644219, 695754};
  transcriptFactorLUT["DLX1"] = {2, 172950207, 172954401};
  transcriptFactorLUT["DLX2"] = {2, 172964165, 172967478};
  transcriptFactorLUT["DLX3"] = {17, 48067368, 48072588};
  transcriptFactorLUT["DLX4"] = {17, 48046561, 48052323};
  transcriptFactorLUT["DLX5"] = {7, 96649701, 96654143};
  transcriptFactorLUT["DLX6"] = {7, 96635289, 96640352};
  transcriptFactorLUT["DMBX1"] = {1, 46972667, 46979886};
  transcriptFactorLUT["DMRT1"] = {9, 841689, 969090};
  transcriptFactorLUT["DMRT2"] = {9, 1050353, 1057554};
  transcriptFactorLUT["DMRT3"] = {9, 976967, 991732};
  transcriptFactorLUT["DMRTA1"] = {9, 22446839, 22452472};
  transcriptFactorLUT["DMRTA2"] = {1, 50883222, 50889119};
  transcriptFactorLUT["DMRTB1"] = {1, 53925071, 53933160};
  transcriptFactorLUT["DMRTC2"] = {19, 42349085, 42356397};
  transcriptFactorLUT["DMTF1"] = {7, 86781676, 86825648};
  transcriptFactorLUT["DNAJC1"] = {10, 22045476, 22292650};
  transcriptFactorLUT["DNAJC2"] = {7, 102952920, 102985320};
  transcriptFactorLUT["DPRX"] = {19, 54135309, 54140263};
  transcriptFactorLUT["DRGX"] = {10, 50574160, 50604062};
  transcriptFactorLUT["DUXA"] = {19, 57663093, 57678856};
  transcriptFactorLUT["E2F1"] = {20, 32263291, 32274210};
  transcriptFactorLUT["E2F2"] = {1, 23832919, 23857712};
  transcriptFactorLUT["E2F3"] = {6, 20403909, 20493945};
  transcriptFactorLUT["E2F4"] = {16, 67226067, 67232821};
  transcriptFactorLUT["E2F5"] = {8, 86099909, 86126753};
  transcriptFactorLUT["E2F6"] = {2, 11584500, 11606303};
  transcriptFactorLUT["E2F7"] = {12, 77415025, 77459360};
  transcriptFactorLUT["E2F8"] = {11, 19245609, 19262507};
  transcriptFactorLUT["E4F1"] = {16, 2273488, 2285743};
  transcriptFactorLUT["EAF2"] = {3, 121554033, 121605373};
  transcriptFactorLUT["EBF1"] = {5, 158122922, 158526788};
  transcriptFactorLUT["EBF2"] = {8, 25699245, 25902640};
  transcriptFactorLUT["EBF3"] = {10, 131633495, 131762091};
  transcriptFactorLUT["EBF4"] = {20, 2673523, 2740754};
  transcriptFactorLUT["EGR1"] = {5, 137801180, 137805004};
  transcriptFactorLUT["EGR2"] = {10, 64571755, 64576126};
  transcriptFactorLUT["EGR3"] = {8, 22545173, 22550815};
  transcriptFactorLUT["EGR4"] = {2, 73518056, 73520829};
  transcriptFactorLUT["EHF"] = {11, 34642587, 34684834};
  transcriptFactorLUT["ELF1"] = {13, 41506054, 41556418};
  transcriptFactorLUT["ELF2"] = {4, 139978870, 140005568};
  transcriptFactorLUT["ELF3"] = {1, 201979689, 201986315};
  transcriptFactorLUT["ELF4"] = {23, 129198894, 129244688};
  transcriptFactorLUT["ELF5"] = {11, 34500341, 34535347};
  transcriptFactorLUT["ELK1"] = {23, 47494918, 47510003};
  transcriptFactorLUT["ELK3"] = {12, 96588159, 96663613};
  transcriptFactorLUT["ELK4"] = {1, 205588395, 205602000};
  transcriptFactorLUT["EMX1"] = {2, 73144603, 73162020};
  transcriptFactorLUT["EMX2"] = {10, 119301955, 119309057};
  transcriptFactorLUT["EN1"] = {2, 119599746, 119605759};
  transcriptFactorLUT["EN2"] = {7, 155250823, 155257526};
  transcriptFactorLUT["EOMES"] = {3, 27757439, 27764206};
  transcriptFactorLUT["ERF"] = {19, 42751712, 42759316};
  transcriptFactorLUT["ERG"] = {21, 39751949, 40033704};
  transcriptFactorLUT["ESR1"] = {6, 152011630, 152424408};
  transcriptFactorLUT["ESR2"] = {14, 64699746, 64761128};
  transcriptFactorLUT["ESRRA"] = {11, 64072999, 64084212};
  transcriptFactorLUT["ESRRB"] = {14, 76837689, 76968180};
  transcriptFactorLUT["ESRRG"] = {1, 216676587, 217113015};
  transcriptFactorLUT["ESX1"] = {23, 103494718, 103499599};
  transcriptFactorLUT["ETS1"] = {11, 128328655, 128457453};
  transcriptFactorLUT["ETS2"] = {21, 40177754, 40196878};
  transcriptFactorLUT["ETV1"] = {7, 13930855, 14026139};
  transcriptFactorLUT["ETV2"] = {19, 36132638, 36135773};
  transcriptFactorLUT["ETV3"] = {1, 157102975, 157108177};
  transcriptFactorLUT["ETV4"] = {17, 41605210, 41623800};
  transcriptFactorLUT["ETV5"] = {3, 185764105, 185826901};
  transcriptFactorLUT["ETV6"] = {12, 11802787, 12048325};
  transcriptFactorLUT["ETV7"] = {6, 36333970, 36355577};
  transcriptFactorLUT["EVX1"] = {7, 27282163, 27287438};
  transcriptFactorLUT["EVX2"] = {2, 176944834, 176948690};
  transcriptFactorLUT["FERD3L"] = {7, 19184404, 19185044};
  transcriptFactorLUT["FEV"] = {2, 219845808, 219850379};
  transcriptFactorLUT["FEZF1"] = {7, 121941362, 121944565};
  transcriptFactorLUT["FEZF2"] = {3, 62355346, 62359190};
  transcriptFactorLUT["FIGLA"] = {2, 71004441, 71017775};
  transcriptFactorLUT["FIZ1"] = {19, 56102736, 56110893};
  transcriptFactorLUT["FLI1"] = {11, 128556429, 128683162};
  transcriptFactorLUT["FOS"] = {14, 75745480, 75748937};
  transcriptFactorLUT["FOSB"] = {19, 45971252, 45978437};
  transcriptFactorLUT["FOSL1"] = {11, 65659606, 65667997};
  transcriptFactorLUT["FOSL2"] = {2, 28615778, 28637516};
  transcriptFactorLUT["FOXA1"] = {14, 38058756, 38064325};
  transcriptFactorLUT["FOXA2"] = {20, 22561641, 22565101};
  transcriptFactorLUT["FOXA3"] = {19, 46367517, 46377055};
  transcriptFactorLUT["FOXB1"] = {15, 60296420, 60298142};
  transcriptFactorLUT["FOXB2"] = {9, 79634570, 79635869};
  transcriptFactorLUT["FOXC1"] = {6, 1610680, 1614129};
  transcriptFactorLUT["FOXC2"] = {16, 86600856, 86602537};
  transcriptFactorLUT["FOXD2"] = {1, 47901688, 47906363};
  transcriptFactorLUT["FOXD3"] = {1, 63788729, 63790797};
  transcriptFactorLUT["FOXD4"] = {9, 116230, 118417};
  transcriptFactorLUT["FOXD4L1"] = {2, 114256660, 114258727};
  transcriptFactorLUT["FOXD4L3"] = {9, 70917782, 70920000};
  transcriptFactorLUT["FOXD4L4"] = {9, 42718065, 42719316};
  transcriptFactorLUT["FOXD4L5"] = {9, 70175706, 70178815};
  transcriptFactorLUT["FOXD4L6"] = {9, 69199479, 69202204};
  transcriptFactorLUT["FOXE1"] = {9, 100615536, 100618997};
  transcriptFactorLUT["FOXE3"] = {1, 47881743, 47883724};
  transcriptFactorLUT["FOXF1"] = {16, 86544132, 86548070};
  transcriptFactorLUT["FOXF2"] = {6, 1390068, 1395832};
  transcriptFactorLUT["FOXG1"] = {14, 29236277, 29239483};
  transcriptFactorLUT["FOXH1"] = {8, 145699114, 145701718};
  transcriptFactorLUT["FOXI1"] = {5, 169532916, 169536729};
  transcriptFactorLUT["FOXI2"] = {10, 129535537, 129539450};
  transcriptFactorLUT["FOXJ1"] = {17, 74132414, 74137380};
  transcriptFactorLUT["FOXJ2"] = {12, 8185358, 8208118};
  transcriptFactorLUT["FOXJ3"] = {1, 42642209, 42800903};
  transcriptFactorLUT["FOXK1"] = {7, 4721929, 4811074};
  transcriptFactorLUT["FOXK2"] = {17, 80477593, 80562483};
  transcriptFactorLUT["FOXL1"] = {16, 86612114, 86615304};
  transcriptFactorLUT["FOXL2"] = {3, 138663065, 138665982};
  transcriptFactorLUT["FOXM1"] = {12, 2966846, 2986321};
  transcriptFactorLUT["FOXN1"] = {17, 26850958, 26865175};
  transcriptFactorLUT["FOXN2"] = {2, 48541794, 48606434};
  transcriptFactorLUT["FOXN3"] = {14, 89622515, 89883454};
  transcriptFactorLUT["FOXN4"] = {12, 109715782, 109747025};
  transcriptFactorLUT["FOXO1"] = {13, 41129800, 41240734};
  transcriptFactorLUT["FOXO3"] = {6, 108881025, 109005971};
  transcriptFactorLUT["FOXO4"] = {23, 70315998, 70323384};
  transcriptFactorLUT["FOXO6"] = {1, 41827602, 41849263};
  transcriptFactorLUT["FOXP1"] = {3, 71247033, 71633140};
  transcriptFactorLUT["FOXP2"] = {7, 114055051, 114333827};
  transcriptFactorLUT["FOXP3"] = {23, 49106896, 49121288};
  transcriptFactorLUT["FOXP4"] = {6, 41514163, 41570122};
  transcriptFactorLUT["FOXQ1"] = {6, 1312674, 1314993};
  transcriptFactorLUT["FOXR1"] = {11, 118842416, 118851995};
  transcriptFactorLUT["FOXR2"] = {23, 55649832, 55652621};
  transcriptFactorLUT["FOXS1"] = {20, 30432102, 30433420};
  transcriptFactorLUT["GABPA"] = {21, 27107257, 27144771};
  transcriptFactorLUT["GATA1"] = {23, 48644981, 48652717};
  transcriptFactorLUT["GATA2"] = {3, 128198264, 128212030};
  transcriptFactorLUT["GATA3"] = {10, 8096666, 8117164};
  transcriptFactorLUT["GATA4"] = {8, 11561716, 11617509};
  transcriptFactorLUT["GATA5"] = {20, 61038552, 61051026};
  transcriptFactorLUT["GATA6"] = {18, 19749397, 19782491};
  transcriptFactorLUT["GATAD1"] = {7, 92076761, 92089381};
  transcriptFactorLUT["GATAD2B"] = {1, 153777202, 153895451};
  transcriptFactorLUT["GBX1"] = {7, 150845675, 150864635};
  transcriptFactorLUT["GBX2"] = {2, 237073878, 237076652};
  transcriptFactorLUT["GCM1"] = {6, 52991759, 53013624};
  transcriptFactorLUT["GCM2"] = {6, 10873455, 10882098};
  transcriptFactorLUT["GFI1"] = {1, 92940317, 92951628};
  transcriptFactorLUT["GFI1B"] = {9, 135854097, 135867084};
  transcriptFactorLUT["GLI1"] = {12, 57853917, 57866047};
  transcriptFactorLUT["GLI2"] = {2, 121554866, 121750229};
  transcriptFactorLUT["GLI3"] = {7, 42000547, 42276618};
  transcriptFactorLUT["GLI4"] = {8, 144349606, 144359101};
  transcriptFactorLUT["GLIS1"] = {1, 53971905, 54199877};
  transcriptFactorLUT["GLIS2"] = {16, 4382224, 4389598};
  transcriptFactorLUT["GLIS3"] = {9, 3824127, 4300035};
  transcriptFactorLUT["GMEB1"] = {1, 28995239, 29042115};
  transcriptFactorLUT["GMEB2"] = {20, 62218954, 62258381};
  transcriptFactorLUT["GRHL1"] = {2, 10091791, 10142412};
  transcriptFactorLUT["GRHL2"] = {8, 102504667, 102681952};
  transcriptFactorLUT["GRHL3"] = {1, 24649529, 24681808};
  transcriptFactorLUT["GSC"] = {14, 95234559, 95236499};
  transcriptFactorLUT["GSC2"] = {22, 19136503, 19137796};
  transcriptFactorLUT["GSX1"] = {13, 28366779, 28368089};
  transcriptFactorLUT["GSX2"] = {4, 54966247, 54968122};
  transcriptFactorLUT["GTF2I"] = {7, 74071990, 74175022};
  transcriptFactorLUT["GTF2IRD1"] = {7, 73868119, 74016920};
  transcriptFactorLUT["GTF2IRD2"] = {7, 74247755, 74267872};
  transcriptFactorLUT["GTF2IRD2B"] = {7, 74508346, 74565623};
  transcriptFactorLUT["GTF3A"] = {13, 27998680, 28009846};
  transcriptFactorLUT["GZF1"] = {20, 23345020, 23353683};
  transcriptFactorLUT["HAND1"] = {5, 153854531, 153857824};
  transcriptFactorLUT["HAND2"] = {4, 174447651, 174451378};
  transcriptFactorLUT["HBP1"] = {7, 106809405, 106842974};
  transcriptFactorLUT["HDX"] = {23, 83572881, 83757487};
  transcriptFactorLUT["HELT"] = {4, 185939994, 185941958};
  transcriptFactorLUT["HES1"] = {3, 193853930, 193856401};
  transcriptFactorLUT["HES2"] = {1, 6475293, 6479979};
  transcriptFactorLUT["HES3"] = {1, 6304251, 6305638};
  transcriptFactorLUT["HES4"] = {1, 934343, 935552};
  transcriptFactorLUT["HES5"] = {1, 2460183, 2461684};
  transcriptFactorLUT["HES6"] = {2, 239146907, 239148765};
  transcriptFactorLUT["HES7"] = {17, 8023907, 8027410};
  transcriptFactorLUT["HESX1"] = {3, 57231943, 57234280};
  transcriptFactorLUT["HEY1"] = {8, 80676244, 80680098};
  transcriptFactorLUT["HEY2"] = {6, 126070731, 126082415};
  transcriptFactorLUT["HEYL"] = {1, 40089102, 40105348};
  transcriptFactorLUT["HHEX"] = {10, 94449680, 94455408};
  transcriptFactorLUT["HIC1"] = {17, 1958392, 1962981};
  transcriptFactorLUT["HIC2"] = {22, 21771692, 21805750};
  transcriptFactorLUT["HIF1A"] = {14, 62162118, 62214977};
  transcriptFactorLUT["HINFP"] = {11, 118992232, 119005765};
  transcriptFactorLUT["HIVEP1"] = {6, 12012723, 12165232};
  transcriptFactorLUT["HIVEP2"] = {6, 143072603, 143266338};
  transcriptFactorLUT["HIVEP3"] = {1, 42312859, 42501596};
  transcriptFactorLUT["HKR1"] = {19, 37825579, 37855357};
  transcriptFactorLUT["HLF"] = {17, 53342320, 53402426};
  transcriptFactorLUT["HLX"] = {1, 221052742, 221058400};
  transcriptFactorLUT["HMBOX1"] = {8, 28747910, 28910242};
  transcriptFactorLUT["HMG20A"] = {15, 77712992, 77777946};
  transcriptFactorLUT["HMG20B"] = {19, 3572942, 3579081};
  transcriptFactorLUT["HMGA1"] = {6, 34204576, 34214008};
  transcriptFactorLUT["HMGA2"] = {12, 66218239, 66346311};
  transcriptFactorLUT["HMGB1"] = {13, 31032878, 31040081};
  transcriptFactorLUT["HMGB2"] = {4, 174252526, 174255595};
  transcriptFactorLUT["HMGB3"] = {23, 150151747, 150159248};
  transcriptFactorLUT["HMGXB3"] = {5, 149380168, 149432706};
  transcriptFactorLUT["HMGXB4"] = {22, 35653444, 35691800};
  transcriptFactorLUT["HMX1"] = {4, 8847801, 8873543};
  transcriptFactorLUT["HMX2"] = {10, 124907637, 124910188};
  transcriptFactorLUT["HMX3"] = {10, 124895566, 124897247};
  transcriptFactorLUT["HNF1A"] = {12, 121416548, 121440314};
  transcriptFactorLUT["HNF1B"] = {17, 36046433, 36105069};
  transcriptFactorLUT["HNF4A"] = {20, 42984440, 43061485};
  transcriptFactorLUT["HNF4G"] = {8, 76452202, 76479061};
  transcriptFactorLUT["HOMEZ"] = {14, 23742843, 23755309};
  transcriptFactorLUT["HOPX"] = {4, 57514153, 57547872};
  transcriptFactorLUT["HOXA1"] = {7, 27132613, 27135625};
  transcriptFactorLUT["HOXA10"] = {7, 27210209, 27213955};
  transcriptFactorLUT["HOXA11"] = {7, 27220775, 27224835};
  transcriptFactorLUT["HOXA13"] = {7, 27236498, 27239725};
  transcriptFactorLUT["HOXA2"] = {7, 27139972, 27142394};
  transcriptFactorLUT["HOXA3"] = {7, 27145808, 27166639};
  transcriptFactorLUT["HOXA4"] = {7, 27168125, 27170399};
  transcriptFactorLUT["HOXA5"] = {7, 27180670, 27183287};
  transcriptFactorLUT["HOXA6"] = {7, 27185201, 27187393};
  transcriptFactorLUT["HOXA7"] = {7, 27193337, 27196296};
  transcriptFactorLUT["HOXA9"] = {7, 27202056, 27205149};
  transcriptFactorLUT["HOXB1"] = {17, 46606806, 46608272};
  transcriptFactorLUT["HOXB13"] = {17, 46802126, 46806111};
  transcriptFactorLUT["HOXB2"] = {17, 46620018, 46622393};
  transcriptFactorLUT["HOXB3"] = {17, 46626231, 46651810};
  transcriptFactorLUT["HOXB4"] = {17, 46652868, 46655743};
  transcriptFactorLUT["HOXB5"] = {17, 46668618, 46671103};
  transcriptFactorLUT["HOXB6"] = {17, 46673098, 46682334};
  transcriptFactorLUT["HOXB7"] = {17, 46684594, 46688383};
  transcriptFactorLUT["HOXB8"] = {17, 46689707, 46692301};
  transcriptFactorLUT["HOXB9"] = {17, 46698518, 46703835};
  transcriptFactorLUT["HOXC10"] = {12, 54378945, 54384062};
  transcriptFactorLUT["HOXC11"] = {12, 54366909, 54370203};
  transcriptFactorLUT["HOXC12"] = {12, 54348713, 54350350};
  transcriptFactorLUT["HOXC13"] = {12, 54332575, 54340328};
  transcriptFactorLUT["HOXC4"] = {12, 54410635, 54449814};
  transcriptFactorLUT["HOXC5"] = {12, 54426831, 54429144};
  transcriptFactorLUT["HOXC6"] = {12, 54422193, 54424607};
  transcriptFactorLUT["HOXC8"] = {12, 54402889, 54406545};
  transcriptFactorLUT["HOXC9"] = {12, 54393876, 54397120};
  transcriptFactorLUT["HOXD1"] = {2, 177053306, 177055635};
  transcriptFactorLUT["HOXD10"] = {2, 176981491, 176984670};
  transcriptFactorLUT["HOXD11"] = {2, 176972083, 176974316};
  transcriptFactorLUT["HOXD12"] = {2, 176964529, 176965488};
  transcriptFactorLUT["HOXD13"] = {2, 176957531, 176960666};
  transcriptFactorLUT["HOXD3"] = {2, 177028804, 177037826};
  transcriptFactorLUT["HOXD4"] = {2, 177016112, 177017949};
  transcriptFactorLUT["HOXD8"] = {2, 176994467, 176997423};
  transcriptFactorLUT["HOXD9"] = {2, 176987412, 176989645};
  transcriptFactorLUT["HSF1"] = {8, 145515269, 145538385};
  transcriptFactorLUT["HSF4"] = {16, 67197287, 67203848};
  transcriptFactorLUT["HSFX1"] = {23, 148674182, 148676974};
  transcriptFactorLUT["HSFX2"] = {23, 148674171, 148676970};
  transcriptFactorLUT["HSFY1"] = {24, 20708576, 20750849};
  transcriptFactorLUT["HSFY2"] = {24, 20708556, 20750849};
  transcriptFactorLUT["ID1"] = {20, 30193085, 30194317};
  transcriptFactorLUT["ID2"] = {2, 8822112, 8824583};
  transcriptFactorLUT["ID3"] = {1, 23884420, 23886285};
  transcriptFactorLUT["ID4"] = {6, 19837600, 19842431};
  transcriptFactorLUT["IKZF1"] = {7, 50344264, 50472798};
  transcriptFactorLUT["IKZF2"] = {2, 213864410, 214016333};
  transcriptFactorLUT["IKZF3"] = {17, 37913967, 38020441};
  transcriptFactorLUT["IKZF4"] = {12, 56414688, 56432219};
  transcriptFactorLUT["IKZF5"] = {10, 124750321, 124768366};
  transcriptFactorLUT["INSM1"] = {20, 20348764, 20351592};
  transcriptFactorLUT["INSM2"] = {14, 36003247, 36006260};
  transcriptFactorLUT["IRF1"] = {5, 131817300, 131826465};
  transcriptFactorLUT["IRF2"] = {4, 185308875, 185395726};
  transcriptFactorLUT["IRF3"] = {19, 50162825, 50169132};
  transcriptFactorLUT["IRF4"] = {6, 391738, 411443};
  transcriptFactorLUT["IRF5"] = {7, 128578270, 128590096};
  transcriptFactorLUT["IRF6"] = {1, 209958967, 209979520};
  transcriptFactorLUT["IRF7"] = {11, 612554, 615999};
  transcriptFactorLUT["IRF8"] = {16, 85932773, 85956211};
  transcriptFactorLUT["IRF9"] = {14, 24630421, 24635774};
  transcriptFactorLUT["IRX1"] = {5, 3596167, 3601517};
  transcriptFactorLUT["IRX2"] = {5, 2746278, 2751769};
  transcriptFactorLUT["IRX3"] = {16, 54317211, 54320378};
  transcriptFactorLUT["IRX4"] = {5, 1877540, 1887098};
  transcriptFactorLUT["IRX5"] = {16, 54965110, 54968395};
  transcriptFactorLUT["IRX6"] = {16, 55358470, 55364672};
  transcriptFactorLUT["ISL1"] = {5, 50678957, 50690563};
  transcriptFactorLUT["ISL2"] = {15, 76629064, 76634816};
  transcriptFactorLUT["ISX"] = {22, 35462128, 35483380};
  transcriptFactorLUT["JARID2"] = {6, 15246205, 15522273};
  transcriptFactorLUT["JDP2"] = {14, 75898836, 75939404};
  transcriptFactorLUT["JUN"] = {1, 59246462, 59249785};
  transcriptFactorLUT["JUNB"] = {19, 12902309, 12904125};
  transcriptFactorLUT["JUND"] = {19, 18390503, 18392466};
  transcriptFactorLUT["KDM5A"] = {12, 389222, 498621};
  transcriptFactorLUT["KDM5B"] = {1, 202696531, 202777549};
  transcriptFactorLUT["KDM5C"] = {23, 53220502, 53254604};
  transcriptFactorLUT["KDM5D"] = {24, 21867300, 21906825};
  transcriptFactorLUT["KIAA2018"] = {3, 113367232, 113415493};
  transcriptFactorLUT["KLF1"] = {19, 12995236, 12998017};
  transcriptFactorLUT["KLF10"] = {8, 103661004, 103666192};
  transcriptFactorLUT["KLF11"] = {2, 10184371, 10194963};
  transcriptFactorLUT["KLF12"] = {13, 74260148, 74708066};
  transcriptFactorLUT["KLF13"] = {15, 31619057, 31727868};
  transcriptFactorLUT["KLF14"] = {7, 130417381, 130418860};
  transcriptFactorLUT["KLF15"] = {3, 126061477, 126076236};
  transcriptFactorLUT["KLF16"] = {19, 1852397, 1863564};
  transcriptFactorLUT["KLF17"] = {1, 44584521, 44600809};
  transcriptFactorLUT["KLF2"] = {19, 16435650, 16438339};
  transcriptFactorLUT["KLF3"] = {4, 38665789, 38703129};
  transcriptFactorLUT["KLF4"] = {9, 110247132, 110252047};
  transcriptFactorLUT["KLF5"] = {13, 73629113, 73651680};
  transcriptFactorLUT["KLF6"] = {10, 3818187, 3827473};
  transcriptFactorLUT["KLF7"] = {2, 207938861, 208031970};
  transcriptFactorLUT["KLF8"] = {23, 56258869, 56314322};
  transcriptFactorLUT["KLF9"] = {9, 72999512, 73029573};
  transcriptFactorLUT["L3MBTL1"] = {20, 42143075, 42170535};
  transcriptFactorLUT["L3MBTL4"] = {18, 5954704, 6414910};
  transcriptFactorLUT["LBX1"] = {10, 102986732, 102988717};
  transcriptFactorLUT["LBX2"] = {2, 74724643, 74726685};
  transcriptFactorLUT["LCOR"] = {10, 98592016, 98724198};
  transcriptFactorLUT["LCORL"] = {4, 17844838, 18023483};
  transcriptFactorLUT["LEF1"] = {4, 108968700, 109087953};
  transcriptFactorLUT["LEUTX"] = {19, 40267233, 40276775};
  transcriptFactorLUT["LHX1"] = {17, 35294771, 35301915};
  transcriptFactorLUT["LHX2"] = {9, 126773888, 126795442};
  transcriptFactorLUT["LHX3"] = {9, 139088095, 139095004};
  transcriptFactorLUT["LHX4"] = {1, 180199432, 180244188};
  transcriptFactorLUT["LHX5"] = {12, 113900693, 113909877};
  transcriptFactorLUT["LHX6"] = {9, 124964855, 124991091};
  transcriptFactorLUT["LHX8"] = {1, 75600566, 75627218};
  transcriptFactorLUT["LHX9"] = {1, 197881634, 197899273};
  transcriptFactorLUT["LIN28A"] = {1, 26737268, 26756219};
  transcriptFactorLUT["LIN28B"] = {6, 105404922, 105531207};
  transcriptFactorLUT["LITAF"] = {16, 11641577, 11680806};
  transcriptFactorLUT["LMX1A"] = {1, 165171103, 165325952};
  transcriptFactorLUT["LMX1B"] = {9, 129376721, 129463311};
  transcriptFactorLUT["LYL1"] = {19, 13209841, 13213974};
  transcriptFactorLUT["MAF"] = {16, 79627744, 79634622};
  transcriptFactorLUT["MAFA"] = {8, 144510229, 144512602};
  transcriptFactorLUT["MAFB"] = {20, 39314487, 39317880};
  transcriptFactorLUT["MAFF"] = {22, 38597938, 38612517};
  transcriptFactorLUT["MAFG"] = {17, 79876144, 79881444};
  transcriptFactorLUT["MAFK"] = {7, 1570367, 1582679};
  transcriptFactorLUT["MAX"] = {14, 65472818, 65569262};
  transcriptFactorLUT["MAZ"] = {16, 29817857, 29822504};
  transcriptFactorLUT["MBD1"] = {18, 47793251, 47808144};
  transcriptFactorLUT["MBD2"] = {18, 51729049, 51751158};
  transcriptFactorLUT["MBD3"] = {19, 1576669, 1592760};
  transcriptFactorLUT["MBD4"] = {3, 129149786, 129159022};
  transcriptFactorLUT["MECOM"] = {3, 168801286, 168865522};
  transcriptFactorLUT["MECP2"] = {23, 153295685, 153363188};
  transcriptFactorLUT["MEF2A"] = {15, 100106132, 100256629};
  transcriptFactorLUT["MEF2B"] = {19, 19256375, 19281098};
  transcriptFactorLUT["MEF2D"] = {1, 156433512, 156460391};
  transcriptFactorLUT["MEIS1"] = {2, 66662531, 66799891};
  transcriptFactorLUT["MEIS2"] = {15, 37183221, 37392341};
  transcriptFactorLUT["MEIS3"] = {19, 47906374, 47922785};
  transcriptFactorLUT["MEOX1"] = {17, 41717757, 41738931};
  transcriptFactorLUT["MEOX2"] = {7, 15650836, 15726308};
  transcriptFactorLUT["MESP1"] = {15, 90293097, 90294540};
  transcriptFactorLUT["MESP2"] = {15, 90319588, 90321982};
  transcriptFactorLUT["MGA"] = {15, 41952609, 42062141};
  transcriptFactorLUT["MIER3"] = {5, 56215428, 56247957};
  transcriptFactorLUT["MITF"] = {3, 69985750, 70017488};
  transcriptFactorLUT["MIXL1"] = {1, 226411318, 226414756};
  transcriptFactorLUT["MLX"] = {17, 40719077, 40725221};
  transcriptFactorLUT["MLXIP"] = {12, 122516633, 122631894};
  transcriptFactorLUT["MLXIPL"] = {7, 73007523, 73038870};
  transcriptFactorLUT["MNT"] = {17, 2287353, 2304258};
  transcriptFactorLUT["MNX1"] = {7, 156797546, 156803347};
  transcriptFactorLUT["MSC"] = {8, 72753776, 72756731};
  transcriptFactorLUT["MSGN1"] = {2, 17997785, 17998367};
  transcriptFactorLUT["MSX1"] = {4, 4861391, 4865660};
  transcriptFactorLUT["MSX2"] = {5, 174151574, 174157902};
  transcriptFactorLUT["MTA1"] = {14, 105886185, 105937057};
  transcriptFactorLUT["MTA2"] = {11, 62360674, 62369312};
  transcriptFactorLUT["MTA3"] = {2, 42795184, 42984087};
  transcriptFactorLUT["MTF1"] = {1, 38275238, 38325292};
  transcriptFactorLUT["MXD1"] = {2, 70142172, 70170076};
  transcriptFactorLUT["MXD3"] = {5, 176734205, 176739292};
  transcriptFactorLUT["MXD4"] = {4, 2249159, 2263739};
  transcriptFactorLUT["MXI1"] = {10, 111985761, 112047123};
  transcriptFactorLUT["MYB"] = {6, 135502452, 135540311};
  transcriptFactorLUT["MYBL1"] = {8, 67474409, 67525484};
  transcriptFactorLUT["MYBL2"] = {20, 42295658, 42345136};
  transcriptFactorLUT["MYC"] = {8, 128748314, 128753680};
  transcriptFactorLUT["MYCN"] = {2, 16080559, 16087129};
  transcriptFactorLUT["MYF5"] = {12, 81110707, 81113447};
  transcriptFactorLUT["MYF6"] = {12, 81101407, 81103256};
  transcriptFactorLUT["MYNN"] = {3, 169490852, 169507504};
  transcriptFactorLUT["MYOD1"] = {11, 17741109, 17743678};
  transcriptFactorLUT["MYOG"] = {1, 203052256, 203055166};
  transcriptFactorLUT["MYSM1"] = {1, 59120410, 59165747};
  transcriptFactorLUT["MYT1"] = {20, 62795826, 62873606};
  transcriptFactorLUT["MYT1L"] = {2, 1792884, 2335085};
  transcriptFactorLUT["MZF1"] = {19, 59073283, 59084942};
  transcriptFactorLUT["NANOG"] = {12, 7941991, 7948657};
  transcriptFactorLUT["NCOR1"] = {17, 15933407, 16118874};
  transcriptFactorLUT["NCOR2"] = {12, 124808956, 125052010};
  transcriptFactorLUT["NEUROD1"] = {2, 182540832, 182545392};
  transcriptFactorLUT["NEUROD2"] = {17, 37760020, 37764175};
  transcriptFactorLUT["NEUROD4"] = {12, 55413728, 55423801};
  transcriptFactorLUT["NEUROD6"] = {7, 31377074, 31380538};
  transcriptFactorLUT["NEUROG1"] = {5, 134869971, 134871639};
  transcriptFactorLUT["NEUROG2"] = {4, 113434671, 113437328};
  transcriptFactorLUT["NEUROG3"] = {10, 71331790, 71333210};
  transcriptFactorLUT["NFAT5"] = {16, 69599868, 69738569};
  transcriptFactorLUT["NFATC1"] = {18, 77155771, 77289323};
  transcriptFactorLUT["NFATC2"] = {20, 50003493, 50159258};
  transcriptFactorLUT["NFATC3"] = {16, 68119268, 68263162};
  transcriptFactorLUT["NFATC4"] = {14, 24836116, 24848811};
  transcriptFactorLUT["NFE2"] = {12, 54685890, 54689563};
  transcriptFactorLUT["NFE2L1"] = {17, 46125685, 46138907};
  transcriptFactorLUT["NFE2L2"] = {2, 178095030, 178129859};
  transcriptFactorLUT["NFE2L3"] = {7, 26191846, 26226756};
  transcriptFactorLUT["NFIA"] = {1, 61547979, 61928460};
  transcriptFactorLUT["NFIB"] = {9, 14081841, 14398982};
  transcriptFactorLUT["NFIC"] = {19, 3359560, 3469215};
  transcriptFactorLUT["NFIL3"] = {9, 94171326, 94186908};
  transcriptFactorLUT["NFIX"] = {19, 13106583, 13209610};
  transcriptFactorLUT["NFKB1"] = {4, 103422485, 103538459};
  transcriptFactorLUT["NFKB2"] = {10, 104154334, 104162286};
  transcriptFactorLUT["NFX1"] = {9, 33290417, 33348721};
  transcriptFactorLUT["NFXL1"] = {4, 47849249, 47916574};
  transcriptFactorLUT["NFYA"] = {6, 41040706, 41070146};
  transcriptFactorLUT["NFYB"] = {12, 104510857, 104532040};
  transcriptFactorLUT["NFYC"] = {1, 41157241, 41237275};
  transcriptFactorLUT["NHLH1"] = {1, 160336860, 160342638};
  transcriptFactorLUT["NHLH2"] = {1, 116378998, 116383333};
  transcriptFactorLUT["NKX1-1"] = {4, 1396719, 1400230};
  transcriptFactorLUT["NKX2-1"] = {14, 36985603, 36988903};
  transcriptFactorLUT["NKX2-2"] = {20, 21491659, 21494664};
  transcriptFactorLUT["NKX2-3"] = {10, 101292689, 101296280};
  transcriptFactorLUT["NKX2-4"] = {20, 21376004, 21378047};
  transcriptFactorLUT["NKX2-5"] = {5, 172659106, 172662315};
  transcriptFactorLUT["NKX2-6"] = {8, 23559963, 23564111};
  transcriptFactorLUT["NKX2-8"] = {14, 37049215, 37051786};
  transcriptFactorLUT["NKX3-1"] = {8, 23536205, 23540434};
  transcriptFactorLUT["NKX3-2"] = {4, 13542453, 13546114};
  transcriptFactorLUT["NKX6-1"] = {4, 85414435, 85419387};
  transcriptFactorLUT["NKX6-2"] = {10, 134598319, 134599537};
  transcriptFactorLUT["NKX6-3"] = {8, 41503828, 41504878};
  transcriptFactorLUT["NOBOX"] = {7, 144094332, 144107320};
  transcriptFactorLUT["NOTO"] = {2, 73429385, 73438340};
  transcriptFactorLUT["NPAS1"] = {19, 47524142, 47549017};
  transcriptFactorLUT["NPAS2"] = {2, 101436612, 101613287};
  transcriptFactorLUT["NPAS3"] = {14, 33408458, 34273382};
  transcriptFactorLUT["NPAS4"] = {11, 66188474, 66194177};
  transcriptFactorLUT["NR0B1"] = {23, 30322538, 30327495};
  transcriptFactorLUT["NR0B2"] = {1, 27237974, 27240567};
  transcriptFactorLUT["NR1D1"] = {17, 38249036, 38256978};
  transcriptFactorLUT["NR1D2"] = {3, 23987611, 24022109};
  transcriptFactorLUT["NR1H2"] = {19, 50879679, 50886285};
  transcriptFactorLUT["NR1H3"] = {11, 47279467, 47290584};
  transcriptFactorLUT["NR1H4"] = {12, 100867550, 100957645};
  transcriptFactorLUT["NR1I2"] = {3, 119501556, 119537332};
  transcriptFactorLUT["NR1I3"] = {1, 161199455, 161208000};
  transcriptFactorLUT["NR2C1"] = {12, 95414004, 95467404};
  transcriptFactorLUT["NR2C2"] = {3, 14989090, 15090786};
  transcriptFactorLUT["NR2E1"] = {6, 108487261, 108510013};
  transcriptFactorLUT["NR2E3"] = {15, 72102887, 72107270};
  transcriptFactorLUT["NR2F1"] = {5, 92919042, 92930315};
  transcriptFactorLUT["NR2F2"] = {15, 96874110, 96883492};
  transcriptFactorLUT["NR2F6"] = {19, 17342693, 17356151};
  transcriptFactorLUT["NR3C1"] = {5, 142657495, 142783254};
  transcriptFactorLUT["NR3C2"] = {4, 148999914, 149363672};
  transcriptFactorLUT["NR4A1"] = {12, 52416615, 52453291};
  transcriptFactorLUT["NR4A2"] = {2, 157180943, 157189287};
  transcriptFactorLUT["NR4A3"] = {9, 102584136, 102596341};
  transcriptFactorLUT["NR5A1"] = {9, 127243514, 127269699};
  transcriptFactorLUT["NR5A2"] = {1, 199996729, 200146550};
  transcriptFactorLUT["NR6A1"] = {9, 127279553, 127533589};
  transcriptFactorLUT["NRF1"] = {7, 129269918, 129396922};
  transcriptFactorLUT["NRL"] = {14, 24549315, 24553832};
  transcriptFactorLUT["OLIG1"] = {21, 34442449, 34444728};
  transcriptFactorLUT["OLIG2"] = {21, 34398215, 34401503};
  transcriptFactorLUT["OLIG3"] = {6, 137813335, 137815531};
  transcriptFactorLUT["ONECUT1"] = {15, 53049159, 53082209};
  transcriptFactorLUT["ONECUT2"] = {18, 55102916, 55158530};
  transcriptFactorLUT["ONECUT3"] = {19, 1753661, 1775444};
  transcriptFactorLUT["OSR1"] = {2, 19551245, 19558372};
  transcriptFactorLUT["OSR2"] = {8, 99956630, 99964332};
  transcriptFactorLUT["OTP"] = {5, 76924536, 76934522};
  transcriptFactorLUT["OTX1"] = {2, 63277936, 63284966};
  transcriptFactorLUT["OTX2"] = {14, 57267424, 57277194};
  transcriptFactorLUT["OVOL1"] = {11, 65554504, 65564690};
  transcriptFactorLUT["OVOL2"] = {20, 18004795, 18039832};
  transcriptFactorLUT["OVOL3"] = {19, 36602104, 36604613};
  transcriptFactorLUT["PATZ1"] = {22, 31721789, 31742249};
  transcriptFactorLUT["PAX1"] = {20, 21686296, 21699124};
  transcriptFactorLUT["PAX2"] = {10, 102495465, 102589698};
  transcriptFactorLUT["PAX3"] = {2, 223064605, 223163715};
  transcriptFactorLUT["PAX4"] = {7, 127250345, 127255780};
  transcriptFactorLUT["PAX5"] = {9, 36833271, 37034476};
  transcriptFactorLUT["PAX6"] = {11, 31806339, 31832901};
  transcriptFactorLUT["PAX7"] = {1, 18957499, 19062632};
  transcriptFactorLUT["PAX8"] = {2, 113973573, 114036498};
  transcriptFactorLUT["PAX9"] = {14, 37126772, 37147011};
  transcriptFactorLUT["PBRM1"] = {3, 52579367, 52719866};
  transcriptFactorLUT["PBX1"] = {1, 164528596, 164821060};
  transcriptFactorLUT["PBX2"] = {6, 32152509, 32157963};
  transcriptFactorLUT["PBX3"] = {9, 128509616, 128729655};
  transcriptFactorLUT["PBX4"] = {19, 19672515, 19729725};
  transcriptFactorLUT["PDX1"] = {13, 28494167, 28500451};
  transcriptFactorLUT["PEG3"] = {19, 57321444, 57352094};
  transcriptFactorLUT["PGR"] = {11, 100900354, 100999794};
  transcriptFactorLUT["PHOX2A"] = {11, 71950120, 71955220};
  transcriptFactorLUT["PHOX2B"] = {4, 41746098, 41750987};
  transcriptFactorLUT["PIAS1"] = {15, 68346571, 68480404};
  transcriptFactorLUT["PIAS2"] = {18, 44390022, 44497495};
  transcriptFactorLUT["PIAS3"] = {1, 145575987, 145586546};
  transcriptFactorLUT["PIAS4"] = {19, 4007595, 4039384};
  transcriptFactorLUT["PINX1"] = {8, 10622470, 10697409};
  transcriptFactorLUT["PITX1"] = {5, 134363423, 134369964};
  transcriptFactorLUT["PITX2"] = {4, 111538579, 111563279};
  transcriptFactorLUT["PITX3"] = {10, 103989945, 104001231};
  transcriptFactorLUT["PKNOX1"] = {21, 44394619, 44454041};
  transcriptFactorLUT["PKNOX2"] = {11, 125034558, 125303285};
  transcriptFactorLUT["PLAG1"] = {8, 57073467, 57123859};
  transcriptFactorLUT["PLAGL1"] = {6, 144261436, 144329541};
  transcriptFactorLUT["PLAGL2"] = {20, 30780306, 30795546};
  transcriptFactorLUT["PMS1"] = {2, 190648810, 190742355};
  transcriptFactorLUT["POU1F1"] = {3, 87308782, 87325737};
  transcriptFactorLUT["POU2F1"] = {1, 167190065, 167396582};
  transcriptFactorLUT["POU2F2"] = {19, 42590261, 42636625};
  transcriptFactorLUT["POU2F3"] = {11, 120107348, 120190653};
  transcriptFactorLUT["POU3F1"] = {1, 38509522, 38512450};
  transcriptFactorLUT["POU3F2"] = {6, 99282579, 99286666};
  transcriptFactorLUT["POU3F3"] = {2, 105470525, 105475031};
  transcriptFactorLUT["POU3F4"] = {23, 82763268, 82764775};
  transcriptFactorLUT["POU4F1"] = {13, 79173229, 79177695};
  transcriptFactorLUT["POU4F2"] = {4, 147560044, 147563623};
  transcriptFactorLUT["POU4F3"] = {5, 145718586, 145720083};
  transcriptFactorLUT["POU5F1"] = {6, 31132113, 31134947};
  transcriptFactorLUT["POU5F1B"] = {8, 128427856, 128429441};
  transcriptFactorLUT["POU6F1"] = {12, 51580718, 51591950};
  transcriptFactorLUT["POU6F2"] = {7, 39017608, 39504390};
  transcriptFactorLUT["PPARA"] = {22, 46546498, 46639653};
  transcriptFactorLUT["PPARD"] = {6, 35310334, 35395968};
  transcriptFactorLUT["PPARG"] = {3, 12329348, 12475855};
  transcriptFactorLUT["PRDM1"] = {6, 106534194, 106557814};
  transcriptFactorLUT["PRDM10"] = {11, 129769600, 129872730};
  transcriptFactorLUT["PRDM12"] = {9, 133539980, 133558384};
  transcriptFactorLUT["PRDM13"] = {6, 100054649, 100063454};
  transcriptFactorLUT["PRDM14"] = {8, 70963885, 70983562};
  transcriptFactorLUT["PRDM15"] = {21, 43218384, 43283411};
  transcriptFactorLUT["PRDM16"] = {1, 2985741, 3355185};
  transcriptFactorLUT["PRDM2"] = {1, 14075875, 14114574};
  transcriptFactorLUT["PRDM4"] = {12, 108126642, 108154914};
  transcriptFactorLUT["PRDM5"] = {4, 121613067, 121844021};
  transcriptFactorLUT["PRDM6"] = {5, 122424840, 122523745};
  transcriptFactorLUT["PRDM8"] = {4, 81106423, 81125482};
  transcriptFactorLUT["PRDM9"] = {5, 23507723, 23528706};
  transcriptFactorLUT["PREB"] = {2, 27353624, 27357542};
  transcriptFactorLUT["PRKRIR"] = {11, 76061000, 76092009};
  transcriptFactorLUT["PROP1"] = {5, 177419235, 177423243};
  transcriptFactorLUT["PROX1"] = {1, 214161843, 214214847};
  transcriptFactorLUT["PROX2"] = {14, 75319735, 75330537};
  transcriptFactorLUT["PRRX1"] = {1, 170633312, 170708541};
  transcriptFactorLUT["PRRX2"] = {9, 132427919, 132484951};
  transcriptFactorLUT["PTF1A"] = {10, 23481459, 23483181};
  transcriptFactorLUT["RARA"] = {17, 38474472, 38513895};
  transcriptFactorLUT["RARB"] = {3, 25469753, 25639422};
  transcriptFactorLUT["RARG"] = {12, 53604349, 53614197};
  transcriptFactorLUT["RAX"] = {18, 56934266, 56940625};
  transcriptFactorLUT["RAX2"] = {19, 3769088, 3772219};
  transcriptFactorLUT["RBAK"] = {7, 5085451, 5109119};
  transcriptFactorLUT["RCOR1"] = {14, 103058995, 103196913};
  transcriptFactorLUT["RCOR2"] = {11, 63678692, 63684316};
  transcriptFactorLUT["RCOR3"] = {1, 211432707, 211486655};
  transcriptFactorLUT["REL"] = {2, 61108629, 61155291};
  transcriptFactorLUT["RELA"] = {11, 65421066, 65430443};
  transcriptFactorLUT["RELB"] = {19, 45504706, 45541456};
  transcriptFactorLUT["REPIN1"] = {7, 150065878, 150071133};
  transcriptFactorLUT["RERE"] = {1, 8412463, 8877699};
  transcriptFactorLUT["REST"] = {4, 57774041, 57802010};
  transcriptFactorLUT["RFX1"] = {19, 14072341, 14117134};
  transcriptFactorLUT["RFX2"] = {19, 5993174, 6110664};
  transcriptFactorLUT["RFX3"] = {9, 3247036, 3395596};
  transcriptFactorLUT["RFX4"] = {12, 106994914, 107156582};
  transcriptFactorLUT["RFX5"] = {1, 151313115, 151319769};
  transcriptFactorLUT["RFX6"] = {6, 117198375, 117253326};
  transcriptFactorLUT["RFX7"] = {15, 56382730, 56535483};
  transcriptFactorLUT["RFX8"] = {2, 102013822, 102091165};
  transcriptFactorLUT["RHOXF1"] = {23, 119243010, 119249847};
  transcriptFactorLUT["RHOXF2"] = {23, 119206240, 119211707};
  transcriptFactorLUT["RHOXF2B"] = {23, 119206228, 119211707};
  transcriptFactorLUT["RNF138"] = {18, 29672568, 29711524};
  transcriptFactorLUT["RORA"] = {15, 60780482, 60884707};
  transcriptFactorLUT["RORB"] = {9, 77112251, 77302117};
  transcriptFactorLUT["RORC"] = {1, 151778546, 151804348};
  transcriptFactorLUT["RREB1"] = {6, 7107829, 7252213};
  transcriptFactorLUT["RUNX1"] = {21, 36160097, 36421595};
  transcriptFactorLUT["RUNX3"] = {1, 25226001, 25291501};
  transcriptFactorLUT["RXRA"] = {9, 137218308, 137332432};
  transcriptFactorLUT["RXRB"] = {6, 33161361, 33168473};
  transcriptFactorLUT["SALL1"] = {16, 51169885, 51185183};
  transcriptFactorLUT["SALL2"] = {14, 21989230, 22005350};
  transcriptFactorLUT["SALL3"] = {18, 76740274, 76758969};
  transcriptFactorLUT["SALL4"] = {20, 50400582, 50419048};
  transcriptFactorLUT["SATB1"] = {3, 18389132, 18466829};
  transcriptFactorLUT["SATB2"] = {2, 200134222, 200322819};
  transcriptFactorLUT["SCRT1"] = {8, 145554227, 145559943};
  transcriptFactorLUT["SCRT2"] = {20, 642239, 656823};
  transcriptFactorLUT["SEBOX"] = {17, 26691289, 26692173};
  transcriptFactorLUT["SETDB1"] = {1, 150898814, 150917797};
  transcriptFactorLUT["SETDB2"] = {13, 50025687, 50069139};
  transcriptFactorLUT["SHOX"] = {23, 585078, 607558};
  transcriptFactorLUT["SHOX2"] = {3, 157813799, 157823952};
  transcriptFactorLUT["SIM1"] = {6, 100836749, 100911551};
  transcriptFactorLUT["SIX1"] = {14, 61111416, 61116155};
  transcriptFactorLUT["SIX2"] = {2, 45232323, 45236542};
  transcriptFactorLUT["SIX3"] = {2, 45169036, 45173216};
  transcriptFactorLUT["SIX4"] = {14, 61176255, 61190852};
  transcriptFactorLUT["SIX5"] = {19, 46268042, 46272497};
  transcriptFactorLUT["SIX6"] = {14, 60975937, 60978525};
  transcriptFactorLUT["SMAD1"] = {4, 146402950, 146480325};
  transcriptFactorLUT["SMAD2"] = {18, 45359465, 45456970};
  transcriptFactorLUT["SMAD3"] = {15, 67358194, 67487533};
  transcriptFactorLUT["SMAD4"] = {18, 48556582, 48611411};
  transcriptFactorLUT["SMAD5"] = {5, 135468533, 135518422};
  transcriptFactorLUT["SMAD6"] = {15, 66994673, 67074337};
  transcriptFactorLUT["SMAD7"] = {18, 46446222, 46477081};
  transcriptFactorLUT["SMAD9"] = {13, 37418967, 37494409};
  transcriptFactorLUT["SMARCC1"] = {3, 47627377, 47823405};
  transcriptFactorLUT["SMARCC2"] = {12, 56555635, 56583351};
  transcriptFactorLUT["SMARCE1"] = {17, 38783975, 38804103};
  transcriptFactorLUT["SNAI1"] = {20, 48599512, 48605420};
  transcriptFactorLUT["SNAI2"] = {8, 49830238, 49833999};
  transcriptFactorLUT["SNAI3"] = {16, 88744089, 88752882};
  transcriptFactorLUT["SNAPC4"] = {9, 139270028, 139292889};
  transcriptFactorLUT["SOHLH1"] = {9, 138585254, 138591374};
  transcriptFactorLUT["SOHLH2"] = {13, 36742344, 36788752};
  transcriptFactorLUT["SOX1"] = {13, 112721912, 112726020};
  transcriptFactorLUT["SOX10"] = {22, 38368318, 38380539};
  transcriptFactorLUT["SOX11"] = {2, 5832798, 5841517};
  transcriptFactorLUT["SOX12"] = {20, 306214, 310872};
  transcriptFactorLUT["SOX13"] = {1, 204042245, 204096871};
  transcriptFactorLUT["SOX14"] = {3, 137483133, 137485172};
  transcriptFactorLUT["SOX15"] = {17, 7491497, 7493488};
  transcriptFactorLUT["SOX17"] = {8, 55370494, 55373456};
  transcriptFactorLUT["SOX18"] = {20, 62679078, 62680979};
  transcriptFactorLUT["SOX2"] = {3, 181429711, 181432223};
  transcriptFactorLUT["SOX21"] = {13, 95361878, 95364797};
  transcriptFactorLUT["SOX3"] = {23, 139585151, 139587225};
  transcriptFactorLUT["SOX30"] = {5, 157052686, 157079428};
  transcriptFactorLUT["SOX4"] = {6, 21593971, 21598849};
  transcriptFactorLUT["SOX5"] = {12, 23685230, 24715383};
  transcriptFactorLUT["SOX6"] = {11, 15987994, 16424413};
  transcriptFactorLUT["SOX7"] = {8, 10581277, 10588084};
  transcriptFactorLUT["SOX8"] = {16, 1031807, 1036979};
  transcriptFactorLUT["SOX9"] = {17, 70117160, 70122560};
  transcriptFactorLUT["SP1"] = {12, 53774427, 53810226};
  transcriptFactorLUT["SP100"] = {2, 231280870, 231372861};
  transcriptFactorLUT["SP110"] = {2, 231041188, 231084827};
  transcriptFactorLUT["SP140"] = {2, 231090444, 231177930};
  transcriptFactorLUT["SP140L"] = {2, 231191893, 231268445};
  transcriptFactorLUT["SP2"] = {17, 45973515, 46006323};
  transcriptFactorLUT["SP3"] = {2, 174771186, 174830430};
  transcriptFactorLUT["SP4"] = {7, 21467688, 21554151};
  transcriptFactorLUT["SP5"] = {2, 171571856, 171574498};
  transcriptFactorLUT["SP6"] = {17, 45922279, 45933240};
  transcriptFactorLUT["SP7"] = {12, 53720359, 53730167};
  transcriptFactorLUT["SP8"] = {7, 20821893, 20826508};
  transcriptFactorLUT["SP9"] = {2, 175199820, 175202268};
  transcriptFactorLUT["SPDEF"] = {6, 34505578, 34524110};
  transcriptFactorLUT["SPI1"] = {11, 47376408, 47400127};
  transcriptFactorLUT["SPIB"] = {19, 50922194, 50934309};
  transcriptFactorLUT["SPIC"] = {12, 101871334, 101880775};
  transcriptFactorLUT["SPZ1"] = {5, 79615789, 79617660};
  transcriptFactorLUT["SREBF1"] = {17, 17714662, 17740325};
  transcriptFactorLUT["SREBF2"] = {22, 42229082, 42303312};
  transcriptFactorLUT["SRF"] = {6, 43139032, 43149244};
  transcriptFactorLUT["SRY"] = {24, 2654895, 2655782};
  transcriptFactorLUT["SSRP1"] = {11, 57093458, 57103351};
  transcriptFactorLUT["ST18"] = {8, 53023391, 53322439};
  transcriptFactorLUT["STAT1"] = {2, 191833761, 191878976};
  transcriptFactorLUT["STAT2"] = {12, 56735381, 56754037};
  transcriptFactorLUT["STAT3"] = {17, 40465342, 40540405};
  transcriptFactorLUT["STAT4"] = {2, 191894301, 192015986};
  transcriptFactorLUT["STAT5A"] = {17, 40439564, 40463960};
  transcriptFactorLUT["STAT5B"] = {17, 40351194, 40428424};
  transcriptFactorLUT["STAT6"] = {12, 57489186, 57505196};
  transcriptFactorLUT["SUB1"] = {5, 32585604, 32604185};
  transcriptFactorLUT["T"] = {6, 166571145, 166582157};
  transcriptFactorLUT["TADA2A"] = {17, 35767311, 35837226};
  transcriptFactorLUT["TADA2B"] = {4, 7045155, 7059677};
  transcriptFactorLUT["TAL1"] = {1, 47681961, 47698007};
  transcriptFactorLUT["TAL2"] = {9, 108424737, 108425385};
  transcriptFactorLUT["TBR1"] = {2, 162272619, 162281573};
  transcriptFactorLUT["TBX1"] = {22, 19744225, 19771112};
  transcriptFactorLUT["TBX10"] = {11, 67398773, 67407031};
  transcriptFactorLUT["TBX15"] = {1, 119425665, 119532179};
  transcriptFactorLUT["TBX18"] = {6, 85442215, 85473954};
  transcriptFactorLUT["TBX19"] = {1, 168250277, 168283664};
  transcriptFactorLUT["TBX2"] = {17, 59477256, 59486827};
  transcriptFactorLUT["TBX20"] = {7, 35242041, 35293711};
  transcriptFactorLUT["TBX21"] = {17, 45810609, 45823485};
  transcriptFactorLUT["TBX22"] = {23, 79270254, 79287268};
  transcriptFactorLUT["TBX3"] = {12, 115108058, 115121969};
  transcriptFactorLUT["TBX4"] = {17, 59533806, 59561664};
  transcriptFactorLUT["TBX5"] = {12, 114791734, 114846247};
  transcriptFactorLUT["TBX6"] = {16, 30097114, 30103205};
  transcriptFactorLUT["TCF12"] = {15, 57210832, 57580714};
  transcriptFactorLUT["TCF15"] = {20, 584636, 590910};
  transcriptFactorLUT["TCF21"] = {6, 134210258, 134216675};
  transcriptFactorLUT["TCF23"] = {2, 27371944, 27375819};
  transcriptFactorLUT["TCF3"] = {19, 1609288, 1652328};
  transcriptFactorLUT["TCF4"] = {18, 52889561, 53255860};
  transcriptFactorLUT["TCF7"] = {5, 133451349, 133483920};
  transcriptFactorLUT["TCF7L1"] = {2, 85360582, 85537511};
  transcriptFactorLUT["TCF7L2"] = {10, 114710008, 114927436};
  transcriptFactorLUT["TCFL5"] = {20, 61472365, 61493115};
  transcriptFactorLUT["TEAD1"] = {11, 12695968, 12966284};
  transcriptFactorLUT["TEAD2"] = {19, 49843852, 49865714};
  transcriptFactorLUT["TEAD3"] = {6, 35441373, 35464861};
  transcriptFactorLUT["TEAD4"] = {12, 3068747, 3149842};
  transcriptFactorLUT["TEF"] = {22, 41763336, 41795332};
  transcriptFactorLUT["TERF1"] = {8, 73921096, 73959987};
  transcriptFactorLUT["TERF2"] = {16, 69389463, 69419891};
  transcriptFactorLUT["TFAM"] = {10, 60144902, 60158990};
  transcriptFactorLUT["TFAP2A"] = {6, 10396915, 10412607};
  transcriptFactorLUT["TFAP2B"] = {6, 50786438, 50815326};
  transcriptFactorLUT["TFAP2C"] = {20, 55204357, 55214338};
  transcriptFactorLUT["TFAP2D"] = {6, 50681256, 50740746};
  transcriptFactorLUT["TFAP2E"] = {1, 36038970, 36060927};
  transcriptFactorLUT["TFAP4"] = {16, 4307186, 4323001};
  transcriptFactorLUT["TFCP2"] = {12, 51487538, 51566926};
  transcriptFactorLUT["TFCP2L1"] = {2, 121974163, 122042778};
  transcriptFactorLUT["TFDP1"] = {13, 114239055, 114295788};
  transcriptFactorLUT["TFDP2"] = {3, 141663269, 141719229};
  transcriptFactorLUT["TFDP3"] = {23, 132350696, 132352376};
  transcriptFactorLUT["TFE3"] = {23, 48886237, 48901043};
  transcriptFactorLUT["TFEB"] = {6, 41651715, 41703997};
  transcriptFactorLUT["TFEC"] = {7, 115575201, 115608367};
  transcriptFactorLUT["TGIF1"] = {18, 3453771, 3458406};
  transcriptFactorLUT["TGIF2"] = {20, 35202956, 35222355};
  transcriptFactorLUT["TGIF2LX"] = {23, 89176939, 89177882};
  transcriptFactorLUT["TGIF2LY"] = {24, 3447125, 3448082};
  transcriptFactorLUT["THAP1"] = {8, 42691816, 42698474};
  transcriptFactorLUT["THAP10"] = {15, 71173680, 71184772};
  transcriptFactorLUT["THAP11"] = {16, 67876212, 67878098};
  transcriptFactorLUT["THAP2"] = {12, 72057676, 72074428};
  transcriptFactorLUT["THAP3"] = {1, 6684924, 6693642};
  transcriptFactorLUT["THAP4"] = {2, 242523819, 242556916};
  transcriptFactorLUT["THAP5"] = {7, 108202587, 108209897};
  transcriptFactorLUT["THAP6"] = {4, 76439653, 76455236};
  transcriptFactorLUT["THAP7"] = {22, 21354060, 21356404};
  transcriptFactorLUT["THAP8"] = {19, 36525886, 36545664};
  transcriptFactorLUT["THAP9"] = {4, 83821836, 83841284};
  transcriptFactorLUT["THRA"] = {17, 38219067, 38250120};
  transcriptFactorLUT["THRB"] = {3, 24158644, 24536313};
  transcriptFactorLUT["TLX1"] = {10, 102891060, 102897546};
  transcriptFactorLUT["TLX2"] = {2, 74741595, 74744275};
  transcriptFactorLUT["TLX3"] = {5, 170736287, 170739138};
  transcriptFactorLUT["TOX"] = {8, 59717976, 60031767};
  transcriptFactorLUT["TOX2"] = {20, 42574535, 42698254};
  transcriptFactorLUT["TOX3"] = {16, 52471917, 52580806};
  transcriptFactorLUT["TOX4"] = {14, 21945334, 21967321};
  transcriptFactorLUT["TP53"] = {17, 7571719, 7590868};
  transcriptFactorLUT["TP63"] = {3, 189349215, 189615068};
  transcriptFactorLUT["TP73"] = {1, 3569128, 3652765};
  transcriptFactorLUT["TRERF1"] = {6, 42192670, 42419789};
  transcriptFactorLUT["TRPS1"] = {8, 116420723, 116681255};
  transcriptFactorLUT["TSC22D1"] = {13, 45006278, 45150701};
  transcriptFactorLUT["TSC22D2"] = {3, 150126121, 150177905};
  transcriptFactorLUT["TSC22D3"] = {23, 106956451, 106960291};
  transcriptFactorLUT["TSC22D4"] = {7, 100064141, 100076902};
  transcriptFactorLUT["TSHZ1"] = {18, 72922709, 73001905};
  transcriptFactorLUT["TSHZ3"] = {19, 31765850, 31840190};
  transcriptFactorLUT["TTF1"] = {9, 135250936, 135282238};
  transcriptFactorLUT["TUB"] = {11, 8060179, 8127654};
  transcriptFactorLUT["TULP1"] = {6, 35465650, 35480679};
  transcriptFactorLUT["TULP2"] = {19, 49384221, 49401996};
  transcriptFactorLUT["TULP3"] = {12, 3000032, 3050306};
  transcriptFactorLUT["TULP4"] = {6, 158733691, 158932856};
  transcriptFactorLUT["TWIST1"] = {7, 19155090, 19157295};
  transcriptFactorLUT["TWIST2"] = {2, 239756672, 239832244};
  transcriptFactorLUT["UBP1"] = {3, 33429828, 33481870};
  transcriptFactorLUT["UBTF"] = {17, 42282400, 42297041};
  transcriptFactorLUT["UBTFL1"] = {11, 89819117, 89820299};
  transcriptFactorLUT["UNCX"] = {7, 1272653, 1276613};
  transcriptFactorLUT["USF1"] = {1, 161009040, 161015769};
  transcriptFactorLUT["USF2"] = {19, 35759895, 35770718};
  transcriptFactorLUT["VAX1"] = {10, 118892800, 118897812};
  transcriptFactorLUT["VAX2"] = {2, 71127719, 71160575};
  transcriptFactorLUT["VDR"] = {12, 48235319, 48298814};
  transcriptFactorLUT["VENTX"] = {10, 135051407, 135055434};
  transcriptFactorLUT["VEZF1"] = {17, 56048909, 56065615};
  transcriptFactorLUT["VSX1"] = {20, 25059178, 25063015};
  transcriptFactorLUT["VSX2"] = {14, 74706174, 74729441};
  transcriptFactorLUT["VTN"] = {17, 26694298, 26697373};
  transcriptFactorLUT["WHSC1"] = {4, 1894508, 1983934};
  transcriptFactorLUT["WIZ"] = {19, 15532317, 15560762};
  transcriptFactorLUT["WT1"] = {11, 32409321, 32452363};
  transcriptFactorLUT["XBP1"] = {22, 29190547, 29196560};
  transcriptFactorLUT["YBX1"] = {1, 43148065, 43168020};
  transcriptFactorLUT["YBX2"] = {17, 7191570, 7197876};
  transcriptFactorLUT["YY1"] = {14, 100705101, 100745371};
  transcriptFactorLUT["YY2"] = {23, 21874104, 21876845};
  transcriptFactorLUT["ZBED1"] = {23, 2404454, 2419008};
  transcriptFactorLUT["ZBED2"] = {3, 111311746, 111314182};
  transcriptFactorLUT["ZBED3"] = {5, 76372531, 76383030};
  transcriptFactorLUT["ZBED4"] = {22, 50247496, 50283726};
  transcriptFactorLUT["ZBTB1"] = {14, 64971291, 65000408};
  transcriptFactorLUT["ZBTB10"] = {8, 81398447, 81438500};
  transcriptFactorLUT["ZBTB11"] = {3, 101368282, 101395988};
  transcriptFactorLUT["ZBTB12"] = {6, 31867393, 31869769};
  transcriptFactorLUT["ZBTB16"] = {11, 113930430, 114121397};
  transcriptFactorLUT["ZBTB17"] = {1, 16268363, 16302627};
  transcriptFactorLUT["ZBTB2"] = {6, 151685249, 151712835};
  transcriptFactorLUT["ZBTB20"] = {3, 114033346, 114866127};
  transcriptFactorLUT["ZBTB24"] = {6, 109802151, 109804440};
  transcriptFactorLUT["ZBTB25"] = {14, 64953554, 64970563};
  transcriptFactorLUT["ZBTB26"] = {9, 125680307, 125693830};
  transcriptFactorLUT["ZBTB3"] = {11, 62518434, 62521656};
  transcriptFactorLUT["ZBTB32"] = {19, 36203829, 36207940};
  transcriptFactorLUT["ZBTB33"] = {23, 119384609, 119392251};
  transcriptFactorLUT["ZBTB34"] = {9, 129622943, 129648156};
  transcriptFactorLUT["ZBTB37"] = {1, 173837492, 173855774};
  transcriptFactorLUT["ZBTB38"] = {3, 141043054, 141168632};
  transcriptFactorLUT["ZBTB4"] = {17, 7362684, 7382944};
  transcriptFactorLUT["ZBTB40"] = {1, 22778343, 22857650};
  transcriptFactorLUT["ZBTB41"] = {1, 197122813, 197169672};
  transcriptFactorLUT["ZBTB42"] = {14, 105267517, 105271048};
  transcriptFactorLUT["ZBTB43"] = {9, 129567284, 129600487};
  transcriptFactorLUT["ZBTB44"] = {11, 130096573, 130184607};
  transcriptFactorLUT["ZBTB45"] = {19, 59024896, 59030921};
  transcriptFactorLUT["ZBTB46"] = {20, 62375020, 62436856};
  transcriptFactorLUT["ZBTB47"] = {3, 42695175, 42709072};
  transcriptFactorLUT["ZBTB48"] = {1, 6640050, 6649340};
  transcriptFactorLUT["ZBTB49"] = {4, 4291923, 4323513};
  transcriptFactorLUT["ZBTB6"] = {9, 125670334, 125675607};
  transcriptFactorLUT["ZBTB7A"] = {19, 4045215, 4066816};
  transcriptFactorLUT["ZBTB7B"] = {1, 154975105, 154991001};
  transcriptFactorLUT["ZBTB7C"] = {18, 45553638, 45663680};
  transcriptFactorLUT["ZBTB8A"] = {1, 33004745, 33071551};
  transcriptFactorLUT["ZEB1"] = {10, 31610063, 31818742};
  transcriptFactorLUT["ZEB2"] = {2, 145141941, 145277958};
  transcriptFactorLUT["ZFAT"] = {8, 135490030, 135708801};
  transcriptFactorLUT["ZFHX3"] = {16, 72816785, 73082274};
  transcriptFactorLUT["ZFHX4"] = {8, 77593514, 77779521};
  transcriptFactorLUT["ZFP1"] = {16, 75182420, 75206132};
  transcriptFactorLUT["ZFP14"] = {19, 36825354, 36870105};
  transcriptFactorLUT["ZFP2"] = {5, 178322915, 178360210};
  transcriptFactorLUT["ZFP28"] = {19, 57050316, 57068170};
  transcriptFactorLUT["ZFP3"] = {17, 4981753, 4999669};
  transcriptFactorLUT["ZFP30"] = {19, 38123388, 38146313};
  transcriptFactorLUT["ZFP37"] = {9, 115804094, 115819071};
  transcriptFactorLUT["ZFP41"] = {8, 144328990, 144335761};
  transcriptFactorLUT["ZFP42"] = {4, 188916924, 188926203};
  transcriptFactorLUT["ZFP57"] = {6, 29640168, 29644931};
  transcriptFactorLUT["ZFP62"] = {5, 180274610, 180288286};
  transcriptFactorLUT["ZFP64"] = {20, 50700549, 50808524};
  transcriptFactorLUT["ZFP82"] = {19, 36882860, 36909550};
  transcriptFactorLUT["ZFP90"] = {16, 68573115, 68609975};
  transcriptFactorLUT["ZFP91"] = {11, 58346586, 58389023};
  transcriptFactorLUT["ZFPM1"] = {16, 88520013, 88601574};
  transcriptFactorLUT["ZFPM2"] = {8, 106331146, 106816767};
  transcriptFactorLUT["ZFX"] = {23, 24167761, 24234372};
  transcriptFactorLUT["ZFY"] = {24, 2803517, 2850547};
  transcriptFactorLUT["ZGLP1"] = {19, 10415478, 10420233};
  transcriptFactorLUT["ZHX1"] = {8, 124260689, 124286727};
  transcriptFactorLUT["ZHX2"] = {8, 123793900, 123986755};
  transcriptFactorLUT["ZHX3"] = {20, 39807088, 39928739};
  transcriptFactorLUT["ZIC1"] = {3, 147127180, 147134506};
  transcriptFactorLUT["ZIC2"] = {13, 100634025, 100639019};
  transcriptFactorLUT["ZIC3"] = {23, 136648345, 136654259};
  transcriptFactorLUT["ZIC4"] = {3, 147103834, 147124596};
  transcriptFactorLUT["ZIC5"] = {13, 100615274, 100624178};
  transcriptFactorLUT["ZIK1"] = {19, 58095627, 58103758};
  transcriptFactorLUT["ZIM2"] = {19, 57285922, 57352097};
  transcriptFactorLUT["ZIM3"] = {19, 57645463, 57656570};
  transcriptFactorLUT["ZKSCAN1"] = {7, 99613194, 99639312};
  transcriptFactorLUT["ZKSCAN2"] = {16, 25247321, 25268855};
  transcriptFactorLUT["ZKSCAN3"] = {6, 28317690, 28336954};
  transcriptFactorLUT["ZKSCAN4"] = {6, 28209482, 28220074};
  transcriptFactorLUT["ZKSCAN5"] = {7, 99102272, 99131445};
  transcriptFactorLUT["ZMIZ1"] = {10, 80828791, 81076285};
  transcriptFactorLUT["ZMIZ2"] = {7, 44795786, 44809479};
  transcriptFactorLUT["ZNF10"] = {12, 133707213, 133736049};
  transcriptFactorLUT["ZNF100"] = {19, 21906842, 21950430};
  transcriptFactorLUT["ZNF101"] = {19, 19779596, 19794315};
  transcriptFactorLUT["ZNF107"] = {7, 64126460, 64171960};
  transcriptFactorLUT["ZNF114"] = {19, 48777006, 48790865};
  transcriptFactorLUT["ZNF117"] = {7, 64434829, 64451414};
  transcriptFactorLUT["ZNF12"] = {7, 6728063, 6746566};
  transcriptFactorLUT["ZNF121"] = {19, 9676291, 9695209};
  transcriptFactorLUT["ZNF124"] = {1, 247318290, 247335319};
  transcriptFactorLUT["ZNF131"] = {5, 43121595, 43176426};
  transcriptFactorLUT["ZNF132"] = {19, 58944180, 58951589};
  transcriptFactorLUT["ZNF133"] = {20, 18268926, 18297640};
  transcriptFactorLUT["ZNF134"] = {19, 58125829, 58133636};
  transcriptFactorLUT["ZNF135"] = {19, 58570606, 58581110};
  transcriptFactorLUT["ZNF136"] = {19, 12273871, 12300064};
  transcriptFactorLUT["ZNF138"] = {7, 64254765, 64294059};
  transcriptFactorLUT["ZNF14"] = {19, 19821280, 19843921};
  transcriptFactorLUT["ZNF140"] = {12, 133656788, 133684258};
  transcriptFactorLUT["ZNF141"] = {4, 331595, 367691};
  transcriptFactorLUT["ZNF142"] = {2, 219502639, 219524355};
  transcriptFactorLUT["ZNF143"] = {11, 9482511, 9550071};
  transcriptFactorLUT["ZNF146"] = {19, 36705503, 36729675};
  transcriptFactorLUT["ZNF148"] = {3, 124944512, 125094198};
  transcriptFactorLUT["ZNF154"] = {19, 58207642, 58220579};
  transcriptFactorLUT["ZNF155"] = {19, 44488321, 44502477};
  transcriptFactorLUT["ZNF157"] = {23, 47229998, 47273098};
  transcriptFactorLUT["ZNF16"] = {8, 146155743, 146176274};
  transcriptFactorLUT["ZNF160"] = {19, 53569866, 53606687};
  transcriptFactorLUT["ZNF165"] = {6, 28048481, 28057340};
  transcriptFactorLUT["ZNF169"] = {9, 97021547, 97064111};
  transcriptFactorLUT["ZNF17"] = {19, 57922528, 57933307};
  transcriptFactorLUT["ZNF174"] = {16, 3451189, 3459364};
  transcriptFactorLUT["ZNF175"] = {19, 52074530, 52092991};
  transcriptFactorLUT["ZNF177"] = {19, 9486991, 9493293};
  transcriptFactorLUT["ZNF18"] = {17, 11880755, 11900827};
  transcriptFactorLUT["ZNF180"] = {19, 44978644, 45004575};
  transcriptFactorLUT["ZNF181"] = {19, 35225479, 35233774};
  transcriptFactorLUT["ZNF182"] = {23, 47834249, 47863377};
  transcriptFactorLUT["ZNF184"] = {6, 27418520, 27440897};
  transcriptFactorLUT["ZNF189"] = {9, 104161135, 104172942};
  transcriptFactorLUT["ZNF19"] = {16, 71507975, 71523254};
  transcriptFactorLUT["ZNF195"] = {11, 3379156, 3400452};
  transcriptFactorLUT["ZNF197"] = {3, 44666510, 44686752};
  transcriptFactorLUT["ZNF2"] = {2, 95831161, 95850064};
  transcriptFactorLUT["ZNF20"] = {19, 12242167, 12251222};
  transcriptFactorLUT["ZNF200"] = {16, 3272324, 3285456};
  transcriptFactorLUT["ZNF202"] = {11, 123594634, 123612391};
  transcriptFactorLUT["ZNF205"] = {16, 3162562, 3170518};
  transcriptFactorLUT["ZNF208"] = {19, 22148896, 22193745};
  transcriptFactorLUT["ZNF211"] = {19, 58144534, 58154147};
  transcriptFactorLUT["ZNF212"] = {7, 148936741, 148952700};
  transcriptFactorLUT["ZNF213"] = {16, 3185056, 3192805};
  transcriptFactorLUT["ZNF214"] = {11, 7020548, 7041586};
  transcriptFactorLUT["ZNF215"] = {11, 6947653, 6979278};
  transcriptFactorLUT["ZNF217"] = {20, 52183609, 52199636};
  transcriptFactorLUT["ZNF219"] = {14, 21558204, 21567173};
  transcriptFactorLUT["ZNF22"] = {10, 45496272, 45500777};
  transcriptFactorLUT["ZNF221"] = {19, 44455396, 44471752};
  transcriptFactorLUT["ZNF222"] = {19, 44529493, 44537262};
  transcriptFactorLUT["ZNF223"] = {19, 44556163, 44572147};
  transcriptFactorLUT["ZNF224"] = {19, 44598481, 44612479};
  transcriptFactorLUT["ZNF225"] = {19, 44617547, 44637255};
  transcriptFactorLUT["ZNF226"] = {19, 44669248, 44679582};
  transcriptFactorLUT["ZNF227"] = {19, 44716680, 44741421};
  transcriptFactorLUT["ZNF229"] = {19, 44930422, 44952665};
  transcriptFactorLUT["ZNF23"] = {16, 71481502, 71496155};
  transcriptFactorLUT["ZNF230"] = {19, 44507076, 44518072};
  transcriptFactorLUT["ZNF232"] = {17, 5009030, 5026397};
  transcriptFactorLUT["ZNF233"] = {19, 44764032, 44779470};
  transcriptFactorLUT["ZNF234"] = {19, 44645709, 44664462};
  transcriptFactorLUT["ZNF235"] = {19, 44790500, 44809178};
  transcriptFactorLUT["ZNF236"] = {18, 74536115, 74682682};
  transcriptFactorLUT["ZNF239"] = {10, 44051792, 44063907};
  transcriptFactorLUT["ZNF24"] = {18, 32912177, 32924426};
  transcriptFactorLUT["ZNF248"] = {10, 38065453, 38146564};
  transcriptFactorLUT["ZNF25"] = {10, 38238794, 38265453};
  transcriptFactorLUT["ZNF250"] = {8, 146102335, 146126846};
  transcriptFactorLUT["ZNF251"] = {8, 145946293, 145980970};
  transcriptFactorLUT["ZNF253"] = {19, 19976713, 20004293};
  transcriptFactorLUT["ZNF254"] = {19, 24216206, 24312769};
  transcriptFactorLUT["ZNF256"] = {19, 58452200, 58459077};
  transcriptFactorLUT["ZNF257"] = {19, 22235265, 22273903};
  transcriptFactorLUT["ZNF26"] = {12, 133562933, 133589154};
  transcriptFactorLUT["ZNF263"] = {16, 3333486, 3341459};
  transcriptFactorLUT["ZNF264"] = {19, 57702867, 57734214};
  transcriptFactorLUT["ZNF266"] = {19, 9523101, 9546254};
  transcriptFactorLUT["ZNF267"] = {16, 31885078, 31928629};
  transcriptFactorLUT["ZNF268"] = {12, 133757994, 133783697};
  transcriptFactorLUT["ZNF273"] = {7, 64363619, 64391955};
  transcriptFactorLUT["ZNF274"] = {19, 58694355, 58724928};
  transcriptFactorLUT["ZNF275"] = {23, 152599612, 152618384};
  transcriptFactorLUT["ZNF276"] = {16, 89787951, 89807332};
  transcriptFactorLUT["ZNF28"] = {19, 53300660, 53324922};
  transcriptFactorLUT["ZNF280D"] = {15, 56922373, 57026284};
  transcriptFactorLUT["ZNF281"] = {1, 200374074, 200379186};
  transcriptFactorLUT["ZNF282"] = {7, 148892553, 148923339};
  transcriptFactorLUT["ZNF283"] = {19, 44331472, 44353050};
  transcriptFactorLUT["ZNF284"] = {19, 44576296, 44591623};
  transcriptFactorLUT["ZNF285"] = {19, 44889807, 44905777};
  transcriptFactorLUT["ZNF286A"] = {17, 15602890, 15624100};
  transcriptFactorLUT["ZNF287"] = {17, 16453630, 16472520};
  transcriptFactorLUT["ZNF292"] = {6, 87865268, 87973406};
  transcriptFactorLUT["ZNF296"] = {19, 45574757, 45579688};
  transcriptFactorLUT["ZNF3"] = {7, 99667593, 99680171};
  transcriptFactorLUT["ZNF30"] = {19, 35417806, 35436076};
  transcriptFactorLUT["ZNF300"] = {5, 150273953, 150284545};
  transcriptFactorLUT["ZNF302"] = {19, 35168543, 35177302};
  transcriptFactorLUT["ZNF304"] = {19, 57862641, 57871265};
  transcriptFactorLUT["ZNF311"] = {6, 28962593, 28973037};
  transcriptFactorLUT["ZNF317"] = {19, 9251055, 9274091};
  transcriptFactorLUT["ZNF319"] = {16, 58028571, 58033762};
  transcriptFactorLUT["ZNF32"] = {10, 44139306, 44144326};
  transcriptFactorLUT["ZNF320"] = {19, 53379424, 53394599};
  transcriptFactorLUT["ZNF324"] = {19, 58978411, 58984945};
  transcriptFactorLUT["ZNF324B"] = {19, 58962970, 58969199};
  transcriptFactorLUT["ZNF329"] = {19, 58637694, 58662148};
  transcriptFactorLUT["ZNF331"] = {19, 54041539, 54083523};
  transcriptFactorLUT["ZNF333"] = {19, 14800869, 14831772};
  transcriptFactorLUT["ZNF334"] = {20, 45128268, 45142198};
  transcriptFactorLUT["ZNF335"] = {20, 44577291, 44600833};
  transcriptFactorLUT["ZNF337"] = {20, 25653830, 25667575};
  transcriptFactorLUT["ZNF33A"] = {10, 38299577, 38348995};
  transcriptFactorLUT["ZNF33B"] = {10, 43084531, 43134285};
  transcriptFactorLUT["ZNF34"] = {8, 145997608, 146012730};
  transcriptFactorLUT["ZNF341"] = {20, 32319565, 32380075};
  transcriptFactorLUT["ZNF343"] = {20, 2462465, 2505168};
  transcriptFactorLUT["ZNF345"] = {19, 37341259, 37370477};
  transcriptFactorLUT["ZNF347"] = {19, 53641956, 53662322};
  transcriptFactorLUT["ZNF35"] = {3, 44690232, 44702283};
  transcriptFactorLUT["ZNF350"] = {19, 52467592, 52490079};
  transcriptFactorLUT["ZNF354A"] = {5, 178138521, 178157703};
  transcriptFactorLUT["ZNF354B"] = {5, 178286953, 178311424};
  transcriptFactorLUT["ZNF354C"] = {5, 178487606, 178507691};
  transcriptFactorLUT["ZNF358"] = {19, 7581003, 7585911};
  transcriptFactorLUT["ZNF362"] = {1, 33722173, 33766320};
  transcriptFactorLUT["ZNF366"] = {5, 71735725, 71803249};
  transcriptFactorLUT["ZNF367"] = {9, 99148224, 99180669};
  transcriptFactorLUT["ZNF37A"] = {10, 38383263, 38412278};
  transcriptFactorLUT["ZNF382"] = {19, 37096206, 37119499};
  transcriptFactorLUT["ZNF383"] = {19, 37717365, 37734574};
  transcriptFactorLUT["ZNF384"] = {12, 6775642, 6798541};
  transcriptFactorLUT["ZNF391"] = {6, 27356523, 27369227};
  transcriptFactorLUT["ZNF394"] = {7, 99090853, 99097877};
  transcriptFactorLUT["ZNF396"] = {18, 32946660, 32957301};
  transcriptFactorLUT["ZNF397"] = {18, 32820993, 32838397};
  transcriptFactorLUT["ZNF398"] = {7, 148844559, 148880118};
  transcriptFactorLUT["ZNF407"] = {18, 72342918, 72516583};
  transcriptFactorLUT["ZNF408"] = {11, 46722316, 46727466};
  transcriptFactorLUT["ZNF41"] = {23, 47305560, 47342345};
  transcriptFactorLUT["ZNF410"] = {14, 74353317, 74398991};
  transcriptFactorLUT["ZNF415"] = {19, 53611131, 53636171};
  transcriptFactorLUT["ZNF416"] = {19, 58082933, 58090243};
  transcriptFactorLUT["ZNF417"] = {19, 58417141, 58427978};
  transcriptFactorLUT["ZNF418"] = {19, 58433251, 58446740};
  transcriptFactorLUT["ZNF419"] = {19, 57999078, 58006048};
  transcriptFactorLUT["ZNF420"] = {19, 37569381, 37620651};
  transcriptFactorLUT["ZNF423"] = {16, 49524514, 49891830};
  transcriptFactorLUT["ZNF425"] = {7, 148799877, 148823438};
  transcriptFactorLUT["ZNF426"] = {19, 9638680, 9649303};
  transcriptFactorLUT["ZNF429"] = {19, 21688436, 21721079};
  transcriptFactorLUT["ZNF43"] = {19, 21987750, 22034870};
  transcriptFactorLUT["ZNF430"] = {19, 21203425, 21242852};
  transcriptFactorLUT["ZNF431"] = {19, 21324839, 21368805};
  transcriptFactorLUT["ZNF432"] = {19, 52536676, 52552073};
  transcriptFactorLUT["ZNF433"] = {19, 12125531, 12146525};
  transcriptFactorLUT["ZNF436"] = {1, 23685940, 23694879};
  transcriptFactorLUT["ZNF439"] = {19, 11976843, 11980306};
  transcriptFactorLUT["ZNF44"] = {19, 12382624, 12405714};
  transcriptFactorLUT["ZNF440"] = {19, 11925106, 11946016};
  transcriptFactorLUT["ZNF441"] = {19, 11877814, 11894893};
  transcriptFactorLUT["ZNF442"] = {19, 12460184, 12476475};
  transcriptFactorLUT["ZNF443"] = {19, 12540519, 12551926};
  transcriptFactorLUT["ZNF444"] = {19, 56652534, 56672262};
  transcriptFactorLUT["ZNF445"] = {3, 44481261, 44519162};
  transcriptFactorLUT["ZNF446"] = {19, 58987530, 58992601};
  transcriptFactorLUT["ZNF449"] = {23, 134478695, 134497338};
  transcriptFactorLUT["ZNF45"] = {19, 44416775, 44439411};
  transcriptFactorLUT["ZNF451"] = {6, 56954827, 57035098};
  transcriptFactorLUT["ZNF454"] = {5, 178368193, 178393218};
  transcriptFactorLUT["ZNF460"] = {19, 57791852, 57805436};
  transcriptFactorLUT["ZNF461"] = {19, 37128282, 37157755};
  transcriptFactorLUT["ZNF467"] = {7, 149461452, 149470295};
  transcriptFactorLUT["ZNF468"] = {19, 53341784, 53360902};
  transcriptFactorLUT["ZNF469"] = {16, 88493878, 88507165};
  transcriptFactorLUT["ZNF470"] = {19, 57078889, 57094262};
  transcriptFactorLUT["ZNF471"] = {19, 57019211, 57040269};
  transcriptFactorLUT["ZNF473"] = {19, 50529211, 50552031};
  transcriptFactorLUT["ZNF479"] = {7, 57187325, 57207571};
  transcriptFactorLUT["ZNF48"] = {16, 30406432, 30411429};
  transcriptFactorLUT["ZNF480"] = {19, 52800421, 52829180};
  transcriptFactorLUT["ZNF483"] = {9, 114287438, 114306712};
  transcriptFactorLUT["ZNF484"] = {9, 95607312, 95640320};
  transcriptFactorLUT["ZNF485"] = {10, 44101854, 44113352};
  transcriptFactorLUT["ZNF486"] = {19, 20278022, 20311299};
  transcriptFactorLUT["ZNF490"] = {19, 12686919, 12721623};
  transcriptFactorLUT["ZNF491"] = {19, 11909390, 11919306};
  transcriptFactorLUT["ZNF492"] = {19, 22817125, 22850472};
  transcriptFactorLUT["ZNF493"] = {19, 21579920, 21610296};
  transcriptFactorLUT["ZNF496"] = {1, 247463621, 247495045};
  transcriptFactorLUT["ZNF497"] = {19, 58865722, 58874214};
  transcriptFactorLUT["ZNF500"] = {16, 4800814, 4817219};
  transcriptFactorLUT["ZNF501"] = {3, 44771097, 44778575};
  transcriptFactorLUT["ZNF502"] = {3, 44754134, 44765323};
  transcriptFactorLUT["ZNF506"] = {19, 19903519, 19932560};
  transcriptFactorLUT["ZNF507"] = {19, 32836513, 32878573};
  transcriptFactorLUT["ZNF510"] = {9, 99518146, 99540328};
  transcriptFactorLUT["ZNF512"] = {2, 27805835, 27846082};
  transcriptFactorLUT["ZNF512B"] = {20, 62588056, 62601223};
  transcriptFactorLUT["ZNF513"] = {2, 27600097, 27603311};
  transcriptFactorLUT["ZNF514"] = {2, 95813399, 95825263};
  transcriptFactorLUT["ZNF516"] = {18, 74069636, 74207146};
  transcriptFactorLUT["ZNF517"] = {8, 146024260, 146034529};
  transcriptFactorLUT["ZNF519"] = {18, 14075988, 14132489};
  transcriptFactorLUT["ZNF521"] = {18, 22641887, 22932214};
  transcriptFactorLUT["ZNF524"] = {19, 56111705, 56114504};
  transcriptFactorLUT["ZNF526"] = {19, 42724491, 42732353};
  transcriptFactorLUT["ZNF527"] = {19, 37862058, 37883966};
  transcriptFactorLUT["ZNF528"] = {19, 52901120, 52921657};
  transcriptFactorLUT["ZNF529"] = {19, 37034516, 37096178};
  transcriptFactorLUT["ZNF530"] = {19, 58111252, 58119637};
  transcriptFactorLUT["ZNF532"] = {18, 56530060, 56653709};
  transcriptFactorLUT["ZNF534"] = {19, 52934666, 52955192};
  transcriptFactorLUT["ZNF536"] = {19, 30863327, 31048965};
  transcriptFactorLUT["ZNF540"] = {19, 38042272, 38105079};
  transcriptFactorLUT["ZNF541"] = {19, 48023941, 48059113};
  transcriptFactorLUT["ZNF543"] = {19, 57831864, 57842144};
  transcriptFactorLUT["ZNF544"] = {19, 58740069, 58775008};
  transcriptFactorLUT["ZNF546"] = {19, 40502942, 40526948};
  transcriptFactorLUT["ZNF547"] = {19, 57874802, 57890925};
  transcriptFactorLUT["ZNF548"] = {19, 57901217, 57913919};
  transcriptFactorLUT["ZNF549"] = {19, 58038692, 58052244};
  transcriptFactorLUT["ZNF550"] = {19, 58053203, 58071231};
  transcriptFactorLUT["ZNF551"] = {19, 58193336, 58201169};
  transcriptFactorLUT["ZNF552"] = {19, 58318449, 58326281};
  transcriptFactorLUT["ZNF554"] = {19, 2819871, 2836733};
  transcriptFactorLUT["ZNF555"] = {19, 2841432, 2860472};
  transcriptFactorLUT["ZNF556"] = {19, 2867332, 2878503};
  transcriptFactorLUT["ZNF557"] = {19, 7069470, 7087978};
  transcriptFactorLUT["ZNF558"] = {19, 8920251, 8942980};
  transcriptFactorLUT["ZNF559"] = {19, 9434447, 9454521};
  transcriptFactorLUT["ZNF560"] = {19, 9577030, 9609279};
  transcriptFactorLUT["ZNF561"] = {19, 9718001, 9731916};
  transcriptFactorLUT["ZNF562"] = {19, 9759337, 9785776};
  transcriptFactorLUT["ZNF563"] = {19, 12428303, 12444534};
  transcriptFactorLUT["ZNF564"] = {19, 12636183, 12662356};
  transcriptFactorLUT["ZNF565"] = {19, 36672961, 36705575};
  transcriptFactorLUT["ZNF566"] = {19, 36936020, 36980463};
  transcriptFactorLUT["ZNF567"] = {19, 37180301, 37212238};
  transcriptFactorLUT["ZNF568"] = {19, 37407230, 37488834};
  transcriptFactorLUT["ZNF569"] = {19, 37902059, 37958339};
  transcriptFactorLUT["ZNF57"] = {19, 2900895, 2918474};
  transcriptFactorLUT["ZNF570"] = {19, 37959959, 37976260};
  transcriptFactorLUT["ZNF571"] = {19, 38055154, 38085693};
  transcriptFactorLUT["ZNF572"] = {8, 125985538, 125991630};
  transcriptFactorLUT["ZNF573"] = {19, 38229202, 38270230};
  transcriptFactorLUT["ZNF574"] = {19, 42580289, 42585720};
  transcriptFactorLUT["ZNF575"] = {19, 44037339, 44040284};
  transcriptFactorLUT["ZNF576"] = {19, 44100543, 44104587};
  transcriptFactorLUT["ZNF577"] = {19, 52374550, 52391229};
  transcriptFactorLUT["ZNF579"] = {19, 56088890, 56092211};
  transcriptFactorLUT["ZNF580"] = {19, 56153462, 56154836};
  transcriptFactorLUT["ZNF581"] = {19, 56154985, 56156989};
  transcriptFactorLUT["ZNF582"] = {19, 56894647, 56904889};
  transcriptFactorLUT["ZNF583"] = {19, 56915717, 56936400};
  transcriptFactorLUT["ZNF584"] = {19, 58920062, 58929692};
  transcriptFactorLUT["ZNF585A"] = {19, 37638339, 37663338};
  transcriptFactorLUT["ZNF585B"] = {19, 37672480, 37701451};
  transcriptFactorLUT["ZNF586"] = {19, 58281019, 58291984};
  transcriptFactorLUT["ZNF587"] = {19, 58361180, 58376491};
  transcriptFactorLUT["ZNF589"] = {3, 48282595, 48312479};
  transcriptFactorLUT["ZNF594"] = {17, 5082830, 5095178};
  transcriptFactorLUT["ZNF596"] = {8, 182396, 197340};
  transcriptFactorLUT["ZNF597"] = {16, 3482421, 3493537};
  transcriptFactorLUT["ZNF599"] = {19, 35248978, 35264134};
  transcriptFactorLUT["ZNF600"] = {19, 53268747, 53290034};
  transcriptFactorLUT["ZNF606"] = {19, 58488440, 58514714};
  transcriptFactorLUT["ZNF607"] = {19, 38187263, 38210089};
  transcriptFactorLUT["ZNF610"] = {19, 52839497, 52870376};
  transcriptFactorLUT["ZNF611"] = {19, 53206065, 53233135};
  transcriptFactorLUT["ZNF613"] = {19, 52430687, 52449011};
  transcriptFactorLUT["ZNF614"] = {19, 52516576, 52531680};
  transcriptFactorLUT["ZNF615"] = {19, 52494586, 52511483};
  transcriptFactorLUT["ZNF616"] = {19, 52617652, 52643191};
  transcriptFactorLUT["ZNF619"] = {3, 40518603, 40531728};
  transcriptFactorLUT["ZNF620"] = {3, 40547601, 40559712};
  transcriptFactorLUT["ZNF621"] = {3, 40566368, 40581285};
  transcriptFactorLUT["ZNF623"] = {8, 144718182, 144735900};
  transcriptFactorLUT["ZNF624"] = {17, 16524047, 16557167};
  transcriptFactorLUT["ZNF625"] = {19, 12255708, 12267546};
  transcriptFactorLUT["ZNF626"] = {19, 20827508, 20844402};
  transcriptFactorLUT["ZNF627"] = {19, 11708234, 11729974};
  transcriptFactorLUT["ZNF628"] = {19, 55987698, 55995854};
  transcriptFactorLUT["ZNF629"] = {16, 30789769, 30798523};
  transcriptFactorLUT["ZNF630"] = {23, 47917566, 47930508};
  transcriptFactorLUT["ZNF639"] = {3, 179040778, 179053323};
  transcriptFactorLUT["ZNF641"] = {12, 48733792, 48744674};
  transcriptFactorLUT["ZNF644"] = {1, 91380856, 91487812};
  transcriptFactorLUT["ZNF646"] = {16, 31085742, 31094833};
  transcriptFactorLUT["ZNF648"] = {1, 182023704, 182030847};
  transcriptFactorLUT["ZNF649"] = {19, 52392487, 52408305};
  transcriptFactorLUT["ZNF652"] = {17, 47366567, 47439476};
  transcriptFactorLUT["ZNF653"] = {19, 11594241, 11616738};
  transcriptFactorLUT["ZNF655"] = {7, 99156447, 99174076};
  transcriptFactorLUT["ZNF658"] = {9, 40771401, 40792112};
  transcriptFactorLUT["ZNF660"] = {3, 44626455, 44637557};
  transcriptFactorLUT["ZNF662"] = {3, 42947657, 42960825};
  transcriptFactorLUT["ZNF664"] = {12, 124457761, 124499986};
  transcriptFactorLUT["ZNF665"] = {19, 53666551, 53696619};
  transcriptFactorLUT["ZNF667"] = {19, 56950692, 56988770};
  transcriptFactorLUT["ZNF668"] = {16, 31072163, 31076409};
  transcriptFactorLUT["ZNF669"] = {1, 247263263, 247267674};
  transcriptFactorLUT["ZNF670"] = {1, 247197939, 247242115};
  transcriptFactorLUT["ZNF671"] = {19, 58231118, 58238995};
  transcriptFactorLUT["ZNF672"] = {1, 249132376, 249143716};
  transcriptFactorLUT["ZNF674"] = {23, 46357159, 46404892};
  transcriptFactorLUT["ZNF675"] = {19, 23835707, 23870017};
  transcriptFactorLUT["ZNF676"] = {19, 22361902, 22379753};
  transcriptFactorLUT["ZNF677"] = {19, 53738637, 53758111};
  transcriptFactorLUT["ZNF678"] = {1, 227751219, 227850164};
  transcriptFactorLUT["ZNF679"] = {7, 63688851, 63727309};
  transcriptFactorLUT["ZNF680"] = {7, 63985033, 64023505};
  transcriptFactorLUT["ZNF681"] = {19, 23921996, 23941693};
  transcriptFactorLUT["ZNF682"] = {19, 20115226, 20150277};
  transcriptFactorLUT["ZNF683"] = {1, 26688124, 26699266};
  transcriptFactorLUT["ZNF684"] = {1, 40997232, 41013841};
  transcriptFactorLUT["ZNF687"] = {1, 151254787, 151264381};
  transcriptFactorLUT["ZNF688"] = {16, 30581018, 30583728};
  transcriptFactorLUT["ZNF689"] = {16, 30613878, 30621754};
  transcriptFactorLUT["ZNF69"] = {19, 11998669, 12025365};
  transcriptFactorLUT["ZNF691"] = {1, 43312243, 43318146};
  transcriptFactorLUT["ZNF692"] = {1, 249144202, 249153315};
  transcriptFactorLUT["ZNF696"] = {8, 144373558, 144382120};
  transcriptFactorLUT["ZNF697"] = {1, 120161999, 120190390};
  transcriptFactorLUT["ZNF699"] = {19, 9405985, 9415795};
  transcriptFactorLUT["ZNF7"] = {8, 146052902, 146068607};
  transcriptFactorLUT["ZNF70"] = {22, 24083770, 24093279};
  transcriptFactorLUT["ZNF700"] = {19, 12035882, 12061588};
  transcriptFactorLUT["ZNF701"] = {19, 53073525, 53090427};
  transcriptFactorLUT["ZNF705A"] = {12, 8325149, 8332642};
  transcriptFactorLUT["ZNF705D"] = {8, 11946846, 11973025};
  transcriptFactorLUT["ZNF705G"] = {8, 7215497, 7220490};
  transcriptFactorLUT["ZNF707"] = {8, 144766621, 144777555};
  transcriptFactorLUT["ZNF708"] = {19, 21473962, 21512212};
  transcriptFactorLUT["ZNF709"] = {19, 12571997, 12595632};
  transcriptFactorLUT["ZNF71"] = {19, 57106663, 57135544};
  transcriptFactorLUT["ZNF710"] = {15, 90544722, 90625432};
  transcriptFactorLUT["ZNF711"] = {23, 84498996, 84528368};
  transcriptFactorLUT["ZNF713"] = {7, 55954969, 56009918};
  transcriptFactorLUT["ZNF714"] = {19, 21264952, 21307883};
  transcriptFactorLUT["ZNF716"] = {7, 57509882, 57533265};
  transcriptFactorLUT["ZNF717"] = {3, 75786028, 75834734};
  transcriptFactorLUT["ZNF732"] = {4, 264463, 289944};
  transcriptFactorLUT["ZNF736"] = {7, 63774250, 63817012};
  transcriptFactorLUT["ZNF737"] = {19, 20720797, 20748626};
  transcriptFactorLUT["ZNF74"] = {22, 20748404, 20762753};
  transcriptFactorLUT["ZNF740"] = {12, 53574534, 53584654};
  transcriptFactorLUT["ZNF746"] = {7, 149169883, 149194898};
  transcriptFactorLUT["ZNF749"] = {19, 57946692, 57957191};
  transcriptFactorLUT["ZNF75A"] = {16, 3355405, 3368576};
  transcriptFactorLUT["ZNF75D"] = {23, 134419718, 134478012};
  transcriptFactorLUT["ZNF76"] = {6, 35227490, 35263764};
  transcriptFactorLUT["ZNF763"] = {19, 12075868, 12091198};
  transcriptFactorLUT["ZNF764"] = {16, 30565084, 30569642};
  transcriptFactorLUT["ZNF765"] = {19, 53898396, 53915262};
  transcriptFactorLUT["ZNF766"] = {19, 52772823, 52795976};
  transcriptFactorLUT["ZNF768"] = {16, 30535321, 30537910};
  transcriptFactorLUT["ZNF77"] = {19, 2933215, 2944969};
  transcriptFactorLUT["ZNF770"] = {15, 35270541, 35280497};
  transcriptFactorLUT["ZNF771"] = {16, 30418734, 30429916};
  transcriptFactorLUT["ZNF772"] = {19, 57980953, 57988938};
  transcriptFactorLUT["ZNF773"] = {19, 58011221, 58024519};
  transcriptFactorLUT["ZNF774"] = {15, 90895476, 90904715};
  transcriptFactorLUT["ZNF775"] = {7, 150076405, 150095719};
  transcriptFactorLUT["ZNF776"] = {19, 58258163, 58269527};
  transcriptFactorLUT["ZNF777"] = {7, 149128453, 149158053};
  transcriptFactorLUT["ZNF778"] = {16, 89284110, 89295965};
  transcriptFactorLUT["ZNF780A"] = {19, 40578898, 40596845};
  transcriptFactorLUT["ZNF780B"] = {19, 40534166, 40562115};
  transcriptFactorLUT["ZNF781"] = {19, 38158649, 38183216};
  transcriptFactorLUT["ZNF782"] = {9, 99579272, 99616389};
  transcriptFactorLUT["ZNF783"] = {7, 148959261, 148982085};
  transcriptFactorLUT["ZNF784"] = {19, 56132106, 56135941};
  transcriptFactorLUT["ZNF785"] = {16, 30591993, 30597092};
  transcriptFactorLUT["ZNF786"] = {7, 148766732, 148787869};
  transcriptFactorLUT["ZNF787"] = {19, 56598728, 56632742};
  transcriptFactorLUT["ZNF789"] = {7, 99070514, 99079948};
  transcriptFactorLUT["ZNF79"] = {9, 130186652, 130207651};
  transcriptFactorLUT["ZNF790"] = {19, 37309223, 37341689};
  transcriptFactorLUT["ZNF791"] = {19, 12721731, 12740676};
  transcriptFactorLUT["ZNF792"] = {19, 35447257, 35454953};
  transcriptFactorLUT["ZNF793"] = {19, 37997840, 38034239};
  transcriptFactorLUT["ZNF799"] = {19, 12500827, 12512088};
  transcriptFactorLUT["ZNF8"] = {19, 58790317, 58807254};
  transcriptFactorLUT["ZNF80"] = {3, 113953479, 113956425};
  transcriptFactorLUT["ZNF800"] = {7, 127010096, 127032778};
  transcriptFactorLUT["ZNF805"] = {19, 57752052, 57774106};
  transcriptFactorLUT["ZNF808"] = {19, 53030908, 53059303};
  transcriptFactorLUT["ZNF81"] = {23, 47696300, 47781655};
  transcriptFactorLUT["ZNF813"] = {19, 53970988, 53997546};
  transcriptFactorLUT["ZNF814"] = {19, 58380746, 58400442};
  transcriptFactorLUT["ZNF816"] = {19, 53452631, 53466164};
  transcriptFactorLUT["ZNF821"] = {16, 71893582, 71917444};
  transcriptFactorLUT["ZNF823"] = {19, 11832079, 11849824};
  transcriptFactorLUT["ZNF827"] = {4, 146681887, 146859607};
  transcriptFactorLUT["ZNF83"] = {19, 53115617, 53141644};
  transcriptFactorLUT["ZNF831"] = {20, 57766074, 57834167};
  transcriptFactorLUT["ZNF836"] = {19, 52658124, 52674896};
  transcriptFactorLUT["ZNF837"] = {19, 58878989, 58892389};
  transcriptFactorLUT["ZNF84"] = {12, 133614167, 133639890};
  transcriptFactorLUT["ZNF841"] = {19, 52567718, 52599018};
  transcriptFactorLUT["ZNF844"] = {19, 12175545, 12188626};
  transcriptFactorLUT["ZNF845"] = {19, 53837001, 53858122};
  transcriptFactorLUT["ZNF846"] = {19, 9868150, 9879410};
  transcriptFactorLUT["ZNF85"] = {19, 21106058, 21133503};
  transcriptFactorLUT["ZNF852"] = {3, 44540461, 44552132};
  transcriptFactorLUT["ZNF860"] = {3, 32023265, 32033228};
  transcriptFactorLUT["ZNF865"] = {19, 56124958, 56129907};
  transcriptFactorLUT["ZNF878"] = {19, 12154619, 12163782};
  transcriptFactorLUT["ZNF879"] = {5, 178450775, 178461388};
  transcriptFactorLUT["ZNF880"] = {19, 52873169, 52889046};
  transcriptFactorLUT["ZNF90"] = {19, 20188802, 20231977};
  transcriptFactorLUT["ZNF91"] = {19, 23540497, 23578362};
  transcriptFactorLUT["ZNF92"] = {7, 64838711, 64866048};
  transcriptFactorLUT["ZNF93"] = {19, 20011721, 20046382};
  transcriptFactorLUT["ZNF98"] = {19, 22573898, 22605148};
  transcriptFactorLUT["ZNF99"] = {19, 22934984, 22966973};
  transcriptFactorLUT["ZSCAN1"] = {19, 58545433, 58565999};
  transcriptFactorLUT["ZSCAN10"] = {16, 3138890, 3149318};
  transcriptFactorLUT["ZSCAN12"] = {6, 28356726, 28367544};
  transcriptFactorLUT["ZSCAN16"] = {6, 28092386, 28097856};
  transcriptFactorLUT["ZSCAN18"] = {19, 58595208, 58609730};
  transcriptFactorLUT["ZSCAN2"] = {15, 85144248, 85159841};
  transcriptFactorLUT["ZSCAN20"] = {1, 33938231, 33961995};
  transcriptFactorLUT["ZSCAN21"] = {7, 99647416, 99662663};
  transcriptFactorLUT["ZSCAN22"] = {19, 58838384, 58853712};
  transcriptFactorLUT["ZSCAN23"] = {6, 28400431, 28411279};
  transcriptFactorLUT["ZSCAN29"] = {15, 43650369, 43662258};
  transcriptFactorLUT["ZSCAN30"] = {18, 32831021, 32870209};
  transcriptFactorLUT["ZSCAN4"] = {19, 58180302, 58190520};
  transcriptFactorLUT["ZSCAN5A"] = {19, 56732678, 56739659};
  transcriptFactorLUT["ZSCAN5B"] = {19, 56701057, 56704421};
  transcriptFactorLUT["ZXDA"] = {23, 57931863, 57937067};
  transcriptFactorLUT["ZXDB"] = {23, 57618268, 57623910};
  transcriptFactorLUT["ZXDC"] = {3, 126177743, 126194762};
  transcriptFactorLUT["ZZZ3"] = {1, 78030189, 78148343};
  return true;
}

