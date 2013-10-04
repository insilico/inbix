/* 
 * File:   EpistasisEQtl.cpp
 * Author: bwhite
 * 
 * Created on October 3, 2013, 11:48 AM
 */

#include <iostream>
#include <string>
#include <vector>

#include "plink.h"
#include "model.h"
#include "linear.h"
#include "stats.h"
#include "helper.h"

#include "EpistasisEQtl.h"

using namespace std;

// Plink object
extern Plink* PP;

EpistasisEQtl::EpistasisEQtl() {
  radius = -1;
  localCis = false;
}

EpistasisEQtl::~EpistasisEQtl() {
}

bool EpistasisEQtl::ReadTranscriptCoordinates(string coordinatesFilename) {
  // open the numeric attributes file if possible
  checkFileExists(coordinatesFilename);
  ifstream coordinatesFile(coordinatesFilename.c_str(), ios::in);
  if(coordinatesFile.fail()) {
    return false;
  }

  int rows = 0;
  while(!coordinatesFile.eof()) {
    char nline[par::MAX_LINE_LENGTH];
    coordinatesFile.getline(nline, par::MAX_LINE_LENGTH, '\n');

    // convert to string
    string sline = nline;
    if(sline == "") continue;

    // read line from text file into a vector of tokens
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
    for(int i=0; i < tokens.size()-1; ++i) {
      int t = 0;
      if(!from_string<int>(t, tokens[i], std::dec)) {
        cerr << "Error parsing transcript info to integer on line " << rows << endl;
        cerr << "token: " << tokens[i] << endl;
        return false;
      }
      coordinates[gene].push_back(t);
    }
  }
  coordinatesFile.close();

  PP->printLOG("Read " + int2str(rows) + 
    " transcript coordinates info from ["  + coordinatesFilename + "]\n");
  
  return true;
}

bool EpistasisEQtl::SetRadius(int newRadius) {
  if(newRadius < 1) {
    return false;
  }
  radius = newRadius;
  return true;
}

bool EpistasisEQtl::SetLocalCis(bool localCisFlag) {
  localCis = localCisFlag;
  return true;
}

bool EpistasisEQtl::Run() {
  // we have SNPs?
  int numSnps = PP->nl_all;
  if(!numSnps) {
    cerr << "Error: no SNPs found " << endl;
    return false;
  }
  
  // we have transcript expression levels?
  int numTranscripts = PP->nlistname.size();
  if(!numTranscripts) {
    cerr << "Error: no transcript values found " << endl;
    return false;
  }
  
  // we have a transcript lookup table that matches the expression data?
  int numTranscriptInfo = coordinates.size();
  if(numTranscriptInfo != numTranscripts) {
    cerr << "Error: number of coordinate file entries does not match "
      "the number of transcript values found " << endl;
    return false;
  }
  vector<string>::const_iterator cit = PP->nlistname.begin();
  for(; cit != PP->nlistname.end(); ++cit) {
    // do all the transcripts in expression values have coordinates info?
    if(coordinates.find(*cit) == coordinates.end()) {
      cerr << "Error: Transcript " << *cit << " not found in coordinate file" << endl;
      return false;
    }
  }
  
  // open output files
  string testnumbersFilename = par::output_file_name + ".testnumbers.txt";
  string eqtlFilename = par::output_file_name + ".eqtl.txt";
  string epiqtlFilename = par::output_file_name + ".epiqtl.txt";

  PP->printLOG("Writing test results to [ " + testnumbersFilename + " ]\n");
  PP->printLOG("Writing eQTL results to [ " + eqtlFilename + " ]\n");
  PP->printLOG("Writing epiQTL results to [ " + epiqtlFilename + " ]\n");

  ofstream TESTNUMBERS;
  TESTNUMBERS.open(testnumbersFilename.c_str(), ios::out);
  ofstream EQTL;
  EQTL.open(eqtlFilename.c_str(), ios::out);
  ofstream EPIQTL;
  EPIQTL.open(epiqtlFilename.c_str(), ios::out);
  
  // for each transcript build main effect and epistasis regression models
  PP->printLOG("epiQTL linear regression loop\n");
	PP->SNP2Ind();
  for(int i=0; i < PP->nlistname.size(); ++i) {
    string thisTranscript = PP->nlistname[i];
    cout << "Transcript: " << thisTranscript << endl;
    
    // get transcript expression vector as phenotype
    PP->setQtlPhenoFromNumericIndex(i);
    
    // get transcript chromosome
    int thisChromosome = coordinates[thisTranscript][0];
    
    // get all SNP indices on the chromosome
    vector<int> thisTranscriptSnpIndices;
    for(int j=0; j < PP->locus.size(); ++j) {
      Locus* thisSnp = PP->locus[j];
      if(thisSnp->chr == thisChromosome) {
        thisTranscriptSnpIndices.push_back(j);
      }
    }
    cout << thisTranscriptSnpIndices.size() << " SNPs found for transcript" << endl;
    TESTNUMBERS << thisTranscript << "\t" << thisTranscriptSnpIndices.size() << endl;
    
    // fit main effect regression model for SNPs
    cout << "Running main effects regression models" << endl;
    pair<double, double> snpResult;
    int j;
#pragma omp parallel for
    for(j=0; j < thisTranscriptSnpIndices.size(); ++j) {
      int thisSnpIndex = thisTranscriptSnpIndices[j];
      string thisSnpName = PP->locus[thisSnpIndex]->name;
      
      Model* mainEffectModel = new LinearModel(PP);
      mainEffectModel->setMissing();
      pair<double, double> snpResult;
      mainEffectModel->addAdditiveSNP(thisSnpIndex);
      mainEffectModel->label.push_back(thisSnpName);
      // add covariates if specified
      if(par::covar_file) {
        for(int k=0; k < par::clist_number; k++) {
          // add covariate to the model
          mainEffectModel->addCovariate(k);
          mainEffectModel->label.push_back(PP->clistname[k]);
        }
      }
      
      // fit the model and write results
      snpResult = fitModel(mainEffectModel);
#pragma omp critical
{
      EQTL 
        << thisSnpName << "\t"
        << thisTranscript << "\t"
        << snpResult.first << "\t"
        << snpResult.second << endl;
}      
      delete mainEffectModel;
    }
    
    // interaction regression model for each pair of SNPs
    cout << "Running interaction regression models" << endl;
    int ii, jj;
#pragma omp parallel for schedule(dynamic, 1) private(ii, jj)
    for(ii=0; ii < thisTranscriptSnpIndices.size(); ++ii) {
      for(jj=ii+1; jj < thisTranscriptSnpIndices.size(); ++jj) {
        //cout << "MODEL" << endl;
        int snpAIndex = thisTranscriptSnpIndices[ii];
        string snpAName = PP->locus[snpAIndex]->name;
        int snpBIndex = thisTranscriptSnpIndices[jj];
        string snpBName = PP->locus[snpBIndex]->name;
        //cout << ii << "\t" << snpAIndex << "\t" << snpAName << endl;
        //cout << jj << "\t" << snpBIndex << "\t" << snpBName << endl;
        
        //cout << "Preparing interaction regression model" << endl;
        Model* interactionModel = new LinearModel(PP);
        interactionModel->setMissing();
        interactionModel->addAdditiveSNP(snpAIndex);
        interactionModel->label.push_back(snpAName);
        interactionModel->addAdditiveSNP(snpBIndex);
        interactionModel->label.push_back(snpBName);
        if(par::covar_file) {
          for(int kk = 0; kk < par::clist_number; kk++) {
            // add covariate to the model
            interactionModel->addCovariate(kk);
            interactionModel->label.push_back(PP->clistname[kk]);
          }
        }
        interactionModel->addInteraction(1, 2);
        interactionModel->label.push_back("EPI");
        
        // fit the model and write results
        //cout << "Fitting interaction regression model" << endl;
        interactionModel->buildDesignMatrix();
        interactionModel->fitLM();
        //cout << "Getting interaction regression model" << endl;
        vector_t betaInteractionCoefs = interactionModel->getCoefs();
        double interactionValue = 
          betaInteractionCoefs[betaInteractionCoefs.size() - 1];
        vector_t betaInteractionCoefPVals = interactionModel->getPVals();
        double interactionPval =
          betaInteractionCoefPVals[betaInteractionCoefPVals.size() - 1];
#pragma omp critical
{
        //cout << "Writing interaction regression model" << endl;
        EPIQTL 
          << snpAName << "\t" << snpBName << "\t"
          << thisTranscript << "\t"
          << interactionValue << "\t"
          << interactionPval << endl;
}
        delete interactionModel;
      }
    }
    
  } // for each transcript loop
  
  // clean up
  TESTNUMBERS.close();
  EQTL.close();
  EPIQTL.close();
  
  return true;
}
