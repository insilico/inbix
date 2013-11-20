/* 
 * File:   EpistasisEQtl.cpp
 * Author: bwhite
 * 
 * Created on October 3, 2013, 11:48 AM
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
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
  
  // for each transcript build main effect and epistasis regression models
  PP->printLOG("epiQTL linear regression loop for all transcripts\n");
  if(localCis) {
    PP->printLOG("epiQTL local cis mode with radius: " + 
      int2str(int(radius / 1000)) + " kilobases\n");
  }
  if(par::epiqtl_interaction_full) {
    PP->printLOG("epiQTL FULL epistatic interaction mode\n");
  } else {
    PP->printLOG("epiQTL cis-cis/cis-trans interaction mode\n");
  }
  
  int transcriptIndex = 0;
  string testnumbersFilename = par::output_file_name + ".testnumbers.txt";
  PP->printLOG("Writing test results to [ " + testnumbersFilename + " ]\n");
  std::ofstream TESTNUMBERS;
  TESTNUMBERS.open(testnumbersFilename.c_str(), ios::out);
  string thisTranscript;
  for(; transcriptIndex < PP->nlistname.size(); ++transcriptIndex) {
    
    thisTranscript = PP->nlistname[transcriptIndex];
    PP->printLOG("Transcript: " + thisTranscript + "\n");
    //cout << "Transcript: " << thisTranscript << endl;

    // get transcript expression vector as phenotype
    PP->setQtlPhenoFromNumericIndex(transcriptIndex);
    
    // EQTL -------------------------------------------------------------------
    // fit main effect regression model for SNPs
    //cout << "Running main effects regression models" << endl;
    string eqtlFilename = par::output_file_name + "." + 
      thisTranscript + ".eqtl.txt";
    PP->printLOG("Writing eQTL results to [ " + eqtlFilename + " ]\n");
    std::ofstream EQTL;
    EQTL.open(eqtlFilename.c_str(), ios::out);
    int thisSnpIndex = 0;
    for(thisSnpIndex=0; thisSnpIndex < PP->nl_all; ++thisSnpIndex) {
      string thisSnpName = PP->locus[thisSnpIndex]->name;
      
      Model* mainEffectModel = new LinearModel(PP);
      mainEffectModel->setMissing();
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
      // Build design matrix
      mainEffectModel->buildDesignMatrix();

      // Fit linear model
      int tp = 1; 
      mainEffectModel->testParameter = tp; // single variable main effect
      mainEffectModel->fitLM();

      // Obtain estimates and statistics
      vector_t betaMainEffectCoefs = mainEffectModel->getCoefs();
      // p-values don't include intercept term
      vector_t betaMainEffectCoefPvals = mainEffectModel->getPVals();
      double mainEffectPval = betaMainEffectCoefPvals[tp - 1];

      // always use first coefficient after intercept as main effect term
      double mainEffectValue = betaMainEffectCoefs[tp];
      double mainEffectPValue = betaMainEffectCoefPvals[tp-1];
      EQTL 
        << thisSnpName << "\t"
        << thisTranscript << "\t"
        << mainEffectValue << "\t"
        << mainEffectPValue << endl;

      delete mainEffectModel;
    }
    EQTL.close();
    
    // EPIQTL -----------------------------------------------------------------
    int nAllSnps = PP->nl_all;
    int nInnerLoop = -1;
    vector<int> thisTranscriptSnpIndices;
    if(par::epiqtl_interaction_full) {
      nInnerLoop = nAllSnps;
    } else {
      GetSnpsForTranscript(thisTranscript, thisTranscriptSnpIndices);
      nInnerLoop = thisTranscriptSnpIndices.size();
    }
    // allocate results matrices
    double** resultsMatrixBetas= new double*[nAllSnps];
    for(int a=0; a < nAllSnps; ++a) {
      resultsMatrixBetas[a] = new double[nInnerLoop];
      for(int b=0; b < nInnerLoop; ++b) {
        resultsMatrixBetas[a][b] = 0.0;
      }
    }
    double** resultsMatrixPvals= new double*[nAllSnps];
    for(int a=0; a < nAllSnps; ++a) {
      resultsMatrixPvals[a] = new double[nInnerLoop];
      for(int b=0; b < nInnerLoop; ++b) {
        resultsMatrixPvals[a][b] = 0.0;
      }
    }
    
    if(par::epiqtl_interaction_full) {
      PP->printLOG("epiQTL linear regression loop: SNP x SNP\n");
#pragma omp parallel for schedule(dynamic, 10) 
      for(int ii=0; ii < nAllSnps; ++ii) {
        if(ii && (ii % 1000 == 0)) {
          cout << ii << "/" << nAllSnps << endl;
        }
        for(int jj=ii+1; jj < nAllSnps; ++jj) {
          int snpAIndex = ii;
          string snpAName = PP->locus[snpAIndex]->name;
          int snpBIndex = jj;
          string snpBName = PP->locus[snpBIndex]->name;
          Model* interactionModel = new LinearModel(PP);
          interactionModel->setMissing();
          interactionModel->addAdditiveSNP(snpAIndex);
          interactionModel->label.push_back(snpAName);
          interactionModel->addAdditiveSNP(snpBIndex);
          interactionModel->label.push_back(snpBName);
          if(par::covar_file) {
            for(int kk = 0; kk < par::clist_number; kk++) {
              interactionModel->addCovariate(kk);
              interactionModel->label.push_back(PP->clistname[kk]);
            }
          }
          interactionModel->addInteraction(1, 2);
          interactionModel->label.push_back("EPI");
          interactionModel->buildDesignMatrix();
          int tp = 3;
          // add # covars to test param to get interaction param
          if(par::covar_file) {
            tp += par::clist_number;
          }
          interactionModel->testParameter = tp; // interaction
          interactionModel->fitLM();
          vector_t betaInteractionCoefs = interactionModel->getCoefs();
          double interactionValue = 
            betaInteractionCoefs[betaInteractionCoefs.size() - 1];
          vector_t betaInteractionCoefPVals = interactionModel->getPVals();
          double interactionPval =
            betaInteractionCoefPVals[betaInteractionCoefPVals.size() - 1];
          resultsMatrixBetas[ii][jj] = interactionValue;
          resultsMatrixPvals[ii][jj] = interactionPval;
//          cout 
//            << (interactionModel->fitConverged() ? "TRUE" : "FALSE")
//            << "\t" << snpAIndex << "\t" << snpAName
//            << "\t" << snpBIndex << "\t" << snpBName 
//            << "\t" << interactionValue << "\t" << interactionPval
//            << "\t" << betaInteractionCoefs[3]
//            << "\t" << betaInteractionCoefPVals[2]
//            << endl;
//          exit(1);
          delete interactionModel;
        }
      }
    } else {
      PP->printLOG("epiQTL linear regression loop: SNP x cis/trans (" + 
        int2str(nInnerLoop) + ")\n");
#pragma omp parallel for schedule(dynamic, 10) 
      for(int ii=0; ii < nAllSnps; ++ii) {
        if(ii && (ii % 1000 == 0)) {
          cout << ii << "/" << nAllSnps << endl;
        }
        for(int jj=0; jj < nInnerLoop; ++jj) {
          int snpAIndex = ii;
          string snpAName = PP->locus[snpAIndex]->name;
          int snpBIndex = thisTranscriptSnpIndices[jj];
          string snpBName = PP->locus[snpBIndex]->name;
          Model* interactionModel = new LinearModel(PP);
          interactionModel->setMissing();
          interactionModel->addAdditiveSNP(snpAIndex);
          interactionModel->label.push_back(snpAName);
          interactionModel->addAdditiveSNP(snpBIndex);
          interactionModel->label.push_back(snpBName);
          if(par::covar_file) {
            for(int kk = 0; kk < par::clist_number; kk++) {
              interactionModel->addCovariate(kk);
              interactionModel->label.push_back(PP->clistname[kk]);
            }
          }
          interactionModel->addInteraction(1, 2);
          interactionModel->label.push_back("EPI");
          interactionModel->buildDesignMatrix();
          int tp = 3;
          // add # covars to test param to get interaction param
          if(par::covar_file) {
            tp += par::clist_number;
          }
          interactionModel->testParameter = tp; // interaction
          interactionModel->fitLM();
          vector_t betaInteractionCoefs = interactionModel->getCoefs();
          double interactionValue = 
            betaInteractionCoefs[betaInteractionCoefs.size() - 1];
          vector_t betaInteractionCoefPVals = interactionModel->getPVals();
          double interactionPval =
            betaInteractionCoefPVals[betaInteractionCoefPVals.size() - 1];
          resultsMatrixBetas[ii][jj] = interactionValue;
          resultsMatrixPvals[ii][jj] = interactionPval;
          delete interactionModel;
        }
      }
    }

    // write regression results
    string epiqtlFilename = par::output_file_name + "." + 
      thisTranscript + ".epiqtl.txt";
    PP->printLOG("Writing epiQTL results to [ " + epiqtlFilename + " ]\n");
    ofstream EPIQTL_OUT;
    EPIQTL_OUT.open(epiqtlFilename.c_str(), ios::out);
    if(par::epiqtl_interaction_full) {
      for(int kk=0; kk < nAllSnps; ++kk) {
        for(int ll=kk+1; ll < nAllSnps; ++ll) {
          int snpAIndex = kk;
          string snpAName = PP->locus[snpAIndex]->name;
          int snpBIndex = ll;
          string snpBName = PP->locus[snpBIndex]->name;
          EPIQTL_OUT
            << snpAName << "\t" << snpBName << "\t"
            << thisTranscript << "\t"
            << resultsMatrixBetas[kk][ll] << "\t"
            << resultsMatrixPvals[kk][ll] << endl;
        }
      }
    }
    else {
      for(int kk=0; kk < nAllSnps; ++kk) {
        for(int ll=0; ll < nInnerLoop; ++ll) {
          int snpAIndex = kk;
          string snpAName = PP->locus[snpAIndex]->name;
          int snpBIndex = thisTranscriptSnpIndices[ll];
          string snpBName = PP->locus[snpBIndex]->name;
          EPIQTL_OUT
            << snpAName << "\t" << snpBName << "\t"
            << thisTranscript << "\t"
            << resultsMatrixBetas[kk][ll] << "\t"
            << resultsMatrixPvals[kk][ll] << endl;
        }
      }
    }
    EPIQTL_OUT.close();

    // release dynamically allocated memory
    for(int a=0; a < nAllSnps; ++a) {
      delete [] resultsMatrixBetas[a];
    }
    delete [] resultsMatrixBetas;
    for(int a=0; a < nAllSnps; ++a) {
      delete [] resultsMatrixPvals[a];
    }
    delete [] resultsMatrixPvals;

    
    TESTNUMBERS << thisTranscript << "\t" 
      << (nAllSnps * nInnerLoop) << endl;
    
  } // END for each transcript loop
  
  TESTNUMBERS.close();
  
  PP->printLOG("epiQTL analysis finished\n");

  return true;
}

bool EpistasisEQtl::GetSnpsForTranscript(string transcript, 
  vector<int>& snpIndices) {

  // get transcript info
  int chromosome = coordinates[transcript][COORD_CHROM];
  int bpStart = coordinates[transcript][COORD_BP_START];
  int bpEnd = coordinates[transcript][COORD_BP_END];
  int lowerThreshold = bpStart - radius;
  int upperThreshold = bpEnd + radius;
  
//  cout 
//    << chromosome  << ", " 
//    << radius  << ", " 
//    << "(" << lowerThreshold << "), "
//    << bpStart  << ", " 
//    << bpEnd << ", "
//    << "(" << upperThreshold << ")"
//    << endl;
  
  // find SNPs matching criteria
  for(int j=0; j < PP->locus.size(); ++j) {
    Locus* thisSnp = PP->locus[j];
    if(thisSnp->chr == chromosome) {
      if(localCis) {
        // on the same chromosome and within radius of transcript
        if(thisSnp->bp >= lowerThreshold && 
          thisSnp->bp <= upperThreshold) {
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
