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
  tfMode = false;
  tfRadius = 0;
  LoadDefaultTranscriptionFactorLUT();
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

  coordinates.clear();
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

bool EpistasisEQtl::ReadTranscriptFactorCoordinates(string coordinatesFilename) {
  // open the numeric attributes file if possible
  checkFileExists(coordinatesFilename);
  ifstream coordinatesFile(coordinatesFilename.c_str(), ios::in);
  if(coordinatesFile.fail()) {
    return false;
  }

  transcriptFactorLUT.clear();
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
    int bpStart;
    int bpEnd;
    if(!from_string<int>(bpStart, tokens[1], std::dec)) {
      cerr << "Error parsing transcription factor info to integer on line " << rows << endl;
      cerr << "token: " << tokens[1] << endl;
    }
    if(!from_string<int>(bpEnd, tokens[2], std::dec)) {
      cerr << "Error parsing transcription factor info to integer on line " << rows << endl;
      cerr << "token: " << tokens[2] << endl;
    }
    transcriptFactorLUT[gene] = make_pair(bpStart, bpEnd);
  }
  coordinatesFile.close();

  PP->printLOG("Read " + int2str(rows) + 
    " transcription factor coordinates info from ["  + coordinatesFilename + "]\n");
  
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


bool EpistasisEQtl::SetTFRadius(int newRadius) {
  if(newRadius < 1) {
    cerr << "Error setting TF radius to: " << newRadius << endl;
    return false;
  }
  // newRadius is in kilobases, but need to store a bases
  tfRadius = newRadius * 1000;
  return true;
}

bool EpistasisEQtl::SetTF(bool tfFlag) {
  tfMode = tfFlag;
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
  
  string testnumbersFilename = par::output_file_name + ".testnumbers.txt";
  PP->printLOG("Writing test results to [ " + testnumbersFilename + " ]\n");
  std::ofstream TESTNUMBERS;
  TESTNUMBERS.open(testnumbersFilename.c_str(), ios::out);
  string thisTranscript;
  int transcriptIndex = 0;
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
      double mainEffectValue = betaMainEffectCoefs[tp];
      // p-values don't include intercept term
      vector_t betaMainEffectCoefPvals = mainEffectModel->getPVals();
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
#pragma omp parallel for
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
#pragma omp critical
{          
          vector_t betaInteractionCoefs = interactionModel->getCoefs();
          double interactionValue = 
            betaInteractionCoefs[betaInteractionCoefs.size() - 1];
          vector_t betaInteractionCoefPVals = interactionModel->getPVals();
          double interactionPval =
            betaInteractionCoefPVals[betaInteractionCoefPVals.size() - 1];
          resultsMatrixBetas[ii][jj] = interactionValue;
          resultsMatrixPvals[ii][jj] = interactionPval;
}
          delete interactionModel;
        }
      }
    } else {
      PP->printLOG("epiQTL linear regression loop: SNP x cis/trans (" + 
        int2str(nInnerLoop) + ")\n");
#pragma omp parallel for
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
#pragma omp critical
{          
          vector_t betaInteractionCoefs = interactionModel->getCoefs();
          double interactionValue = 
            betaInteractionCoefs[betaInteractionCoefs.size() - 1];
          vector_t betaInteractionCoefPVals = interactionModel->getPVals();
          double interactionPval =
            betaInteractionCoefPVals[betaInteractionCoefPVals.size() - 1];
           // if((interactionPval < 0) || (interactionPval > 1)) {
           //   cout 
           //    << "!!!!! DANGER !!!!!\t"
           //     << (interactionModel->fitConverged() ? "TRUE" : "FALSE")
           //     << "\t" << snpAIndex << "\t" << snpAName
           //     << "\t" << snpBIndex << "\t" << snpBName 
           //     << "\t" << interactionValue << "\t" << interactionPval
           //     << "\t" << betaInteractionCoefs[3]
           //     << "\t" << betaInteractionCoefPVals[2]
           //     << endl;
           //    exit(1);
           //  }
          resultsMatrixBetas[ii][jj] = interactionValue;
          resultsMatrixPvals[ii][jj] = interactionPval;
}
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

bool EpistasisEQtl::LoadDefaultTranscriptionFactorLUT() {
  transcriptFactorLUT["ADNP"] = make_pair(49505454, 49547527);
  transcriptFactorLUT["AFF1"] = make_pair(87856153, 88062206);
  transcriptFactorLUT["AFF2"] = make_pair(147582138, 148082193);
  transcriptFactorLUT["AFF3"] = make_pair(100163715, 100759037);
  transcriptFactorLUT["AFF4"] = make_pair(132211070, 132299354);
  transcriptFactorLUT["AHR"] = make_pair(17338275, 17385775);
  transcriptFactorLUT["AHRR"] = make_pair(304290, 438405);
  transcriptFactorLUT["AIRE"] = make_pair(45705720, 45718102);
  transcriptFactorLUT["ALX1"] = make_pair(85674035, 85695561);
  transcriptFactorLUT["ALX3"] = make_pair(110602996, 110613322);
  transcriptFactorLUT["ALX4"] = make_pair(44282277, 44331716);
  transcriptFactorLUT["AR"] = make_pair(66763873, 66950461);
  transcriptFactorLUT["ARGFX"] = make_pair(121286777, 121309469);
  transcriptFactorLUT["ARID1A"] = make_pair(27022521, 27108601);
  transcriptFactorLUT["ARID1B"] = make_pair(157099063, 157531913);
  transcriptFactorLUT["ARID2"] = make_pair(46123619, 46301819);
  transcriptFactorLUT["ARID3A"] = make_pair(926036, 972803);
  transcriptFactorLUT["ARID3B"] = make_pair(74833547, 74890472);
  transcriptFactorLUT["ARID3C"] = make_pair(34621454, 34628011);
  transcriptFactorLUT["ARID4A"] = make_pair(58765221, 58840451);
  transcriptFactorLUT["ARID4B"] = make_pair(235330209, 235491532);
  transcriptFactorLUT["ARID5A"] = make_pair(97202463, 97218371);
  transcriptFactorLUT["ARID5B"] = make_pair(63661012, 63856707);
  transcriptFactorLUT["ARNT"] = make_pair(150782180, 150849244);
  transcriptFactorLUT["ARNT2"] = make_pair(80696691, 80890277);
  transcriptFactorLUT["ARNTL"] = make_pair(13299273, 13408812);
  transcriptFactorLUT["ARNTL2"] = make_pair(27485786, 27578746);
  transcriptFactorLUT["ARX"] = make_pair(25021812, 25034065);
  transcriptFactorLUT["ASCL1"] = make_pair(103351451, 103354294);
  transcriptFactorLUT["ASCL2"] = make_pair(2289727, 2292182);
  transcriptFactorLUT["ASCL3"] = make_pair(8959118, 8964580);
  transcriptFactorLUT["ASCL4"] = make_pair(108168161, 108170421);
  transcriptFactorLUT["ATF1"] = make_pair(51157788, 51214943);
  transcriptFactorLUT["ATF2"] = make_pair(175936977, 176032934);
  transcriptFactorLUT["ATF3"] = make_pair(212738675, 212794119);
  transcriptFactorLUT["ATF4"] = make_pair(39916568, 39918691);
  transcriptFactorLUT["ATF5"] = make_pair(50431958, 50437193);
  transcriptFactorLUT["ATF6"] = make_pair(161736033, 161933860);
  transcriptFactorLUT["ATF6B"] = make_pair(32083044, 32096017);
  transcriptFactorLUT["ATF7"] = make_pair(53905842, 54020199);
  transcriptFactorLUT["ATOH1"] = make_pair(94750077, 94751142);
  transcriptFactorLUT["ATOH7"] = make_pair(69990351, 69991870);
  transcriptFactorLUT["ATOH8"] = make_pair(85980908, 86018506);
  transcriptFactorLUT["BACH1"] = make_pair(30671115, 30718469);
  transcriptFactorLUT["BACH2"] = make_pair(90636246, 91006627);
  transcriptFactorLUT["BARHL1"] = make_pair(135457992, 135465640);
  transcriptFactorLUT["BARHL2"] = make_pair(91177578, 91182794);
  transcriptFactorLUT["BARX1"] = make_pair(96713908, 96717608);
  transcriptFactorLUT["BARX2"] = make_pair(129245880, 129322174);
  transcriptFactorLUT["BATF"] = make_pair(75988783, 76013334);
  transcriptFactorLUT["BATF2"] = make_pair(64755416, 64757749);
  transcriptFactorLUT["BATF3"] = make_pair(212859758, 212873327);
  transcriptFactorLUT["BAZ2A"] = make_pair(56989379, 57030163);
  transcriptFactorLUT["BAZ2B"] = make_pair(160175489, 160473112);
  transcriptFactorLUT["BBX"] = make_pair(107241782, 107530176);
  transcriptFactorLUT["BCL11A"] = make_pair(60684328, 60780633);
  transcriptFactorLUT["BCL11B"] = make_pair(99635624, 99738050);
  transcriptFactorLUT["BCL6"] = make_pair(187439164, 187454285);
  transcriptFactorLUT["BCL6B"] = make_pair(6926368, 6932961);
  transcriptFactorLUT["BHLHA15"] = make_pair(97841565, 97842271);
  transcriptFactorLUT["BHLHA9"] = make_pair(1173857, 1174565);
  transcriptFactorLUT["BHLHE22"] = make_pair(65492794, 65496191);
  transcriptFactorLUT["BHLHE23"] = make_pair(61637330, 61638387);
  transcriptFactorLUT["BHLHE40"] = make_pair(5021096, 5026865);
  transcriptFactorLUT["BHLHE41"] = make_pair(26272958, 26278003);
  transcriptFactorLUT["BNC2"] = make_pair(16409500, 16870786);
  transcriptFactorLUT["BSX"] = make_pair(122848356, 122852379);
  transcriptFactorLUT["C11orf95"] = make_pair(63527363, 63536113);
  transcriptFactorLUT["CAMTA1"] = make_pair(6845383, 7829766);
  transcriptFactorLUT["CAMTA2"] = make_pair(4871286, 4890960);
  transcriptFactorLUT["CARHSP1"] = make_pair(8946798, 8962869);
  transcriptFactorLUT["CBFB"] = make_pair(67063049, 67134958);
  transcriptFactorLUT["CCDC79"] = make_pair(66788878, 66835523);
  transcriptFactorLUT["CDC5L"] = make_pair(44355250, 44418161);
  transcriptFactorLUT["CDX1"] = make_pair(149546343, 149564121);
  transcriptFactorLUT["CDX2"] = make_pair(28536204, 28543505);
  transcriptFactorLUT["CDX4"] = make_pair(72667089, 72674421);
  transcriptFactorLUT["CEBPA"] = make_pair(33790839, 33793470);
  transcriptFactorLUT["CEBPB"] = make_pair(48807119, 48809227);
  transcriptFactorLUT["CEBPD"] = make_pair(48649475, 48650726);
  transcriptFactorLUT["CEBPE"] = make_pair(23586514, 23588820);
  transcriptFactorLUT["CEBPG"] = make_pair(33865429, 33873592);
  transcriptFactorLUT["CIC"] = make_pair(42788816, 42799948);
  transcriptFactorLUT["CLOCK"] = make_pair(56294067, 56413076);
  transcriptFactorLUT["CREB1"] = make_pair(208394615, 208470284);
  transcriptFactorLUT["CREB3"] = make_pair(35732316, 35737005);
  transcriptFactorLUT["CREB3L1"] = make_pair(46299188, 46342972);
  transcriptFactorLUT["CREB3L2"] = make_pair(137597556, 137686847);
  transcriptFactorLUT["CREB3L3"] = make_pair(4153597, 4173051);
  transcriptFactorLUT["CREB3L4"] = make_pair(153940314, 153946840);
  transcriptFactorLUT["CREB5"] = make_pair(28475233, 28865511);
  transcriptFactorLUT["CREBL2"] = make_pair(12764766, 12798042);
  transcriptFactorLUT["CREM"] = make_pair(35456466, 35501886);
  transcriptFactorLUT["CRX"] = make_pair(48325098, 48346586);
  transcriptFactorLUT["CSDC2"] = make_pair(41957013, 41972670);
  transcriptFactorLUT["CSDE1"] = make_pair(115259533, 115300671);
  transcriptFactorLUT["CTCF"] = make_pair(67596309, 67673088);
  transcriptFactorLUT["CTCFL"] = make_pair(56081797, 56100163);
  transcriptFactorLUT["CUX1"] = make_pair(101459183, 101927250);
  transcriptFactorLUT["CUX2"] = make_pair(111471827, 111788358);
  transcriptFactorLUT["DBP"] = make_pair(49133816, 49140807);
  transcriptFactorLUT["DBX1"] = make_pair(20177759, 20181870);
  transcriptFactorLUT["DBX2"] = make_pair(45408538, 45444882);
  transcriptFactorLUT["DDIT3"] = make_pair(57910370, 57914300);
  transcriptFactorLUT["DEAF1"] = make_pair(644219, 695754);
  transcriptFactorLUT["DLX1"] = make_pair(172950207, 172954401);
  transcriptFactorLUT["DLX2"] = make_pair(172964165, 172967478);
  transcriptFactorLUT["DLX3"] = make_pair(48067368, 48072588);
  transcriptFactorLUT["DLX4"] = make_pair(48046561, 48052323);
  transcriptFactorLUT["DLX5"] = make_pair(96649701, 96654143);
  transcriptFactorLUT["DLX6"] = make_pair(96635289, 96640352);
  transcriptFactorLUT["DMBX1"] = make_pair(46972667, 46979886);
  transcriptFactorLUT["DMRT1"] = make_pair(841689, 969090);
  transcriptFactorLUT["DMRT2"] = make_pair(1050353, 1057554);
  transcriptFactorLUT["DMRT3"] = make_pair(976967, 991732);
  transcriptFactorLUT["DMRTA1"] = make_pair(22446839, 22452472);
  transcriptFactorLUT["DMRTA2"] = make_pair(50883222, 50889119);
  transcriptFactorLUT["DMRTB1"] = make_pair(53925071, 53933160);
  transcriptFactorLUT["DMRTC2"] = make_pair(42349085, 42356397);
  transcriptFactorLUT["DMTF1"] = make_pair(86781676, 86825648);
  transcriptFactorLUT["DNAJC1"] = make_pair(22045476, 22292650);
  transcriptFactorLUT["DNAJC2"] = make_pair(102952920, 102985320);
  transcriptFactorLUT["DPRX"] = make_pair(54135309, 54140263);
  transcriptFactorLUT["DRGX"] = make_pair(50574160, 50604062);
  transcriptFactorLUT["DUXA"] = make_pair(57663093, 57678856);
  transcriptFactorLUT["E2F1"] = make_pair(32263291, 32274210);
  transcriptFactorLUT["E2F2"] = make_pair(23832919, 23857712);
  transcriptFactorLUT["E2F3"] = make_pair(20403909, 20493945);
  transcriptFactorLUT["E2F4"] = make_pair(67226067, 67232821);
  transcriptFactorLUT["E2F5"] = make_pair(86099909, 86126753);
  transcriptFactorLUT["E2F6"] = make_pair(11584500, 11606303);
  transcriptFactorLUT["E2F7"] = make_pair(77415025, 77459360);
  transcriptFactorLUT["E2F8"] = make_pair(19245609, 19262507);
  transcriptFactorLUT["E4F1"] = make_pair(2273488, 2285743);
  transcriptFactorLUT["EAF2"] = make_pair(121554033, 121605373);
  transcriptFactorLUT["EBF1"] = make_pair(158122922, 158526788);
  transcriptFactorLUT["EBF2"] = make_pair(25699245, 25902640);
  transcriptFactorLUT["EBF3"] = make_pair(131633495, 131762091);
  transcriptFactorLUT["EBF4"] = make_pair(2673523, 2740754);
  transcriptFactorLUT["EGR1"] = make_pair(137801180, 137805004);
  transcriptFactorLUT["EGR2"] = make_pair(64571755, 64576126);
  transcriptFactorLUT["EGR3"] = make_pair(22545173, 22550815);
  transcriptFactorLUT["EGR4"] = make_pair(73518056, 73520829);
  transcriptFactorLUT["EHF"] = make_pair(34642587, 34684834);
  transcriptFactorLUT["ELF1"] = make_pair(41506054, 41556418);
  transcriptFactorLUT["ELF2"] = make_pair(139978870, 140005568);
  transcriptFactorLUT["ELF3"] = make_pair(201979689, 201986315);
  transcriptFactorLUT["ELF4"] = make_pair(129198894, 129244688);
  transcriptFactorLUT["ELF5"] = make_pair(34500341, 34535347);
  transcriptFactorLUT["ELK1"] = make_pair(47494918, 47510003);
  transcriptFactorLUT["ELK3"] = make_pair(96588159, 96663613);
  transcriptFactorLUT["ELK4"] = make_pair(205588395, 205602000);
  transcriptFactorLUT["EMX1"] = make_pair(73144603, 73162020);
  transcriptFactorLUT["EMX2"] = make_pair(119301955, 119309057);
  transcriptFactorLUT["EN1"] = make_pair(119599746, 119605759);
  transcriptFactorLUT["EN2"] = make_pair(155250823, 155257526);
  transcriptFactorLUT["EOMES"] = make_pair(27757439, 27764206);
  transcriptFactorLUT["ERF"] = make_pair(42751712, 42759316);
  transcriptFactorLUT["ERG"] = make_pair(39751949, 40033704);
  transcriptFactorLUT["ESR1"] = make_pair(152011630, 152424408);
  transcriptFactorLUT["ESR2"] = make_pair(64699746, 64761128);
  transcriptFactorLUT["ESRRA"] = make_pair(64072999, 64084212);
  transcriptFactorLUT["ESRRB"] = make_pair(76837689, 76968180);
  transcriptFactorLUT["ESRRG"] = make_pair(216676587, 217113015);
  transcriptFactorLUT["ESX1"] = make_pair(103494718, 103499599);
  transcriptFactorLUT["ETS1"] = make_pair(128328655, 128457453);
  transcriptFactorLUT["ETS2"] = make_pair(40177754, 40196878);
  transcriptFactorLUT["ETV1"] = make_pair(13930855, 14026139);
  transcriptFactorLUT["ETV2"] = make_pair(36132638, 36135773);
  transcriptFactorLUT["ETV3"] = make_pair(157102975, 157108177);
  transcriptFactorLUT["ETV4"] = make_pair(41605210, 41623800);
  transcriptFactorLUT["ETV5"] = make_pair(185764105, 185826901);
  transcriptFactorLUT["ETV6"] = make_pair(11802787, 12048325);
  transcriptFactorLUT["ETV7"] = make_pair(36333970, 36355577);
  transcriptFactorLUT["EVX1"] = make_pair(27282163, 27287438);
  transcriptFactorLUT["EVX2"] = make_pair(176944834, 176948690);
  transcriptFactorLUT["FERD3L"] = make_pair(19184404, 19185044);
  transcriptFactorLUT["FEV"] = make_pair(219845808, 219850379);
  transcriptFactorLUT["FEZF1"] = make_pair(121941362, 121944565);
  transcriptFactorLUT["FEZF2"] = make_pair(62355346, 62359190);
  transcriptFactorLUT["FIGLA"] = make_pair(71004441, 71017775);
  transcriptFactorLUT["FIZ1"] = make_pair(56102736, 56110893);
  transcriptFactorLUT["FLI1"] = make_pair(128556429, 128683162);
  transcriptFactorLUT["FOS"] = make_pair(75745480, 75748937);
  transcriptFactorLUT["FOSB"] = make_pair(45971252, 45978437);
  transcriptFactorLUT["FOSL1"] = make_pair(65659606, 65667997);
  transcriptFactorLUT["FOSL2"] = make_pair(28615778, 28637516);
  transcriptFactorLUT["FOXA1"] = make_pair(38058756, 38064325);
  transcriptFactorLUT["FOXA2"] = make_pair(22561641, 22565101);
  transcriptFactorLUT["FOXA3"] = make_pair(46367517, 46377055);
  transcriptFactorLUT["FOXB1"] = make_pair(60296420, 60298142);
  transcriptFactorLUT["FOXB2"] = make_pair(79634570, 79635869);
  transcriptFactorLUT["FOXC1"] = make_pair(1610680, 1614129);
  transcriptFactorLUT["FOXC2"] = make_pair(86600856, 86602537);
  transcriptFactorLUT["FOXD2"] = make_pair(47901688, 47906363);
  transcriptFactorLUT["FOXD3"] = make_pair(63788729, 63790797);
  transcriptFactorLUT["FOXD4"] = make_pair(116230, 118417);
  transcriptFactorLUT["FOXD4L1"] = make_pair(114256660, 114258727);
  transcriptFactorLUT["FOXD4L3"] = make_pair(70917782, 70920000);
  transcriptFactorLUT["FOXD4L4"] = make_pair(42718065, 42719316);
  transcriptFactorLUT["FOXD4L5"] = make_pair(70175706, 70178815);
  transcriptFactorLUT["FOXD4L6"] = make_pair(69199479, 69202204);
  transcriptFactorLUT["FOXE1"] = make_pair(100615536, 100618997);
  transcriptFactorLUT["FOXE3"] = make_pair(47881743, 47883724);
  transcriptFactorLUT["FOXF1"] = make_pair(86544132, 86548070);
  transcriptFactorLUT["FOXF2"] = make_pair(1390068, 1395832);
  transcriptFactorLUT["FOXG1"] = make_pair(29236277, 29239483);
  transcriptFactorLUT["FOXH1"] = make_pair(145699114, 145701718);
  transcriptFactorLUT["FOXI1"] = make_pair(169532916, 169536729);
  transcriptFactorLUT["FOXI2"] = make_pair(129535537, 129539450);
  transcriptFactorLUT["FOXJ1"] = make_pair(74132414, 74137380);
  transcriptFactorLUT["FOXJ2"] = make_pair(8185358, 8208118);
  transcriptFactorLUT["FOXJ3"] = make_pair(42642209, 42800903);
  transcriptFactorLUT["FOXK1"] = make_pair(4721929, 4811074);
  transcriptFactorLUT["FOXK2"] = make_pair(80477593, 80562483);
  transcriptFactorLUT["FOXL1"] = make_pair(86612114, 86615304);
  transcriptFactorLUT["FOXL2"] = make_pair(138663065, 138665982);
  transcriptFactorLUT["FOXM1"] = make_pair(2966846, 2986321);
  transcriptFactorLUT["FOXN1"] = make_pair(26850958, 26865175);
  transcriptFactorLUT["FOXN2"] = make_pair(48541794, 48606434);
  transcriptFactorLUT["FOXN3"] = make_pair(89622515, 89883454);
  transcriptFactorLUT["FOXN4"] = make_pair(109715782, 109747025);
  transcriptFactorLUT["FOXO1"] = make_pair(41129800, 41240734);
  transcriptFactorLUT["FOXO3"] = make_pair(108881025, 109005971);
  transcriptFactorLUT["FOXO4"] = make_pair(70315998, 70323384);
  transcriptFactorLUT["FOXO6"] = make_pair(41827602, 41849263);
  transcriptFactorLUT["FOXP1"] = make_pair(71247033, 71633140);
  transcriptFactorLUT["FOXP2"] = make_pair(114055051, 114333827);
  transcriptFactorLUT["FOXP3"] = make_pair(49106896, 49121288);
  transcriptFactorLUT["FOXP4"] = make_pair(41514163, 41570122);
  transcriptFactorLUT["FOXQ1"] = make_pair(1312674, 1314993);
  transcriptFactorLUT["FOXR1"] = make_pair(118842416, 118851995);
  transcriptFactorLUT["FOXR2"] = make_pair(55649832, 55652621);
  transcriptFactorLUT["FOXS1"] = make_pair(30432102, 30433420);
  transcriptFactorLUT["GABPA"] = make_pair(27107257, 27144771);
  transcriptFactorLUT["GATA1"] = make_pair(48644981, 48652717);
  transcriptFactorLUT["GATA2"] = make_pair(128198264, 128212030);
  transcriptFactorLUT["GATA3"] = make_pair(8096666, 8117164);
  transcriptFactorLUT["GATA4"] = make_pair(11561716, 11617509);
  transcriptFactorLUT["GATA5"] = make_pair(61038552, 61051026);
  transcriptFactorLUT["GATA6"] = make_pair(19749397, 19782491);
  transcriptFactorLUT["GATAD1"] = make_pair(92076761, 92089381);
  transcriptFactorLUT["GATAD2B"] = make_pair(153777202, 153895451);
  transcriptFactorLUT["GBX1"] = make_pair(150845675, 150864635);
  transcriptFactorLUT["GBX2"] = make_pair(237073878, 237076652);
  transcriptFactorLUT["GCM1"] = make_pair(52991759, 53013624);
  transcriptFactorLUT["GCM2"] = make_pair(10873455, 10882098);
  transcriptFactorLUT["GFI1"] = make_pair(92940317, 92951628);
  transcriptFactorLUT["GFI1B"] = make_pair(135854097, 135867084);
  transcriptFactorLUT["GLI1"] = make_pair(57853917, 57866047);
  transcriptFactorLUT["GLI2"] = make_pair(121554866, 121750229);
  transcriptFactorLUT["GLI3"] = make_pair(42000547, 42276618);
  transcriptFactorLUT["GLI4"] = make_pair(144349606, 144359101);
  transcriptFactorLUT["GLIS1"] = make_pair(53971905, 54199877);
  transcriptFactorLUT["GLIS2"] = make_pair(4382224, 4389598);
  transcriptFactorLUT["GLIS3"] = make_pair(3824127, 4300035);
  transcriptFactorLUT["GMEB1"] = make_pair(28995239, 29042115);
  transcriptFactorLUT["GMEB2"] = make_pair(62218954, 62258381);
  transcriptFactorLUT["GRHL1"] = make_pair(10091791, 10142412);
  transcriptFactorLUT["GRHL2"] = make_pair(102504667, 102681952);
  transcriptFactorLUT["GRHL3"] = make_pair(24649529, 24681808);
  transcriptFactorLUT["GSC"] = make_pair(95234559, 95236499);
  transcriptFactorLUT["GSC2"] = make_pair(19136503, 19137796);
  transcriptFactorLUT["GSX1"] = make_pair(28366779, 28368089);
  transcriptFactorLUT["GSX2"] = make_pair(54966247, 54968122);
  transcriptFactorLUT["GTF2I"] = make_pair(74071990, 74175022);
  transcriptFactorLUT["GTF2IRD1"] = make_pair(73868119, 74016920);
  transcriptFactorLUT["GTF2IRD2"] = make_pair(74247755, 74267872);
  transcriptFactorLUT["GTF2IRD2B"] = make_pair(74508346, 74565623);
  transcriptFactorLUT["GTF3A"] = make_pair(27998680, 28009846);
  transcriptFactorLUT["GZF1"] = make_pair(23345020, 23353683);
  transcriptFactorLUT["HAND1"] = make_pair(153854531, 153857824);
  transcriptFactorLUT["HAND2"] = make_pair(174447651, 174451378);
  transcriptFactorLUT["HBP1"] = make_pair(106809405, 106842974);
  transcriptFactorLUT["HDX"] = make_pair(83572881, 83757487);
  transcriptFactorLUT["HELT"] = make_pair(185939994, 185941958);
  transcriptFactorLUT["HES1"] = make_pair(193853930, 193856401);
  transcriptFactorLUT["HES2"] = make_pair(6475293, 6479979);
  transcriptFactorLUT["HES3"] = make_pair(6304251, 6305638);
  transcriptFactorLUT["HES4"] = make_pair(934343, 935552);
  transcriptFactorLUT["HES5"] = make_pair(2460183, 2461684);
  transcriptFactorLUT["HES6"] = make_pair(239146907, 239148765);
  transcriptFactorLUT["HES7"] = make_pair(8023907, 8027410);
  transcriptFactorLUT["HESX1"] = make_pair(57231943, 57234280);
  transcriptFactorLUT["HEY1"] = make_pair(80676244, 80680098);
  transcriptFactorLUT["HEY2"] = make_pair(126070731, 126082415);
  transcriptFactorLUT["HEYL"] = make_pair(40089102, 40105348);
  transcriptFactorLUT["HHEX"] = make_pair(94449680, 94455408);
  transcriptFactorLUT["HIC1"] = make_pair(1958392, 1962981);
  transcriptFactorLUT["HIC2"] = make_pair(21771692, 21805750);
  transcriptFactorLUT["HIF1A"] = make_pair(62162118, 62214977);
  transcriptFactorLUT["HINFP"] = make_pair(118992232, 119005765);
  transcriptFactorLUT["HIVEP1"] = make_pair(12012723, 12165232);
  transcriptFactorLUT["HIVEP2"] = make_pair(143072603, 143266338);
  transcriptFactorLUT["HIVEP3"] = make_pair(42312859, 42501596);
  transcriptFactorLUT["HKR1"] = make_pair(37825579, 37855357);
  transcriptFactorLUT["HLF"] = make_pair(53342320, 53402426);
  transcriptFactorLUT["HLX"] = make_pair(221052742, 221058400);
  transcriptFactorLUT["HMBOX1"] = make_pair(28747910, 28910242);
  transcriptFactorLUT["HMG20A"] = make_pair(77712992, 77777946);
  transcriptFactorLUT["HMG20B"] = make_pair(3572942, 3579081);
  transcriptFactorLUT["HMGA1"] = make_pair(34204576, 34214008);
  transcriptFactorLUT["HMGA2"] = make_pair(66218239, 66346311);
  transcriptFactorLUT["HMGB1"] = make_pair(31032878, 31040081);
  transcriptFactorLUT["HMGB2"] = make_pair(174252526, 174255595);
  transcriptFactorLUT["HMGB3"] = make_pair(150151747, 150159248);
  transcriptFactorLUT["HMGXB3"] = make_pair(149380168, 149432706);
  transcriptFactorLUT["HMGXB4"] = make_pair(35653444, 35691800);
  transcriptFactorLUT["HMX1"] = make_pair(8847801, 8873543);
  transcriptFactorLUT["HMX2"] = make_pair(124907637, 124910188);
  transcriptFactorLUT["HMX3"] = make_pair(124895566, 124897247);
  transcriptFactorLUT["HNF1A"] = make_pair(121416548, 121440314);
  transcriptFactorLUT["HNF1B"] = make_pair(36046433, 36105069);
  transcriptFactorLUT["HNF4A"] = make_pair(42984440, 43061485);
  transcriptFactorLUT["HNF4G"] = make_pair(76452202, 76479061);
  transcriptFactorLUT["HOMEZ"] = make_pair(23742843, 23755309);
  transcriptFactorLUT["HOPX"] = make_pair(57514153, 57547872);
  transcriptFactorLUT["HOXA1"] = make_pair(27132613, 27135625);
  transcriptFactorLUT["HOXA10"] = make_pair(27210209, 27213955);
  transcriptFactorLUT["HOXA11"] = make_pair(27220775, 27224835);
  transcriptFactorLUT["HOXA13"] = make_pair(27236498, 27239725);
  transcriptFactorLUT["HOXA2"] = make_pair(27139972, 27142394);
  transcriptFactorLUT["HOXA3"] = make_pair(27145808, 27166639);
  transcriptFactorLUT["HOXA4"] = make_pair(27168125, 27170399);
  transcriptFactorLUT["HOXA5"] = make_pair(27180670, 27183287);
  transcriptFactorLUT["HOXA6"] = make_pair(27185201, 27187393);
  transcriptFactorLUT["HOXA7"] = make_pair(27193337, 27196296);
  transcriptFactorLUT["HOXA9"] = make_pair(27202056, 27205149);
  transcriptFactorLUT["HOXB1"] = make_pair(46606806, 46608272);
  transcriptFactorLUT["HOXB13"] = make_pair(46802126, 46806111);
  transcriptFactorLUT["HOXB2"] = make_pair(46620018, 46622393);
  transcriptFactorLUT["HOXB3"] = make_pair(46626231, 46651810);
  transcriptFactorLUT["HOXB4"] = make_pair(46652868, 46655743);
  transcriptFactorLUT["HOXB5"] = make_pair(46668618, 46671103);
  transcriptFactorLUT["HOXB6"] = make_pair(46673098, 46682334);
  transcriptFactorLUT["HOXB7"] = make_pair(46684594, 46688383);
  transcriptFactorLUT["HOXB8"] = make_pair(46689707, 46692301);
  transcriptFactorLUT["HOXB9"] = make_pair(46698518, 46703835);
  transcriptFactorLUT["HOXC10"] = make_pair(54378945, 54384062);
  transcriptFactorLUT["HOXC11"] = make_pair(54366909, 54370203);
  transcriptFactorLUT["HOXC12"] = make_pair(54348713, 54350350);
  transcriptFactorLUT["HOXC13"] = make_pair(54332575, 54340328);
  transcriptFactorLUT["HOXC4"] = make_pair(54410635, 54449814);
  transcriptFactorLUT["HOXC5"] = make_pair(54426831, 54429144);
  transcriptFactorLUT["HOXC6"] = make_pair(54422193, 54424607);
  transcriptFactorLUT["HOXC8"] = make_pair(54402889, 54406545);
  transcriptFactorLUT["HOXC9"] = make_pair(54393876, 54397120);
  transcriptFactorLUT["HOXD1"] = make_pair(177053306, 177055635);
  transcriptFactorLUT["HOXD10"] = make_pair(176981491, 176984670);
  transcriptFactorLUT["HOXD11"] = make_pair(176972083, 176974316);
  transcriptFactorLUT["HOXD12"] = make_pair(176964529, 176965488);
  transcriptFactorLUT["HOXD13"] = make_pair(176957531, 176960666);
  transcriptFactorLUT["HOXD3"] = make_pair(177028804, 177037826);
  transcriptFactorLUT["HOXD4"] = make_pair(177016112, 177017949);
  transcriptFactorLUT["HOXD8"] = make_pair(176994467, 176997423);
  transcriptFactorLUT["HOXD9"] = make_pair(176987412, 176989645);
  transcriptFactorLUT["HSF1"] = make_pair(145515269, 145538385);
  transcriptFactorLUT["HSF4"] = make_pair(67197287, 67203848);
  transcriptFactorLUT["HSFX1"] = make_pair(148674182, 148676974);
  transcriptFactorLUT["HSFX2"] = make_pair(148674171, 148676970);
  transcriptFactorLUT["HSFY1"] = make_pair(20708576, 20750849);
  transcriptFactorLUT["HSFY2"] = make_pair(20708556, 20750849);
  transcriptFactorLUT["ID1"] = make_pair(30193085, 30194317);
  transcriptFactorLUT["ID2"] = make_pair(8822112, 8824583);
  transcriptFactorLUT["ID3"] = make_pair(23884420, 23886285);
  transcriptFactorLUT["ID4"] = make_pair(19837600, 19842431);
  transcriptFactorLUT["IKZF1"] = make_pair(50344264, 50472798);
  transcriptFactorLUT["IKZF2"] = make_pair(213864410, 214016333);
  transcriptFactorLUT["IKZF3"] = make_pair(37913967, 38020441);
  transcriptFactorLUT["IKZF4"] = make_pair(56414688, 56432219);
  transcriptFactorLUT["IKZF5"] = make_pair(124750321, 124768366);
  transcriptFactorLUT["INSM1"] = make_pair(20348764, 20351592);
  transcriptFactorLUT["INSM2"] = make_pair(36003247, 36006260);
  transcriptFactorLUT["IRF1"] = make_pair(131817300, 131826465);
  transcriptFactorLUT["IRF2"] = make_pair(185308875, 185395726);
  transcriptFactorLUT["IRF3"] = make_pair(50162825, 50169132);
  transcriptFactorLUT["IRF4"] = make_pair(391738, 411443);
  transcriptFactorLUT["IRF5"] = make_pair(128578270, 128590096);
  transcriptFactorLUT["IRF6"] = make_pair(209958967, 209979520);
  transcriptFactorLUT["IRF7"] = make_pair(612554, 615999);
  transcriptFactorLUT["IRF8"] = make_pair(85932773, 85956211);
  transcriptFactorLUT["IRF9"] = make_pair(24630421, 24635774);
  transcriptFactorLUT["IRX1"] = make_pair(3596167, 3601517);
  transcriptFactorLUT["IRX2"] = make_pair(2746278, 2751769);
  transcriptFactorLUT["IRX3"] = make_pair(54317211, 54320378);
  transcriptFactorLUT["IRX4"] = make_pair(1877540, 1887098);
  transcriptFactorLUT["IRX5"] = make_pair(54965110, 54968395);
  transcriptFactorLUT["IRX6"] = make_pair(55358470, 55364672);
  transcriptFactorLUT["ISL1"] = make_pair(50678957, 50690563);
  transcriptFactorLUT["ISL2"] = make_pair(76629064, 76634816);
  transcriptFactorLUT["ISX"] = make_pair(35462128, 35483380);
  transcriptFactorLUT["JARID2"] = make_pair(15246205, 15522273);
  transcriptFactorLUT["JDP2"] = make_pair(75898836, 75939404);
  transcriptFactorLUT["JUN"] = make_pair(59246462, 59249785);
  transcriptFactorLUT["JUNB"] = make_pair(12902309, 12904125);
  transcriptFactorLUT["JUND"] = make_pair(18390503, 18392466);
  transcriptFactorLUT["KDM5A"] = make_pair(389222, 498621);
  transcriptFactorLUT["KDM5B"] = make_pair(202696531, 202777549);
  transcriptFactorLUT["KDM5C"] = make_pair(53220502, 53254604);
  transcriptFactorLUT["KDM5D"] = make_pair(21867300, 21906825);
  transcriptFactorLUT["KIAA2018"] = make_pair(113367232, 113415493);
  transcriptFactorLUT["KLF1"] = make_pair(12995236, 12998017);
  transcriptFactorLUT["KLF10"] = make_pair(103661004, 103666192);
  transcriptFactorLUT["KLF11"] = make_pair(10184371, 10194963);
  transcriptFactorLUT["KLF12"] = make_pair(74260148, 74708066);
  transcriptFactorLUT["KLF13"] = make_pair(31619057, 31727868);
  transcriptFactorLUT["KLF14"] = make_pair(130417381, 130418860);
  transcriptFactorLUT["KLF15"] = make_pair(126061477, 126076236);
  transcriptFactorLUT["KLF16"] = make_pair(1852397, 1863564);
  transcriptFactorLUT["KLF17"] = make_pair(44584521, 44600809);
  transcriptFactorLUT["KLF2"] = make_pair(16435650, 16438339);
  transcriptFactorLUT["KLF3"] = make_pair(38665789, 38703129);
  transcriptFactorLUT["KLF4"] = make_pair(110247132, 110252047);
  transcriptFactorLUT["KLF5"] = make_pair(73629113, 73651680);
  transcriptFactorLUT["KLF6"] = make_pair(3818187, 3827473);
  transcriptFactorLUT["KLF7"] = make_pair(207938861, 208031970);
  transcriptFactorLUT["KLF8"] = make_pair(56258869, 56314322);
  transcriptFactorLUT["KLF9"] = make_pair(72999512, 73029573);
  transcriptFactorLUT["L3MBTL1"] = make_pair(42143075, 42170535);
  transcriptFactorLUT["L3MBTL4"] = make_pair(5954704, 6414910);
  transcriptFactorLUT["LBX1"] = make_pair(102986732, 102988717);
  transcriptFactorLUT["LBX2"] = make_pair(74724643, 74726685);
  transcriptFactorLUT["LCOR"] = make_pair(98592016, 98724198);
  transcriptFactorLUT["LCORL"] = make_pair(17844838, 18023483);
  transcriptFactorLUT["LEF1"] = make_pair(108968700, 109087953);
  transcriptFactorLUT["LEUTX"] = make_pair(40267233, 40276775);
  transcriptFactorLUT["LHX1"] = make_pair(35294771, 35301915);
  transcriptFactorLUT["LHX2"] = make_pair(126773888, 126795442);
  transcriptFactorLUT["LHX3"] = make_pair(139088095, 139095004);
  transcriptFactorLUT["LHX4"] = make_pair(180199432, 180244188);
  transcriptFactorLUT["LHX5"] = make_pair(113900693, 113909877);
  transcriptFactorLUT["LHX6"] = make_pair(124964855, 124991091);
  transcriptFactorLUT["LHX8"] = make_pair(75600566, 75627218);
  transcriptFactorLUT["LHX9"] = make_pair(197881634, 197899273);
  transcriptFactorLUT["LIN28A"] = make_pair(26737268, 26756219);
  transcriptFactorLUT["LIN28B"] = make_pair(105404922, 105531207);
  transcriptFactorLUT["LITAF"] = make_pair(11641577, 11680806);
  transcriptFactorLUT["LMX1A"] = make_pair(165171103, 165325952);
  transcriptFactorLUT["LMX1B"] = make_pair(129376721, 129463311);
  transcriptFactorLUT["LYL1"] = make_pair(13209841, 13213974);
  transcriptFactorLUT["MAF"] = make_pair(79627744, 79634622);
  transcriptFactorLUT["MAFA"] = make_pair(144510229, 144512602);
  transcriptFactorLUT["MAFB"] = make_pair(39314487, 39317880);
  transcriptFactorLUT["MAFF"] = make_pair(38597938, 38612517);
  transcriptFactorLUT["MAFG"] = make_pair(79876144, 79881444);
  transcriptFactorLUT["MAFK"] = make_pair(1570367, 1582679);
  transcriptFactorLUT["MAX"] = make_pair(65472818, 65569262);
  transcriptFactorLUT["MAZ"] = make_pair(29817857, 29822504);
  transcriptFactorLUT["MBD1"] = make_pair(47793251, 47808144);
  transcriptFactorLUT["MBD2"] = make_pair(51729049, 51751158);
  transcriptFactorLUT["MBD3"] = make_pair(1576669, 1592760);
  transcriptFactorLUT["MBD4"] = make_pair(129149786, 129159022);
  transcriptFactorLUT["MECOM"] = make_pair(168801286, 168865522);
  transcriptFactorLUT["MECP2"] = make_pair(153295685, 153363188);
  transcriptFactorLUT["MEF2A"] = make_pair(100106132, 100256629);
  transcriptFactorLUT["MEF2B"] = make_pair(19256375, 19281098);
  transcriptFactorLUT["MEF2D"] = make_pair(156433512, 156460391);
  transcriptFactorLUT["MEIS1"] = make_pair(66662531, 66799891);
  transcriptFactorLUT["MEIS2"] = make_pair(37183221, 37392341);
  transcriptFactorLUT["MEIS3"] = make_pair(47906374, 47922785);
  transcriptFactorLUT["MEOX1"] = make_pair(41717757, 41738931);
  transcriptFactorLUT["MEOX2"] = make_pair(15650836, 15726308);
  transcriptFactorLUT["MESP1"] = make_pair(90293097, 90294540);
  transcriptFactorLUT["MESP2"] = make_pair(90319588, 90321982);
  transcriptFactorLUT["MGA"] = make_pair(41952609, 42062141);
  transcriptFactorLUT["MIER3"] = make_pair(56215428, 56247957);
  transcriptFactorLUT["MITF"] = make_pair(69985750, 70017488);
  transcriptFactorLUT["MIXL1"] = make_pair(226411318, 226414756);
  transcriptFactorLUT["MLX"] = make_pair(40719077, 40725221);
  transcriptFactorLUT["MLXIP"] = make_pair(122516633, 122631894);
  transcriptFactorLUT["MLXIPL"] = make_pair(73007523, 73038870);
  transcriptFactorLUT["MNT"] = make_pair(2287353, 2304258);
  transcriptFactorLUT["MNX1"] = make_pair(156797546, 156803347);
  transcriptFactorLUT["MSC"] = make_pair(72753776, 72756731);
  transcriptFactorLUT["MSGN1"] = make_pair(17997785, 17998367);
  transcriptFactorLUT["MSX1"] = make_pair(4861391, 4865660);
  transcriptFactorLUT["MSX2"] = make_pair(174151574, 174157902);
  transcriptFactorLUT["MTA1"] = make_pair(105886185, 105937057);
  transcriptFactorLUT["MTA2"] = make_pair(62360674, 62369312);
  transcriptFactorLUT["MTA3"] = make_pair(42795184, 42984087);
  transcriptFactorLUT["MTF1"] = make_pair(38275238, 38325292);
  transcriptFactorLUT["MXD1"] = make_pair(70142172, 70170076);
  transcriptFactorLUT["MXD3"] = make_pair(176734205, 176739292);
  transcriptFactorLUT["MXD4"] = make_pair(2249159, 2263739);
  transcriptFactorLUT["MXI1"] = make_pair(111985761, 112047123);
  transcriptFactorLUT["MYB"] = make_pair(135502452, 135540311);
  transcriptFactorLUT["MYBL1"] = make_pair(67474409, 67525484);
  transcriptFactorLUT["MYBL2"] = make_pair(42295658, 42345136);
  transcriptFactorLUT["MYC"] = make_pair(128748314, 128753680);
  transcriptFactorLUT["MYCN"] = make_pair(16080559, 16087129);
  transcriptFactorLUT["MYF5"] = make_pair(81110707, 81113447);
  transcriptFactorLUT["MYF6"] = make_pair(81101407, 81103256);
  transcriptFactorLUT["MYNN"] = make_pair(169490852, 169507504);
  transcriptFactorLUT["MYOD1"] = make_pair(17741109, 17743678);
  transcriptFactorLUT["MYOG"] = make_pair(203052256, 203055166);
  transcriptFactorLUT["MYSM1"] = make_pair(59120410, 59165747);
  transcriptFactorLUT["MYT1"] = make_pair(62795826, 62873606);
  transcriptFactorLUT["MYT1L"] = make_pair(1792884, 2335085);
  transcriptFactorLUT["MZF1"] = make_pair(59073283, 59084942);
  transcriptFactorLUT["NANOG"] = make_pair(7941991, 7948657);
  transcriptFactorLUT["NCOR1"] = make_pair(15933407, 16118874);
  transcriptFactorLUT["NCOR2"] = make_pair(124808956, 125052010);
  transcriptFactorLUT["NEUROD1"] = make_pair(182540832, 182545392);
  transcriptFactorLUT["NEUROD2"] = make_pair(37760020, 37764175);
  transcriptFactorLUT["NEUROD4"] = make_pair(55413728, 55423801);
  transcriptFactorLUT["NEUROD6"] = make_pair(31377074, 31380538);
  transcriptFactorLUT["NEUROG1"] = make_pair(134869971, 134871639);
  transcriptFactorLUT["NEUROG2"] = make_pair(113434671, 113437328);
  transcriptFactorLUT["NEUROG3"] = make_pair(71331790, 71333210);
  transcriptFactorLUT["NFAT5"] = make_pair(69599868, 69738569);
  transcriptFactorLUT["NFATC1"] = make_pair(77155771, 77289323);
  transcriptFactorLUT["NFATC2"] = make_pair(50003493, 50159258);
  transcriptFactorLUT["NFATC3"] = make_pair(68119268, 68263162);
  transcriptFactorLUT["NFATC4"] = make_pair(24836116, 24848811);
  transcriptFactorLUT["NFE2"] = make_pair(54685890, 54689563);
  transcriptFactorLUT["NFE2L1"] = make_pair(46125685, 46138907);
  transcriptFactorLUT["NFE2L2"] = make_pair(178095030, 178129859);
  transcriptFactorLUT["NFE2L3"] = make_pair(26191846, 26226756);
  transcriptFactorLUT["NFIA"] = make_pair(61547979, 61928460);
  transcriptFactorLUT["NFIB"] = make_pair(14081841, 14398982);
  transcriptFactorLUT["NFIC"] = make_pair(3359560, 3469215);
  transcriptFactorLUT["NFIL3"] = make_pair(94171326, 94186908);
  transcriptFactorLUT["NFIX"] = make_pair(13106583, 13209610);
  transcriptFactorLUT["NFKB1"] = make_pair(103422485, 103538459);
  transcriptFactorLUT["NFKB2"] = make_pair(104154334, 104162286);
  transcriptFactorLUT["NFX1"] = make_pair(33290417, 33348721);
  transcriptFactorLUT["NFXL1"] = make_pair(47849249, 47916574);
  transcriptFactorLUT["NFYA"] = make_pair(41040706, 41070146);
  transcriptFactorLUT["NFYB"] = make_pair(104510857, 104532040);
  transcriptFactorLUT["NFYC"] = make_pair(41157241, 41237275);
  transcriptFactorLUT["NHLH1"] = make_pair(160336860, 160342638);
  transcriptFactorLUT["NHLH2"] = make_pair(116378998, 116383333);
  transcriptFactorLUT["NKX1-1"] = make_pair(1396719, 1400230);
  transcriptFactorLUT["NKX2-1"] = make_pair(36985603, 36988903);
  transcriptFactorLUT["NKX2-2"] = make_pair(21491659, 21494664);
  transcriptFactorLUT["NKX2-3"] = make_pair(101292689, 101296280);
  transcriptFactorLUT["NKX2-4"] = make_pair(21376004, 21378047);
  transcriptFactorLUT["NKX2-5"] = make_pair(172659106, 172662315);
  transcriptFactorLUT["NKX2-6"] = make_pair(23559963, 23564111);
  transcriptFactorLUT["NKX2-8"] = make_pair(37049215, 37051786);
  transcriptFactorLUT["NKX3-1"] = make_pair(23536205, 23540434);
  transcriptFactorLUT["NKX3-2"] = make_pair(13542453, 13546114);
  transcriptFactorLUT["NKX6-1"] = make_pair(85414435, 85419387);
  transcriptFactorLUT["NKX6-2"] = make_pair(134598319, 134599537);
  transcriptFactorLUT["NKX6-3"] = make_pair(41503828, 41504878);
  transcriptFactorLUT["NOBOX"] = make_pair(144094332, 144107320);
  transcriptFactorLUT["NOTO"] = make_pair(73429385, 73438340);
  transcriptFactorLUT["NPAS1"] = make_pair(47524142, 47549017);
  transcriptFactorLUT["NPAS2"] = make_pair(101436612, 101613287);
  transcriptFactorLUT["NPAS3"] = make_pair(33408458, 34273382);
  transcriptFactorLUT["NPAS4"] = make_pair(66188474, 66194177);
  transcriptFactorLUT["NR0B1"] = make_pair(30322538, 30327495);
  transcriptFactorLUT["NR0B2"] = make_pair(27237974, 27240567);
  transcriptFactorLUT["NR1D1"] = make_pair(38249036, 38256978);
  transcriptFactorLUT["NR1D2"] = make_pair(23987611, 24022109);
  transcriptFactorLUT["NR1H2"] = make_pair(50879679, 50886285);
  transcriptFactorLUT["NR1H3"] = make_pair(47279467, 47290584);
  transcriptFactorLUT["NR1H4"] = make_pair(100867550, 100957645);
  transcriptFactorLUT["NR1I2"] = make_pair(119501556, 119537332);
  transcriptFactorLUT["NR1I3"] = make_pair(161199455, 161208000);
  transcriptFactorLUT["NR2C1"] = make_pair(95414004, 95467404);
  transcriptFactorLUT["NR2C2"] = make_pair(14989090, 15090786);
  transcriptFactorLUT["NR2E1"] = make_pair(108487261, 108510013);
  transcriptFactorLUT["NR2E3"] = make_pair(72102887, 72107270);
  transcriptFactorLUT["NR2F1"] = make_pair(92919042, 92930315);
  transcriptFactorLUT["NR2F2"] = make_pair(96874110, 96883492);
  transcriptFactorLUT["NR2F6"] = make_pair(17342693, 17356151);
  transcriptFactorLUT["NR3C1"] = make_pair(142657495, 142783254);
  transcriptFactorLUT["NR3C2"] = make_pair(148999914, 149363672);
  transcriptFactorLUT["NR4A1"] = make_pair(52416615, 52453291);
  transcriptFactorLUT["NR4A2"] = make_pair(157180943, 157189287);
  transcriptFactorLUT["NR4A3"] = make_pair(102584136, 102596341);
  transcriptFactorLUT["NR5A1"] = make_pair(127243514, 127269699);
  transcriptFactorLUT["NR5A2"] = make_pair(199996729, 200146550);
  transcriptFactorLUT["NR6A1"] = make_pair(127279553, 127533589);
  transcriptFactorLUT["NRF1"] = make_pair(129269918, 129396922);
  transcriptFactorLUT["NRL"] = make_pair(24549315, 24553832);
  transcriptFactorLUT["OLIG1"] = make_pair(34442449, 34444728);
  transcriptFactorLUT["OLIG2"] = make_pair(34398215, 34401503);
  transcriptFactorLUT["OLIG3"] = make_pair(137813335, 137815531);
  transcriptFactorLUT["ONECUT1"] = make_pair(53049159, 53082209);
  transcriptFactorLUT["ONECUT2"] = make_pair(55102916, 55158530);
  transcriptFactorLUT["ONECUT3"] = make_pair(1753661, 1775444);
  transcriptFactorLUT["OSR1"] = make_pair(19551245, 19558372);
  transcriptFactorLUT["OSR2"] = make_pair(99956630, 99964332);
  transcriptFactorLUT["OTP"] = make_pair(76924536, 76934522);
  transcriptFactorLUT["OTX1"] = make_pair(63277936, 63284966);
  transcriptFactorLUT["OTX2"] = make_pair(57267424, 57277194);
  transcriptFactorLUT["OVOL1"] = make_pair(65554504, 65564690);
  transcriptFactorLUT["OVOL2"] = make_pair(18004795, 18039832);
  transcriptFactorLUT["OVOL3"] = make_pair(36602104, 36604613);
  transcriptFactorLUT["PATZ1"] = make_pair(31721789, 31742249);
  transcriptFactorLUT["PAX1"] = make_pair(21686296, 21699124);
  transcriptFactorLUT["PAX2"] = make_pair(102495465, 102589698);
  transcriptFactorLUT["PAX3"] = make_pair(223064605, 223163715);
  transcriptFactorLUT["PAX4"] = make_pair(127250345, 127255780);
  transcriptFactorLUT["PAX5"] = make_pair(36833271, 37034476);
  transcriptFactorLUT["PAX6"] = make_pair(31806339, 31832901);
  transcriptFactorLUT["PAX7"] = make_pair(18957499, 19062632);
  transcriptFactorLUT["PAX8"] = make_pair(113973573, 114036498);
  transcriptFactorLUT["PAX9"] = make_pair(37126772, 37147011);
  transcriptFactorLUT["PBRM1"] = make_pair(52579367, 52719866);
  transcriptFactorLUT["PBX1"] = make_pair(164528596, 164821060);
  transcriptFactorLUT["PBX2"] = make_pair(32152509, 32157963);
  transcriptFactorLUT["PBX3"] = make_pair(128509616, 128729655);
  transcriptFactorLUT["PBX4"] = make_pair(19672515, 19729725);
  transcriptFactorLUT["PDX1"] = make_pair(28494167, 28500451);
  transcriptFactorLUT["PEG3"] = make_pair(57321444, 57352094);
  transcriptFactorLUT["PGR"] = make_pair(100900354, 100999794);
  transcriptFactorLUT["PHOX2A"] = make_pair(71950120, 71955220);
  transcriptFactorLUT["PHOX2B"] = make_pair(41746098, 41750987);
  transcriptFactorLUT["PIAS1"] = make_pair(68346571, 68480404);
  transcriptFactorLUT["PIAS2"] = make_pair(44390022, 44497495);
  transcriptFactorLUT["PIAS3"] = make_pair(145575987, 145586546);
  transcriptFactorLUT["PIAS4"] = make_pair(4007595, 4039384);
  transcriptFactorLUT["PINX1"] = make_pair(10622470, 10697409);
  transcriptFactorLUT["PITX1"] = make_pair(134363423, 134369964);
  transcriptFactorLUT["PITX2"] = make_pair(111538579, 111563279);
  transcriptFactorLUT["PITX3"] = make_pair(103989945, 104001231);
  transcriptFactorLUT["PKNOX1"] = make_pair(44394619, 44454041);
  transcriptFactorLUT["PKNOX2"] = make_pair(125034558, 125303285);
  transcriptFactorLUT["PLAG1"] = make_pair(57073467, 57123859);
  transcriptFactorLUT["PLAGL1"] = make_pair(144261436, 144329541);
  transcriptFactorLUT["PLAGL2"] = make_pair(30780306, 30795546);
  transcriptFactorLUT["PMS1"] = make_pair(190648810, 190742355);
  transcriptFactorLUT["POU1F1"] = make_pair(87308782, 87325737);
  transcriptFactorLUT["POU2F1"] = make_pair(167190065, 167396582);
  transcriptFactorLUT["POU2F2"] = make_pair(42590261, 42636625);
  transcriptFactorLUT["POU2F3"] = make_pair(120107348, 120190653);
  transcriptFactorLUT["POU3F1"] = make_pair(38509522, 38512450);
  transcriptFactorLUT["POU3F2"] = make_pair(99282579, 99286666);
  transcriptFactorLUT["POU3F3"] = make_pair(105470525, 105475031);
  transcriptFactorLUT["POU3F4"] = make_pair(82763268, 82764775);
  transcriptFactorLUT["POU4F1"] = make_pair(79173229, 79177695);
  transcriptFactorLUT["POU4F2"] = make_pair(147560044, 147563623);
  transcriptFactorLUT["POU4F3"] = make_pair(145718586, 145720083);
  transcriptFactorLUT["POU5F1"] = make_pair(31132113, 31134947);
  transcriptFactorLUT["POU5F1B"] = make_pair(128427856, 128429441);
  transcriptFactorLUT["POU6F1"] = make_pair(51580718, 51591950);
  transcriptFactorLUT["POU6F2"] = make_pair(39017608, 39504390);
  transcriptFactorLUT["PPARA"] = make_pair(46546498, 46639653);
  transcriptFactorLUT["PPARD"] = make_pair(35310334, 35395968);
  transcriptFactorLUT["PPARG"] = make_pair(12329348, 12475855);
  transcriptFactorLUT["PRDM1"] = make_pair(106534194, 106557814);
  transcriptFactorLUT["PRDM10"] = make_pair(129769600, 129872730);
  transcriptFactorLUT["PRDM12"] = make_pair(133539980, 133558384);
  transcriptFactorLUT["PRDM13"] = make_pair(100054649, 100063454);
  transcriptFactorLUT["PRDM14"] = make_pair(70963885, 70983562);
  transcriptFactorLUT["PRDM15"] = make_pair(43218384, 43283411);
  transcriptFactorLUT["PRDM16"] = make_pair(2985741, 3355185);
  transcriptFactorLUT["PRDM2"] = make_pair(14075875, 14114574);
  transcriptFactorLUT["PRDM4"] = make_pair(108126642, 108154914);
  transcriptFactorLUT["PRDM5"] = make_pair(121613067, 121844021);
  transcriptFactorLUT["PRDM6"] = make_pair(122424840, 122523745);
  transcriptFactorLUT["PRDM8"] = make_pair(81106423, 81125482);
  transcriptFactorLUT["PRDM9"] = make_pair(23507723, 23528706);
  transcriptFactorLUT["PREB"] = make_pair(27353624, 27357542);
  transcriptFactorLUT["PRKRIR"] = make_pair(76061000, 76092009);
  transcriptFactorLUT["PROP1"] = make_pair(177419235, 177423243);
  transcriptFactorLUT["PROX1"] = make_pair(214161843, 214214847);
  transcriptFactorLUT["PROX2"] = make_pair(75319735, 75330537);
  transcriptFactorLUT["PRRX1"] = make_pair(170633312, 170708541);
  transcriptFactorLUT["PRRX2"] = make_pair(132427919, 132484951);
  transcriptFactorLUT["PTF1A"] = make_pair(23481459, 23483181);
  transcriptFactorLUT["RARA"] = make_pair(38474472, 38513895);
  transcriptFactorLUT["RARB"] = make_pair(25469753, 25639422);
  transcriptFactorLUT["RARG"] = make_pair(53604349, 53614197);
  transcriptFactorLUT["RAX"] = make_pair(56934266, 56940625);
  transcriptFactorLUT["RAX2"] = make_pair(3769088, 3772219);
  transcriptFactorLUT["RBAK"] = make_pair(5085451, 5109119);
  transcriptFactorLUT["RCOR1"] = make_pair(103058995, 103196913);
  transcriptFactorLUT["RCOR2"] = make_pair(63678692, 63684316);
  transcriptFactorLUT["RCOR3"] = make_pair(211432707, 211486655);
  transcriptFactorLUT["REL"] = make_pair(61108629, 61155291);
  transcriptFactorLUT["RELA"] = make_pair(65421066, 65430443);
  transcriptFactorLUT["RELB"] = make_pair(45504706, 45541456);
  transcriptFactorLUT["REPIN1"] = make_pair(150065878, 150071133);
  transcriptFactorLUT["RERE"] = make_pair(8412463, 8877699);
  transcriptFactorLUT["REST"] = make_pair(57774041, 57802010);
  transcriptFactorLUT["RFX1"] = make_pair(14072341, 14117134);
  transcriptFactorLUT["RFX2"] = make_pair(5993174, 6110664);
  transcriptFactorLUT["RFX3"] = make_pair(3247036, 3395596);
  transcriptFactorLUT["RFX4"] = make_pair(106994914, 107156582);
  transcriptFactorLUT["RFX5"] = make_pair(151313115, 151319769);
  transcriptFactorLUT["RFX6"] = make_pair(117198375, 117253326);
  transcriptFactorLUT["RFX7"] = make_pair(56382730, 56535483);
  transcriptFactorLUT["RFX8"] = make_pair(102013822, 102091165);
  transcriptFactorLUT["RHOXF1"] = make_pair(119243010, 119249847);
  transcriptFactorLUT["RHOXF2"] = make_pair(119206240, 119211707);
  transcriptFactorLUT["RHOXF2B"] = make_pair(119206228, 119211707);
  transcriptFactorLUT["RNF138"] = make_pair(29672568, 29711524);
  transcriptFactorLUT["RORA"] = make_pair(60780482, 60884707);
  transcriptFactorLUT["RORB"] = make_pair(77112251, 77302117);
  transcriptFactorLUT["RORC"] = make_pair(151778546, 151804348);
  transcriptFactorLUT["RREB1"] = make_pair(7107829, 7252213);
  transcriptFactorLUT["RUNX1"] = make_pair(36160097, 36421595);
  transcriptFactorLUT["RUNX3"] = make_pair(25226001, 25291501);
  transcriptFactorLUT["RXRA"] = make_pair(137218308, 137332432);
  transcriptFactorLUT["RXRB"] = make_pair(33161361, 33168473);
  transcriptFactorLUT["SALL1"] = make_pair(51169885, 51185183);
  transcriptFactorLUT["SALL2"] = make_pair(21989230, 22005350);
  transcriptFactorLUT["SALL3"] = make_pair(76740274, 76758969);
  transcriptFactorLUT["SALL4"] = make_pair(50400582, 50419048);
  transcriptFactorLUT["SATB1"] = make_pair(18389132, 18466829);
  transcriptFactorLUT["SATB2"] = make_pair(200134222, 200322819);
  transcriptFactorLUT["SCRT1"] = make_pair(145554227, 145559943);
  transcriptFactorLUT["SCRT2"] = make_pair(642239, 656823);
  transcriptFactorLUT["SEBOX"] = make_pair(26691289, 26692173);
  transcriptFactorLUT["SETDB1"] = make_pair(150898814, 150917797);
  transcriptFactorLUT["SETDB2"] = make_pair(50025687, 50069139);
  transcriptFactorLUT["SHOX"] = make_pair(585078, 607558);
  transcriptFactorLUT["SHOX2"] = make_pair(157813799, 157823952);
  transcriptFactorLUT["SIM1"] = make_pair(100836749, 100911551);
  transcriptFactorLUT["SIX1"] = make_pair(61111416, 61116155);
  transcriptFactorLUT["SIX2"] = make_pair(45232323, 45236542);
  transcriptFactorLUT["SIX3"] = make_pair(45169036, 45173216);
  transcriptFactorLUT["SIX4"] = make_pair(61176255, 61190852);
  transcriptFactorLUT["SIX5"] = make_pair(46268042, 46272497);
  transcriptFactorLUT["SIX6"] = make_pair(60975937, 60978525);
  transcriptFactorLUT["SMAD1"] = make_pair(146402950, 146480325);
  transcriptFactorLUT["SMAD2"] = make_pair(45359465, 45456970);
  transcriptFactorLUT["SMAD3"] = make_pair(67358194, 67487533);
  transcriptFactorLUT["SMAD4"] = make_pair(48556582, 48611411);
  transcriptFactorLUT["SMAD5"] = make_pair(135468533, 135518422);
  transcriptFactorLUT["SMAD6"] = make_pair(66994673, 67074337);
  transcriptFactorLUT["SMAD7"] = make_pair(46446222, 46477081);
  transcriptFactorLUT["SMAD9"] = make_pair(37418967, 37494409);
  transcriptFactorLUT["SMARCC1"] = make_pair(47627377, 47823405);
  transcriptFactorLUT["SMARCC2"] = make_pair(56555635, 56583351);
  transcriptFactorLUT["SMARCE1"] = make_pair(38783975, 38804103);
  transcriptFactorLUT["SNAI1"] = make_pair(48599512, 48605420);
  transcriptFactorLUT["SNAI2"] = make_pair(49830238, 49833999);
  transcriptFactorLUT["SNAI3"] = make_pair(88744089, 88752882);
  transcriptFactorLUT["SNAPC4"] = make_pair(139270028, 139292889);
  transcriptFactorLUT["SOHLH1"] = make_pair(138585254, 138591374);
  transcriptFactorLUT["SOHLH2"] = make_pair(36742344, 36788752);
  transcriptFactorLUT["SOX1"] = make_pair(112721912, 112726020);
  transcriptFactorLUT["SOX10"] = make_pair(38368318, 38380539);
  transcriptFactorLUT["SOX11"] = make_pair(5832798, 5841517);
  transcriptFactorLUT["SOX12"] = make_pair(306214, 310872);
  transcriptFactorLUT["SOX13"] = make_pair(204042245, 204096871);
  transcriptFactorLUT["SOX14"] = make_pair(137483133, 137485172);
  transcriptFactorLUT["SOX15"] = make_pair(7491497, 7493488);
  transcriptFactorLUT["SOX17"] = make_pair(55370494, 55373456);
  transcriptFactorLUT["SOX18"] = make_pair(62679078, 62680979);
  transcriptFactorLUT["SOX2"] = make_pair(181429711, 181432223);
  transcriptFactorLUT["SOX21"] = make_pair(95361878, 95364797);
  transcriptFactorLUT["SOX3"] = make_pair(139585151, 139587225);
  transcriptFactorLUT["SOX30"] = make_pair(157052686, 157079428);
  transcriptFactorLUT["SOX4"] = make_pair(21593971, 21598849);
  transcriptFactorLUT["SOX5"] = make_pair(23685230, 24715383);
  transcriptFactorLUT["SOX6"] = make_pair(15987994, 16424413);
  transcriptFactorLUT["SOX7"] = make_pair(10581277, 10588084);
  transcriptFactorLUT["SOX8"] = make_pair(1031807, 1036979);
  transcriptFactorLUT["SOX9"] = make_pair(70117160, 70122560);
  transcriptFactorLUT["SP1"] = make_pair(53774427, 53810226);
  transcriptFactorLUT["SP100"] = make_pair(231280870, 231372861);
  transcriptFactorLUT["SP110"] = make_pair(231041188, 231084827);
  transcriptFactorLUT["SP140"] = make_pair(231090444, 231177930);
  transcriptFactorLUT["SP140L"] = make_pair(231191893, 231268445);
  transcriptFactorLUT["SP2"] = make_pair(45973515, 46006323);
  transcriptFactorLUT["SP3"] = make_pair(174771186, 174830430);
  transcriptFactorLUT["SP4"] = make_pair(21467688, 21554151);
  transcriptFactorLUT["SP5"] = make_pair(171571856, 171574498);
  transcriptFactorLUT["SP6"] = make_pair(45922279, 45933240);
  transcriptFactorLUT["SP7"] = make_pair(53720359, 53730167);
  transcriptFactorLUT["SP8"] = make_pair(20821893, 20826508);
  transcriptFactorLUT["SP9"] = make_pair(175199820, 175202268);
  transcriptFactorLUT["SPDEF"] = make_pair(34505578, 34524110);
  transcriptFactorLUT["SPI1"] = make_pair(47376408, 47400127);
  transcriptFactorLUT["SPIB"] = make_pair(50922194, 50934309);
  transcriptFactorLUT["SPIC"] = make_pair(101871334, 101880775);
  transcriptFactorLUT["SPZ1"] = make_pair(79615789, 79617660);
  transcriptFactorLUT["SREBF1"] = make_pair(17714662, 17740325);
  transcriptFactorLUT["SREBF2"] = make_pair(42229082, 42303312);
  transcriptFactorLUT["SRF"] = make_pair(43139032, 43149244);
  transcriptFactorLUT["SRY"] = make_pair(2654895, 2655782);
  transcriptFactorLUT["SSRP1"] = make_pair(57093458, 57103351);
  transcriptFactorLUT["ST18"] = make_pair(53023391, 53322439);
  transcriptFactorLUT["STAT1"] = make_pair(191833761, 191878976);
  transcriptFactorLUT["STAT2"] = make_pair(56735381, 56754037);
  transcriptFactorLUT["STAT3"] = make_pair(40465342, 40540405);
  transcriptFactorLUT["STAT4"] = make_pair(191894301, 192015986);
  transcriptFactorLUT["STAT5A"] = make_pair(40439564, 40463960);
  transcriptFactorLUT["STAT5B"] = make_pair(40351194, 40428424);
  transcriptFactorLUT["STAT6"] = make_pair(57489186, 57505196);
  transcriptFactorLUT["SUB1"] = make_pair(32585604, 32604185);
  transcriptFactorLUT["T"] = make_pair(166571145, 166582157);
  transcriptFactorLUT["TADA2A"] = make_pair(35767311, 35837226);
  transcriptFactorLUT["TADA2B"] = make_pair(7045155, 7059677);
  transcriptFactorLUT["TAL1"] = make_pair(47681961, 47698007);
  transcriptFactorLUT["TAL2"] = make_pair(108424737, 108425385);
  transcriptFactorLUT["TBR1"] = make_pair(162272619, 162281573);
  transcriptFactorLUT["TBX1"] = make_pair(19744225, 19771112);
  transcriptFactorLUT["TBX10"] = make_pair(67398773, 67407031);
  transcriptFactorLUT["TBX15"] = make_pair(119425665, 119532179);
  transcriptFactorLUT["TBX18"] = make_pair(85442215, 85473954);
  transcriptFactorLUT["TBX19"] = make_pair(168250277, 168283664);
  transcriptFactorLUT["TBX2"] = make_pair(59477256, 59486827);
  transcriptFactorLUT["TBX20"] = make_pair(35242041, 35293711);
  transcriptFactorLUT["TBX21"] = make_pair(45810609, 45823485);
  transcriptFactorLUT["TBX22"] = make_pair(79270254, 79287268);
  transcriptFactorLUT["TBX3"] = make_pair(115108058, 115121969);
  transcriptFactorLUT["TBX4"] = make_pair(59533806, 59561664);
  transcriptFactorLUT["TBX5"] = make_pair(114791734, 114846247);
  transcriptFactorLUT["TBX6"] = make_pair(30097114, 30103205);
  transcriptFactorLUT["TCF12"] = make_pair(57210832, 57580714);
  transcriptFactorLUT["TCF15"] = make_pair(584636, 590910);
  transcriptFactorLUT["TCF21"] = make_pair(134210258, 134216675);
  transcriptFactorLUT["TCF23"] = make_pair(27371944, 27375819);
  transcriptFactorLUT["TCF3"] = make_pair(1609288, 1652328);
  transcriptFactorLUT["TCF4"] = make_pair(52889561, 53255860);
  transcriptFactorLUT["TCF7"] = make_pair(133451349, 133483920);
  transcriptFactorLUT["TCF7L1"] = make_pair(85360582, 85537511);
  transcriptFactorLUT["TCF7L2"] = make_pair(114710008, 114927436);
  transcriptFactorLUT["TCFL5"] = make_pair(61472365, 61493115);
  transcriptFactorLUT["TEAD1"] = make_pair(12695968, 12966284);
  transcriptFactorLUT["TEAD2"] = make_pair(49843852, 49865714);
  transcriptFactorLUT["TEAD3"] = make_pair(35441373, 35464861);
  transcriptFactorLUT["TEAD4"] = make_pair(3068747, 3149842);
  transcriptFactorLUT["TEF"] = make_pair(41763336, 41795332);
  transcriptFactorLUT["TERF1"] = make_pair(73921096, 73959987);
  transcriptFactorLUT["TERF2"] = make_pair(69389463, 69419891);
  transcriptFactorLUT["TFAM"] = make_pair(60144902, 60158990);
  transcriptFactorLUT["TFAP2A"] = make_pair(10396915, 10412607);
  transcriptFactorLUT["TFAP2B"] = make_pair(50786438, 50815326);
  transcriptFactorLUT["TFAP2C"] = make_pair(55204357, 55214338);
  transcriptFactorLUT["TFAP2D"] = make_pair(50681256, 50740746);
  transcriptFactorLUT["TFAP2E"] = make_pair(36038970, 36060927);
  transcriptFactorLUT["TFAP4"] = make_pair(4307186, 4323001);
  transcriptFactorLUT["TFCP2"] = make_pair(51487538, 51566926);
  transcriptFactorLUT["TFCP2L1"] = make_pair(121974163, 122042778);
  transcriptFactorLUT["TFDP1"] = make_pair(114239055, 114295788);
  transcriptFactorLUT["TFDP2"] = make_pair(141663269, 141719229);
  transcriptFactorLUT["TFDP3"] = make_pair(132350696, 132352376);
  transcriptFactorLUT["TFE3"] = make_pair(48886237, 48901043);
  transcriptFactorLUT["TFEB"] = make_pair(41651715, 41703997);
  transcriptFactorLUT["TFEC"] = make_pair(115575201, 115608367);
  transcriptFactorLUT["TGIF1"] = make_pair(3453771, 3458406);
  transcriptFactorLUT["TGIF2"] = make_pair(35202956, 35222355);
  transcriptFactorLUT["TGIF2LX"] = make_pair(89176939, 89177882);
  transcriptFactorLUT["TGIF2LY"] = make_pair(3447125, 3448082);
  transcriptFactorLUT["THAP1"] = make_pair(42691816, 42698474);
  transcriptFactorLUT["THAP10"] = make_pair(71173680, 71184772);
  transcriptFactorLUT["THAP11"] = make_pair(67876212, 67878098);
  transcriptFactorLUT["THAP2"] = make_pair(72057676, 72074428);
  transcriptFactorLUT["THAP3"] = make_pair(6684924, 6693642);
  transcriptFactorLUT["THAP4"] = make_pair(242523819, 242556916);
  transcriptFactorLUT["THAP5"] = make_pair(108202587, 108209897);
  transcriptFactorLUT["THAP6"] = make_pair(76439653, 76455236);
  transcriptFactorLUT["THAP7"] = make_pair(21354060, 21356404);
  transcriptFactorLUT["THAP8"] = make_pair(36525886, 36545664);
  transcriptFactorLUT["THAP9"] = make_pair(83821836, 83841284);
  transcriptFactorLUT["THRA"] = make_pair(38219067, 38250120);
  transcriptFactorLUT["THRB"] = make_pair(24158644, 24536313);
  transcriptFactorLUT["TLX1"] = make_pair(102891060, 102897546);
  transcriptFactorLUT["TLX2"] = make_pair(74741595, 74744275);
  transcriptFactorLUT["TLX3"] = make_pair(170736287, 170739138);
  transcriptFactorLUT["TOX"] = make_pair(59717976, 60031767);
  transcriptFactorLUT["TOX2"] = make_pair(42574535, 42698254);
  transcriptFactorLUT["TOX3"] = make_pair(52471917, 52580806);
  transcriptFactorLUT["TOX4"] = make_pair(21945334, 21967321);
  transcriptFactorLUT["TP53"] = make_pair(7571719, 7590868);
  transcriptFactorLUT["TP63"] = make_pair(189349215, 189615068);
  transcriptFactorLUT["TP73"] = make_pair(3569128, 3652765);
  transcriptFactorLUT["TRERF1"] = make_pair(42192670, 42419789);
  transcriptFactorLUT["TRPS1"] = make_pair(116420723, 116681255);
  transcriptFactorLUT["TSC22D1"] = make_pair(45006278, 45150701);
  transcriptFactorLUT["TSC22D2"] = make_pair(150126121, 150177905);
  transcriptFactorLUT["TSC22D3"] = make_pair(106956451, 106960291);
  transcriptFactorLUT["TSC22D4"] = make_pair(100064141, 100076902);
  transcriptFactorLUT["TSHZ1"] = make_pair(72922709, 73001905);
  transcriptFactorLUT["TSHZ3"] = make_pair(31765850, 31840190);
  transcriptFactorLUT["TTF1"] = make_pair(135250936, 135282238);
  transcriptFactorLUT["TUB"] = make_pair(8060179, 8127654);
  transcriptFactorLUT["TULP1"] = make_pair(35465650, 35480679);
  transcriptFactorLUT["TULP2"] = make_pair(49384221, 49401996);
  transcriptFactorLUT["TULP3"] = make_pair(3000032, 3050306);
  transcriptFactorLUT["TULP4"] = make_pair(158733691, 158932856);
  transcriptFactorLUT["TWIST1"] = make_pair(19155090, 19157295);
  transcriptFactorLUT["TWIST2"] = make_pair(239756672, 239832244);
  transcriptFactorLUT["UBP1"] = make_pair(33429828, 33481870);
  transcriptFactorLUT["UBTF"] = make_pair(42282400, 42297041);
  transcriptFactorLUT["UBTFL1"] = make_pair(89819117, 89820299);
  transcriptFactorLUT["UNCX"] = make_pair(1272653, 1276613);
  transcriptFactorLUT["USF1"] = make_pair(161009040, 161015769);
  transcriptFactorLUT["USF2"] = make_pair(35759895, 35770718);
  transcriptFactorLUT["VAX1"] = make_pair(118892800, 118897812);
  transcriptFactorLUT["VAX2"] = make_pair(71127719, 71160575);
  transcriptFactorLUT["VDR"] = make_pair(48235319, 48298814);
  transcriptFactorLUT["VENTX"] = make_pair(135051407, 135055434);
  transcriptFactorLUT["VEZF1"] = make_pair(56048909, 56065615);
  transcriptFactorLUT["VSX1"] = make_pair(25059178, 25063015);
  transcriptFactorLUT["VSX2"] = make_pair(74706174, 74729441);
  transcriptFactorLUT["VTN"] = make_pair(26694298, 26697373);
  transcriptFactorLUT["WHSC1"] = make_pair(1894508, 1983934);
  transcriptFactorLUT["WIZ"] = make_pair(15532317, 15560762);
  transcriptFactorLUT["WT1"] = make_pair(32409321, 32452363);
  transcriptFactorLUT["XBP1"] = make_pair(29190547, 29196560);
  transcriptFactorLUT["YBX1"] = make_pair(43148065, 43168020);
  transcriptFactorLUT["YBX2"] = make_pair(7191570, 7197876);
  transcriptFactorLUT["YY1"] = make_pair(100705101, 100745371);
  transcriptFactorLUT["YY2"] = make_pair(21874104, 21876845);
  transcriptFactorLUT["ZBED1"] = make_pair(2404454, 2419008);
  transcriptFactorLUT["ZBED2"] = make_pair(111311746, 111314182);
  transcriptFactorLUT["ZBED3"] = make_pair(76372531, 76383030);
  transcriptFactorLUT["ZBED4"] = make_pair(50247496, 50283726);
  transcriptFactorLUT["ZBTB1"] = make_pair(64971291, 65000408);
  transcriptFactorLUT["ZBTB10"] = make_pair(81398447, 81438500);
  transcriptFactorLUT["ZBTB11"] = make_pair(101368282, 101395988);
  transcriptFactorLUT["ZBTB12"] = make_pair(31867393, 31869769);
  transcriptFactorLUT["ZBTB16"] = make_pair(113930430, 114121397);
  transcriptFactorLUT["ZBTB17"] = make_pair(16268363, 16302627);
  transcriptFactorLUT["ZBTB2"] = make_pair(151685249, 151712835);
  transcriptFactorLUT["ZBTB20"] = make_pair(114033346, 114866127);
  transcriptFactorLUT["ZBTB24"] = make_pair(109802151, 109804440);
  transcriptFactorLUT["ZBTB25"] = make_pair(64953554, 64970563);
  transcriptFactorLUT["ZBTB26"] = make_pair(125680307, 125693830);
  transcriptFactorLUT["ZBTB3"] = make_pair(62518434, 62521656);
  transcriptFactorLUT["ZBTB32"] = make_pair(36203829, 36207940);
  transcriptFactorLUT["ZBTB33"] = make_pair(119384609, 119392251);
  transcriptFactorLUT["ZBTB34"] = make_pair(129622943, 129648156);
  transcriptFactorLUT["ZBTB37"] = make_pair(173837492, 173855774);
  transcriptFactorLUT["ZBTB38"] = make_pair(141043054, 141168632);
  transcriptFactorLUT["ZBTB4"] = make_pair(7362684, 7382944);
  transcriptFactorLUT["ZBTB40"] = make_pair(22778343, 22857650);
  transcriptFactorLUT["ZBTB41"] = make_pair(197122813, 197169672);
  transcriptFactorLUT["ZBTB42"] = make_pair(105267517, 105271048);
  transcriptFactorLUT["ZBTB43"] = make_pair(129567284, 129600487);
  transcriptFactorLUT["ZBTB44"] = make_pair(130096573, 130184607);
  transcriptFactorLUT["ZBTB45"] = make_pair(59024896, 59030921);
  transcriptFactorLUT["ZBTB46"] = make_pair(62375020, 62436856);
  transcriptFactorLUT["ZBTB47"] = make_pair(42695175, 42709072);
  transcriptFactorLUT["ZBTB48"] = make_pair(6640050, 6649340);
  transcriptFactorLUT["ZBTB49"] = make_pair(4291923, 4323513);
  transcriptFactorLUT["ZBTB6"] = make_pair(125670334, 125675607);
  transcriptFactorLUT["ZBTB7A"] = make_pair(4045215, 4066816);
  transcriptFactorLUT["ZBTB7B"] = make_pair(154975105, 154991001);
  transcriptFactorLUT["ZBTB7C"] = make_pair(45553638, 45663680);
  transcriptFactorLUT["ZBTB8A"] = make_pair(33004745, 33071551);
  transcriptFactorLUT["ZEB1"] = make_pair(31610063, 31818742);
  transcriptFactorLUT["ZEB2"] = make_pair(145141941, 145277958);
  transcriptFactorLUT["ZFAT"] = make_pair(135490030, 135708801);
  transcriptFactorLUT["ZFHX3"] = make_pair(72816785, 73082274);
  transcriptFactorLUT["ZFHX4"] = make_pair(77593514, 77779521);
  transcriptFactorLUT["ZFP1"] = make_pair(75182420, 75206132);
  transcriptFactorLUT["ZFP14"] = make_pair(36825354, 36870105);
  transcriptFactorLUT["ZFP2"] = make_pair(178322915, 178360210);
  transcriptFactorLUT["ZFP28"] = make_pair(57050316, 57068170);
  transcriptFactorLUT["ZFP3"] = make_pair(4981753, 4999669);
  transcriptFactorLUT["ZFP30"] = make_pair(38123388, 38146313);
  transcriptFactorLUT["ZFP37"] = make_pair(115804094, 115819071);
  transcriptFactorLUT["ZFP41"] = make_pair(144328990, 144335761);
  transcriptFactorLUT["ZFP42"] = make_pair(188916924, 188926203);
  transcriptFactorLUT["ZFP57"] = make_pair(29640168, 29644931);
  transcriptFactorLUT["ZFP62"] = make_pair(180274610, 180288286);
  transcriptFactorLUT["ZFP64"] = make_pair(50700549, 50808524);
  transcriptFactorLUT["ZFP82"] = make_pair(36882860, 36909550);
  transcriptFactorLUT["ZFP90"] = make_pair(68573115, 68609975);
  transcriptFactorLUT["ZFP91"] = make_pair(58346586, 58389023);
  transcriptFactorLUT["ZFPM1"] = make_pair(88520013, 88601574);
  transcriptFactorLUT["ZFPM2"] = make_pair(106331146, 106816767);
  transcriptFactorLUT["ZFX"] = make_pair(24167761, 24234372);
  transcriptFactorLUT["ZFY"] = make_pair(2803517, 2850547);
  transcriptFactorLUT["ZGLP1"] = make_pair(10415478, 10420233);
  transcriptFactorLUT["ZHX1"] = make_pair(124260689, 124286727);
  transcriptFactorLUT["ZHX2"] = make_pair(123793900, 123986755);
  transcriptFactorLUT["ZHX3"] = make_pair(39807088, 39928739);
  transcriptFactorLUT["ZIC1"] = make_pair(147127180, 147134506);
  transcriptFactorLUT["ZIC2"] = make_pair(100634025, 100639019);
  transcriptFactorLUT["ZIC3"] = make_pair(136648345, 136654259);
  transcriptFactorLUT["ZIC4"] = make_pair(147103834, 147124596);
  transcriptFactorLUT["ZIC5"] = make_pair(100615274, 100624178);
  transcriptFactorLUT["ZIK1"] = make_pair(58095627, 58103758);
  transcriptFactorLUT["ZIM2"] = make_pair(57285922, 57352097);
  transcriptFactorLUT["ZIM3"] = make_pair(57645463, 57656570);
  transcriptFactorLUT["ZKSCAN1"] = make_pair(99613194, 99639312);
  transcriptFactorLUT["ZKSCAN2"] = make_pair(25247321, 25268855);
  transcriptFactorLUT["ZKSCAN3"] = make_pair(28317690, 28336954);
  transcriptFactorLUT["ZKSCAN4"] = make_pair(28209482, 28220074);
  transcriptFactorLUT["ZKSCAN5"] = make_pair(99102272, 99131445);
  transcriptFactorLUT["ZMIZ1"] = make_pair(80828791, 81076285);
  transcriptFactorLUT["ZMIZ2"] = make_pair(44795786, 44809479);
  transcriptFactorLUT["ZNF10"] = make_pair(133707213, 133736049);
  transcriptFactorLUT["ZNF100"] = make_pair(21906842, 21950430);
  transcriptFactorLUT["ZNF101"] = make_pair(19779596, 19794315);
  transcriptFactorLUT["ZNF107"] = make_pair(64126460, 64171960);
  transcriptFactorLUT["ZNF114"] = make_pair(48777006, 48790865);
  transcriptFactorLUT["ZNF117"] = make_pair(64434829, 64451414);
  transcriptFactorLUT["ZNF12"] = make_pair(6728063, 6746566);
  transcriptFactorLUT["ZNF121"] = make_pair(9676291, 9695209);
  transcriptFactorLUT["ZNF124"] = make_pair(247318290, 247335319);
  transcriptFactorLUT["ZNF131"] = make_pair(43121595, 43176426);
  transcriptFactorLUT["ZNF132"] = make_pair(58944180, 58951589);
  transcriptFactorLUT["ZNF133"] = make_pair(18268926, 18297640);
  transcriptFactorLUT["ZNF134"] = make_pair(58125829, 58133636);
  transcriptFactorLUT["ZNF135"] = make_pair(58570606, 58581110);
  transcriptFactorLUT["ZNF136"] = make_pair(12273871, 12300064);
  transcriptFactorLUT["ZNF138"] = make_pair(64254765, 64294059);
  transcriptFactorLUT["ZNF14"] = make_pair(19821280, 19843921);
  transcriptFactorLUT["ZNF140"] = make_pair(133656788, 133684258);
  transcriptFactorLUT["ZNF141"] = make_pair(331595, 367691);
  transcriptFactorLUT["ZNF142"] = make_pair(219502639, 219524355);
  transcriptFactorLUT["ZNF143"] = make_pair(9482511, 9550071);
  transcriptFactorLUT["ZNF146"] = make_pair(36705503, 36729675);
  transcriptFactorLUT["ZNF148"] = make_pair(124944512, 125094198);
  transcriptFactorLUT["ZNF154"] = make_pair(58207642, 58220579);
  transcriptFactorLUT["ZNF155"] = make_pair(44488321, 44502477);
  transcriptFactorLUT["ZNF157"] = make_pair(47229998, 47273098);
  transcriptFactorLUT["ZNF16"] = make_pair(146155743, 146176274);
  transcriptFactorLUT["ZNF160"] = make_pair(53569866, 53606687);
  transcriptFactorLUT["ZNF165"] = make_pair(28048481, 28057340);
  transcriptFactorLUT["ZNF169"] = make_pair(97021547, 97064111);
  transcriptFactorLUT["ZNF17"] = make_pair(57922528, 57933307);
  transcriptFactorLUT["ZNF174"] = make_pair(3451189, 3459364);
  transcriptFactorLUT["ZNF175"] = make_pair(52074530, 52092991);
  transcriptFactorLUT["ZNF177"] = make_pair(9486991, 9493293);
  transcriptFactorLUT["ZNF18"] = make_pair(11880755, 11900827);
  transcriptFactorLUT["ZNF180"] = make_pair(44978644, 45004575);
  transcriptFactorLUT["ZNF181"] = make_pair(35225479, 35233774);
  transcriptFactorLUT["ZNF182"] = make_pair(47834249, 47863377);
  transcriptFactorLUT["ZNF184"] = make_pair(27418520, 27440897);
  transcriptFactorLUT["ZNF189"] = make_pair(104161135, 104172942);
  transcriptFactorLUT["ZNF19"] = make_pair(71507975, 71523254);
  transcriptFactorLUT["ZNF195"] = make_pair(3379156, 3400452);
  transcriptFactorLUT["ZNF197"] = make_pair(44666510, 44686752);
  transcriptFactorLUT["ZNF2"] = make_pair(95831161, 95850064);
  transcriptFactorLUT["ZNF20"] = make_pair(12242167, 12251222);
  transcriptFactorLUT["ZNF200"] = make_pair(3272324, 3285456);
  transcriptFactorLUT["ZNF202"] = make_pair(123594634, 123612391);
  transcriptFactorLUT["ZNF205"] = make_pair(3162562, 3170518);
  transcriptFactorLUT["ZNF208"] = make_pair(22148896, 22193745);
  transcriptFactorLUT["ZNF211"] = make_pair(58144534, 58154147);
  transcriptFactorLUT["ZNF212"] = make_pair(148936741, 148952700);
  transcriptFactorLUT["ZNF213"] = make_pair(3185056, 3192805);
  transcriptFactorLUT["ZNF214"] = make_pair(7020548, 7041586);
  transcriptFactorLUT["ZNF215"] = make_pair(6947653, 6979278);
  transcriptFactorLUT["ZNF217"] = make_pair(52183609, 52199636);
  transcriptFactorLUT["ZNF219"] = make_pair(21558204, 21567173);
  transcriptFactorLUT["ZNF22"] = make_pair(45496272, 45500777);
  transcriptFactorLUT["ZNF221"] = make_pair(44455396, 44471752);
  transcriptFactorLUT["ZNF222"] = make_pair(44529493, 44537262);
  transcriptFactorLUT["ZNF223"] = make_pair(44556163, 44572147);
  transcriptFactorLUT["ZNF224"] = make_pair(44598481, 44612479);
  transcriptFactorLUT["ZNF225"] = make_pair(44617547, 44637255);
  transcriptFactorLUT["ZNF226"] = make_pair(44669248, 44679582);
  transcriptFactorLUT["ZNF227"] = make_pair(44716680, 44741421);
  transcriptFactorLUT["ZNF229"] = make_pair(44930422, 44952665);
  transcriptFactorLUT["ZNF23"] = make_pair(71481502, 71496155);
  transcriptFactorLUT["ZNF230"] = make_pair(44507076, 44518072);
  transcriptFactorLUT["ZNF232"] = make_pair(5009030, 5026397);
  transcriptFactorLUT["ZNF233"] = make_pair(44764032, 44779470);
  transcriptFactorLUT["ZNF234"] = make_pair(44645709, 44664462);
  transcriptFactorLUT["ZNF235"] = make_pair(44790500, 44809178);
  transcriptFactorLUT["ZNF236"] = make_pair(74536115, 74682682);
  transcriptFactorLUT["ZNF239"] = make_pair(44051792, 44063907);
  transcriptFactorLUT["ZNF24"] = make_pair(32912177, 32924426);
  transcriptFactorLUT["ZNF248"] = make_pair(38065453, 38146564);
  transcriptFactorLUT["ZNF25"] = make_pair(38238794, 38265453);
  transcriptFactorLUT["ZNF250"] = make_pair(146102335, 146126846);
  transcriptFactorLUT["ZNF251"] = make_pair(145946293, 145980970);
  transcriptFactorLUT["ZNF253"] = make_pair(19976713, 20004293);
  transcriptFactorLUT["ZNF254"] = make_pair(24216206, 24312769);
  transcriptFactorLUT["ZNF256"] = make_pair(58452200, 58459077);
  transcriptFactorLUT["ZNF257"] = make_pair(22235265, 22273903);
  transcriptFactorLUT["ZNF26"] = make_pair(133562933, 133589154);
  transcriptFactorLUT["ZNF263"] = make_pair(3333486, 3341459);
  transcriptFactorLUT["ZNF264"] = make_pair(57702867, 57734214);
  transcriptFactorLUT["ZNF266"] = make_pair(9523101, 9546254);
  transcriptFactorLUT["ZNF267"] = make_pair(31885078, 31928629);
  transcriptFactorLUT["ZNF268"] = make_pair(133757994, 133783697);
  transcriptFactorLUT["ZNF273"] = make_pair(64363619, 64391955);
  transcriptFactorLUT["ZNF274"] = make_pair(58694355, 58724928);
  transcriptFactorLUT["ZNF275"] = make_pair(152599612, 152618384);
  transcriptFactorLUT["ZNF276"] = make_pair(89787951, 89807332);
  transcriptFactorLUT["ZNF28"] = make_pair(53300660, 53324922);
  transcriptFactorLUT["ZNF280D"] = make_pair(56922373, 57026284);
  transcriptFactorLUT["ZNF281"] = make_pair(200374074, 200379186);
  transcriptFactorLUT["ZNF282"] = make_pair(148892553, 148923339);
  transcriptFactorLUT["ZNF283"] = make_pair(44331472, 44353050);
  transcriptFactorLUT["ZNF284"] = make_pair(44576296, 44591623);
  transcriptFactorLUT["ZNF285"] = make_pair(44889807, 44905777);
  transcriptFactorLUT["ZNF286A"] = make_pair(15602890, 15624100);
  transcriptFactorLUT["ZNF287"] = make_pair(16453630, 16472520);
  transcriptFactorLUT["ZNF292"] = make_pair(87865268, 87973406);
  transcriptFactorLUT["ZNF296"] = make_pair(45574757, 45579688);
  transcriptFactorLUT["ZNF3"] = make_pair(99667593, 99680171);
  transcriptFactorLUT["ZNF30"] = make_pair(35417806, 35436076);
  transcriptFactorLUT["ZNF300"] = make_pair(150273953, 150284545);
  transcriptFactorLUT["ZNF302"] = make_pair(35168543, 35177302);
  transcriptFactorLUT["ZNF304"] = make_pair(57862641, 57871265);
  transcriptFactorLUT["ZNF311"] = make_pair(28962593, 28973037);
  transcriptFactorLUT["ZNF317"] = make_pair(9251055, 9274091);
  transcriptFactorLUT["ZNF319"] = make_pair(58028571, 58033762);
  transcriptFactorLUT["ZNF32"] = make_pair(44139306, 44144326);
  transcriptFactorLUT["ZNF320"] = make_pair(53379424, 53394599);
  transcriptFactorLUT["ZNF324"] = make_pair(58978411, 58984945);
  transcriptFactorLUT["ZNF324B"] = make_pair(58962970, 58969199);
  transcriptFactorLUT["ZNF329"] = make_pair(58637694, 58662148);
  transcriptFactorLUT["ZNF331"] = make_pair(54041539, 54083523);
  transcriptFactorLUT["ZNF333"] = make_pair(14800869, 14831772);
  transcriptFactorLUT["ZNF334"] = make_pair(45128268, 45142198);
  transcriptFactorLUT["ZNF335"] = make_pair(44577291, 44600833);
  transcriptFactorLUT["ZNF337"] = make_pair(25653830, 25667575);
  transcriptFactorLUT["ZNF33A"] = make_pair(38299577, 38348995);
  transcriptFactorLUT["ZNF33B"] = make_pair(43084531, 43134285);
  transcriptFactorLUT["ZNF34"] = make_pair(145997608, 146012730);
  transcriptFactorLUT["ZNF341"] = make_pair(32319565, 32380075);
  transcriptFactorLUT["ZNF343"] = make_pair(2462465, 2505168);
  transcriptFactorLUT["ZNF345"] = make_pair(37341259, 37370477);
  transcriptFactorLUT["ZNF347"] = make_pair(53641956, 53662322);
  transcriptFactorLUT["ZNF35"] = make_pair(44690232, 44702283);
  transcriptFactorLUT["ZNF350"] = make_pair(52467592, 52490079);
  transcriptFactorLUT["ZNF354A"] = make_pair(178138521, 178157703);
  transcriptFactorLUT["ZNF354B"] = make_pair(178286953, 178311424);
  transcriptFactorLUT["ZNF354C"] = make_pair(178487606, 178507691);
  transcriptFactorLUT["ZNF358"] = make_pair(7581003, 7585911);
  transcriptFactorLUT["ZNF362"] = make_pair(33722173, 33766320);
  transcriptFactorLUT["ZNF366"] = make_pair(71735725, 71803249);
  transcriptFactorLUT["ZNF367"] = make_pair(99148224, 99180669);
  transcriptFactorLUT["ZNF37A"] = make_pair(38383263, 38412278);
  transcriptFactorLUT["ZNF382"] = make_pair(37096206, 37119499);
  transcriptFactorLUT["ZNF383"] = make_pair(37717365, 37734574);
  transcriptFactorLUT["ZNF384"] = make_pair(6775642, 6798541);
  transcriptFactorLUT["ZNF391"] = make_pair(27356523, 27369227);
  transcriptFactorLUT["ZNF394"] = make_pair(99090853, 99097877);
  transcriptFactorLUT["ZNF396"] = make_pair(32946660, 32957301);
  transcriptFactorLUT["ZNF397"] = make_pair(32820993, 32838397);
  transcriptFactorLUT["ZNF398"] = make_pair(148844559, 148880118);
  transcriptFactorLUT["ZNF407"] = make_pair(72342918, 72516583);
  transcriptFactorLUT["ZNF408"] = make_pair(46722316, 46727466);
  transcriptFactorLUT["ZNF41"] = make_pair(47305560, 47342345);
  transcriptFactorLUT["ZNF410"] = make_pair(74353317, 74398991);
  transcriptFactorLUT["ZNF415"] = make_pair(53611131, 53636171);
  transcriptFactorLUT["ZNF416"] = make_pair(58082933, 58090243);
  transcriptFactorLUT["ZNF417"] = make_pair(58417141, 58427978);
  transcriptFactorLUT["ZNF418"] = make_pair(58433251, 58446740);
  transcriptFactorLUT["ZNF419"] = make_pair(57999078, 58006048);
  transcriptFactorLUT["ZNF420"] = make_pair(37569381, 37620651);
  transcriptFactorLUT["ZNF423"] = make_pair(49524514, 49891830);
  transcriptFactorLUT["ZNF425"] = make_pair(148799877, 148823438);
  transcriptFactorLUT["ZNF426"] = make_pair(9638680, 9649303);
  transcriptFactorLUT["ZNF429"] = make_pair(21688436, 21721079);
  transcriptFactorLUT["ZNF43"] = make_pair(21987750, 22034870);
  transcriptFactorLUT["ZNF430"] = make_pair(21203425, 21242852);
  transcriptFactorLUT["ZNF431"] = make_pair(21324839, 21368805);
  transcriptFactorLUT["ZNF432"] = make_pair(52536676, 52552073);
  transcriptFactorLUT["ZNF433"] = make_pair(12125531, 12146525);
  transcriptFactorLUT["ZNF436"] = make_pair(23685940, 23694879);
  transcriptFactorLUT["ZNF439"] = make_pair(11976843, 11980306);
  transcriptFactorLUT["ZNF44"] = make_pair(12382624, 12405714);
  transcriptFactorLUT["ZNF440"] = make_pair(11925106, 11946016);
  transcriptFactorLUT["ZNF441"] = make_pair(11877814, 11894893);
  transcriptFactorLUT["ZNF442"] = make_pair(12460184, 12476475);
  transcriptFactorLUT["ZNF443"] = make_pair(12540519, 12551926);
  transcriptFactorLUT["ZNF444"] = make_pair(56652534, 56672262);
  transcriptFactorLUT["ZNF445"] = make_pair(44481261, 44519162);
  transcriptFactorLUT["ZNF446"] = make_pair(58987530, 58992601);
  transcriptFactorLUT["ZNF449"] = make_pair(134478695, 134497338);
  transcriptFactorLUT["ZNF45"] = make_pair(44416775, 44439411);
  transcriptFactorLUT["ZNF451"] = make_pair(56954827, 57035098);
  transcriptFactorLUT["ZNF454"] = make_pair(178368193, 178393218);
  transcriptFactorLUT["ZNF460"] = make_pair(57791852, 57805436);
  transcriptFactorLUT["ZNF461"] = make_pair(37128282, 37157755);
  transcriptFactorLUT["ZNF467"] = make_pair(149461452, 149470295);
  transcriptFactorLUT["ZNF468"] = make_pair(53341784, 53360902);
  transcriptFactorLUT["ZNF469"] = make_pair(88493878, 88507165);
  transcriptFactorLUT["ZNF470"] = make_pair(57078889, 57094262);
  transcriptFactorLUT["ZNF471"] = make_pair(57019211, 57040269);
  transcriptFactorLUT["ZNF473"] = make_pair(50529211, 50552031);
  transcriptFactorLUT["ZNF479"] = make_pair(57187325, 57207571);
  transcriptFactorLUT["ZNF48"] = make_pair(30406432, 30411429);
  transcriptFactorLUT["ZNF480"] = make_pair(52800421, 52829180);
  transcriptFactorLUT["ZNF483"] = make_pair(114287438, 114306712);
  transcriptFactorLUT["ZNF484"] = make_pair(95607312, 95640320);
  transcriptFactorLUT["ZNF485"] = make_pair(44101854, 44113352);
  transcriptFactorLUT["ZNF486"] = make_pair(20278022, 20311299);
  transcriptFactorLUT["ZNF490"] = make_pair(12686919, 12721623);
  transcriptFactorLUT["ZNF491"] = make_pair(11909390, 11919306);
  transcriptFactorLUT["ZNF492"] = make_pair(22817125, 22850472);
  transcriptFactorLUT["ZNF493"] = make_pair(21579920, 21610296);
  transcriptFactorLUT["ZNF496"] = make_pair(247463621, 247495045);
  transcriptFactorLUT["ZNF497"] = make_pair(58865722, 58874214);
  transcriptFactorLUT["ZNF500"] = make_pair(4800814, 4817219);
  transcriptFactorLUT["ZNF501"] = make_pair(44771097, 44778575);
  transcriptFactorLUT["ZNF502"] = make_pair(44754134, 44765323);
  transcriptFactorLUT["ZNF506"] = make_pair(19903519, 19932560);
  transcriptFactorLUT["ZNF507"] = make_pair(32836513, 32878573);
  transcriptFactorLUT["ZNF510"] = make_pair(99518146, 99540328);
  transcriptFactorLUT["ZNF512"] = make_pair(27805835, 27846082);
  transcriptFactorLUT["ZNF512B"] = make_pair(62588056, 62601223);
  transcriptFactorLUT["ZNF513"] = make_pair(27600097, 27603311);
  transcriptFactorLUT["ZNF514"] = make_pair(95813399, 95825263);
  transcriptFactorLUT["ZNF516"] = make_pair(74069636, 74207146);
  transcriptFactorLUT["ZNF517"] = make_pair(146024260, 146034529);
  transcriptFactorLUT["ZNF519"] = make_pair(14075988, 14132489);
  transcriptFactorLUT["ZNF521"] = make_pair(22641887, 22932214);
  transcriptFactorLUT["ZNF524"] = make_pair(56111705, 56114504);
  transcriptFactorLUT["ZNF526"] = make_pair(42724491, 42732353);
  transcriptFactorLUT["ZNF527"] = make_pair(37862058, 37883966);
  transcriptFactorLUT["ZNF528"] = make_pair(52901120, 52921657);
  transcriptFactorLUT["ZNF529"] = make_pair(37034516, 37096178);
  transcriptFactorLUT["ZNF530"] = make_pair(58111252, 58119637);
  transcriptFactorLUT["ZNF532"] = make_pair(56530060, 56653709);
  transcriptFactorLUT["ZNF534"] = make_pair(52934666, 52955192);
  transcriptFactorLUT["ZNF536"] = make_pair(30863327, 31048965);
  transcriptFactorLUT["ZNF540"] = make_pair(38042272, 38105079);
  transcriptFactorLUT["ZNF541"] = make_pair(48023941, 48059113);
  transcriptFactorLUT["ZNF543"] = make_pair(57831864, 57842144);
  transcriptFactorLUT["ZNF544"] = make_pair(58740069, 58775008);
  transcriptFactorLUT["ZNF546"] = make_pair(40502942, 40526948);
  transcriptFactorLUT["ZNF547"] = make_pair(57874802, 57890925);
  transcriptFactorLUT["ZNF548"] = make_pair(57901217, 57913919);
  transcriptFactorLUT["ZNF549"] = make_pair(58038692, 58052244);
  transcriptFactorLUT["ZNF550"] = make_pair(58053203, 58071231);
  transcriptFactorLUT["ZNF551"] = make_pair(58193336, 58201169);
  transcriptFactorLUT["ZNF552"] = make_pair(58318449, 58326281);
  transcriptFactorLUT["ZNF554"] = make_pair(2819871, 2836733);
  transcriptFactorLUT["ZNF555"] = make_pair(2841432, 2860472);
  transcriptFactorLUT["ZNF556"] = make_pair(2867332, 2878503);
  transcriptFactorLUT["ZNF557"] = make_pair(7069470, 7087978);
  transcriptFactorLUT["ZNF558"] = make_pair(8920251, 8942980);
  transcriptFactorLUT["ZNF559"] = make_pair(9434447, 9454521);
  transcriptFactorLUT["ZNF560"] = make_pair(9577030, 9609279);
  transcriptFactorLUT["ZNF561"] = make_pair(9718001, 9731916);
  transcriptFactorLUT["ZNF562"] = make_pair(9759337, 9785776);
  transcriptFactorLUT["ZNF563"] = make_pair(12428303, 12444534);
  transcriptFactorLUT["ZNF564"] = make_pair(12636183, 12662356);
  transcriptFactorLUT["ZNF565"] = make_pair(36672961, 36705575);
  transcriptFactorLUT["ZNF566"] = make_pair(36936020, 36980463);
  transcriptFactorLUT["ZNF567"] = make_pair(37180301, 37212238);
  transcriptFactorLUT["ZNF568"] = make_pair(37407230, 37488834);
  transcriptFactorLUT["ZNF569"] = make_pair(37902059, 37958339);
  transcriptFactorLUT["ZNF57"] = make_pair(2900895, 2918474);
  transcriptFactorLUT["ZNF570"] = make_pair(37959959, 37976260);
  transcriptFactorLUT["ZNF571"] = make_pair(38055154, 38085693);
  transcriptFactorLUT["ZNF572"] = make_pair(125985538, 125991630);
  transcriptFactorLUT["ZNF573"] = make_pair(38229202, 38270230);
  transcriptFactorLUT["ZNF574"] = make_pair(42580289, 42585720);
  transcriptFactorLUT["ZNF575"] = make_pair(44037339, 44040284);
  transcriptFactorLUT["ZNF576"] = make_pair(44100543, 44104587);
  transcriptFactorLUT["ZNF577"] = make_pair(52374550, 52391229);
  transcriptFactorLUT["ZNF579"] = make_pair(56088890, 56092211);
  transcriptFactorLUT["ZNF580"] = make_pair(56153462, 56154836);
  transcriptFactorLUT["ZNF581"] = make_pair(56154985, 56156989);
  transcriptFactorLUT["ZNF582"] = make_pair(56894647, 56904889);
  transcriptFactorLUT["ZNF583"] = make_pair(56915717, 56936400);
  transcriptFactorLUT["ZNF584"] = make_pair(58920062, 58929692);
  transcriptFactorLUT["ZNF585A"] = make_pair(37638339, 37663338);
  transcriptFactorLUT["ZNF585B"] = make_pair(37672480, 37701451);
  transcriptFactorLUT["ZNF586"] = make_pair(58281019, 58291984);
  transcriptFactorLUT["ZNF587"] = make_pair(58361180, 58376491);
  transcriptFactorLUT["ZNF589"] = make_pair(48282595, 48312479);
  transcriptFactorLUT["ZNF594"] = make_pair(5082830, 5095178);
  transcriptFactorLUT["ZNF596"] = make_pair(182396, 197340);
  transcriptFactorLUT["ZNF597"] = make_pair(3482421, 3493537);
  transcriptFactorLUT["ZNF599"] = make_pair(35248978, 35264134);
  transcriptFactorLUT["ZNF600"] = make_pair(53268747, 53290034);
  transcriptFactorLUT["ZNF606"] = make_pair(58488440, 58514714);
  transcriptFactorLUT["ZNF607"] = make_pair(38187263, 38210089);
  transcriptFactorLUT["ZNF610"] = make_pair(52839497, 52870376);
  transcriptFactorLUT["ZNF611"] = make_pair(53206065, 53233135);
  transcriptFactorLUT["ZNF613"] = make_pair(52430687, 52449011);
  transcriptFactorLUT["ZNF614"] = make_pair(52516576, 52531680);
  transcriptFactorLUT["ZNF615"] = make_pair(52494586, 52511483);
  transcriptFactorLUT["ZNF616"] = make_pair(52617652, 52643191);
  transcriptFactorLUT["ZNF619"] = make_pair(40518603, 40531728);
  transcriptFactorLUT["ZNF620"] = make_pair(40547601, 40559712);
  transcriptFactorLUT["ZNF621"] = make_pair(40566368, 40581285);
  transcriptFactorLUT["ZNF623"] = make_pair(144718182, 144735900);
  transcriptFactorLUT["ZNF624"] = make_pair(16524047, 16557167);
  transcriptFactorLUT["ZNF625"] = make_pair(12255708, 12267546);
  transcriptFactorLUT["ZNF626"] = make_pair(20827508, 20844402);
  transcriptFactorLUT["ZNF627"] = make_pair(11708234, 11729974);
  transcriptFactorLUT["ZNF628"] = make_pair(55987698, 55995854);
  transcriptFactorLUT["ZNF629"] = make_pair(30789769, 30798523);
  transcriptFactorLUT["ZNF630"] = make_pair(47917566, 47930508);
  transcriptFactorLUT["ZNF639"] = make_pair(179040778, 179053323);
  transcriptFactorLUT["ZNF641"] = make_pair(48733792, 48744674);
  transcriptFactorLUT["ZNF644"] = make_pair(91380856, 91487812);
  transcriptFactorLUT["ZNF646"] = make_pair(31085742, 31094833);
  transcriptFactorLUT["ZNF648"] = make_pair(182023704, 182030847);
  transcriptFactorLUT["ZNF649"] = make_pair(52392487, 52408305);
  transcriptFactorLUT["ZNF652"] = make_pair(47366567, 47439476);
  transcriptFactorLUT["ZNF653"] = make_pair(11594241, 11616738);
  transcriptFactorLUT["ZNF655"] = make_pair(99156447, 99174076);
  transcriptFactorLUT["ZNF658"] = make_pair(40771401, 40792112);
  transcriptFactorLUT["ZNF660"] = make_pair(44626455, 44637557);
  transcriptFactorLUT["ZNF662"] = make_pair(42947657, 42960825);
  transcriptFactorLUT["ZNF664"] = make_pair(124457761, 124499986);
  transcriptFactorLUT["ZNF665"] = make_pair(53666551, 53696619);
  transcriptFactorLUT["ZNF667"] = make_pair(56950692, 56988770);
  transcriptFactorLUT["ZNF668"] = make_pair(31072163, 31076409);
  transcriptFactorLUT["ZNF669"] = make_pair(247263263, 247267674);
  transcriptFactorLUT["ZNF670"] = make_pair(247197939, 247242115);
  transcriptFactorLUT["ZNF671"] = make_pair(58231118, 58238995);
  transcriptFactorLUT["ZNF672"] = make_pair(249132376, 249143716);
  transcriptFactorLUT["ZNF674"] = make_pair(46357159, 46404892);
  transcriptFactorLUT["ZNF675"] = make_pair(23835707, 23870017);
  transcriptFactorLUT["ZNF676"] = make_pair(22361902, 22379753);
  transcriptFactorLUT["ZNF677"] = make_pair(53738637, 53758111);
  transcriptFactorLUT["ZNF678"] = make_pair(227751219, 227850164);
  transcriptFactorLUT["ZNF679"] = make_pair(63688851, 63727309);
  transcriptFactorLUT["ZNF680"] = make_pair(63985033, 64023505);
  transcriptFactorLUT["ZNF681"] = make_pair(23921996, 23941693);
  transcriptFactorLUT["ZNF682"] = make_pair(20115226, 20150277);
  transcriptFactorLUT["ZNF683"] = make_pair(26688124, 26699266);
  transcriptFactorLUT["ZNF684"] = make_pair(40997232, 41013841);
  transcriptFactorLUT["ZNF687"] = make_pair(151254787, 151264381);
  transcriptFactorLUT["ZNF688"] = make_pair(30581018, 30583728);
  transcriptFactorLUT["ZNF689"] = make_pair(30613878, 30621754);
  transcriptFactorLUT["ZNF69"] = make_pair(11998669, 12025365);
  transcriptFactorLUT["ZNF691"] = make_pair(43312243, 43318146);
  transcriptFactorLUT["ZNF692"] = make_pair(249144202, 249153315);
  transcriptFactorLUT["ZNF696"] = make_pair(144373558, 144382120);
  transcriptFactorLUT["ZNF697"] = make_pair(120161999, 120190390);
  transcriptFactorLUT["ZNF699"] = make_pair(9405985, 9415795);
  transcriptFactorLUT["ZNF7"] = make_pair(146052902, 146068607);
  transcriptFactorLUT["ZNF70"] = make_pair(24083770, 24093279);
  transcriptFactorLUT["ZNF700"] = make_pair(12035882, 12061588);
  transcriptFactorLUT["ZNF701"] = make_pair(53073525, 53090427);
  transcriptFactorLUT["ZNF705A"] = make_pair(8325149, 8332642);
  transcriptFactorLUT["ZNF705D"] = make_pair(11946846, 11973025);
  transcriptFactorLUT["ZNF705G"] = make_pair(7215497, 7220490);
  transcriptFactorLUT["ZNF707"] = make_pair(144766621, 144777555);
  transcriptFactorLUT["ZNF708"] = make_pair(21473962, 21512212);
  transcriptFactorLUT["ZNF709"] = make_pair(12571997, 12595632);
  transcriptFactorLUT["ZNF71"] = make_pair(57106663, 57135544);
  transcriptFactorLUT["ZNF710"] = make_pair(90544722, 90625432);
  transcriptFactorLUT["ZNF711"] = make_pair(84498996, 84528368);
  transcriptFactorLUT["ZNF713"] = make_pair(55954969, 56009918);
  transcriptFactorLUT["ZNF714"] = make_pair(21264952, 21307883);
  transcriptFactorLUT["ZNF716"] = make_pair(57509882, 57533265);
  transcriptFactorLUT["ZNF717"] = make_pair(75786028, 75834734);
  transcriptFactorLUT["ZNF732"] = make_pair(264463, 289944);
  transcriptFactorLUT["ZNF736"] = make_pair(63774250, 63817012);
  transcriptFactorLUT["ZNF737"] = make_pair(20720797, 20748626);
  transcriptFactorLUT["ZNF74"] = make_pair(20748404, 20762753);
  transcriptFactorLUT["ZNF740"] = make_pair(53574534, 53584654);
  transcriptFactorLUT["ZNF746"] = make_pair(149169883, 149194898);
  transcriptFactorLUT["ZNF749"] = make_pair(57946692, 57957191);
  transcriptFactorLUT["ZNF75A"] = make_pair(3355405, 3368576);
  transcriptFactorLUT["ZNF75D"] = make_pair(134419718, 134478012);
  transcriptFactorLUT["ZNF76"] = make_pair(35227490, 35263764);
  transcriptFactorLUT["ZNF763"] = make_pair(12075868, 12091198);
  transcriptFactorLUT["ZNF764"] = make_pair(30565084, 30569642);
  transcriptFactorLUT["ZNF765"] = make_pair(53898396, 53915262);
  transcriptFactorLUT["ZNF766"] = make_pair(52772823, 52795976);
  transcriptFactorLUT["ZNF768"] = make_pair(30535321, 30537910);
  transcriptFactorLUT["ZNF77"] = make_pair(2933215, 2944969);
  transcriptFactorLUT["ZNF770"] = make_pair(35270541, 35280497);
  transcriptFactorLUT["ZNF771"] = make_pair(30418734, 30429916);
  transcriptFactorLUT["ZNF772"] = make_pair(57980953, 57988938);
  transcriptFactorLUT["ZNF773"] = make_pair(58011221, 58024519);
  transcriptFactorLUT["ZNF774"] = make_pair(90895476, 90904715);
  transcriptFactorLUT["ZNF775"] = make_pair(150076405, 150095719);
  transcriptFactorLUT["ZNF776"] = make_pair(58258163, 58269527);
  transcriptFactorLUT["ZNF777"] = make_pair(149128453, 149158053);
  transcriptFactorLUT["ZNF778"] = make_pair(89284110, 89295965);
  transcriptFactorLUT["ZNF780A"] = make_pair(40578898, 40596845);
  transcriptFactorLUT["ZNF780B"] = make_pair(40534166, 40562115);
  transcriptFactorLUT["ZNF781"] = make_pair(38158649, 38183216);
  transcriptFactorLUT["ZNF782"] = make_pair(99579272, 99616389);
  transcriptFactorLUT["ZNF783"] = make_pair(148959261, 148982085);
  transcriptFactorLUT["ZNF784"] = make_pair(56132106, 56135941);
  transcriptFactorLUT["ZNF785"] = make_pair(30591993, 30597092);
  transcriptFactorLUT["ZNF786"] = make_pair(148766732, 148787869);
  transcriptFactorLUT["ZNF787"] = make_pair(56598728, 56632742);
  transcriptFactorLUT["ZNF789"] = make_pair(99070514, 99079948);
  transcriptFactorLUT["ZNF79"] = make_pair(130186652, 130207651);
  transcriptFactorLUT["ZNF790"] = make_pair(37309223, 37341689);
  transcriptFactorLUT["ZNF791"] = make_pair(12721731, 12740676);
  transcriptFactorLUT["ZNF792"] = make_pair(35447257, 35454953);
  transcriptFactorLUT["ZNF793"] = make_pair(37997840, 38034239);
  transcriptFactorLUT["ZNF799"] = make_pair(12500827, 12512088);
  transcriptFactorLUT["ZNF8"] = make_pair(58790317, 58807254);
  transcriptFactorLUT["ZNF80"] = make_pair(113953479, 113956425);
  transcriptFactorLUT["ZNF800"] = make_pair(127010096, 127032778);
  transcriptFactorLUT["ZNF805"] = make_pair(57752052, 57774106);
  transcriptFactorLUT["ZNF808"] = make_pair(53030908, 53059303);
  transcriptFactorLUT["ZNF81"] = make_pair(47696300, 47781655);
  transcriptFactorLUT["ZNF813"] = make_pair(53970988, 53997546);
  transcriptFactorLUT["ZNF814"] = make_pair(58380746, 58400442);
  transcriptFactorLUT["ZNF816"] = make_pair(53452631, 53466164);
  transcriptFactorLUT["ZNF821"] = make_pair(71893582, 71917444);
  transcriptFactorLUT["ZNF823"] = make_pair(11832079, 11849824);
  transcriptFactorLUT["ZNF827"] = make_pair(146681887, 146859607);
  transcriptFactorLUT["ZNF83"] = make_pair(53115617, 53141644);
  transcriptFactorLUT["ZNF831"] = make_pair(57766074, 57834167);
  transcriptFactorLUT["ZNF836"] = make_pair(52658124, 52674896);
  transcriptFactorLUT["ZNF837"] = make_pair(58878989, 58892389);
  transcriptFactorLUT["ZNF84"] = make_pair(133614167, 133639890);
  transcriptFactorLUT["ZNF841"] = make_pair(52567718, 52599018);
  transcriptFactorLUT["ZNF844"] = make_pair(12175545, 12188626);
  transcriptFactorLUT["ZNF845"] = make_pair(53837001, 53858122);
  transcriptFactorLUT["ZNF846"] = make_pair(9868150, 9879410);
  transcriptFactorLUT["ZNF85"] = make_pair(21106058, 21133503);
  transcriptFactorLUT["ZNF852"] = make_pair(44540461, 44552132);
  transcriptFactorLUT["ZNF860"] = make_pair(32023265, 32033228);
  transcriptFactorLUT["ZNF865"] = make_pair(56124958, 56129907);
  transcriptFactorLUT["ZNF878"] = make_pair(12154619, 12163782);
  transcriptFactorLUT["ZNF879"] = make_pair(178450775, 178461388);
  transcriptFactorLUT["ZNF880"] = make_pair(52873169, 52889046);
  transcriptFactorLUT["ZNF90"] = make_pair(20188802, 20231977);
  transcriptFactorLUT["ZNF91"] = make_pair(23540497, 23578362);
  transcriptFactorLUT["ZNF92"] = make_pair(64838711, 64866048);
  transcriptFactorLUT["ZNF93"] = make_pair(20011721, 20046382);
  transcriptFactorLUT["ZNF98"] = make_pair(22573898, 22605148);
  transcriptFactorLUT["ZNF99"] = make_pair(22934984, 22966973);
  transcriptFactorLUT["ZSCAN1"] = make_pair(58545433, 58565999);
  transcriptFactorLUT["ZSCAN10"] = make_pair(3138890, 3149318);
  transcriptFactorLUT["ZSCAN12"] = make_pair(28356726, 28367544);
  transcriptFactorLUT["ZSCAN16"] = make_pair(28092386, 28097856);
  transcriptFactorLUT["ZSCAN18"] = make_pair(58595208, 58609730);
  transcriptFactorLUT["ZSCAN2"] = make_pair(85144248, 85159841);
  transcriptFactorLUT["ZSCAN20"] = make_pair(33938231, 33961995);
  transcriptFactorLUT["ZSCAN21"] = make_pair(99647416, 99662663);
  transcriptFactorLUT["ZSCAN22"] = make_pair(58838384, 58853712);
  transcriptFactorLUT["ZSCAN23"] = make_pair(28400431, 28411279);
  transcriptFactorLUT["ZSCAN29"] = make_pair(43650369, 43662258);
  transcriptFactorLUT["ZSCAN30"] = make_pair(32831021, 32870209);
  transcriptFactorLUT["ZSCAN4"] = make_pair(58180302, 58190520);
  transcriptFactorLUT["ZSCAN5A"] = make_pair(56732678, 56739659);
  transcriptFactorLUT["ZSCAN5B"] = make_pair(56701057, 56704421);
  transcriptFactorLUT["ZXDA"] = make_pair(57931863, 57937067);
  transcriptFactorLUT["ZXDB"] = make_pair(57618268, 57623910);
  transcriptFactorLUT["ZXDC"] = make_pair(126177743, 126194762);
  transcriptFactorLUT["ZZZ3"] = make_pair(78030189, 78148343);

  return true;
}
