/*
 * PlinkInternalsDatasetInstance.cpp - Bill White
 * 
 * Class to hold dataset instances (rows) - 6/14/05
 * Reworked entirely for McKinney Lab work - 2/28/11
 * Modified for inclusion in inbix - 8/3/16
*/

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "plink.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "PlinkInternalsDatasetInstance.h"
#include "StringUtils.h"
#include "BestN.h"
#include "Insilico.h"

using namespace std;
using namespace insilico;

PlinkInternalsDatasetInstance::
PlinkInternalsDatasetInstance(Dataset* ds, string instanceID, 
        Plink* plinkPtr, Individual* plinkInd): DatasetInstance(ds) {
  PP = plinkPtr;
  ID = instanceID;
  individual = plinkInd;
  // switch from SNP-major to individual-major data orientation!
  // PP->SNP2Ind();
}

PlinkInternalsDatasetInstance::~PlinkInternalsDatasetInstance() {
}

double PlinkInternalsDatasetInstance::GetSimpleSNPValue(int snp) {
  bool i1 = individual->one[snp];
	bool i2 = individual->two[snp];
  ///////////////////////
  // Autosomal coding: counts minor alleles: aka dosage
  double retVal = 0;
  if(i1) {
    if(!i2) {
      retVal = 0;
    } else {
      retVal = 0;
    }
  } else {
    if(i2) {
      retVal = 1; // het 
    }
    else {
      retVal = 2; // hom
    }
  }
  return retVal;
}

unsigned int PlinkInternalsDatasetInstance::NumAttributes() {
  return PP->nl_all;
}

AttributeLevel PlinkInternalsDatasetInstance::GetAttribute(unsigned int index) {
  if(dataset->HasGenotypes()) {
    if(index < dataset->NumAttributes()) {
      return static_cast<AttributeLevel>(GetSimpleSNPValue(index));
    } else {
      cerr << "ERROR: Attribute index is out of range: " << index << endl;
      exit(1);
    }
  } else {
    cerr << "ERROR: Attempting to access SNP value when none "
            << "have been loaded" << endl;
    exit(1);
  }
}

unsigned int PlinkInternalsDatasetInstance::NumNumerics() {
  return(PP->nlistname.size());
}

double PlinkInternalsDatasetInstance::GetNumeric(unsigned int index) {
  if(dataset->HasNumerics()) {
    if(index < dataset->NumNumerics()) {
      return individual->nlist[index];
    } else {
      cerr << "ERROR: Numeric index out of range: " << index << endl;
      exit(1);
    }
  } else {
    cerr << "ERROR: Attempting to access numeric value when none "
            << "have been loaded" << endl;
    exit(1);
  }
}

void PlinkInternalsDatasetInstance::Print() {
  cout << "Instance ID [" << ID << "]" << endl;
  cout << "PLINK pointer [" << PP << "]" << endl;
  cout << "PLINK individual pointer [" << individual << endl;
  // for RReliefF - bcw - 9/30/11
  if(dataset->HasContinuousPhenotypes()) {
    cout << "=> [" << predictedValueTau << "]" << endl;
  } else {
    cout << "=> [" << classLabel << "]" << endl;
  }
}
