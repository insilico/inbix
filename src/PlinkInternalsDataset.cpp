/*
 * PlinkInternalsDataset.cpp - Bill White - 2/24/11
 *
 * Collection class holding PlinkInternalDatasetInstance from 
 * Plink internal data structures accessed through a pointer.
 * 
 * Modified for inclusion in inbix 8/3/16
 */

#include "plink.h"
#include "options.h"
#include "PlinkInternalsDataset.h"
#include "PlinkInternalsDatasetInstance.h"
#include "helper.h"

using namespace std;

PlinkInternalsDataset::PlinkInternalsDataset(Plink* plinkPtr): 
  Dataset::Dataset() {
  plinkInternalsPtr = plinkPtr;
}

PlinkInternalsDataset::~PlinkInternalsDataset() {
}

bool PlinkInternalsDataset::LoadDatasetFromPlink() {
  plinkInternalsPtr->printLOG(Timestamp() + "Adapting PLINK data structures to Dataset object\n");
  // --------------------------------------------------------------------------
  // get variable metadata
  unsigned int numAttributes = plinkInternalsPtr->nl_all;
  if(numAttributes) {
    plinkInternalsPtr->printLOG(Timestamp() + "Get PLINK SNP variable metadata for " + 
                 int2str(numAttributes) + " SNPs\n");
    hasGenotypes = true;
    attributeAlleles.resize(numAttributes);
    attributeMinorAllele.resize(numAttributes);
    attributeMutationTypes.resize(numAttributes);
    for(int i=0; i < numAttributes; i++) {
      Locus* locus = plinkInternalsPtr->locus[i];
      string locusName = locus->name;
      attributeNames.push_back(locusName);
      attributesMask[locusName] = i;
      char a1 = locus->allele1[0];
      char a2 = locus->allele2[0];
      pair<char, char> alleles = make_pair(a1, a2);
      attributeAlleles[i] = alleles;
      attributeMinorAllele[i] = make_pair(a1, locus->freq);
      attributeMutationTypes[i] = attributeMutationMap[alleles];
    }
  }
  unsigned int numNumerics = plinkInternalsPtr->nlistname.size();
  if(numNumerics) {
    plinkInternalsPtr->printLOG(Timestamp() + "Get PLINK Numeric variable metadata for " + 
                 int2str(numNumerics) + " numerics\n");
    hasNumerics = true;
    for(int i=0; i < plinkInternalsPtr->nlistname.size(); i++) {
      string numericName = plinkInternalsPtr->nlistname[i];
      numericsNames.push_back(numericName);
      numericsMask[numericName] = i;
    }
  }
  
  if((numAttributes + numNumerics) < 1) {
    error("No SNPs or numeric attributes found");
  }
  // --------------------------------------------------------------------------
  // look at all SNP and numeric values for all subjects
  plinkInternalsPtr->printLOG(Timestamp() + "Get PLINK variable data for all subjects\n");
  
	attributeLevelsSeen.resize(numAttributes);
	genotypeCounts.resize(numAttributes);
	attributeAlleleCounts.resize(numAttributes);
	attributeMutationTypes.resize(numAttributes);
  levelCounts.resize(numAttributes);

  // binary or continuous trait/phenotype
  if(par::bt) {
    plinkInternalsPtr->printLOG(Timestamp() + "Detected BINARY phenotype\n");
    hasContinuousPhenotypes = false;
    levelCountsByClass.resize(numAttributes);
  } else {
    plinkInternalsPtr->printLOG(Timestamp() + "Detected CONTINUOUS phenotype\n");
    hasContinuousPhenotypes = true;
  }
  hasPhenotypes = true;  
  
  unsigned int nextInstanceIdx = 0;
  instances.clear();
  for(int sampleIdx=0; sampleIdx < plinkInternalsPtr->sample.size(); sampleIdx++) {
    string ID = plinkInternalsPtr->sample[sampleIdx]->fid + plinkInternalsPtr->sample[sampleIdx]->iid;
//    if(par::verbose) {
//      PP->printLOG("[ DEBUG ] individual i: " + int2str(i)
//              + ", ID: " + ID 
//              + ", next index: " + int2str(nextInstanceIdx) + "\n");
//    }
    if(plinkInternalsPtr->sample[sampleIdx]->missing) {
      // missing phenotype so skip this individual
      if(par::verbose) {
        plinkInternalsPtr->printLOG(Timestamp() + "ID: " + ID + " missing, skipping\n");
      }
      continue;
    }
    // this individual is a "go"!
    PlinkInternalsDatasetInstance* tmpInd = 
      new PlinkInternalsDatasetInstance(this, ID, plinkInternalsPtr, plinkInternalsPtr->sample[sampleIdx]);
    double phenotype = -9;
    if(par::bt) {
      // encode 0/1 for my Dataset class assumptions
      uint intPheno = plinkInternalsPtr->sample[sampleIdx]->aff? 1: 0;
      phenotype = static_cast<double>(intPheno);
      classIndexes[intPheno].push_back(sampleIdx);
      tmpInd->SetClass(phenotype);
    } else {
      phenotype = plinkInternalsPtr->sample[sampleIdx]->phenotype;
      if(std::isnan(phenotype) || std::isinf(phenotype)) {
        error("Invalid continuous phenotype value: [ " + dbl2str(phenotype) + " ]");
      }
      if(sampleIdx == 0) {
        continuousPhenotypeMinMax.first = phenotype;
        continuousPhenotypeMinMax.second = phenotype;
      } else {
        if(phenotype < continuousPhenotypeMinMax.first) {
          continuousPhenotypeMinMax.first = phenotype;
        }
        if(phenotype > continuousPhenotypeMinMax.second) {
          continuousPhenotypeMinMax.second = phenotype;
        }
      }
      tmpInd->SetPredictedValueTau(phenotype);
    }
    instanceIds.push_back(ID);
    instanceIdsToLoad.push_back(ID);
    instancesMask[ID] = nextInstanceIdx;
    nextInstanceIdx++;
    vector<AttributeLevel> tmpSnps;
    if(hasGenotypes) {
      for(int attrIdx=0; attrIdx < numAttributes; attrIdx++) {
        AttributeLevel attr = -9;
        attr = static_cast<AttributeLevel>(tmpInd->GetSimpleSNPValue(attrIdx));
        tmpSnps.push_back(attr);
        if(attr == -9) {
          missingValues[ID].push_back(attrIdx);
        }
        ++levelCounts[attrIdx][attr];
        if(!HasContinuousPhenotypes()) {
          ClassLevel classLevel = static_cast<ClassLevel>(tmpInd->GetClass());
          ++levelCountsByClass[attrIdx][make_pair(attr, classLevel)];
        }
        pair<char, char> alleles = attributeAlleles[attrIdx];
        ++attributeAlleleCounts[attrIdx][alleles.first];
        ++attributeAlleleCounts[attrIdx][alleles.second];
        string genotype = "  ";
        genotype[0] = alleles.first;
        genotype[1] = alleles.second;
        ++genotypeCounts[attrIdx][genotype];
        attributeLevelsSeen[attrIdx].insert(genotype);
      }
      tmpInd->LoadInstanceFromVector(tmpSnps);
      MaskIncludeAllAttributes(DISCRETE_TYPE);
    }
    if(hasNumerics) {
      numericsSums.resize(numNumerics);
      numericsMinMax.resize(numNumerics);
      for(int numericIdx=0; numericIdx < numNumerics; numericIdx++) {
        NumericLevel numericValue = plinkInternalsPtr->sample[sampleIdx]->nlist[numericIdx];
        if(numericValue == -9) {
          missingNumericValues[ID].push_back(numericIdx);
          continue;
        }
        numericsSums[numericIdx] += numericValue;
        if(numericValue < numericsMinMax[numericIdx].first) {
          numericsMinMax[numericIdx].first = numericValue;
        }
        if(numericValue > numericsMinMax[numericIdx].second) {
          numericsMinMax[numericIdx].second = numericValue;
        }
        tmpInd->AddNumeric(numericValue);
      }
      MaskIncludeAllAttributes(NUMERIC_TYPE);
    }
    instances.push_back(tmpInd);
    instancesMask[ID] = sampleIdx;
  }

  if(hasGenotypes) {
    hasAllelicInfo = true;
  } else {
    hasAllelicInfo = false;
  }

  plinkInternalsPtr->printLOG(Timestamp() + "Final number of subjects loaded: " + int2str(instances.size()) + "\n");
  plinkInternalsPtr->printLOG(Timestamp() + "Final number of SNPs loaded: " + int2str(numAttributes) + "\n");
  plinkInternalsPtr->printLOG(Timestamp() + "Final number of numerics loaded: " + int2str(numNumerics) + "\n");

  return true;
}
