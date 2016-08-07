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

using namespace std;

PlinkInternalsDataset::PlinkInternalsDataset(Plink* plinkPtr): 
  Dataset::Dataset() {
  PP = plinkPtr;
}

bool PlinkInternalsDataset::LoadDataset() {
  PP->printLOG("Adapting PLINK data structures to Dataset object\n");
  // --------------------------------------------------------------------------
  // get variable metadata
  PP->printLOG("Get PLINK variable metadata\n");
  unsigned int numAttributes = PP->nl_all;
  if(numAttributes) {
    hasGenotypes = true;
    attributeAlleles.resize(numAttributes);
    attributeMinorAllele.resize(numAttributes);
    attributeMutationTypes.resize(numAttributes);
    for(int i=0; i < numAttributes; i++) {
      Locus* locus = PP->locus[i];
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
  unsigned int numNumerics = PP->nlistname.size();
  if(numNumerics) {
    hasNumerics = true;
    for(int i=0; i < PP->nlistname.size(); i++) {
      string numericName = PP->nlistname[i];
      numericsNames.push_back(numericName);
      numericsMask[numericName] = i;
    }
  }

  // --------------------------------------------------------------------------
  // look at all SNP and numeric values for all subjects
  PP->printLOG("Get PLINK variable data for all subject\n");
	attributeLevelsSeen.resize(numAttributes);
	genotypeCounts.resize(numAttributes);
	attributeAlleleCounts.resize(numAttributes);
	attributeMutationTypes.resize(numAttributes);
  for(int i=0; i < PP->sample.size(); i++) {
    string ID = PP->sample[i]->fid + PP->sample[i]->iid;
    instanceIds.push_back(ID);
    if(PP->sample[i]->phenotype == -9) {
      // missing phenotype so skip this individual
      PP->printLOG("Missing phenotype for ID: " + ID + ", skipping\n");
      continue;
    }
    PlinkInternalsDatasetInstance* tmpInd = 0;
    tmpInd = new PlinkInternalsDatasetInstance(this, ID, PP, PP->sample[i]);
    double phenotype = -9;
    if(par::bt) {
      hasContinuousPhenotypes = true;
      phenotype = PP->sample[i]->phenotype;
      tmpInd->SetPredictedValueTau(phenotype);
      instanceIdsToLoad.push_back(ID);
    } else {
      hasContinuousPhenotypes = false;
      int intPheno = PP->sample[i]->aff? 2: 1;
      phenotype = static_cast<double>(PP->sample[i]->aff? 2: 1);
      classIndexes[intPheno].push_back(i);
      tmpInd->SetClass(phenotype);
    }
    instances.push_back(tmpInd);
    instancesMask[ID] = i;
    if(hasGenotypes) {
      for(int j=0; j < numAttributes; j++) {
        AttributeLevel attr;
        attr = static_cast<AttributeLevel>(tmpInd->GetSimpleSNPValue(j));
        if(attr == -9) {
          missingValues[ID].push_back(j);
        }
        ++levelCounts[j][attr];
        ++levelCountsByClass[j][make_pair(attr, tmpInd->GetClass())];
        pair<char, char> alleles = attributeAlleles[j];
        ++attributeAlleleCounts[j][alleles.first];
        ++attributeAlleleCounts[j][alleles.second];
        string genotype = "" + alleles.first + alleles.second;
        ++genotypeCounts[j][genotype];
        attributeLevelsSeen[j].insert(genotype);
      }
      MaskIncludeAllAttributes(DISCRETE_TYPE);
    }
    if(hasNumerics) {
      numericsSums.resize(numNumerics);
      numericsMinMax.resize(numNumerics);
      for(int j=0; j < numNumerics; j++) {
        NumericLevel numeric = PP->sample[i]->nlist[j];
        if(numeric == -9) {
          missingNumericValues[ID].push_back(j);
        }
        numericsSums[j] += numeric;
        if(numeric < numericsMinMax[j].first) {
          numericsMinMax[j].first = numeric;
        }
        if(numeric > numericsMinMax[j].second) {
          numericsMinMax[j].second = numeric;
        }
      }
      MaskIncludeAllAttributes(NUMERIC_TYPE);
    }
  }
  hasPhenotypes = true;  
  hasAllelicInfo = true;

  // this has already been done above; called by iterative Relief-F to update
  // UpdateAllLevelCounts();

  return true;
}
