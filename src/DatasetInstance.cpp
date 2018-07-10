/*
 * DatasetInstance.cpp - Bill White - 6/14/05
 * 
 * Class to hold dataset instances (rows)
 * Reworked entirely for McKinney Lab work - 2/28/11
 */

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "Dataset.h"
#include "DatasetInstance.h"
#include "StringUtils.h"
#include "BestN.h"
#include "helper.h"

using namespace std;
using namespace insilico;

/// functor for T comparison
typedef DistancePair T;

class deref_less_bcw : public std::binary_function<T, T, bool>
{
public:

  bool operator()(const T a, const T b) const {
    return(a.first < b.first);
  }
};

DatasetInstance::DatasetInstance(Dataset* ds, string newId) {
  dataset = ds;
  id = newId;
  classLabel = MISSING_DISCRETE_CLASS_VALUE;
  predictedValueTau = MISSING_DISCRETE_CLASS_VALUE;
  cvSetType = CV_NONE;
}

DatasetInstance::~DatasetInstance() {
}

Dataset* DatasetInstance::GetDatasetPtr() {
  return dataset;
}

bool
DatasetInstance::LoadInstanceFromVector(vector<AttributeLevel> newAttributes) {
  if(!newAttributes.size()) {
  	cerr << "ERROR: LoadInstanceFromVector: vector is empty" << endl;
    return false;
  }
  attributes.clear();
  vector<AttributeLevel>::const_iterator it;
  for(it = newAttributes.begin(); it != newAttributes.end(); it++) {
    attributes.push_back(*it);
  }
  return true;
}

bool 
DatasetInstance::LoadInstanceFromInstancePtr(Dataset* srcDs, 
                                             DatasetInstance* srcInstance) {
  dataset = srcDs;
  classLabel = srcInstance->GetClass();
  attributes = srcInstance->GetAttributes();
  numerics = srcInstance->GetNumerics();
  predictedValueTau = srcInstance->GetPredictedValueTau();
  id = srcInstance->GetId();
  
  return true;
}

bool DatasetInstance::SetId(std::string newId) {
  id = newId;
  return true;
}

unsigned int DatasetInstance::NumAttributes() {
  return(attributes.size());
}

AttributeLevel DatasetInstance::GetAttribute(unsigned int index) {
  if(attributes.size()) {
    if(index < attributes.size()) {
      return attributes[index];
    } else {
      cerr << "ERROR: Attribute index is out of range: " << index << endl;
      exit(1);
    }
  } else {
    cerr << "ERROR: Attempting to access attribute value when none "
            << "have been loaded" << endl;
    exit(1);
  }

}

unsigned int DatasetInstance::NumNumerics() {
  return(numerics.size());
}

double DatasetInstance::GetNumeric(unsigned int index) {
  if(numerics.size()) {
    if(index < numerics.size()) {
      return numerics[index];
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

bool DatasetInstance::AddNumeric(NumericLevel newNum) {
  //  cout << "DatasetInstance::AddNumeric(double newNum): " << newNum << endl;
  numerics.push_back(newNum);
  return true;
}

ClassLevel DatasetInstance::GetClass() {
  return classLabel;
}

void DatasetInstance::SetClass(ClassLevel classValue) {
  classLabel = classValue;
}

double DatasetInstance::GetPredictedValueTau() {
  return predictedValueTau;
}

void DatasetInstance::SetPredictedValueTau(double newValue) {
  predictedValueTau = newValue;
}

double DatasetInstance::GetInfluenceFactorD(unsigned int neighborIndex) {
  return neighborInfluenceFactorDs[neighborIndex];
}

void DatasetInstance::ClearInfluenceFactors() {
  neighborInfluenceFactorDs.clear();
}

bool DatasetInstance::AddInfluenceFactorD(double factor) {
  neighborInfluenceFactorDs.push_back(factor);
  return true;
}

void DatasetInstance::Print() {
  vector<AttributeLevel>::const_iterator it = attributes.begin();
  for(; it != attributes.end(); ++it) {
    cout << *it << " ";
  }
  if(numerics.size()) {
    cout << " | ";
    vector<double>::const_iterator dit = numerics.begin();
    for(; dit != numerics.end(); ++dit) {
      cout << *dit << " ";
    }
  }
  // for RReliefF - bcw - 9/30/11
  if(dataset->HasContinuousPhenotypes()) {
    cout << "=> [" << predictedValueTau << "]" << endl;
  } else {
    cout << "=> [" << classLabel << "]" << endl;
  }
}

bool DatasetInstance::SwapAttributes(unsigned int a1, unsigned int a2) {
  if(a1 >= attributes.size()) {
    return false;
  }
  if(a2 >= attributes.size()) {
    return false;
  }
  // hahaha
  if(a1 == a2) {
    return true;
  }

  AttributeLevel temp = attributes[a1];
  attributes[a1] = attributes[a2];
  attributes[a2] = temp;

  return true;
}

void DatasetInstance::SetDistanceSums(unsigned int kNearestNeighbors,
                                      DistancePairs& sameClassSums,
                                      map<ClassLevel, DistancePairs>& diffClassSums) {
  // added 9/22/11 for iterative Relief-F
  bestNeighborIdsSameClass.clear();
  bestNeighborIdsDiffClass.clear();
  
  // use Nate's best_n.h algorithm
  //PrintDistancePairs(sameClassSums);
  DistancePairs bestInstancesHits;
  best_n(sameClassSums.begin(), sameClassSums.end(),
         back_insert_iterator<DistancePairs>(bestInstancesHits),
         kNearestNeighbors, deref_less_bcw());
  DistancePairsIt hit;
  for(hit = bestInstancesHits.begin(); hit != bestInstancesHits.end(); ++hit) {
    DistancePair thisHit = *hit;
    bestNeighborIdsSameClass.push_back(thisHit.second);
  }
  map<ClassLevel, DistancePairs>::const_iterator it = diffClassSums.begin();
  for(; it != diffClassSums.end(); ++it) {
    ClassLevel thisClass = it->first;
    DistancePairs thisDiffSums = it->second;
    DistancePairs bestInstancesMisses;
    best_n(thisDiffSums.begin(), thisDiffSums.end(),
           back_insert_iterator<DistancePairs > (bestInstancesMisses),
           kNearestNeighbors, deref_less_bcw());
    DistancePairsIt mit;
    //PrintDistancePairs(bestInstancesMisses);
    for(mit = bestInstancesMisses.begin(); mit != bestInstancesMisses.end(); ++mit) {
      DistancePair thisMiss = *mit;
      bestNeighborIdsDiffClass[thisClass].push_back(thisMiss.second);
    }
  }
}

void DatasetInstance::SetDistanceSums(unsigned int kNearestNeighbors,
                                      DistancePairs instanceSums) {
  bestNeighborIds.clear();

  // cout << "Instance sums:" << endl;
  // PrintDistancePairs(instanceSums);

  // use Nate's best_n.h algorithm to select the nearest neighbors
  DistancePairs bestInstances;
  best_n(instanceSums.begin(), instanceSums.end(),
         back_insert_iterator<DistancePairs > (bestInstances),
         kNearestNeighbors, deref_less_bcw());
  sort(bestInstances.begin(), bestInstances.end());
  // cout << "Best instances:" << endl;
  // PrintDistancePairs(bestInstances);

  DistancePairsIt it = bestInstances.begin();
  for(; it != bestInstances.end(); ++it) {
    //    cout << it->first << " => " << it->second << endl;
    bestNeighborIds.push_back(it->second);
  }
  //  PrintVector(bestNeighborIds, "Best neighbor IDs");
  //  cout << "----------------------------------------------------------" << endl;
}

void DatasetInstance::PrintDistancePairs(const DistancePairs& distPairs) {
  for(DistancePairsIt dpit = distPairs.begin(); dpit != distPairs.end(); ++dpit) {
    cout << (*dpit).first << "\t" << (*dpit).second << endl;
  }
}

bool DatasetInstance::GetNNearestInstances(unsigned int n,
                                           vector<unsigned int>& sameClassInstances,
                                           vector<unsigned int>& diffClassInstances) {
  if((bestNeighborIdsSameClass.size() < n) || (bestNeighborIdsDiffClass.size() < n)) {
    cerr << endl << "ERROR: GetNNearestInstances: N: [" << n
            << "] is larger than the number of neighbors: "
            << "Same: " << bestNeighborIdsSameClass.size()
            << ", Different: " << bestNeighborIdsDiffClass.size() << endl;
    return false;
  }

  sameClassInstances.clear();
  diffClassInstances.clear();
  for(unsigned int i = 0; i < n; ++i) {
    unsigned int sameIdx;
    dataset->GetInstanceIndexForID(bestNeighborIdsSameClass[i], sameIdx);
    sameClassInstances.push_back(sameIdx);
    unsigned int diffIdx;
    dataset->GetInstanceIndexForID(bestNeighborIdsSameClass[i], diffIdx);
    sameClassInstances.push_back(diffIdx);
  }

  return true;
}

bool DatasetInstance::GetNNearestInstances(unsigned int n,
        vector<unsigned int>& sameClassInstances,
        map<ClassLevel, vector<unsigned int> >& diffClassInstances) {

  if(bestNeighborIdsSameClass.size() < n) {
    error("GetNNearestInstances: N: [" + int2str(n) + 
          "] is larger than the number of neighbors in same class: [" + 
          int2str(bestNeighborIdsSameClass.size()) +"]\n");
    return false;
  }
  map<string, uint> instMap = dataset->MaskGetInstanceMask();
  sameClassInstances.clear();
  for(uint i = 0; i < n; ++i) {
    uint sameIdx = instMap[bestNeighborIdsSameClass[i]];
    //dataset->GetInstanceIndexForID(bestNeighborIdsSameClass[i], sameIdx);
    sameClassInstances.push_back(sameIdx);
  }
  diffClassInstances.clear();
  map<ClassLevel, vector<string> >::const_iterator it;
  for(it = bestNeighborIdsDiffClass.begin();
      it != bestNeighborIdsDiffClass.end(); ++it) {
    ClassLevel thisClass = it->first;
    vector<string> thisClassIds = it->second;
    if(thisClassIds.size() < n) {
      error("GetNNearestInstances: N: [" + int2str(n) + 
             "] is larger than the number of neighbors for class " +
             int2str(thisClass) + ": [" + 
             int2str(bestNeighborIdsDiffClass.size()) + "]\n");
    }
    for(unsigned int i = 0; i < n; ++i) {
      uint diffIdx = instMap[thisClassIds[i]];
      //dataset->GetInstanceIndexForID(thisClassIds[i], diffIdx);
      diffClassInstances[thisClass].push_back(diffIdx);
    }
  }
  
  return true;
}

bool DatasetInstance::GetNNearestInstances(unsigned int n,
                                           vector<unsigned int>& closestInstances) {

  if(bestNeighborIds.size() < n) {
    cerr << "ERROR: GetNNearestInstances: k: [" << n
            << "] is larger than the number of neighbors" << endl;
    return false;
  }

  //  cout << "Same sums (" << sameSums.size() << ")" << endl;
  //  copy(neighborSums.begin(), neighborSums.end(), ostream_iterator<double>(cout, "\n"));

  closestInstances.clear();
  for(unsigned int i = 0; i < n; ++i) {
    unsigned int instanceIndex;
    dataset->GetInstanceIndexForID(bestNeighborIds[i], instanceIndex);
    closestInstances.push_back(instanceIndex);
  }

  return true;
}

bool DatasetInstance::ResetNearestNeighbors() {
	bestNeighborIdsSameClass.clear();
	bestNeighborIdsDiffClass.clear();
	bestNeighborIds.clear();
	
	return true;
}

bool DatasetInstance::SetCvSetType(CvSetType newType) {
  bool success = true;
  switch(newType) {
    case CV_TRAIN:
    case CV_HOLDOUT:
    case CV_TEST:
    case CV_NONE:
      cvSetType = newType;
    default:
      success = false;
  }
  return success;
}

CvSetType DatasetInstance::GetCvSetType() {
  return cvSetType;
}
