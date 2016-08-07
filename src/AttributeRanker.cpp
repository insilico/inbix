/*
 * AttributeRanker.cpp
 *
 *  Created on: Aug 13, 2012
 *      Author: bwhite
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "AttributeRanker.h"
#include "Dataset.h"
#include "Insilico.h"

using namespace std;

AttributeRanker::AttributeRanker(Dataset* ds) {
	dataset = ds;
	classificationAccuracy = 1.0;
  k = 0;
  normalizeScores = false;
}

AttributeRanker::~AttributeRanker() {
}

AttributeScores AttributeRanker::GetScores() {
	return scores;
}

void AttributeRanker::WriteScores(string baseFilename) {
	string resultsFilename = baseFilename + ".ranks";
	ofstream outFile;
	outFile.open(resultsFilename.c_str());
	if (outFile.bad()) {
		cerr << "ERROR: Could not open scores file " << resultsFilename
				<< "for writing" << endl;
		exit(1);
	}
	PrintScores(outFile);
	outFile.close();
}

void AttributeRanker::PrintScores(ofstream& outStream) {
	for (AttributeScoresCIt scoresIt = scores.begin(); scoresIt != scores.end();
			++scoresIt) {
		outStream << scoresIt->first << "\t" << scoresIt->second << endl;
	}
}

double AttributeRanker::GetClassificationError() {
	return classificationAccuracy;
}

bool AttributeRanker::DoNormalize() {
	return normalizeScores;
}

bool AttributeRanker::SetK(unsigned int newK) {

  if(dataset->HasContinuousPhenotypes()) {
    if(newK && (newK < dataset->NumInstances())) {
      k = newK;
      return true;
    }
    return false;
  }

  map<ClassLevel, vector<unsigned int> > classLevels =
          dataset->GetClassIndexes();
  map<ClassLevel, vector<unsigned int> >::const_iterator ciIt;
  unsigned int minSampleCount = dataset->NumInstances();
  unsigned int maxSampleCount = 0;
  for(ciIt = classLevels.begin(); ciIt != classLevels.end(); ++ciIt) {
    if(ciIt->second.size() < minSampleCount) {
      minSampleCount = ciIt->second.size();
    }
    if(ciIt->second.size() > maxSampleCount) {
      maxSampleCount = ciIt->second.size();
    }
  }

  if(minSampleCount < newK) {
    cout << Timestamp() << "WARNING: Minimum class size " << minSampleCount
            << " is less than k " << newK << endl;
    cout << Timestamp() << "WARNING: Setting k nearest neighbors to " <<
            minSampleCount << endl;
    k = minSampleCount;
  } else {
    if(newK > maxSampleCount) {
      cout << Timestamp() << "WARNING: Maximum class size " << maxSampleCount
              << " is less than k " << newK << endl;
      cout << Timestamp() << "WARNING: Setting k nearest neighbors to " <<
              maxSampleCount - 1 << endl;
      k = maxSampleCount - 1;
    }
    else {
      k = newK;
    }
  }

  return true;
}
