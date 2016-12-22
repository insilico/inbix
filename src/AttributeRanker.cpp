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
#include "options.h"
#include "helper.h"

using namespace std;

AttributeRanker::AttributeRanker(Dataset* ds) {
	dataset = ds;
	classificationAccuracy = 1.0;
  normalizeScores = false;
  k = 10;
}

AttributeRanker::~AttributeRanker() {
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

AttributeScores AttributeRanker::GetScores() {
//  AttributeScores returnScores;
////  unsigned int nameIdx = 0;
////  vector<string> maskNames = dataset->MaskGetAllVariableNames();
//  AttributeScoresCIt scoresIt = scores.begin();
//  for(; scoresIt != scores.end(); ++scoresIt) {
//    returnScores.push_back(make_pair((*scoresIt).first, (*scoresIt).second));
//  }
  return scores;
}


void AttributeRanker::WriteScores(string baseFilename) {
	string resultsFilename = baseFilename + ".relieff.tab";
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

void AttributeRanker::SetNormalize(bool switchTF) {
  normalizeScores = switchTF;
}

bool AttributeRanker::GetNormalizeFlag() {
	return normalizeScores;
}

bool AttributeRanker::NormalizeScores() {
  PP->printLOG(Timestamp() + "Normalizing ReliefF scores to 0-1\n");
  if(!scores.size()) {
    error("AttributeRanker::NormalizeScores() Cannot normalize zero-length scores\n");
  }
  ScoreVarPair firstScore = scores[0];
  double minScore = firstScore.first;
  double maxScore = firstScore.first;
  AttributeScoresCIt scoresIt = scores.begin();
  for(; scoresIt != scores.end(); ++scoresIt) {
    ScoreVarPair thisScore = *scoresIt;
    if(thisScore.first < minScore) {
      minScore = thisScore.first;
    }
    if(thisScore.first > maxScore) {
      maxScore = thisScore.first;
    }
  }
  // with min and max scores, scale values to min=0, max=1
  AttributeScores newScores;
  double scoreRange = maxScore - minScore;
  if(scoreRange < par::epsilon) {
    scoreRange = par::epsilon;
  }
  for(AttributeScoresIt it = scores.begin(); it != scores.end(); ++it) {
    ScoreVarPair thisScore = *it;
    double key = thisScore.first;
    string val = thisScore.second;
    newScores.push_back(make_pair((key - minScore) / scoreRange, val));
  }
  scores.clear();
  scores = newScores;
  
  return true;
}
