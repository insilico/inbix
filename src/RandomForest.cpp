/* 
 * RandomForest.cpp
 * 
 * Adapter class to map EC call for Random Forest importance scores
 * to Ranger Random Forest library functions.
 */

#include <iostream>
#include <algorithm>    // std::copy
#include <vector>       // std::vector
#include <iomanip>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sstream>

#include "Insilico.h"

#include "plink.h"
#include "options.h"
#include "helper.h"

#include "RandomForest.h"
#include "AttributeRanker.h"
#include "Dataset.h"
// Ranger random forest project integration - bcw - 9/26/16
#include "DataDouble.h"
#include "Forest.h"
#include "ForestClassification.h"
#include "ForestRegression.h"

using namespace std;

RandomForest::RandomForest(Dataset* ds, Plink* plinkPtr):
	AttributeRanker::AttributeRanker(ds) {
	dataset = ds;
	PP->printLOG(Timestamp() + "Random Forest constructor\n");
  forest = 0;
  PP = plinkPtr;
  try {
    CreateDefaultForestForPheno();
  } catch (std::exception& e) {
    cerr << "RandomForest constructor Error: " << e.what() 
         << " inbix will EXIT now." << std::endl;
    if(forest) {
      delete forest;
      forest = 0;
    }
    shutdown();
  }

  InitializeData(false);
}

RandomForest::RandomForest(Dataset* ds, vector<string> bestAttributeNames):
	AttributeRanker::AttributeRanker(ds) {
	dataset = ds;
	PP->printLOG(Timestamp() + "Random Forest constructor\n");
  forest = 0;
  try {
    CreateDefaultForestForPheno();
  } catch (std::exception& e) {
    cerr << "RandomForest constructor Error: " << e.what() 
         << " inbix will EXIT now." << std::endl;
    if(forest) {
      delete forest;
      forest = 0;
    }
    shutdown();
  }
  
  InitializeData(true);
}

bool RandomForest::InitializeData(bool useMask) {
  // set class variables for reserving memory and other operations
  PP->printLOG(Timestamp() + "Loading data object from inbix internal data structures\n");
  if(!dataset->HasNumerics()) {
    //data = new DataChar();
    error("RandomForest only supports numeric attributes at this time");
  }
  CreateDefaultForestForPheno();
  DataDouble* data = new DataDouble();
  if(useMask) {
    if(data->loadFromDatasetMask(dataset, dataset->GetVariableNames())) {
      error("RandomForest loadFromDatasetMask(&P)");
    }
  } else {
    if(data->loadFromPlink(PP)) {
      error("RandomForest loadFromPlink(&P)");
    }
  }
  try {
    forest->init(par::depvarname, 
            par::memmode, 
            data, 
            par::mtry, 
            par::output_file_name,
            par::ntree, 
            0,
            par::nrfthreads, 
            par::impmeasure,
            minNodeSize,
            par::statusvarname,
            par::do_rfprobability, 
            par::rfreplace,
            par::catvars, 
            par::savemem, 
            par::splitrule,
            par::predall, 
            par::fraction, 
            par::alpha, 
            par::minprop, 
            par::holdout);
  } catch (std::exception& e) {
    stringstream msg;
    msg << "RandomForest InitializeData exception: "
            << e.what() << " inbix will EXIT now." << endl;
    error(msg.str());
    if(forest) {
      delete forest;
      forest = 0;
    }
    shutdown();
  }
  
  return true;
}

RandomForest::~RandomForest() {
  if(forest) delete forest;
}

AttributeScores RandomForest::ComputeScores() {
  PP->printLOG(Timestamp() + "Computing Random Forest variable importance scores\n");
  if(forest) {
    PP->printLOG(Timestamp() + "Running random forest algorithm\n");
    forest->run(par::verbose);
    const vector<double>& rfScores = forest->getVariableImportance();
    const vector<string>& varNames = dataset->GetVariableNames();
    scores.clear();
    for (size_t i=0; i < varNames.size(); ++i) {
      if(varNames[i] != "Class") {
        scores.push_back(make_pair(rfScores[i], varNames[i]));
      }
    }
  } else {
    error("RandomForest::ComputeScores object has not been constructed");
  }

  return scores;
}

double RandomForest::Predict(Dataset* testData) {
  double error = -1;
  
  return error;
}

double RandomForest::GetClassificationError() {
  if(forest) {
    return forest->getOverallPredictionError();
  } else {
    error("RandomForest::GetClassificationError object has not been constructed");
  }
}

void RandomForest::WriteScores(string baseFilename) {
  // safely ignore baseFilename in interface for AttributeRanker
  if(forest) {
    PP->printLOG("Writing output files\n");
    forest->writeOutput();
  } else {
    error("RandomForest::WriteScores object has not been constructed");
  }
}

void RandomForest::WriteScoresInternal() {
  // safely ignore baseFilename in interface for AttributeRanker
  if(forest) {
    PP->printLOG("Writing forest in Ranger internal format\n");
    forest->writeOutputInternal();
  } else {
    error("RandomForest::WriteScoresInternal object has not been constructed");
  }
}

bool RandomForest::CreateDefaultForestForPheno() {
  if(forest) delete forest;
  if(par::bt) {
    PP->printLOG(Timestamp() + "Creating ForestClassification\n");
    forest = new ForestClassification;
  } else {
    PP->printLOG(Timestamp() + "Creating ForestRegression\n");
    forest = new ForestRegression;
  }
  PP->printLOG(Timestamp() + "Initializing forest with inbix data\n");
  if(par::bt) {
    minNodeSize = DEFAULT_MIN_NODE_SIZE_CLASSIFICATION;
  } else {
    minNodeSize = DEFAULT_MIN_NODE_SIZE_REGRESSION;
  }
  forest->setVerboseOutput(&cout);
  
  return true;
}
