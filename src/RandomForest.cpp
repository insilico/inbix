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

#include "Insilico.h"

#include "plink.h"
#include "options.h"
#include "helper.h"

#include "RandomForest.h"
#include "AttributeRanker.h"
// Ranger random forest project integration - bcw - 9/26/16
#include "Data.h"
#include "DataDouble.h"
#include "Forest.h"
#include "ForestClassification.h"
#include "ForestRegression.h"

using namespace std;

RandomForest::RandomForest(Dataset* ds, Plink* plinkPtr):
	AttributeRanker::AttributeRanker(ds) {
	dataset = ds;
	cout << Timestamp()	<< "Random Forest constructor" << endl;
  data = 0;
  forest = 0;
  PP = plinkPtr;
  try {
    PP->printLOG("Creating RandomForest object\n");
    if(par::bt) {
      PP->printLOG("Creating ForestClassification\n");
      forest = new ForestClassification;
    } else {
      PP->printLOG("Creating ForestRegression\n");
      forest = new ForestRegression;
    }
    forest->setVerboseOutput(&cout);
    PP->printLOG("Loading data object from inbix internal data structures\n");
    if(ds->HasNumerics()) {
      data = new DataDouble();
    } else {
      //data = new DataChar();
      error("RandomForest only  supports numeric attributes");
    }
    if(data->loadFromPlink(PP)) {
      error("loadFromPlink(&P)");
    }
    PP->printLOG("Initializing forest with inbix data\n");
    unsigned int minNodeSize = 0;
    if(par::bt) {
      minNodeSize = DEFAULT_MIN_NODE_SIZE_CLASSIFICATION;
    } else {
      minNodeSize = DEFAULT_MIN_NODE_SIZE_REGRESSION;
    }
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
    std::cerr << "RandomForest constructor Error: " << e.what() 
            << " inbix will EXIT now." << std::endl;
    if(forest) {
      delete forest;
      forest = 0;
    }
  }
}

RandomForest::~RandomForest() {
  if(forest) delete forest;
  if(data) delete data;
}

AttributeScores RandomForest::ComputeScores() {
	cout << Timestamp() << "Computing Random Forest variable importance scores"
			<< endl;
  if(forest) {
    PP->printLOG("Running random forest algorithm\n");
    forest->run(par::verbose);
    const vector<double>& rfScores = forest->getVariableImportance();
    const vector<string>& varNames = data->getVariableNames();
    scores.clear();
    for (size_t i=0; i < varNames.size(); ++i) {
      scores.push_back(make_pair(rfScores[i], varNames[i]));
    }
  } else {
    error("RandomForest::ComputeScores object has not been constructed");
  }

  return scores;
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

