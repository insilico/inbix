/* 
 * File:   EvaporativeCooling.cpp
 * Author: billwhite
 * 
 * Created on July 14, 2011, 9:25 PM
 *
 * Implements the Evaporative Cooling algorithm in:
 * McKinney, et. al. "Capturing the Spectrum of Interaction Effects in Genetic
 * Association Studies by Simulated Evaporative Cooling Network Analysis."
 * PLoS Genetics, Vol 5, Issue 3, 2009.
  *
 * Made even more generic with main effects and interaction effects algorithms
 * in a class hierarchy from a AttributeRanker base. 8/12/12
 * 
 * Modified for inbix. 9/29/16
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <omp.h>

#include <gsl/gsl_rng.h>
#include <boost/lexical_cast.hpp>

#include "plink.h"
#include "options.h"

// EC project
#include "Insilico.h"
#include "Dataset.h"
#include "Statistics.h"
#include "StringUtils.h"
#include "EvaporativeCooling.h"
#include "Deseq.h"
#include "Edger.h"
#include "RandomForest.h"
#include "ReliefF.h"
#include "RReliefF.h"
#include "ReliefFSeq.h"
#include "PlinkInternalsDataset.h"
#include "helper.h"

using namespace std;
using namespace insilico;

EvaporativeCooling::EvaporativeCooling(Dataset* ds, Plink* plinkPtr,
        AnalysisType anaType) {
	PP->printLOG(Timestamp() + "Evaporative Cooling constructor\n");
	if (ds) {
		dataset = ds;
	} else {
		error("EvaporativeCooling constructor: PLINK data set is not initialized");
	}
  PP = plinkPtr;
	analysisType = anaType;
	if(par::ecOptimizeTemp) {
		optimizeTemperature = true;
	}
	optimalTemperature = 1.0;
	bestClassificationError = 1.0;
	// set the number of target attributes
	numTargetAttributes = par::ecNumTarget;
	if (numTargetAttributes == 0) {
		numTargetAttributes = ds->NumVariables();
		numToRemovePerIteration = 0;
	}
	if (numTargetAttributes > dataset->NumVariables()) {
		error("--ec-num-target must be less than or equal to the number of attributes in the data set\n");
		exit(EXIT_FAILURE);
	}
	PP->printLOG(Timestamp() + "EC is removing attributes until best " +
			int2str(numTargetAttributes) + " remain\n");
	/// set the EC steps to perform
  string ecAlgParam = par::ecAlgorithmSteps;
  string meAlgorithmParam = par::ecMeAlgorithm;
  string itAlgorithmParam = par::ecItAlgorithm;
	algorithmType = EC_ALG_ME_IT;
  if (ecAlgParam == "all") {
    algorithmType = EC_ALG_ME_IT;
    PP->printLOG(Timestamp() + "Running EC in standard mode: Main effects + interaction effects\n");
  } else {
    if (ecAlgParam == "me") {
      algorithmType = EC_ALG_ME_ONLY;
      PP->printLOG(Timestamp() + "Running EC in main effects only mode\n");
    } else {
      if (ecAlgParam == "it") {
        algorithmType = EC_ALG_IT_ONLY;
        PP->printLOG(Timestamp() + "Running EC in interactions effects only mode\n");
      } else {
        error("-ec-algorithm-steps must be one of: all, me or it\n");
      }
    }
  }
	/// set the main effects algorithm
	meAlgorithmType = EC_ME_ALG_RJ;
  string ecMeAlgParam = par::ecMeAlgorithm;
  if (ecMeAlgParam == "randomforest") {
    meAlgorithmType = EC_ME_ALG_RJ;
    PP->printLOG(Timestamp() + "EC main effects algorithm set to: Random Forest\n");
  } else {
    if (ecMeAlgParam == "deseq") {
      meAlgorithmType = EC_ME_ALG_DESEQ;
      PP->printLOG(Timestamp() + "Running EC in main effects algorithm set to: DESeq\n");
    } else {
      error("--ec-me-algorithm must be one of: rrelieff or deseq\n");
    }
  }
	maineffectAlgorithm = NULL;
	switch(meAlgorithmType) {
	case EC_ME_ALG_RJ:
    PP->printLOG(Timestamp() + "Creating RandomForest main effects object\n");
		maineffectAlgorithm = new RandomForest(ds, PP);
		break;
	case EC_ME_ALG_DESEQ:
    PP->printLOG(Timestamp() + "Creating DESeq main effects object\n");
		maineffectAlgorithm = new Deseq(ds);
		break;
	case EC_ME_ALG_EDGER:
    PP->printLOG(Timestamp() + "Creating edgeR main effects object\n");
		maineffectAlgorithm = new Edger(ds);
		break;
	}
	/// set the interaction algorithm
	itAlgorithmType = EC_IT_ALG_RF;
  string ecItAlgParam = par::ecItAlgorithm;
  if (ecItAlgParam == "relieff") {
    itAlgorithmType = EC_IT_ALG_RF;
    PP->printLOG(Timestamp() + "EC interaction effects algorithm set to: ReliefF\n");
  } else {
    if (ecItAlgParam == "reliefseq") {
      itAlgorithmType = EC_IT_ALG_RFSEQ;
      PP->printLOG(Timestamp() + "EC interaction effects algorithm set to: ReliefFSeq\n");
    } else {
      error("--ec-it-algorithm must be one of: relieff or reliefseq");
    }
  }
	interactionAlgorithm = NULL;
	switch(itAlgorithmType) {
	case EC_IT_ALG_RF:
    if(ds->HasContinuousPhenotypes()) {
      PP->printLOG(Timestamp() + "Creating Regression ReliefF object" );
      interactionAlgorithm = new RReliefF(ds, PP);
    } else {
      PP->printLOG(Timestamp() + "Creating Standard ReliefF object\n");
      interactionAlgorithm = new ReliefF(ds, PP, anaType);
    }
		break;
	case EC_IT_ALG_RFSEQ:
    PP->printLOG(Timestamp() + "Creating ReliefFSeq object\n");
		interactionAlgorithm = new ReliefFSeq(ds, PP);
		break;
	}
	// set the number of attributes to remove per iteration
  if(par::do_iterative_removal_ec) {
    numToRemovePerIteration = par::ecIterNumToRemove;
    if(!numToRemovePerIteration) {
      unsigned int iterPercentToRemove = par::ecIterPercentToRemove;
      numToRemovePerIteration = (unsigned int) (((double) iterPercentToRemove
        / 100.0) * dataset->NumAttributes());
    }
    PP->printLOG(Timestamp() + "EC will remove " + int2str(numToRemovePerIteration) + 
        " attributes on first iteration\n");
  }

	// openMP multi-threading setup
	numRFThreads = omp_get_num_procs();
	PP->printLOG(Timestamp() + "EC has " + int2str(numRFThreads) + " threads\n");
} // end of constructor

EvaporativeCooling::~EvaporativeCooling() {
	if (interactionAlgorithm) {
		delete interactionAlgorithm;
	}
	if (maineffectAlgorithm) {
		delete maineffectAlgorithm;
	}
}

bool EvaporativeCooling::ComputeECScores() {
  stringstream msg;
	unsigned int numWorkingAttributes = dataset->NumVariables();
	if (numWorkingAttributes < numTargetAttributes) {
		msg << "ERROR: The number of attributes in the data set "
				<< numWorkingAttributes
				<< " is less than the number of target attributes "
				<< numTargetAttributes << endl;
		error(msg.str());
	}

	// EC algorithm as in Figure 5, page 10 of the paper referenced
	// at top of this file. Modified per Brett's email to not do the
	// varying temperature and classifier accuracy optimization steps.
	// Added temperature optimization - 4/9/12
	unsigned int iteration = 1;
	optimalTemperature = 1.0;
	while (numWorkingAttributes >= numTargetAttributes) {
		pair<unsigned int, unsigned int> titvCounts;
    double titvRatio = 0;
    if(dataset->HasGenotypes()) {
      titvCounts = dataset->GetAttributeTiTvCounts();
      titvRatio = titvCounts.first;
      if(titvCounts.second) {
        titvRatio = (double) titvCounts.first / (double) titvCounts.second;
      }
    }

		msg << Timestamp()
				<< "----------------------------------------------------"
				<< "-------------------------" << endl;
		msg << Timestamp() << "EC algorithm...iteration: " << iteration
				<< ", working attributes: " << numWorkingAttributes
				<< ", target attributes: " << numTargetAttributes
				<< ", temperature: " << optimalTemperature
				<< ", best classification error: " << bestClassificationError << endl;
    if(dataset->HasGenotypes()) {
      msg << Timestamp()
          << "Ti/Tv: transitions: " << titvCounts.first
          << " transversions: " << titvCounts.second
          << " ratio: " << titvRatio
          << endl;
    }
		PP->printLOG(msg.str());
    msg.clear();
    
		// -------------------------------------------------------------------------
		// run main effects algorithm and get the normalized scores for use in EC
		if ((algorithmType == EC_ALG_ME_IT) || (algorithmType == EC_ALG_ME_ONLY)) {
      if(!RunRandomForest()) {
				error("In EC algorithm: Random forest failed");
			}
        
			PP->printLOG(Timestamp() + "Running Random Forest\n");
			double classificationError = maineffectAlgorithm->GetClassificationError();
      PP->printLOG(Timestamp() + "Classifier error: " + dbl2str(classificationError) + "\n");
			if(classificationError < bestClassificationError) {
				bestClassificationError = classificationError;
			}
			classificationErrors.push_back(classificationError);
			// Random forest standalone runs
			if ((algorithmType == EC_ALG_ME_ONLY)
					&& (numWorkingAttributes == numTargetAttributes)) {
				sort(maineffectScores.begin(), maineffectScores.end(), scoresSortDesc);
				ecScores.resize(numTargetAttributes);
				copy(maineffectScores.begin(),
						maineffectScores.begin() + numTargetAttributes,
						ecScores.begin());
				return true;
			}
		}

		// -------------------------------------------------------------------------
		// run interaction effects algorithm and get normalized score for use in EC
		if ((algorithmType == EC_ALG_ME_IT) || (algorithmType == EC_ALG_IT_ONLY)) {
			PP->printLOG(Timestamp() + "Running ReliefF\n");
			if (!RunReliefF()) {
				error("In EC algorithm: ReliefF failed");
			}
			// ReliefF standalone runs
			if ((algorithmType == EC_ALG_IT_ONLY)
					&& (numWorkingAttributes == numTargetAttributes)) {
				sort(interactionScores.begin(), interactionScores.end(), scoresSortDesc);
				ecScores.resize(numTargetAttributes);
				copy(interactionScores.begin(), interactionScores.begin() + numTargetAttributes,
						ecScores.begin());
				return true;
			}
		}

		// -------------------------------------------------------------------------
		// compute free energy for all attributes
		PP->printLOG(Timestamp() + "Computing free energy\n");
		if (!ComputeFreeEnergy(optimalTemperature)) {
			error("In EC algorithm: ComputeFreeEnergy failed\n");
		}
		// PrintAllScoresTabular();
		// PrintKendallTaus();

		// -------------------------------------------------------------------------
		// optimize the temperature by sampling a set of delta values around T
		if(optimizeTemperature) {
			PP->printLOG(Timestamp() + "Optimizing coupling temperature T\n");
			vector<double> temperatureDeltas;
			temperatureDeltas.push_back(-0.2);
			temperatureDeltas.push_back(0.2);
			optimalTemperature = OptimizeTemperature(temperatureDeltas);
			temperatures.push_back(optimalTemperature);
		}

		// write scores for each iteration
		stringstream scoreFilename;
		scoreFilename << "ec." << iteration << ".scores.dat";
		ofstream outFile;
		outFile.open(scoreFilename.str().c_str());
		if(outFile.bad()) {
			error("Could not open scores file " + scoreFilename.str() + "for writing");
		}
		PP->printLOG(Timestamp() + "Writing ALL EC scores to [" + scoreFilename.str() + "]\n");
		PrintAllAttributeScores(outFile);
		outFile.close();

		// -------------------------------------------------------------------------
		// remove the worst attributes and iterate
		PP->printLOG(Timestamp() + "Removing the worst attributes\n");
		unsigned int numToRemove = numToRemovePerIteration;
		numToRemoveNextIteration = numToRemove - numToRemovePerIteration;
		if (par::ecIterPercentToRemove) {
			unsigned int iterPercentToRemove = par::ecIterPercentToRemove;
			numToRemove = (int) (((double) iterPercentToRemove / 100.0)
					* dataset->NumVariables());
			numToRemoveNextIteration = (int) (((double) iterPercentToRemove / 100.0)
					* numToRemove);
			numToRemovePerIteration = numToRemove;
		}
		if ((numWorkingAttributes - numToRemove) < numTargetAttributes) {
			numToRemove = numWorkingAttributes - numTargetAttributes;
		}
		if (numToRemove < 1) {
//      cerr << "ERROR: Number of attributes to remove is less than one." << endl;
//      return false;
			break;
		}
		PP->printLOG(Timestamp() + 
      "EvaporativeCooling::ComputeECScores removing the worst " + 
      int2str(numToRemove) + " attributes\n");
		if (!RemoveWorstAttributes(numToRemove)) {
			error("EvaporativeCooling::ComputeECScores: RemoveWorstAttribute failed");
		}
    cout << "DEBUG: numToRemove: " << numToRemove << endl;
		numWorkingAttributes -= numToRemove;

		++iteration;
    maineffectAlgorithm->InitializeData(true);
    interactionAlgorithm->ResetForNextIteration();
	}

	PP->printLOG(Timestamp() + "EC algorithm ran for " + int2str(iteration) + " iterations\n");

	// remaining free energy attributes are the ones we want to write as a
	// new dataset to be analyzed with (re)GAIN + SNPrank
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);
	ecScores.resize(numTargetAttributes);
	copy(freeEnergyScores.begin(),
			freeEnergyScores.begin() + numTargetAttributes,
			ecScores.begin());

	return true;
}

AttributeScores& EvaporativeCooling::GetMaineffectScores() {
	return maineffectScores;
}

AttributeScores& EvaporativeCooling::GetInteractionScores() {
	return interactionScores;
}

AttributeScores& EvaporativeCooling::GetECScores() {
	return ecScores;
}

EcAlgorithmType EvaporativeCooling::GetAlgorithmType() {
	return algorithmType;
}

void EvaporativeCooling::PrintAttributeScores(ofstream& outStream) {
	for(AttributeScoresCIt ecScoresIt = ecScores.begin(); ecScoresIt != ecScores.end();
			++ecScoresIt) {
		outStream << fixed << setprecision(8) << (*ecScoresIt).first << "\t"
				<< (*ecScoresIt).second << endl;
	}
}

void EvaporativeCooling::PrintAllAttributeScores(ofstream& outStream) {
	for(AttributeScoresCIt ecScoresIt = freeEnergyScores.begin();
			ecScoresIt != freeEnergyScores.end();	++ecScoresIt) {
		outStream << fixed << setprecision(8) << (*ecScoresIt).first << "\t"
				<< (*ecScoresIt).second << endl;
	}
	for (AttributeScoresCIt ecScoresIt = evaporatedAttributes.begin();
			ecScoresIt != evaporatedAttributes.end();	++ecScoresIt) {
		outStream << fixed << setprecision(8) << 0 << "\t"
				<< (*ecScoresIt).second << endl;
	}
}

void EvaporativeCooling::PrintMainEffectAttributeScores(ofstream& outStream) {
	sort(maineffectScores.begin(), maineffectScores.end(), scoresSortDesc);
	for(AttributeScoresCIt rjScoresIt = maineffectScores.begin();
			rjScoresIt != maineffectScores.end(); ++rjScoresIt) {
		outStream << fixed << setprecision(8) << (*rjScoresIt).first << "\t"
				<< (*rjScoresIt).second << endl;
	}
}

void EvaporativeCooling::PrintInteractionAttributeScores(ofstream& outStream) {
	sort(interactionScores.begin(), interactionScores.end(), scoresSortDesc);
	for(AttributeScoresCIt rfScoresIt = interactionScores.begin();
			rfScoresIt != interactionScores.end(); ++rfScoresIt) {
		outStream << fixed << setprecision(8) << (*rfScoresIt).first << "\t"
				<< (*rfScoresIt).second << endl;
	}
}

void EvaporativeCooling::WriteAttributeScores(string baseFilename) {
  stringstream msg;
	string resultsFilename = baseFilename;
	ofstream outFile;
	// added 9/26/11 for reflecting the fact that only parts of the
	// complete EC algorithm were performed
	switch (algorithmType) {
	case EC_ALG_ME_IT:
		resultsFilename = baseFilename + ".ec.tab";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			error("ERROR: Could not open scores file " + resultsFilename +
					  "for writing");
		}
		msg << Timestamp()
				<< "Writing EC scores to [" + resultsFilename + "]" << endl;
    PP->printLOG(msg.str());
    msg.str("");
		PrintAttributeScores(outFile);
		outFile.close();

		resultsFilename = baseFilename + ".ec.randomforest.tab";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			error("ERROR: Could not open scores file " + resultsFilename + "for writing");
		}
		msg << Timestamp()
				<< "Writing EC main effects scores to [" + resultsFilename + "]" << endl;
    PP->printLOG(msg.str());
    msg.str("");
		PrintMainEffectAttributeScores(outFile);
		outFile.close();

		resultsFilename = baseFilename + ".ec.relieff.tab";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			error("ERROR: Could not open scores file " + resultsFilename + "for writing");
		}
		msg << Timestamp()
				<< "Writing EC interaction effects scores to [" + resultsFilename + "]"
				<< endl;
    PP->printLOG(msg.str());
    msg.str("");
		PrintInteractionAttributeScores(outFile);
		outFile.close();
		break;
	case EC_ALG_ME_ONLY:
		resultsFilename += ".randomforest.tab";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			error("ERROR: Could not open scores file " + resultsFilename + "for writing");
		}
		msg << Timestamp()
				<< "Writing EC main effects scores to [" + resultsFilename + "]" << endl;
    PP->printLOG(msg.str());
    msg.str("");
		PrintAttributeScores(outFile);
		outFile.close();
		break;
	case EC_ALG_IT_ONLY:
		if(itAlgorithmType == EC_IT_ALG_RF) {
			resultsFilename += ".relieff.tab";
		}
		if(itAlgorithmType == EC_IT_ALG_RFSEQ) {
			resultsFilename += ".reliefseq.tab";
		}
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			error("ERROR: Could not open scores file " + resultsFilename + "for writing");
		}
		msg << Timestamp()
				<< "Writing EC interaction effects scores to [" + resultsFilename + "]"
				<< endl;
    PP->printLOG(msg.str());
    msg.str("");
		PrintAttributeScores(outFile);
		outFile.close();
		break;
	default:
		// we should not get here by the CLI front end but it is possible to call
		// this from other programs in the future or when used as a library
		// TODO: better message
		error("ERROR: Attempting to write attribute scores before the analysis type was determined.");
	}
}

void EvaporativeCooling::WriteClassificationErrors(string filename) {
	ofstream outFile;
	outFile.open(filename.c_str());
	if (outFile.bad()) {
		cerr << "ERROR: Could not open classification errors file "
				<< filename	<< "for writing" << endl;
		exit(1);
	}
	cout << Timestamp()
			<< "Writing EC classification errors to [" + filename + "]" << endl;
	vector<double>::const_iterator ceIt = classificationErrors.begin();
	for(; ceIt != classificationErrors.end(); ++ceIt) {
		outFile << *ceIt << endl;
	}
	outFile.close();
}

void EvaporativeCooling::WriteTemperatures(string filename) {
	ofstream outFile;
	outFile.open(filename.c_str());
	if (outFile.bad()) {
		cerr << "ERROR: Could not open temperatures file "
				<< filename	<< "for writing" << endl;
		exit(1);
	}
	cout << Timestamp()
			<< "Writing EC temperatures to [" + filename + "]" << endl;
	vector<double>::const_iterator tempsIt = temperatures.begin();
	for(; tempsIt != temperatures.end(); ++tempsIt) {
		outFile << *tempsIt << endl;
	}
	outFile.close();
}

bool EvaporativeCooling::PrintAllScoresTabular() {
	// sanity checks
	if (maineffectScores.size() != interactionScores.size()) {
		cerr
				<< "ERROR: Random Forest and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}
	if (freeEnergyScores.size() != interactionScores.size()) {
		cerr
				<< "ERROR: Random Forest and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}

	sort(maineffectScores.begin(), maineffectScores.end(), scoresSortDesc);
	sort(interactionScores.begin(), interactionScores.end(), scoresSortDesc);
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);

	cout << "\t\t\tE (RF)\t\tS (RJ)\t\tF (free energy)\n";
	unsigned int numScores = freeEnergyScores.size();
	for (unsigned int i = 0; i < numScores; ++i) {
		pair<double, string> thisRJScores = maineffectScores[i];
		pair<double, string> thisRFScores = interactionScores[i];
		pair<double, string> thisFEScores = freeEnergyScores[i];
		printf("\t\t\t%s\t%6.4f\t%s\t%6.4f\t%s\t%6.4f\n",
				thisRFScores.second.c_str(), thisRFScores.first,
				thisRJScores.second.c_str(), thisRJScores.first,
				thisFEScores.second.c_str(), thisFEScores.first);
	}

	return true;
}

bool EvaporativeCooling::PrintKendallTaus() {
	// sanity checks
	if (maineffectScores.size() != interactionScores.size()) {
		cerr
				<< "ERROR: Random Forest and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}
	if (freeEnergyScores.size() != interactionScores.size()) {
		cerr
				<< "ERROR: Random Forest and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}

	sort(maineffectScores.begin(), maineffectScores.end(), scoresSortDesc);
	sort(interactionScores.begin(), interactionScores.end(), scoresSortDesc);
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);

	vector<string> rjNames;
	vector<string> rfNames;
	vector<string> feNames;
	unsigned int numScores = freeEnergyScores.size();
	for (unsigned int i = 0; i < numScores; ++i) {
		pair<double, string> thisRJScores = maineffectScores[i];
		pair<double, string> thisRFScores = interactionScores[i];
		pair<double, string> thisFEScores = freeEnergyScores[i];
		rjNames.push_back(thisRJScores.second);
		rfNames.push_back(thisRFScores.second);
		feNames.push_back(thisFEScores.second);
	}

	double tauRJRF = KendallTau(rjNames, rfNames);
	double tauRJFE = KendallTau(rjNames, feNames);
	double tauRFFE = KendallTau(rfNames, feNames);

	cout << "\t\t\tKendall tau's: " << "RJvRF: " << tauRJRF << ", RJvFE: "
			<< tauRJFE << ", RFvFE: " << tauRFFE << endl;

	return true;
}


bool EvaporativeCooling::RunRandomForest() {
	maineffectScores = maineffectAlgorithm->ComputeScores();
	if(maineffectScores.size() == 0) {
		error("ERROR: RunRandomForest: No scores computed");
	}
  
  if(par::do_normalize_scores) {
    cout << Timestamp() << "Normalizing ReliefF scores to 0-1" << endl;
    pair<double, string> firstScore = maineffectScores[0];
    double minRFScore = firstScore.first;
    double maxRFScore = firstScore.first;
    AttributeScoresCIt rfScoresIt = maineffectScores.begin();
    for (; rfScoresIt != maineffectScores.end(); ++rfScoresIt) {
      pair<double, string> thisScore = *rfScoresIt;
      if (thisScore.first < minRFScore) {
        minRFScore = thisScore.first;
      }
      if (thisScore.first > maxRFScore) {
        maxRFScore = thisScore.first;
      }
    }
    // normalize attribute scores if necessary
    if (minRFScore == maxRFScore) {
      cout << Timestamp() << "WARNING: random forest min and max scores are the same. "
          << "No normalization necessary" << endl;
      return true;
    }
    AttributeScores newRFScores;
    double rfRange = maxRFScore - minRFScore;
    for (AttributeScoresIt it = maineffectScores.begin(); it != maineffectScores.end(); ++it) {
      pair<double, string> thisScore = *it;
      double key = thisScore.first;
      string val = thisScore.second;
      newRFScores.push_back(make_pair((key - minRFScore) / rfRange, val));
    }

    maineffectScores.clear();
    maineffectScores = newRFScores;
  }
  return true;
}

bool EvaporativeCooling::RunReliefF() {
	/// postcondition: interactionScores contains the newly-computed scores
	interactionScores = interactionAlgorithm->ComputeScores();
	if(interactionScores.size() == 0) {
		error("ERROR: RunReliefF: No scores computed");
	}

	// added for rnaSeq analysis - bcw - 8/21/12
	// normalizing "stretches" the distribution to many zeroes, thus
	// complicating ranking of many ties; don't do it!
	if(itAlgorithmType == EC_IT_ALG_RFSEQ) {
		cout << Timestamp() << "rnaSeq skipping normalization."	<< endl;
		return true;
	}

  if(par::do_normalize_scores) {
    cout << Timestamp() << "Normalizing ReliefF scores to 0-1" << endl;
    pair<double, string> firstScore = interactionScores[0];
    double minRFScore = firstScore.first;
    double maxRFScore = firstScore.first;
    AttributeScoresCIt rfScoresIt = interactionScores.begin();
    for (; rfScoresIt != interactionScores.end(); ++rfScoresIt) {
      pair<double, string> thisScore = *rfScoresIt;
      if (thisScore.first < minRFScore) {
        minRFScore = thisScore.first;
      }
      if (thisScore.first > maxRFScore) {
        maxRFScore = thisScore.first;
      }
    }
    // normalize attribute scores if necessary
    if (minRFScore == maxRFScore) {
      cout << Timestamp() << "WARNING: Relief-F min and max scores are the same. "
          << "No normalization necessary" << endl;
      return true;
    }
    AttributeScores newRFScores;
    double rfRange = maxRFScore - minRFScore;
    for (AttributeScoresIt it = interactionScores.begin(); it != interactionScores.end(); ++it) {
      pair<double, string> thisScore = *it;
      double key = thisScore.first;
      string val = thisScore.second;
      newRFScores.push_back(make_pair((key - minRFScore) / rfRange, val));
    }

    interactionScores.clear();
    interactionScores = newRFScores;
  }
  
	return true;
}

bool EvaporativeCooling::ComputeFreeEnergy(double temperature) {
  stringstream msg;
	if (algorithmType == EC_ALG_ME_IT) {
		if (maineffectScores.size() != interactionScores.size()) {
			msg << "EvaporativeCooling::ComputeFreeEnergy scores lists are "
              "unequal. RandomForest: " << maineffectScores.size() 
              << " vs. Relief-F: "  << interactionScores.size();
			error(msg.str());
		}
	}

	freeEnergyScores.clear();
	AttributeScoresCIt rjIt = maineffectScores.begin();
	AttributeScoresCIt rfIt = interactionScores.begin();
	switch (algorithmType) {
	case EC_ALG_ME_IT:
		sort(maineffectScores.begin(), maineffectScores.end(), scoresSortAscByName);
		sort(interactionScores.begin(), interactionScores.end(), scoresSortAscByName);
		for (; rjIt != maineffectScores.end(); ++rjIt, ++rfIt) {
			string val = rjIt->second;
			double key = rjIt->first;
			freeEnergyScores.push_back(
					make_pair((*rfIt).first + (temperature * key), val));
		}
		break;
	case EC_ALG_ME_ONLY:
		for (; rjIt != maineffectScores.end(); ++rjIt) {
			freeEnergyScores.push_back(make_pair(rjIt->first, rjIt->second));
		}
		break;
	case EC_ALG_IT_ONLY:
		for (; rfIt != interactionScores.end(); ++rfIt) {
			freeEnergyScores.push_back(make_pair(rfIt->first, rfIt->second));
		}
		break;
	default:
		msg << "EvaporativeCooling::ComputeFreeEnergy: "
				<< "could not determine EC algorithm type";
		error(msg.str());
	}

	return true;
}

bool EvaporativeCooling::RemoveWorstAttributes(unsigned int numToRemove) {
  if(par::verbose) PP->printLOG("ENTERING EvaporativeCooling::RemoveWorstAttributes\n");
  stringstream msg;
	unsigned int numToRemoveAdj = numToRemove;
	unsigned int numAttr = dataset->NumAttributes();
	if ((numAttr - numToRemove) < numTargetAttributes) {
		msg << Timestamp() << "WARNING: attempt to remove " << numToRemove
				<< " attributes which will remove more than target "
				<< "number of attributes " << numTargetAttributes << ". Adjusting"
				<< endl;
    PP->printLOG(msg.str());
    msg.str("");
		numToRemoveAdj = numAttr - numTargetAttributes;
	}
	msg << Timestamp() << "EvaporativeCooling::RemoveWorstAttributes " 
          << numToRemoveAdj << " attributes" << endl;
  PP->printLOG(msg.str());
  msg.str("");
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortAsc);
	for (unsigned int i=0; i < numToRemoveAdj; ++i) {
		// worst score and attribute name
		pair<double, string> worst = freeEnergyScores[i];
    if(par::verbose) {
      msg << Timestamp() << "current worst attribute" << i << ": " 
              << worst.second << " (" << worst.first << ")" << endl;
      PP->printLOG(msg.str());
      msg.str("");
    }
		// save worst
    if(par::verbose) PP->printLOG("saving worst attributes to evaporatedAttributes\n");
		evaporatedAttributes.push_back(worst);
		// remove the attribute from those under consideration
    if(par::verbose) PP->printLOG("calling dataset->MaskRemoveVariable\n");
		if (!dataset->MaskRemoveVariable(worst.second)) {
			error("Could not remove worst attribute: " + worst.second + "\n");
		}
	}
  if(par::verbose) PP->printLOG("RETURNING FROM EvaporativeCooling::RemoveWorstAttributes\n");

	return true;
}

double EvaporativeCooling::OptimizeTemperature(vector<double> deltas) {
  stringstream msg;
  msg << Timestamp() << "--- OPTIMIZER BEGIN: Classification error to beat: "
          << dbl2str(bestClassificationError) << endl;
	PP->printLOG(msg.str());
  msg.clear();
	/// for each delta, run a classifier on the best attributes according
	/// to the free energy and update best temperature
	vector<double>::const_iterator deltaIt = deltas.begin();
	AttributeScores bestFreeEnergyScores = freeEnergyScores;
	for(; deltaIt != deltas.end(); ++deltaIt) {
		double thisTemp = optimalTemperature + *deltaIt;
		ComputeFreeEnergy(thisTemp);
		double thisClassificationError = ComputeClassificationErrorRandomForest();
		msg << Timestamp()
				<< "OPTIMIZER: Trying temperature: " << thisTemp
				<< " => Classification Error: " << thisClassificationError << endl;
    PP->printLOG(msg.str());
    msg.clear();
		/// if classification error is lower at this delta, update best temperature
		/// and best classification error
		if(thisClassificationError < bestClassificationError) {
			msg << Timestamp() << "--- OPTIMIZER: found better temperature: "
					<< thisTemp << endl;
      PP->printLOG(msg.str());
      msg.clear();
			bestClassificationError = thisClassificationError;
			optimalTemperature = thisTemp;
			bestFreeEnergyScores = freeEnergyScores;
		}
	}
	msg << Timestamp() << "--- OPTIMIZER RESULTS: "
			<< "Temperature: " << optimalTemperature
			<< ", Classification error: " << bestClassificationError << endl;
  PP->printLOG(msg.str());
	freeEnergyScores = bestFreeEnergyScores;

	return optimalTemperature;
}

double EvaporativeCooling::ComputeClassificationErrorRandomForest() {
  stringstream msg;
	/// get the best attribute names based on free energy score
	unsigned int numBest = freeEnergyScores.size() - numToRemoveNextIteration;
	if(!numBest) {
		msg << "ERROR: Best results calculation results in zero attributes" << endl;
		msg << "Number of best to use in classifier: " << numBest << endl;
		msg << "Free energy scores: " << freeEnergyScores.size() << endl;
		msg << "Number to remove next iteration: " << numToRemoveNextIteration << endl;
		error(msg.str());
	}
  
  // ---------------------------------------------------------------------------
  PP->printLOG(Timestamp() + "Getting best free energy attribute scores for RandomForest\n");
  vector<string> bestAttributes;
	unsigned int numCopied = 0;
	AttributeScoresCIt scoreIt = freeEnergyScores.begin();
	for(; numCopied < numBest && scoreIt != freeEnergyScores.end();
			++numCopied, ++scoreIt) {
		bestAttributes.push_back(scoreIt->second);
	}
	if(!bestAttributes.size() || (bestAttributes.size() != numBest)) {
		msg << "ERROR: could not get " << numBest << " attributes, got: "
				<< bestAttributes.size();
		error(msg.str());
	}

  // ---------------------------------------------------------------------------
  PP->printLOG(Timestamp() + "Running RandomForest on [ " + int2str(bestAttributes.size()) + " ] attributes\n");
	double classifierError = 1.0;
  RandomForest* forest = new RandomForest(dataset, bestAttributes);
  forest->ComputeScores();
  classifierError = forest->GetClassificationError();
  if(par::verbose) {
    PP->printLOG(Timestamp() + "EvaporativeCooling::ComputeClassificationErrorRandomForest => " + dbl2str(classifierError) + "\n");
  }
  delete forest;
 
	/// return the classification error on this data
	return classifierError;
}
