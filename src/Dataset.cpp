/*
 * Dataset.cpp - Bill White - 6/14/05
 *
 * Collection class holding DatasetInstances
 * Added interaction information week of 4/18-26/06
 * Reworked entirely for McKinney Lab work - 2/28/11
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <utility> 
#include <random>

#include <limits.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>

#include <boost/lexical_cast.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <armadillo>

//#include <R.h>
//#include <Rcpp.h>
//#include <RInside.h>

#include "ChiSquared.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "StringUtils.h"
#include "Statistics.h"
#include "Insilico.h"
#include "DgeData.h"
#include "BirdseedData.h"
#include "DistanceMetrics.h"
#include "ReliefF.h"

#include "Insilico.h"
#include "helper.h"

extern Plink* PP;

using namespace std;
using namespace insilico;
using namespace boost;
using namespace boost::math;
using namespace arma;

Dataset::Dataset() {
	PP->printLOG(Timestamp() +  "Dataset: Default constructor\n");
	/// Set defaults.
	snpsFilename = "";
	hasGenotypes = false;
	hasAllelicInfo = false;

	numericsFilename = "";
	hasNumerics = false;

	alternatePhenotypesFilename = "";
	hasAlternatePhenotypes = false;
	hasPhenotypes = true;
	hasContinuousPhenotypes = false;
	classColumn = 0;

	maskIsPushed = false;

	/// Load attribute mutation map for transitions/transversions.
	attributeMutationMap[make_pair('A', 'G')] = TRANSITION_MUTATION;
	attributeMutationMap[make_pair('G', 'A')] = TRANSITION_MUTATION;
	attributeMutationMap[make_pair('C', 'T')] = TRANSITION_MUTATION;
	attributeMutationMap[make_pair('T', 'C')] = TRANSITION_MUTATION;
	attributeMutationMap[make_pair('A', 'C')] = TRANSVERSION_MUTATION;
	attributeMutationMap[make_pair('C', 'A')] = TRANSVERSION_MUTATION;
	attributeMutationMap[make_pair('G', 'T')] = TRANSVERSION_MUTATION;
	attributeMutationMap[make_pair('T', 'G')] = TRANSVERSION_MUTATION;
	attributeMutationMap[make_pair('A', 'T')] = TRANSVERSION_MUTATION;
	attributeMutationMap[make_pair('T', 'A')] = TRANSVERSION_MUTATION;
	attributeMutationMap[make_pair('C', 'G')] = TRANSVERSION_MUTATION;
	attributeMutationMap[make_pair('G', 'C')] = TRANSVERSION_MUTATION;
	attributeMutationMap[make_pair('X', 'Y')] = UNKNOWN_MUTATION;
	attributeMutationMap[make_pair('Y', 'X')] = UNKNOWN_MUTATION;

	/// metric defaults
	snpDiffMetricName = "gm";
	snpDiffFuncPtr = diffGMM;
	snpNearestNeighborMetricName = "gm";
	snpNearestNeighborFuncPtr = diffGMM;
	numDiffMetricName = "manhattan";
	numDiffFuncPtr = diffManhattan;

  if(par::verbose) {
    PP->printLOG(Timestamp() + 
                 "Default SNP diff metric: " + snpDiffMetricName + "\n");
    PP->printLOG(Timestamp() + 
                 "Default SNP nearest neighbors distance metric: " + 
                 snpNearestNeighborMetricName + "\n");
    PP->printLOG(Timestamp() + 
                 "Default continuous distance metric: " + 
                 numDiffMetricName + "\n");
  }
	PP->printLOG(Timestamp() +  "Dataset constructor done\n");
}

Dataset::~Dataset() {
	vector<DatasetInstance*>::const_iterator it;
	for (it = instances.begin(); it != instances.end(); it++) {
		if (*it) {
			delete *it;
		}
	}
}

bool Dataset::LoadDataset(vector<vector<int> >& dataMatrix,
		vector<int>& classLabels, vector<string>& attrNames) {

	cout << Timestamp() << "Loading Dataset from raw data vectors" << endl;

	uint numInstances = dataMatrix.size();
	uint numAttributes = dataMatrix[0].size();
	cout << Timestamp() << "Source data has " << numInstances << " rows and "
			<< numAttributes << " columns" << endl;

	if (!dataMatrix.size()) {
		cerr << "ERROR: LoadDataset from raw data. Data matrix is empty. "
				<< endl;
		return false;
	}
	if (!dataMatrix[0].size()) {
		cerr
				<< "ERROR: LoadDataset from raw data. Data matrix first row is empty. "
				<< endl;
		return false;
	}
	if (dataMatrix[0].size() != numAttributes) {
		cerr << "ERROR: LoadDataset from raw data. Number of attributes in the "
				<< "data matrix " << dataMatrix[0].size()
				<< " is not equal to the " << "number of attributes names"
				<< numAttributes << endl;
		return false;
	}
	if (classLabels.size() != numInstances) {
		cerr << "ERROR: LoadDataset from raw data. Number of class labels "
				<< classLabels.size()
				<< " is not equal to the number of instances " << numInstances
				<< endl;
		return false;
	}

	for (uint i = 0; i < numAttributes; ++i) {
		attributeNames.push_back(attrNames[i]);
		attributesMask[attrNames[i]] = i;
	}

	cout << Timestamp() << attributeNames.size() << " attribute names read"
			<< endl;

	levelCounts.resize(numAttributes);
	levelCountsByClass.resize(numAttributes);
	attributeLevelsSeen.resize(numAttributes);
	genotypeCounts.resize(numAttributes);

	uint rowIndex = 0;
	vector<vector<int> >::const_iterator rowIt = dataMatrix.begin();
	for (; rowIt != dataMatrix.end(); ++rowIt, ++rowIndex) {
		vector<int> row(numAttributes);
		copy(rowIt->begin(), rowIt->end(), row.begin());
		string ID = zeroPadNumber(rowIndex, 8) + zeroPadNumber(rowIndex, 8);
		DatasetInstance* dsi = new DatasetInstance(this, ID);
		dsi->LoadInstanceFromVector(row);
		dsi->SetClass(classLabels[rowIndex]);
		classIndexes[classLabels[rowIndex]].push_back(rowIndex);
		instances.push_back(dsi);
		instanceIds.push_back(ID);
		instancesMask[ID] = rowIndex;
	}

	hasGenotypes = true;
	hasNumerics = false;
	hasContinuousPhenotypes = false;
	hasPhenotypes = true;

	UpdateAllLevelCounts();
	CreateDummyAlleles();
	PrintStatsSimple();
	// PrintLevelCounts();

	return true;
}

bool Dataset::LoadDataset(string snpsFilename, string numericsFilename,
		string altPhenoFilename, vector<string> ids) {

	// instance IDs to use when loading from data, numeric and
	// alternate phenotype files
	if (ids.size()) {
		instanceIdsToLoad.resize(ids.size());
		copy(ids.begin(), ids.end(), instanceIdsToLoad.begin());
		if (numericsFilename != "") {
			hasNumerics = true;
		}
		if (altPhenoFilename != "") {
			hasAlternatePhenotypes = true;
		}
	}

	// load SNPs
	if (snpsFilename != "") {
		if (!LoadSnps(snpsFilename)) {
			cerr << "ERROR in LoadDataset. Could not load SNPs file" << endl;
			return false;
		}
		if (instancesMask.size() == 0) {
			cerr << "ERROR: No instances for analysis" << endl;
			return false;
		}
		if ((attributesMask.size() == 0) && (numericsMask.size() == 0)) {
			cerr << "ERROR: No variables for analysis" << endl;
			return false;
		}
		if (attributesMask.size()) {
			hasGenotypes = true;
		}
		if (numericsMask.size()) {
			hasNumerics = true;
		}

		// remove the instances that don't match the covariate and/or phenotype files
		uint nextInstanceIndex = 0;
		if (instanceIdsToLoad.size()) {
			cout << Timestamp() << "Finding matching IDs" << endl;
			vector<DatasetInstance*> tempInstances;
			vector<string> tempInstanceIds;
			map<string, uint> tempInstancesMask;
			for (uint i = 0; i < instanceIdsToLoad.size(); ++i) {
				uint instanceIndex = 0;
				string ID = instanceIdsToLoad[i];
				if (!GetInstanceIndexForID(ID, instanceIndex)) {
					cerr << "ERROR: Could not find ID in data set: " << ID
							<< endl;
					return false;
				}
				tempInstances.push_back(instances[instanceIndex]);
				tempInstanceIds.push_back(ID);
				tempInstancesMask[ID] = nextInstanceIndex;
				++nextInstanceIndex;
			}
			instances = tempInstances;
			instanceIds = tempInstanceIds;
			instancesMask = tempInstancesMask;
		}
		cout << Timestamp() << NumInstances()
				<< " instances remain after covariate/phenotype matching"
				<< endl;
	}

	// load numerics
	if (numericsFilename != "") {
		if (!LoadNumerics(numericsFilename)) {
			cerr << "ERROR in LoadDataset. Could not load numerics file"
					<< endl;
			return false;
		}
	}

	// load alternate phenotypes
	if (altPhenoFilename != "") {
		if (!LoadAlternatePhenotypes(altPhenoFilename)) {
			cerr
					<< "ERROR in LoadDataset. Could not load alternate phenotypes file"
					<< endl;
			return false;
		}
	}

	hasPhenotypes = true;

	return true;
}

bool Dataset::LoadDataset(DgeData* dgeData) {
	// TODO: check for genotypes already loaded; assume DGE only for now

	numericsFilename = dgeData->GetCountsFilename();
	cout << Timestamp() << "Reading numerics from DgeData class whose "
			<< "original counts are from ["	<< numericsFilename << "]" << endl;

	// populate numericNames, numericsMinMax and numericsMask
	vector<string> geneNames = dgeData->GetGeneNames();
	for (int geneIndex = 0; geneIndex < (int) geneNames.size(); ++geneIndex) {
		numericsNames.push_back(geneNames[geneIndex]);
		numericsMinMax.push_back(dgeData->GetGeneMinMax(geneIndex));
		numericsSums.push_back(dgeData->GetGeneCountsSum(geneIndex));
		numericsMask[geneNames[geneIndex]] = geneIndex;
	}

	// load the data set instances: set the instance numerics,
	// instance IDs, instance mask and phenotype
	vector<string> sampleNames = dgeData->GetSampleNames();
	for (int instanceIndex = 0; instanceIndex < (int) sampleNames.size();
			++instanceIndex) {
		string ID = sampleNames[instanceIndex] + sampleNames[instanceIndex];
		vector<double> sampleValues = dgeData->GetSampleCounts(instanceIndex);
		DatasetInstance* dsi = new DatasetInstance(this, ID);
		for (int numericIndex = 0; numericIndex < (int) geneNames.size();
				++numericIndex) {
			dsi->AddNumeric(sampleValues[numericIndex]);
		}
		instances.push_back(dsi);
		instanceIds.push_back(ID);
		numericsIds.push_back(ID);
		instancesMask[ID] = instanceIndex;
		ClassLevel thisClass = dgeData->GetSamplePhenotype(instanceIndex);
		dsi->SetClass(thisClass);
		classIndexes[thisClass].push_back(instanceIndex);
	}
	hasPhenotypes = true;

	cout << Timestamp() << "Read " << NumNumerics() << " numeric attributes"
			<< " from DGE counts data" << endl;

	hasNumerics = true;
	hasAllelicInfo = false;

	return true;
}

bool Dataset::LoadDataset(BirdseedData* birdseedData) {
	// TODO: check for genotypes already loaded; assume birdseed only for now

	snpsFilename = "BirdseedData CLASS";

	cout << Timestamp() << "Reading attributes from " << snpsFilename << endl;

	// populate attributeNames and attributesMask
	vector<string> snpNames = birdseedData->GetSNPNames();
	int numAttributes = snpNames.size();
	for (int i = 0; i < numAttributes; ++i) {
		attributeNames.push_back(snpNames[i]);
		attributesMask[snpNames[i]] = i;
	}
	levelCounts.resize(numAttributes);
	levelCountsByClass.resize(numAttributes);
	attributeLevelsSeen.resize(numAttributes);
	genotypeCounts.resize(numAttributes);
	attributeAlleles.resize(numAttributes);
	attributeAlleleCounts.resize(numAttributes);
	attributeMinorAllele.resize(numAttributes);
	attributeMutationTypes.resize(numAttributes);

	// load the data set instances: set the instance attributes,
	// instance IDs, instance mask and phenotype
	vector<string> sampleNames;
	if (birdseedData->HasSubjectLabels()) {
		sampleNames = birdseedData->GetSubjectLabels();
	} else {
		sampleNames = birdseedData->GetSubjectNames();
	}
	int numSamples = sampleNames.size();
	for (int instanceIndex = 0; instanceIndex < numSamples; ++instanceIndex) {
		string ID = sampleNames[instanceIndex] + sampleNames[instanceIndex];
		vector<AttributeLevel> sampleValues = birdseedData->GetSubjectGenotypes(
				instanceIndex);
		DatasetInstance* dsi = new DatasetInstance(this, ID);
		vector<AttributeLevel> thisSnps(numAttributes);
		for (int snpIndex = 0; snpIndex < numAttributes; ++snpIndex) {
			AttributeLevel thisSnp = sampleValues[snpIndex];
			// attributeLevelsSeen[snpIndex].insert(thisSnpString);
			thisSnps[snpIndex] = thisSnp;
		}
    dsi->LoadInstanceFromVector(thisSnps);
		instances.push_back(dsi);

		instanceIds.push_back(ID);
//			attributeIds.push_back(ID);
		instancesMask[ID] = instanceIndex;

		vector<uint> bsMissingValues;
		bool hasMissingValues = birdseedData->GetMissingValues(ID,
				bsMissingValues);
		if (hasMissingValues) {
			missingValues[ID] = bsMissingValues;
		}

		if (birdseedData->HasPhenotypes()) {
			ClassLevel thisClass = birdseedData->GetSamplePhenotype(
					instanceIndex);
			dsi->SetClass(thisClass);
			classIndexes[thisClass].push_back(instanceIndex);
		} else {
			dsi->SetClass(MISSING_DISCRETE_CLASS_VALUE);
		}
	}

	for (int snpIndex = 0; snpIndex < numAttributes; ++snpIndex) {
		/// add allelic info
		pair<char, char> thisSnpAlleles = birdseedData->GetMajorMinorAlleles(
				snpIndex);
		double thisSnpMajAlleleFreq = birdseedData->GetMajorAlleleFrequency(
				snpIndex);
		attributeAlleles[snpIndex] = thisSnpAlleles;
		attributeMinorAllele[snpIndex] = make_pair(thisSnpAlleles.second,
				1.0 - thisSnpMajAlleleFreq);
		++attributeAlleleCounts[snpIndex][thisSnpAlleles.first];
		++attributeAlleleCounts[snpIndex][thisSnpAlleles.second];
		genotypeCounts[snpIndex] = birdseedData->GetGenotypeCounts(snpIndex);
		attributeMutationTypes[snpIndex] = attributeMutationMap[make_pair(
				thisSnpAlleles.second, thisSnpAlleles.first)];
	}

	cout << Timestamp() << "Dataset: Read " << NumInstances() << " instances, "
			<< NumAttributes() << " SNP attributes" << " from birdseed data"
			<< endl;

	if (!birdseedData->HasPhenotypes()) {
		cout << Timestamp() << "Birdseed data DOES NOT have phenotypes."
				<< endl;
		hasPhenotypes = false;
	} else {
		hasPhenotypes = true;
	}

	hasGenotypes = true;
	hasAllelicInfo = true;

	// birdseedData->PrintAlleleCounts();
	UpdateAllLevelCounts();
	// PrintLevelCounts();

	return true;
}

bool Dataset::LoadOtherDatasetInstances(Dataset* otherDs, 
                                        vector<uint> instIdx) {
  // --------------------------------------------------------------------------
  // Create new instance meta data in this Dataset from meta data and instances 
  // from the other data instances in passed vector instIdx as indexes
  instances.clear();
  instanceIds.clear();
  instanceIdsToLoad.clear();
  instancesMask.clear();
  for(uint newInstanceIdx=0; newInstanceIdx < instIdx.size(); ++newInstanceIdx) {
    uint otherInstanceIdx = instIdx[newInstanceIdx];
    DatasetInstance* srcInstance = otherDs->GetInstance(otherInstanceIdx);
    string destInstanceId = srcInstance->GetId();
    DatasetInstance* dstInstance = new DatasetInstance(this, destInstanceId);
    dstInstance->LoadInstanceFromInstancePtr(this, srcInstance);
    instanceIds.push_back(destInstanceId);
    classIndexes[dstInstance->GetClass()].push_back(newInstanceIdx);
    instanceIdsToLoad.push_back(destInstanceId);
    instances.push_back(dstInstance);
    instancesMask[destInstanceId] = newInstanceIdx;
  }
  hasPhenotypes = otherDs->HasPhenotypes();
  hasAlternatePhenotypes = true;
  hasContinuousPhenotypes = otherDs->HasContinuousPhenotypes();
  if(!hasContinuousPhenotypes) {
    this->PrintClassIndexInfo(cout);
  }
  // --------------------------------------------------------------------------
	hasGenotypes = otherDs->HasGenotypes();
  if(hasGenotypes) {
    UpdateAllLevelCounts();
  	hasAllelicInfo = otherDs->HasAllelicInfo();
    // TODO: how to copy allelic info???
    MaskIncludeAllAttributes(DISCRETE_TYPE);
  }
  hasNumerics = otherDs->HasNumerics();
  if(hasNumerics) {
    numericsNames = otherDs->GetNumericsNames();
    MaskIncludeAllAttributes(NUMERIC_TYPE);
    // find the min and max values for each numeric attribute
    // used in diff/distance calculation metrics
    vector<NumericLevel> numericColumn;
    for (uint i = 0; i < NumNumerics(); ++i) {
      double columnSum = 0.0;
      GetNumericValues(i, numericColumn);
      double minElement = *numericColumn.begin();
      double maxElement = *numericColumn.begin();
      for (vector<NumericLevel>::const_iterator it = numericColumn.begin();
          it != numericColumn.end(); ++it) {
        if ((*it != MISSING_NUMERIC_VALUE) && (*it < minElement)) {
          minElement = *it;
        }
        if ((*it != MISSING_NUMERIC_VALUE) && (*it > maxElement)) {
          maxElement = *it;
        }
        columnSum += *it;
      }
      numericsMinMax.push_back(make_pair(minElement, maxElement));
      numericsSums.push_back(columnSum);
    }
  }
    
  return true;
}

bool Dataset::GetAttributeRowCol(uint row, uint col,
		AttributeLevel& attrVal) {
	unsigned long numInstances = instances.size();
	unsigned long numAttributes = instances[0]->NumAttributes();

	if ((row < numInstances) && (col < numAttributes)) {
		attrVal = instances[row]->GetAttribute(col);
		return true;
	}

	return false;
}

bool Dataset::GetNumericRowCol(uint row, uint col,
		NumericLevel& numVal) {
	unsigned long numInstances = instances.size();
	unsigned long numNumerics = instances[0]->NumNumerics();

	if ((row < numInstances) && (col < numNumerics)) {
		numVal = instances[row]->GetNumeric(col);
		return true;
	}

	return false;
}

bool Dataset::WriteNewDataset(string newDatasetFilename,
		OutputDatasetType outputDatasetType) {

	switch (outputDatasetType) {
	case PLINK_PED_DATASET:
		if (!WriteNewPlinkPedDataset(GetFileBasename(newDatasetFilename))) {
			cerr << "ERROR: Failed to write new PLINK ped data set";
			return false;
		} else {
			return true;
		}
	case PLINK_COVAR_DATASET:
		if (!WriteNewPlinkCovarDataset(GetFileBasename(newDatasetFilename))) {
			cerr << "ERROR: Failed to write new PLINK ped data set";
			return false;
		} else {
			return true;
		}
	case TAB_DELIMITED_DATASET:
	case CSV_DELIMITED_DATASET:
	case ARFF_DATASET:
	case PLINK_BED_DATASET:
	case NO_OUTPUT_DATASET:
		break;
	default:
		cerr << "ERROR: Output data set file type not recognized." << endl;
		return false;
	}

	ofstream newDatasetStream(newDatasetFilename.c_str());
	if (!newDatasetStream.is_open()) {
		cerr << "ERROR: Could not open new dataset file: " << newDatasetFilename
				<< endl;
		return false;
	}

	/// write the attribute names header
	if (outputDatasetType == ARFF_DATASET) {
		newDatasetStream << "@RELATION dataset" << endl << endl;
	}
	map<string, uint>::const_iterator ait = attributesMask.begin();
	for (; ait != attributesMask.end(); ++ait) {
		switch (outputDatasetType) {
		case TAB_DELIMITED_DATASET:
			newDatasetStream << ait->first << "\t";
			break;
		case CSV_DELIMITED_DATASET:
			newDatasetStream << ait->first << ",";
			break;
		case ARFF_DATASET:
			newDatasetStream << "@ATTRIBUTE " << ait->first << " {0,1,2}"
					<< endl;
			break;
		case PLINK_PED_DATASET:
		case PLINK_BED_DATASET:
		case NO_OUTPUT_DATASET:
		default:
			cerr << "ERROR: Unrecognized output data set type: "
					<< outputDatasetType << endl;
			return false;
		}
	}
	map<string, uint>::const_iterator nit = numericsMask.begin();
	for (; nit != numericsMask.end(); ++nit) {
		switch (outputDatasetType) {
		case TAB_DELIMITED_DATASET:
			newDatasetStream << nit->first << "\t";
			break;
		case CSV_DELIMITED_DATASET:
			newDatasetStream << nit->first << ",";
			break;
		case ARFF_DATASET:
			newDatasetStream << "@ATTRIBUTE " << nit->first << " numeric"
					<< endl;
			break;
		case PLINK_PED_DATASET:
		case PLINK_BED_DATASET:
			break;
		case NO_OUTPUT_DATASET:
		default:
			cerr << "ERROR: Unrecognized output data set type: "
					<< outputDatasetType << endl;
			return false;
		}
	}
	switch (outputDatasetType) {
	case TAB_DELIMITED_DATASET:
	case CSV_DELIMITED_DATASET:
		newDatasetStream << "Class" << endl;
		break;
	case ARFF_DATASET:
		if (hasContinuousPhenotypes) {
			newDatasetStream << "@ATTRIBUTE Class numeric" << endl;
		} else {
			newDatasetStream << "@ATTRIBUTE Class {0,1}" << endl;
		}
		break;
	case PLINK_PED_DATASET:
	case PLINK_BED_DATASET:
	case NO_OUTPUT_DATASET:
	default:
		cerr << "ERROR: Unrecognized output data set type: "
				<< outputDatasetType << endl;
		return false;
	}

	/// write the data, respecting the masked attributes, numerics
	/// and masked instances - 10/28/11
	/// write the attribute names header
	if (outputDatasetType == ARFF_DATASET) {
		newDatasetStream << endl << "@DATA" << endl;
	}
	vector<string> instanceIds = GetInstanceIds();
	for (uint iIdx = 0; iIdx < NumInstances(); iIdx++) {
		unsigned instanceIndex = 0;
		GetInstanceIndexForID(instanceIds[iIdx], instanceIndex);
		// write discrete attribute values
		vector<uint> attrIndices = MaskGetAttributeIndices(
				DISCRETE_TYPE);
		for (uint aIdx = 0; aIdx < attrIndices.size(); aIdx++) {
			AttributeLevel A = instances[instanceIndex]->GetAttribute(
					attrIndices[aIdx]);
			switch (outputDatasetType) {
			case TAB_DELIMITED_DATASET:
				if (A == MISSING_ATTRIBUTE_VALUE) {
					newDatasetStream << "?" << "\t";
				} else {
					newDatasetStream << A << "\t";
				}
				break;
			case CSV_DELIMITED_DATASET:
			case ARFF_DATASET:
				if (A == MISSING_ATTRIBUTE_VALUE) {
					newDatasetStream << "?" << ",";
				} else {
					newDatasetStream << A << ",";
				}
				break;
			case PLINK_PED_DATASET:
			case PLINK_BED_DATASET:
			case NO_OUTPUT_DATASET:
			default:
				cerr << "ERROR: Unrecognized output data set type: "
						<< outputDatasetType << endl;
				return false;
			}
		}
		/// write continuous attribute values
		vector<uint> numIndices = MaskGetAttributeIndices(NUMERIC_TYPE);
		for (uint nIdx = 0; nIdx < numIndices.size(); nIdx++) {
			NumericLevel N = instances[instanceIndex]->GetNumeric(
					numIndices[nIdx]);
			switch (outputDatasetType) {
			case TAB_DELIMITED_DATASET:
				if (N == MISSING_NUMERIC_VALUE) {
					newDatasetStream << "?" << "\t";
				} else {
					newDatasetStream << N << "\t";
				}
				break;
			case CSV_DELIMITED_DATASET:
			case ARFF_DATASET:
				if (N == MISSING_NUMERIC_VALUE) {
					newDatasetStream << "?" << ",";
				} else {
					newDatasetStream << N << ",";
				}
				break;
			case PLINK_PED_DATASET:
			case PLINK_BED_DATASET:
			case NO_OUTPUT_DATASET:
			default:
				cerr << "ERROR: Unrecognized output data set type: "
						<< outputDatasetType << endl;
				return false;
			}
		}
		// class/phenotype
		ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
		NumericLevel continuousClassLevel = MISSING_NUMERIC_CLASS_VALUE;
		if (hasContinuousPhenotypes) {
			continuousClassLevel =
					instances[instanceIndex]->GetPredictedValueTau();
			newDatasetStream << continuousClassLevel << endl;
		} else {
			discreteClassLevel = instances[instanceIndex]->GetClass();
			newDatasetStream << discreteClassLevel << endl;
		}
	}

	newDatasetStream.close();

	return true;
}

bool Dataset::WriteNewDataset(string newDatasetFilename,
		vector<string> attributes, OutputDatasetType outputDatasetType) {

	cout << Timestamp() << "Dataset::WriteNewDataset with " << attributes.size()
			<< " named attributes" << endl;
	PrintStats();

	if (!attributes.size()) {
		error("Dataset::WriteNewDataset: no named attributes to write\n");
	}

	switch (outputDatasetType) {
	case PLINK_PED_DATASET:
		error("This version of WriteNewDataset does not support PED\n");
	case TAB_DELIMITED_DATASET:
	case CSV_DELIMITED_DATASET:
	case ARFF_DATASET:
	case PLINK_BED_DATASET:
	case NO_OUTPUT_DATASET:
		break;
	default:
		error("Output data set file type not recognized.\n");
		return false;
	}

	ofstream newDatasetStream(newDatasetFilename.c_str());
	if (!newDatasetStream.is_open()) {
		error("Could not open new dataset file: " + newDatasetFilename + "\n");
	}

	/// write the attribute names header
	if (outputDatasetType == ARFF_DATASET) {
		newDatasetStream << "@RELATION dataset" << endl << endl;
	}
	map<string, uint>::const_iterator ait = attributesMask.begin();
	for (; ait != attributesMask.end(); ++ait) {
		/// is this attribute in the list passed in as a parameter
		if (find(attributes.begin(), attributes.end(), ait->first)
				== attributes.end()) {
			continue;
		}
		switch (outputDatasetType) {
		case TAB_DELIMITED_DATASET:
			newDatasetStream << ait->first << "\t";
			break;
		case CSV_DELIMITED_DATASET:
			newDatasetStream << ait->first << ",";
			break;
		case ARFF_DATASET:
			newDatasetStream << "@ATTRIBUTE " << ait->first << " {0,1,2}"
					<< endl;
			break;
		case PLINK_PED_DATASET:
		case PLINK_BED_DATASET:
		case NO_OUTPUT_DATASET:
		default:
			error("Unrecognized output data set type: " + 
            int2str(outputDatasetType) + "\n");
		}
	}
	map<string, uint>::const_iterator nit = numericsMask.begin();
	for (; nit != numericsMask.end(); ++nit) {
		if (find(attributes.begin(), attributes.end(), nit->first)
				== attributes.end()) {
			continue;
		}
  	switch (outputDatasetType) {
      case TAB_DELIMITED_DATASET:
        newDatasetStream << nit->first << "\t";
        break;
      case CSV_DELIMITED_DATASET:
        newDatasetStream << nit->first << ",";
        break;
      case ARFF_DATASET:
        newDatasetStream << "@ATTRIBUTE " << nit->first << " numeric" << endl;
        break;
      case PLINK_PED_DATASET:
      case PLINK_BED_DATASET:
      case NO_OUTPUT_DATASET:
      default:
        error("Unrecognized output data set type: " + 
              int2str(outputDatasetType) + "\n");
    }
	}
	switch (outputDatasetType) {
    case TAB_DELIMITED_DATASET:
    case CSV_DELIMITED_DATASET:
      newDatasetStream << "Class" << endl;
      break;
    case ARFF_DATASET:
      if (hasContinuousPhenotypes) {
        newDatasetStream << "@ATTRIBUTE Class numeric" << endl;
      } else {
        newDatasetStream << "@ATTRIBUTE Class {0,1}" << endl;
      }
      break;
    case PLINK_PED_DATASET:
    case PLINK_BED_DATASET:
    case NO_OUTPUT_DATASET:
    default:
      error("Unrecognized output data set type: " + 
            int2str(outputDatasetType) + "\n");
	}

	/// write the data, respecting the masked attributes, numerics
	/// and masked instances - 10/28/11
	/// write the attribute names header
	if (outputDatasetType == ARFF_DATASET) {
		newDatasetStream << endl << "@DATA" << endl;
	}
	vector<string> instanceIds = GetInstanceIds();
	for (uint iIdx = 0; iIdx < NumInstances(); iIdx++) {
		unsigned instanceIndex = 0;
		GetInstanceIndexForID(instanceIds[iIdx], instanceIndex);
		// write discrete attribute values
		vector<uint> attrIndices = MaskGetAttributeIndices(
				DISCRETE_TYPE);
		for (uint aIdx = 0; aIdx < attrIndices.size(); aIdx++) {
			AttributeLevel A = instances[instanceIndex]->GetAttribute(
					attrIndices[aIdx]);
			if (find(attributes.begin(), attributes.end(),
					attributeNames[attrIndices[aIdx]]) == attributes.end()) {
				continue;
			}
			switch (outputDatasetType) {
			case TAB_DELIMITED_DATASET:
				if (A == MISSING_ATTRIBUTE_VALUE) {
					newDatasetStream << "?" << "\t";
				} else {
					newDatasetStream << A << "\t";
				}
				break;
			case CSV_DELIMITED_DATASET:
			case ARFF_DATASET:
				if (A == MISSING_ATTRIBUTE_VALUE) {
					newDatasetStream << "?" << ",";
				} else {
					newDatasetStream << A << ",";
				}
				break;
			case PLINK_PED_DATASET:
			case PLINK_BED_DATASET:
			case NO_OUTPUT_DATASET:
			default:
				error("Unrecognized output data set type: " + 
              int2str(outputDatasetType) + "\n");
			}
		}
		/// write continuous attribute values
		vector<uint> numIndices = MaskGetAttributeIndices(NUMERIC_TYPE);
		for (uint nIdx = 0; nIdx < numIndices.size(); nIdx++) {
			NumericLevel N = instances[instanceIndex]->GetNumeric(
					numIndices[nIdx]);
			if (find(attributes.begin(), attributes.end(),
					numericsNames[numIndices[nIdx]]) == attributes.end()) {
				continue;
			}
			switch (outputDatasetType) {
			case TAB_DELIMITED_DATASET:
				if (N == MISSING_NUMERIC_VALUE) {
					newDatasetStream << "?" << "\t";
				} else {
					newDatasetStream << N << "\t";
				}
				break;
			case CSV_DELIMITED_DATASET:
			case ARFF_DATASET:
				if (N == MISSING_NUMERIC_VALUE) {
					newDatasetStream << "?" << ",";
				} else {
					newDatasetStream << N << ",";
				}
				break;
			case PLINK_PED_DATASET:
			case PLINK_BED_DATASET:
			case NO_OUTPUT_DATASET:
			default:
				error("ERROR: Unrecognized output data set type: " +
              int2str(outputDatasetType) + "\n");
			}
		}
		// class/phenotype
		ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
		NumericLevel continuousClassLevel = MISSING_NUMERIC_CLASS_VALUE;
		if (hasContinuousPhenotypes) {
			continuousClassLevel =
					instances[instanceIndex]->GetPredictedValueTau();
			newDatasetStream << continuousClassLevel << endl;
		} else {
			discreteClassLevel = instances[instanceIndex]->GetClass();
			newDatasetStream << discreteClassLevel << endl;
		}
	}

	newDatasetStream.close();

	return true;
}

bool Dataset::ExtractAttributes(string scoresFilename, uint topN,
		string newDatasetFilename) {
	// read attribute scores from file
	ifstream scoresStream(scoresFilename.c_str());
	if (!scoresStream.is_open()) {
		error("Could not open scores file: " + scoresFilename + "\n");
	}
	// build a "rank map": pair<ranker score, attribute index>
	// for all scores in the score file
	string line;
	vector<pair<double, string> > rankMap;
	vector<string> tokens;
	while (getline(scoresStream, line)) {
//    cout << line << endl;
		tokens.clear();
		split(tokens, line);
//    cout << "t0: " << tokens[0] << ", t1: " << tokens[1] << endl;
		rankMap.push_back(
				make_pair(lexical_cast<double>(tokens[0]), tokens[1]));
	}
	scoresStream.close();

	// sort the rank map in ascending score order
	sort(rankMap.begin(), rankMap.end());

	// print the top N sorted index
//  vector<pair<double, string> >::const_iterator rankMapIt;
//  for(rankMapIt = rankMap.end() - 1;
//      rankMapIt != (rankMap.end() - topN - 1);
//      rankMapIt--) {
//    cout << rankMapIt->first << " => " << rankMapIt->second << endl;
//  }

	// write top N ranked attributes to new file
	// honor MDR format by using header row information - bcw - 9/13/05
	ofstream newDatasetStream(newDatasetFilename.c_str());
	if (!newDatasetStream.is_open()) {
		error("Could not open new data set file: " + newDatasetFilename + "\n");
	}

	// write MDR header - bcw - 9/13/05
	vector<pair<double, string> >::const_iterator rankHdrIt;
	for (rankHdrIt = rankMap.end() - 1; rankHdrIt != rankMap.end() - topN - 1;
			rankHdrIt--) {
		newDatasetStream << rankHdrIt->second << "\t";
	}
	newDatasetStream << "Class" << endl;

	// write dataset instances
	// MDR format -> tab-delimited attributes + class - bcw - 9/13/05
	// adjusted for integrated data and continuous class
	vector<DatasetInstance*>::const_iterator dsiIt;
	vector<pair<double, string> >::const_iterator rankIt;
	for (dsiIt = instances.begin(); dsiIt != instances.end(); dsiIt++) {
		for (rankIt = rankMap.end() - 1; rankIt != rankMap.end() - topN - 1;
				rankIt--) {
			if (MaskSearchVariableType(rankIt->second, DISCRETE_TYPE)) {
				newDatasetStream
						<< (*dsiIt)->GetAttribute(
								GetAttributeIndexFromName(rankIt->second))
						<< "\t";
			} else {
				if (MaskSearchVariableType(rankIt->second, NUMERIC_TYPE)) {
					newDatasetStream
							<< (*dsiIt)->GetNumeric(
									GetNumericIndexFromName(rankIt->second))
							<< "\t";
				} else {
					error("Dataset::ExtractAttributes: requested variable name [" +
                rankIt->second + "] is not in the data set" + "\n");
				}
			}
		}
		if (hasContinuousPhenotypes) {
			newDatasetStream << (*dsiIt)->GetPredictedValueTau() << endl;
		} else {
			newDatasetStream << (*dsiIt)->GetClass() << endl;
		}
	}

	newDatasetStream.close();

	return true;
}

bool Dataset::SwapAttributes(uint a1, uint a2) {
	vector<DatasetInstance*>::const_iterator it;
	for (it = instances.begin(); it != instances.end(); it++) {
		(*it)->SwapAttributes(a1, a2);
	}

	return true;
}

uint Dataset::NumVariables() {
	return (NumAttributes() + NumNumerics());
}

vector<string> Dataset::GetVariableNames() {
	vector<string> variableNames(attributesMask.size() + numericsMask.size());
	vector<string> attrNames = GetAttributeNames();
	vector<string> numNames = GetNumericsNames();

	copy(attrNames.begin(), attrNames.end(), variableNames.begin());
	copy(numNames.begin(), numNames.end(),
			variableNames.begin() + attrNames.size());

	return variableNames;
}

uint Dataset::NumInstances() {
	return instancesMask.size();
}

DatasetInstance* Dataset::GetInstance(uint index) {
	if (index < instances.size()) {
		return instances[index];
	}
	error("Instance index out of range: " + int2str(index) + "\n");
}

DatasetInstance* Dataset::GetRandomInstance() {
  uniform_int_distribution<int> runif(0, NumInstances() - 1);
	return instances[runif(engine)];
}

vector<string> Dataset::GetInstanceIds() {
	vector<string> idsToReturn;
	map<string, uint>::const_iterator it = instancesMask.begin();
	for (; it != instancesMask.end(); ++it) {
		idsToReturn.push_back(it->first);
	}
	return idsToReturn;
}

bool Dataset::GetInstanceIndexForID(string ID, uint& instanceIndex) {

	if (instancesMask.find(ID) != instancesMask.end()) {
		instanceIndex = instancesMask[ID];
		return true;
	}
	error("Could not find instance ID: " + ID + "\n");

	return false;
}

uint Dataset::NumAttributes() {
	return attributesMask.size();
}

vector<string> Dataset::GetAttributeNames() {
	vector<string> names;
	map<string, uint>::const_iterator it = attributesMask.begin();
	for (; it != attributesMask.end(); ++it) {
		names.push_back(it->first);
	}
	return names;
}

vector<string> Dataset::GetFileAttributeNames() {
	vector<string> names;
	vector<string>::const_iterator it = attributeNames.begin();
	for (; it != attributeNames.end(); ++it) {
		if (MaskSearchVariableType(*it, DISCRETE_TYPE)) {
			names.push_back(*it);
		}
	}
	return names;
}

bool Dataset::GetAttributeValues(uint attributeIndex,
		vector<AttributeLevel>& attributeValues) {
	if (attributeIndex > attributeNames.size()) {
		error("Dataset::GetAttributeValues: attribute index out of range: "
				  + int2str(attributeIndex) + "\n");
		return false;
	}
	if (hasGenotypes) {
		attributeValues.clear();
		map<string, uint>::const_iterator it;
		for (it = instancesMask.begin(); it != instancesMask.end(); it++) {
			AttributeLevel thisAttribute = instances[it->second]->GetAttribute(
					attributeIndex);
			attributeValues.push_back(thisAttribute);
		}
	} else {
		error("Attempting to access SNP data when none have been loaded\n");
		return false;
	}

	return true;
}

bool Dataset::GetAttributeValues(string attributeName,
		vector<AttributeLevel>& attributeValues) {
	if (MaskSearchVariableType(attributeName, DISCRETE_TYPE)) {
		GetAttributeValues(GetAttributeIndexFromName(attributeName),
				attributeValues);
	} else {
		error("Dataset::GetAttributeValues cannot get attribute values for: " + 
          attributeName + ". Either doesn't exist or is excluded\n");
		return false;
	}
	return true;
}

std::string Dataset::GetSnpsFilename() {
	return snpsFilename;
}

bool Dataset::HasGenotypes() {
	return hasGenotypes;
}

bool Dataset::HasAllelicInfo() {
	return hasAllelicInfo;
}

AttributeLevel Dataset::GetAttribute(unsigned instanceIndex, string name) {
	if (instanceIndex >= instances.size()) {
		error("Dataset::GetAttribute: instance index " + 
          int2str(instanceIndex) + " out of range\n");
	}
	map<string, uint>::iterator pos = attributesMask.find(name);
	if (pos != attributesMask.end()) {
		return instances[instanceIndex]->GetAttribute(pos->second);
	} else {
		error("Dataset::GetAttribute: " + name + " at instance index: " + 
          int2str(instanceIndex) + " not found\n");
	}
}

pair<char, char> Dataset::GetAttributeAlleles(uint attributeIndex) {
	if (!hasAllelicInfo) {
		error("GetAttributeAlleles: This data set does not have "
				  "allelic information.\n");
	}
	pair<char, char> returnPair = make_pair(' ', ' ');
	if (attributeIndex < attributeAlleles.size()) {
		returnPair = attributeAlleles[attributeIndex];
	} else {
		error("GetAttributeAlleles: attribute index out of range: " +
				  int2str(attributeIndex) + "\n");
	}
	return returnPair;
}

pair<char, double> Dataset::GetAttributeMAF(uint attributeIndex) {
	/// An Introduction to Genetic Analysis by Griffiths, Miller, Suzuki,
	/// Lewontin and Gelbart, 2000, page 715.
	pair<char, double> returnPair = make_pair(' ', 0.0);
	double p = 0.0;
	double q = 0.0;

	map<AttributeLevel, int> genoCounts;
	if (attributeIndex < NumAttributes()) {
		vector<AttributeLevel> thisAttrCol;
		GetAttributeValues(attributeIndex, thisAttrCol);
		for (uint i = 0; i < thisAttrCol.size(); ++i) {
			++genoCounts[thisAttrCol[i]];
		}
		if (genoCounts.size() == 3) {
			double f_AA = ((double) genoCounts[0]) / NumInstances();
			double f_Aa = ((double) genoCounts[1]) / NumInstances();
			double f_aa = ((double) genoCounts[2]) / NumInstances();
			p = f_AA + 0.5 * f_Aa;
			q = f_aa + 0.5 * f_Aa;
			if (p < q) {
				returnPair = make_pair('p', p);
			} else {
				returnPair = make_pair('q', q);
			}
		}
	}
	return returnPair;
}

bool Dataset::ProcessExclusionFile(string exclusionFilename) {
	ifstream dataStream(exclusionFilename.c_str());
	if (!dataStream.is_open()) {
		error("Could not open exclusion file: " + exclusionFilename + "\n");
	}

	// temporary string for reading file lines
	string line;
	uint lineNumber = 0;
	while (getline(dataStream, line)) {
		++lineNumber;
		string attributeName = trim(line);
		if (!MaskRemoveVariable(attributeName)) {
			error("Attribute to exclude [" + attributeName + "] on line [" +
            int2str(lineNumber) + 
            "] in the exclusion file. It is not in the data set\n");
		}
	}
	dataStream.close();

	return true;
}

AttributeMutationType Dataset::GetAttributeMutationType(
		uint attributeIndex) {
	if (attributeIndex < attributeMutationTypes.size()) {
		return attributeMutationTypes[attributeIndex];
	}
	return UNKNOWN_MUTATION;
}

double Dataset::GetJukesCantorDistance(DatasetInstance* dsi1,
		DatasetInstance* dsi2) {
	double returnValue = 0.0;

	double d = 0.0;
	vector<uint> attributeIndicies = MaskGetAttributeIndices(
			DISCRETE_TYPE);
	for (uint attrIdx = 0; attrIdx < attributeIndicies.size();
			++attrIdx) {
		AttributeLevel dsi1Al = dsi1->GetAttribute(attributeIndicies[attrIdx]);
		AttributeLevel dsi2Al = dsi2->GetAttribute(attributeIndicies[attrIdx]);
		d += ((dsi1Al == dsi2Al) ? 0 : 1);
	}

	returnValue = (-3.0 / 4.0) * log(1.0 - (4.0 / 3.0) * ((double) d));

	return returnValue;
}

double Dataset::GetKimuraDistance(DatasetInstance* dsi1,
		DatasetInstance* dsi2) {
	vector<uint> attributeIndicies = MaskGetAttributeIndices(
			DISCRETE_TYPE);
	double p = 0, q = 0;
	for (uint attrIdx = 0; attrIdx < attributeIndicies.size();
			++attrIdx) {
		AttributeLevel dsi1Al = dsi1->GetAttribute(attributeIndicies[attrIdx]);
		AttributeLevel dsi2Al = dsi2->GetAttribute(attributeIndicies[attrIdx]);
		if (dsi1Al != dsi2Al) {
			if (GetAttributeMutationType(attributeIndicies[attrIdx])
					== TRANSITION_MUTATION) {
				++p;
			}
			if (GetAttributeMutationType(attributeIndicies[attrIdx])
					== TRANSVERSION_MUTATION) {
				++q;
			}
		}
	}
	double npos = (double) attributeIndicies.size();
	double P = p / npos;
	double Q = q / npos;
	double w1 = 1.0 - 2.0 * P - Q;
	double w2 = 1.0 - 2.0 * Q;
	double kimuraDistance = -0.5 * ProtectedLog(w1) - 0.25 * ProtectedLog(w2);

	return kimuraDistance;
}

uint Dataset::NumLevels(uint index) {
	if (index < attributeLevelsSeen.size()) {
		return attributeLevelsSeen[index].size();
	}
	cerr << "ERROR: Attempt to access number of levels of an attribute "
			<< "index that is out of range: " << index << " out of "
			<< attributeLevelsSeen.size() << std::endl;
	return 0;
}

uint Dataset::GetAttributeIndexFromName(string attributeName) {
	for (uint i = 0; i < attributeNames.size(); i++) {
		if (attributeNames[i] == attributeName) {
			return i;
		}
	}
	return INVALID_INDEX;
}

uint Dataset::NumNumerics() {
	return numericsMask.size();
}

vector<string> Dataset::GetNumericsNames() {
	vector<string> names;
	map<string, uint>::const_iterator it = numericsMask.begin();
	for (; it != numericsMask.end(); ++it) {
		names.push_back(it->first);
	}
	return names;
}

vector<string> Dataset::GetFileNumericsNames() {
	vector<string> names;
	vector<string>::const_iterator it = numericsNames.begin();
	for (; it != numericsNames.end(); ++it) {
		if (MaskSearchVariableType(*it, NUMERIC_TYPE)) {
			names.push_back(*it);
		}
	}
	return names;
}

pair<NumericLevel, NumericLevel> Dataset::GetMinMaxForNumeric(
		uint numericIdx) {
	return numericsMinMax[numericIdx];
}

double Dataset::GetMeanForNumeric(uint numericIdx) {
	double sum = 0.0;
	vector<uint> instanceIndicies = MaskGetInstanceIndices();
	for (uint i = 0; i < instanceIndicies.size(); ++i) {
		sum += instances[instanceIndicies[i]]->GetNumeric(numericIdx);
	}

	return (sum / ((double) instanceIndicies.size()));
}

bool Dataset::HasNumerics() {
	return hasNumerics;
}

void Dataset::HasNumerics(bool setHasNumerics) {
	hasNumerics = setHasNumerics;
}

NumericLevel Dataset::GetNumeric(unsigned instanceIndex, string name) {
	if (instanceIndex >= instances.size()) {
		cerr << "ERROR: Dataset::GetNumeric: instance index " << instanceIndex
				<< " out of range" << endl;
		exit(1);
	}
	map<string, uint>::iterator pos = numericsMask.find(name);
	if (pos != numericsMask.end()) {
		return instances[instanceIndex]->GetNumeric(pos->second);
	} else {
		cerr << "ERROR: Dataset::GetNumeric" << endl;
		exit(1);
	}
}

bool Dataset::GetNumericValues(string numericName,
		vector<NumericLevel>& numericValues) {
	if (MaskSearchVariableType(numericName, NUMERIC_TYPE)) {
		GetNumericValues(GetNumericIndexFromName(numericName), numericValues);
	} else {
		cerr
				<< "ERROR: Dataset::GetNumericValues cannot get numeric values for: "
				<< numericName << ". Either doesn't exist or is excluded"
				<< endl;
		return false;
	}
	return true;
}

std::string Dataset::GetNumericsFilename() {
	return numericsFilename;
}

uint Dataset::GetNumericIndexFromName(string numericName) {
  return numericsMask[numericName];
//	for (uint i = 0; i < numericsNames.size(); i++) {
//		if (numericsNames[i] == numericName) {
//			return i;
//		}
//	}
//	return INVALID_INDEX;
}

mat Dataset::GetNumericMatrix() {
  mat returnMatrix(NumInstances(), NumNumerics());
  returnMatrix.zeros();
  	vector<string> instanceIds = GetInstanceIds();
	for (uint iIdx = 0; iIdx < NumInstances(); iIdx++) {
		unsigned instanceIndex = 0;
		GetInstanceIndexForID(instanceIds[iIdx], instanceIndex);
		/// write continuous attribute values
		vector<uint> numIndices = MaskGetAttributeIndices(NUMERIC_TYPE);
		for (uint nIdx = 0; nIdx < numIndices.size(); nIdx++) {
			NumericLevel N = instances[instanceIndex]->GetNumeric(numIndices[nIdx]);
      returnMatrix(iIdx, nIdx) = N;
    }
  }
  return returnMatrix;
}

uint Dataset::NumClasses() {
	if (hasPhenotypes) {
		return classIndexes.size();
	}
	return INVALID_INT_VALUE;
}

uint Dataset::GetClassColumn() {
	if (hasPhenotypes) {
		return classColumn;
	} else {
		return INVALID_INT_VALUE;
	}
}

bool Dataset::GetClassValues(vector<ClassLevel>& classValues) {
	if (hasPhenotypes) {
		classValues.clear();
		map<string, uint>::const_iterator it;
		for (it = instancesMask.begin(); it != instancesMask.end(); it++) {
			classValues.push_back(instances[it->second]->GetClass());
		}
		return true;
	} else {
		return false;
	}
}

const std::map<ClassLevel, std::vector<uint> >&
Dataset::GetClassIndexes() {
	return classIndexes;
}

bool Dataset::HasAlternatePhenotypes() {
	return hasAlternatePhenotypes;
}

void Dataset::HasAlternatePhenotypes(bool setHasAlternatePhenotypes) {
	hasAlternatePhenotypes = setHasAlternatePhenotypes;
}

string Dataset::GetAlternatePhenotypesFilename() {
	return alternatePhenotypesFilename;
}

bool Dataset::HasContinuousPhenotypes() {
	return hasContinuousPhenotypes;
}

bool Dataset::HasPhenotypes() {
	return hasPhenotypes;
}

pair<double, double> Dataset::GetMinMaxForContinuousPhenotype() {
	return continuousPhenotypeMinMax;
}

void Dataset::Print() {
	PrintStats();
	cout << Timestamp() << "Data set values:" << endl << endl;
	map<string, uint>::const_iterator it;
	for (it = instancesMask.begin(); it != instancesMask.end(); it++) {
		DatasetInstance* dsi = instances[it->second];
		cout << instanceIds[it->second] << "\t";
		map<string, uint>::const_iterator ait = attributesMask.begin();
		for (; ait != attributesMask.end(); ++ait) {
			cout << dsi->GetAttribute(ait->second) << "\t";
		}
		map<string, uint>::const_iterator nit = numericsMask.begin();
		for (; nit != numericsMask.end(); ++nit) {
			cout << dsi->GetNumeric(nit->second) << "\t";
		}
		if (hasPhenotypes) {
			if (hasContinuousPhenotypes) {
				cout << dsi->GetPredictedValueTau() << endl;
			} else {
				cout << dsi->GetClass() << endl;
			}
		} else {
			cout << endl;
		}
	}
	cout << endl;
}

void Dataset::PrintStats() {
	uint numInstances = NumInstances();
	uint numClasses = NumClasses();
	uint numAttributes = NumAttributes();
	uint numNumerics = NumNumerics();
	uint numElements = (numInstances * (numAttributes + numNumerics))
			+ numInstances;

	// TODO: add class stats
	cout << Timestamp() << "Dataset has:" << endl << Timestamp()
			<< "instances:      " << numInstances << endl;
	if (hasGenotypes) {
		cout << Timestamp() << "SNPs:           " << numAttributes << endl;
	}
	if (hasNumerics) {
		cout << Timestamp() << "numerics:       " << numNumerics << endl;
	}
	if (hasPhenotypes) {
		if (hasContinuousPhenotypes) {
			cout << Timestamp() << "continuous phenotype, min: "
					<< continuousPhenotypeMinMax.first << ", max: "
					<< continuousPhenotypeMinMax.second << endl;
		} else {
			cout << Timestamp() << "classes:        " << numClasses << endl;
			PrintClassIndexInfo();
		}
	} else {
		cout << Timestamp() << "This data set has not phenotypic information"
				<< endl;
	}

	cout << Timestamp() << "total elements: " << numElements << endl;

	//PrintIntToStringMap();
	//PrintStringToIntMap();
	//  PrintLevelCounts();
	// PrintAttributeLevels();
	PrintMissingValuesStats();
}

void Dataset::PrintNumericsStats() {
	uint numInstances = NumInstances();
	uint numClasses = NumClasses();
	uint numAttributes = NumAttributes();
	uint numNumerics = NumNumerics();
	uint numElements = (numInstances * (numAttributes + numNumerics))
			+ numInstances;

	// TODO: add class stats
	cout << Timestamp() << "Dataset has:" << endl << Timestamp()
			<< "instances:      " << numInstances << endl;
	if (hasGenotypes) {
		cout << Timestamp() << "SNPs:           " << numAttributes << endl;
	}
	if (hasNumerics) {
		cout << Timestamp() << "numerics:       " << numNumerics << endl;
	}
//  vector< pair<double, double> >::const_iterator minMaxIt = numericsMinMax.begin();
//  for(uint i = 0; minMaxIt != numericsMinMax.end(); ++minMaxIt, ++i) {
//    cout << Timestamp()
//            << (*minMaxIt).first << " <= "
//            << numericsNames[i] << " <= "
//            << (*minMaxIt).second << endl;
//  }
	if (hasPhenotypes) {
		if (hasContinuousPhenotypes) {
			cout << Timestamp() << "continuous phenotype, min: "
					<< continuousPhenotypeMinMax.first << ", max: "
					<< continuousPhenotypeMinMax.second << endl;
		} else {
			cout << Timestamp() << "classes:        " << numClasses << endl;
			PrintClassIndexInfo();
		}
	} else {
		cout << Timestamp() << "This data set has not phenotypic information"
				<< endl;
	}
	cout << Timestamp() << "total elements: " << numElements << endl;

	PrintMissingValuesStats();
}

void Dataset::PrintStatsSimple(ostream& outStream) {
	uint numInstances = NumInstances();
	uint numClasses = NumClasses();
	uint numAttributes = NumAttributes();
	uint numNumerics = NumNumerics();
	uint numElements = (numInstances * (numAttributes + numNumerics))
			+ numInstances;

	outStream << Timestamp() << "Dataset has:" << endl << Timestamp()
			<< "instances:      " << numInstances << endl;
	if (hasGenotypes) {
		outStream << Timestamp() << "SNPs:           " << numAttributes << endl;
	}
	if (hasNumerics) {
		outStream << Timestamp() << "numerics:       " << numNumerics << endl;
	}

	if (hasPhenotypes) {
		if (hasContinuousPhenotypes) {
			outStream << Timestamp() << "continuous phenotype: " << endl;
			outStream << Timestamp() << "continuous phenotype, min: "
					<< continuousPhenotypeMinMax.first << ", max: "
					<< continuousPhenotypeMinMax.second << endl;
		} else {
			outStream << Timestamp() << "classes:        " << numClasses
					<< endl;
			PrintClassIndexInfo(outStream);
		}
	} else {
		outStream << Timestamp()
				<< "This data set has not phenotypic information" << endl;
	}

	outStream << Timestamp() << "total elements: " << numElements << endl;
}

void Dataset::PrintClassIndexInfo(ostream& outStream) {
	if (hasPhenotypes) {
		outStream << Timestamp() << "Data Set Class Index" << endl;
		outStream << Timestamp() << "Index has [" << classIndexes.size()
				<< "] entries:" << endl;
		map<ClassLevel, vector<uint> >::const_iterator mit =
				classIndexes.begin();
		for (; mit != classIndexes.end(); ++mit) {
			outStream << Timestamp() << (*mit).first << ": "
					<< (*mit).second.size() << endl;
		}
	} else {
		outStream << Timestamp() << "No index info since no phenotypes" << endl;
	}
}

void Dataset::PrintMissingValuesStats() {
	unsigned attrsMissing = 0;
	if (missingValues.size()) {
		cout << Timestamp() << "Missing Attributes Values Detected" << endl;
		cout << Timestamp() << "Instance # Missing" << endl;
		map<string, vector<uint> >::const_iterator mit;
		for (mit = missingValues.begin(); mit != missingValues.end(); ++mit) {
			cout << Timestamp() << mit->first << setw(10) << mit->second.size()
					<< endl;
			attrsMissing += mit->second.size();
		}
	} else {
		cout << Timestamp() << "0 missing attribute values detected" << endl;
	}
	double genotypingRate = 1.0
			- (attrsMissing
					/ ((double) instances.size() * attributeNames.size()));
	cout << Timestamp() << "Total genotyping rate: " << genotypingRate << endl;

	if (missingNumericValues.size()) {
		cout << Timestamp() << "Missing Numeric Values Detected" << endl;
		cout << Timestamp() << "Instance # Missing" << endl;
		map<string, vector<uint> >::const_iterator mit;
		for (mit = missingNumericValues.begin();
				mit != missingNumericValues.end(); ++mit) {
			cout << Timestamp() << mit->first << setw(10) << mit->second.size()
					<< endl;
		}
	} else {
		cout << Timestamp() << "0 missing numeric values detected" << endl;
	}
}

void Dataset::PrintLevelCounts() {
	cout << Timestamp() << "Data set attribute level counts:" << endl;
	vector<map<AttributeLevel, uint> >::const_iterator levelCountsIt =
			levelCounts.begin();
	for (uint attrIdx = 0; levelCountsIt != levelCounts.end();
			++levelCountsIt, ++attrIdx) {
		cout << Timestamp() << "Attribute [" << attrIdx << "]" << endl;
		map<AttributeLevel, uint>::const_iterator itsIt =
				(*levelCountsIt).begin();
		for (; itsIt != (*levelCountsIt).end(); ++itsIt) {
			cout << Timestamp() << (*itsIt).first << "->" << (*itsIt).second
					<< endl;
		}
	}
}

void Dataset::WriteLevelCounts(std::string levelsFilename) {
	if (levelCounts.size() == 0) {
		cout << "WARNING: No level counts to write to file" << endl;
	}

	ofstream outFile;
	outFile.open(levelsFilename.c_str());
	if (outFile.bad()) {
		cerr << "ERROR: Could not open levels counts file for writing" << endl;
		exit(EXIT_FAILURE);
	}

	cout << Timestamp() << "Writing level counts to [" << levelsFilename << "]"
			<< endl;
	vector<map<AttributeLevel, uint> >::const_iterator levelCountsIt =
			levelCounts.begin();
	for (uint attrIdx = 0; levelCountsIt != levelCounts.end();
			++levelCountsIt, ++attrIdx) {
		outFile << attributeNames[attrIdx] << "\t";

		pair<char, double> thisMinorAllele;
		pair<char, char> thisAttrAlleles;
		if (hasAllelicInfo) {
			thisMinorAllele = attributeMinorAllele[attrIdx];
			thisAttrAlleles = attributeAlleles[attrIdx];
			outFile << thisAttrAlleles.first << "\t" << thisAttrAlleles.second
					<< "\t" << thisMinorAllele.second;
		}

		map<AttributeLevel, uint>::const_iterator itsIt =
				(*levelCountsIt).begin();
		outFile << "\t[";
		bool first = true;
		for (; itsIt != (*levelCountsIt).end(); ++itsIt) {
			if (first) {
				outFile << itsIt->first << ":" << setw(4) << itsIt->second;
				first = false;
			} else {
				outFile << " / " << itsIt->first << ":" << setw(4)
						<< itsIt->second;
			}
		}
		outFile << " ]";

		if (hasAllelicInfo) {
			outFile << "\t" << thisMinorAllele.first << "\t"
					<< thisMinorAllele.second;
		} else {
			pair<char, double> minorAllele = GetAttributeMAF(attrIdx);
			outFile << "\t" << minorAllele.first << "\t" << minorAllele.second;
		}

		AttributeMutationType attributeMutationType = GetAttributeMutationType(
				attrIdx);
		string mutationTypeDesc = "Unknown";
		if (thisAttrAlleles.first == thisAttrAlleles.second) {
			mutationTypeDesc = "Monomorphic";
		} else {
			if (attributeMutationType == TRANSITION_MUTATION) {
				mutationTypeDesc = "Transition";
			}
			if (attributeMutationType == TRANSVERSION_MUTATION) {
				mutationTypeDesc = "Transversion";
			}
		}
		outFile << "\t" << mutationTypeDesc << endl;
	}

	outFile.close();
}

void Dataset::PrintAttributeLevelsSeen() {
	cout << Timestamp() << "Dataset attribute levels seen:" << endl;
	vector<set<string> >::const_iterator attributeIt =
			attributeLevelsSeen.begin();
	for (uint i = 0; attributeIt != attributeLevelsSeen.end();
			++attributeIt, ++i) {
		cout << Timestamp() << "Attribute [" << i << "] => [ ";
		set<string>::const_iterator sIt = (*attributeIt).begin();
		for (; sIt != (*attributeIt).end(); ++sIt) {
			cout << *sIt << " ";
		}
		cout << "]" << endl;
	}
}

bool Dataset::MaskRemoveVariable(string variableName) {
	if (MaskSearchVariableType(variableName, DISCRETE_TYPE)) {
		return MaskRemoveVariableType(variableName, DISCRETE_TYPE);
	}
	if (MaskSearchVariableType(variableName, NUMERIC_TYPE)) {
		return MaskRemoveVariableType(variableName, NUMERIC_TYPE);
	}
  error("Dataset::MaskRemoveVariable variableName: " + variableName + "\n");
}

bool Dataset::MaskRemoveVariableType(string variableName,
		AttributeType varType) {
	map<string, uint>::iterator pos;
	if (varType == DISCRETE_TYPE) {
		pos = attributesMask.find(variableName);
		if (pos != attributesMask.end()) {
			attributesMask.erase(pos);
		} else {
			error("Dataset::MaskRemoveVariable failed for SNP attribute name: "
					  + variableName + " name not found\n");
			return false;
		}
	} else {
		pos = numericsMask.find(variableName);
		if (pos != numericsMask.end()) {
			numericsMask.erase(pos);
		} else {
			error("Dataset::MaskRemoveVariable failed for numerics attribute name: " 
            + variableName + " name not found\n");
			return false;
		}
	}
	return true;
}

bool Dataset::MaskSearchVariableType(string variableName,
		AttributeType varType) {
	map<string, uint>::iterator pos;
	if (varType == DISCRETE_TYPE) {
		pos = attributesMask.find(variableName);
		if (pos != attributesMask.end()) {
			return true;
		} else {
			return false;
		}
	} else {
		pos = numericsMask.find(variableName);
		if (pos != numericsMask.end()) {
			return true;
		} else {
			return false;
		}
	}
}

bool Dataset::MaskIncludeAllAttributes(AttributeType attrType) {
	if (attrType == DISCRETE_TYPE) {
		attributesMask.clear();
		uint attributeIndex = 0;
		vector<string>::const_iterator it = attributeNames.begin();
		for (; it != attributeNames.end(); ++it) {
			attributesMask[*it] = attributeIndex;
			++attributeIndex;
		}
		return true;
	} else {
		numericsMask.clear();
		uint attributeIndex = 0;
		vector<string>::const_iterator it = numericsNames.begin();
		for (; it != numericsNames.end(); ++it) {
			numericsMask[*it] = attributeIndex;
			++attributeIndex;
		}
		return true;
	}
}

vector<uint> Dataset::MaskGetAttributeIndices(AttributeType attrType) {
	vector<uint> indices;
	if (attrType == DISCRETE_TYPE) {
		map<string, uint>::const_iterator it = attributesMask.begin();
		for (; it != attributesMask.end(); ++it) {
			indices.push_back(it->second);
		}
	} else {
		map<string, uint>::const_iterator it = numericsMask.begin();
		for (; it != numericsMask.end(); ++it) {
			indices.push_back(it->second);
		}
	}
	return indices;
}

const map<string, uint>&
Dataset::MaskGetAttributeMask(AttributeType attrType) {
	if (attrType == DISCRETE_TYPE) {
		return attributesMask;
	} else {
		return numericsMask;
	}
}

vector<string> Dataset::MaskGetAllVariableNames() {
	vector<string> names;
	map<string, uint>::const_iterator ait = attributesMask.begin();
	for (; ait != attributesMask.end(); ++ait) {
		names.push_back(ait->first);
	}
	map<string, uint>::const_iterator nit = numericsMask.begin();
	for (; nit != numericsMask.end(); ++nit) {
		names.push_back(nit->first);
	}
	return names;
}

bool Dataset::MaskRemoveInstance(std::string instanceId) {
	map<string, uint>::iterator pos = instancesMask.find(instanceId);
	if (pos != instancesMask.end()) {
		instancesMask.erase(pos);
	} else {
		cerr << "ERROR: Dataset::MaskRemoveInstance failed for instance ID: "
				<< instanceId << endl;
		return false;
	}
	return true;
}

bool Dataset::MaskSearchInstance(string instanceId) {
	map<string, uint>::iterator pos = instancesMask.find(instanceId);
	if (pos != instancesMask.end()) {
		return true;
	} else {
		return false;
	}
}

bool Dataset::MaskIncludeAllInstances() {
  PP->printLOG("Clearing the instance mask\n");
	instancesMask.clear();
  PP->printLOG("Adding all instance IDs to the instance mask\n");
	vector<string>::const_iterator it = instanceIds.begin();
	for(uint instanceIndex=0; it != instanceIds.end(); ++it, ++instanceIndex) {
    instancesMask[*it] = instanceIndex;
    cout << instanceIndex << "\t" << *it << endl;
	}
  cout << "Final mask size: " << instancesMask.size() << endl;
  
	return true;
}

vector<uint> Dataset::MaskGetInstanceIndices() {
	vector<uint> indices;
	map<string, uint>::const_iterator it = instancesMask.begin();
	for (; it != instancesMask.end(); ++it) {
		indices.push_back(it->second);
	}
	return indices;
}

vector<string> Dataset::MaskGetInstanceIds() {
	vector<string> ids;
	map<string, uint>::const_iterator it = instancesMask.begin();
	for (; it != instancesMask.end(); ++it) {
		ids.push_back(it->first);
	}
	return ids;
}

const map<string, uint>& Dataset::MaskGetInstanceMask() {
	return instancesMask;
}

bool Dataset::MaskPushAll() {
	if (!maskIsPushed) {
		attributesMaskPushed = attributesMask;
		numericsMaskPushed = numericsMask;
		instancesMaskPushed = instancesMask;
		maskIsPushed = true;
		return true;
	} else {
		cerr << "ERROR: only one level of pushing is allowed for "
				<< "attribute masks" << endl;
	}
	return false;
}

bool Dataset::MaskPopAll() {
	if (maskIsPushed) {
		attributesMask = attributesMaskPushed;
		numericsMask = numericsMaskPushed;
		instancesMask = instancesMaskPushed;
		maskIsPushed = false;
		return true;
	} else {
		cerr << "ERROR: attempt to pop an un-pushed attribute mask" << endl;
		return false;
	}
}

bool Dataset::MaskWriteNewDataset(string newDatasetFilename) {

	ofstream outFile;
	outFile.open(newDatasetFilename.c_str());
	if (outFile.bad()) {
		cerr << "ERROR: Could not open new data set file ["
				<< newDatasetFilename << "] for writing" << endl;
		return false;
	}

	if (HasNumerics()) {
		vector<string> numericNames = GetNumericsNames();
		outFile << join(numericNames.begin(), numericNames.end(), "\t");
		outFile << "Class" << endl;
		for (uint i = 0; i < NumInstances(); ++i) {
			ostringstream newInstanceString;
			for (uint j = 0; j < numericNames.size() - 1; ++j) {
				newInstanceString << GetNumeric(i, numericNames[j]) << "\t";
			}
			if (hasContinuousPhenotypes) {
				outFile << newInstanceString.str()
						<< instances[i]->GetPredictedValueTau() << endl;
			} else {
				outFile << newInstanceString.str() << instances[i]->GetClass()
						<< endl;
			}
		}
	} else {
		vector<string> attributeNames = GetAttributeNames();
		outFile << join(attributeNames.begin(), attributeNames.end(), "\t");
		outFile << "Class" << endl;
		for (uint i = 0; i < NumInstances(); ++i) {
			ostringstream newInstanceString;
			for (uint j = 0; j < attributeNames.size() - 1; ++j) {
				newInstanceString << GetAttribute(i, attributeNames[j]) << "\t";
			}
			if (hasContinuousPhenotypes) {
				outFile << newInstanceString.str()
						<< instances[i]->GetPredictedValueTau() << endl;
			} else {
				outFile << newInstanceString.str() << instances[i]->GetClass()
						<< endl;
			}
		}
	}
	outFile.close();

	return true;
}

void Dataset::PrintMaskStats() {
	PP->printLOG(Timestamp() + "Dataset Mask Statistics\n");
	PP->printLOG(Timestamp() + "Attributes mask size: " + int2str(attributesMask.size()) + "\n");
	PP->printLOG(Timestamp() + "Numerics mask size: " + int2str(numericsMask.size()) + "\n");
	PP->printLOG(Timestamp() + "Instances mask size: " + int2str(instancesMask.size()) + "\n");
	PP->printLOG(Timestamp() + "Mask is pushed?\n");
	PP->printLOG(Timestamp() + (maskIsPushed ? "true" : "false") + "\n");
}

void Dataset::RunSnpDiagnosticTests(string logFilename,
		double globalGenotypeThreshold, uint cellThreshold) {
	/// open the diagnostic log file
	ofstream outFile;
	outFile.open(logFilename.c_str());
	if (outFile.bad()) {
		cerr << "Could not open diagnostic log file [" << logFilename
				<< "] for writing" << endl;
		exit(EXIT_FAILURE);
	}
	cout << Timestamp() << "Writing diagnostics log to file [" << logFilename
			<< "]" << endl;

	if (!hasGenotypes) {
		cout << Timestamp() << "ERROR: This data set does not have any SNPs - "
				<< "Cannot perform diagnostic tests" << endl;
		return;
	}

	cout << Timestamp() << "BEGIN diagnostic tests BEGIN" << endl;
	cout << Timestamp() << "Diagnostic Information for: [" << snpsFilename
			<< "]" << endl;
	outFile << Timestamp() << "Diagnostic Information for: [" << snpsFilename
			<< "]" << endl;
	cout << Timestamp() << "Global genotype threshold: "
			<< globalGenotypeThreshold << endl;
	cout << Timestamp() << "X^2 cell threshold: " << cellThreshold << endl;
	outFile << Timestamp() << "Global genotype threshold: "
			<< globalGenotypeThreshold << endl;
	outFile << Timestamp() << "X^2 cell threshold: " << cellThreshold << endl;
	PrintStatsSimple();
	PrintStatsSimple(outFile);

	map<string, vector<string> > screwySnps;
	vector<string> badSnps;

	// ---------------------------------------------------------------------------
	// check for missing data in instances (subjects) - bcw - 7/12/11
	cout << Timestamp() << "missing values check" << endl;
	outFile << Timestamp() << "missing values check" << endl;
	uint totalMissing = 0;
	map<string, vector<uint> >::const_iterator mit =
			missingValues.begin();
	for (; mit != missingValues.end(); ++mit) {
		outFile << mit->first << endl;
		vector<uint>::const_iterator aiIt = (mit->second).begin();
		for (; aiIt != (mit->second).end(); ++aiIt) {
			outFile << *aiIt << " ";
		}
		outFile << endl;
		totalMissing += (*mit).second.size();
	}
	cout << Timestamp() << missingValues.size()
			<< " individuals had a total of " << totalMissing
			<< " missing values" << endl;
	outFile << Timestamp() << missingValues.size()
			<< " individuals had a total of " << totalMissing
			<< " missing values" << endl;

	// ---------------------------------------------------------------------------
	cout << Timestamp() << "global frequency count threshold check" << endl;
	outFile << Timestamp() << "global frequency count threshold check" << endl;
	uint attributeIndex = 0;
	uint instanceThreshold = (uint) (NumInstances()
			* globalGenotypeThreshold);
	vector<map<AttributeLevel, uint> >::const_iterator lcIt;
	// PrintLevelCounts();
	uint freqCountBad = 0;
	for (lcIt = levelCounts.begin(); lcIt != levelCounts.end();
			++lcIt, ++attributeIndex) {
		vector<uint> ftGenotypeCounts;
		map<AttributeLevel, uint>::const_iterator countsIt =
				(*lcIt).begin();
		for (uint levelIdx = 0; countsIt != (*lcIt).end();
				++countsIt, ++levelIdx) {
			string attributeName = attributeNames[attributeIndex];
			if ((*countsIt).second < instanceThreshold) {
				//        cout << "Attribute [" << attributeName
				//                << "] does not pass global threshold check of [" << instanceThreshold << "]. "
				//                << "Level: [" << (*countsIt).first
				//                << "], count : [" << (*countsIt).second
				//                << "]" << endl;
				ostringstream gtmsg;
				gtmsg << "does not pass global threshold check of ["
						<< instanceThreshold << "]. " << "Level: ["
						<< (*countsIt).first << "], count : ["
						<< (*countsIt).second << "]";
				screwySnps[attributeName].push_back(gtmsg.str());
				++freqCountBad;
			}
			if ((*countsIt).second == 0) {
				//        cout << "Attribute [" << attributeName
				//                << "] Level: [" << levelIdx
				//                << "] is not represented in the dataset" << endl << endl;
				ostringstream mgmsg;
				mgmsg << "Level: [" << levelIdx
						<< "] is not represented in the dataset";
				screwySnps[attributeName].push_back(mgmsg.str());
				++freqCountBad;
			}
			ftGenotypeCounts.push_back((*countsIt).second);
		}

		//    if(!CheckHardyWeinbergEquilibrium(genotypeCounts)) {
		//      cerr << "Attribute [" << attributeNames[attributeIndex]
		//              << "] fails the HWE test" << endl << endl;
		//      screwySnps[attributeNames[attributeIndex]].push_back("HWE");
		//    }
	}
	cout << Timestamp() << freqCountBad << " bad frequency counts" << endl;
	outFile << Timestamp() << freqCountBad << " bad frequency counts" << endl;

	// ---------------------------------------------------------------------------
	if (hasPhenotypes) {
		cout << Timestamp() << "chi-squared cell count threshold check" << endl;
		outFile << Timestamp() << "chi-squared cell count threshold check"
				<< endl;
		ChiSquared chiSq(this);
		uint chiSqBad = 0;
		for (uint attributeIndex = 0; attributeIndex < NumAttributes();
				attributeIndex++) {
			chiSq.ComputeScore(attributeIndex);
			vector<vector<double> > frequencyCounts =
					chiSq.GetFrequencyCounts();
			for (uint curClass = 0; curClass < NumClasses();
					curClass++) {
				for (uint curLevel = 0;
						curLevel < NumLevels(attributeIndex); curLevel++) {
					if (frequencyCounts[curClass][curLevel] < cellThreshold) {
						//          cout << "Attribute [" << attributeNames[attributeIndex]
						//                  << "] does not pass x^2 cell threshold check of [" << cellThreshold << "]. "
						//                  << "Level: [" << curLevel << "], Class: [" << curClass << "] => "
						//                  << "count: [" << frequencyCounts[curClass][curLevel]
						//                  << "]" << endl;
						ostringstream chimsg;
						chimsg << "does not pass x^2 cell threshold check of ["
								<< cellThreshold << "]. " << "Level: ["
								<< curLevel << "], Class: [" << curClass
								<< "] => " << "count: ["
								<< frequencyCounts[curClass][curLevel] << "]";
						screwySnps[attributeNames[attributeIndex]].push_back(
								chimsg.str());
						++chiSqBad;
					}
				}
			}
		}
		cout << Timestamp() << chiSqBad << " bad x^2 threshold checks" << endl;
		outFile << Timestamp() << chiSqBad << " bad x^2 threshold checks"
				<< endl;
	} else {
		cout << Timestamp()
				<< "WARNNING: x^2 test skipped - no phenotypes in this data set"
				<< endl;
		outFile << Timestamp()
				<< "WARNNING: x^2 test skipped - no phenotypes in this data set"
				<< endl;
	}

	// ---------------------------------------------------------------------------
	cout << Timestamp() << "Hardy-Weinberg Equilibrium (HWE) check" << endl;
	outFile << Timestamp() << "Hardy-Weinberg Equilibrium (HWE) check" << endl;
	uint badHWE = 0;
	for (uint attributeIndex = 0; attributeIndex < NumAttributes();
			attributeIndex++) {
		// get the vector of genotype counts
		map<AttributeLevel, uint> gcounts = levelCounts[attributeIndex];
		vector<uint> counts;
		map<AttributeLevel, uint>::const_iterator it = gcounts.begin();
		for (; it != gcounts.end(); ++it) {
			// cout << it->first << " " << it->second << endl;
			counts.push_back(it->second);
		}
		if (!CheckHardyWeinbergEquilibrium(counts)) {
			cout << Timestamp() << attributeNames[attributeIndex]
					<< " fails Hardy-Weinberg Equilibrium test" << endl;
			outFile << Timestamp() << attributeNames[attributeIndex]
					<< " fails Hardy-Weinberg Equilibrium test" << endl;
			++badHWE;
		}
	}
	cout << Timestamp() << badHWE << " SNPs failed HWE checks" << endl;
	outFile << Timestamp() << badHWE << " SNPs failed HWE checks" << endl;

	/// write diagnostic log information collected in screwySnps to file
	cout << Timestamp() << "Writing bad SNP information to log file" << endl;
	outFile << Timestamp() << "Writing bad SNP information to log file" << endl;
	if (screwySnps.size()) {
		map<string, vector<string> >::const_iterator it = screwySnps.begin();
		for (; it != screwySnps.end(); ++it) {
			outFile << (*it).first << "\n" << "\t"
					<< join((*it).second.begin(), (*it).second.end(), "\n\t")
					<< endl;
		}
	} else {
		cout << Timestamp() << "No failed diagnostics to write to log file"
				<< endl;
		outFile << Timestamp() << "No failed diagnostics to write to log file"
				<< endl;
	}

	cout << Timestamp() << "END diagnostic tests" << endl;
	outFile << Timestamp() << "END diagnostic tests" << endl;

	outFile.close();
}

bool Dataset::CheckHardyWeinbergEquilibrium(
		vector<uint>& chkGenotypeCounts) {
	/// observed counts
	uint AA = chkGenotypeCounts[0];
	uint Aa = chkGenotypeCounts[1];
	uint aa = chkGenotypeCounts[2];
	uint sum = AA + Aa + aa;
	if (sum < 1) {
		cout << Timestamp()
				<< "WARNING: Hardy-Weinberg test failed to find any observed values"
				<< endl;
		return false;
	}
	if (sum != NumInstances()) {
		cerr << Timestamp() << "HWE FAIL: genotype counts sum [" << sum
				<< "] does not add up to number of instances ["
				<< NumInstances() << "]" << endl;
		PrintLevelCounts();
		return false;
	}

	/// HWE probabilities
	double fAA = (double) AA / sum;
	double fAa = (double) Aa / sum;
	double p = fAA + 0.5 * fAa;
	double q = 1.0 - p;
	//  cout << p << " " << q << " " << endl;

	/// expected values
	double eAA = p * p * sum;
	double eAa = 2.0 * p * q * sum;
	double eaa = q * q * sum;
	if (eAA == 0 || eAa == 0 || eaa == 0) {
		cerr << Timestamp() << "HWE FAIL: HWE chi-square could not be "
				<< "computed due to zeroes in expected counts" << endl;
		return false;
	}
	uint expectedSum = (uint) (eAA + eAa + eaa + 0.5);
	if (expectedSum != NumInstances()) {
		cerr << Timestamp() << "HWE FAIL: HWE expected genotype counts sum ["
				<< expectedSum << "] does not add up to number of instances ["
				<< NumInstances() << "]" << endl;
		return false;
	}

	/// perform Pearson's chi-squared test
	double chiSquaredSum = ((AA - eAA) / eAA) + ((Aa - eAa) / eAa)
			+ ((aa - eaa) / eaa);

	//  cout << eAA << " " << eAa << " " << eaa << endl;

	/// one degree of freedom (# genotypes - # alleles), 5% significance level
	//double pValue = 1.0 - gsl_cdf_chisq_Q(chiSquaredSum, 1);
  chi_squared chiDist(5);
  double pValue = cdf(chiDist, chiSquaredSum);
	if (pValue > 0.1) {
		cerr << Timestamp() << "HWE FAIL: x^2 p-value too high: [" << pValue
				<< "]" << endl;
		cerr << Timestamp() << "Genotype counts: " << AA << " / " << Aa << " / "
				<< aa << endl;
		return false;
	}

	return true;
}

double Dataset::SNPHWE(int obs_hets, int obs_hom1, int obs_hom2) {
	if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) {
		printf("FATAL ERROR - SNP-HWE: Current genotype configuration "
				"(%d  %d %d ) includes a negative count", obs_hets, obs_hom1,
				obs_hom2);
		exit(EXIT_FAILURE);
	}

	int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

	int rare_copies = 2 * obs_homr + obs_hets;
	int genotypes = obs_hets + obs_homc + obs_homr;

	double * het_probs = (double *) malloc(
			(size_t) (rare_copies + 1) * sizeof(double));
	if (het_probs == NULL) {
		printf("FATAL ERROR - SNP-HWE: Unable to allocate array for "
				"heterozygote probabilities");
		exit(EXIT_FAILURE);
	}

	int i;
	for (i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;

	/* start at midpoint */
	int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1))
		mid++;

	int curr_hets = mid;
	int curr_homr = (rare_copies - mid) / 2;
	int curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
		het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets
				* (curr_hets - 1.0)
				/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		sum += het_probs[curr_hets - 2];

		/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
		curr_homr++;
		curr_homc++;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr
				* curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));
		sum += het_probs[curr_hets + 2];

		/* add 2 heterozygotes for next iteration -> subtract one rare,
		 * one common homozygote */
		curr_homr--;
		curr_homc--;
	}

	for (i = 0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	/* alternate p-value calculation for p_hi/p_lo
	 double p_hi = het_probs[obs_hets];
	 for (i = obs_hets + 1; i <= rare_copies; i++)
	 p_hi += het_probs[i];

	 double p_lo = het_probs[obs_hets];
	 for (i = obs_hets - 1; i >= 0; i--)
	 p_lo += het_probs[i];


	 double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
	 */

	double p_hwe = 0.0;
	/*  p-value calculation for p_hwe  */
	for (i = 0; i <= rare_copies; i++) {
		if (het_probs[i] > het_probs[obs_hets])
			continue;
		p_hwe += het_probs[i];
	}

	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

	free(het_probs);

	return p_hwe;
}

double Dataset::GetClassProbability(ClassLevel thisClass) {
	if (NumInstances()) {
		return ((double) classIndexes[thisClass].size()
				/ (double) NumInstances());
	}
	return 0.0;
}

double Dataset::GetProbabilityValueGivenClass(uint attributeIndex,
		AttributeLevel A, ClassLevel classValue) {
	map<pair<AttributeLevel, ClassLevel>, uint> thisAttrLevelCounts =
			levelCountsByClass[attributeIndex];
	uint instancesInThisClass = classIndexes[classValue].size();
	pair<AttributeLevel, ClassLevel> key = make_pair(A, classValue);
	if (thisAttrLevelCounts.find(key) != thisAttrLevelCounts.end()) {
		//    cout << attributeIndex << " [" << A << "|" << classValue << "]: "
		//            << thisAttrLevelCounts[key] << "/" << instancesInThisClass << endl;
		return ((double) thisAttrLevelCounts[key]
				/ (double) instancesInThisClass);
	} else {
		return 0.0;
	}
}

void Dataset::AttributeInteractionInformation() {
	/// call CalculateInteractionInformation
	map<pair<int, int>, map<string, double> > results;
	CalculateInteractionInformation(results);
	if (!results.size()) {
		return;
	}

	/// display results of interaction calculations
	cout << endl << "---------------------------------" << endl
			<< "Attribute Interaction Information" << endl
			<< "---------------------------------" << endl << endl;

	/// display results header line
	cout << "A\tB\t";
	map<pair<int, int>, map<string, double> >::const_iterator rIt =
			results.begin();
	map<string, double> hdrResults = (*rIt).second;
	map<string, double>::const_iterator rdIt;
	for (rdIt = hdrResults.begin(); rdIt != hdrResults.end(); rdIt++) {
		if (rdIt->first == "I(A;B|C)") {
			cout << fixed << setw(10) << "I(A;B|C)" << setw(18) << "( % )";
		} else {
			cout << setw(10) << rdIt->first;
		}
	}
	cout << endl;

	/// get the column sum
	double metricSum = 0.0;
	for (rIt = results.begin(); rIt != results.end(); rIt++) {
		map<string, double> thisResults = (*rIt).second;
		for (rdIt = thisResults.begin(); rdIt != thisResults.end(); rdIt++) {
			string metric = rdIt->first;
			if (metric == "I(A;B|C)") {
				metricSum += rdIt->second;
			}
		}
	}

	/// display results detail; I(A;B|C) column as percentage
	for (rIt = results.begin(); rIt != results.end(); rIt++) {
		pair<uint, uint> thisCombo = (*rIt).first;
		map<string, double> thisResults = (*rIt).second;
		cout << thisCombo.first << "\t" << thisCombo.second << "\t";
		for (rdIt = thisResults.begin(); rdIt != thisResults.end(); rdIt++) {
			string metric = rdIt->first;
			if (metric == "I(A;B|C)") {
				cout << fixed << setw(10) << rdIt->second
						<< setw(10) << "( "
						<< ((rdIt->second / metricSum) * 100) << "% )";
			} else {
				cout << fixed << setw(10) << rdIt->second;
			}
		}
		cout << endl;
	}
	cout << endl;
}

void Dataset::CalculateInteractionInformation(
		map<pair<int, int>, map<string, double> >& results) {

	/// so many way to fail before getting started
	if (!hasPhenotypes) {
		cerr << "ERROR: CalculateInteractionInformation requires discrete, "
				<< "case-control phenotypes: no phenotypes found" << endl;
		return;
	}
	if (hasContinuousPhenotypes) {
		cerr << "ERROR: CalculateInteractionInformation requires discrete, "
				<< "case-control phenotypes: continuous phenotypes found"
				<< endl;
		return;
	}
	if (hasNumerics) {
		cerr << "ERROR: CalculateInteractionInformation requires discrete, "
				<< "geontypic data: continuous data found" << endl;
		return;
	}

	/// vectors to hold sequences for attributes a and b with class c
	/// ab is an attribute constructed from the cartesian product of a and b
	vector<AttributeLevel> a;
	vector<AttributeLevel> b;
	vector<ClassLevel> c;
	vector<AttributeLevel> ab;

	/// Get the class values once
	GetClassValues(c);
	//copy(c.begin(), c.end(), ostream_iterator<InstanceClass>(cout, "\n"));

	/// for all possible (unique) interactions, ie nCk
	cout << Timestamp()
			<< "Constructing attribute interaction matrix in parallel:" << endl
			<< Timestamp();
	int numAttributes = NumAttributes();
	vector<string> attrNames = GetAttributeNames();
	/// use OpenMP to run in parallel the construction of the
	/// attribute interaction matrix
#pragma omp for schedule(dynamic, 1)
	for (int i = 0; i < numAttributes; i++) {
		/// THREAD STARTS
		for (int j = i + 1; j < numAttributes; j++) {

			/// load attribute values (columns) into vectors for entropy routines
			GetAttributeValues(attrNames[i], a);
			GetAttributeValues(attrNames[j], b);

			/// construct a new attribute that is the cartesian product of a and b
			ConstructAttributeCart(a, b, ab);

			/// compute all entropies, calculate information theoretic quantities
			/// and save the results
			map<string, double> r;

			r["H(A)"] = Entropy(a);
			r["H(B)"] = Entropy(b);
			r["H(C)"] = Entropy(c);
			r["H(AB)"] = Entropy(ab);

			r["H(A|C)"] = ConditionalEntropy(a, c);
			r["H(B|C)"] = ConditionalEntropy(b, c);
			r["H(AB|C)"] = ConditionalEntropy(ab, c);
			r["H(C|A)"] = ConditionalEntropy(c, a);
			r["H(C|B)"] = ConditionalEntropy(c, b);

			r["I(A;C)"] = r["H(C)"] - r["H(C|A)"];
			r["I(B;C)"] = r["H(C)"] - r["H(C|B)"];
			r["I(A;B)"] = r["H(A)"] + r["H(B)"] - r["H(AB)"];
			r["I(A;B|C)"] = r["H(A|C)"] + r["H(B|C)"] - r["H(AB|C)"];
			r["I(A;B;C)"] = r["I(A;B|C)"] - r["I(A;B)"];
			r["I(AB;C)"] = r["I(A;B;C)"] + r["I(A;C)"] + r["I(B;C)"];

			// from gain.py:
			//      def interaction_information(self, attrA, attrB):
			//        # I(A;B;C)=I(A;B|C)-I(A;B)
			//        # I(A;B;C)=H(AB)+H(BC)+H(AC)-H(A)-H(B)-H(C)-H(ABC)
			//        # where C is the class
			//
			//        H_ABC   = self.entropy(attrA, attrB, self.class_idx)
			//        H_AB    = self.entropy(attrA, attrB)
			//        H_AC    = self.entropy(attrA, self.class_idx)
			//        H_BC    = self.entropy(attrB, self.class_idx)
			//        H_A     = self.entropy(attrA)
			//        H_B     = self.entropy(attrB)
			//        H_C     = self.entropy(self.class_idx)
			//        return H_AB + H_BC + H_AC - H_A - H_B - H_C - H_ABC

			r["I_3_pgain"] = r["H(AB)"] + r["H(B|C)"] + r["H(A|C)"] - r["H(A)"]
					- r["H(B)"] - r["H(C)"] - r["H(AB|C)"];
			// I_3_paper from:
			// McKinney BA, Crowe JE Jr, Guo J, Tian D (2009) Capturing the Spectrum
			// of Interaction Effects in Genetic Association Studies by Simulated
			// Evaporative Cooling Network Analysis. PLoS Genet 5(3): e1000432.
			// doi:10.1371/journal.pgen.1000432
			r["I_3_paper"] = r["I(AB;C)"] - r["I(A;C)"] - r["I(B;C)"];
			// are the two I_3 values not close?
			//			if(fabs(I_3_pygain - I_3_paper) > 1e-6) {
			//				cerr << "WARNING: I_3 values are not the same. "
			//								<< "I_3 from pygain: " << I_3_pygain
			//								<< " I_3 from paper:  " << I_3_paper << endl;
			//			}

			//			cout << "Adding pair: " << i << ", " << j << endl;
			//      cout << attributeNames[i] << "\t" << attributeNames[j]
			//              << "\t" << r["I(A;B|C)"] << "\t" << r["I(A;C)"] << endl;

			//      map<string, double>::const_iterator it = r.begin();
			//      for(; it != r.end(); ++it) {
			//        cout << (*it).first << " = " << (*it).second << endl;
			//      }
			//      exit(1);

			results[make_pair(i, j)] = r;
		} /// inner loop over j ends
		  /// THREAD ENDS

		  // happy lights
		if (i && (i % 100 == 0)) {
			cout << i << "/" << numAttributes << " ";
			cout.flush();
		}
		if (i && ((i % 1000) == 0)) {
			cout << endl << Timestamp();
		}
	} /// outer loop over i ends
	cout << numAttributes << "/" << numAttributes << " done" << endl;

	/// end - no return value
}

bool Dataset::CalculateGainMatrix(double** gainMatrix, string matrixFilename) {
	if (!HasGenotypes() || HasNumerics()) {
		cerr << "Dataset::CalculateGainMatrix only works on SNP data" << endl;
		return false;
	}

	if ((!HasPhenotypes()) || (NumClasses() != 2)) {
		cerr
				<< "Dataset::CalculateGainMatrix only works with case-control phenotypes"
				<< endl;
		return false;
	}

	map<string, uint> attributeMask = MaskGetAttributeMask(
			DISCRETE_TYPE);
	vector<string> attributeNames = MaskGetAllVariableNames();
	int numAttributes = attributeNames.size();

	/// Calculate the interaction information from entropies
	map<pair<int, int>, map<string, double> > results;
	CalculateInteractionInformation(results);
	/// Populate the GAIN matrix
	vector<AttributeLevel> c;
	vector<AttributeLevel> a;
	GetClassValues(c);
	for (int i = 0; i < numAttributes; i++) {
		GetAttributeValues(attributeNames[i], a);
		gainMatrix[i][i] = SelfEntropy(a, c);
		for (int j = i + 1; j < numAttributes; j++) {
			map<string, double> r = results[make_pair(i, j)];
			gainMatrix[i][j] = gainMatrix[j][i] = r["I_3_paper"];
		}
	}

	if (matrixFilename != "") {
		cout << Timestamp() << "Writing GAIN matrix to file [" << matrixFilename
				<< "]" << endl;
		ofstream outFile(matrixFilename.c_str());
		/// write header
		for (int i = 0; i < numAttributes; ++i) {
			if (i) {
				outFile << "\t" << attributeNames[i];
			} else {
				outFile << attributeNames[i];
			}
		}
		outFile << endl;
		/// write all m-by-m matrix entries
		for (int i = 0; i < numAttributes; ++i) {
			for (int j = i; j < numAttributes; ++j) {
				if (j == i) { // fill in symmetric entries, replacing j < i with tabs					string tabs = "";
					string tabs = "";
					for (int k = 0; k < j; k++)
						tabs += "\t";
					outFile << tabs << gainMatrix[i][j];
				} else
					outFile << "\t" << gainMatrix[i][j];
			}
			outFile << endl;
		}
		outFile.close();
	}

	return true;
}

bool Dataset::CalculateRegainMatrix(double** gainMatrix,
		string matrixFilename) {
	if (HasGenotypes() || !HasNumerics()) {
		cerr << "Dataset::CalculateRegainMatrix only works on numeric data" << endl;
		return false;
	}

	if ((!HasPhenotypes()) || (NumClasses() != 2)) {
		cerr
				<< "Dataset::CalculateGainMatrix only works with case-control phenotypes"
				<< endl;
		return false;
	}

	map<string, uint> attributeMask = MaskGetAttributeMask(
			NUMERIC_TYPE);
	vector<string> attributeNames = MaskGetAllVariableNames();
	int numAttributes = attributeNames.size();

	// do linear or logistic regression here
	// THIS IS IMPLEMENTED IN THE inbix PROJECT.

	if (matrixFilename != "") {
		cout << Timestamp() << "Writing GAIN matrix to file [" << matrixFilename
				<< "]" << endl;
		ofstream outFile(matrixFilename.c_str());
		/// write header
		for (int i = 0; i < numAttributes; ++i) {
			if (i) {
				outFile << "\t" << attributeNames[i];
			} else {
				outFile << attributeNames[i];
			}
		}
		outFile << endl;
		/// write all m-by-m matrix entries
		for (int i = 0; i < numAttributes; ++i) {
			for (int j = i; j < numAttributes; ++j) {
				if (j == i) { // fill in symmetric entries, replacing j < i with tabs					string tabs = "";
					string tabs = "";
					for (int k = 0; k < j; k++)
						tabs += "\t";
					outFile << tabs << gainMatrix[i][j];
				} else
					outFile << "\t" << gainMatrix[i][j];
			}
			outFile << endl;
		}
		outFile.close();
	}

	return true;
}

vector<double> Dataset::GetMAFs() {
  vector<double> retVals(NumAttributes());
  for(int i=0; i < attributeMinorAllele.size(); ++i) {
    retVals[i] = attributeMinorAllele[i].second;
  }
  return retVals;
}

double Dataset::ComputeInstanceToInstanceDistance(DatasetInstance* dsi1,
		DatasetInstance* dsi2) {
	double distance = 0;

	if (HasGenotypes()) {
		if (snpNearestNeighborMetricName == "KM") {
			distance = GetKimuraDistance(dsi1, dsi2);
		} else {
			if (snpNearestNeighborMetricName == "JC") {
				distance = GetJukesCantorDistance(dsi1, dsi2);
			} else {
				vector<uint> attributeIndices = MaskGetAttributeIndices(DISCRETE_TYPE);
				for (uint i = 0; i < attributeIndices.size(); ++i) {
					// DEBUG
			  	pair<char, char> alleles =
		  			dsi1->GetDatasetPtr()->GetAttributeAlleles(attributeIndices[i]);
			  	string a1 = " ";
			  	a1[0] = alleles.first;
			  	string a2 = " ";
			  	a2[0] = alleles.second;
			  	map<AttributeLevel, string> genotypeMap;
			  	genotypeMap[0] = a1 + a1;
			  	genotypeMap[1] = a1 + a2;
			  	genotypeMap[2] = a2 + a2;
			  	AttributeLevel attrLevel1 = dsi1->GetAttribute(attributeIndices[i]);
			  	AttributeLevel attrLevel2 = dsi2->GetAttribute(attributeIndices[i]);
			  	string genotype1 = genotypeMap[attrLevel1];
			  	string genotype2 = genotypeMap[attrLevel2];
			  	double tempDistance = snpNearestNeighborFuncPtr(attributeIndices[i], dsi1, dsi2);
					distance += tempDistance;
					// cout 
					// 	<< a1 << "\t" << a2 << "\t" 
					// 	<< attrLevel1 << "\t" << attrLevel2 << "\t"
					// 	<< genotype1 << "\t" << genotype2 << "\t" 
					// 	<< tempDistance << "\t" << distance << endl;
				}
			}
		}
		// cout << "SNP distance = " << distance << endl;
	}

	// added 6/16/11
	// compute numeric distances
	if (HasNumerics()) {
		// cout << "Computing numeric instance-to-instance distance..." << endl;
		vector<uint> numericIndices = MaskGetAttributeIndices(NUMERIC_TYPE);
		// cout << "\tNumber of numerics: " << numericIndices.size() << endl;
		for (uint i = 0; i < numericIndices.size(); ++i) {
			// cout << "\t\tNumeric index: " << numericIndices[i] << endl;
			double numDistance = numDiffFuncPtr(numericIndices[i], dsi1, dsi2);
			// cout << "Numeric distance " << i << " => " << numDistance << endl;
			distance += numDistance;
		}
	}

	return distance;
}

bool Dataset::SetDistanceMetrics(string newSnpDiffMetric, string newSnpNNMetric, 
	string newNumMetric) {
	/// set the SNP metric function pointer
	bool snpDiffMetricFunctionUnset = true;
	if (snpDiffMetricFunctionUnset && to_upper(newSnpDiffMetric) == "GM") {
		snpDiffFuncPtr = diffGMM;
		snpDiffMetricFunctionUnset = false;
	}
	if (snpDiffMetricFunctionUnset && to_upper(newSnpDiffMetric) == "AM") {
		snpDiffFuncPtr = diffAMM;
		snpDiffMetricFunctionUnset = false;
	}
	if (snpDiffMetricFunctionUnset && to_upper(newSnpDiffMetric) == "NCA") {
		snpDiffFuncPtr = diffNCA;
		snpDiffMetricFunctionUnset = false;
	}
	if (snpDiffMetricFunctionUnset && to_upper(newSnpDiffMetric) == "NCA6") {
		snpDiffFuncPtr = diffNCA6;
		snpDiffMetricFunctionUnset = false;
	}
	if (snpDiffMetricFunctionUnset && to_upper(newSnpDiffMetric) == "KM") {
		snpDiffFuncPtr = diffKM;
		snpDiffMetricFunctionUnset = false;
	}
	if (snpDiffMetricFunctionUnset) {
		cerr << "ERROR: Cannot set SNP diff metric to [ "
				<< newSnpDiffMetric << " ]" << endl;
		return false;
	}
	snpDiffMetricName = newSnpDiffMetric;

	/// set the nearest neighbors metric
	bool snpNNMetricFunctionUnset = true;
	if (snpNNMetricFunctionUnset && to_upper(newSnpNNMetric) == "GM") {
		snpNearestNeighborFuncPtr = diffGMM;
		snpNNMetricFunctionUnset = false;
	}
	if (snpNNMetricFunctionUnset && to_upper(newSnpNNMetric) == "AM") {
		snpNearestNeighborFuncPtr = diffAMM;
		snpNNMetricFunctionUnset = false;
	}
	if (snpNNMetricFunctionUnset && to_upper(newSnpNNMetric) == "NCA") {
		snpNearestNeighborFuncPtr = diffNCA;
		snpNNMetricFunctionUnset = false;
	}
	if (snpNNMetricFunctionUnset && to_upper(newSnpNNMetric) == "NCA6") {
		snpNearestNeighborFuncPtr = diffNCA6;
		snpNNMetricFunctionUnset = false;
	}
	if (snpNNMetricFunctionUnset && to_upper(newSnpNNMetric) == "KM") {
		snpNearestNeighborFuncPtr = diffKM;
		snpNNMetricFunctionUnset = false;
	}
	if (snpNNMetricFunctionUnset && to_upper(newSnpNNMetric) == "GRM") {
    // no need to set a function pointer for GRM, matrix computed all at once
		snpNearestNeighborFuncPtr = 0;
		snpNNMetricFunctionUnset = false;
	}
	if (snpNNMetricFunctionUnset) {
		cerr << "ERROR: Cannot set SNP nearest neighbors metric to [ "
				<< newSnpNNMetric << " ]" << endl;
		return false;
	}
	snpNearestNeighborMetricName = newSnpNNMetric;

	if (to_upper(newNumMetric) == "MANHATTAN") {
		numDiffFuncPtr = diffManhattan;
	} else {
		if (to_upper(newNumMetric) == "EUCLIDEAN") {
			numDiffFuncPtr = diffEuclidean;
		}
		else {
			cerr << "ERROR: [ " << newNumMetric
					<< " ] is not a valid numeric metric type" << endl;
			return false;
		}
	}
	numDiffMetricName = newNumMetric;

	cout << Timestamp() << "New SNP distance diff metric: "
			<< snpDiffMetricName << " " << (void*) snpDiffFuncPtr << endl;
	cout << Timestamp() << "New SNP distance nearest neighbor metric: "
			<< snpNearestNeighborMetricName << " " 
          << (void *) snpNearestNeighborFuncPtr << endl;
	cout << Timestamp() << "New continuous distance metric: " 
      << numDiffMetricName << " " << (void *) numDiffFuncPtr << endl;

	return true;
}

vector<string> Dataset::GetDistanceMetrics() {
	vector<string> metrics;
	metrics.push_back(snpDiffMetricName);
	metrics.push_back(snpNearestNeighborMetricName);
	metrics.push_back(numDiffMetricName);

	return metrics;
}

pair<uint, uint> Dataset::GetAttributeTiTvCounts() {
	vector<uint> attrIndices = MaskGetAttributeIndices(DISCRETE_TYPE);
	uint tiCount = 0, tvCount = 0;
	for (uint aIdx = 0; aIdx < attrIndices.size(); aIdx++) {
		if (attributeMutationTypes[aIdx] == TRANSITION_MUTATION) {
			++tiCount;
		} else {
			++tvCount;
		}
	}

	return make_pair(tiCount, tvCount);
}

bool Dataset::WriteSnpTiTvInfo(string titvFilename) {
	if (!hasAllelicInfo) {
		cerr
				<< "Dataset::WriteSnpTiTvInfo Cannot write transition/transversion "
						"for data with no allelic info" << endl;
    
		return false;
	}

	ofstream outFile(titvFilename.c_str());
	vector<AttributeMutationType>::const_iterator it =
			attributeMutationTypes.begin();
	vector<std::string>::const_iterator nit = attributeNames.begin();
	for (; it != attributeMutationTypes.end(); ++it, ++nit) {
		outFile << *nit << "\t" << *it << endl;
	}
	outFile.close();

	return true;
}

bool Dataset::ResetNearestNeighbors() {
//	cout << Timestamp() 
//					<< "INFO: Dataset is clearing instance nearest neighbor
//         	<< "information" << endl;
	map<string, uint>::const_iterator it = instancesMask.begin();
	for (; it != instancesMask.end(); ++it) {
		instances[it->second]->ResetNearestNeighbors();
	}	
	
	return true;
}


bool Dataset::CalculateDistanceMatrix(double** distanceMatrix,
		string matrixFilename) {
	cout << Timestamp() << "Calculating distance matrix" << endl;
	map<string, uint> instanceMask = MaskGetInstanceMask();
	vector<string> instanceIds = MaskGetInstanceIds();
	int numInstances = instanceIds.size();

	// populate the matrix - upper triangular
	// NOTE: make complete symmetric matrix for neighbor-to-neighbor sums
	cout << Timestamp() 
          <<  "Computing instance-to-instance nearest neighbor distances with " 
          << snpNearestNeighborMetricName << "... " << endl;
	//  omp_set_nested(1);
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < numInstances; ++i) {
		// cout << "Computing instance to instance distances. Row: " << i << endl;
		// #pragma omp parallel for
		for (int j = i + 1; j < numInstances; ++j) {
			uint dsi1Index = 0;
			GetInstanceIndexForID(instanceIds[i], dsi1Index);
			uint dsi2Index = 0;
			GetInstanceIndexForID(instanceIds[j], dsi2Index);
#pragma omp critical
      {
			distanceMatrix[i][j] = 
      distanceMatrix[j][i] =
					ComputeInstanceToInstanceDistance(
              GetInstance(dsi1Index),
							GetInstance(dsi2Index));
			// cout << i << ", " << j << " => " << distanceMatrix[i][j] << endl;
      }
    }
    #pragma omp critical 
{
    if (i && (i % 100 == 0)) {
      cout << Timestamp() << i << "/" << numInstances << endl;
    }
    distanceMatrix[i][i] = 0.0;
}
  }
  // end parallel openMP section


	cout << Timestamp() << numInstances << "/" << numInstances << " done"
			<< endl;

	if (matrixFilename != "") {
		string outFilename = matrixFilename + ".mat";
		cout << Timestamp() << "Writing distance matrix to file ["
				<< outFilename << "]" << endl;
		ofstream outFile(outFilename.c_str());
		/// write header
		for (int i = 0; i < numInstances; ++i) {
			if (i) {
				outFile << "\t" << instanceIds[i];
			} else {
				outFile << instanceIds[i];
			}
		}
		outFile << endl;
		/// write all n-by-n matrix entries
		for (int i = 0; i < numInstances; ++i) {
			for (int j = 0; j < numInstances; ++j) {
				if (j)
					outFile << "\t" << distanceMatrix[i][j];
				else
					outFile << distanceMatrix[i][j];
			}
			outFile << endl;
		}
		outFile.close();

		/// write phenotypes for the instances
		string phenoFilename = matrixFilename + ".matpheno";
		cout << Timestamp() << "Writing phenotypes to file [" << phenoFilename
				<< "]" << endl;
		ofstream phenoFile(phenoFilename.c_str());
		for (int i = 0; i < numInstances; ++i) {
			uint dsiIndex = 0;
			GetInstanceIndexForID(instanceIds[i], dsiIndex);
			if (hasContinuousPhenotypes) {
				phenoFile << instances[dsiIndex]->GetPredictedValueTau()
						<< endl;
			} else {
				phenoFile << instances[dsiIndex]->GetClass() << endl;
			}
		}
		phenoFile.close();
	}

	return true;
}

bool Dataset::CalculateDistanceMatrix(vector<vector<double> >& distanceMatrix) {
	map<string, uint> instanceMask = MaskGetInstanceMask();
	vector<string> instanceIds = MaskGetInstanceIds();
	int numInstances = instanceIds.size();

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < numInstances; ++i) {
		for (int j = i + 1; j < numInstances; ++j) {
			uint dsi1Index = 0;
			GetInstanceIndexForID(instanceIds[i], dsi1Index);
			uint dsi2Index = 0;
			GetInstanceIndexForID(instanceIds[j], dsi2Index);
			distanceMatrix[i][j] = distanceMatrix[j][i] =
					ComputeInstanceToInstanceDistance(GetInstance(dsi1Index),
							GetInstance(dsi2Index));
		}
		distanceMatrix[i][i] = 0.0;
	}

	return true;
}

/// ------------ Beginning of private methods ------------------

bool Dataset::LoadSnps(std::string filename) {

	/// Open the data file and read line-by-line
	snpsFilename = filename;
	ifstream dataStream(snpsFilename.c_str());
	if (!dataStream.is_open()) {
		cerr << "ERROR: Could not open SNP data set: " << snpsFilename << endl;
		return false;
	}
	cout << Timestamp()
			<< "Reading whitespace-delimited SNP data set lines from "
			<< snpsFilename << ":" << endl;

	// temporary string for reading file lines
	string line;

	// read the header row - whitespace delimited attribute names
	// special attribute named "class" can be in any position
	// keep attribute names
	getline(dataStream, line);
	vector<string> tokens;
	split(tokens, line);
	vector<string>::const_iterator it;
	uint numAttributes = 0;
	uint classIndex = 0;
	for (it = tokens.begin(); it != tokens.end(); ++it) {
		if (to_upper(*it) == "CLASS") {
			cout << Timestamp() << "Class column detect at " << classIndex
					<< endl;
			classColumn = classIndex;
		} else {
			attributeNames.push_back(*it);
			attributesMask[*it] = numAttributes;
			++numAttributes;
		}
		++classIndex;
	}

	/// Detect the class type
	bool classDetected = false;
	switch (DetectClassType(filename, classColumn + 1, true)) {
	case CASE_CONTROL_CLASS_TYPE:
		cout << Timestamp() << "Case-control phenotypes detected" << endl;
		hasContinuousPhenotypes = false;
		classDetected = true;
		break;
	case CONTINUOUS_CLASS_TYPE:
		cout << Timestamp() << "Continuous phenotypes detected" << endl;
		hasContinuousPhenotypes = true;
		classDetected = true;
		break;
	case MULTI_CLASS_TYPE:
		cout << "ERROR: more than two discrete phenotypes detected" << endl;
		break;
	case NO_CLASS_TYPE:
		cout << "ERROR: phenotypes could not be detected" << endl;
		break;
	}
	if (!classDetected) {
		return false;
	}

	levelCounts.resize(numAttributes);
	levelCountsByClass.resize(numAttributes);
	attributeLevelsSeen.resize(numAttributes);
	genotypeCounts.resize(numAttributes);

	// read instance attributes from whitespace-delimited lines
	uint instanceIndex = 0;
	double minPheno = 0.0, maxPheno = 0.0;
	uint lineNumber = 0;
	while (getline(dataStream, line)) {
		++lineNumber;
		// only load matching IDs from numerics and phenotype files
		// the delimited text file uses line numbers for instance IDs
		string trimmedLine = trim(line);
		// skip blank lines in the data section (usually end of file)
		if (!trimmedLine.size()) {
			continue;
		}

		string ID = zeroPadNumber(lineNumber, 8) + zeroPadNumber(lineNumber, 8);
		if (!IsLoadableInstanceID(ID)) {
			cout << Timestamp() << "WARNING: Dataset ID [" << ID
					<< "] skipped. " << " Line number: " << lineNumber
					<< ". Not found in numerics and/or phenotype file(s)"
					<< endl;
			continue;
		}

		vector<string> attributesStringVector;
		split(attributesStringVector, trimmedLine);
		uint numAttributesRead = attributesStringVector.size() - 1;
		if ((numAttributesRead == 0) || (numAttributesRead != numAttributes)) {
			cout << Timestamp() << "WARNING: Skipping line " << lineNumber
					<< " instance has " << numAttributesRead << " should have "
					<< numAttributes << endl;
			continue;
		}

		vector<AttributeLevel> attributeVector;
		uint attrIdx = 0;
		vector<string>::const_iterator it = attributesStringVector.begin();
		ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
		NumericLevel numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
		for (; it != attributesStringVector.end(); ++it, ++attrIdx) {
			string thisAttr = *it;
			if (attrIdx == classColumn) {
				discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
				numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
				if (hasContinuousPhenotypes) {
					if (thisAttr != "-9") {
						numericClassLevel = lexical_cast<NumericLevel>(
								thisAttr);
						if (lineNumber == 1) {
							minPheno = maxPheno = numericClassLevel;
						} else {
							if (numericClassLevel < minPheno) {
								minPheno = numericClassLevel;
							}
							if (numericClassLevel > maxPheno) {
								maxPheno = numericClassLevel;
							}
						}
					} else {
						if (!hasAlternatePhenotypes) {
							cout << Timestamp() << "Instance ID " << ID
									<< " filtered out by missing value" << endl;
							continue;
						}
					}
				} else {
					if (thisAttr != "-9") {
						discreteClassLevel = lexical_cast<ClassLevel>(thisAttr);
					} else {
						if (!hasAlternatePhenotypes) {
							cout << Timestamp() << "Instance ID " << ID
									<< " filtered out by missing value" << endl;
							continue;
						}
					}
				}
			} else {
				// attributes
				AttributeLevel thisAttrLevel = MISSING_ATTRIBUTE_VALUE;
				if (thisAttr == "?") {
					missingValues[ID].push_back(attrIdx);
				} else {
					thisAttrLevel = lexical_cast<AttributeLevel>(thisAttr);
					attributeLevelsSeen[attrIdx].insert(thisAttr);
				}
				attributeVector.push_back(thisAttrLevel);
			}
		}

		// create an instance from the vector of attribute and class values
		if (attributeVector.size() != numAttributes) {
			cerr << "ERROR: Number of attributes parsed on line " << lineNumber
					<< ": " << attributesStringVector.size()
					<< " is not equal to the number of attributes"
					<< " read from the data file header: " << numAttributes
					<< endl;
			return false;
		}

		DatasetInstance* newInst = new DatasetInstance(this, ID);
		if (newInst) {
			if (hasContinuousPhenotypes) {
				newInst->SetPredictedValueTau(numericClassLevel);
			} else {
				newInst->SetClass(discreteClassLevel);
				classIndexes[discreteClassLevel].push_back(instanceIndex);
			}
			newInst->LoadInstanceFromVector(attributeVector);
			instances.push_back(newInst);
			instanceIds.push_back(ID);
			instancesMask[ID] = instanceIndex;
		} else {
			cerr << "ERROR: loading tab-delimited data set. "
					<< "Could not create dataset instance for line number "
					<< lineNumber << endl;
			return false;
		}
		++instanceIndex;

		// happy lights
		if ((lineNumber - 1) && ((lineNumber % 100) == 0)) {
			cout << Timestamp() << lineNumber << endl;
		}
	}
	cout << Timestamp() << lineNumber << " lines read" << endl;

	dataStream.close();

	cout << Timestamp() << "There are " << NumInstances()
			<< " instances in the data set" << endl;
	cout << Timestamp() << "There are " << instancesMask.size()
			<< " instances in the instance mask" << endl;
	if (instancesMask.size() == 0) {
		cerr << "ERROR: no instances in the instance mask" << endl;
		return false;
	}

	if (hasContinuousPhenotypes) {
		continuousPhenotypeMinMax = make_pair(minPheno, maxPheno);
		cout << Timestamp() << "Continuous phenotypes." << endl;
	} else {
		cout << Timestamp() << "There are " << classIndexes.size()
				<< " classes in the data set" << endl;
	}

	hasGenotypes = true;
	UpdateAllLevelCounts();

	CreateDummyAlleles();

	return true;
}

void Dataset::UpdateAllLevelCounts() {
	cout << Timestamp() << "Updating all level counts:" << endl;
  if(NumAttributes() < 1) {
    error("Dataset::UpdateAllLevelCounts: No attributes to update\n");
  }
  if(!hasGenotypes) {
    error("Dataset::UpdateAllLevelCounts: Level counts only works on genotype data\n");
  }
	levelCounts.clear();
	levelCounts.resize(NumAttributes());
	if (hasPhenotypes && !hasContinuousPhenotypes) {
		levelCountsByClass.clear();
		levelCountsByClass.resize(NumAttributes());
	}
	/// initialize level count maps to contain at least three levels
	for (uint i = 0; i < NumAttributes(); ++i) {
		levelCounts[i][0] = 0;
		levelCounts[i][1] = 0;
		levelCounts[i][2] = 0;
	}
	uint instanceCount = 0;
	map<string, uint>::const_iterator it = instancesMask.begin();
	for (; it != instancesMask.end(); ++it) {
		UpdateLevelCounts(instances[it->second]);
		if (instanceCount && ((instanceCount % 100) == 0)) {
			cout << Timestamp() << instanceCount << "/" << instancesMask.size()
					<< endl;
			cout.flush();
		}
		++instanceCount;
	}
	cout << Timestamp() << instanceCount << "/" << instancesMask.size()
			<< " done" << endl;
}

void Dataset::UpdateLevelCounts(DatasetInstance* dsi) {
	ClassLevel thisClassLevel = dsi->GetClass();
	map<string, uint>::const_iterator it = attributesMask.begin();
	for (; it != attributesMask.end(); ++it) {
		uint attributeIndex = it->second;
		AttributeLevel thisAttributeLevel = dsi->GetAttribute(attributeIndex);
		if (thisAttributeLevel != MISSING_ATTRIBUTE_VALUE) {
			++levelCounts[attributeIndex][thisAttributeLevel];
			if (hasPhenotypes && !hasContinuousPhenotypes) {
				++levelCountsByClass[attributeIndex][make_pair(
						thisAttributeLevel, thisClassLevel)];
			}
		}
	}
}

void Dataset::ExcludeMonomorphs() {
	uint attrIdx = 0;
	uint attrsExcluded = 0;
	vector<map<AttributeLevel, uint> >::const_iterator it;
	for (it = levelCounts.begin(); it != levelCounts.end(); ++it, ++attrIdx) {
		if (it->size() == 1) {
			cout << Timestamp() << "WARNING: attribute "
					<< attributeNames[attrIdx]
					<< " is monomorphic and being marked as excluded from the analysis"
					<< endl;
			MaskRemoveVariable(attributeNames[attrIdx]);
			++attrsExcluded;
		}
	}
	cout << Timestamp() << attrsExcluded << " SNPs excluded as monomorphic"
			<< endl;
}

void Dataset::CreateDummyAlleles() {
	attributeAlleles.clear();
	attributeAlleleCounts.clear();
	attributeMinorAllele.clear();
	map<string, uint>::const_iterator ait = attributesMask.begin();
	for (; ait != attributesMask.end(); ++ait) {
		map<string, uint>::const_iterator it = instancesMask.begin();
		map<char, uint> alleleCounts;
		for (; it != instancesMask.end(); ++it) {
			DatasetInstance* dsi = instances[it->second];
			uint attributeIndex = ait->second;
			AttributeLevel thisAttributeLevel = dsi->GetAttribute(
					attributeIndex);
			if (thisAttributeLevel == MISSING_ATTRIBUTE_VALUE) {
				continue;
			} else {
				switch (thisAttributeLevel) {
				case 0:
					alleleCounts['X'] += 2;
					break;
				case 1:
					++alleleCounts['X'];
					++alleleCounts['Y'];
					break;
				case 2:
					alleleCounts['Y'] += 2;
					break;
				}
			} // end count alleles
		} // end all instances
		  // cout << alleleCounts['X'] << " => " << alleleCounts['Y'] << endl;
		attributeAlleleCounts.push_back(alleleCounts);
		/// assign major and minor alleles
		double maf = 0;
		char minorAllele = 'Y';
		if (alleleCounts['X'] > alleleCounts['Y']) {
			attributeAlleles.push_back(make_pair('X', 'Y'));
			maf = (double) alleleCounts['Y'] / ((double) NumInstances() * 2.0);
			minorAllele = 'Y';
		} else {
			attributeAlleles.push_back(make_pair('Y', 'X'));
			maf = (double) alleleCounts['X'] / ((double) NumInstances() * 2.0);
			minorAllele = 'X';
		}
		attributeMinorAllele.push_back(make_pair(minorAllele, maf));
		attributeMutationTypes.push_back(UNKNOWN_MUTATION);
	} // end all attributes

	hasAllelicInfo = true;
}

bool Dataset::LoadNumerics(string filename) {
	numericsFilename = filename;
	ifstream dataStream(numericsFilename.c_str());
	if (!dataStream.is_open()) {
		cerr << "ERROR: Could not open numerics file: " << numericsFilename
				<< endl;
		return false;
	}
	cout << Timestamp() << "Reading numerics from " << numericsFilename << endl;
	//  PrintStats();

	// temporary string for reading file lines
	uint lineNumber = 0;
	string line;
	// read the header for covariate names
	getline(dataStream, line);
	++lineNumber;
	vector<string> numNames;
	split(numNames, line);
	if (numNames.size() < 3) {
		cerr << "ERROR: Covariate file must have at least three columns: "
				<< "FID IID COV1 ... COVN" << endl;
		return false;
	}
	numericsNames.resize(numNames.size() - 2);
	copy(numNames.begin() + 2, numNames.end(), numericsNames.begin());
	vector<string>::const_iterator it = numericsNames.begin();
	uint numIdx = 0;
	for (; it != numericsNames.end(); ++it) {
		numericsMask[*it] = numIdx;
		++numIdx;
	}

	// if no snp data then need to create instances in this loop - 6/19/11
	// read each new set of numerics
	map<string, bool> idsSeen;
	uint newInstanceIdx = 0;
	DatasetInstance* tempInstance = 0;
	uint instanceIndex = 0;
	while (getline(dataStream, line)) {
		++lineNumber;

		// cout << lineNumber << ": " << line << endl;
		vector<string> numericsStringVector;
		split(numericsStringVector, line);
		// cout << line << endl;
		if (numericsStringVector.size() < 3) {
			cerr << "ERROR: Covariate file must have at least three columns: "
					<< "FID IID COV1 ... COVN" << endl;
			return false;
		}
		if (numericsNames.size() != (numericsStringVector.size() - 2)) {
			cerr
					<< "ERROR: Number of numeric values read from the covariate file header: ["
					<< numericsNames.size()
					<< "] is not equal to the number of numeric"
					<< " values read [" << (numericsStringVector.size() - 2)
					<< "] on line: " << lineNumber << " of " << numericsFilename
					<< endl;
			return false;
		}

		// first column is the ID for matching data set rows - bcw - 8/3/11
		// use both FID and IID so all PLINK files will work - 4/10/13
		string ID = numericsStringVector[0] + numericsStringVector[1];
		if (!IsLoadableInstanceID(ID)) {
			cout << Timestamp() << "WARNING: Skipping Numeric ID [" << ID
					<< "]. "
					<< "It does not match the data set and/or phenotype file"
					<< endl;
			continue;
		}

		if (idsSeen.find(ID) == idsSeen.end()) {
			idsSeen[ID] = true;
		} else {
			cout << Timestamp() << "WARNING: Duplicate ID [" << ID
					<< "] detected and " << "skipped on line [" << lineNumber
					<< "]" << endl;
			continue;
		}

		// cout << "Numerics ID string from file: " << thisID << endl;
		numericsIds.push_back(ID);
		if (!hasGenotypes) {
			instancesMask[ID] = newInstanceIdx++;
			tempInstance = new DatasetInstance(this, ID);
		}
		// skip the first two columns: family and individual IDs
		vector<string>::const_iterator it = numericsStringVector.begin() + 2;
		// add new numeric columns to this instance
		uint numericsIndex = 0;
		for (; it != numericsStringVector.end(); it++) {
			//      cout << "instance " << instanceIndex << ", read from file: " << *it
			//              << ", numeric value: " << strtod((*it).c_str(), NULL) << endl;
			NumericLevel thisValue = 0.0;
			if ((*it == "-9") || (*it == "?")) {
				thisValue = MISSING_NUMERIC_VALUE;
				missingNumericValues[ID].push_back(numericsIndex);
			} else {
				thisValue = lexical_cast<NumericLevel>(*it);
			}
			if (!hasGenotypes) {
				tempInstance->AddNumeric(thisValue);
			} else {
				uint lookupIdIndex = 0;
				if (GetInstanceIndexForID(ID, lookupIdIndex)) {
					// cout << "ID: " << thisID << ", Lookup index: " << lookupIdIndex << ", numInstances = " << instances.size() << endl;
					instances[lookupIdIndex]->AddNumeric(thisValue);
				}
			}
			++numericsIndex;
			// cout << "-----------------------------------------------------" << endl;
		}
		if (!hasGenotypes) {
			instances.push_back(tempInstance);
		}
		++instanceIndex;
	}

	hasNumerics = true;

	// find the min and max values for each numeric attribute
	// used in diff/distance calculation metrics
	vector<NumericLevel> numericColumn;
	for (uint i = 0; i < NumNumerics(); ++i) {
		double columnSum = 0.0;
		GetNumericValues(i, numericColumn);
		double minElement = *numericColumn.begin();
		double maxElement = *numericColumn.begin();
		for (vector<NumericLevel>::const_iterator it = numericColumn.begin();
				it != numericColumn.end(); ++it) {
			if ((*it != MISSING_NUMERIC_VALUE) && (*it < minElement)) {
				minElement = *it;
			}
			if ((*it != MISSING_NUMERIC_VALUE) && (*it > maxElement)) {
				maxElement = *it;
			}
			columnSum += *it;
		}
		numericsMinMax.push_back(make_pair(minElement, maxElement));
		numericsSums.push_back(columnSum);
	}

  if(par::verbose) PrintStatsSimple();

	return true;
}

bool Dataset::LoadPrivacySim(string filename) {
  // --------------------------------------------------------------------------
	numericsFilename = filename;
	ifstream dataStream(numericsFilename.c_str());
	if (!dataStream.is_open()) {
		error("ERROR: Could not open sim data file: " + numericsFilename + "\n");
		return false;
	}
	PP->printLOG(Timestamp() + "Reading sim data from " + numericsFilename + "\n");

  // --------------------------------------------------------------------------
	// temporary string for reading file lines
	string line;
	uint lineNumber = 0;
  uint newInstanceIndex = 0;

  // --------------------------------------------------------------------------
	// read the header for variable names
	getline(dataStream, line);
	++lineNumber;
	vector<string> parsedHeaderLine;
	split(parsedHeaderLine, line);
	numericsNames.resize(parsedHeaderLine.size() - 1);
	copy(parsedHeaderLine.begin(), parsedHeaderLine.end() - 1, numericsNames.begin());
	vector<string>::const_iterator it = numericsNames.begin();
	for (uint numIdx = 0; it != numericsNames.end(); ++it, ++numIdx) {
		numericsMask[*it] = numIdx;
	}
  PP->printLOG(Timestamp() + "Read " + int2str(numericsNames.size()) + 
               " numeric names from the file header\n");
  
  // --------------------------------------------------------------------------
	// no snp data in privacy sims, so create instances in this loop - 12/23/16
  hasGenotypes = false;
  hasNumerics = true;
	map<string, bool> idsSeen;
	DatasetInstance* tempInstance = 0;
  uint phenoIdx = numericsNames.size();
  hasPhenotypes = true;
  hasContinuousPhenotypes = false;
	while (getline(dataStream, line)) {
    // parse line by line the matrix of expression and phenotype
		++lineNumber;
		vector<string> parsedNumericsFileLine;
		split(parsedNumericsFileLine, line);
		if (parsedNumericsFileLine.size() < 2) {
      error("ERROR: simulated data file must have at least two columns: "
            "VAR1 ... PHENO");
		}
		if (numericsNames.size() != (parsedNumericsFileLine.size() - 1)) {
			error("Number of numeric values read from the simulated data file header: [" +
            int2str(numericsNames.size()) + "] is not equal to the number of numeric" +
            " values read [" + int2str(parsedNumericsFileLine.size() - 1) +
            "] on line: " + int2str(lineNumber) + " of " + numericsFilename + "\n");
		}

    // new instance information
    stringstream instIDss;
    instIDss << "ind" << (lineNumber - 1);
    string instID = instIDss.str();
    instanceIds.push_back(instID);
    instanceIdsToLoad.push_back(instID);
    
		if (!IsLoadableInstanceID(instID)) {
			PP->printLOG(Timestamp() + "WARNING: Skipping Numeric ID [" + instID	+ 
                   "]. It does not match the data set and/or phenotype file\n");
			continue;
		}
		if (idsSeen.find(instID) == idsSeen.end()) {
			idsSeen[instID] = true;
		} else {
			PP->printLOG(Timestamp() + "WARNING: Duplicate ID [" + instID +
                   "] detected and skipped on line [" + int2str(lineNumber) +
                   + "]\n");
			continue;
		}

    tempInstance = new DatasetInstance(this, instID);
		vector<string>::const_iterator it = parsedNumericsFileLine.begin();
    numericsIds.clear();
		for (uint numericsIndex = 0; 
            it != parsedNumericsFileLine.end(); 
            it++, numericsIndex++) {
  
      if(numericsIndex != phenoIdx) {
        string varName = numericsNames[numericsIndex];
//        cout << "instance index|numerics index|numerics name|numerics value: " 
//                << newInstanceIndex << "|" 
//                << numericsIndex << "|"
//                << varName << "|"
//                << *it << endl;
        numericsIds.push_back(varName);
        NumericLevel thisValue = 0.0;
        if ((*it == "-9") || (*it == "?")) {
          thisValue = MISSING_NUMERIC_VALUE;
          missingNumericValues[instID].push_back(numericsIndex);
        } else {
          thisValue = lexical_cast<NumericLevel>(*it);
        }
        tempInstance->AddNumeric(thisValue);
      } else {
        ClassLevel thisClass = lexical_cast<ClassLevel>(*it);
        // convert -1, 1 encoding to 0/1 phenotype for inbix expectations
        thisClass = (thisClass == -1)? 0: 1;
        tempInstance->SetClass(thisClass);
        classIndexes[thisClass].push_back(newInstanceIndex);
      }
    
    }
		instances.push_back(tempInstance);
    instancesMask[instID] = newInstanceIndex;
    ++newInstanceIndex;
	} // end read line by line from file
	PP->printLOG(Timestamp() + "Read " + int2str(NumNumerics()) + 
          " simulated data attributes\n");

	// find the min and max values for each numeric attribute
	// used in diff/distance calculation metrics
	vector<NumericLevel> numericColumn;
	for (uint i = 0; i < NumNumerics(); ++i) {
		double columnSum = 0.0;
		GetNumericValues(i, numericColumn);
		double minElement = *numericColumn.begin();
		double maxElement = *numericColumn.begin();
		for (vector<NumericLevel>::const_iterator it = numericColumn.begin();
				it != numericColumn.end(); ++it) {
			if ((*it != MISSING_NUMERIC_VALUE) && (*it < minElement)) {
				minElement = *it;
			}
			if ((*it != MISSING_NUMERIC_VALUE) && (*it > maxElement)) {
				maxElement = *it;
			}
			columnSum += *it;
		}
		numericsMinMax.push_back(make_pair(minElement, maxElement));
		numericsSums.push_back(columnSum);
	}
  
  if(par::verbose) PrintStatsSimple();
  
	return true;
}

bool Dataset::GetNumericValues(uint numericIndex,	vector<double>& numericValues) {
	if (!hasNumerics) {
    cout << Timestamp()
    << "WARNING: Dataset::GetNumericValues attempting to access "
    << "numeric data when none have been loaded" << endl;
		return false;
	}
  
  if (numericIndex > numericsNames.size()) {
    cerr
        << "ERROR: Dataset::GetNumericValues: numeric index out of range: "
        << numericIndex << endl;
    return false;
  }
  if (!NumInstances()) {
    cerr << "ERROR: Dataset::GetNumericValues: no instances/subjects found" << endl;
    return false;
  }

  // for numericIndex variable in the variables mask
  string thisVarName = numericsNames[numericIndex];
  uint thisVarIdx = numericsMask[thisVarName];

  if(!MaskSearchVariableType(thisVarName, NUMERIC_TYPE)) {
    cout << "WARNING SKIPPING NUMERIC VARIABLE: " 
            << numericIndex << ", " 
            << thisVarName << endl;
  }

    // for all instances in the instance mask
  numericValues.clear();
  for (uint instOrdIdx=0; instOrdIdx < instanceIds.size(); ++instOrdIdx) {
    string thisInstanceID = instanceIds[instOrdIdx];
    uint thisInstanceIdx = instancesMask[thisInstanceID];
    NumericLevel thisNumeric = instances[thisInstanceIdx]->GetNumeric(thisVarIdx);
//    cout << "INSTANCE ORDIDX INSTIDX INSTID: " 
//            << instOrdIdx << ", " 
//            << thisInstanceIdx << ", "
//            << thisInstanceID << ", "
//            << thisNumeric
//            << endl;
    numericValues.push_back(thisNumeric);
  }

	return true;
}

bool Dataset::LoadAlternatePhenotypes(string phenotypesFilename) {
	if ((!hasGenotypes) && (!hasNumerics)) {
		cerr
				<< "ERROR: Dataset::LoadAlternatePhenotypes: SNP and/or numeric data "
				<< "must be loaded before alternate phenotypes" << endl;
		return false;
	}

	/// Detect the class type
	bool classDetected = false;
	switch (DetectClassType(phenotypesFilename, 3, false)) {
	case CASE_CONTROL_CLASS_TYPE:
		cout << Timestamp() << "Case-control phenotypes detected" << endl;
		hasContinuousPhenotypes = false;
		classDetected = true;
		break;
	case CONTINUOUS_CLASS_TYPE:
		cout << Timestamp() << "Continuous phenotypes detected" << endl;
		hasContinuousPhenotypes = true;
		classDetected = true;
		break;
	case MULTI_CLASS_TYPE:
		cout << Timestamp() << "Multiclass phenotypes detected" << endl;
		hasContinuousPhenotypes = false;
		classDetected = true;
		break;
	case NO_CLASS_TYPE:
		cout << "ERROR: phenotypes could not be detected" << endl;
		break;
	}
	if (!classDetected) {
		return false;
	}

	alternatePhenotypesFilename = phenotypesFilename;
	ifstream dataStream(alternatePhenotypesFilename.c_str());
	if (!dataStream.is_open()) {
		cerr << endl << "ERROR: Could not open alternate phenotype file: "
				<< alternatePhenotypesFilename << endl;
		return false;
	}
	cout << Timestamp() << "Reading alternate phenotypes from "
			<< alternatePhenotypesFilename << "... " << endl;

	// temporary string for reading file lines
	string line;
	uint lineNumber = 0;

	// read each new phenotype value
	// remove header? decided we will not use a header
	// getline(dataStream, line);

	// clear any existing class info for reading alternate phenotype file
	// PrintClassIndexInfo();
	if (!hasContinuousPhenotypes) {
		classIndexes.clear();
	}

	uint instancesRead = 0;
	map<string, bool> idsSeen;
	double minPheno = 0.0, maxPheno = 0.0;
	vector<string> idsToDelete;
	while (getline(dataStream, line)) {
		++lineNumber;
		string trimmedLine = trim(line);
		if (trimmedLine == "") {
			cout << Timestamp() << "WARNING: Line [" << lineNumber
					<< "] in phenotype file is blank. Skipping" << endl;
			continue;
		}

		vector<string> lineParts;
		split(lineParts, trimmedLine);

		// use both FID and IID so all PLINK files will work - 4/10/13
		// string ID = lineParts[0];
		string ID = lineParts[0] + lineParts[1];
		// skip IDs that don't match common IDs between numerics, phenotypes and data
		if (!IsLoadableInstanceID(ID)) {
			cout << Timestamp() << "WARNING: Skipping alternate phenotype ID ["
					<< ID << "]. "
					<< "It does not match the data set and/or numerics file(s)"
					<< endl;
			continue;
		}
		// skip duplicate IDs
		if (idsSeen.find(ID) == idsSeen.end()) {
			idsSeen[ID] = true;
		} else {
			cout << Timestamp() << "WARNING: Duplicate ID [" << ID
					<< "] detected and " << "skipped on line [" << lineNumber
					<< "]" << endl;
			continue;
		}

		// lookup instance index for this ID; cannot assume in same order
		uint instanceIndex = 0;
		if (!GetInstanceIndexForID(ID, instanceIndex)) {
			cerr << "ERROR: on lookup in GetInstanceIndexForID: " << ID
					<< " on line: " << lineNumber << endl;
			return false;
		}

		ClassLevel classValue = 0;
		NumericLevel predictedValue = 0.0;
		string classString = lineParts[2];
		if ((classString == "-9") || (classString == "?")) {
			cout << Timestamp() << "WARNING: missing phenotype value read from "
					<< "alternate phenotype file line: " << lineNumber << endl;
			cout << Timestamp() << "\tID: " << ID << " will be removed from the"
					" analysis." << endl;
			idsToDelete.push_back(ID);
			continue;
		} else {
			if (hasContinuousPhenotypes) {
				predictedValue = lexical_cast<NumericLevel>(classString);
				instances[instanceIndex]->SetPredictedValueTau(predictedValue);
				if (lineNumber == 1) {
					minPheno = maxPheno = predictedValue;
				} else {
					if (predictedValue < minPheno) {
						minPheno = predictedValue;
					}
					if (predictedValue > maxPheno) {
						maxPheno = predictedValue;
					}
				}
			} else {
				classValue = lexical_cast<ClassLevel>(classString);
				instances[instanceIndex]->SetClass(classValue);
				classIndexes[classValue].push_back(instanceIndex);
			}
		}

		++instancesRead;
	}

	// if we found missing genotypes in the alternate phenotype file load,
	// we need to remove the already loaded instances with the missing phenotype
	// instance IDs - also update meta data - 10/27/11
	if (idsToDelete.size()) {
		cout << Timestamp() << "Removing " << idsToDelete.size()
				<< " instances with missing phenotypes" << endl;
		for (vector<string>::const_iterator delIt = idsToDelete.begin();
				delIt != idsToDelete.end(); ++delIt) {

			string delId = *delIt;
			uint delIdIndex = instancesMask[delId];
			ClassLevel delClass = instances[delIdIndex]->GetClass();

			// remove instanceIndex from classIndexes
			vector<uint>::iterator ciIt = find(
					classIndexes[delClass].begin(),
					classIndexes[delClass].end(), delIdIndex);
			if (ciIt != classIndexes[delClass].end()) {
				classIndexes[delClass].erase(ciIt);
			}

			// remove ID from instanceIds
			vector<string>::iterator iiIt = find(instanceIds.begin(),
					instanceIds.end(), delId);
			if (iiIt != instanceIds.end()) {
				instanceIds.erase(iiIt);
			}

			// remove from instancesMask
			map<string, uint>::iterator imIt = instancesMask.find(
					delId);
			if (imIt != instancesMask.end()) {
				instancesMask.erase(imIt);
			}
		}
		classIndexes.erase(MISSING_DISCRETE_CLASS_VALUE);
	}

	// update all counts
	UpdateAllLevelCounts();

	if (hasContinuousPhenotypes) {
		continuousPhenotypeMinMax = make_pair(minPheno, maxPheno);
	}

	hasAlternatePhenotypes = true;

	cout << Timestamp() << "Read " << instancesRead << " phenotypes. "
			<< NumInstances() << " instances in the data set" << endl;

	return true;
}

bool Dataset::IsLoadableInstanceID(std::string ID) {
	if (!instanceIdsToLoad.size()) {
		return true;
	}
	if (std::find(instanceIdsToLoad.begin(), instanceIdsToLoad.end(), ID)
			== instanceIdsToLoad.end()) {
		return false;
	} else {
		return true;
	}
}

bool Dataset::WriteNewPlinkPedDataset(string baseDatasetFilename) {
	// only write ped format if it makes sense
	if (!hasAllelicInfo) {
		cerr << "ERROR: WriteNewPlinkPedDataset: "
				<< "No allelic information to write new data set" << endl;
		return false;
	}
	if (hasNumerics) {
		cerr << "ERROR: WriteNewPlinkPedDataset: "
				<< "This data set type does not handle numeric attributes"
				<< endl;
		return false;
	}

	// write map file
	string mapFilename = baseDatasetFilename + ".map";
	ofstream newMapStream(mapFilename.c_str());
	if (!newMapStream.is_open()) {
		cerr << "ERROR: Could not open new map file: " << baseDatasetFilename
				<< ".map" << endl;
		return false;
	}
	map<string, uint>::const_iterator ait = attributesMask.begin();
	for (; ait != attributesMask.end(); ++ait) {
		newMapStream << "0 " << ait->first << " 0 0" << endl;
	}
	newMapStream.close();

	// write ped file
	string pedFilename = baseDatasetFilename + ".ped";
	ofstream newPedStream(pedFilename.c_str());
	if (!newPedStream.is_open()) {
		cerr << "ERROR: Could not open new ped file: " << baseDatasetFilename
				<< ".ped" << endl;
		return false;
	}
	vector<string> instanceIds = GetInstanceIds();
	for (uint iIdx = 0; iIdx < NumInstances(); iIdx++) {
		unsigned instanceIndex = 0;
		string instanceId = instanceIds[iIdx];
		GetInstanceIndexForID(instanceId, instanceIndex);
		// write first six columns:
		/*
		 * Family ID
		 * Individual ID
		 * Paternal ID
		 * Maternal ID
		 * Sex (1=male; 2=female; other=unknown)
		 * Phenotype
		 */
		newPedStream << instanceId << " " << instanceId << " 0 0 0 ";
		// TODO: check for missing phenotype?
		if (hasContinuousPhenotypes) {
			newPedStream << instances[instanceIndex]->GetPredictedValueTau()
					<< " ";
		} else {
			if (instances[instanceIndex]->GetClass()
					== MISSING_DISCRETE_CLASS_VALUE) {
				newPedStream << MISSING_DISCRETE_CLASS_VALUE;
			} else {
				newPedStream << (instances[instanceIndex]->GetClass() + 1);
			}
		}
		// write discrete attribute SNP values as pairs of alleles
		vector<uint> attrIndices = MaskGetAttributeIndices(
				DISCRETE_TYPE);
		for (uint aIdx = 0; aIdx < attrIndices.size(); aIdx++) {
			AttributeLevel A = instances[instanceIndex]->GetAttribute(
					attrIndices[aIdx]);
			// determine alleles and genotype
			pair<char, char> alleles = GetAttributeAlleles(attrIndices[aIdx]);
			stringstream genotype;
			switch (A) {
			case 0:
				genotype << alleles.first << " " << alleles.first;
				break;
			case 1:
				genotype << alleles.second << " " << alleles.first;
				break;
			case 2:
				genotype << alleles.second << " " << alleles.second;
				break;
			case MISSING_ATTRIBUTE_VALUE:
				genotype << "0 0";
				break;
			}
			newPedStream << " " << genotype.str();
		}
		newPedStream << endl;
	}
	newPedStream.close();

	return true;
}

bool Dataset::WriteNewPlinkCovarDataset(string baseDatasetFilename) {
	// not sure I need this
	return true;
}
