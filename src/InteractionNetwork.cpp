/*
 * InteractionNetwork.cpp - Bill White - 12/1/12
 *
 * In Silico Lab interaction network class.
 * 
 * Imported into inbix and removed external library dependencies: 
 * boost, armadillo and igraph. - bcw - 5/13/13
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <cmath>
#include <numeric>

#include "plink.h"
#include "helper.h"
#include "stats.h"
#include "StringUtils.h"

#include "InteractionNetwork.h"

using namespace std;
using namespace Insilico;

InteractionNetwork::InteractionNetwork(string matrixFileParam,
		 MatrixFileType fileType, bool isUpperTriangular, Plink* pp)
{
  // provide access to the inbix (PLINK) environment)
  inbixEnv = pp;
  
	// read the matrix file and store in adjacencyMatrix member variable
	switch(fileType) {
	case REGAIN_FILE:
		if(!ReadGainFile(matrixFileParam, isUpperTriangular)) {
			error("FATAL ERROR: Reading (re)GAIN file: " + matrixFileParam + "\n");
		}
		break;
	case CORR_1D_FILE:
		if(!ReadBrainCorr1DFile(matrixFileParam)) {
			error("FATAL ERROR: Reading matrix file: " + 
							matrixFileParam + "\n");
		}
		break;
	case CSV_FILE:
		if(!ReadCsvFile(matrixFileParam)) {
			error("FATAL ERROR: Reading matrix file: " + matrixFileParam + "\n");
		}
		break;
	case SIF_FILE:
		if(!ReadSifFile(matrixFileParam)) {
			error("FATAL ERROR: Reading SIF file: " + matrixFileParam + "\n");
		}
		break;
	default:
		error("Could not determine the matrix file type: " + 
						int2str(fileType) + "\n");
	}
	networkFile = matrixFileParam;
  
  edgeThreshold = DEFAULT_EDGE_THRESHOLD;
}

InteractionNetwork::InteractionNetwork(double** variablesMatrix,
		unsigned int dim, vector<string>& variableNames, Plink* pp)
{
  // provide access to the inbix (PLINK) environment)
  inbixEnv = pp;

	// setup G matrix for SNPrank algorithm
	sizeMatrix(adjMatrix, dim, dim);
	// copy variable matrix values into G
	for(unsigned int i=0; i < dim; ++i) {
		for(unsigned int j=0; j < dim; ++j) {
			adjMatrix[i][j] = adjMatrix[j][i] = variablesMatrix[i][j];
		}
	}
	for(unsigned int i=0; i < dim; ++i) {
		nodeNames.push_back(variableNames[i]);
	}
	// set default values
	networkFile = "";
  edgeThreshold = DEFAULT_EDGE_THRESHOLD;
}

InteractionNetwork::~InteractionNetwork()
{}

bool InteractionNetwork::SetThreshold(double threshold) {
  edgeThreshold = threshold;
  
  return true;
}

unsigned int InteractionNetwork::NumNodes()
{
	return adjMatrix.size();
}

matrix_t InteractionNetwork::GetAdjacencyMatrix() {
	return adjMatrix;
}

vector<string> InteractionNetwork::GetNodeNames() {
	return nodeNames;
}

void InteractionNetwork::PrintAdjacencyMatrix() {
	for(unsigned int i=0; i < nodeNames.size(); ++i) {
		cout << setw(12) << nodeNames[i];
	}
	cout << endl;
	for(unsigned int i=0; i < adjMatrix.size(); ++i) {
		for(unsigned int j=0; j < adjMatrix.size(); ++j) {
			if(j <= i) {
				printf("%8.6f\t", adjMatrix[i][j]);
			}
		}
		cout << endl;
	}
}

void InteractionNetwork::PrintSummary()
{
	unsigned int n = adjMatrix.size();
	inbixEnv->printLOG("Adjacency Matrix: "	+ int2str(n) + " x " + 
    int2str(n) + "\n");
	double minElement = adjMatrix[0][0];
	double maxElement = adjMatrix[0][0];
  for(int i=0; i < adjMatrix.size(); ++i) {
    for(int j=0; j < adjMatrix.size(); ++j) {
      if(adjMatrix[i][j] < minElement) {
        minElement = adjMatrix[i][j];
      }
      if(adjMatrix[i][j] > maxElement) {
        maxElement = adjMatrix[i][j];
      }
    }
  }
  inbixEnv->printLOG("Minimum: " + dbl2str(minElement) + "\n");
  inbixEnv->printLOG("Maximum: " + dbl2str(maxElement) + "\n");
  if(edgeThreshold > 0) {
    inbixEnv->printLOG("Edge Threshold: " + dbl2str(edgeThreshold) + "\n");
  } else {
    inbixEnv->printLOG("Edge Threshold: NO THRESHOLD\n");
  }
}

bool InteractionNetwork::WriteToFile(string outFile, MatrixFileType fileType)
{
	bool success = false;
  
	switch(fileType) {
	case CSV_FILE:
		success = WriteDelimitedFile(outFile, ",");
		break;
	case REGAIN_FILE:
		success = WriteDelimitedFile(outFile, "\t");
		break;
	case SIF_FILE:
		success = WriteSifFile(outFile);
		break;
	default:
		cerr << "InteractionNetwork::WriteToFile: "
				<< "ERROR: Unknown file type: " << fileType << endl;
	}

	return success;
}

bool InteractionNetwork::WriteDelimitedFile(string outFilename, string delimiter)
{
	ofstream outputFileHandle(outFilename.c_str());
	for(unsigned int i=0; i < nodeNames.size(); ++i) {
		if(i) {
			outputFileHandle << delimiter << nodeNames[i];
		}
		else {
			outputFileHandle << nodeNames[i];
		}
	}
	outputFileHandle << endl << setiosflags(ios::fixed) << setprecision(8);
	for(unsigned int i=0; i < adjMatrix.size(); ++i) {
		for(unsigned int j=0; j < adjMatrix.size(); ++j) {
			if(j) {
				outputFileHandle << delimiter << adjMatrix[i][j];
			}
			else {
				outputFileHandle << adjMatrix[i][j];
			}
		}
		outputFileHandle << endl;
	}
	outputFileHandle.close();

	return true;
}

bool InteractionNetwork::WriteSifFile(string outFilename)
{
	ofstream outputFileHandle(outFilename.c_str());
	for(unsigned int i=0; i < adjMatrix.size(); ++i) {
		for(unsigned int j=i+1; j < adjMatrix.size(); ++j) {
			if(adjMatrix[i][j]) {
				outputFileHandle
					<< nodeNames[i] << "\t" << adjMatrix[i][j] << nodeNames[j] << endl;
			}
		}
	}
	outputFileHandle.close();

	return true;
}

bool InteractionNetwork::Merge(
		InteractionNetwork& toMerge,
		double priorProbEdges,
		double alpha,
		double omega,
		double threshold
		)
{
	if(toMerge.NumNodes() != adjMatrix.size()) {
		cerr << "ERROR: Cannot merge networks of different sizes." << endl;
		return false;
	}
	matrix_t otherAdjacencyMatrix = toMerge.GetAdjacencyMatrix();

	double posteriorProb = 0.0;
	for(unsigned int i=0; i < adjMatrix.size(); ++i) {
		for(unsigned int j=i; j < adjMatrix.size(); ++j) {
//			if(i == j) {
//				adjMatrix[i][j] = 0;
//				continue;
//			}
			double beta_ij_1 = adjMatrix[i][j];
			double beta_ij_2 = otherAdjacencyMatrix[i][j];
			double probWgE1 = alpha * (1.0 - exp(-omega * beta_ij_1));
			double probWgE2 = alpha * (1.0 - exp(-omega * beta_ij_2));
			posteriorProb = probWgE1 * probWgE2 * priorProbEdges;
//			cout
//				<< "B_ij_1: " << beta_ij_1 << " "
//				<< "B_ij_2: " << beta_ij_2 << " "
//				<< "P(W1|E1): " << probWgE1 << " "
//				<< "P(W1|E1): " << probWgE2 << " "
//				<< "posterior: " << posteriorProb
//				<< endl;
			if(posteriorProb > threshold) {
				adjMatrix[i][j] = posteriorProb;
				adjMatrix[j][i] = posteriorProb;
			}
			else {
				adjMatrix[i][j] = 0;
				adjMatrix[j][i] = 0;
			}
		}
	}

	return true;
}
// ------------------ P R I V A T E   M E T H O D S --------------------------

bool InteractionNetwork::ReadCsvFile(string matrixFilename)
{
  ifstream matrixFileHandle(matrixFilename.c_str());
  if(!matrixFileHandle.is_open()) {
      cerr << "ERROR: Could not open matrix file" << matrixFilename << endl;
      return false;
  }

  string line;
  string delimiter = ",";

  // read first line - header with node names
  getline(matrixFileHandle, line);
  string trimmedLine = trim(line);
  vector<string> lineParts;
  split(lineParts, trimmedLine, delimiter);
  for(unsigned int nn=0; nn < lineParts.size(); ++nn) {
  	nodeNames.push_back(lineParts[nn]);
  	nodeNameIndex[lineParts[nn]] = nn;
  }

  // initialize dimensions of the adjacency matrix
  size_t adjDim = lineParts.size();
  if(!adjDim) {
      cerr << "ERROR: Could not parse header values" << endl;
      return false;
  }
  sizeMatrix(adjMatrix, adjDim, adjDim);

  // read the rest of the file as adjacency matrix values
  unsigned int row = 0;
  unsigned int tokensExpected = adjDim;
  while(getline(matrixFileHandle, line)) {
    string trimmedLine = trim(line);
    lineParts.clear();
    split(lineParts, trimmedLine, delimiter);
    if(lineParts.size() != tokensExpected) {
      cerr << "ERROR line:" << endl << endl << line << endl << endl
      		<< "ERROR: Could not parse file row: " << (row+1) << endl
      		<< "Expecting " << tokensExpected << " values, got "
      		<< lineParts.size() << endl;
      return false;
    }
    // set the matrix values to the parsed row values
    unsigned int col = 0;
      for(; col < lineParts.size(); ++col) {
        string trimmedPart = trim(lineParts[col]);
				double t;
				if(!from_string<double>(t, lineParts[col], std::dec)) {
					error("Parsing CSV line " + line);
				}
        adjMatrix[row][col] = t;
      }
    row++;
  }
  matrixFileHandle.close();

  return true;
}

bool InteractionNetwork::ReadGainFile(string gainFilename, bool isUpperTriangular)
{
	ifstream gainFileHandle(gainFilename.c_str());
	if(!gainFileHandle.is_open()) {
		cerr << "ERROR: Could not open (re)GAIN file: " << gainFilename << endl;
		return false;
	}

	string line;

	// read first line (header)
	getline(gainFileHandle, line);
	trim(line);

	// tokenize header
  vector<string> lineParts;
  split(lineParts, line);
  for(unsigned int nn=0; nn < lineParts.size(); ++nn) {
  	nodeNames.push_back(lineParts[nn]);
  	nodeNameIndex[lineParts[nn]] = nn;
  }

	// initialize dimensions of data matrix
	size_t numVars = nodeNames.size();
	if(!numVars) {
		cerr << "ERROR: Could not parse SNP names from (re)GAIN file header" << endl;
		return false;
	}
	sizeMatrix(adjMatrix, numVars, numVars);

	// read numeric data into G
	size_t row = 0;
	vector<string> lineTokens;
	unsigned int tokensExpected = numVars;
	while(getline(gainFileHandle, line)) {
		trim(line);
		lineTokens.clear();
		split(lineTokens, line);
		if(lineTokens.size() != tokensExpected) {
			cerr << "ERROR line:" << endl << endl << line << endl << endl;
			cerr << "ERROR: Could not parse (re)GAIN file row: " << (row+2) << endl
					<< "Expecting " << tokensExpected << " values, got "
					<< lineTokens.size() << endl;
			return false;
		}
		size_t startIndex = numVars - tokensExpected;
		for (size_t tokenIndex=0, col=startIndex; col < numVars; ++tokenIndex, ++col) {
			string token = lineTokens[tokenIndex];
			trim(token);
			if(token == "") {
				cerr << "ERROR: parsing line " << row << " col " << col
						<< " of (re)Gain file, token [" << token << "]" << endl;
				return false;
			}
			double t;
			if(!from_string<double>(t, lineParts[col], std::dec)) {
				error("Parsing REGAIN line " + line);
			}
			adjMatrix[row][col] = t;
			if((row != col) && isUpperTriangular) {
				if(!from_string<double>(t, lineParts[col], std::dec)) {
					error("Parsing REGAIN line " + line);
				}
				adjMatrix[col][row] = t;
			}
		}
		row++;
		if(isUpperTriangular) {
			--tokensExpected;
		}
	}

	gainFileHandle.close();

	return true;
}

bool InteractionNetwork::ReadSifFile(string sifFilename)
{
	// parse the graph data structures

	// unique node names
  set<string> nodeNameSet;
  // edges between nodes
  vector<pair<pair<string, string>, double> > edges;

  // parse the SIF text file
  ifstream sifFileHandle(sifFilename.c_str());
  if(!sifFileHandle.is_open()) {
      cerr << "ERROR: Could not open SIF file" << sifFilename << endl;
      return false;
  }

  string line;
  string delimiter = "\t";
  while(getline(sifFileHandle, line)) {
    trim(line);
    if(line == "") {
    	cout << "WARNING: Blank line skipped" << endl;
    	continue;
    }
    vector<string> sifValues;
    split(sifValues, line, delimiter);
    string node1 = sifValues[0];
		double t;
		if(!from_string<double>(t, sifValues[1], std::dec)) {
			error("Parsing SIF line " + line);
		}
    double weight = t;
    string node2 = sifValues[2];
    nodeNameSet.insert(node1);
    nodeNameSet.insert(node2);
    edges.push_back(make_pair(make_pair(node1, node2), weight));
  }
  sifFileHandle.close();

  // assign each node name a unique index
  map<string, unsigned int> nodeNameMap;
  set<string>::const_iterator nnsIt = nodeNameSet.begin();
  for(unsigned int nnIndex=0;
  		nnsIt != nodeNameSet.end();
  		++nnsIt, ++nnIndex) {
  	nodeNameMap[*nnsIt] = nnIndex;
  	nodeNames.push_back(*nnsIt);
  	nodeNameIndex[*nnsIt] = nnIndex;
  }

  // set symmetric adjacency matrix for the edges
  sizeMatrix(adjMatrix, nodeNameSet.size(), nodeNameSet.size());
  matrixFill(adjMatrix, 0.0);
  vector<pair<pair<string, string>, double> >::const_iterator edgeIt;
  for(edgeIt = edges.begin(); edgeIt != edges.end(); ++edgeIt) {
  	pair<string, string> nodeNames = edgeIt->first;
  	double weight = edgeIt->second;
  	unsigned int node1Index = nodeNameMap[nodeNames.first];
  	unsigned int node2Index = nodeNameMap[nodeNames.second];
  	adjMatrix[node1Index][node2Index] = weight;
  	adjMatrix[node2Index][node1Index] = weight;
  }

  return true;
}

bool InteractionNetwork::ReadBrainCorr1DFile(string corr1dFilename) {

	ifstream corr1dFileHandle(corr1dFilename.c_str());
  if(!corr1dFileHandle.is_open()) {
      cerr << "ERROR: ReadBrainCorr1DFile: Could not open "
      		<< corr1dFilename << endl;
      return false;
  }

  string line;
  string delimiter = " ";
  unsigned int row = 0;

  // read header line
  getline(corr1dFileHandle, line);
  trim(line);
  // strip off first the two characters: #<space>
  string header = line.substr(2, line.size()-2);
  vector<string> headerValues;
  split(headerValues, header);
  unsigned int adjDim = headerValues.size();
  if(!adjDim) {
		cerr << "ERROR: ReadBrainCorr1DFile: Could not parse 1D correlation values"
				<< endl;
		return false;
  }
  sizeMatrix(adjMatrix, adjDim, adjDim);

  vector<string>::const_iterator hIt = headerValues.begin();
  unsigned int hIndex = 0;
  for(; hIt != headerValues.end(); ++hIt, ++hIndex) {
  	nodeNames.push_back(*hIt);
  	nodeNameIndex[*hIt] = hIndex;
  }

  // read the rest of the matrix file
  vector<string> corr1dValues;
  unsigned int tokensExpected = adjDim;
  while(getline(corr1dFileHandle, line)) {
    trim(line);
    corr1dValues.clear();
    split(corr1dValues, line, delimiter);
    if(corr1dValues.size() != tokensExpected) {
      cerr << "ERROR line:" << endl << endl << line << endl << endl;
      cerr << "ERROR: Could not parse 1D correlation file row: "
        << (row+1) << endl << "Expecting " << tokensExpected
        << " values, got " << corr1dValues.size() << endl;
      return false;
    }
    // set the matrix values to the parsed row values
  	unsigned int col=0;
		for(; col < corr1dValues.size(); ++col) {
			string trimmedCor = trim(corr1dValues[col]);
			double t;
			if(!from_string<double>(t, trimmedCor, std::dec)) {
				error("Parsing Corr1D line " + line);
			}
			adjMatrix[row][col] = t;
		}
    row++;
  }
  corr1dFileHandle.close();

	return true;
}

pair<double, vector<vector<unsigned int> > >
	InteractionNetwork::ModularityLeadingEigenvector() {

	// keep original adjacency matrix
	matrix_t A(adjMatrix);
	unsigned int n = A.size();

	// remove the diagonal
  vector_t diagZ(n, 0);
	matrixSetDiag(A, diagZ);

	// make the adjacency matrix binary using a threshold (from user params)
  if(edgeThreshold > 0) {
    for(unsigned int i=0; i < n; ++i) {
      for(unsigned int j=0; j < n; ++j) {
        if(A[i][j] < edgeThreshold) {
          A[i][j] = 0.0;
        }
        else {
          A[i][j] = 1.0;
        }
      }
    }
  }
  
	// get column sums k_i, which correspond to number of edges * 2
  vector_t k;
  matrixSums(A, k, 1);
  double m = 0;
  for(int i=0; i < k.size(); ++i) {
    m += k[i];
  }
	m *= 0.5;

	// real symmetric modularity matrix B
  // B = A - k_vec * k_vec.t() / (2.0 * m);
	matrix_t B;
  sizeMatrix(B, n, n);
  double scaleFactor = 1.0 / (2.0 * m);
  for(int i=0; i < n; ++i) {
    for(int j=0; j < n; ++j) {
      B[i][j] = A[i][j] - k[i] * k[j] * scaleFactor;
    }
  }
  
	// ------------------------- I T E R A T I O N ------------------------------

	stack<vector<unsigned int> > processStack;

	// the starting module is the entire network
	vector<unsigned int> firstModule;
	for(unsigned int i=0; i < n; ++i) {
		firstModule.push_back(i);
	}
	processStack.push(firstModule);

	// iterate until stack is empty
	unsigned int iteration = 1;
	while(!processStack.empty()) {
		++iteration;

		vector<unsigned int> thisModule = processStack.top();
		processStack.pop();
		unsigned int newDim = thisModule.size();

		// get the submatrix Bg defined by the indices of this module (Eqn 6)
		matrix_t Bg;
    sizeMatrix(Bg, newDim, newDim);
		for(unsigned int l1=0; l1 < newDim; ++l1) {
			for(unsigned int l2=0; l2 < newDim; ++l2) {
				Bg[l1][l2] = B[thisModule[l1]][thisModule[l2]];
			}
		}

		// adjust the diagonal
		matrix_t BgRowSumDiag;
    sizeMatrix(BgRowSumDiag, Bg.size(), Bg.size());
    vector_t rowsums;
    matrixSums(Bg, rowsums, 0);
		for(unsigned int i=0; i < rowsums.size(); ++i) {
			Bg[i][i] = Bg[i][i] - rowsums[i];
		}

		// call the community finding/modularity function
		pair<double, vector_t> sub_modules = ModularityBestSplit(Bg, m);
		double deltaQ = sub_modules.first;
		vector_t s = sub_modules.second;

		// find the split indices
		vector<unsigned int> s1;
		vector<unsigned int> s2;
		for(unsigned int mi=0; mi < s.size(); ++mi) {
			if(s[mi] > 0) {
				s1.push_back(thisModule[mi]);
			}
			else {
				s2.push_back(thisModule[mi]);
			}
		}

		// have we hit any stopping criteria?
		if((s1.size() == 0) || (s2.size() == 0)) {
			modules.push_back(thisModule);
			if(Q == 0) {
				Q = deltaQ;
			}
		}
		else {
			if(deltaQ <= MODULARITY_THRESHOLD) {
				modules.push_back(thisModule);
			} else {
				// add the splits to the processing stack and recurse
				processStack.push(s1);
				processStack.push(s2);
				// accumulate global Q
				Q += deltaQ;
			}
		}
	}

	return make_pair(Q, modules);
}

pair<double, vector<double> >	InteractionNetwork::Homophily() {

  if(!modules.size()) {
    error("Cannot compute homphily: no modules exist");
  }
  
	double globalHomophily = 0.0;
	vector<double> localHomophilies;

	unsigned int totalNodes = adjMatrix.size();
//	cout << "Total nodes: " << totalNodes << endl;

	// for each module in the modules list
//	cout << "Number of module: " << modules.size() << endl;
	for(unsigned int curModule=0; curModule < modules.size(); ++curModule) {

		// get the indices of the nodes in the module
		unsigned int modSize = modules[curModule].size();
//		cout << "Module size: " << modSize << endl;

		intvec_t modIndices(modSize);
		for(unsigned int mi=0; mi < modSize; ++mi) {
			modIndices[mi] = modules[curModule][mi];
		}

		// get the indices of the nodes not in the module
		intvec_t notIndices(totalNodes - modSize);
		unsigned int notIndex = 0;
		for(unsigned int j=0; j < modules.size(); ++j) {
			if(j != curModule) {
				for(unsigned int k=0; k < modules[j].size(); ++k) {
					notIndices[notIndex] = modules[j][k];
					++notIndex;
				}
			}
		}
//    cout << "Module indices:" << endl;
//    display(modIndices);
//    cout << "Not in module indices:" << endl;
//    display(notIndices);

		// get the number of internal connections
		// matrix_t modMatrix = adjMatrix(modIndices, modIndices);
    matrix_t modMatrix;
    matrixExtractRowColIdx(adjMatrix, modIndices, modIndices, modMatrix);
		// double internalConnections = sum(sum(trimatu(modMatrix)));
    // sum the column sums of the upper triangular
    double internalConnections = 0;
    for(int j=0; j < modMatrix[0].size(); ++j) {
      for(int i=0; i < j; ++i) {
        internalConnections += modMatrix[i][j];
      }
    }
    
		// get the number of external connections
		matrix_t notMatrix;
    matrixExtractRowColIdx(adjMatrix, modIndices, notIndices, notMatrix);
    vector_t notSums;
    matrixSums(notMatrix, notSums, 1);
		//double externalConnections = sum(sum(notMatrix));
    double externalConnections = 0;
    for(int i=0; i < notSums.size(); ++i) {
      externalConnections += notSums[i];
    }

//		cout << "int: " << internalConnections
//				<< ", ext: " << externalConnections << endl;

		// calculate and save local homophily
		double modHomophily = 0.0;
		if((internalConnections != 0) && (externalConnections != 0)) {
			modHomophily = (internalConnections - externalConnections) /
					(internalConnections + externalConnections);
		}
//		cout << "Module homophily: " << modHomophily << endl;
		double localHomophily = modSize * modHomophily / totalNodes;
//		cout << "Module frac: " << localHomophily << endl;
		localHomophilies.push_back(localHomophily);

		// update global homophily
		globalHomophily += localHomophily;
	}

  pair<double, vector<double> > results;
	results.first = globalHomophily;
	results.second.resize(localHomophilies.size());
	copy(localHomophilies.begin(), localHomophilies.end(), results.second.begin());

	return results;
}

void InteractionNetwork::ShowHomophily() {
  inbixEnv->printLOG("Q from existing modules: " + dbl2str(ComputeQ()) + "\n");
  pair<double, vector<double> > homophily = Homophily();
  inbixEnv->printLOG("Total homophily: " + dbl2str(homophily.first) + "\n");
  vector<double>::const_iterator hIt = homophily.second.begin();
  unsigned int modIdx = 0;
  for(; hIt != homophily.second.end(); ++hIt, ++modIdx) {
    inbixEnv->printLOG("Homophily for module " + int2str(modIdx+1) +
             ": " + dbl2str(*hIt) + "\n");
  }
}

double InteractionNetwork::ComputeQ() {

	intvec_t allModules = FlattenModules();
  // double m = sum(sum(adjMatrix)) / 2.0;
  vector_t rowSums;
  matrixSums(adjMatrix, rowSums, 1);
  double m = 0;
  for(int i=0; i < rowSums.size(); ++i) {
    m += rowSums[i];
  }
  m /= 2.0;

  double q = 0.0;
  vector_t k;
  matrixSums(adjMatrix, k, 1);
  for(int i=0; i < adjMatrix.size(); ++i) {
    for(int j=0; j < adjMatrix.size(); ++j) {
			q += (adjMatrix[i][j] - k[i] * k[j] / (2.0 * m)) *
					((double) (allModules[i] == allModules[j]) - 0.5) * 2.0;
    }
  }
  q = q / (4.0 * m);
  
  return q;
}

bool InteractionNetwork::SetModulesFromFile(string modulesFilename) {
  // parse the module text file
  ifstream modFileHandle(modulesFilename.c_str());
  if(!modFileHandle.is_open()) {
      cerr << "ERROR: Could not open modules file" << modulesFilename << endl;
      return false;
  }

  string line;
  string delimiter = "\t";
  map<string, unsigned int> modMap;
  set<unsigned int> moduleNumbers;
  while(getline(modFileHandle, line)) {
    trim(line);
    if(line == "") {
    	cout << "WARNING: Blank line skipped" << endl;
    	continue;
    }
    vector<string> modValues;
    split(modValues, line);
    // TODO: use inbix conversion here
    int moduleNumber = atoi(modValues[1].c_str());
    moduleNumbers.insert(moduleNumber);
    modMap[modValues[0]] = moduleNumber;
  }
  modFileHandle.close();

  modules.resize(moduleNumbers.size());
  map<string, unsigned int>::const_iterator modIt = modMap.begin();
  for(; modIt != modMap.end(); ++modIt) {
  	string nodeName = modIt->first;
  	unsigned int nodeModule = modIt->second;
  	// get index of node name
  	modules[nodeModule].push_back(nodeNameIndex[nodeName]);
  }

	return true;
}
void InteractionNetwork::ShowModules() {
	inbixEnv->printLOG("Modules:\n");
	for(unsigned int moduleIdx=0; moduleIdx < modules.size(); ++moduleIdx) {
		inbixEnv->printLOG("Nodes in module " + int2str(moduleIdx+1) + ": ");
		for(unsigned int memberIdx=0; memberIdx < modules[moduleIdx].size();
				++memberIdx) {
			inbixEnv->printLOG(nodeNames[modules[moduleIdx][memberIdx]] + " ");
		}
    inbixEnv->printLOG("\n");
	}
}

void InteractionNetwork::SaveModules(string saveFilename) {
	ofstream outputFileHandle(saveFilename.c_str());
  if(outputFileHandle.fail()) {
    error("Could not open network modules file for saving: " + saveFilename + "\n");
  }
  inbixEnv->printLOG("Saving network modules to [" + saveFilename + "]\n");
	for(unsigned int moduleIdx=0; moduleIdx < modules.size(); ++moduleIdx) {
		for(unsigned int memberIdx=0; memberIdx < modules[moduleIdx].size();
				++memberIdx) {
			unsigned int nodeIndex = modules[moduleIdx][memberIdx];
			outputFileHandle << nodeNames[nodeIndex] << "\t" << (moduleIdx + 1) << endl;
		}
	}
	outputFileHandle.close();
}

pair<double, vector_t> 
  InteractionNetwork::ModularityBestSplit(matrix_t& B, double m) {

  // Armadillo version:
  //	eig_sym(eigval, eigvec, B);
  //	uword  maxeig_idx;
  //	eigval.max(maxeig_idx);
  //	colvec s_out = eigvec.col(maxeig_idx);

  // eigenvectors changes B!!!
  matrix_t tempB(B);
  
  // compute eigenvectors and eigenvalues
  Eigen eigen = eigenvectors(tempB);
  vector_t eigval = eigen.d;
  matrix_t eigvec = eigen.z;
  //display(eigvec);
  //display(eigval);
  
  int maxEigIdx = 0;
  double maxEigvalue = eigval[maxEigIdx];
  for(int i=1; i < eigval.size(); ++i) {
    if(eigval[i] > maxEigvalue) {
      maxEigvalue = eigval[i];
      maxEigIdx = i;
    }
  }
  
  // get the max eigenvector into s_out
  vector_t s_out;
  for(int i=0; i < eigvec[0].size(); ++i) {
    s_out.push_back(eigvec[i][maxEigIdx]);
  }
  //display(s_out);
  //exit(1);

  // transform s_out into -1 or 1 values
	for(unsigned int i=0; i < s_out.size(); ++i) {
		if(s_out[i] < 0) {
			s_out[i] = -1;
		}
		else {
			s_out[i] = 1;
		}
	}
  //display(s_out);
  
  //                    (n x 1^T) (n x n) (n x 1)
	//matrix_t Q_matrix_t = s_out.t() * B * s_out;
  vector_t s_outTB;
  s_outTB.resize(B.size(), 0);
  for(int i=0; i < s_out.size(); ++i) {
    for(int j=0; j < B.size(); ++j) {
      s_outTB[i] += s_out[j] * B[j][i];
    }
  }
  double Q = 0.0;
  for(int i=0; i < s_outTB.size(); ++i) {
    Q += s_outTB[i] * s_out[i];
  }
	Q *= (1.0 / (m * 4.0));

	return(make_pair(Q, s_out));
}

intvec_t InteractionNetwork::FlattenModules() {
	intvec_t flatModules(adjMatrix.size());
	for(unsigned int i=0; i < modules.size(); ++i) {
		for(unsigned int j=0; j < modules[i].size(); ++j) {
			flatModules[modules[i][j]] = i;
		}
	}
	return flatModules;
}
