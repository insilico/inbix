/*
 * ImpSnpRank.cpp - Bill White - 11/5/12
 *
 * SNPrank using the IMP database interactions table.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <mysql++/mysql++.h>
#include <mysql++/ssqls.h>

#include "SnpRankNextGen.h"

using namespace std;
namespace po = boost::program_options;

// ----------------------------------------------------------------------------
// The Specialized SQL Structure (SSQLS) feature lets you easily define
// C++ structures that match the form of your SQL tables.

// define the struct the 'interactions' table queries
sql_create_3(interactions,
    1, 3,
    string, gene1,
    string, gene2,
    mysqlpp::sql_double, weight
);
// some shorthand for the interaction records
typedef vector<interactions> ImpList;
typedef vector<interactions>::const_iterator ImpListCIt;

// define the struct the 'count' queries
sql_create_2(counts,
    1, 2,
    string, gene,
    mysqlpp::sql_int, count
);
// some shorthand for the counts records
typedef vector<counts> CountsList;
typedef vector<counts>::const_iterator CountsListCIt;

// ----------------------------------------------------------------------------
// function prototypes
bool ReadGeneList(
		string geneListFilename,
		bool hasMainEffects,
		vector<string>& geneNames,
		vector<double>& geneMainEffects
);
bool FindGeneListInteractions(
		mysqlpp::Connection& conn,
		vector<string> geneNames,
		double confidenceThreshold,
		ImpList& results
);
bool WriteImpListToSifFile(
		string sifFilename,
		ImpList& impList
);
void PrintImpList(
		ImpList& impList
);
void PrintVariablesMatrix(
		double ** variablesMatrix,
		unsigned int dim,
		vector<string>& snpNames
);
bool WriteInteractionsMatrix(
		string matrixFilename,
		double** matrix,
		unsigned int dim,
		vector<string> snpNames);
bool GrowGeneList(
		mysqlpp::Connection& conn,
		vector<string>& geneNames,
		unsigned int maxSize,
		double confidenceThreshold
);
vector<string> FindTopConnectedGenes(
		mysqlpp::Connection& conn,
		string thisGene,
		int topN,
		double confidenceThreshold
);

// ---------------------- M A I N  P R O G R A M ------------------------------

int main(int argc, char *argv[]) {
	// --------------------------------------------------------------------------
	// process the command line arguments
	string geneFile;
	string mainEffectsSource = "graduated";
	unsigned int targetNetworkSize = 0;
	double confidenceThreshold = 0.05;
	string outFilePrefix = "impsnprank";
	po::options_description desc("impsnprank");
	desc.add_options()
		("gene-input-file,i", po::value<string>(&geneFile),
		 "input gene list filename (tab-delimited), REQUIRED")
		("main-effects-source,m", po::value<string>(&mainEffectsSource)->default_value(mainEffectsSource),
		 "main effects come from this source (zero|one_over_n|one_over_rank|graduated|file)")
		("target-network-size,n", po::value<unsigned int>(&targetNetworkSize)->default_value(targetNetworkSize),
			 "target network size (default = size of gene list)")
		("confidence-threshold,t", po::value<double>(&confidenceThreshold)->default_value(confidenceThreshold),
			"confidence threshold/cutoff")
		("output-prefix,o", po::value<string>(&outFilePrefix)->default_value(outFilePrefix),
		 "file prefix for all output")
		("write-interactions-matrix,w",
				"write interactions matrix to <output-prefix>.matrix")
		("verbose,v", "verbose output")
		("help,h", "display this help screen");
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		cout << desc << endl;
		return 1;
	}
	else {
		if(vm.count("gene-input-file")) {
			cout << "Command line options:" << endl;
			cout << "Gene list filename: " << geneFile << endl;
			cout << "Target network size: " << targetNetworkSize << endl;
			cout << "Confidence threshold: " << confidenceThreshold << endl;
			cout << "Output File: " << outFilePrefix << endl << endl;
		}
		else {
			cerr << endl << "ERROR: Gene input file is a REQUIRED parameter" << endl;
			cout << endl << desc << endl;
			return 1;
		}
	}

	// --------------------------------------------------------------------------
	// read in the gene list
	cout << "Reading the gene list file" << endl;
	vector<string> geneNames;
	vector<double> mainEffects;
	bool geneFileHasMainEffects = (mainEffectsSource == "file")? true: false;
	if(!ReadGeneList(geneFile, geneFileHasMainEffects,
			geneNames, mainEffects)) {
		cerr << "ERROR: failed to read gene list file: "
				<< geneFile << endl;
		return 1;
	}
	unsigned int numGenes = geneNames.size();
	if(numGenes < 2) {
		cerr << "ERROR: At least two genes are needed, " << numGenes
				<< " read from gene list file" << endl;
		return 1;
	}
	cout << "Read " << numGenes << " genes from gene list file" << endl;
	if(vm.count("verbose")) {
		for(unsigned int i=0; i < geneNames.size(); ++i) {
			cout << i << "\t" << geneNames[i] << endl;
		}
	}

	// --------------------------------------------------------------------------
	// connect to the IMP data base
	cout << "Connecting to the IMP MySQL database" << endl;
	mysqlpp::Connection conn(false);
	// connect to the database
	if(!conn.connect("imp", "localhost", "bwhite", "b314w159")) {
		cerr << "ERROR: DB connection failed: "
				<< conn.error() << endl;
		return 1;
	}

	// --------------------------------------------------------------------------
	// if necessary "grow" the gene list to targetNetworkSize
	if(targetNetworkSize > numGenes) {
		cout << "Growing the gene list to size: " << targetNetworkSize << endl;
		if(!GrowGeneList(conn, geneNames, targetNetworkSize, confidenceThreshold)) {
			cerr << "ERROR: Could not grow network" << endl;
			return 1;
		}
		else {
			numGenes = geneNames.size();
			cout << "New gene list has " << numGenes << " genes" << endl;
			if(vm.count("verbose")) {
				for(unsigned int i=0; i < numGenes; ++i) {
					cout << i << "\t" << geneNames[i] << endl;
				}
			}
		}
	}
	else {
		cout << "WARNING: Target network size is less than gene list size. "
				<< " Network will be created for the gene list as-is." << endl;
	}

	// --------------------------------------------------------------------------
	// find the connections among the gene list from IMP
	cout << "Finding gene list interactions from IMP" << endl;
	ImpList impList;
	if(!FindGeneListInteractions(conn, geneNames, confidenceThreshold, impList)) {
		cerr << "ERROR: failed to get gene list interactions" << endl;
		return 1;
	}
	if(impList.size() < 1) {
		cerr << "ERROR: SQL query for interactions returned no records." << endl;
		return 1;
	}
	else {
		cout << "SQL query returned " << impList.size()
				<< " interaction records" << endl;
	}
	if(vm.count("verbose")) {
		PrintImpList(impList);
	}

	// --------------------------------------------------------------------------
	// set main effects if not from gene list file
	if(mainEffectsSource == "zero") {
		cout << "Setting main effects to zero" << endl;
		mainEffects.clear();
		for(unsigned int i=0; i < numGenes; ++i) {
			mainEffects.push_back(0.0);
		}
	}
	if(mainEffectsSource == "one_over_n") {
		cout << "Setting main effects to 1/n" << endl;
		mainEffects.clear();
		for(unsigned int i=0; i < numGenes; ++i) {
			mainEffects.push_back(1.0/((double) numGenes));
		}
	}
	if(mainEffectsSource == "one_over_rank") {
		cout << "Setting main effects to 1/rank" << endl;
		mainEffects.clear();
		for(unsigned int i=0; i < numGenes; ++i) {
			mainEffects.push_back(1.0/((double) i + 1.0));
		}
	}
	if(mainEffectsSource == "graduated") {
		cout << "Setting main effects to graduated" << endl;
		mainEffects.clear();
		double n = (double) numGenes;
		double sumN = (n * (n + 1.0)) / 2.0;
		for(int i=(int) numGenes; i >= 1; --i) {
			// cout << i << "\t" << ((double) i / sumN) << endl;
			mainEffects.push_back((double) i / sumN);
		}
	}
	if(vm.count("verbose")) {
		for(unsigned int i=0; i < mainEffects.size(); ++i) {
			cout << i << "\t" << mainEffects[i] << endl;
		}
	}

	// --------------------------------------------------------------------------
	cout << "Creating gene lookup table" << endl;
	map<string, unsigned int> nameLookupTable;
	for(unsigned int i=0; i < numGenes; ++i) {
		nameLookupTable[geneNames[i]] = i;
	}

	// --------------------------------------------------------------------------
	// construct a variable interactions matrix for SNPrank
	cout << "Creating variable interactions matrix" << endl;
	double** variablesMatrix;
	variablesMatrix = new double*[numGenes];
	for(unsigned int i=0; i < numGenes; ++i) {
		variablesMatrix[i] = new double[numGenes];
		for(unsigned int j=0; j < numGenes; ++j) {
			variablesMatrix[i][j] = 0.0;
		}
	}
	for(ImpListCIt it=impList.begin(); it != impList.end(); ++it) {
		unsigned int gene1Index = nameLookupTable[it->gene1];
		unsigned int gene2Index = nameLookupTable[it->gene2];
		variablesMatrix[gene1Index][gene2Index] = 1.0;
		variablesMatrix[gene2Index][gene1Index] = 1.0;
	}
	// add main effects diagonal
	for(unsigned int i=0; i < numGenes; ++i) {
		variablesMatrix[i][i] = mainEffects[i];
	}
	if(vm.count("verbose")) {
		PrintVariablesMatrix(variablesMatrix, numGenes, geneNames);
	}
	if(vm.count("write-interactions-matrix")) {
		cout << "Writing interactions matrix to ["
				<< outFilePrefix + ".matrix" << "]" << endl;
		WriteInteractionsMatrix(outFilePrefix + ".matrix", variablesMatrix,
				numGenes, geneNames);
	}

	// --------------------------------------------------------------------------
	cout << "Calling SNPrankNextGen" << endl;
	SnpRankNextGen snprank(variablesMatrix, numGenes, geneNames);
	if(!snprank.CalculateSNPrank(GAUSS_ELIMINATION)) {
		cerr << "ERROR: CalculateSNPrank() failed." << endl;
		return 1;
	}

	// output verbose results
	if(vm.count("verbose")) {
		cout << "Verbose output requested. SNPrank results:" << endl;
		snprank.WriteToConsole();
	}

	// output snprank results to file
	string snpRankFilename = outFilePrefix + ".snprank";
	cout << "Writing SNPrank results [" << snpRankFilename << "]" << endl;
	snprank.WriteToFile(snpRankFilename);

	// write the interaction information to a Cytoscape SIF file
	string sifFilename = outFilePrefix + ".sif";
	cout << "Writing SIF file [" << sifFilename << "]" << endl;
	if(!WriteImpListToSifFile(sifFilename, impList)) {
		cerr << "ERROR: failed to write SIF file"
				<< sifFilename << endl;
		return 1;
	}

	// --------------------------------------------------------------------------
	// clean up and exit gracefully
	cout << "Cleaning up dynamically allocated interactions matrix" << endl;
	for(unsigned int i=0; i < numGenes; ++i) {
		delete [] variablesMatrix[i];
	}
	delete [] variablesMatrix;

	cout << "Done." << endl << endl;

	return 0;
}

// ---------------------- F U N C T I O N S -----------------------------------

bool ReadGeneList(string geneListFilename, bool hasMainEffects,
		vector<string>& geneNames, vector<double>& geneMainEffects)
{
	ifstream geneListFileHandle(geneListFilename.c_str());
	if(!geneListFileHandle.is_open()) {
		cerr << "ERROR: Could not open gene list file: "
				<< geneListFilename << endl;
		return false;
	}

	string line;
	string delimiter = "\t";
	int row = 0;
	vector<string> lineTokens;
	while(getline(geneListFileHandle, line)) {
		boost::trim(line);
		if(hasMainEffects) {
			lineTokens.clear();
			boost::split(lineTokens, line, boost::is_any_of(delimiter));
			if(lineTokens.size() != 2) {
				cerr << "ERROR line:" << endl << endl << line << endl << endl;
				cerr << "ERROR: Could not parse gene list: " << (row+2) << endl
						<< "Expecting 2 values, got "	<< lineTokens.size() << endl;
				return false;
			}
			geneNames.push_back(lineTokens[0]);
			geneMainEffects.push_back(boost::lexical_cast<double>(lineTokens[1]));
		}
		else {
			geneNames.push_back(line);
			geneMainEffects.push_back(0.0);
		}
	}

	geneListFileHandle.close();

	return true;
}

bool FindGeneListInteractions(mysqlpp::Connection& conn,
		vector<string> geneNames, double confidenceThreshold,
		ImpList& results)
{
	// make gene names list into an SQL IN list
	unsigned int nameCounter = 0;
	string inPredicate = "IN (";
	for(vector<string>::const_iterator gnIt = geneNames.begin();
			gnIt != geneNames.end(); ++gnIt) {
		if(nameCounter) {
			inPredicate += ", '" + *gnIt + "'";
		}
		else {
			inPredicate += "'" + *gnIt + "'";
		}
		++nameCounter;
	}
	inPredicate += ")";

	// select intersection of gene1 and gene2 names for the interactions
	string sql = "SELECT * FROM interactions WHERE ";
	sql += "gene1 " + inPredicate +	" AND gene2 " + inPredicate;
	// and cut off any interactions below the threshold
	sql += " AND weight > ";
	stringstream ss;
	ss << sql << confidenceThreshold << ";";
	string completeSql = ss.str();
#ifdef DEBUG_IMPSNPRANK
	cout << "Executing SQL query: " << endl << completeSql << endl;
#endif

	// execute the query
	try {
		mysqlpp::Query query = conn.query(completeSql);
		query.storein(results);
	}
	catch (const mysqlpp::BadQuery& er) {
		// Handle any query errors
		cerr << "Query error: " << er.what() << endl;
		return false;
	}
	catch (const mysqlpp::Exception& er) {
		// Catch-all for any other MySQL++ exceptions
		cerr << "Error: " << er.what() << endl;
		return false;
	}

	return true;
}

bool WriteImpListToSifFile(string sifFilename, ImpList& impList)
{
	ofstream outFile;
	outFile.open(sifFilename.c_str());
	if(outFile.bad()) {
		cerr << "ERROR: Could not open SIF file " << sifFilename
				<< "for writing" << endl;
		exit(1);
	}
	for(ImpListCIt it=impList.begin(); it != impList.end(); ++it) {
		outFile << it->gene1 << "\t" << it->weight << "\t" << it->gene2 << endl;
	}
	outFile.close();

	return true;
}

void PrintImpList(ImpList& impList)
{
	for(ImpListCIt it=impList.begin(); it != impList.end(); ++it) {
		cout << it->gene1 << "\t" << it->gene2 << "\t" << it->weight << endl;
	}
}

void PrintVariablesMatrix(double ** variablesMatrix, unsigned int dim,
		vector<string>& snpNames)
{
	cout << "Variables matrix:" << endl;
	for(unsigned int i=0; i < snpNames.size(); ++i) {
		if(i ==0) {
			cout << snpNames[i];
		}
		else {
			cout << "\t" << snpNames[i];
		}
	}
	cout << endl;
	for(unsigned int i=0; i < dim; ++i) {
		for(unsigned int j=0; j < dim; ++j) {
			if(j) {
				cout << "\t" << variablesMatrix[i][j];
			}
			else {
				cout << variablesMatrix[i][j];
			}
		}
		cout << endl;
	}
}

bool WriteInteractionsMatrix(string matrixFilename, double** matrix,
		unsigned int dim, vector<string> snpNames)
{
	ofstream outputFileHandle(matrixFilename.c_str());
	if(!outputFileHandle.is_open()) {
		cerr << "ERROR: Could not open matrix file: " << matrixFilename << endl;
		return false;
	}
	for(unsigned int i=0; i < snpNames.size(); ++i) {
		if(i ==0) {
			outputFileHandle << snpNames[i];
		}
		else {
			outputFileHandle << "\t" << snpNames[i];
		}
	}
	outputFileHandle << endl;
	for(unsigned int i=0; i < dim; ++i) {
		for(unsigned int j=0; j < dim; ++j) {
			if(j==0) {
				outputFileHandle << matrix[i][j];
			}
			else {
				outputFileHandle << "\t" << matrix[i][j];
			}
		}
		outputFileHandle << endl;
	}
	outputFileHandle.close();

	return true;
}

bool GrowGeneList(mysqlpp::Connection& conn,	vector<string>& geneNames,
		unsigned int maxSize,	double confidenceThreshold)
{
	// make gene names list into an SQL 'IN' list
	unsigned int nameCounter = 0;
	string inPredicate = "IN (";
	for(vector<string>::const_iterator gnIt = geneNames.begin();
			gnIt != geneNames.end(); ++gnIt) {
		if(nameCounter) {
			inPredicate += ", '" + *gnIt + "'";
		}
		else {
			inPredicate += "'" + *gnIt + "'";
		}
		++nameCounter;
	}
	inPredicate += ")";

	// get the counts for each gene in the gene list, in both directions
	stringstream sql;
	sql << "SELECT gene1 as gene, COUNT(*) AS count FROM interactions WHERE "
			<< "gene1 " + inPredicate + " AND weight > "
			<< confidenceThreshold << " GROUP BY gene1;";
	string sql1 = sql.str();
	sql.str("");
	sql << "SELECT gene2 as gene, COUNT(*) AS count FROM interactions WHERE "
			<< "gene2 " + inPredicate + " AND weight > "
			<< confidenceThreshold << " GROUP BY gene2;";
	string sql2 = sql.str();

	// execute the queries
	vector<counts> counts1;
	vector<counts> counts2;
	try {
#ifdef DEBUG_IMPSNPRANK
		cout << "DEBUG:" << endl << sql1 << endl;
#endif
		mysqlpp::Query query1 = conn.query(sql1);
		query1.storein(counts1);
#ifdef DEBUG_IMPSNPRANK
		cout << "DEBUG:" << endl << sql2 << endl;
#endif
		mysqlpp::Query query2 = conn.query(sql2);
		query2.storein(counts2);
	}
	catch (const mysqlpp::BadQuery& er) {
		// Handle any query errors
		cerr << "Query error: " << er.what() << endl;
		return false;
	}
	catch (const mysqlpp::Exception& er) {
		// Catch-all for any other MySQL++ exceptions
		cerr << "Error: " << er.what() << endl;
		return false;
	}

	// combine the counts from the two gene counts results sets
#ifdef DEBUG_IMPSNPRANK
	cout << endl << "DEBUG: combine the counts" << endl;
#endif
	map<string, unsigned int> combinedCounts;
	unsigned int combinedCountsSum = 0;
	for(CountsListCIt it=counts1.begin(); it != counts1.end(); ++it) {
		combinedCounts[it->gene] += it->count;
		combinedCountsSum += it->count;
#ifdef DEBUG_IMPSNPRANK
		cout << "DEBUG: gene1 query record counts for " << it->gene
				<< "\t" << it->count << endl;
#endif
	}
	for(CountsListCIt it=counts2.begin(); it != counts2.end(); ++it) {
		combinedCounts[it->gene] += it->count;
		combinedCountsSum += it->count;
#ifdef DEBUG_IMPSNPRANK
		cout << "DEBUG: gene2 query record counts for " << it->gene
				<< "\t" << it->count << endl;
#endif
	}
#ifdef DEBUG_IMPSNPRANK
	cout << "DEBUG: combined counts sum: " << combinedCountsSum << endl;
#endif

	// determine each gene's proportions in the IMP database
#ifdef DEBUG_IMPSNPRANK
	cout << endl << "DEBUG: combined count sums (degree)" << endl;
#endif
	map<string, double> geneProportions;
	map<string, unsigned int>::const_iterator sumsIt = combinedCounts.begin();
	for(; sumsIt != combinedCounts.end(); ++sumsIt) {
		geneProportions[sumsIt->first] =
				((double) sumsIt->second) / ((double) combinedCountsSum);
#ifdef DEBUG_IMPSNPRANK
		cout << "DEBUG: " << sumsIt->first << "\t" << sumsIt->second << endl;
#endif
	}

#ifdef DEBUG_IMPSNPRANK
	cout << endl << "DEBUG: SQL SELECT LIMIT values for each gene" << endl;
	cout << "gene\tproportion\tnodes to add" << endl;
#endif
	double pSum = 0.0;
	map<string, int> geneLimits;
	map<string, double>::const_iterator pIt = geneProportions.begin();
	int numAdded = 0;
	for(; pIt != geneProportions.end(); ++pIt) {
		geneLimits[pIt->first] =
				static_cast<int>(round(pIt->second * maxSize)-1);
#ifdef DEBUG_IMPSNPRANK
		cout << pIt->first << "\t" << pIt->second
				<< "\t" << geneLimits[pIt->first] << endl;
#endif
		pSum += pIt->second;
		numAdded += geneLimits[pIt->first];
	}
#ifdef DEBUG_IMPSNPRANK
	cout << "\tpSum = " << pSum << "\tadding " << numAdded << " nodes" << endl;
#endif

	// use the proportions to select the top connected genes to each gene
	// and grow the gene list to the target size
	cout << endl << "DEBUG: adding new genes/nodes" << endl;
	map<string, int>::const_iterator glIt = geneLimits.begin();
	for(; glIt != geneLimits.end(); ++glIt) {
		string thisGene = glIt->first;
		int limit = glIt->second;
		if(limit) {
#ifdef DEBUG_IMPSNPRANK
			cout << "DEBUG: Getting " << limit << " new genes" << endl;
#endif
			vector<string> newGenes = FindTopConnectedGenes(conn, thisGene,
					limit, confidenceThreshold);
#ifdef DEBUG_IMPSNPRANK
			cout << "DEBUG: " << newGenes.size() << " found" << endl;
#endif
			vector<string>::const_iterator ngIt = newGenes.begin();
			for(; ngIt != newGenes.end(); ++ngIt) {
#ifdef DEBUG_IMPSNPRANK
				cout << "DEBUG: Adding gene for " << thisGene << ": " << *ngIt << endl;
#endif
				geneNames.push_back(*ngIt);
			}
		}
		else {
#ifdef DEBUG_IMPSNPRANK
			cout << "DEBUG: no new genes required, skipping " << thisGene << endl;
#endif
		}
	}

	return true;
}

vector<string> FindTopConnectedGenes(mysqlpp::Connection& conn,
		string thisGene, int topN, double confidenceThreshold)
{
	vector<string> topGenes;

	// get the interactions for the gene, in both directions
	stringstream sql;
	sql << "SELECT * FROM interactions WHERE "
			<< "(gene1='" << thisGene << "' OR gene2='" << thisGene << "')"
			<< " AND weight > "	<< confidenceThreshold << " ORDER BY weight DESC;";
	string sql1 = sql.str();

	// execute the query
	vector<interactions> results;
	try {
#ifdef DEBUG_IMPSNPRANK
		cout << "DEBUG:" << sql1 << endl;
#endif
		mysqlpp::Query query1 = conn.query(sql1);
		query1.storein(results);
	}
	catch (const mysqlpp::BadQuery& er) {
		// Handle any query errors
		cerr << "Query error: " << er.what() << endl;
		return topGenes;
	}
	catch (const mysqlpp::Exception& er) {
		// Catch-all for any other MySQL++ exceptions
		cerr << "Error: " << er.what() << endl;
		return topGenes;
	}

#ifdef DEBUG_IMPSNPRANK
	cout << "DEBUG: query returned " << results.size() << " records" << endl;
#endif
	vector<interactions>::const_iterator resIt = results.begin();
	for(; (resIt != results.end()) && (topGenes.size() < (unsigned int) topN);
			++resIt) {
		string gene1 = resIt->gene1;
		string gene2 = resIt->gene2;
		// cout << "DEBUG: " << gene1 << "\t" << gene2 << endl;
		if(gene1 == thisGene) {
			if(find(topGenes.begin(), topGenes.end(), gene2) == topGenes.end()) {
				topGenes.push_back(gene2);
#ifdef DEBUG_IMPSNPRANK
				cout << "DEBUG: Adding to topGenes: " << gene2 << endl;
#endif
			}
		}
		if(gene2 == thisGene) {
			if(find(topGenes.begin(), topGenes.end(), gene1) == topGenes.end()) {
				topGenes.push_back(gene1);
#ifdef DEBUG_IMPSNPRANK
				cout << "DEBUG: Adding to topGenes: " << gene1 << endl;
#endif
			}
		}
	}

	return topGenes;
}
