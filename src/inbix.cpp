// inbix.cpp - Insilico Bioinformatics - Bill White
//
// (c) 2014 McKinney Insilico Bioinformatics Lab
// The University of Tulsa
//
// Code borrowed heavily from the PLINK project
// PLINK (c) 2005-2009 Shaun Purcell

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <cassert>

#include <armadillo>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "crandom.h"
#include "perm.h"
#include "sets.h"
#include "linear.h"
#include "logistic.h"
#include "phase.h"
#include "clumpld.h"
#include "nlist.h"
#include "sets.h"
#include "stats.h"
#include "idhelp.h"
#include "zed.h"

#include <omp.h>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>

#include "regain.h"
#include "InteractionNetwork.h"
#include "CentralityRanker.h"
#include "ArmadilloFuncs.h"
#include "EpistasisEQtl.h"

// ReliefSeq project integration - bcw - 8/7/16
#include "Dataset.h"
#include "PlinkInternalsDataset.h"
#include "AttributeRanker.h"
#include "ReliefSeqController.h"
#include "ReliefF.h"
#include "RReliefF.h"

// Ranger random forest project integration - bcw - 9/26/16
#include "Data.h"
#include "DataDouble.h"
#include "Forest.h"
#include "ForestClassification.h"
#include "ForestRegression.h"

// Evaporative Cooling (EC) - bcw - 9/29/16
#include "EvaporativeCooling.h"
#include "Deseq.h"
#include "Edger.h"

using namespace std;
using namespace arma;
using namespace boost;
namespace po = boost::program_options;

ofstream LOG;
string PVERSION;
string PDATE;
string PREL;
Plink * PP;
map<string, int> Range::groupNames;

bool pvalComparator(const matrixElement &l, const matrixElement &r) {
  return l.first < r.first;
}

int main(int argc, char* argv[]) {
	/////////////////////////
	// Setup, display title

	cout.setf(ios::fixed);
	cout.precision(8);

	set_new_handler(NoMem);

	PVERSION = "0.04"; // 4 chars
	PREL = " "; // space or p (full, or prelease) 
	PDATE = "2016       "; // 11 chars

	//////////////////
	// The major class
	Plink P;
	PP = &P;

	/////////////////////////////////////////////////////
	// General class for all haplotype-related functions
	P.haplo = new HaploPhase(P);

	//////////////////////////
	// Command line arguments
	CArgs a(argc, argv);
	getOutputFilename(a);

	//////////////////////////
	// Start logging, title
	LOG.open(string(par::output_file_name + ".log").c_str());

#if !defined(IQTL) && !defined(DCVAR)
	P.printLOG("\n"
					"@----------------------------------------------------------@\n"
					"|        inbix        |     v" + PVERSION + PREL + "     |   " + PDATE + "     |\n"
					"|----------------------------------------------------------|\n"
					"|  (C) 2016 Bill White, GNU General Public License, v2     |\n"
					"@----------------------------------------------------------@\n"
					"\n");
#endif
  
	//////////////////////////
	// Fully parse command line
	setOptions(a);

	/////////////////////
	// Permutation class
	if(par::random_seed == 0)
		CRandom::srand(time(0));
	else
		CRandom::srand(par::random_seed);

	Perm perm(P);
	P.pperm = &perm;

	/////////////
	// Time stamp
	P.printLOG("Writing console output to log file [ " +
					par::output_file_name + ".log ]\n");

	time_t curr = time(0);
	string tdstamp = (string) ctime(&curr);
	P.printLOG("Analysis started: " + tdstamp + "\n");

	/////////////////////////////////////
	// Validate and record all arguments
	a.check_unused_options(P);

	if(par::output_file_name.find(".", 0) != string::npos)
		P.printLOG("** For gPLINK compatibility, do not use '.' in --out **\n");

	//////////////////////////
	// Some basic definitions
	if(par::species_dog) defineDogChromosomes();
	else if(par::species_sheep) defineSheepChromosomes();
	else if(par::species_cow) defineCowChromosomes();
	else if(par::species_horse) defineHorseChromosomes();
	else if(par::species_rice) defineRiceChromosomes();
	else if(par::species_mouse) defineMouseChromosomes();
	else defineHumanChromosomes();

	///////////////////////////////
	// Web-based SNPServer lookup?
	if(par::lookup) {
		P.lookup();
		shutdown();
	}

	if(par::lookup2) {
		P.lookup2();
		shutdown();
	}

	/////////////////////////
	// ID helper?
	if(par::idhelp) {
		IDHelper ID;
		ID.idHelp();
		shutdown();
	}

	/////////////////////////
	// File compression utility
	if(par::compress_file) {
		fileCompress();
		shutdown();
	}

	if(par::uncompress_file) {
		fileUncompress();
		shutdown();
	}

	//////////////////////////////////////////////////
	// Main Input files

	// Simulate or read in data:
	// Binary or ASCII format; transposed/long/generic
	if(par::dummy) P.dummyLoader();
	else if(par::greport) P.displayGeneReport();
	else if(par::annot_file) P.annotateFile();
	else if(par::meta_analysis) P.metaAnalysis();
	else if(par::rare_test_score_range) P.displayRareRange();
	else if(par::simul) {
		if(par::simul_qt)
			P.simulateSNPs_QT();
		else
			P.simulateSNPs();
	} else if(par::cnv_list) P.setUpForCNVList();
	else if(par::read_bitfile) P.readBinData();
	else if(par::lfile_input) P.readDataLongFormat();
	else if(par::tfile_input) P.readTransposedData();
	else if(par::read_ped) P.readData();
	else if(par::gvar) {
		par::load_gvar = true;
		P.readGenericVariantData();
	} else if(par::dosage_assoc) {
		P.readFamFile(par::famfile);
		if(par::dosage_hasMap) {
			checkFileExists(par::mapfile);
			vector<bool> include;
			vector<int> include_pos(0);
			int nvar = 0;
			P.readMapFile(par::mapfile,
							include,
							include_pos,
							nvar);
		}
	}

	////////////////////////////////////////////
	// network deconvolution - bcw - 10/22/13
  if(par::do_deconvolution) {
    P.printLOG("Performing network deconvolution\n");
    string matrixFile = "";
    MatrixFileType fileType = INVALID_FILE;
    bool isUpperTriangular = false;
    if(par::do_regain_post) {
      matrixFile = par::regainFile;
      fileType = REGAIN_FILE;
      P.printLOG("Reading correlation from reGAIN format\n");
    }
    if(par::sifNetwork) {
      matrixFile = par::sifFile;
      fileType = SIF_FILE;
      P.printLOG("Reading correlation SIF format\n");
    }
    if(fileType == INVALID_FILE) {
      error("Error running network deconvolution: no valid network file type");
    }
    P.printLOG("Creating network\n");
    InteractionNetwork* network = 
      new InteractionNetwork(matrixFile, fileType, isUpperTriangular, &P);
    mat nd;
    P.printLOG("Running deconvolve\n");
    if(!network->Deconvolve(nd, par::deconvolutionAlpha, 
      par::deconvolutionBeta, par::deconvolutionControl)) {
      error("Deconvolution failed");
    }
    // write resulting matrix
    string outMatrixFilename = par::output_file_name + ".deconvolved";
		armaWriteMatrix(nd, outMatrixFilename, network->GetNodeNames());
    shutdown();
  }
		
	////////////////////////////////////////////
	// A SIF file specified? - bcw - 5/23/13
	if(par::sifNetwork && par::sifToGain) {
		P.printLOG("Converting SIF network to reGAIN\n");
		P.outputSifToGain(par::sifFile);
		shutdown();
	}

	////////////////////////////////////////////
	// A numeric file specified? - bcw - 4/20/13
	if(par::numeric_file) {
		// check for SNPs already loaded
		par::have_snps = par::read_bitfile || par::read_ped ||
						par::tfile_input || par::lfile_input;

		// compute a variance-covariance matrix of the numeric data
		if(par::do_covariance_matrix) {
			if(par::have_snps) {
				error("Covariance matrix is only supported in numeric data");
			}
//			matrix_t testX;
//			matrix_t testCovMatrix;
//			matrix_t testCorMatrix;
//			vector<string> variableNames;
//			matrixRead(par::numeric_filename, testX, variableNames);
//			matrixComputeCovariance(testX, testCovMatrix, testCorMatrix);
//			matrixWrite(testCorMatrix, "foobar.txt", variableNames);
			
			// read the data matrix
			mat X;
			vector<string> variableNames;
			if(!armaReadMatrix(par::numeric_filename, X, variableNames)) {
				error("Cannot read matrix file: " + par::numeric_filename);
			}
			// compute covariances/correlations
			mat covMatrix;
			mat corMatrix;
			if(armaComputeCovariance(X, covMatrix, corMatrix)) {
				// write results
				string covFilename = par::output_file_name + ".covariance";
				string corFilename = par::output_file_name + ".correlation";
				armaWriteMatrix(covMatrix, covFilename, variableNames);
				armaWriteMatrix(corMatrix, corFilename, variableNames);
			} else {
				error("Could not compute covariance/correlation matrices");
			}
			shutdown();
		}

		// prepare numeric attributes for further processing
		if(!P.readNumericFile()) {
			error("Problem reading the numeric file: [" + par::numeric_filename + "]");
		}
		par::have_numerics = true;

		if(par::do_numeric_summary) {
			P.printLOG("Reporting numeric file summary statistics.\n");
			reportNumericSummaryStats();
			shutdown();
		}
    
    // data set transforms prior to analysis - bcw - 10/30/13
    if(par::do_numeric_standardize) {
      P.printLOG("Standardizing numeric variables.\n");
      if(!numericMeanCenter()) {
        error("Mean centering numerics failed.");
      }
      if(!numericStandardize()) {
        error("Standardizing numerics failed.");
      }
      if(par::exportDelimited) {
        P.printLOG("Performing data set export to delimited format\n");
        // switch from SNP-major to individual-major data orientation!
        P.SNP2Ind();
        string fileExtension = ".delim";
        if(par::exportDelimiter == "\t") {
          fileExtension = ".tab";
        }
        if(par::exportDelimiter == ",") {
          fileExtension = ".csv";
        }
        string delimitedFilename = par::output_file_name + fileExtension;
        P.outputDelimitedFile(delimitedFilename, par::exportDelimiter);
        shutdown();
      }
    }

		// extract attributes listed in user file to new file
		if(par::do_numeric_extract) {
			P.printLOG("Extracting numeric attributes to a new numeric file.\n");
			P.outputNumericExtract(par::numeric_extract_file);
			shutdown();
		}
    
    // low value filter
    if(par::do_numeric_lowval_filter) {
      P.printLOG("Numeric filter mode\n");
      boolvec_t varIndicesThatPass(P.nlistname.size(), true);
      P.printLOG("Running numeric low value filter\n");
      numericLowValueFilter(par::numeric_lowval_percentile, varIndicesThatPass);
      string filteredFilename = par::output_file_name + ".lowvalfilter.num";
      P.printLOG("Writing filtered numeric file [ " + filteredFilename + " ]\n");
      P.outputNumericFiltered(filteredFilename, varIndicesThatPass);
      shutdown();
    }
    
    // low variance filter
    if(par::do_numeric_lowvar_filter) {
      P.printLOG("Running numeric low variance filter\n");
      boolvec_t varIndicesThatPass(P.nlistname.size(), true);
      numericVarianceFilter(par::numeric_lowvar_percentile, varIndicesThatPass);
      string filteredFilename = par::output_file_name + ".lowvarfilter.num";
      P.printLOG("Writing filtered numeric file [ " + filteredFilename + " ]\n");
      P.outputNumericFiltered(filteredFilename, varIndicesThatPass);
      shutdown();
    }

	}

	// Set number of individuals
	P.n = P.sample.size();

	// Set number of pairs
	P.np = (int) ((double) (P.n * (P.n - 1)) / (double) 2);

	// Total number of all (test+background) loci
	P.nl_all = P.locus.size();

	// Number of permutations to store
	P.maxr2.resize(par::replicates);

	// Check for duplicate individual or SNP names
	checkDupes(P);

	/////////////////////////////////////
	// Merge with a secondary data file 
	// Standard (non-list) mode
	if(par::merge_data && !par::merge_list) {

		if(par::merge_binary)
			P.mergeBinaryData();
		else
			P.mergeData();

		// Reset number of individuals
		P.n = P.sample.size();

		// Set number of pairs
		P.np = (int) ((double) (P.n * (P.n - 1)) / (double) 2);

		// Total number of all (test+background) loci
		P.nl_all = P.locus.size();
	}

	/////////////////////////////////////
	// Merge with a secondary data file 
	// List mode
	if(par::merge_list)
		P.mergeList();

	//////////////////////////////////////////  
	// A different phenotype file specified?
	if(par::pheno_file) P.readPhenoFile();
	else if(par::make_pheno) P.makePhenotype();
	else if(par::multiple_phenotypes) P.readMultiplePhenoFile();

	////////////////////////////////
	// Remove any individuals with 
	// missing phenotypes

	if(!par::ignore_phenotypes)
		removeMissingPhenotypes(P);

	//////////////////////////////////
	// Binary affection status coding
	if(par::bt)
		affCoding(P);

	/////////////////////////////////
	// Update MAP file information?
	if(par::update_map)
		P.updateMapFile();

	/////////////////////////////////
	// Update FAM information?
	if(par::update_ids || par::update_parents || par::update_sex || par::update_pheno)
		P.updateFamFile();

	/////////////////////////////////
	// Update allele file information?
	if(par::update_alleles)
		P.updateAlleles();

	/////////////////////////////////
	// Flip DNA strand for any SNPs? 
	if(par::flip_strand)
		P.flipStrand();

	/////////////////////////////////
	// Recode any alleles? 
	if(par::recode_ACGT || par::recode_1234)
		P.alleleRecoding();

	//////////////////////////////////////////////////////////
	// Output a specific set of SNPs (--extract or --exclude)
	if(par::extract_before_exclude) {
		if(par::extract_set)
			P.extractExcludeSet(false);
		if(par::exclude_set)
			P.extractExcludeSet(true);
	} else {
		if(par::exclude_set)
			P.extractExcludeSet(true);
		if(par::extract_set)
			P.extractExcludeSet(false);
	}

	/////////////////////////////////////////////////////////////
	// Output a specific set of individuals --remove or --keep
	if(par::remove_before_keep) {
		if(par::remove_indiv)
			P.removeIndividuals(false);
		if(par::keep_indiv)
			P.removeIndividuals(true);
	} else {
		if(par::keep_indiv)
			P.removeIndividuals(true);
		if(par::remove_indiv)
			P.removeIndividuals(false);
	}

	///////////////////////////////////////////////
	// Filter based on attribute files
	if(par::snp_attrib_filter)
		P.attribFilterSNP();

	if(par::ind_attrib_filter)
		P.attribFilterInd();

	///////////////////////////////////////////////
	// Filter based on qualiy scores
	if(par::read_snp_qual)
		P.filterQualSNPs();

	if(par::read_geno_qual)
		P.filterQualGenotypes();

	//////////////////////////////////////////////////
	// Pull a random subset of SNPs?
	if(par::thin_snps)
		P.thinSNPs();

	/////////////////////////////////////////////////////////////
	// If in --genome list mode, keep the two lists of individuals
	if(par::genome_2sets) {
		P.keep2SetsForGenome();
	}

	///////////////////////////////////////////////
	// Read a list of obligatory missing genotypes?
	if(par::oblig_missing)
		P.setObligMissing();

	//////////////////////////////////////////////////
	// Filter individuals based on external covariate? 
	if(par::filter_on_covar) {
		P.filterOnCovariate();
		// Reset number of individuals
		P.n = P.sample.size();
		P.np = (int) ((double) (P.n * (P.n - 1)) / (double) 2);
	}

	////////////////////////////
	// Any simple preset filters
	if(par::filter_males)
		P.filterOnMale();
	else if(par::filter_females)
		P.filterOnFemale();

	if(par::filter_cases)
		P.filterOnCase();
	else if(par::filter_controls)
		P.filterOnControl();

	if(par::filter_founders)
		P.filterOnFounder();
	else if(par::filter_nonfounders)
		P.filterOnNonFounder();

	//////////////////////////////// 
	// A covariate file specified?
	if(par::covar_file) {
		// Multiple covariates?
		if(par::clist) {
			if(!P.readCovListFile())
				error("Problem reading the covariates");
		} else // a single covariate
		{
			if(!P.readCovariateFile())
				error("Problem reading the specified covariate from the covariate file");
		}
	}

	//////////////////////////////////////
	// Assign cluster solution from file 
	if(par::include_cluster_from_file) {
		P.printLOG("Reading clusters from [ " +
						par::include_cluster_filename + " ]\n");
		if(!P.readClusterFile())
			error("Problem reading from [ " + par::include_cluster_filename + " ]");
	} else if(par::sol_family) {
		P.printLOG("Setting clusters based on family IDs\n");
		vector<string> famlist;
		P.kname.resize(0);
		for(int i = 0; i < P.n; i++) {
			Individual * person = P.sample[i];
			bool match = false;
			for(unsigned int j = 0; j < famlist.size(); j++)
				if(person->fid == famlist[j]) {
					match = true;
					person->sol = j;
				}
			if(!match) {
				famlist.push_back(person->fid);
				person->sol = famlist.size() - 1;
				P.kname.push_back(person->fid);
			}
		}

		// Set number of clusters/families
		P.nk = famlist.size();

		// Set klist variable
		P.klist.clear();
		for(int j = 0; j < P.nk; j++)
			P.klist.push_back(new Cluster);

		for(int i = 0; i < P.n; i++)
			if(P.sample[i]->sol > -1)
				P.klist[P.sample[i]->sol]->person.push_back(P.sample[i]);

	} else {
		P.klist.clear();
		P.klist.push_back(new Cluster);
		for(int i = 0; i < P.n; i++)
			P.klist[0]->person.push_back(P.sample[i]);
	}

	/////////////////////////////////////////
	// Zero-out specific sets of genotypes?
	if(par::zero_cluster)
		P.zeroOnCluster();

	/////////////////////////////////
	// Fix reference allele? 

	if(par::set_reference_allele)
		P.setReferenceAllele();

	//////////////////////////////////
	// Determine formats for printing
	P.prettyPrintLengths();

	//////////////////////////////////////////////////
	//                                              //
	// Process a dosage file                        //
	//                                              //
	//////////////////////////////////////////////////
	if(par::dosage_assoc) {
		// Normal behavior is to load data, and perform 
		// analysis; if the hard-call option is specified, 
		// then this will generate a dataset, that we can
		// subsequent filter and save, etc, as usual, i.e.
		// in that case, do not halt
		P.processDosageFile();

		if(!par::dosage_hard_call)
			shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// Handle CNV segments separately               //
	//                                              //
	//////////////////////////////////////////////////
	if(par::cnv_list) {
		P.readCNVList();
		P.processCNVList();
		shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// Handle non-SNP data separately               //
	//                                              //
	//////////////////////////////////////////////////
	if(par::gvar || par::gvar_write) {

		// We might want to load generic variants on top
		// of existing SNP data; or afresh if none of the 
		// above have been specified

		if(!par::load_gvar)
			P.readGenericVariantData();

		if(par::gvar_write) {
			P.outputGenericVariantFile();
			shutdown();
		}

		P.processGVAR();
		shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// Misc. .genome grouper utility                //
	//                                              //
	//////////////////////////////////////////////////
	if(par::genome_groups) {
		P.groupGenome();
		shutdown();
	}

	//////////////////////////////////
	// Missing code

	//   if ( par::bt && ! par::missing_phenotype_explicit ) 
	//     par::missing_phenotype = "0";

	//////////////////////////////////////////////////
	//                                              //
	// Basic MAF, genotyping filters & HWE/ME       //
	//                                              //
	//////////////////////////////////////////////////
	P.printLOG("Before frequency and genotyping pruning, there are "
					+ int2str(P.nl_all) + " SNPs\n");

	if(!par::FIXED_p) {
		P.filterSNPs();
	} else
		for(int i = 0; i < P.nl_all; i++) {
			if(P.locus[i]->allele1 == "1")
				P.locus[i]->freq = par::FIX_p;
			else
				P.locus[i]->freq = 1 - par::FIX_p;
		}

	P.printLOG("After frequency and genotyping pruning, there are "
					+ int2str(P.nl_all) + " SNPs\n");

	// NOTE:: some methods do not require SNP and/or numeric data to be loaded,
	// so do not exit with an error when these conditions are detected.
	if((P.nl_all == 0) && (P.nlistname.size() == 0) &&
					(!par::do_modularity) && (!par::do_ranking) &&
					(!par::do_regain_post) && (!par::sifNetwork)) {
		error("Stopping as there are no SNPs or numeric attributes "
						"left for analysis\n");
	}

	if((P.n == 0) && (!par::do_modularity) && (!par::do_ranking)
					&& (!par::do_regain_post) && (!par::sifNetwork)) {
		error("Stopping as there are no individuals left for analysis\n");
	}

	//////////////////////////////////////////////////
	// Re-report basic case/control counts, etc
	summaryBasics(P);

	//////////////////////////////////////////////////
	// Any null allele codes (monomorhpics)?
	for(int l = 0; l < P.nl_all; l++) {
		if(P.locus[l]->allele1 == "")
			P.locus[l]->allele1 = "0";
	}

	/////////////////////////////////////////////////////////////////////////////
	// perform epistatic eQTL analysis
  if(par::do_iqtl) {
		P.printLOG("Performing iQTL analysis\n");

    if(par::iqtl_expression_file == "") {
      error("Transcript expression file is required. Use --transcript-matrix");
    }
    if(par::iqtl_coord_file == "") {
      error("Transcript coordinate file is required. Use --coordinates");
    }
    
    // read the expression data as a numeric file in PLINK format
		P.printLOG("Reading transcripts from [" + par::iqtl_expression_file + "]\n");
    par::numeric_filename = par::iqtl_expression_file;
		if(!P.readNumericFile()) {
			error("Cannot read eQTL expression file: " + par::iqtl_expression_file);
		}
    
    EpistasisEQtl* iqtl = new EpistasisEQtl();
    
    // read transcript coordinate information from file
		P.printLOG("Reading transcript coordinates from [" + par::iqtl_coord_file + "]\n");
    if(!iqtl->ReadTranscriptCoordinates(par::iqtl_coord_file)) {
      error("Cannot read coordinates file: " + par::iqtl_coord_file);
    }

    // set parameters
    iqtl->SetLocalCis(par::iqtl_local_cis);
    iqtl->SetRadius(par::iqtl_radius);
    // added 4/21/15
	  if(par::do_iqtl_tf) {
	    iqtl->SetTF(par::do_iqtl_tf);
	    iqtl->SetTFRadius(par::iqtl_tf_radius);
	    if(par::iqtl_tf_coord_file != "") {
		    // read transcription factor coordinate information from file
				P.printLOG("Reading TF coordinates from [" + par::iqtl_tf_coord_file + "]\n");
		    if(!iqtl->ReadTranscriptFactorCoordinates(par::iqtl_tf_coord_file)) {
		      error("Cannot read TF coordinates file: " + par::iqtl_tf_coord_file);
		    }
	    }
    }

    // run the analysis
  	P.SNP2Ind();
    if(!iqtl->Run(true)) {
      error("iQTL analysis failed");
    }

    // clean up and exit gracefully
    delete iqtl;
    shutdown();
  }
  
	/////////////////////////////////////////////////////////////////////////////
	// compute a coexpression matrix of the numeric data
	if(par::do_coexpression_all || par::do_coexpression_casecontrol) {
		if(!par::numeric_file) {
			error("Co-expression matrix requires numeric data");
		}
		if(par::have_snps) {
			error("Co-expression matrix is only supported in numeric data");
		}
		if(!par::bt) {
			error("Co-expression matrix is only supported in case-control data");
		}
		if(par::do_coexpression_all) {
			P.printLOG("Computing coexpression for ALL variables.\n");
			mat X;
			if(!armaGetPlinkNumericToMatrixAll(X)) {
				error("Cannot read numeric data into matrix");
			}
			// compute covariances/correlations
			mat covMatrix;
			mat corMatrix;
			if(armaComputeCovariance(X, covMatrix, corMatrix)) {
				// write results
				string coexpFilename = par::output_file_name + ".coexpression";
				armaWriteMatrix(corMatrix, coexpFilename, P.nlistname);
			} else {
				error("Could not compute coexpression matrix");
			}
		} else {
			P.printLOG("Computing coexpression for CASES and CONTROLS.\n");
			mat X;
			mat Y;
			if(!armaGetPlinkNumericToMatrixCaseControl(X, Y)) {
				error("Cannot read numeric data into case-control matrices");
			}
			// compute covariances/correlations
			mat covMatrixX;
			mat corMatrixX;
			if(armaComputeCovariance(X, covMatrixX, corMatrixX)) {
				// write results
				string coexpFilename = par::output_file_name + ".coexpression.cases";
				armaWriteMatrix(corMatrixX, coexpFilename, P.nlistname);
			} else {
				error("Could not compute coexpression matrix for cases");
			}
			mat covMatrixY;
			mat corMatrixY;
			if(armaComputeCovariance(Y, covMatrixY, corMatrixY)) {
				// write results
				string coexpFilename = par::output_file_name + ".coexpression.controls";
				armaWriteMatrix(corMatrixY, coexpFilename, P.nlistname);
			} else {
				error("Could not compute coexpression matrix for controls");
			}
      // write difference matrix - bcw - 10/18/13
      string coexpFilename = par::output_file_name + ".coexpression.ccdiff";
      mat diffMatrix = corMatrixX - corMatrixY; 
      armaWriteMatrix(diffMatrix, coexpFilename, P.nlistname);
		}
		shutdown();
	}

	/////////////////////////////////////////////////////////////////////////////
	// data set export requested - bcw - 5/21/13
	if(par::exportArff) {
		P.printLOG("Performing data set export to Weka ARFF format\n");
		// switch from SNP-major to individual-major data orientation!
		P.SNP2Ind();
		string arffFilename = par::output_file_name + ".arff";
		P.outputArffFile(arffFilename);
		shutdown();
	}

 	/////////////////////////////////////////////////////////////////////////////
	// delimited format for Excel, R, etc - bcw - 5/22/13
	if(par::exportDelimited) {
		P.printLOG("Performing data set export to delimited format\n");
		// switch from SNP-major to individual-major data orientation!
		P.SNP2Ind();
		string fileExtension = ".delim";
		if(par::exportDelimiter == "\t") {
			fileExtension = ".tab";
		}
		if(par::exportDelimiter == ",") {
			fileExtension = ".csv";
		}
		string delimitedFilename = par::output_file_name + fileExtension;
		P.outputDelimitedFile(delimitedFilename, par::exportDelimiter);
		shutdown();
	}

	/////////////////////////////////////////////////////////////////////////////
	// variable ranking requested - bcw - 5/16/13
	if(par::do_ranking) {
		P.printLOG("Performing variable ranking\n");
		P.SNP2Ind();
		if(par::ranker_method == "centrality" ||
						par::ranker_method == "centrality_power" ||
						par::ranker_method == "centrality_gauss") {
			if(!par::do_regain_post) {
				error("Centrality ranking requires a reGAIN file");
			}
			P.printLOG("Ranking by network centrality: " + par::ranker_method + "\n");
			CentralityRanker cr(par::regainFile);
			if(par::ranker_centrality_gamma > 0) {
				P.printLOG("Network centrality gamma set to: " +
								dbl2str(par::ranker_centrality_gamma) + "\n");
				cr.SetGlobalGamma(par::ranker_centrality_gamma);
			}
			if(par::ranker_method == "centrality_power") {
				P.printLOG("Centrality using fixed gamma 0.85, power method\n");
				cr.SetGlobalGamma(0.85);
				if(!cr.CalculateCentrality(POWER_METHOD)) {
					error("Centrality ranking failed");
				}
			}
			if(par::ranker_method == "centrality" ||
							par::ranker_method == "centrality_gauss") {
				P.printLOG("Centrality using adaptive gamma\n");
				if(!cr.CalculateCentrality(GAUSS_ELIMINATION)) {
					error("Centrality ranking failed");
				}
			}
			if(par::verbose) {
				cr.WriteToConsole(par::ranker_top_n);
			}
			string saveFilename = par::output_file_name + ".ranks";
			cr.WriteToFile(saveFilename, par::ranker_top_n);
		}
		if(par::ranker_method == "regressions" ||
						par::ranker_method == "regressionb" ||
						par::ranker_method == "regressionp") {
			if(par::bt) {
				P.printLOG("Ranking by logistic regression\n");
			} else {
				P.printLOG("Ranking by linear regression\n");
			}
			rankedlist_t ranks;
			RegressionRankResults results;
			if(par::ranker_method == "regressions") {
				P.printLOG("Using regression standardized coefficient values\n");
				rankByRegression(REGRESSION_RANK_STAT, ranks, results);
			}
			if(par::ranker_method == "regressionb") {
				P.printLOG("Using regression beta coefficient values\n");
				rankByRegression(REGRESSION_RANK_BETA, ranks, results);
			}
			if(par::ranker_method == "regressionp") {
				P.printLOG("Using regression p-values\n");
				rankByRegression(REGRESSION_RANK_PVAL, ranks, results);
			}
			string outFile = par::output_file_name + ".ranks";
			P.printLOG("Writing scores to [" + outFile + "]\n");
			ofstream outputFileHandle(outFile.c_str());
			int numToWrite = ranks.size();
			int topN = par::ranker_top_n;
			if((topN > 0) && (topN <= ranks.size())) {
				numToWrite = topN;
			} else {
				if(topN != -1) {
					cout << "WARNING: Attempting to use top N outside valid range: "
									<< topN << ". Using all ranks." << endl;
				}
			}
			for(int i = 0; i < numToWrite; i++) {
				double score = ranks[i].first;
				string name = ranks[i].second;
				outputFileHandle << name << "\t" << score << endl;
			}
			outputFileHandle.close();

			string outFileDetail = par::output_file_name + ".ranks.detail";
			P.printLOG("Writing detailed results to [" + outFileDetail + "]\n");
			ofstream outputFileDetailHandle(outFileDetail.c_str());
			outputFileDetailHandle << "var\tcoef\tpval\tstat"	<< endl;
			for(int i = 0; i < results.vars.size(); i++) {
				outputFileDetailHandle 
								<< results.vars[i] << "\t"
								<< results.coefs[i] << "\t"
								<< results.pvals[i] << "\t"
								<< results.stats[i]
								<< endl;								
			}
			outputFileDetailHandle.close();
		}
		if(par::ranker_method == "random") {
			P.printLOG("Ranking Random\n");
			P.printLOG("***** Random ranking not implemented *****\n");
		}
		shutdown();
	}

	// --------------------------------------------------------------------------
	// permuted GAIN - bcw - 5/23/14
	if(par::do_ranker_permutation) {
	
		P.printLOG("Performing GAIN permutation analysis\n");
		P.SNP2Ind();

		int M = P.nlistname.size();
		int N = P.sample.size();
		mat permResults;
		permResults.resize(par::rankerPermNum, M);
		int numPerms = par::rankerPermNum;
		for(int perm = 0; perm < numPerms; ++perm)	{

			// permute phenotype labels - careful!
			int n1 = 0, n2 = 0;
			vector<int> pIdx(N);
			permute(pIdx); // this might not do what you expect
			vector<int> newPhenos;
			for(int i=0; i < N; i++) {
				newPhenos.push_back(P.sample[pIdx[i]]->phenotype);
			}
			for(int i=0; i < N; i++) {
				if(newPhenos[i] == 1) {
					P.sample[i]->pperson->aff = false;
					P.sample[i]->aff = false;
					P.sample[i]->phenotype = 1;
					++n1;
				} else {
					P.sample[i]->pperson->aff = true;
					P.sample[i]->aff = true;
					P.sample[i]->phenotype = 2;
					++n2;
				}
			}

			CentralityRanker* cr = 0;
			if(par::rankerPermMethod == "regain") {
				// run reGAIN 
				Regain* regain = new Regain(
								par::regainCompress,
								par::regainSifThreshold,
								par::have_numerics,
								par::regainComponents,
								par::regainFdrPrune,
								true);
				// reGAIN transform options
				if(par::regainMatrixTransform == "none") {
					regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_NONE);
				} else {
					if(par::regainMatrixTransform == "threshold") {
						regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_THRESH);
					} else {
						if(par::regainMatrixTransform == "abs") {
							regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_ABS);
						} else {
							error("reGAIN output transform allowed options: {none, threshold, abs}");
						}
					}
				}
				regain->performPureInteraction(false);
				regain->run();
				// SNPrank
				cr = new CentralityRanker(regain->getRawMatrix(), M, P.nlistname);

				delete regain;
			}
			
			if(par::rankerPermMethod == "dcgain") {
				mat dcgain(M, M);
		    // t-test for diagonal
		    double df = N - 2;
		    for(int i=0; i < M; ++i) {
		      double t;
		      tTest(i, t);
		      dcgain(i, i) = t;
		    }

		    // z-test for off-diagonal elements
		    mat X;
		    mat Y;
		    if(!armaGetPlinkNumericToMatrixCaseControl(X, Y)) {
		      error("Cannot read numeric data into case-control matrices");
		    }

		    // compute covariances/correlations
		    mat covMatrixX;
		    mat corMatrixX;
		    if(!armaComputeCovariance(X, covMatrixX, corMatrixX)) {
		      error("Could not compute coexpression matrix for cases");
		    }
		    mat covMatrixY;
		    mat corMatrixY;
		    if(!armaComputeCovariance(Y, covMatrixY, corMatrixY)) {
		      error("Could not compute coexpression matrix for controls");
		    }

		    // algorithm from R script z_test.R
		    int i, j;
#pragma omp parallel for schedule(dynamic, 1) private(i, j)
  	    for(i=0; i < M; ++i) {
		      for(j=0; j < M; ++j) {
		        if(j <= i) {
		          continue;
		        }
		        double r_ij_1 = corMatrixX(i, j);
		        double r_ij_2 = corMatrixY(i, j);
		        double z_ij_1 = 0.5 * log((abs((1 + r_ij_1) / (1 - r_ij_1))));
		        double z_ij_2 = 0.5 * log((abs((1 + r_ij_2) / (1 - r_ij_2))));
		        double Z_ij = abs(z_ij_1 - z_ij_2) / sqrt((1 / ((double) n1 - 3) + 1 / ((double) n2 - 3)));
		        dcgain(i, j) = Z_ij;
		        dcgain(j, i) = Z_ij;
		      }
		    }

		    if(par::do_dcgain_abs) {
		      for(int i=0; i < dcgain.n_rows; ++i) {
		        for(int j=0; j < dcgain.n_cols; ++j) {
		          dcgain(i, j) = abs(dcgain(i, j));
		        }
		      }
		    }

				cr = new CentralityRanker(dcgain, P.nlistname);
			}			

			// run snprank
			assert(cr);
			cr->SetGlobalGamma(0.85);
			if(!cr->CalculateCentrality(GAUSS_ELIMINATION)) {
				error("SNPrank failed");
			}
			
			// save scores to dcgain matrix
			vec snprankResults = cr->GetResultsByVariable();
			//cout << snprankResults << endl;
			permResults.row(perm) = snprankResults.t();

			delete cr;
			cr = 0;
		}
		// display(permResults);

		// calculate variable thresholds
		int thresholdIndex = (int) floor(numPerms * (1.0 - par::rankerPermThreshold)) - 1;
	 	P.printLOG("\nUsing permutation threshold [" + dbl2str(par::rankerPermThreshold) + "]\n");
		// cout 
		// 	<< "M: " << M
		// 	<< " N: " << N 
		// 	<< " perm: " << perm
		// 	<< " threshold index: " << thresholdIndex 
		// 	<< endl;
	 	string saveFilename = par::output_file_name + "_thresholds.txt";
	 	P.printLOG("Writing permutation thresholds to [" + saveFilename + "]\n");
		ofstream outputFileHandle(saveFilename.c_str());
		for(int col=0; col < M; ++col) {
			colvec colScores = permResults.col(col);
			// sort the scores
			sort(colScores.begin(), colScores.end());
			// get the threshold value
			outputFileHandle << P.nlistname[col] << "\t" << colScores(thresholdIndex) << endl;
		}
		outputFileHandle.close();

	 	saveFilename = par::output_file_name + ".perm";
		armaWriteMatrix(permResults, saveFilename, P.nlistname);

		shutdown();
	}

	/////////////////////////////////////////////////////////////////////////////
	// modularity analysis requested - bcw - 5/13/13
	// NOTE: if regain file specified AND modularity assume 
	// no transform option.
	if(par::do_modularity) {
		P.printLOG("Performing network modularity analysis\n");
		InteractionNetwork* network = 0;
		// check for input file type, and construct a new network
		if(par::sifNetwork) {
			P.printLOG("Reading network from SIF file [" + par::sifFile + "]\n");
			network = new InteractionNetwork(par::sifFile, SIF_FILE, false, &P);
		}
		if(par::afniNetwork) {
			P.printLOG("Reading network from corr.1D file [" + par::afni1dFile + "]\n");
			network = new InteractionNetwork(par::afni1dFile, CORR_1D_FILE, false, &P);
		}
		if(par::do_regain_post) {
			P.printLOG("Reading network from reGAIN file [" + par::regainFile + "]\n");
			network = new InteractionNetwork(par::regainFile, REGAIN_FILE, false, &P);
		}
		if(!network) {
			error("Network construction for modularity analysis failed for the "
							"given inbix options");
		}
		P.printLOG("--- Network loaded\n");
		network->PrintSummary();

		// preprocessing transformations
		if(par::modPowerTransform) {
			P.printLOG("Transforming adjacency matrix connectivity using "
				"power with exponent " + dbl2str(par::modPowerTransformExponent) + "\n");
			network->ApplyPowerTransform(par::modPowerTransformExponent);
			P.printLOG("--- Power transformed\n");
			network->PrintSummary();
		}
		if(par::modFisherTransform) {
			P.printLOG("Transforming adjacency matrix connectivity using "
				"Fisher correlation transform\n");
			network->ApplyFisherTransform();
			P.printLOG("--- Fisher transformed\n");
			network->PrintSummary();
		}
		
		if(par::modEnableConnectivityThreshold) {
			P.printLOG("Thresholding adjacency matrix connectivity to 0 if <= " +
							dbl2str(par::modConnectivityThreshold) + "\n");
			network->SetConnectivityThresholding(par::modEnableConnectivityThreshold);
			network->SetConnectivityThreshold(par::modConnectivityThreshold);
			if(par::modUseBinaryThreshold) {
				P.printLOG("Using binary thresholding to 1 if > threshold\n");
				network->SetBinaryThresholding(par::modUseBinaryThreshold);
			}
		}

		// compute modularity
		ModularityResult modules = network->ModularityLeadingEigenvector();
		P.printLOG("Total modularity Q = " + dbl2str(modules.first) + "\n");
		P.printLOG("There are " + int2str(modules.second.size()) +
						" modules" + "\n");
		network->ShowModuleSizes();
		if(par::modComputeHomophily) {
			network->ShowHomophily();
		}
		// save modules
		network->SaveModules(par::output_file_name + ".modules");

		// clean up and shut down
		delete network;
		shutdown();
	}

	/////////////////////////////////////////////////////////////////////////////
	// recursive indirect paths modularity analysis requested - bcw - 5/31/16
	if(par::do_randomforest) {
		P.printLOG("\nPerforming random forest analysis\n");
    P.SNP2Ind();
    Forest* forest = 0;
    try {
      P.printLOG("Creating random forest object\n");
      if(par::bt) {
        P.printLOG("Creating ForestClassification\n");
        forest = new ForestClassification;
      } else {
        P.printLOG("Creating ForestRegression\n");
        forest = new ForestRegression;
      }
      forest->setVerboseOutput(&cout);
      P.printLOG("Loading data object from inbix internal data structures\n");
      Data* data = new DataDouble();
      if(data->loadFromPlink(&P)) {
        error("loadFromPlink(&P)");
      }
      P.printLOG("Initializing forest with inbix data\n");
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
      P.printLOG("Running random forest algorithm\n");
      forest->run(par::verbose);
      P.printLOG("Writing output files\n");
      forest->writeOutput();
    } catch (std::exception& e) {
      std::cerr << "Error: " << e.what() << " inbix will EXIT now." << std::endl;
      if(forest) {
        delete forest;
        forest = 0;
      }
    }

    // write the forest in Ranger internal format for prediction
    if(par::writeforest) {
      forest->writeOutputInternal();
    }

    // end gracefully
    if(forest) {
      delete forest;
      forest = 0;
    }
		P.printLOG("Random Forest analysis complete\n");
		shutdown();
  }

	/////////////////////////////////////////////////////////////////////////////
	// Relief-F analysis requested - bcw - 5/31/16
	if(par::do_relieff) {

    AnalysisType analysisType = NO_ANALYSIS;
    Dataset* ds = new PlinkInternalsDataset(&P);
    // individual-major mode for SNP bit vectors
		P.printLOG("\nLoading data set for Relief-F analysis\n");
    P.SNP2Ind();
    if(!((PlinkInternalsDataset*)ds)->LoadDataset()) {
      error("Could not load data set from PLINK internal data structures");
    }
    cout << "***** PlinkInternalsDataset*)ds)->LoadDataset() *****" << endl;

    ds->SetDistanceMetrics(par::snpMetricWeights, par::snpMetricNN, par::numMetric);
    
    if(ds->HasGenotypes() && ds->HasNumerics()) {
      analysisType = INTEGRATED_ANALYSIS;
  		P.printLOG("Performing INTEGRATED analysis\n");
    }
    if(ds->HasGenotypes()) {
      analysisType = SNP_ONLY_ANALYSIS;
  		P.printLOG("Performing SNP analysis\n");
    }
    if(ds->HasNumerics()) {
      analysisType = NUMERIC_ONLY_ANALYSIS;
  		P.printLOG("Performing NUMERIC analysis\n");
    }
    
		P.printLOG("\nPerforming Relief-F analysis\n");
    ReliefSeqController rsc(ds, &P, analysisType);
    if(par::k == 0) {
      if(!rsc.ComputeScoresKopt()) {
        error("ERROR: Failed to calculate optimum k ReliefSeq scores");
      }
    }
    else {
      if(!rsc.ComputeScores()) {
        error("ERROR: Failed to calculate ReliefSeq scores");
      }
      //rsc.PrintScores();
    }
    cout << Timestamp() << "Relief-F algorithm done" << endl;
    
    // ---------------------------------------------------------------------------
    // write results files
    string resultsFile = par::output_file_name;
		P.printLOG("\nWriting Relief-F results to: " + resultsFile + ".relief.tab\n");
    rsc.WriteAttributeScores(resultsFile);
    
		shutdown();
  }
  
	/////////////////////////////////////////////////////////////////////////////
	// Evaporative Cooling analysis requested - bcw - 5/31/16
	if(par::do_ec) {

    AnalysisType analysisType = NO_ANALYSIS;
    Dataset* ds = new PlinkInternalsDataset(&P);
    // individual-major mode for SNP bit vectors
		P.printLOG("\nLoading data set for Evaporative Cooling analysis\n");
    P.SNP2Ind();
    if(!((PlinkInternalsDataset*)ds)->LoadDataset()) {
      error("Could not load data set from PLINK internal data structures");
    }
    cout << "***** PlinkInternalsDataset*)ds)->LoadDataset() *****" << endl;

    ds->SetDistanceMetrics(par::snpMetricWeights, par::snpMetricNN, par::numMetric);
    
    if(ds->HasGenotypes() && ds->HasNumerics()) {
      analysisType = INTEGRATED_ANALYSIS;
  		P.printLOG("Performing INTEGRATED analysis\n");
    }
    if(ds->HasGenotypes()) {
      analysisType = SNP_ONLY_ANALYSIS;
  		P.printLOG("Performing SNP analysis\n");
    }
    if(ds->HasNumerics()) {
      analysisType = NUMERIC_ONLY_ANALYSIS;
  		P.printLOG("Performing NUMERIC analysis\n");
    }
    
		P.printLOG("\nPerforming Evaporative Cooling analysis\n");
    EvaporativeCooling ec(ds, &P, analysisType);
    ec.ComputeECScores();
    cout << Timestamp() << "Evaporative Cooling algorithm done" << endl;
    
    // ---------------------------------------------------------------------------
    // write results files
    string resultsFile = par::output_file_name;
		P.printLOG("\nWriting Evaporative Cooling results to: " + resultsFile + ".relief.tab\n");
    ec.WriteAttributeScores(resultsFile);
    
		shutdown();
  }
  
	/////////////////////////////////////////////////////////////////////////////
	// recursive indirect paths modularity analysis requested - bcw - 5/31/16
	if(par::do_ripm) {
		P.printLOG("\nPerforming rip-M analysis\n");
		InteractionNetwork* network = 0;
		// check for input file type, and construct a new network
		if(par::sifNetwork) {
			P.printLOG("Reading network from SIF file [" + par::sifFile + "]\n");
			network = new InteractionNetwork(par::sifFile, SIF_FILE, false, &P);
		}
		if(par::afniNetwork) {
			P.printLOG("Reading network from corr.1D file [" + par::afni1dFile + "]\n");
			network = new InteractionNetwork(par::afni1dFile, CORR_1D_FILE, false, &P);
		}
		if(par::do_regain_post) {
			P.printLOG("Reading network from reGAIN file [" + par::regainFile + "]\n");
			network = new InteractionNetwork(par::regainFile, REGAIN_FILE, false, &P);
		}
		if(!network) {
			error("Network construction for modularity analysis failed for the "
							"given inbix options");
		}
		P.printLOG("--- Network loaded\n");
    network->PrintSummary();
    
		// preprocessing transformations
		if(par::modFisherTransform) {
			P.printLOG("Transforming adjacency matrix connectivity using "
				"Fisher correlation transform\n");
			network->ApplyFisherTransform();
			P.printLOG("--- Fisher transformed\n");
			network->PrintSummary();
		}
		
		if(par::thresholdType == "hard") {
			P.printLOG("Thresholding adjacency matrix connectivity to 0 if <= " +
							dbl2str(par::modConnectivityThreshold) + "\n");
			network->SetConnectivityThresholding(par::thresholdValue);
			network->SetConnectivityThreshold(par::thresholdValue);
		} else {
			if(par::thresholdType == "soft") {
				P.printLOG("Transforming adjacency matrix connectivity using "
					"power with exponent (soft threshold): " + dbl2str(par::thresholdValue) + "\n");
				network->ApplyPowerTransform(par::thresholdValue);
				P.printLOG("--- Power transformed\n");
				network->PrintSummary();
			}
		}
		if(par::useWeighted) {
			P.printLOG("Using weighted edges\n");
			network->SetBinaryThresholding(false);
		} else {
			P.printLOG("Using binary edges thresholding to 1 if edge > threshold, else 0\n");
			network->SetBinaryThresholding(true);
		}
  	network->PrintSummary();

		// compute modularity recursively and merge small modules with rip-M
		network->SetDebugMode(par::verbose);
		network->ripM(par::startMergeOrder, par::maxMergeOrder,
		              par::minModuleSize, par::maxModuleSize);

		// save modules to tab delimited file
		if(par::verbose) network->ShowModuleIndices();
		network->SaveModules(par::output_file_name + ".ripm.modules");

    // save adjacency matrix for R analysis
		network->WriteToFile(par::output_file_name + ".adjacency.tab", 
                         REGAIN_FILE, 
                         NET_MATRIX_ADJ);
    // save connectivity matrix for R analysis
		network->WriteToFile(par::output_file_name + ".connectivity.tab", 
                         REGAIN_FILE, 
                         NET_MATRIX_CON);

		// clean up and shut down
		delete network;
		shutdown();
	}

	/////////////////////////////////////////////////////////////////////////////
	// dcVar analysis requested - bcw - 2/22/15
	if(par::do_dcvar) {
		P.printLOG("Performing dcVar analysis\n");
		// NOTE: THE SNP2Ind() CALL IS CRITICAL!!! 2/24/15
		P.SNP2Ind();
    int numVariants = P.nl_all;
    int numGenes = P.nlistname.size();
		// make sure we have variants
		if(numVariants < 1) {
      error("Variants file must specified at least one variant for this analysis!");
    }
		// make sure we have genes
		if(numGenes < 2) {
      error("Gene expression file must specified for this analysis!");
    }
		P.printLOG(int2str(numVariants) + " variants, and " + int2str(numGenes) + " genes\n");

    // for all variants
    unsigned int variantIdx;
    for(variantIdx=0; variantIdx < numVariants; ++variantIdx) {
			string variantName = P.locus[variantIdx]->name;
			// P.printLOG("Variant: " + variantName + "\n");
			// get variant info as case-control phenotype based on variant model
			// cout << "P.n: " << P.n << endl;
			// cout << "sample size: " << P.sample.size() << endl;
			for(int sampleIdx=0; sampleIdx < P.n; sampleIdx++) {
				Individual* person = P.sample[sampleIdx];
				// cout << "variantIdx: " << variantIdx << endl;
				// cout << "sampleIdx: " << sampleIdx << endl;
				// cout << "phenotype: " << person->phenotype << endl;
				// cout << "aff: " << person->aff << endl;
				// cout << "locus size: " << P.locus.size() << endl;
				// cout << "SNP allele1 size: " << person->one.size() << endl;
				// cout << "SNP allele2 size: " << person->two.size() << endl;
			  bool i1 = person->one[variantIdx];
				bool i2 = person->two[variantIdx];
				// cout << "i1: "<< i1 << ", i2:  " << i2 << endl;
				// see Caleb's email of 2/24/15 for phenotype assignment based on var model param
				// and bit-wise genotype encoding
				double thisPheno = -9;
				bool thisAff = false;
			  if(i1) {
			    if(!i2) {
				    // 10 het
			      thisPheno = 1;
						thisAff = true;
			    } else {
			      // 11
			      thisPheno = 1;
						thisAff = true;
			    }
			  } else {
			    // 01 // het 
			    if(i2) {
			    	if(par::dcvar_var_model == "rec") {
			    		thisPheno = 1;
							thisAff = true;
			    	}
			    	else {
			    		if(par::dcvar_var_model == "dom") {
			    			thisPheno = 0;
								thisAff = false;
		    			}
			    		// else "hom" missing pheno = -9
			    	}
			    }
			    // 00
			    else {
			      thisPheno = 0; // hom
 						thisAff = false;
 					}
			  }
				// cout 
				// 	<< "Variant index: " << variantIdx << "\t" 
				// 	<< "[ " << i1 << ", " << i2 << " ]" << "\t"
				// 	<< "Sample: " << sampleIdx << "\t" 
				// 	<< "Phenotype: " << thisPheno << "\t" 
				// 	<< "Affected: " << thisAff
				// 	<< endl;
				person->phenotype = thisPheno;
				person->aff = thisAff;
			}

	    // run dcGAIN for this variant phenotype
	    mat results(numGenes, numGenes);
	    mat pvals(numGenes, numGenes);
	    armaDcgain(results, pvals);
		  // DEBUG
		  // cout << "results" << endl << results.submat(0,0,4,4) << endl;
  		// cout << "pvals" << endl << pvals.submat(0,0,4,4) << endl;
		  // armaWriteMatrix(results, "DEBUG.dcgain", PP->nlistname);
		  // armaWriteMatrix(pvals, "DEBUG.pvals", PP->nlistname);

      // save p-values that pass BH rejection threshold
  		// setup output file
      string dcvarFilename = variantName + ".dcVarTest.txt";
      ofstream dcvarFile;
      P.printLOG("Writing results to [ " + dcvarFilename + " ]\n");
      dcvarFile.open(dcvarFilename.c_str());
      if(dcvarFile.fail()) {
        error("Cannot open dcVar test results file for writing.");
      }
      dcvarFile.precision(6);
      dcvarFile.fixed;

	    if(par::do_dcvar_pfilter) {
  		  double nVars = (double) numGenes;
  	    double nCombs = (nVars * (nVars - 1.0)) / 2.0;
		    double minP = 1.0;
		    double maxP = 1.0;
		    int goodPvalCount = 0;
	    	if(par::dcvar_pfilter_type == "fdr") {
			    // ------------------------------------------------------------------------
			    P.printLOG("Filtering using Benjamini-Hochberg FDR threshold\n");
			    // get all p-values
			    vector<matrixElement> test_pvals;
		 	    for(int i=0; i < pvals.n_rows; ++i) {
			      for(int j=0; j < pvals.n_cols; ++j) {
			      	if(j <= i) { continue; }
			        test_pvals.push_back(make_pair(pvals(i, j), make_pair(i, j)));
			      }
			    }
			    // sort them
			    sort(test_pvals.begin(), test_pvals.end(), pvalComparator);
			    // use rough FDR (RFDR) to estimate alpha based on input FDR
			    int num_pvals = test_pvals.size();
			    double m = (double) num_pvals * (double) numVariants;
			    double alpha = 2 * m * par::dcvar_pfilter_value  / (m + 1);
			    int threshold_index = -1;
			    // BH method
			    for(int i = 0; i < num_pvals; i++) {
			      double l = (i + 1) * alpha / (double) num_pvals;
			      // test whether current p-value < current l
			      if(test_pvals[i].first < l) {
			        threshold_index = i;
			      } else {
			        break;
			      }
			    }
			    // BH threshold condition not met with any p-values, so break out of this iteration
			    if(threshold_index == -1) {
			      P.printLOG("No p-value meets BH threshold criteria, so nothing saved\n");
			    } else {
				    // BH rejection threshold
				    double T = test_pvals[threshold_index].first;
				    P.printLOG("BH rejection threshold T = " + dbl2str(T) + ", R = " +
				      int2str(threshold_index) + "\n");
				    goodPvalCount = threshold_index + 1;
				    minP = test_pvals[0].first;
				    maxP = test_pvals[test_pvals.size()-1].first;
				    for(int i=0; i < threshold_index; i++) {
				    	double p = test_pvals[i].first;
				      pair<int, int> idx = test_pvals[i].second;
			        string gene1 = P.nlistname[idx.first];
			        string gene2 = P.nlistname[idx.second];
			        dcvarFile << gene1 << "\t" << gene2 << "\t" << p << endl;
				    }
					}
				} else {
			    P.printLOG("Filtering using Bonferroni threshold\n");
	    		// insure doubles used in all intermediate calculations
	        double correctedP = par::dcvar_pfilter_value / (nCombs * numVariants);
	        // cout 
	        // 	<< nVars << "\t" 
	        // 	<< nCombs <<  "\t"
	        // 	<< numVariants << "\t"
	        // 	<< correctedP
	        // 	<< endl;
	        // printf("FDR Corrected p-value: %g\n", correctedP);
			    minP = pvals(0, 0);
			    maxP = pvals(0, 0);
			    for(int i=0; i < pvals.n_rows; ++i) {
			      for(int j=0; j < pvals.n_cols; ++j) {
			      	if(j <= i) { continue; }
			        string gene1 = P.nlistname[i];
			        string gene2 = P.nlistname[j];
			        double p = pvals(i, j);
					    // printf("p-value [%g] < [%g] ?\n", p, correctedP);
			        // cout << gene1 << ", " << gene2 << " => p=" << p << " corrected=" 
			        //  << correctedP << " Passed FDR test!" << endl;
			        if(p < minP) minP = p;
			        if(p > maxP) maxP = p;
			        if(p <  correctedP) {
			          ++goodPvalCount;
			          // cout << gene1 << ", " << gene2 << " => p=" << p << " corrected=" 
			          //  << correctedP << " Passed FDR test!" << endl;
						    //printf("p-value [%g] < [%g] PASSED!\n", p, correctedP);
			          dcvarFile << gene1 << "\t" << gene2 << "\t" << p << endl;
			      	}
			    	} // end pvals cols
			    } // end pvals rows
				}
		    P.printLOG("Found [" + int2str(goodPvalCount) + "] tested p-values, min/max: " + 
		    	dbl2str(minP) + " / " + dbl2str(maxP) + "\n");
		  } else {
		  	// no p-value filtering
			  P.printLOG("Saving ALL p-values\n");
				for(int i=0; i < pvals.n_rows; ++i) {
				  for(int j=0; j < pvals.n_cols; ++j) {
				  	if(j <= i) { continue; }
				    string gene1 = P.nlistname[i];
				    string gene2 = P.nlistname[j];
				    double p = pvals(i, j);
	          dcvarFile << gene1 << "\t" << gene2 << "\t" << p << endl;
				  }
				}
		  }
	    
	    dcvarFile.close();

    } // END all variants loop

    shutdown();
	}

	/////////////////////////////////////////////////////////////////////////////
	// dcGAIN analysis requested - bcw - 10/30/13
	// moved algorithm to armaDcgain function - bcw - 3/12/15
	if(par::do_differential_coexpression) {
		P.printLOG("Performing dcGAIN analysis\n");
    int numVars = P.nlistname.size();
    mat results(numVars, numVars);
    mat pvals(numVars, numVars);

    armaDcgain(results, pvals);
    
    // write results
    if(par::do_dcgain_abs) {
      for(int i=0; i < results.n_rows; ++i) {
        for(int j=0; j < results.n_cols; ++j) {
          results(i, j) = abs(results(i, j));
        }
      }
    }
    string dcgainFilename = par::output_file_name + ".dcgain";
    armaWriteMatrix(results, dcgainFilename, P.nlistname);
    string dcgainPvalsFilename = par::output_file_name + ".pvals.dcgain";
    armaWriteMatrix(pvals, dcgainPvalsFilename, P.nlistname);
    shutdown();
  }
  
	////////////////////////////////////////////
	// dmGAIN analysis requested - bcw - 7/31/14
	// from bam email modification of dcGAIN - 7/29/14
	if(par::do_differential_modularity) {
		P.printLOG("Performing dmGAIN analysis\n");
    int numVars = P.nlistname.size();
    mat results(numVars, numVars);
    mat pvals(numVars, numVars);

    // t-test for diagonal
    int nAff = 0;
    int nUnaff = 0;
    for(int i=0; i < PP->sample.size(); i++) {
      if(PP->sample[i]->aff) {
        ++nAff;
      }
      else {
        ++nUnaff;
      }
    }
    double df = nAff + nUnaff - 2;
		P.printLOG("Performing z-tests\n");
    for(int i=0; i < numVars; ++i) {
      // double t;
      // tTest(i, t);
      // double p = pT(t, df);
      // results(i, i) = t;
      double z;
      zTest(i, z);
      double p = 1.0;
      results(i, i) = z;
      if(par::do_regain_pvalue_threshold) {
        if(p > par::regainPvalueThreshold) {
          results(i, i) = 0;
        }
      }
      pvals(i, i) = p;
    }

    // z-test for off-diagonal elements
    P.printLOG("Computing coexpression for CASES and CONTROLS.\n");
    mat X;
    mat Y;
    if(!armaGetPlinkNumericToMatrixCaseControl(X, Y)) {
      error("Cannot read numeric data into case-control matrices");
    }
    // compute covariances/correlations
    mat covMatrixX;
    mat corMatrixX;
    if(!armaComputeCovariance(X, covMatrixX, corMatrixX)) {
      error("Could not compute coexpression matrix for cases");
    }
    mat covMatrixY;
    mat corMatrixY;
    if(!armaComputeCovariance(Y, covMatrixY, corMatrixY)) {
      error("Could not compute coexpression matrix for controls");
    }

    // algorithm from R script z_test.R, modified by bam email 7/29/14
    colvec k_1 = sum(corMatrixX, 1);
		double two_m_1 = sum(k_1);
		colvec k_2 = sum(corMatrixY, 1);
		double two_m_2 = sum(k_2);

		P.printLOG("Performing Z-tests for interactions\n");
    double n1 = nAff;
    double n2 = nUnaff;
    for(int i=0; i < numVars; ++i) {
      for(int j=0; j < numVars; ++j) {
        if(j <= i) {
          continue;
        }
        double r_ij_1 = corMatrixX(i, j);
        double r_ij_2 = corMatrixY(i, j);
        double z_ij_1 = corMatrixX(i, j) - k_1(i) * k_1(j) / two_m_1;
        double z_ij_2 = corMatrixY(i, j) - k_2(i) * k_2(j) / two_m_2;
        double Z_ij = abs(z_ij_1 - z_ij_2) / sqrt((1/(n1 - 3) + 1 / (n2 - 3)));
        double p = 2 * normdist(-abs(Z_ij));
        results(i, j) = Z_ij;
        results(j, i) = Z_ij;
        if(par::do_regain_pvalue_threshold) {
          if(p > par::regainPvalueThreshold) {
            results(i, j) = 0;
            results(j, i) = 0;
          }
        }
        pvals(i, j) = p;
        pvals(j, i) = p;
      }
    }
    // write results
    if(par::do_dmgain_abs) {
      for(int i=0; i < results.n_rows; ++i) {
        for(int j=0; j < results.n_cols; ++j) {
          results(i, j) = abs(results(i, j));
        }
      }
    }
    string dmgainFilename = par::output_file_name + ".dmgain";
    armaWriteMatrix(results, dmgainFilename, P.nlistname);
    string dmgainPvalsFilename = par::output_file_name + ".pvals.dmgain";
    armaWriteMatrix(pvals, dmgainPvalsFilename, P.nlistname);
    shutdown();
  }

	/////////////////////////////////////////////////////////////////////////////
	// reGAIN analysis requested - bcw - 4/22/13
	if(par::do_regain) {
		P.printLOG("Performing reGAIN analysis\n");
		P.SNP2Ind();
		Regain* regain = new Regain(
						par::regainCompress,
						par::regainSifThreshold,
						par::have_numerics,
						par::regainComponents,
						par::regainFdrPrune,
						true);
		// reGAIN output options - bcw - 4/30/13
		if(par::regainMatrixThreshold) {
			regain->setOutputThreshold(par::regainMatrixThresholdValue);
			regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_THRESH);
		}
		if(par::regainMatrixFormat == "upper") {
			regain->setOutputFormat(REGAIN_OUTPUT_FORMAT_UPPER);
		} else {
			if(par::regainMatrixFormat == "full") {
				regain->setOutputFormat(REGAIN_OUTPUT_FORMAT_FULL);
			} else {
				error("reGAIN output format allowed options: {upper, full}");
			}
		}
		if(par::regainMatrixTransform == "none") {
			regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_NONE);
		} else {
			if(par::regainMatrixTransform == "threshold") {
				regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_THRESH);
			} else {
				if(par::regainMatrixTransform == "abs") {
					regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_ABS);
				} else {
					error("reGAIN output transform allowed options: {none, threshold, abs}");
				}
			}
		}
		if(par::regainPureInteractions) {
			regain->performPureInteraction(true);
		} else {
			regain->performPureInteraction(false);
		}
		regain->run();
		if(par::regainFdrPrune) {
			regain->writeRegain(false);
			regain->fdrPrune(par::regainFdr);
		}
		// write output options to stdout and log file - bcw - 5/1/13
		regain->logOutputOptions();
		regain->logMatrixStats();
		regain->writeRegain(false, par::regainFdrPrune);
		regain->writeRegain(true);
		delete regain;
		// stop inbix processing
		shutdown();
	}

	/////////////////////////////////////////////////////////////////////////////
	// reGAIN post processing requested - bcw - 5/3/13
	if(par::do_regain_post) {
		P.printLOG("Performing reGAIN file post processing\n");
		P.SNP2Ind();
		Regain* regain = new Regain(
						par::regainCompress,
						par::regainSifThreshold,
						par::regainComponents
						);
		regain->readRegainFromFile(par::regainFile);
		// reGAIN output options - bcw - 4/30/13
		if(par::regainMatrixThreshold) {
			regain->setOutputThreshold(par::regainMatrixThresholdValue);
			regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_THRESH);
		}
		if(par::regainMatrixFormat == "upper") {
			P.printLOG("Reformatting reGAIN matrix to upper triangular\n");
			regain->setOutputFormat(REGAIN_OUTPUT_FORMAT_UPPER);
		} else {
			if(par::regainMatrixFormat == "full") {
				P.printLOG("Reformatting reGAIN matrix to full matrix\n");
				regain->setOutputFormat(REGAIN_OUTPUT_FORMAT_FULL);
			} else {
				error("reGAIN output format allowed options: {upper, full}");
			}
		}
		if(par::regainMatrixTransform == "none") {
			regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_NONE);
		} else {
			if(par::regainMatrixTransform == "threshold") {
				P.printLOG("Thresholding reGAIN matrix\n");
				regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_THRESH);
			} else {
				if(par::regainMatrixTransform == "abs") {
					P.printLOG("Absolute value reGAIN matrix\n");
					regain->setOutputTransform(REGAIN_OUTPUT_TRANSFORM_ABS);
				} else {
					error("reGAIN output transform allowed options: {none, threshold, abs}");
				}
			}
		}
		if(par::regainMatrixToSif) {
			regain->writeRegainToSifFile(par::output_file_name + ".sif");
		} else {
			regain->writeRegainToFile(par::regainFile + ".postproc");
		}

    delete regain;
		shutdown();
	}

  // we're done if we have no SNPs at this point
  // (though we might have done some numeric file reading)
  if(P.nl_all == 0) {
    shutdown();
  }
  
	/////////////////////////////////////////
	// SET statistics?
	if(par::read_set)
		P.readSet();
	else if(par::make_set)
		P.outputSetFile();

	Set S(P.snpset);
	P.pS = &S;

	// Remove any SNPs not in a set
	// unless using particular commands
	// (set-by-all epistasis, set-table)
	if(par::read_set || par::make_set) {
		if(par::drop_sets)
			P.pS->dropNotSet(P);
	}

	//////////////////////////////////////////////////
	// Build final marker scaffold 
	makeScaffold(P);

	//////////////////////////////////////////////////
	//                                              //
	// Create family units?                         //
	//                                              //
	//////////////////////////////////////////////////

	///////////////////////
	// Create family units?
	if(par::MENDEL_test ||
					par::MENDEL_report ||
					par::TDT_test ||
					par::QTDT_test ||
					par::make_founders &&
					!par::built_families) {
		map<string, Individual*> fnd;
		map<Individual*, int> idmap;
		P.linkRelateds(idmap, fnd);
		P.parseTrios();
		par::built_families = true;

		// Perform now, so that the user has an option to 
		// save a new fileset with mendel errors removed
		if(par::MENDEL_report || par::MENDEL_test)
			P.checkMendel();
	}

	////////////////////////////////////////////////
	// Reset PAT/MAT codes of any non- nonfounders?
	// i.e. if parents not actually present in sample?
	if(par::make_founders) {
		P.makeFounders();
	}

	//////////////////////////////////////
	// Sex check
	if(par::check_sex) {
		P.sexCheck();
	}

	//////////////////////////////////////
	// Split TDT units to case/controls
	if(par::tucc) {
		if(!par::built_families) {
			map<string, Individual*> fnd;
			map<Individual*, int> idmap;
			P.linkRelateds(idmap, fnd);
			P.parseTrios();
			par::built_families = true;
		}

		P.checkMendel();
		P.pseudoCaseControl();
	}

	//////////////////////////////////////////////////
	//                                              //
	// Haplotype imputation methods                 //
	//                                              //
	//////////////////////////////////////////////////

	// Do not use this old IMPUTATION method
	// Restrict to --proxy-impute, or original
	// --hap-impute (i.e. based on multi-marker list)

	if(par::meta_large_phase) {

		// Automatically try to impute all one window per chromosome
		// We can put in some other restraints here if need be

		if(par::has_nonfounders && !par::built_families) {
			map<string, Individual*> fnd;
			map<Individual*, int> idmap;
			P.linkRelateds(idmap, fnd);
			P.parseTrios();
			P.checkMendel();
			par::built_families = true;
		}

		P.printLOG("Estimating haplotype frequencies/phases ( MHF >= "
						+ dbl2str(par::min_hf) + " )\n");
		P.printLOG("Considering phases P(H|G) >= "
						+ dbl2str(par::hap_min_phase_prob) + "\n");
		P.printLOG("Requiring per individual per haplotype missingness < "
						+ dbl2str(par::hap_missing_geno) + " \n");

		P.printLOG("Initial EM window size " + int2str(par::haplo_plem_window)
						+ " SNPs with " + int2str(par::haplo_plem_overlap) + " SNP overlap\n");


		// Count number of founders

		P.haplo->cnt_f = 0;
		vector<Individual*>::iterator person = P.sample.begin();
		while(person != P.sample.end()) {
			if((*person)->founder) P.haplo->cnt_f++;
			person++;
		}

		if(P.haplo->cnt_f < P.n) {
			P.haplo->nonfounders = true;
			P.printLOG("Initial phasing based on " +
							int2str(P.haplo->cnt_f) + " founders (" +
							int2str(P.n - P.haplo->cnt_f) +
							" non-founders)\n");
		}

		// Start off just with the autosomes
		// We assume that "--chr" has been specified on the command line, 
		// and so we are only dealing with a single chromosome here
		if(par::impute_verbose) {
			P.printLOG("Writing verbose imputation output to [ "
							+ par::output_file_name + ".phased.out ]\n");
			P.haplo->HIMPUTE.open((par::output_file_name + ".phased.out").c_str(),
							ios::out);
			P.haplo->HIMPUTE.setf(ios::fixed);
			P.haplo->HIMPUTE.precision(2);
		}

		// Run imputation in blocks of up to 1000 SNPs
		P.haplo->makeSlidingWindow("20+20");

		P.haplo->phaseAllHaplotypes(true, perm);

		if(par::impute_verbose)
			P.haplo->HIMPUTE.close();

	}

	////////////////////////////////////
	// Proxy-based haplotype imputation 
	if(par::proxy_impute) {
		P.proxyWrapper();
		// Do not shut down: we assume a --make-bed will
		// be called below
	}

	//////////////////////////////////////////////////
	//                                              //
	// Generate dummy permuted phenotype file       //
	//                                              //
	//////////////////////////////////////////////////
	if(par::output_pheno_perm) {
		P.outputPermedPhenotypes(perm);
		shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// Output formats and transformations           //
	//                                              //
	//////////////////////////////////////////////////

	// Covariate files can also be output (--covar) 
	// for the major options: --make-bed, --recode* 
	// and also just --write-covar option

	if(par::set_table) {
		P.setTable();
		shutdown();
	}

	if(par::write_set) {
		P.writeSetFile();
		shutdown();
	}

	if(par::dump_covar) {
		P.write_covariates();
		shutdown();
	}

	if(par::dump_clst) {
		P.write_clusters();
		shutdown();
	}

	if(par::write_snplist) {
		P.write_snplist();
		shutdown();
	}

	if(par::write_bitfile) {
		P.write_BITFILE();
		if(par::clist)
			P.write_covariates();
		shutdown();
	}

	if(par::recode_fastphase) {
		P.output_fastphase_format();
		shutdown();
	}

	if(par::recode_bimbam) {
		P.output_bimbam_format();
		shutdown();
	}

	if(par::recode_structure) {
		P.output_structure_format();
		shutdown();
	}

	if(par::recode || par::recode_HV || par::recode_12 || par::recode_whap) {
		if(!par::recode_transpose)
			P.display_recoded_PEDFILE();
		else
			P.display_recoded_PEDFILE_transpose();
		if(par::clist)
			P.write_covariates();
		shutdown();
	}

	if(par::recode_AD) {
		P.display_recoded_PEDFILE_AD();
		if(par::clist)
			P.write_covariates();
		shutdown();
	}

	if(par::recode_long) {
		P.display_recoded_LONG();
		if(par::clist)
			P.write_covariates();
		shutdown();
	}

	if(par::recode_mutlist) {
		P.display_recoded_MUTLIST();
		if(par::clist)
			P.write_covariates();
		shutdown();
	}


	if(par::list_by_allele) {
		P.display_listByAllele();
		shutdown();
	}

	if(par::plist) {
		P.display_pairList();
	}

	if(par::indiv_report) {
		P.display_indivReport();
		shutdown();
	}

	if(par::list_twolocus) {
		P.display_twolocus();
		shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// LD-based lookups                             //
	//                                              //
	//////////////////////////////////////////////////

	////////////////////////////////////////////
	// Set summary statistics
	if(par::set_screen) {
		P.setAssocSummary();
		shutdown();
	}

	////////////////////////////////////////////
	// LD-based clumping
	if(par::clumpld) {
		clump_LD cld(&P, P.haplo,
						par::clumpld_p1,
						par::clumpld_kb,
						par::clumpld_p2,
						par::clumpld_r2);
		cld.clump();
		shutdown();
	}

	////////////////////////////////////////////
	// Show tags
	if(par::gettag_mode) {
		P.tagMode();
		shutdown();
	}

	////////////////////////////////////////////
	// Haplotype block action
	if(par::make_blocks) {
		P.mkBlks(0, P.nl_all - 1);
		shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// Main set of whole-genome tests               //
	//                                              //
	//////////////////////////////////////////////////

	/////////////////////////////////
	// Some initial set-up work here

	/////////////////////
	// Conditioning SNPs
	if(par::conditioning_snps) {
		if(par::conditioning_snp_single) {

			// ** todo ** change this to allow a NList

			int x = getMarkerNumber(P, par::conditioning_snp_name);

			if(x < 0) error("Marker "
							+ par::conditioning_snp_name
							+ " does not exist in filtered data\n");

			P.conditioner.push_back(x);
			P.conditioner_mask.push_back(false);
		} else
			P.readConditioningList();
	}

	//////////////////////////////////////////
	// Warn if not enough markers in analysis
	if(par::plink
					|| par::cluster
					|| par::cluster_plot
					|| par::outlier_detection
					|| par::genome_output
					|| par::inbreeding) {
		if(P.nl_all < 10000)
			P.printLOG("\n **Warning** this analysis typically requires whole-genome level data\n"
						"             to give accurate results \n\n");
	}

	//////////////////////////////////////////
	// Arbitrary external functions
	if(par::myfunction) {
		if(1) {
			if(par::has_nonfounders && !par::built_families) {
				map<string, Individual*> fnd;
				map<Individual*, int> idmap;
				P.linkRelateds(idmap, fnd);
				P.parseTrios();
				P.checkMendel();
				par::built_families = true;
			}
		}

		// P.callMe();
		shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// IBS and IBD genome-wide analyses             //
	//                                              //
	//////////////////////////////////////////////////

	//////////////////////////////////////////////
	// Perform a cluster analysis and/or MDS plot 
	if(par::cluster || par::cluster_plot || par::outlier_detection) {
		P.buildCluster();
		shutdown();
	}

	///////////////////////////////////
	// Permutation test between groups
	// based on IBS diffeences
	if(par::ibs_test) {
		P.permutationIBSTest(perm);
		shutdown();
	}

	////////////////////////////////////////////////////
	// Precalculate frequency-averaged P(IBD|IBS) table
	if(par::plink || par::genome_output) {
		if(par::has_nonfounders && !par::built_families) {
			map<string, Individual*> fnd;
			map<Individual*, int> idmap;
			P.linkRelateds(idmap, fnd);
			P.parseTrios();
			// P.checkMendel(); // skip this when in --rel-check mode
			par::built_families = true;
		}

		// So that correct IBD expectation is calculated, 
		// we need to fill in empty slots for missing parents
		P.makeMissingParents();

		P.preCalcGenomeIBD();
	}

	//////////////////////////////
	// Genome-wide output only
	if(par::genome_output) {
		P.displayGenomeWideInfo();

		if(par::genome_test)
			P.testGenomeIBDByCovariate(perm);

		shutdown();
	}

	//////////////////////////////////////
	// Genome-wide inbreeding output only
	if(par::inbreeding) {

		if(par::SNP_major)
			P.SNP2Ind();

		ofstream HET;
		string f = par::output_file_name + ".het";
		HET.open(f.c_str(), ios::out);
		HET.precision(4);

		P.printLOG("Writing individual heterozygosity information to [ " + f + " ] \n");
		HET << setw(par::pp_maxfid) << "FID" << " "
						<< setw(par::pp_maxiid) << "IID" << " "
						<< setw(12) << "O(HOM)" << " "
						<< setw(12) << "E(HOM)" << " "
						<< setw(12) << "N(NM)" << " "
						<< setw(12) << "F" << "\n";

		for(int i1 = 0; i1 < P.n; i1++)
			P.calcInbreeding(P.sample[i1], 0, P.nl_all - 1, HET);

		HET.close();

	}

	///////////////////////
	// Runs of homozygosity
	if(par::homo_run) {

		P.findAllHomozygousRuns(perm);

		if(par::segment_test_individual)
			P.segmentIndividualTest(perm);

		shutdown();
	}

	///////////////////////////////////
	// Runs of missingness (deletions)
	if(par::miss_run) {

		if(par::SNP_major) P.SNP2Ind();

		ofstream RUN;
		string f = par::output_file_name + ".rum";
		RUN.open(f.c_str(), ios::out);

		P.printLOG("Writing run-of-missings information to [ " + f + " ] \n");

		string msg = "Run defined as " + int2str(par::miss_run_length);
		if(par::miss_run_length_kb) msg += " kb\n";
		else msg += " SNPs\n";
		P.printLOG(msg);

		stringstream s2;
		s2 << "With at least " << par::miss_run_level << " missingness\n";
		P.printLOG(s2.str());

		for(int i1 = 0; i1 < P.n; i1++) {
			if(!par::silent)
				cout << i1 + 1 << " of " << P.n << " individuals      \r";
			P.findMissRuns(P.sample[i1], RUN);
		}

		if(!par::silent)
			cout << "\n\n";

		RUN.close();

		shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// LD and haplotype-based analyses              //
	//                                              //
	//////////////////////////////////////////////////

	///////////////////
	// LD-based pruning
	if(par::prune_ld) {
		P.pruneLD();
		shutdown();
	}

	//////////////////////////////
	// Flip-scan
	if(par::flip_scan) {
		P.calcFlipScan();
		shutdown();
	}

	//////////////////////////////
	// LD statistics
	if(par::calc_SNPSNP_LD) {
		P.calcPairwiseLD();
		shutdown();
	}

	if(par::disp_r1 || par::disp_r2) {
		P.calcLDStatistics();
		shutdown();
	}

	///////////////////////////////////////////////////
	// General class for haplotype phasing and testing
	if(par::test_hap_CC && par::qt) {
		par::test_hap_CC = false;
		par::test_hap_QTL = true;
	}

	// In case families are included, build family structure if not
	// already done
	if(par::phase_snps || par::mishap_test || par::proxy_assoc) {

		// Read in list of tests, or make sliding window?
		if(par::phase_snps) {
			if(par::sliding_window)
				P.haplo->makeSlidingWindow(par::sliding_window_size);
			else if(par::hap_specific_snps)
				P.haplo->setSpecificSNPs(par::hap_specific_snps_list);
			else
				P.haplo->readTagFile();
		}

		if(par::has_nonfounders && !par::built_families) {
			map<string, Individual*> fnd;
			map<Individual*, int> idmap;
			P.linkRelateds(idmap, fnd);
			P.parseTrios();
			P.checkMendel();
			par::built_families = true;
		}

		P.printLOG("Estimating haplotype frequencies/phases ( MHF >= "
						+ dbl2str(par::min_hf) + " )\n");
		P.printLOG("Considering phases P(H|G) >= "
						+ dbl2str(par::hap_min_phase_prob) + "\n");
		P.printLOG("Requiring per individual per haplotype missingness < "
						+ dbl2str(par::hap_missing_geno) + " \n");

		// Count number of founders
		P.haplo->cnt_f = 0;
		vector<Individual*>::iterator person = P.sample.begin();
		while(person != P.sample.end()) {
			if((*person)->founder) P.haplo->cnt_f++;
			person++;
		}

		if(P.haplo->cnt_f < P.n) {

			P.haplo->nonfounders = true;
			P.printLOG("Initial phasing based on " +
							int2str(P.haplo->cnt_f) + " founders (" +
							int2str(P.n - P.haplo->cnt_f) +
							" non-founders)\n");
		}

		if(P.n == P.haplo->cnt_f &&
						(par::test_hap_TDT || par::proxy_TDT))
			error("Can not perform TDT in sample with no non-founders");

	}

	/////////////////////////
	// Haplotype frequencies
	if(par::phase_snps && par::display_hap_freqs) {
		P.haplo->calculateHaplotypeFrequencies();
		shutdown();
	}


	////////////////////////////////
	// Haplotype phase probabilities
	if(par::phase_snps && par::display_phase_probs) {
		P.haplo->calculateHaplotypeFrequencies();
		shutdown();
	}

	/////////////////////////////////////////////
	// Haplotypic test of non-random missing data
	if(par::mishap_test) {
		P.performMisHapTests();
		shutdown();
	}

	////////////////////////////////////////////////////
	// Haplotype tracking of an extended region, for an 
	// individual or pair
	if(par::phase_snps && par::segment_haplotrack) {
		P.haplo->trackSharedHaplotypes();
		shutdown();
	}

	////////////////////////////////////////////////////
	// Haplotype tracking of an extended region, for an 
	// individual or pair
	if(par::phase_snps && par::impute_tags) {
		P.haplo->imputeAllHaplotypes();
		shutdown();
	}

	///////////////////////////////////////////////////////
	// Haplotypic test of SNP proxy (convenience function)
	// (we've already done the imputation step above)
	if(par::proxy_assoc && !par::proxy_impute) {
		P.proxyWrapper();
		shutdown();
	}

	//////////////////////////////////////////////////
	//                                              //
	// Misc tests that do not fall within the       //
	// main phenotype loop                          //
	//                                              //
	//////////////////////////////////////////////////

	//////////////////////////////
	// Genome-wide IBS sharing test
	if(par::ibs_sharing_test) {
		P.perm_sharingIBSTest(perm);
		shutdown();
	}

	//////////////////////////////
	// Gene-based test of epistasis
	if(par::epi_genebased) {
		P.driverSCREEPI();
		shutdown();
	}

	//////////////////////////////
	// Genome-wide epistasis tests
	if(par::epistasis) {
		P.calcEpistasis();
		shutdown();
	}

	/////////////////////////////////////////////
	// Determine per-individual risk profiles
	if(par::score_risk) {
		P.scoreIndividuals();
		shutdown();
	}

	//////////////////////////////////
	// Apply an R-script to the data?
	if(par::run_R_script) {

#ifdef WITH_R_PLUGINS
		P.Rfunc();
		shutdown();
#else
		error("R plugin support has not been compiled in");
#endif  
	}

	//////////////////////////////////////////////////
	//                                              //
	// Genome-wide association tests                //
	//                                              //
	//////////////////////////////////////////////////

	// Allow for the fact that we might be iterating 
	// over multiple phenotypes
	string original_file_root = par::output_file_name;

	if(!par::plink)
		while(1) {

			if(par::all_pheno) {
				if(par::loop_over) {
					P.phenoLabel = P.kname[ par::loop_counter ];
					par::output_file_name = original_file_root + "." + P.phenoLabel;
					par::bt = true;
					par::qt = false;

					for(int i = 0; i < P.n; i++) {
						// Include all samples
						P.sample[i]->missing = false;
						P.sample[i]->aff = P.sample[i]->sol == par::loop_counter ? true : false;
					}
				} else {
					if(P.phenotype_name == "")
						P.phenoLabel = "P" + int2str(par::mult_pheno);
					else
						P.phenoLabel = P.phenotype_name;

					par::output_file_name = original_file_root + "." + P.phenoLabel;
				}
			}

			if(par::assoc_test) {

				if(par::CMH_test_2)
					P.calcMH();
				else if(par::OR_homog_test)
					P.calcHomog();
				else if(par::QTDT_test) {
					// Force a Mendel error check
					if(!(par::MENDEL_report || par::MENDEL_test))
						P.checkMendel();
					P.perm_testQTDT(perm);
				} else if(par::boot) {
					// Redundant
					error("Bootstrap option is no longer supported\n");

					P.calcAssociationWithBootstrap();
				} else {

					// Includes 
					//        basic allelic test
					//        model-based tests
					//        linear & logistic models
					//        2x2xK Cochran-Mantel-Haenszel

					P.calcAssociationWithPermutation(perm);
				}

				if(!par::all_pheno)
					shutdown();

			}

			/////////////////////////////////
			// Haplotype association analysis
			if(par::phase_snps && (par::test_hap_CC ||
							par::test_hap_GLM ||
							par::test_hap_QTL ||
							par::test_hap_TDT)) {

				// This is done separaytely, via the main
				// assoc. loop

				if(par::test_hap_GLM)
					P.calcAssociationWithPermutation(perm);
				else {

					////////////////////////////////////////////////
					// Perform omnibus and haplotype-specific tests
					string f;

					if(par::test_hap_CC)
						f = par::output_file_name + ".assoc.hap";
					else if(par::test_hap_QTL)
						f = par::output_file_name + ".qassoc.hap";
					else if(par::test_hap_TDT)
						f = par::output_file_name + ".tdt.hap";

					if(par::test_hap_CC) {

						P.printLOG("Writing haplotype association statistics to [ " + f + " ]\n");
						P.haplo->HTEST.open(f.c_str(), ios::out);
						P.haplo->HTEST.precision(4);

						P.haplo->HTEST << setw(10) << "LOCUS" << " "
										<< setw(12) << "HAPLOTYPE" << " "
										<< setw(10) << "F_A" << " "
										<< setw(10) << "F_U" << " "
										<< setw(10) << "CHISQ" << " "
										<< setw(4) << "DF" << " "
										<< setw(10) << "P" << " "
										<< "SNPS" << "\n";
					}


					if(par::test_hap_QTL) {
						P.printLOG("Writing haplotype association statistics to [ " + f + " ]\n");
						P.haplo->HTEST.open(f.c_str(), ios::out);
						P.haplo->HTEST.precision(4);

						P.haplo->HTEST << setw(10) << "LOCUS" << " "
										<< setw(12) << "HAPLOTYPE" << " "
										<< setw(8) << "NANAL" << " "
										<< setw(10) << "BETA" << " "
										<< setw(10) << "R2" << " "
										<< setw(8) << "STAT" << " "
										<< setw(10) << "P" << " "
										<< "SNPS" << "\n";

					}

					if(par::test_hap_TDT) {
						P.printLOG("Writing haplotype TDT statistics to [ " + f + " ]\n");
						P.haplo->HTEST.open(f.c_str(), ios::out);
						P.haplo->HTEST.precision(4);

						P.haplo->HTEST << setw(10) << "LOCUS" << " "
										<< setw(12) << "HAPLOTYPE" << " "
										<< setw(10) << "T" << " "
										<< setw(10) << "U" << " "
										<< setw(10) << "CHISQ" << " "
										<< setw(10) << "P" << " "
										<< "SNPS" << "\n";
					}


					P.haplo->phaseAllHaplotypes(true, perm);

					P.haplo->HTEST.close();
				}

				if(!par::all_pheno)
					shutdown();

			}

			//////////////////////////////////////////////////////
			// Haplotypic conditional tests (WHAP implementation, 
			// now called CHAP, for conditional haplotype
			if(par::chap_test) {

				P.conditionalHaplotypeTest(true, perm);

				if(!par::all_pheno)
					shutdown();
			}

			//////////////////////////////
			// QTL interaction test
			if(par::assoc_gxe) {
				P.perm_testGXE2(perm);
				if(!par::all_pheno)
					shutdown();
			}

			/////////////////////////////
			// Rare allele test 
			if(par::elf_baseline) {
				P.elfBaseline();
				shutdown();
			}

			if(par::rare_test) {
				P.permTestRareDistribution(perm);
				if(!par::all_pheno)
					shutdown();
			}

			/////////////////////////
			// Hotelling's T^2 test
			if(par::hotel) {
				P.perm_testHotel(perm);
				if(!par::all_pheno)
					shutdown();
			}

			///////////////////////////////////
			// Test difference in missing rates
			if(par::test_missing) {
				P.calcAssociationWithPermutation(perm);
				if(!par::all_pheno)
					shutdown();
			}

			//////////////////////////////////
			// Genome-wide family-based (TDT)
			// and Parent-of-origin analysis
			if(par::TDT_test) {

				// Force a Mendel error check, if we have not 
				// already

				if(!(par::MENDEL_report || par::MENDEL_test))
					P.checkMendel();

				// Either basic TDT or Parent-Of-Origin analysis

				if(par::parent_of_origin)
					P.perm_testTDT_POO(perm);
				else if(par::sibTDT_test)
					P.perm_testTDT(perm);
				else
					P.perm_testTDT(perm);

				if(!par::all_pheno)
					shutdown();

			}

			// Read next phenotype: repeat, or shutdown
			if(par::all_pheno) {

				if(par::loop_over) {
					// Construct next phenotype from cluster file

					par::loop_counter++;
					if(par::loop_counter == P.nk)
						shutdown();

				} else {
					// Read next phenotype from file

					par::mult_pheno++;
					if(!P.readPhenoFile())
						shutdown();

					// and recode, if a binary affection status coding

					if(par::bt)
						affCoding(P);

				}
			}

			if(!par::all_pheno)
				shutdown();

		} // Next potential phenotype

	//////////////////////////////////////////////////
	//                                              //
	// PLINK segmental sharing analyses             //
	//                                              //
	//////////////////////////////////////////////////

	// Stop now, unless a plink analysis is specified
	if(!par::plink)
		shutdown();

	if(par::SNP_major)
		P.SNP2Ind();

	//////////////////////////////////////////////
	// Read pre-computed segment list and perform 
	// segmental tests?
	if(par::read_segment_file) {

		ifstream SEG;
		SEG.open(par::read_segment_filename.c_str(), ios::in);
		P.printLOG("Reading IBD-segment information from [ "
						+ par::read_segment_filename + " ]\n");
		checkFileExists(par::read_segment_filename);

		if(par::segment_minimal)
			P.readSegmentFileMinimal(SEG);
		else
			P.readSegmentFile(SEG);
		SEG.close();

		// IBS validation of segments (i.e. possibly in a larger
		// datafile? but one that must be a superset of all SNPs in
		// segment file)

		if(false) {
			P.validateSegments();
			shutdown();
		}

		// Find overlap in segments?
		if(par::segment_overlap)
			P.summariseHomoRuns();

		// Per-individual summary/test?
		if(par::segment_test_individual) {
			P.segmentIndividualTest(perm);
			shutdown();
		}

		// Perform pairwise summary/analysis of segments?
		P.summaryIBSsegments(perm);

		P.printLOG("Writing segment summary to [ " + par::output_file_name
						+ ".segment.indiv ]\n\n");

		P.indivSegmentSummary();

		shutdown();
	}

	//////////////////////////////
	// Pair inclusion/exclusion

	// Number of informative pairs
	int c = 0;

	// Read or calculate informative pairs?
	if(par::ibd_read)
		c = P.readInformative();
	else
		c = P.calcInformative();

	////////////////////////////////////////////////
	// Test of genome-wide relatedness by covariate
	if(par::genome_test) {
		P.testGenomeIBDByCovariate(perm);
		shutdown();
	}

	///////////////////////////////////
	// Save pairs to be included? i.e. 
	// after removing all pairs for 

	// a) low IBD
	// b) not being an affected pair
	// c) being a concordant unaffected pair

	if(par::inc_write)
		P.writeInformative();

		/////////////////////////////////////////////////////////////////
		// Get and display information on chromosomal range to be tested

		//  else if (par::singlepoint)
		//    P.printLOG("Using singlepoint analysis mode\n");
	else if(par::inter_grid > 0) {
		stringstream s2;
		s2 << "Using multipoint analysis: step = "
						<< par::inter_grid
						<< " and fringe = "
						<< par::fringe << " cM\n";
		P.printLOG(s2.str());
	} else {
		stringstream s2;
		s2 << "Using multipoint analysis: grid = "
						<< par::grid
						<< " and fringe = "
						<< par::fringe << " cM\n";
		P.printLOG(s2.str());
	}

	vector<int> chrs;
	if(par::run_chr == 0) {

		vector<int> r = getChromosomeRange(P);

		P.printLOG("\nScanning from autosomes from chromosome " +
						chromosomeName(r[0]) + " to " +
						chromosomeName(r[1]) + "\n\n");
		for(int i = r[0]; i <= r[1]; i++)
			if((!par::chr_haploid[i]) &&
							(!par::chr_sex[i]))
				chrs.push_back(i);
	} else chrs.push_back(par::run_chr);


	ofstream SEG;
	if(par::segment_output) {
		string f = par::output_file_name + ".segment";
		SEG.open(f.c_str(), ios::out);
		P.printLOG("Writing IBD-segment information to [ " + f + " ]\n");
		if(par::segment_minimal) P.printLOG("Minimal segment file format\n");

		// Header row for non-minimal format
		if(!par::segment_minimal) {
			SEG << setw(par::pp_maxfid) << "FID1" << " "
							<< setw(par::pp_maxiid) << "IID1" << " "
							<< setw(par::pp_maxfid) << "FID2" << " "
							<< setw(par::pp_maxiid) << "IID2" << " ";

			if(par::bt) SEG << setw(4) << "PHE" << " ";

			SEG << setw(4) << "CHR" << " "
							<< setw(10) << "BP1" << " "
							<< setw(10) << "BP2" << " "
							<< setw(par::pp_maxsnp) << "SNP1" << " "
							<< setw(par::pp_maxsnp) << "SNP2" << " "
							<< setw(6) << "NSNP" << " "
							<< setw(10) << "KB" << "\n";
		}

		f = par::output_file_name + ".segment.summary";
		P.printLOG("Writing IBD-segment summary to [ " + f + " ]\n\n");

		P.printLOG("Minimum segment length is "
						+ dbl2str((double) par::segment_length / (double) 1000)
						+ " kb and " + int2str(par::segment_snp) + " SNPs\n");
		P.printLOG("Segment thresholds are " + dbl2str(par::segment_threshold_start)
						+ " and " + dbl2str(par::segment_threshold_finish) + "\n");
		P.printLOG("Maximum intra-segment inter-SNP distance is "
						+ int2str(par::segment_inter_snp_distance)
						+ "\n");
	}

	ofstream MP;
	if(par::multi_output) {
		string f = par::output_file_name + ".multi";
		MP.open(f.c_str(), ios::out);
		MP.setf(ios::fixed);
		MP.precision(5);
		P.printLOG("Writing multipoint IBD estimates to [ " + f + " ]\n");
	}

	ofstream GMULTI;
	if(par::gmulti_output) {
		string f = par::output_file_name + ".gmulti";
		GMULTI.open(f.c_str(), ios::out);
		GMULTI.precision(4);
		P.printLOG("Writing genotype/multipoint IBD estimates to [ " + f + " ]\n");
	}

	//////////////////////////////
	// Consider each chromosome
	for(int i = 0; i < chrs.size(); i++) {

		//////////////////////////
		// Reset main variables

		P.phenotype.resize(0);
		P.pihat.resize(0);
		P.Zlocus.resize(0);
		P.m_pihat.resize(0);
		P.v_pihat.resize(0);
		P.pair1.resize(0);
		P.pair2.resize(0);

		// Set chromosome
		par::run_chr = chrs[i];

		// Find scan range
		P.setMarkerRange();

		// Total number of all (test+background) loci
		P.nl_all = P.locus.size();

		// Number of loci for test loci
		P.nl = par::run_end - par::run_start + 1;
		P.printLOG(int2str(P.nl) + " markers in this scan\n");

		// Set up (single) multipoint marker map
		//  if (par::singlepoint)
		//   P.preCalcSinglePoint();

		P.preCalcMultiPoint();

		/////////////////////////////////
		//  For each pair of individuals
		//  calculate mutlipoint pihats

		// reset counter
		int c1 = 0;
		int c2 = 0;

		for(int i1 = 0; i1 < P.n - 1; i1++)
			for(int i2 = i1 + 1; i2 < P.n; i2++) {

				if(!P.skip_pair[c1++]) {
					Individual * p1 = P.sample[i1];
					Individual * p2 = P.sample[i2];

					/////////////////////////////////
					// Skip if genome sets specified

					if(par::genome_2sets) {
						if(!((P.gset1.find(p1) != P.gset1.end() &&
										P.gset2.find(p2) != P.gset2.end()) ||
										(P.gset1.find(p2) != P.gset1.end() &&
										P.gset2.find(p1) != P.gset2.end())))
							continue;
					}

					/////////////////////////////////
					// Skip if within-cluster analysis
					// specified
					if(par::IBD_within) {
						if(p1->sol != p2->sol)
							continue;
					}

					/////////////////////////////////
					//  1. Calculate IBD(g) | IBS(g) 
					Z IBDg = P.saved_IBDg[c2];

					if(!par::silent) {
						cout << "IBD calculation: "
										<< ++c2
										<< " of "
										<< c
										<< "                  \r";
						cout.flush();
					}

					/////////////////////////////////
					//  2. Calculate IBD(l) - IBD(g)
					vector<Z> IBDl = P.calcLocusIBD(p1, p2, IBDg);


					/////////////////////////////////
					//  3. Multipoint calculation
					P.pairid = itoa((int) p1->phenotype, 10)
									+ " " + itoa((int) p2->phenotype, 10) + " ";
					P.pairid += itoa((int) c2, 10) + " ";
					P.pairid += p1->fid + "_" + p1->iid + "_ ";
					P.pairid += p2->fid + "_" + p2->iid + "_";

					vector_t p;

					//////////////////////////////////
					// Perform either using
					//     Singlepoint analysis
					//     Multipoint analysis (default)

					// if (par::singlepoint) 
					//  p = P.calcSinglePoint(IBDl,IBDg);

					p = P.calcMultiPoint(IBDl, IBDg, MP);

					////////////////////////////////
					//   3b. Verbose output:
					//       genotypes for each pair
					if(par::gmulti_output) {
						for(int l = par::run_start; l <= par::run_end; l++)
							P.displayGMULTI(p1, p2, l, GMULTI);
					}

					///////////////////////////////
					//  4. Scan for segments of IBD
					P.findSegments(i1, i2, p, SEG);

					/////////////////////////////
					//  5. Add to list 

					// Do not bother saving for now...
					// only save segments...
					if(false) P.pihat.push_back(p);

					// And (A,B) pair to list
					P.pair1.push_back(i1);
					P.pair2.push_back(i2);

				}

			}

		if(!par::silent)
			cout << "\n";

		/////////////////////////////////
		//  Make list of unique individuals

		// copy first set of individuals
		P.in_anal = P.pair1;
		for(unsigned int ind = 0; ind < P.pair2.size(); ind++)
			P.in_anal.push_back(P.pair2[ind]);
		sort(P.in_anal.begin(), P.in_anal.end());
		vector<int>::iterator new_end =
						unique(P.in_anal.begin(), P.in_anal.end());
		// delete all elements past new_end 
		P.in_anal.erase(new_end, P.in_anal.end());

		P.printLOG(int2str(P.in_anal.size()) +
						" unique, informative individuals in analysis\n");

		if(P.in_anal.size() == 0) {
			error("No individuals left in analysis: halting");
		}

		/////////////////////////////////
		//  Verbose output: summarise IBD

		if(par::segment_output) {
			// P.summaryIBDsegments(perm);
		} else
			if(par::summary_ibd_output)
			P.summaryIBD();

		/////////////////////////////////
		//  Next chromosome
		par::done_global_pihat = true;

		if(!par::silent)
			cout << "\n";

	}

	// Now do IBD segment (as IBS...)
	if(par::segment_output) {
		P.summaryIBSsegments(perm);
		P.indivSegmentSummary();
	}

	//////////////////////////////
	// Find overlap in segments?
	if(par::segment_overlap)
		P.summariseHomoRuns();

	//////////////////////////////////////
	// Shut segment and multipoint files
	if(par::segment_output)
		SEG.close();
	if(par::multi_output)
		MP.close();
	if(par::gmulti_output)
		GMULTI.close();

	////////////////////////////////
	//  Output genome-wide p-values
	if(par::permute)
		if(chrs.size() >= 1 && (!par::ignore_phenotypes))
			P.displayGenomePV();

	////////////////////////////////
	//  We're definitely done now
	shutdown();
}
