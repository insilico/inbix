/* =============================================================================
 *
 * Filename: regain.cpp - Bill White - 4/23/13
 *
 * Description:  Regression GAIN calculation
 *
 * Created:  06/20/2011
 * Original Author:  Nick Davis, nick-davis@utulsa.edu
 * =============================================================================
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iterator>

#include "gsl/gsl_cdf.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "plink.h"
#include "options.h"
#include "model.h"
#include "logistic.h"
#include "linear.h"
#include "stats.h"
#include "helper.h"

#include "regain.h"

// Plink object
extern Plink* PP;

Regain::Regain(bool compr, double sifthr, bool integrative, bool compo, 
				bool fdrpr) {
	// set class vars to passed args
	compressed = compr;
	sif_thresh = sifthr;
	intregain = integrative;
	component = compo;
	fdrprune = fdrpr;

	// set integrative/normal regain vars
	// additional ext for integrative
	string ext = intregain ? ".block" : "";
	// header in betas files
	string hdr = intregain ? "attr" : "SNP";
	
	// total number of attributes
	numattr = intregain ? PP->nl_all + PP->nlistname.size() : PP->nl_all;
	PP->printLOG("Total number of attributes: " + int2str(numattr) + "\n");
	
	// initialize matrices and open output files
	string beta_f = par::output_file_name + ext + ".betas";
	string mebeta_f = par::output_file_name + ext + ".mebetas";
	BETAS.open(beta_f.c_str(), ios::out);
	MEBETAS.open(mebeta_f.c_str(), ios::out);
	PP->printLOG("Writing interaction beta values to [ " + beta_f + " ]\n");
	PP->printLOG("Writing main effect beta values to [ " + mebeta_f + " ]\n");
	BETAS.precision(6);
	MEBETAS.precision(6);
	// print header
	BETAS << hdr << "1\t" << hdr << "2\tB_0\tB_1\tB_1 P-VAL\tB_2\tB_2 P-VAL";
	if (par::covar_file) {
		for (int i = 0; i < par::clist_number; i++) {
			BETAS << "\t" << PP->clistname[i] << "\t" << PP->clistname[i] << " P-VAL";
		}
	}
	BETAS << "\tB_3\tB_3 P-VAL" << endl;

	MEBETAS << hdr << "\tB_0\tB_1\tB_1 P-VAL";
	if (par::covar_file) {
		for (int i = 0; i < par::clist_number; i++) {
			MEBETAS << "\t" << PP->clistname[i] << "\t" << PP->clistname[i] 
							<< " P-VAL";
		}
	}
	MEBETAS << endl;

	string sif_f = par::output_file_name + ext + ".sif";
	SIF.open(sif_f.c_str(), ios::out);
	PP->printLOG("Writing Cytoscape network file (SIF) to [ " + sif_f + " ]\n");
	SIF.precision(6);
	if (component) {
		string snp_sif_f = par::output_file_name + ".snp.sif";
		string num_sif_f = par::output_file_name + ".num.sif";
		string int_sif_f = par::output_file_name + ".int.sif";
		SNP_SIF.open(snp_sif_f.c_str(), ios::out);
		PP->printLOG("Writing SNP Cytoscape network file (SIF) to [ " + 
		snp_sif_f + " ]\n");
		SNP_SIF.precision(6);
		NUM_SIF.open(num_sif_f.c_str(), ios::out);
		PP->printLOG("Writing numeric Cytoscape network file (SIF) to [ " + 
		num_sif_f + " ]\n");
		NUM_SIF.precision(6);
		INT_SIF.open(int_sif_f.c_str(), ios::out);
		PP->printLOG("Writing integrative Cytoscape network file (SIF) to [ " + 
		int_sif_f + " ]\n");
		INT_SIF.precision(6);
	}

	regainMatrix = new double*[numattr];
	regainPMatrix = new double*[numattr];
	// allocate reGAIN matrix
	for (int i = 0; i < numattr; ++i) {
		regainMatrix[i] = new double[numattr];
	}
	// allocate reGAIN p-value matrix
	for (int i = 0; i < numattr; ++i) {
		regainPMatrix[i] = new double[numattr];
	}

}

Regain::~Regain() {
	// close BETAS and SIF ofstreams
	BETAS.close();
	SIF.close();

	// free regain matrix memory
	for (int i = 0; i < numattr; ++i) {
		delete [] regainMatrix[i];
	}
	delete [] regainMatrix;

	// free regain p-value matrix memory
	for (int i = 0; i < numattr; ++i) {
		delete [] regainPMatrix[i];
	}
	delete [] regainPMatrix;
}

void Regain::run() {
	int e1, e2;
#ifdef _OPENMP
	// OpenMP parallelization of this outer loop
	int numThreads = omp_get_num_threads();
	int numProcs = omp_get_num_procs();
	cout << "\t\t" << numThreads << " OpenMP threads available" << endl;
	cout << "\t\t" << numProcs << " OpenMP processors available" << endl;
#pragma omp parallel for schedule(dynamic, 1) private(e1, e2)
#endif
	for (e1 = 0; e1 < numattr; e1++) {
		for (e2 = 0; e2 < numattr; e2++) {
			// We've already performed this test, since the matrix is symmetric
			if (e1 > e2) continue;

			// main effect of SNP/numeric attribute 1 - diagonal of the reGAIN matrix
			if (e1 == e2) {
#ifdef _OPENMP
#pragma omp critical
#endif
				mainEffect(e1, e1 >= PP->nl_all);
			}	else {
				interactionEffect(e1, e1 >= PP->nl_all, e2, e2 >= PP->nl_all);
			}
		}
	} // Next pair of SNPs/numeric attributes
}

void Regain::mainEffect(int e1, bool numeric) {
	Model *lm_main_effect;

	// logisic regression for binary phenotypes (traits), linear otherwise
	if (par::bt) {
		LogisticModel * m = new LogisticModel(PP);
		lm_main_effect = m;
	} else {
		LinearModel * m = new LinearModel(PP);
		lm_main_effect = m;
	}

	// Set missing data
	lm_main_effect->setMissing();

	// label in regression model
	string label = numeric ? PP->nlistname[e1 - PP->nl_all] : PP->locus[e1]->name;

	// Main effect of SNP/numeric attribute
	if (numeric) {
		lm_main_effect->addNumeric(e1 - PP->nl_all);
	}
	else {
		lm_main_effect->addAdditiveSNP(e1);
	}
	lm_main_effect->label.push_back(label);

	// add covariates if specified
	if (par::covar_file) {
		addCovariates(*lm_main_effect);
	}

	// Build design matrix
	lm_main_effect->buildDesignMatrix();

	// Fit linear model
	lm_main_effect->fitLM();

	// Did model fit okay?
	if(!lm_main_effect->isValid()) {
		PP->printLOG("\nERROR: Invalid main effect regression fit variable " + 
		label + ")\n");
		shutdown();
	}

	lm_main_effect->validParameters();

	// Obtain estimates and statistic
	int tp = 1; // always use first coefficient after intercept as main effect term
	lm_main_effect->testParameter = tp; // single variable main effect
	vector_t b_main_effect = lm_main_effect->getCoefs();
	vector_t b_p_values = lm_main_effect->getPVals();

	// set main effect (diagonal) beta coefficient and corresponding p-value
	regainMatrix[e1][e1] = b_main_effect[tp];
	// p-values don't include intercept term
	regainPMatrix[e1][e1] = b_p_values[tp - 1];

	// update main effect betas file
	if (numeric) {
		MEBETAS << PP->nlistname[e1 - PP->nl_all];
	}
	else {
		MEBETAS << PP->locus[e1]->name;
	}
	for (unsigned int i = 0; i < b_main_effect.size(); ++i) {
		// B0 coefficient doesn't have pval
		if (i == 0) {
			MEBETAS << "\t" << b_main_effect[i];
		}
		else {
			// adjust pvals index since there's no B0 pval
			MEBETAS << "\t" << b_main_effect[i] << "\t" << b_p_values[i - 1];
		}
	}
	MEBETAS << endl;

	// free model memory
	delete lm_main_effect;
}

void Regain::addCovariates(Model &m) {
	for (int i = 0; i < par::clist_number; i++) {
		// add covariate to the model
		m.addCovariate(i);
		m.label.push_back(PP->clistname[i]);
	}
}

void Regain::interactionEffect(int e1, bool numeric1, int e2, bool numeric2) {
	Model * lm;

	// logistic regression for binary phenotypes (traits), linear otherwise
	if (par::bt) {
		LogisticModel * m = new LogisticModel(PP);
		lm = m;
	} else {
		LinearModel * m = new LinearModel(PP);
		lm = m;
	}

	// Set missing data
	lm->setMissing();

	// labels in regression model
	string label1 = numeric1? PP->nlistname[e1 - PP->nl_all]: PP->locus[e1]->name;
	string label2 = numeric2? PP->nlistname[e2 - PP->nl_all]: PP->locus[e2]->name;

	// Main effect of SNP/numeric attribute 1
	if (numeric1) {
		lm->addNumeric(e1 - PP->nl_all);
	}	else {
		lm->addAdditiveSNP(e1);
	}
	lm->label.push_back(label1);

	// Main effect of SNP/numeric attribute 2
	if (numeric2) {
		lm->addNumeric(e2 - PP->nl_all);
	}
	else {
		lm->addAdditiveSNP(e2);
	}
	lm->label.push_back(label2);

	// add covariates if specified
	if (par::covar_file) addCovariates(*lm);

	// interaction
	lm->addInteraction(1, 2);
	lm->label.push_back("EPI");

	// Build design matrix
	lm->buildDesignMatrix();

	// Fit linear model
	lm->fitLM();
	
	// Did model fit okay?
//	if(!lm->isValid()) {
//		PP->printLOG("\nWARNING: Invalid regression fit for interaction variables (" + 
//						label1 + ", " + label2 + ")\n");
//		vector<bool> vp = lm->validParameters();
//		cout << "valid parameters?" << endl;
//		copy(vp.begin(), vp.end(), ostream_iterator<bool>(cout, "\n"));
//		cout << endl;
//		// shutdown();
//	}
	
	// interaction
	// TODO: change this when considering main effects in model or not
	int tp = 3;

	// add # covars to test param to get interaction param
	if (par::covar_file) tp += par::clist_number;

	lm->testParameter = tp; // interaction

#ifdef _OPENMP
#pragma omp critical
	{
#endif
		vector_t b = lm->getCoefs();
//		cout << "beta vector:" << endl;
//		copy(b.begin(), b.end(), ostream_iterator<double>(cout, "\n"));

		vector_t beta_p = lm->getPVals();
		//	cout << "pvals vector:" << endl;
		//	copy(beta_p.begin(), beta_p.end(), ostream_iterator<double>(cout, "\n"));

		regainMatrix[e1][e2] = b[b.size() - 1];
		regainMatrix[e2][e1] = b[b.size() - 1];
		regainPMatrix[e1][e2] = beta_p[beta_p.size() - 1];
		regainPMatrix[e2][e1] = beta_p[beta_p.size() - 1];

		double lm_stat = lm->getStatistic();
		vector_t lm_se = lm->getSE();
		//	cout << "SE vector:" << endl;
		//	copy(lm_se.begin(), lm_se.end(), ostream_iterator<double>(cout, "\n"));

		// create a new value for the reGAIN matrix from beta/SE
		vector_t::const_iterator bIt = b.begin();
		vector_t::const_iterator sIt = lm_se.begin();
		vector_t tTestValues;
		for (; bIt != b.end(); ++bIt, ++sIt) {
			tTestValues.push_back(*bIt / *sIt);
		}
		//	cout << "new betas:" << endl;
		//	copy(tTestValues.begin(), tTestValues.end(), 
		//					ostream_iterator<double>(cout, "\n"));

//		double interactionValue = tTestValues[b.size() - 1];
//		regainMatrix[e1][e2] = interactionValue;
//		regainMatrix[e2][e1] = interactionValue;
		// degrees or freedom = numSamples - 3
//		unsigned int numSamples = PP->sample.size();
		// cout << numSamples << " samples for " << (numSamples-3) 
		// << " degrees of freedom" << endl;
//		double gslPval = gsl_cdf_tdist_P(abs(interactionValue), numSamples - 3);
//		double pValue = 2 * (1 - gslPval);
//		regainPMatrix[e1][e2] = pValue;
//		regainPMatrix[e2][e1] = pValue;

		//	cout << "DEBUG info (" << e1 << ", " << e2 << ")" 
		//					<< ", beta = " << setw(10) << b[b.size() - 1] 
		//					<< ", p = " << setw(10) << beta_p[beta_p.size() - 1] 
		//					<< ", x^2 = " << setw(10) << lm_stat
		//					<< ", is valid? " << (lm->isValid()? "TRUE": "FALSE")
		//					<< ", VIF check? " << (lm->checkVIF()? "TRUE": "FALSE")
		//					<< " => New interaction value: " << setw(10) << tTestValues[3] 
		//					<< " p = " << setw(10) << pValue 
		//					<< endl;
//		cout << setprecision(3);
//		cout << fixed;
//		cout << (lm->isValid() ? "TRUE" : "FALSE")
//						<< "\t(" << label1 << ", " << label2 << ")"
//						<< ", beta = " << setw(12) << b[3]
//						<< ", p = " << setw(6) << beta_p[2]
//						<< ", S.E. = " << setw(12) << lm_se[3]
//						<< "\t => " << setw(6) << tTestValues[3]
//						<< " p = " << setw(6) << pValue
//						<< endl;
		cout << (lm->fitConverged() ? "TRUE" : "FALSE")
						<< "\t" << label1 << "\t" << label2 
						<< "\t" << setw(12) << b[3]
						<< "\t" << setw(6) << beta_p[2]
						<< "\t" << setw(12) << lm_se[3]
						<< "\t" << setw(12) << b[3] / lm_se[3]
						<< "\t" << setw(12) << lm_stat
						<< endl;

		// store p-value along with (e1, e2) location of
		// item.  This is used later for FDR pruning
		if (fdrprune) {
			pair<int, int> p = make_pair(e1, e2);
			mat_el Pint = make_pair(beta_p[beta_p.size() - 1], p);
			// mat_el Pint = make_pair(pValue, p);
			gainPint.push_back(Pint);
		}

		// update BETAS and SIF files

		// numeric attributes or SNPs
		if (numeric1) {
			BETAS << PP->nlistname[e1 - PP->nl_all] << "\t";
		}
		else {
			BETAS << PP->locus[e1]->name << "\t";
		}
		if (numeric2) {
			BETAS << PP->nlistname[e2 - PP->nl_all];
		}
		else {
			BETAS << PP->locus[e2]->name;
		}

		for (unsigned int i = 0; i < b.size(); ++i) {
			// B0 coefficient doesn't have pval
			if (i == 0) {
				BETAS << "\t" << b[i];
			}
			else {
				// adjust pvals index since there's no B0 pval
				BETAS << "\t" << b[i] << "\t" << beta_p[i - 1];
			}
		}
		BETAS << endl;

		// add to SIF if interaction >= threshold
		if (abs(b[b.size() - 1]) >= sif_thresh) {
			if (numeric1) {
				SIF << PP->nlistname[e1 - PP->nl_all] << "\t" << abs(b[b.size() - 1]) 
								<< "\t";
			}
			else {
				SIF << PP->locus[e1]->name << "\t" << abs(b[b.size() - 1]) << "\t";
			}
			if (numeric2) {
				SIF << PP->nlistname[e2 - PP->nl_all] << endl;
			}
			else {
				SIF << PP->locus[e2]->name << endl;
			}

			if (component) {
				// numeric
				if (numeric1 && numeric2) {
					NUM_SIF << PP->nlistname[e1 - PP->nl_all] << "\t" 
									<< abs(b[b.size() - 1]) << "\t";
					NUM_SIF << PP->nlistname[e2 - PP->nl_all] << endl;
				}					// integrative
				else if (numeric1 && !numeric2) {
					INT_SIF << PP->nlistname[e1 - PP->nl_all] << "\t" 
									<< abs(b[b.size() - 1]) << "\t";
					INT_SIF << PP->locus[e2]->name << endl;
				}					// integrative
				else if (!numeric1 && numeric2) {
					INT_SIF << PP->locus[e1]->name << "\t" << abs(b[b.size() - 1]) << "\t";
					INT_SIF << PP->nlistname[e2 - PP->nl_all] << endl;
				}					// SNP
				else {
					SNP_SIF << PP->locus[e1]->name << "\t" << abs(b[b.size() - 1]) << "\t";
					SNP_SIF << PP->locus[e2]->name << endl;
				}
			}

		}
#ifdef _OPENMP
	}
#endif
	// free model memory
	delete lm;
}

void Regain::writeRegain(bool pvals, bool fdrprune) {
	
	double** regainMat;
	if (pvals) {
		regainMat = regainPMatrix;
	}	else {
		regainMat = regainMatrix;
	}
	
	// write the reGAIN matrix to file named <output_file_name>.regain
	string snp_f = par::output_file_name;
	string num_f = par::output_file_name;
	string int_f = par::output_file_name; 
	string regain_matrix_f = par::output_file_name;
	
	// additional prefixes/extension for output filename 
	// FDR-pruned
	string prnpre = fdrprune ? ".pruned" : "";
	// p-values file
	string pvpre = pvals ? ".pvals" : "";
	// integrative
	string intpre = intregain ? ".block" : "";
	// compressed/binary file
	string tail = compressed ? ".gz" : "";

	// additional output text	
	string pvtext = pvals ? "p-value " : "";
	string fdrtext = fdrprune ? "FDR-pruned " : "";

	regain_matrix_f += intpre + pvpre + prnpre + ".regain" + tail;

	PP->printLOG("Writing " + fdrtext + " REGAIN " + pvtext + 
	"matrix [ " + regain_matrix_f + " ]\n");
	REGAIN_MATRIX.open(regain_matrix_f.c_str(), compressed);
	if (component) {
		snp_f += ".snp" + pvpre + prnpre + ".regain" + tail;
		PP->printLOG("Writing " + fdrtext + "SNP REGAIN " + pvtext + 
		"matrix [ " + snp_f + " ]\n");
		SNP_MATRIX.open(snp_f.c_str(), compressed);

		num_f += ".num" + pvpre + prnpre + ".regain" + tail;
		PP->printLOG("Writing " + fdrtext + "numeric REGAIN " + pvtext + 
		"matrix [ " + num_f + " ]\n");
		NUM_MATRIX.open(num_f.c_str(), compressed);

		int_f += ".int" + pvpre + prnpre + ".regain" + tail;
		PP->printLOG("Writing " + fdrtext + "integrative REGAIN " + 
		pvtext + "matrix [ " + int_f + " ]\n");
		INT_MATRIX.open(int_f.c_str(), compressed);
	}
	// write SNP column names
	for (int cn = 0; cn < PP->nl_all; ++cn) {
		if (cn) {
			REGAIN_MATRIX << "\t" << PP->locus[cn]->name;
			if (component) SNP_MATRIX << "\t" << PP->locus[cn]->name;
		} else {
			REGAIN_MATRIX << PP->locus[cn]->name;
			if (component) SNP_MATRIX << PP->locus[cn]->name;
		}
	}
	// write numeric attribute column names
	for (int cn = 0; cn < PP->nlistname.size(); ++cn) {
		if (!cn && !PP->nl_all) REGAIN_MATRIX << PP->nlistname[cn];
		else REGAIN_MATRIX << "\t" << PP->nlistname[cn];
		if (component) {
			if (cn) {
				NUM_MATRIX << "\t" << PP->nlistname[cn];
				INT_MATRIX << "\t" << PP->nlistname[cn];
			} else {
				NUM_MATRIX << PP->nlistname[cn];
				INT_MATRIX << PP->nlistname[cn];
			}
		}
	}
	REGAIN_MATRIX << "\n";
	if (component) {
		NUM_MATRIX << "\n";
		INT_MATRIX << "\n";
		SNP_MATRIX << "\n";
	}
	// write matrix entries
	for (int i = 0; i < numattr; ++i) {
		for (int j = i; j < numattr; ++j) {
			// use absolute value (magnitude) of betas 
			regainMat[i][j] = abs(regainMat[i][j]);
			// regainMat[i][j] = fabs(regainMat[i][j]);
			if (j == i) {// fill in symmetric entries, replacing j < i with tabs
				string tabs = "";
				for (int k = 0; k < j; k++)
					tabs += "\t";
				REGAIN_MATRIX << tabs << dbl2str_fixed(regainMat[i][j], 6);
				if (component) {
					if (i < PP->nl_all) 
						SNP_MATRIX << tabs << dbl2str_fixed(regainMat[i][j], 6);
					else {
						tabs = "";
						for (int k = PP->nl_all; k < j; k++)
							tabs += "\t";
						NUM_MATRIX << tabs << dbl2str_fixed(regainMat[i][j], 6);
					}
				}
			} else {
				REGAIN_MATRIX << "\t" << dbl2str_fixed(regainMat[i][j], 6);
				if (component) {
					if (i < PP->nl_all) {
						if (j < PP->nl_all) 
							SNP_MATRIX << "\t" << dbl2str_fixed(regainMat[i][j], 6);
						else {
							if (j == PP->nl_all)
								INT_MATRIX << dbl2str_fixed(regainMat[i][j], 6);
							else INT_MATRIX << "\t" << dbl2str_fixed(regainMat[i][j], 6);
						}
					} else NUM_MATRIX << "\t" << dbl2str_fixed(regainMat[i][j], 6);
				}
			}
		}
		REGAIN_MATRIX << "\n";
		if (component) {
			if (i < PP->nl_all) {
				SNP_MATRIX << "\n";
				INT_MATRIX << "\n";
			} else NUM_MATRIX << "\n";
		}
	}

	// close output stream
	REGAIN_MATRIX.close();
	if (component) {
		SNP_MATRIX.close();
		NUM_MATRIX.close();
		INT_MATRIX.close();
	}
}

void Regain::fdrPrune(double fdr) {
	PP->printLOG("Calculating Benjamini Hochberg FDR for pruning\n");
	int m = gainPint.size();
	// sort gain interaction mal_el type by p-value, maintaining
	// gainPMatrix location (row, col) with sorted values
	sort(gainPint.begin(), gainPint.end(), Regain::mecomp);

	// use rough FDR (RFDR) to estimate alpha based on input FDR
	double alpha = 2 * m * fdr / (m + 1);
	int R = -1;
	// BH method
	for (int i = 0; i < m; i++) {
		double l = (i + 1) * alpha / (double) m;
		// test whether current p-value < current l
		if (gainPint[i].first < l) {
			R = i;
		} else break;
	}

	// BH threshold condition not met with any p-values, so exit
	if (R == -1) {
		PP->printLOG("No p-value meets BH threshold criteria, so nothing pruned\n");
		return;
	}

	// BH rejection threshold
	double T = gainPint[R].first;
	PP->printLOG("BH rejection threshold: T = " + dbl2str(T) + ", R = " + 
	int2str(R) + "\n");
	PP->printLOG("Pruning reGAIN interaction terms with p-values > T (" + 
	dbl2str(T) + ")\n");

	// now prune (set to 0.0) all values greater than R index
	for (int i = R + 1; i < m; i++) {
		pair<int, int> p = gainPint[i].second;
		// symmetric matrix, so set [e1][e2] and [e2][e1]
		regainMatrix[p.first][p.second] = 0.0;
		regainMatrix[p.second][p.first] = 0.0;
	}
	PP->printLOG("Pruned " + int2str(m - (R + 1)) + 
	" values from reGAIN interaction terms\n");
	// use threshold to write R commands to generate FDR plot 
	writeRcomm(T, fdr);
}

void Regain::writeRcomm(double T, double fdr) {
	ofstream RCOMM;
	RCOMM.precision(6);
	string fdr_r_file = par::output_file_name + ".R";
	string betas_file = par::output_file_name + ".betas";
	PP->printLOG("Writing R commands to generate FDR plot [" + fdr_r_file + "]\n");

	RCOMM.open(fdr_r_file.c_str(), ios::out);
	RCOMM << "fdrvars <- read.delim(\"" << betas_file << "\")" << endl;
	RCOMM << "library(calibrate)" << endl;
	RCOMM << "betas <- fdrvars$B_3" << endl;
	RCOMM << "pvals <- fdrvars$B_3.P.VAL" << endl;
	RCOMM << "betas <- abs(betas)" << endl;
	RCOMM << "T <- " << T << endl;
	RCOMM << "partition <- " << fdr << endl;
	RCOMM << "plot(betas, -log10(pvals), type=\"n\")" << endl;
	RCOMM << "abline(h=-log10(T), col=\"green4\", lwd=3)" << endl;
	RCOMM << "accept <- which(-log10(pvals) >= -log10(T))" << endl;
	RCOMM << "reject <- which(-log10(pvals) < -log10(T))" << endl;
	RCOMM << "prnidx <- partition * length(betas[accept])" << endl;
	RCOMM << "srtaccbetas <- sort(betas[accept])" << endl;
	RCOMM << "prnval <- srtaccbetas[prnidx]" << endl;
	RCOMM << "if(prnidx%%1!=0){" << endl;
	RCOMM << "prnval <- (srtaccbetas[floor(prnidx)] + srtaccbetas[ceiling(prnidx)]) / 2" << endl;
	RCOMM << "}" << endl;
	RCOMM << "prunex <- which(betas <= prnval)" << endl;
	RCOMM << "pruney <- which(-log10(pvals) >= -log10(T))" << endl;
	RCOMM << "prune <- intersect(prunex, pruney)" << endl;
	RCOMM << "accept <- setdiff(accept, prune)" << endl;
	RCOMM << "points(betas[accept], -log10(pvals[accept]), bg=\"green4\", pch=21)" << endl;
	RCOMM << "snp1 <- fdrvars$SNP1" << endl;
	RCOMM << "snp2 <- fdrvars$SNP2" << endl;
	RCOMM << "textxy(betas[accept], -log10(pvals[accept]), paste(snp1, snp2, sep=\",\")[accept])" << endl;
	RCOMM << "points(betas[reject], -log10(pvals[reject]), bg=\"blue\", pch=22)" << endl;
	RCOMM << "points(betas[prune], -log10(pvals[prune]), bg=\"red\", pch=24)" << endl;
	RCOMM << "abline(v=prnval, col=\"red\", lwd=3)" << endl;
	RCOMM << "title(\"Scatter plot of -log10 transformed p-values vs. regression betas\")" << endl;
	RCOMM << "legend(\"topleft\", inset=.05, title=\"Type\", c(\"Accepted\", \"Rejected\", \"Pruned\"), pch=c(21,22,24), pt.bg=c(\"green4\", \"blue\", \"red\"))" << endl;
	RCOMM.close();
}

bool Regain::mecomp(const mat_el &l, const mat_el &r) {
	return l.first < r.first;
}
