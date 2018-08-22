/* 
 * File:   DcVar.h
 * Author: bwhite 10/19/17
 */

#ifndef DCVAR_H
#define	DCVAR_H

#include <string>
#include <map>
#include <vector>

#include <armadillo>

#include "plink.h"
// PLINK interface pointer PP defined in Insilico.h
#include "Insilico.h"
// bcw - 1/3/18/ - for CoordinateTable, TransciptFactorTable
#include "EpistasisEQtl.h"

// handle both PLINK (Caleb, et al paper) BED/BIM/BAM and 
// OMRF (Courtney Montgomery) separate files
enum SNP_INPUT_TYPE {
  SNP_SRC_PLINK, SNP_SRC_FILE
};

// SNPs/genotypes
// ID	#CHROM	POS	REF	REAL.REF
enum SNP_EXTRACT_FIELD_IDX {
  SNP_ID=0, SNP_CHROM=1, SNP_POS=2, SNP_REF_ALLELE=3, SNP_DBSNP_ALLELE=4
};
struct SNP_INFO {
  std::string chrom;
  uint position;
  std::string rsnum;
  char refAllele;
};
typedef std::vector<SNP_INFO> SNP_INFO_LIST;
typedef std::vector<SNP_INFO>::const_iterator SNP_INFO_LIST_IT;

// ChIP-Seq histone modification site record fields of interest
// ChIP-seq expression = histone modification site reads
// chromosome	position	null	full	likelihood_ratio	p_value	alpha_dispersion
// 0          1         2     3     4                 5       6
// beta_dispersioneta	r_dispersion	total.allele.specific	ref.allele.specific
// 7                  8             9                     10
// alt.allele.specific	total_region_reads	negLog10_p_value	BH_FDR	rsID
// 11                   12                  13                14      15
enum CHIP_SEQ_EXTRACT_FIELD_IDX {
  CHIP_SEQ_CHROM=0, CHIP_SEQ_POS=1, CHIP_SEQ_EXPR=12, CHIP_SEQ_SNP=15
};
struct CHIP_SEQ_INFO {
  std::string chrom;
  uint position;
  std::string rsnum;
  uint totalRegionReads;
};
typedef std::vector<CHIP_SEQ_INFO> CHIP_SEQ_INFO_LIST;

const string CHECKPOINT_FILENAME = "dcvar.chk";

class DcVar {
public:
  DcVar(SNP_INPUT_TYPE snpInputTypeParam=SNP_SRC_PLINK);
  bool Run();
  bool ReadTranscriptCoordinates(std::string coordinatesFile);
  bool ReadTranscriptFactorCoordinates(std::string coordinatesFile);
  bool SetRadius(int newRadius);
  int GetRadius() { return radius; }
  bool SetLocalCis(bool localCisFlag);
  bool GetLocalCis() { return localCis; }
  // ripped 1/7/18 from EpistasisEQtl for transcription factor considerations
  bool SetTFRadius(int newRadius);
  int GetTFRadius() { return tfRadius; }
  bool SetTF(bool tfFlag);
  bool GetTF() { return tfMode; }
  bool GetTFInfo(std::string tf, std::vector<int>& tfInfo);
  virtual ~DcVar();
private:
  bool RunPlink();
  bool RunOMRF();
  bool RunOMRFChipSeq();
  void PrintState();
  bool ReadGenotypesFile(uint chrom);
  bool ReadSnpLocationsFile(uint chrom);
  bool ReadGeneExpressionFile();
  bool ReadChipSeqFile();
  bool FindSnps(uint pos, std::vector<uint>& inRadius);
  bool MapPhenosToModel(std::vector<uint> phenos, std::string varModel,
                        std::vector<uint>& mappedPhenos);
  std::pair<uint, uint> MapSnpIndexToPlinkPhenos(uint snpIndex, 
                                                 std::string varModel);
  bool SplitExpressionCaseControl(arma::mat& caseMatrix, arma::mat& ctrlMatrix);
  bool ComputeDifferentialCorrelationZsparse(std::string snp, 
                                            arma::mat& cases, 
                                            arma::mat& ctrls);
  bool ComputeDifferentialCorrelationZ(std::string snp, 
                                       arma::mat& cases, 
                                       arma::mat& ctrls, 
                                       double correctedP);
  bool ComputeDifferentialCorrelationZvals(std::string snp, 
                                           arma::mat& cases, 
                                           arma::mat& ctrls);
  bool FlattenPvals(vector_t& retPvals);
  bool FilterPvalues(uint& numFiltered);
  double CalculateFdrBHThreshold();
  bool WriteCheckpoint(uint snpIndex, std::string snpName);
  bool ReadCheckpoint(std::pair<uint, string>& lastSnp);
  bool WriteResults(std::string filename, std::string curSnp);
  // INPUTS
  // ChIP-Seq
  bool chipSeqMode;
  CHIP_SEQ_INFO_LIST chipSeqExpression;
  // SNP meta data
  SNP_INPUT_TYPE snpInputType;
  std::vector<std::string> snpNames;
  SNP_INFO_LIST snpLocations;
  // SNP genotypes
  std::vector<std::string> genotypeSubjects;
  matrix_t genotypeMatrix;
  // RNA-Seq
  std::vector<std::string> geneExprNames;
  std::vector<std::string> geneExprSubjects;
  matrix_t expressionMatrix;
  // ALGORITHM
  double numCombs;
  uint totalTests;
  std::vector<uint> caseIdxCol;
  std::vector<uint> ctrlIdxCol;
  // OUTPUTS
  arma::mat zVals;
  arma::mat pVals;
  // added from EpistasisEQtl (iQTL) - bcw - 1/3/18
  bool CheckInputs();
  bool GetSnpsForTranscript(std::string transcript, 
    std::vector<int>& snpIndices);
  bool GetSnpsForTFs(std::vector<int>& snpIndices, std::vector<std::string>& tfs);
  bool LoadDefaultTranscriptionFactorLUT();
  bool IsSnpInTFs(int chr, int bp, std::string& tf);
  // search radius in kilobases
  uint radius;
  bool localCis;
  CoordinateTable coordinates;
  // added 4/21/15
  bool tfTableLoaded;
  bool tfMode;
  int tfRadius;
  TransciptFactorTable transcriptFactorLUT;
  std::vector<int> thisTranscriptSnpIndices;
  uint nOuterLoop;
  uint nInnerLoop;
  std::vector<int> thisTFSnpIndices;
};

#endif	/* DCVAR_H */
