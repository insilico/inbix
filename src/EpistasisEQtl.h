/* 
 * File:   EpistasisEQtl.h
 * Author: bwhite
 *
 * Created on October 3, 2013, 11:48 AM
 */

#ifndef EPISTASISEQTL_H
#define	EPISTASISEQTL_H

#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>

#include <armadillo>

typedef std::map<std::string, std::vector<uint> > CoordinateTable;
typedef std::map<std::string, std::vector<uint> >::const_iterator CoordinateTableCIt;

typedef std::map<std::string, std::vector<uint> > TransciptFactorTable;
typedef std::map<std::string, std::vector<uint> >::const_iterator TranscriptFactorTableCIt;

enum COORD_FIELDS {
  COORD_CHROM, COORD_BP_START, COORD_BP_END
};

class EpistasisEQtl {
public:
  EpistasisEQtl();
  virtual ~EpistasisEQtl();
  void PrintState();
  bool Run(bool debug=false);
  bool RunEqtl(std::string transcript);
  bool RunIqtlFull();
  bool RunIqtlCisTrans();
  bool ReadTranscriptCoordinates(std::string coordinatesFile);
  bool ReadTranscriptFactorCoordinates(std::string coordinatesFile);
  bool SetDebugMode(bool debugFlag=true);
  bool SetRadius(uint newRadius);
  uint GetRadius() { return radius; }
  bool SetLocalCis(bool localCisFlag);
  bool GetLocalCis() { return localCis; }
  // added 4/21/15 for transcription factor considerations
  bool SetTFRadius(uint newRadius);
  uint GetTFRadius() { return tfRadius; }
  bool SetTF(bool tfFlag);
  bool GetTF() { return tfMode; }
  bool GetTFInfo(std::string tf, std::vector<uint>& tfInfo);
private:
  bool CheckInputs();
  bool GetSnpsForTranscript(std::string transcript, 
    std::vector<uint>& snpIndices);
  bool GetSnpsForTFs(std::vector<uint>& snpIndices, std::vector<std::string>& tfs);
  bool LoadDefaultTranscriptionFactorLUT();
  bool IsSnpInTFs(uint chr, uint bp, std::string& tf);
  bool debugMode;
  uint radius;
  bool localCis;
  CoordinateTable coordinates;
  // added 4/21/15
  bool tfTableLoaded;
  bool tfMode;
  uint tfRadius;
  TransciptFactorTable transcriptFactorLUT;
  std::vector<uint> thisTranscriptSnpIndices;
  uint nOuterLoop;
  uint nInnerLoop;
  std::vector<uint> thisTFSnpIndices;
  uint goodModels;
  uint badModels;
  std::set<uint> outerLoopSnps;
  arma::mat resultsMatrixBetas;
  arma::mat resultsMatrixPvals;
};

#endif	/* EPISTASISEQTL_H */

