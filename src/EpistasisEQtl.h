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

typedef std::map<std::string, std::vector<int> > CoordinateTable;
typedef std::map<std::string, std::vector<int> >::const_iterator CoordinateTableCIt;

typedef std::map<std::string, std::vector<int> > TransciptFactorTable;
typedef std::map<std::string, std::vector<int> >::const_iterator TranscriptFactorTableCIt;

enum COORD_FIELDS {
  COORD_CHROM, COORD_BP_START, COORD_BP_END
};

class EpistasisEQtl {
public:
  EpistasisEQtl();
  virtual ~EpistasisEQtl();
  bool ReadTranscriptCoordinates(std::string coordinatesFile);
  bool ReadTranscriptFactorCoordinates(std::string coordinatesFile);
  bool SetRadius(int newRadius);
  int GetRadius() { return radius; }
  bool SetLocalCis(bool localCisFlag);
  bool GetLocalCis() { return localCis; }
  bool Run(bool debug=false);
  // added 4/21/15
  bool SetTFRadius(int newRadius);
  int GetTFRadius() { return tfRadius; }
  bool SetTF(bool tfFlag);
  bool GetTF() { return tfMode; }
  bool GetTFInfo(std::string tf, std::vector<int>& tfInfo);
private:
  bool GetSnpsForTranscript(std::string transcript, 
    std::vector<int>& snpIndices);
  bool GetSnpsForTFs(std::vector<int>& snpIndices, std::vector<std::string>& tfs);
  bool LoadDefaultTranscriptionFactorLUT();
  bool IsSnpInTFs(int chr, int bp, std::string& tf);
  int radius;
  bool localCis;
  CoordinateTable coordinates;
  // added 4/21/15
  bool tfMode;
  int tfRadius;
  TransciptFactorTable transcriptFactorLUT;
};

#endif	/* EPISTASISEQTL_H */

