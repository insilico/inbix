/* 
 * File:   EpistasisEQtl.h
 * Author: bwhite
 *
 * Created on October 3, 2013, 11:48 AM
 */

#ifndef EPISTASISEQTL_H
#define	EPISTASISEQTL_H

#include <vector>
#include <string>
#include <map>

typedef std::map<std::string, std::vector<int> > CoordinateTable;
typedef std::map<std::string, std::vector<int> >::const_iterator CoordinateTableCIt;

class EpistasisEQtl {
public:
  EpistasisEQtl();
  virtual ~EpistasisEQtl();
  bool ReadTranscriptCoordinates(std::string coordinatesFile);
  bool SetRadius(int newRadius);
  int GetRadius() { return radius; }
  bool SetLocalCis(bool localCisFlag);
  bool GetLocalCis() { return localCis; }
  bool Run();
private:
  int radius;
  bool localCis;
  CoordinateTable coordinates;
};

#endif	/* EPISTASISEQTL_H */

