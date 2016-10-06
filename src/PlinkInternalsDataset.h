/**
 * \class PlinkInternalsDataset
 *
 * \brief Plink internals adapter class.
 *
 * \sa Dataset
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 2/24/11
 * 
 * Modified for inclusion in inbix 8/3/16
 */

#ifndef PLINKINTERNALSDATASET_H
#define	PLINKINTERNALSDATASET_H

#include <string>
#include <vector>

#include "plink.h"
#include "Dataset.h"

class PlinkInternalsDataset : public Dataset
{
public:
  PlinkInternalsDataset(Plink* plinkPtr);
  bool LoadDatasetPP();
  ~PlinkInternalsDataset();
private:
  Plink* PP;
};

#endif	/* PLINKINTERNALSDATASET_H */
