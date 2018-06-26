/**
 * \class PlinkInternalsDatasetInstance
 *
 * \brief Class to hold dataset instances (rows of attributes).
 *
 * Reworked entirely for McKinney Lab work - 2/28/11
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 6/14/05
 * Modified for inclusion in inbix 8/3/16
*/

#ifndef PLINK_INTERNALS_DATASET_INSTANCE_H
#define PLINK_INTERNALS_DATASET_INSTANCE_H

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "plink.h"
#include "Insilico.h"

/// forward reference to avoid circular include problems
class Dataset;
class DatasetInstance;

class PlinkInternalsDatasetInstance : public DatasetInstance
{
public:
  /*************************************************************************//**
   * Construct an data set instance object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  PlinkInternalsDatasetInstance(Dataset* ds, std::string instanceID,
                                Plink* plinkPtr, Individual* plinkInd);
  ~PlinkInternalsDatasetInstance();
  double GetSimpleSNPValue(int snp);
  /// return the number of discrete attributes
  unsigned int NumAttributes() override;
  /*************************************************************************//**
   * Get and return an attribute value at index.
   * \param [in] index attribute index
   * \return attribute value at index
   ****************************************************************************/
  AttributeLevel GetAttribute(unsigned int index) override;
  /// return the number of continuous attributes
  unsigned int NumNumerics() override;
  /*************************************************************************//**
   * Get and return numeric value at index.
   * \param [in] index numeric index
   * \return numeric value at index
   ****************************************************************************/
  NumericLevel GetNumeric(unsigned int index) override;
  void Print() override;
private:
  Plink* plinkInternalsPtr;
  string ID;
  Individual* individual;
};

#endif
