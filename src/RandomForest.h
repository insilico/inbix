/**
 * \class RandomForest
 *
 * \brief RandomForest attribute ranking algorithm.
 *
 * Adapter class to map EC call for Random Jungle importance scores
 * to Random Jungle library functions.
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 10/16/11
 */

#ifndef RANDOMFOREST_H
#define	RANDOMFOREST_H

#include <fstream>

#include "Insilico.h"
#include "Dataset.h"
#include "AttributeRanker.h"

#include "plink.h"
#include "Forest.h"

class RandomForest: public AttributeRanker
{
public:
  /*************************************************************************//**
   * Construct an RandomForest algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] plinkPtr pointer to a PLINK environment
   ****************************************************************************/
  RandomForest(Dataset* ds, Plink* plinkPtr);
  /*************************************************************************//**
   * Construct an RandomForest algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] bestAttributeNames best attribute names in the data set to use
   ****************************************************************************/
  RandomForest(Dataset* ds, std::vector<std::string> bestAttributeNames);
  /*************************************************************************//**
   * Deconstruct an RandomForest algorithm object.
   ****************************************************************************/
  virtual ~RandomForest();
  AttributeScores ComputeScores() override;
  double GetClassificationError() override;
  void WriteScores(std::string baseFilename) override;
  void WriteScoresInternal();
  // bcw 10/12/16
  bool InitializeData(bool useMask=false) override;
  double Predict(Dataset* testData);
private:
  bool CreateDefaultForestForPheno();
  Forest* forest;
  unsigned int minNodeSize;
  Plink* PP;
};

#endif	/* RANDOMFOREST_H */
