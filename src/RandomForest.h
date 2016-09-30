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
#include "Data.h"
#include "Forest.h"

class RandomForest: public AttributeRanker
{
public:
  /*************************************************************************//**
   * Construct an RandomForest algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] vm reference to a Boost map of command line options
   ****************************************************************************/
  RandomForest(Dataset* ds, Plink* plinkPtr);
  /*************************************************************************//**
   * Construct an RandomForest algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] configMap reference ConfigMap (map<string, string>)
   ****************************************************************************/
  virtual ~RandomForest();
  AttributeScores ComputeScores() override;
  double GetClassificationError() override;
  void WriteScores(std::string baseFilename) override;
  void WriteScoresInternal();
private:
  Data* data;
  Forest* forest;
  Plink* PP;
};

#endif	/* RANDOMFOREST_H */

