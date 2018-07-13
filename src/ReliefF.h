/**
 * \class ReliefF
 *
 * \brief ReliefF attribute ranking algorithm.
 *
 * Totally redone for the McKinney in silico lab in 2011.
 * Large refactoring to move all attribute elimination handling to the
 * Dataset and its subclasses. 9/11/11
 *
 * \sa RReliefF
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 7/16/05
 */

#ifndef RELIEFF_H
#define RELIEFF_H

#include <vector>
#include <fstream>

#include "plink.h"

#include "AttributeRanker.h"
#include "Dataset.h"
#include "Insilico.h"

class ReliefF : public AttributeRanker
{
public:
  /*************************************************************************//**
   * Construct an ReliefF algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] plinkPtr pointer to a PLINK object   
   * \param [in] anaType analysis type \see AnalysisType
   ****************************************************************************/
  ReliefF(Dataset* ds, Plink* plinkPtr, AnalysisType anaType);
  virtual ~ReliefF();
  /**
   * Compute the ReliefF scores for the current set of attributes.
   * Implements ReliefF algorithm:
   * Marko Robnik-Sikonja, Igor Kononenko: Theoretical and Empirical Analysis of
   * ReliefF and RReliefF. Machine Learning Journal, 53:23-69, 2003
   * http://lkm.fri.uni-lj.si/rmarko/papers/robnik03-mlj.pdf
   */
  virtual bool ComputeAttributeScores();
  /// Compute the ReliefF scores by iteratively removing worst attributes.
  bool ComputeAttributeScoresIteratively();
  /// Resets some data structures for the next iteration of ReliefF
  bool ResetForNextIteration() override;
  /*************************************************************************//**
   * Write the scores and attribute names to stream.
   * \param [in] outStream stream to write score-attribute name pairs
   ****************************************************************************/
  void PrintAttributeScores(std::ofstream& outStream);
  /*************************************************************************//**
   * Write the scores and attribute names to file.
   * \param [in] baseFilename filename to write score-attribute name pairs
   ****************************************************************************/
  void WriteAttributeScores(std::string baseFilename);
  /// Pre-compute all pairwise instance-to-instance distances.
  bool PreComputeDistances();
  /// Implements AttributeRanker interface.
  AttributeScores ComputeScores() override;
  void PrintBestKs();
   /*************************************************************************//**
	 * Write the best k-nearest neighbors best k and attribute names to file.
	 * \param [in] baseFilename filename to write best-k-attribute name pairs
	 ****************************************************************************/
  void WriteBestKs(std::string baseFilename);
   /// Compute scores based on optimum k
  bool ComputeAttributeScoresKopt();
  /// Print the ID to instance index map to stdout
  void PrintInstancesMask();
private:
  /// no default constructor
  ReliefF();
  
protected:
  // pointer to the PLINK ecosystem
  Plink* PP;
  // distance matrix
  matrix_t distanceMatrix;

  /// Compute the const AttributeScores& ComputeScores(); weight by distance factors for nearest neighbors.
  bool ComputeWeightByDistanceFactors();
  /// type of analysis to perform
  AnalysisType analysisType;
  /*************************************************************************//**
   * Compute the discrete difference in an attribute between two instances.
   * \param [in] attributeIndex index into vector of all attributes
   * \param [in] dsi1 pointer to DatasetInstance 1
   * \param [in] dsi2 pointer to DatasetInstance 2
   * \return diff(erence)
   ****************************************************************************/
  double (*snpDiffFuncPtr)(unsigned int attributeIndex,
                    DatasetInstance* dsi1,
                    DatasetInstance* dsi2);
  /*************************************************************************//**
   * Compute the continuous difference in an attribute between two instances.
   * \param [in] attributeIndex index into vector of all attributes
   * \param [in] dsi1 pointer to DatasetInstance 1
   * \param [in] dsi2 pointer to DatasetInstance 2
   * \return diff(erence)
   ****************************************************************************/
  double (*numDiffFuncPtr)(unsigned int attributeIndex,
                    DatasetInstance* dsi1,
                    DatasetInstance* dsi2);
  /// the name of discrete diff(erence) function
  std::string snpDiffMetricName;
  /// the name of continuous diff(erence) function
  std::string numDiffMetricName;
  /// number of instances to sample
  unsigned int m;
  /// are instances being randomly selected?
  bool randomlySelect;
  /// number of attributes to remove each iteration if running iteratively
  unsigned int removePerIteration;
  /// are we removing a percentage per iteration?
  bool doRemovePercent;
  /// percentage of attributes to remove per iteration if running iteratively
  double removePercentage;
  /// number of target attributes
  unsigned int numTarget;
  /// name of the weight-by-distance method
  std::string weightByDistanceMethod;
  /// sigma value used in exponential decay weight-by-distance
  double weightByDistanceSigma;
  /// RAW attribute scores/weights - ie, no normalization
  std::vector<double> W;

  // --------------------------------------------------------------------------
  // kopt/best k algorithm from ReliefSeq - bcw - 10/7/16
  /*************************************************************************//**
	 * Remove the worst attribute based on free energy scores.
	 * \param [in] numToRemove number of attributes to remove/evaporate
	 * \return distance
	 ****************************************************************************/
  bool RemoveWorstAttributes(unsigned int numToRemove=1);
  bool SetKoptParameters();
  /// Determine the maximum k value for optimization.
  unsigned int GetKmax();
  
  /// attributes that have been evaporated so far
  AttributeScores removedAttributes;
  /// optimize k begin value
  unsigned int koptBegin;
  /// optimize k end value
  unsigned int koptEnd;
  /// optimize k step value
  unsigned int koptStep;
  /// best k by attribute
  std::map<std::string, unsigned int> bestKs;
  
  bool ComputeGRM();
};

#endif
