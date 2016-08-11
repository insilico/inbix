/**
 * \class ReliefFSeq
 *
 * \brief ReliefFSeq attribute ranking algorithm.
 *
 * Designed to handle digital gene expression (DGE) data sets, particularly
 * RNA-Seq high-throughput count data, by accounting for variable-specific
 * variance in counts. Data is known to follow a Poisson or negative binomial
 * distribution. This algorithm is a more computationally practical approach
 * than others that use more sophisticated statistical methods and models.
 * Our approach keeps the ReliefF algorithm general while addressing the
 * variance "dispersion" problem as a special case.
 *
 * \sa ReliefF
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 7/23/12
 */

#ifndef RELIEFFSEQ_H
#define	RELIEFFSEQ_H

#include <vector>
#include <map>
#include <string>
#include <fstream>

#include "plink.h"

#include "ReliefF.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "Insilico.h"

class ReliefFSeq : public ReliefF
{
public:
  /*************************************************************************//**
   * Construct an ReliefFSeq algorithm object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  ReliefFSeq(Dataset* ds, Plink* plinkPtr);
  bool ComputeAttributeScores();
  AttributeScores GetScores();
  // average hit and miss diffs for gene alpha
  std::pair<double, double> MuDeltaAlphas(unsigned int alpha);
  /// standard deviations of hit and miss diffs for gene alpha
  std::pair<double, double> SigmaDeltaAlphas(unsigned int alpha,
  		double muDeltaHit, double muDeltaMiss);
  virtual ~ReliefFSeq();
private:
	/// ReliefSeq mode: signal-to-noise ratio (snr) or t-statistic (tstat)
  std::string mode;
	/// ReliefSeq signal-to-noise ratio mode: signal-to-noise ratio (snr) or
	/// ReliefF (relieff)
	std::string snrMode;
	/// ReliefFSeq t-statistic mode: 1-pvalue (pval) or the absolute value
	/// of the t-statistic (abst)
	std::string tstatMode;
	/// variance denominator adjustment s0
  double s0;
};

#endif	/* SNReliefF_H */
