
/**
 * \class RReliefF
 *
 * \brief Regression ReliefF attribute ranking algorithm.
 *
 * Totally redone for the McKinney insilico lab in 2011.
 * Large refactoring to move all attribute elimination handling to the
 * Dataset and its subclasses. 9/11/11
 *
 * \sa ReliefF
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 9/27/11
 */

#ifndef RRELIEFF_H
#define	RRELIEFF_H

#include <vector>
#include <fstream>

#include "plink.h"

#include "ReliefF.h"
#include "Dataset.h"
#include "Insilico.h"

class RReliefF : public ReliefF
{
public:
  /*************************************************************************//**
   * Construct an ReliefF algorithm object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  RReliefF(Dataset* ds, Plink* plinkPtr);
  bool ComputeAttributeScores();
  virtual ~RReliefF();
private:  
};

#endif	/* RRELIEFF_H */
