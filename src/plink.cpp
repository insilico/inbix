//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include "plink.h"
#include "options.h"

using namespace std;

////////////////////////////////
// Clean-up

void Plink::cleanUp()
{

  for (int i=0; i<sample.size(); i++)
    delete sample[i];
  
  for (int l=0; l<locus.size(); l++)
    delete locus[l];      
  
  if (par::SNP_major)
    for (int l=0; l<SNP.size(); l++)
      delete SNP[l];      

}

