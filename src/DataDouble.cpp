/*-------------------------------------------------------------------------------
This file is part of Ranger.
    
Ranger is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Ranger is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Ranger. If not, see <http://www.gnu.org/licenses/>.

Written by: 

Marvin N. Wright
Institut für Medizinische Biometrie und Statistik
Universität zu Lübeck
Ratzeburger Allee 160
23562 Lübeck 

http://www.imbs-luebeck.de
wright@imbs.uni-luebeck.de
#-------------------------------------------------------------------------------*/

#include <string>

#include "DataDouble.h"

// inbix headers
#include "Insilico.h"
#include "Dataset.h"
#include "InteractionNetwork.h"
#include "plink.h"
#include "options.h"
#include "helper.h"

DataDouble::DataDouble() :
    data(0) {
}

DataDouble::~DataDouble() {
  if (!externalData) {
    delete[] data;
  }
}

bool DataDouble::loadFromPlink(Plink* PP) {
  /* load a Ranger data set from pointer to a PLINK environment 
   * bcw - October 2016
   */
  num_cols = PP->nlistname.size() + 1;
  num_rows = PP->sample.size();
  num_cols_no_sparse = num_cols;
  variable_names.clear();
  for(int i=0; i < PP->nlistname.size(); i++) {
    string numericName = PP->nlistname[i];
    variable_names.push_back(numericName);
  }
  variable_names.push_back("Class");
  reserveMemory();
  
  bool error = false;
  for(int i=0; i < PP->sample.size(); i++) {
    string ID = PP->sample[i]->fid + PP->sample[i]->iid;
    if(PP->sample[i]->missing) {
      // missing phenotype so skip this individual
      if(par::verbose) {
        PP->printLOG("ID: " + ID + " missing, skipping\n");
      }
      continue;
    }
    unsigned int numNumerics = PP->nlistname.size();
    for(int j=0; j < numNumerics; j++) {
      NumericLevel numeric = PP->sample[i]->nlist[j];
      set(j, i, numeric, error);
    }
    double phenotype = -9;
    if(par::bt) {
      int intPheno = PP->sample[i]->aff? 2: 1;
      phenotype = static_cast<double>(intPheno);
    } else {
      phenotype = PP->sample[i]->phenotype;
    }
    set(variable_names.size() - 1, i, phenotype, error);
  }

  externalData = false;

  return error;
}

bool DataDouble::loadFromDatasetMask(Dataset* ds, vector<string> bestAttributes) {
  // load a Ranger data set from my Dataset class, honoring inclusion masks
  if(ds->HasGenotypes()) {
    error("Genotype attributes not supported yet");
  }
  PP->printLOG(Timestamp() + "Copying *masked* PLINK data structures to this Data object\n");

  // set class variables for reserving memory and other operations
  num_cols = bestAttributes.size() + 1;
  num_rows = ds->NumInstances();
  num_cols_no_sparse = num_cols;
  variable_names.resize(num_cols);
  copy(bestAttributes.begin(), bestAttributes.end(), variable_names.begin());
  variable_names[bestAttributes.size()] = (par::depvarname);
  if(data) delete [] data;
  reserveMemory(); // allocates num_cols * num_rows double matrix

  // loop over all instances
  bool error = false;
  vector<string> instanceIds = ds->GetInstanceIds();
  vector<Indices> numInstanceIdx = ds->MaskGetInstanceIndices();
  map<string, unsigned int> varMap = ds->MaskGetAttributeMask(NUMERIC_TYPE);
  for(int i=0; i < numInstanceIdx.size(); i++) {
    string ID = instanceIds[i];
    Indices thisInstanceIndex = numInstanceIdx[i];
    // loop over masked numeric attributes
    for(int j=0; j < bestAttributes.size(); j++) {
      string thisVarName = bestAttributes[j];
      Indices thisAttributeIndex = varMap[thisVarName];
      NumericLevel numeric = 
              ds->GetInstance(thisInstanceIndex)->GetNumeric(thisAttributeIndex);
      set(j, i, numeric, error);
    }
    // phenotype is last column
    double phenotype = -9; 
    if(par::bt) {
      int intPheno = (ds->GetInstance(i)->GetClass()) == 0? -1: 1;
      phenotype = static_cast<double>(intPheno);
    } else {
      phenotype = ds->GetInstance(i)->GetPredictedValueTau();
    }
    set(bestAttributes.size(), i, phenotype, error);
  }

  externalData = false;

  return error;
}
