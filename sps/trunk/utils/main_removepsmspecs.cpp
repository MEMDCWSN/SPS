//
//  main_removepsmspecs - removes spectra that have a PSM
//
#include "abruijn.h"
#include "ClusterData.h"
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ReportTableClusterConsensus.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include "tuple.h"

#include <stdlib.h>

#define DEBUG_ 0

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;


// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 4) {
    cerr << "Usage: main_removepsmspecs in_spec_pklbin in_psm_file out_spec_file" << endl;
    return -1;
  }

  SpecSet inSpectra;
  SpecSet outSpectra;

  inSpectra.Load(argv[1]);
  size_t size1 = inSpectra.size();
  DEBUG_VAR(size1);
  if (size1 == 0) {
    ERROR_MSG("Loading spectra from [" << argv[1] << "]");
    return -2;
  }
  
  PeptideSpectrumMatchSet	psmSet;
  if (!psmSet.loadFromFile(argv[2])) {
    ERROR_MSG("Loading PSM file [" << argv[2] << "]");
    return -3;
  }
  DEBUG_VAR(psmSet.size());
  
  outSpectra.resize(size1);

  //----------------------------------------------------------
  // Create a PSM map
  //----------------------------------------------------------
  map<int, psmPtr> mapPsm;
  for (int iPSM = 0; iPSM < psmSet.size(); iPSM++ ) {
    psmPtr psmSpec = psmSet[iPSM];
    mapPsm[psmSpec->m_scanNum] = psmSpec;
  }
  DEBUG_VAR(mapPsm.size());

   for (size_t i = 0; i < size1; i++) {
     int scanNum = inSpectra[i].scan;
     if (mapPsm.find(scanNum) == mapPsm.end()) {
       outSpectra.push_back(inSpectra[i]);
     }
  }

  if (outSpectra.SaveSpecSet(argv[3]) == 0) {
    ERROR_MSG("Saving spectra to [" << argv[3] << "]");
    return -4
;  }

  
  return 0;
}
