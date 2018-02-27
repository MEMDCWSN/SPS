//
// util_parsimony - Create parsimonious PSM set
//
#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"
#include <stdlib.h>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 4)
  {
    cerr << "Usage: convertpsm input_type PSM_file_in PSM_file_out" << endl;
    cerr << "       Valid types are: inspect msgfdb msgfplus moda msplit msplit2 multipass specnets specres specrep" << endl << endl;
    cerr << "       Converts various PSM file types to the Specnets standard format" << endl;
    return -1;
  }
  
  PeptideSpectrumMatchSet psmSet;

  string type(argv[1]);
  DEBUG_VAR(type);

  DEBUG_MSG("Loading file [" << argv[2] << "] as type [" << type << "]");
  if (type == "inspect") {
    if (!psmSet.loadInspectResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [inspect]");
      return -1;
    }
  } else if (type == "msgfdb") {
    if (!psmSet.loadMSGFDBResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [msgfdb]");
      return -1;
    }
  } else if (type == "msgfplus") {
    if (!psmSet.loadMSGFPlusResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [msgfplus]");
      return -1;
    }
  } else if (type == "moda") {
    if (!psmSet.loadModaResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [moda]");
      return -1;
    }
  } else if (type == "msplit") {
    if (!psmSet.loadMsplitResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [msplit]");
      return -1;
    }
  } else if (type == "msplit2") {
    if (!psmSet.loadMsplitResultsFile2(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [msplit2]");
      return -1;
    }
  } else if (type == "multipass") {
    if (!psmSet.loadMultipassResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [multipass]");
      return -1;
    }
  } else if (type == "specnets") {
    if (!psmSet.loadFromFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [specnets]");
      return -1;
    }
  } else if (type == "specrep") {
    if (!psmSet.loadSpecnetsReportFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [specnets report]");
      return -1;
    }
  } else if (type == "specres") {
    if (!psmSet.loadSpecnetsResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [specnets result]");
      return -1;
    }
  } else {
      ERROR_MSG("Unknown type [" << argv[2] << "]");
      return -1;
  }
  //DEBUG_VAR(psmSet.size());
  
  // Add information is possible
  for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
    vector<float> modifications;
    vector<unsigned int> positions;
    vector<unsigned int> lengths;
    psmSet[iPsm]->getModificationsAndPositions(modifications, positions, lengths, false);
    psmSet[iPsm]->m_numMods = modifications.size();
  }

  DEBUG_MSG("Saving file [" << argv[3] << "]");
  psmSet.saveToFile(argv[3]);

  return 0;
}
