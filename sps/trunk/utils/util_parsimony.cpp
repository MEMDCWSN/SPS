//
// util_setdbindexes - Lookup protein names and find their indexes in the DB
//
#include "Logger.h"
#include "db_fasta.h"
#include <stdlib.h>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 5)
  {
    cerr << "Usage: util_parsimony input_type PSM_file_in PSM_file_out db_file" << endl;
    cerr << "       Valid types are: inspect msgfdb moda psm" << endl;
    return -1;
  }
  
  PeptideSpectrumMatchSet psmSet;

  string type(argv[1]);
  DEBUG_VAR(type);

  //DEBUG_MSG("Loading file [" << argv[2] << "] as type [" << type << "]");
  if (type == "inspect") {
    if (!psmSet.loadInspectResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [inspect]");
      return -1;
    }
  } else if (type == "msfgdb") {
    if (!psmSet.loadMSGFDBResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [msgfdb]");
      return -1;
    }
  } else if (type == "moda") {
    if (!psmSet.loadModaResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [moda]");
      return -1;
    }
  } else if (type == "psm") {
    if (!psmSet.loadFromFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [psm]");
      return -1;
    }
  } else {
      ERROR_MSG("Unknown type [" << argv[2] << "] as type [moda]");
      return -1;
  }
  DEBUG_VAR(psmSet.size());

  DB_fasta db;
  unsigned int dbSize = db.Load(argv[4]);
  if (dbSize == 0) {
      ERROR_MSG("Error loading db file [" << argv[4] << "]");
      return -1;
  }

  db.replaceAA('L', 'I');
  db.replaceAA('K', 'Q');

  PeptideSpectrumMatchSet psmSetOut;

  for (int i = 0; i < psmSet.size(); i++) {
    psmPtr p = psmSet[i];
    string peptide(p->m_annotation);
    size_t dotpos = peptide.find('.',0);
    if (dotpos != string::npos) {
      peptide.erase(0,dotpos+1);
      size_t dotpos = peptide.find('.',0);
      if (dotpos != string::npos) {
        peptide.erase(dotpos+1,peptide.length());
      }
    }
    for (int i = 0; i < peptide.length(); i++) {
      if (peptide[i] == 'L') peptide[i] = 'I';
      if (peptide[i] == 'K') peptide[i] = 'Q';
    }

    list<sps::tuple<int,float,string> > matches;
    unsigned int numResults = db.find(peptide, matches);
    if (numResults == 0) {
      WARN_MSG(peptide << " not found in database!");
    }

    list<sps::tuple<int,float,string> >::iterator itr = matches.begin();
    list<sps::tuple<int,float,string> >::iterator itr_end = matches.end();
    for ( ; itr != itr_end; itr++) {
      psmPtr p2(new PeptideSpectrumMatch(*p));
      int dbIndex = itr->m0;
      p2->m_protein = db.IDs[dbIndex];
      p2->m_dbIndex = dbIndex;
      psmSetOut.push_back(p2);
    }
  }

  DEBUG_VAR(psmSetOut.size());
  psmSetOut.maximumParsimony();
  DEBUG_VAR(psmSetOut.size());

  DEBUG_MSG("Saving file [" << argv[3] << "]");
  psmSetOut.saveToFile(argv[3]);

  return 0;
}
