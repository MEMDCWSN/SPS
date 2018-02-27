//
//  main_convert_psm_to_proteosafe - Change PSMs to form suitable for proteosafe
//
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"
#include "SpectrumPairSet.h"

#include <stdlib.h>

#define DEBUG_ 0

using namespace std;
using namespace specnets;

// -------------------------------------------------------------------------
void addFlankingAa(psmPtr psm,
                   DB_fasta & dbFasta)
{
  string cleanAnno;
  PeptideSpectrumMatch::getUnmodifiedPeptide(psm->m_annotation, cleanAnno);
  string dbSequence = dbFasta.getSequence(psm->m_dbIndex);
  string prevAA = "_";
  string postAA = "_";

  int start = dbSequence.find(cleanAnno);
  if (start != string::npos) {
    if (start != 0) {
      prevAA = dbSequence[start - 1];
    }
    int end = start + cleanAnno.length() - 1;
    if (end < dbSequence.length() - 1) {
      postAA = dbSequence[end + 1];
    }
  }

  string newAnnotation  = prevAA + "." + psm->m_annotation + "." + postAA;
  //DEBUG_VAR(newAnnotation);
  psm->m_annotation = newAnnotation;
  //DEBUG_VAR(psm->m_annotation);
}


// -------------------------------------------------------------------------
int roundedMass(float mass) 
{
  return int(abs(mass) + 0.5) * (mass < 0 ? -1 : 1);
}


void displayUsage(void)
{
  cerr << "Usage: main_convert_psm_to_proteosafe <psm_file> <out_psm_file> [options]" << endl;
  cerr << "   -addflank <db_file>   Will add flanking AAs to the PSMs (requires DB file)" << endl;
  cerr << "   -fixc <aa_file>       Add fixed mods to cysteines (requires AA file)" << endl;
  cerr << "   -outdir <directoryu>  Output directory for output files (if none specified current directory will be used)" << endl;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 3) {
    displayUsage();
    return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("addflank", "ADD_FLANKING_AAS", 1));
  listOptions.push_back(CommandLineParser::Option("fixc", "FIX_CYSTEINES", 1));
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));  

  CommandLineParser clp(argc, argv, 2, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    displayUsage();
    cerr << "Invalid options" << endl;
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);
  
  PeptideSpectrumMatchSet	psmSet;
  DEBUG_MSG("Loading PSM file [" << argv[1] << "]");
  if (!psmSet.loadFromFile(argv[1])) {
    ERROR_MSG("Loading PSM file [" << argv[1] << "]");
    return -2;
  }

  bool fixCysteines = commandLineParams.exists("FIX_CYSTEINES");
  DEBUG_VAR(fixCysteines);
  string aaFile = commandLineParams.getValue("FIX_CYSTEINES");
  DEBUG_VAR(aaFile);

  // Load amino acid masses
  DEBUG_MSG("Loading Database ...");
  AAJumps jumps(1);
  if (!aaFile.empty()) {
    jumps.loadJumps(aaFile.c_str());
  }

  for (int i = 0; i < psmSet.size(); i++) {
    psmSet[i]->changeGapAnnosToSingle();
    if (fixCysteines) {
      psmSet[i]->addFixedCysteineMods();
    }
  }

  bool addFlankingAas = commandLineParams.exists("ADD_FLANKING_AAS");
  DEBUG_VAR(addFlankingAas);
  
  if (addFlankingAas) {
  
    DB_fasta dbFasta;
    string dbFile = commandLineParams.getValue("ADD_FLANKING_AAS");
    DEBUG_MSG("Loading Database ...");
    if (dbFasta.Load(dbFile.c_str()) <= 0) {
      ERROR_MSG("Loading FASTA file [" << dbFile.c_str() << "]");
      return -3;
    }
  
    for (int i = 0; i < psmSet.size(); i++) {
    
      // Can't do decoys - they are from different DB
      if (psmSet[i]->m_isDecoy) {
        continue;
      }
      addFlankingAa(psmSet[i], dbFasta);
    }
  }
  
  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[2];
  DEBUG_MSG("Saving file [" << outputFile << "]");
  psmSet.saveToFile(outputFile.c_str(), true, true);

  return 0;
}

