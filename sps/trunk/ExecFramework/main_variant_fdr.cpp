//
// util_parsimony - Create parsimonious PSM set
//
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "ExecFdrPeptide.h"
#include "Logger.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"

#include <stdlib.h>

using namespace specnets;
using namespace std;
using namespace sps;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 3) {
    cerr << "Usage: main_variant_fdr psm_file pmtol out_psm_file [options]" << endl;
    cerr << "     -decoy filename  Input decoy psm file for 2 file [target is psm_file]" << endl;
    cerr << "     -msgfdb          Input PSMs are in MSGFDB style file" << endl;
    cerr << "     -threshold N     FDR cutoff threshold (default is 0.1)" << endl;
    cerr << "     -pmtol N         Parent mass tolerance for variant group calculation" << endl;
    cerr << "     -type T          Type of FDR: concatenated, merged, split, separate)" << endl;
    cerr << "                                    (default is concatenated)" << endl;
    cerr << "     -usepvalue       Use pvalue for comparison instead of score" << endl;
	  return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("decoy", "INPUT_DECOY_PSMS", 1));
  listOptions.push_back(CommandLineParser::Option("threshold", "PEPTIDE_FDR_CUTOFF", 1));
  listOptions.push_back(CommandLineParser::Option("type", "TDA_TYPE", 1));
  listOptions.push_back(CommandLineParser::Option("usepvalue", "USE_PVALUE_COMPARE", 0));
  listOptions.push_back(CommandLineParser::Option("msgfdb", "LOAD_MSGFDB_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));

  CommandLineParser clp(argc, argv, 3, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "Usage: main_variant_fdr psm_file pmtol out_psm_file [options]" << endl;
    cerr << "     -decoy filename  Input decoy psm file for 2 file [target is psm_file]" << endl;
    cerr << "     -msgfdb          Input PSMs are in MSGFDB style file" << endl;
    cerr << "     -threshold N     FDR cutoff threshold (default is 0.1)" << endl;
    cerr << "     -pmtol N         Parent mass tolerance for variant group calculation" << endl;
    cerr << "     -type T          Type of FDR: concatenated, merged, split, separate)" << endl;
    cerr << "                                    (default is concatenated)" << endl;
    cerr << "     -usepvalue       Use pvalue for comparison instead of score" << endl;
    cerr << parserError << endl;
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  PeptideSpectrumMatchSet	psmSet;
  DEBUG_MSG("Loading target PSM file [" << argv[1] << "]");
  if (commandLineParams.exists("LOAD_MSGFDB_STYLE_FILE")) {
    if (!psmSet.loadMSGFDBResultsFile(argv[1])) {
      ERROR_MSG("Loading target PSM file [" << argv[1] << "]");
      return -2;
    }
    // This should be looked at, but MScore and p-value are reversed here?
    for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
      float temp = psmSet[iPsm]->m_pValue;
      psmSet[iPsm]->m_pValue = psmSet[iPsm]->m_score;
      psmSet[iPsm]->m_score = temp;
    }
  } else {
    if (!psmSet.loadFromFile(argv[1])) {
      ERROR_MSG("Loading target PSM file [" << argv[1] << "]");
      return -2;
    }
  }

  // Make sure the isDecoy fields are set on any decoy PSMs
  for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
    if (psmSet[iPsm]->m_protein.find("XXX") != string::npos ||
        psmSet[iPsm]->m_protein.find("REV") != string::npos) {
      psmSet[iPsm]->m_isDecoy = true;
    }
  }
  PeptideSpectrumMatchSet	psmSetDecoy;
  if (commandLineParams.exists("INPUT_DECOY_PSMS")) {
    string decoyFile = commandLineParams.getValue("INPUT_DECOY_PSMS");
    DEBUG_MSG("Loading decoy PSM file [" << decoyFile << "]");
    if (!psmSetDecoy.loadFromFile(decoyFile.c_str())) {
      ERROR_MSG("Loading decoy PSM file [" << decoyFile << "]");
      return -2;
    }
    
    // Make sure the isDecoy fields are set on the decoy PSMs
    for (int iDecoy = 0; iDecoy < psmSetDecoy.size(); iDecoy++) {
      psmSetDecoy[iDecoy]->m_isDecoy = true;
    }
    psmSetDecoy.saveToFile("debug_psm_decoy.txt");
  }

  if (!commandLineParams.exists("TDA_TYPE")) {
    commandLineParams.setValue("TDA_TYPE", "concatenated");
  }
  
  if (commandLineParams.exists("USE_PVALUE_COMPARE")) {
    commandLineParams.setValue("PEPTIDE_FDR_USE_PVALUE_COMPARE", "1");
  }
  
  if (!commandLineParams.exists("PEPTIDE_FDR_CUTOFF")) {
    commandLineParams.setValue("PEPTIDE_FDR_CUTOFF", "0.01");
  }
    
  commandLineParams.setValue("PEPTIDE_FDR_REMOVE_DECOYS", "1");
  commandLineParams.setValue("USE_VARIANT_GROUPS", "1");
  commandLineParams.setValue("TOLERANCE_PM", argv[2]);

  PeptideSpectrumMatchSet psmSetFdr;
  if (psmSetDecoy.size() != 0) {
    ExecFdrPeptide moduleFdrPeptide(commandLineParams,
                                    &psmSet,
                                    &psmSetDecoy,
                                    &psmSetFdr);

    bool returnStatus = moduleFdrPeptide.invoke();
    if (!returnStatus) {
      ERROR_MSG("Problem invoking ExecFdrPeptide!");
      return -1;
    }
  } else {
    ExecFdrPeptide moduleFdrPeptide(commandLineParams,
                                    &psmSet,
                                    &psmSetFdr);

    bool returnStatus = moduleFdrPeptide.invoke();
    if (!returnStatus) {
      ERROR_MSG("Problem invoking ExecFdrPeptide!");
      return -1;
    }
  }

  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[3];
  DEBUG_MSG("Saving file [" << outputFile << "]");
  psmSetFdr.saveToFile(outputFile.c_str(), true, true);

  return 0;
}

