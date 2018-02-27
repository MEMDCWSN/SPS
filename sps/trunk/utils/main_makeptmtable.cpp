//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "CommandLineParser.h"
#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"

#include <stdlib.h>

#define DEBUG_ 0

using namespace std;
using namespace specnets;



// -------------------------------------------------------------------------
void displayUsage(void)
{
  cerr << "Usage: main_makeptmtable <psm_file> <output_ptmtable>" << endl;
  cerr << "   -recurse              Compute recursive localization" << endl;
  cerr << "   -newpsms <filename>   Output new PSM with adjusted mod locations from recursive PTM table" << endl;
  cerr << "   -outdir <directory>   Output directory for output files (if none specified current directory will be used)" << endl;
  cerr << "   -useorig              Use original annotation" << endl;
  cerr << "   -excludeC <mass>      exclude single AA PTMs on C with the given mass" << endl;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 3) {
    displayUsage();
    return -1;
  }

  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("recurse", "RECURSIVE_PTM_TABLE", 0));
  listOptions.push_back(CommandLineParser::Option("useorig", "USE_ORIGINAL_ANNO", 0));
  listOptions.push_back(CommandLineParser::Option("newpsms", "OUTPUT_NEW_PSMS", 1));
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));
  listOptions.push_back(CommandLineParser::Option("excludeC", "EXCLUDE_C_MASS", 1));

  CommandLineParser clp(argc, argv, 2, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    displayUsage();
    cerr << parserError << endl;
    return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  PeptideSpectrumMatchSet	psmSet;
  if (!psmSet.loadFromFile(argv[1])) {
    ERROR_MSG("Loading PSM file [" << argv[1] << "]");
    return -1;
  }

  bool useOrig = commandLineParams.exists("USE_ORIGINAL_ANNO");
  DEBUG_VAR(useOrig);

  bool recurse = commandLineParams.exists("RECURSIVE_PTM_TABLE");
  DEBUG_VAR(recurse);

  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[2];

  vector<pair<string,float> > exclusionList;
  if (commandLineParams.exists("EXCLUDE_C_MASS")) {
    float mass = commandLineParams.getValueFloat("EXCLUDE_C_MASS");
    pair<string,float> excludePair = make_pair<string,float>("C", mass);
    exclusionList.push_back(excludePair);
  }
  listOptions.push_back(CommandLineParser::Option("excludeC", "EXCLUDE_C_MASS", 1));

  DEBUG_MSG("Saving PTM table to file [" << outputFile << "]");
  if (!psmSet.saveModMatrix(outputFile.c_str(), useOrig, recurse, &exclusionList)) {
    ERROR_MSG("Error saving PTM table to file [ " << outputFile << "]");
    return -1;
  }

  if (commandLineParams.exists("OUTPUT_NEW_PSMS")) {
    map<string, map<float, float> > mapModCount;   // AA = mod, count
    PeptideSpectrumMatchSet returnSet;
    psmSet.getRecursiveModCounts(mapModCount, returnSet, useOrig, &exclusionList);

    string outputFile2;
    if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
      outputFile2 += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
      outputFile2 += "/";
    }
    outputFile2 += commandLineParams.getValue("OUTPUT_NEW_PSMS");
    returnSet.saveToFile(outputFile2.c_str(), true, true);
  }

  return 0;
}
