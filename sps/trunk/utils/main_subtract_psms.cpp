//
//  main_subtract_psms - subtract one set of psms from another
//
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"

#include <stdlib.h>

#define DEBUG_ 0


using namespace std;
using namespace sps;
using namespace specnets;

// -------------------------------------------------------------------------
void displayUsage(void)
{
  cerr << "Usage: main_subtract_psms <psm_file> <subtract_psm_file> <out_file> [options]" << endl;
  cerr << "   -msgfdb               Input PSMs are in MSGFDB format file" << endl;
  cerr << "   -moda                 Input IDs are in MODa (multiple results) style file" << endl;
  cerr << "   -msplit               Input IDs are in MSplit format file" << endl;
  cerr << "   -outdir <directory>   Output directory for output files (if none specified current directory will be used)" << endl;
  cerr << "   -invert               Invert the subtraction list (only keep PSMs in the subtraction list" << endl;
}


// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 4) {
    displayUsage();
	  return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("msgfdb", "LOAD_MSGFDB_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("moda", "LOAD_MODA_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("msplit", "LOAD_MSPLIT_STYLE_FILE", 0));
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));  
  listOptions.push_back(CommandLineParser::Option("invert", "INVERT_SUBTRACTION_LIST", 0));
  
  CommandLineParser clp(argc, argv, 3, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    displayUsage();
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);
  
  PeptideSpectrumMatchSet	psmSet1;
  if (!psmSet1.loadFromFile(argv[1])) {
    ERROR_MSG("Loading PSM file [" << argv[1] << "]");
    return -2;
  }
  DEBUG_VAR(psmSet1.size());
  
  PeptideSpectrumMatchSet	psmSet2;
  if (!psmSet2.loadFromFile(argv[2])) {
    ERROR_MSG("Loading PSM file [" << argv[2] << "]");
    return -3;
  }
  DEBUG_VAR(psmSet2.size());

  //----------------------------------------------------------
  // Create the PSM map
  //----------------------------------------------------------
  map<int, psmPtr> mapPsm2;
  for (int iPSM = 0; iPSM < psmSet2.size(); iPSM++ ) {
    psmPtr psmSpec = psmSet2[iPSM];
    mapPsm2[psmSpec->m_scanNum] = psmSpec;
    //DEBUG_VAR(psmSpec->m_scanNum);
  }
  DEBUG_VAR(mapPsm2.size());

  bool invert = commandLineParams.exists("INVERT_SUBTRACTION_LIST");
  DEBUG_VAR(invert);

  //----------------------------------------------------------
  // Subtract the PSMs
  //----------------------------------------------------------

  PeptideSpectrumMatchSet	psmSetOutput;
  for (int iPsm = 0; iPsm < psmSet1.size(); iPsm++) {
    psmPtr psmSpec = psmSet1[iPsm];
    //DEBUG_VAR(psmSpec->m_scanNum);
    bool found = mapPsm2.find(psmSpec->m_scanNum) != mapPsm2.end();
    if (invert && found) {
      psmSetOutput.push_back(psmSpec);
    } else if (!invert && !found) {
      psmSetOutput.push_back(psmSpec);
    }
  }  


  //----------------------------------------------------------
  // Output results
  //----------------------------------------------------------
  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[3];

  DEBUG_MSG("Saving to file [" << outputFile << "]");
  psmSetOutput.saveToFile(outputFile.c_str(), true, true);
  
  return 0;
}
