//
//  main_dump_pairs_histogram - output the histogram of mass differences from pairs
//
#include "CommandLineParser.h"
#include "Logger.h"
#include "SpectrumPairSet.h"

#include <stdlib.h>

#define DEBUG_ 0

const int TOP_N_DEFAULT = 1000;

using namespace std;
using namespace specnets;

// -------------------------------------------------------------------------
int roundedMass(float mass) 
{
  return int(abs(mass) + 0.5) * (mass < 0 ? -1 : 1);
}

// -------------------------------------------------------------------------
void displayUsage(void)
{
  cerr << "Usage: main_dump_pairs_histogram <pair_file> <out_file> [options]" << endl;
  cerr << "   -maxdiff N            Maximum mass difference to list (default is 200)" << endl;
  cerr << "   -outdir <directoryu>  Output directory for output files (if none specified current directory will be used)" << endl;
  cerr << "   -res N                Resolution for mass bins (default is 1)" << endl;
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
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));
  listOptions.push_back(CommandLineParser::Option("res", "RESOLUTION", 1));
  listOptions.push_back(CommandLineParser::Option("maxdiff", "MAX_MASS_DIFF", 1));
  

  CommandLineParser clp(argc, argv, 2, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    displayUsage();
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);
  
  SpectrumPairSet pairs;
  string pairsFileName = argv[1];
  DEBUG_MSG("Loading pairs file [" << pairsFileName << "]...");
  if (!pairs.loadFromBinaryFile(pairsFileName.c_str())) {
    ERROR_MSG("Could not load " << pairsFileName);
    return false;
  }

  float resolution = commandLineParams.getValueFloat("PARENT_MASS_TOLERANCE", 1.0);
  DEBUG_VAR(resolution);

  float maxDiffMass = commandLineParams.getValueFloat("MAX_MASS_DIFF", 200.0);
  DEBUG_VAR(maxDiffMass);

  map<float, float> modFreqs;
  if (!pairs.getModificationFrequencies(resolution, maxDiffMass, modFreqs)) {
    ERROR_MSG("Could not get mod frequencies!");
    return false;
  }
  
  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[2];

  DEBUG_TRACE;
  ofstream ofs(outputFile.c_str());
  if (!ofs) {
    ERROR_MSG("Could not write to " << outputFile);
  }
  ofs << "mass\tfreq" << endl;

  map<float, float>::iterator itr = modFreqs.begin();
  map<float, float>::iterator itr_end = modFreqs.end();
  for(; itr != itr_end; itr++) {
    float mass = itr->first;
    float freq = itr->second;
    if (mass > maxDiffMass) continue;
    ofs << mass << "\t"  << freq  << endl;
  }
  ofs.close();

  return 0;
}

