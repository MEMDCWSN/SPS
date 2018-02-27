//
//  main_create_psm_specific_penalties - Find the PTMs that are supported by pair data
//
#include "CommandLineParser.h"
#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"
#include "SpectrumPairSet.h"

#include <stdlib.h>

#define DEBUG_ 0
#define DEBUG_PTM_TABLE 0
#define DEBUG_SCAN -1 
#define DEBUG_MASS 100000

const int TOP_N_DEFAULT = 1000;

using namespace std;
using namespace specnets;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope


  if (argc < 4) {
    cerr << "Usage: main_create_psm_specific_penalties <pair_file> <star_file> <out_file> [options]" << endl;
    cerr << "   -maxmass N            Maximum mass tolerance to consider" << endl;
    cerr << "   -minmass N            Minimum mass tolerance to consider" << endl;
    cerr << "   -outdir <directoryu>  Output directory for output files (if none specified current directory will be used)" << endl;
    cerr << "   -pmtol N              Parent mass tolerance for variant group calculation" << endl;
    return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("maxmass", "MAXIMUM_MASS", 1));
  listOptions.push_back(CommandLineParser::Option("minmass", "MINIMUM_MASS", 1));
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));
  listOptions.push_back(CommandLineParser::Option("pmtol", "PARENT_MASS_TOLERANCE", 1));

  CommandLineParser clp(argc, argv, 3, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "Usage: main_create_psm_specific_penalties <pair_file> <star_file> <out_file> [options]" << endl;
    cerr << "   -maxmass N            Maximum mass tolerance to consider" << endl;
    cerr << "   -minmass N            Minimum mass tolerance to consider" << endl;
    cerr << "   -outdir <directoryu>  Output directory for output files (if none specified current directory will be used)" << endl;
    cerr << "   -pmtol N              Parent mass tolerance for variant group calculation" << endl;
    cerr << "Invalid options" << endl;
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);
  
  float parentMassTol = commandLineParams.getValueFloat("PARENT_MASS_TOLERANCE", 1.0);
  DEBUG_VAR(parentMassTol);

  float maximumMass = commandLineParams.getValueFloat("MAXIMUM_MASS", 200.0);
  DEBUG_VAR(maximumMass);
  float minimumMass = commandLineParams.getValueFloat("MINIMUM_MASS", -200.0);
  DEBUG_VAR(minimumMass);

  SpectrumPairSet pairs;
  SpecSet starSpectra;

  string pairsFileName = argv[1];
  DEBUG_MSG("Loading pairs file [" << pairsFileName << "]...");
  if (!pairs.loadFromBinaryFile(pairsFileName.c_str())) {
    ERROR_MSG("Could not load " << pairsFileName);
    return false;
  }

  string starsFileName = argv[2];
  DEBUG_MSG("Loading stars file [" << starsFileName << "]...");
  if (!starSpectra.loadPklBin(starsFileName.c_str())) {
    ERROR_MSG("Could not load " << starsFileName);
    return false;
  }

  map<int, float> mapSpecToMass;
  for (int iStar = 0; iStar < starSpectra.size(); iStar++) {
    mapSpecToMass[starSpectra[iStar].scan] = starSpectra[iStar].parentMass;
  }

  //--------------------------------------------------------------------------
  map<int, map<int, float> > mapSpecToMassToScore;
  DEBUG_VAR(pairs.size()); 
  for (int iPair = 0; iPair < pairs.size(); iPair++) {
    int spec1 = pairs[iPair].spec1 + 1;  // These are 0 based
    int spec2 = pairs[iPair].spec2 + 1;
    
    float rawDiff1 = mapSpecToMass[spec1] - mapSpecToMass[spec2];
    if (DEBUG_) DEBUG_VAR(rawDiff1);

    // not interested in small mass diffs
    if (abs(rawDiff1) < parentMassTol) {
      continue;
    }
    
    int   sign1    = rawDiff1 < 0 ? -1 : 1;
    if (DEBUG_) DEBUG_VAR(sign1);
    float rawDiff2 = mapSpecToMass[spec2] - mapSpecToMass[spec1];
    if (DEBUG_) DEBUG_VAR(rawDiff2);
    int   sign2    = rawDiff2 < 0 ? -1 : 1;
    if (DEBUG_) DEBUG_VAR(sign2);
    
    int massDiff1 = (int)(abs(rawDiff1) + 0.5) * sign1;
    if (DEBUG_) DEBUG_VAR(massDiff1);
    int massDiff2 = (int)(abs(rawDiff2) + 0.5) * sign2;
    if (DEBUG_) DEBUG_VAR(massDiff2);
    
    if (spec1 == DEBUG_SCAN || spec2 == DEBUG_SCAN ||
        massDiff1 == DEBUG_MASS || massDiff2 == DEBUG_MASS) {
      DEBUG_MSG(spec1 << "  " << spec2 << "  " << 
                massDiff1 << "  " << massDiff2 << "  " << 
                pairs[iPair].score1 << "  " << pairs[iPair].score2);
    }

    if (massDiff1 >= minimumMass && massDiff1 <= maximumMass) {
      mapSpecToMassToScore[spec1][massDiff1] += pairs[iPair].score1;
    }
    if (massDiff2 >= minimumMass && massDiff2 <= maximumMass) {
      mapSpecToMassToScore[spec2][massDiff2] += pairs[iPair].score2;
    }
  }

  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[3];

  DEBUG_MSG("Saving file [" << outputFile << "]");
  ofstream ofs(outputFile.c_str());
  if (!ofs) {
    ERROR_MSG("Could not write to " << outputFile);
  }
  ofs << "scan\tmass\tscore" << endl;

  map<int, map<int, float> >::iterator itr1 = mapSpecToMassToScore.begin();
  map<int, map<int, float> >::iterator itrEnd1 = mapSpecToMassToScore.end();
  for (; itr1 != itrEnd1; itr1++) {
    int scan = itr1->first;

    map<int, float>::iterator itr2 = itr1->second.begin();
    map<int, float>::iterator itrEnd2 = itr1->second.end();
    for (; itr2 != itrEnd2; itr2++) {
      int mass = itr2->first;
      float score = itr2->second;
      ofs << scan << "\t"  << mass << "\t" << score << endl;
    }
  }

  return 0;
}

