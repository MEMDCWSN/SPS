//
//  main_find_supported_mods - Find the PTMs that are supported by pair data
//
#include "CommandLineParser.h"
#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"
#include "SpectrumPairSet.h"

#include <stdlib.h>

#define DEBUG_ 0


inline int convertStringToInt(const char * stringValue)
{
  int returnValue = -1;
  returnValue = atoi(stringValue);
  return returnValue;
}


using namespace std;
using namespace specnets;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope


  if (argc < 3) {
    cerr << "Usage: main_split_db <db_file> <num_splits> [options]" << endl;
    cerr << "     -count          only count (do not actually split)" << endl;
    return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("count", "COUNT_ONLY", 0));

  CommandLineParser clp(argc, argv, 2, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "Usage: main_split_db <db_file> <num_splits> [options]" << endl;
    cerr << "     -count          only count (do not actually split)" << endl;
    cerr << "Invalid options" << endl;
	  return -1;
  }

  ifstream ifs(argv[1]);

  int numProteins = 0;
  string lineBuff;
  while (ifs.good() && !ifs.eof()) {
    if (!getline(ifs, lineBuff)) { 
      break;
    }
    if (lineBuff[0] == '>') {
      numProteins++;
    }
  }
  ifs.close();
  DEBUG_VAR(numProteins);
  
  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  if (commandLineParams.exists("COUNT_ONLY")) {
    return 0;
  }

  int numSplits = convertStringToInt(argv[2]);
  DEBUG_VAR(numSplits);
  if (numSplits < 2) {
    ERROR_MSG("Number of files to split into must be > 1");
    return -2;
  }

  int numProtsPerFile = numProteins / numSplits;
  if (numProtsPerFile * numSplits < numProteins) {
    numProtsPerFile++;
  }
  DEBUG_VAR(numProtsPerFile);

  ifstream ifs2(argv[1]);
  
  string inputFile(argv[1]);
  string outputFileBase; 
  string outputFile; 
  size_t posDot = inputFile.find_last_of(".");
  if (posDot == string::npos) {
    outputFileBase = inputFile;
  } else {
    outputFileBase = inputFile.substr(0, posDot + 1);
  }
  
  int protCount = 1;
  int splitNum = 1;
  ofstream ofs;

  outputFile = outputFileBase +  + "_1.fasta";
  DEBUG_MSG("Opening file [" << outputFile << "]");
  ofs.open(outputFile.c_str());
  if (!ofs) {
    ERROR_MSG("Opening file [" << outputFile << "]");
    return -4;
  }

  while (ifs2.good() && !ifs2.eof()) {
    if (!getline(ifs2, lineBuff)) { 
      break;
    }
    if (lineBuff[0] == '>') {
      protCount++;
    }
    if (protCount > numProtsPerFile) {
      DEBUG_VAR(protCount);
      ofs.close();
      splitNum++;
      DEBUG_VAR(splitNum);
      protCount = 1;
      DEBUG_VAR(protCount);
      outputFile = outputFileBase + "_" + intToString(splitNum) + ".fasta";
      DEBUG_MSG("Opening file [" << outputFile << "]");
      ofs.open(outputFile.c_str());
      if (!ofs) {
        ERROR_MSG("Opening file [" << outputFile << "]");
        return -4;
      }
    }
    
    ofs << lineBuff << endl;

  }
  ofs.close();
  ifs.close();
  

  return 0;
}

