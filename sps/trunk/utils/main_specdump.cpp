//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "CommandLineParser.h"
#include "Logger.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include <stdlib.h>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
void displayUsage(void)
{
  cerr << "Usage: main_specdump <type> <specfile> [options]" << endl;
  cerr << "       Valid types are: mgf, prms, specset, specpairset, pklbin" << endl;
  cerr << "   -size               dump only the header info and size" << endl;
  cerr << "   -index N            dump only the single spectrum with index = N" << endl;
  cerr << "   -scan N             dump only the single spectrum with scan number = N" << endl;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  DEBUG_TRACE;
  if (argc < 3) {
    displayUsage();
    return -1;
  }

  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("scan", "SINGLE_SCAN", 1));
  listOptions.push_back(CommandLineParser::Option("index", "SINGLE_INDEX", 1));
  listOptions.push_back(CommandLineParser::Option("size", "SIZE_ONLY", 0));

  CommandLineParser clp(argc, argv, 2, listOptions);
  string parserError = "";
  if (!clp.validate(parserError)) {
    displayUsage();
    cerr << "Invalid options" << endl;
    return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  
  string type = argv[1];

  if (type == "mgf" || type == "prms" || type == "specset" || type == "pklbin")
  {
    SpecSet spectra1;
    size_t size1 = 0;


    if (type == "mgf")
    {
      DEBUG_TRACE;
      spectra1.LoadSpecSet_mgf(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
    }
    else if (type == "prms")
    {
      DEBUG_TRACE;
      spectra1.LoadSpecSet_prmsv3(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
    }
    else if (type == "specset" || type == "pklbin")
    {
      DEBUG_TRACE;
      spectra1.loadPklBin(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
    }

    if (commandLineParams.exists("SIZE_ONLY")) {
      cout << "Total Spectra = " << size1 << endl;
      return 0;
    }

    for (size_t i = 0; i < size1; i++)
    {
      if (commandLineParams.exists("SINGLE_INDEX") &&
          commandLineParams.getValueInt("SINGLE_INDEX", -1) != i) {
        continue;
      }
      
      if (commandLineParams.exists("SINGLE_SCAN") &&
          commandLineParams.getValueInt("SINGLE_SCAN", -1) != spectra1[i].scan) {
        continue;
      }
      
      cout << "i = " << i << endl;
      
      //DEBUG_VAR(i);
      size_t peakSize1 = spectra1[i].size();
      //DEBUG_VAR(peakSize1);

      cout << "Parent Mass = " << spectra1[i].parentMass << endl;
      cout << "Scan Number = " << spectra1[i].scan << endl;
      cout << "MS Level = " << spectra1[i].msLevel << endl;
      cout << "Charge = " << spectra1[i].parentCharge << endl;

      for (size_t j = 0; j < peakSize1; j++)
      {
        TwoValues<float> peak1 = spectra1[i][j];
        //DEBUG_VAR(peak1.values[0]);
        //DEBUG_VAR(peak1.values[1]);
        cout << peak1[0] << ", " << peak1[1] << endl;
      }  
      // If we got here and only doing one spectrum we can exit the loop now
      if (commandLineParams.exists("SINGLE_INDEX") ||
          commandLineParams.exists("SINGLE_SCAN")) {
        break;
      }
    }
  }
  else if (type == "specpairset")
  {
    SpectrumPairSet specpairset1;
    specpairset1.loadFromBinaryFile(argv[2]);
    size_t size1 = specpairset1.size();
    DEBUG_VAR(size1);

    if (commandLineParams.exists("SIZE_ONLY")) {
      cout << "Total Spectra = " << size1 << endl;
      return 0;
    }

    for (size_t i = 0; i < size1; i++)
    {
      SpectrumPair pair1 = specpairset1[i];

      if (commandLineParams.exists("SINGLE_INDEX") &&
          commandLineParams.getValueInt("SINGLE_INDEX", -1) != i) {
        continue;
      }
      // In this case only show those pairs that contain the specified scan
      if (commandLineParams.exists("SINGLE_SCAN") &&
          commandLineParams.getValueInt("SINGLE_SCAN", -1) != pair1.spec1 &&
          commandLineParams.getValueInt("SINGLE_SCAN", -1) != pair1.spec2) {
        continue;
      }
      cout << "i =" << i << endl;
      cout << "pair1.spec1 = " << pair1.spec1 << endl;
      cout << "pair1.spec2 = " << pair1.spec2 << endl;
      cout << "pair1.score1 = " << pair1.score1 << endl;
      cout << "pair1.score2 = " << pair1.score2 << endl;
      cout << "pair1.shift1 = " << pair1.shift1 << endl;
      cout << "pair1.shift2 = " << pair1.shift2 << endl;
      cout << "pair1.specC = " << pair1.specC << endl;
      cout << "pair1.spec2rev = " << pair1.spec2rev << endl;

      // If we got here and only doing one index we can exit the loop now
      if (commandLineParams.exists("SINGLE_INDEX")) {
        break;
      }
    }
  }
  else
  {
    displayUsage();
    cerr << "Unknown type" << endl;
    return -1;
  }  
  
  return 0;
}

