//
//  main_find_supported_mods - Find the PTMs that are supported by pair data
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
int roundedMass(float mass) 
{
  return int(abs(mass) + 0.5) * (mass < 0 ? -1 : 1);
}


// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope


  if (argc < 5) {
    cerr << "Usage: main_find_supported_mods <psm_file> <pair_file> <star_file> <out_file> [options]" << endl;
    cerr << "   -outdir <directoryu>  Output directory for output files (if none specified current directory will be used)" << endl;
    cerr << "   -useorig              Use original annotation" << endl;
    return -1;
  }
  
  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("outdir", "RESULTS_OUTPUT_DIR", 1));
  listOptions.push_back(CommandLineParser::Option("useorig", "USE_ORIGINAL_ANNO", 0));

  CommandLineParser clp(argc, argv, 4, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "Usage: main_find_supported_mods <psm_file> <pair_file> <star_file> <out_file> [options]" << endl;
    cerr << "   -outdir <directoryu>  Output directory for output files (if none specified current directory will be used)" << endl;
    cerr << "   -useorig              Use original annotation" << endl;
    cerr << "Invalid options" << endl;
	  return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);
  
  PeptideSpectrumMatchSet	psmSet;
  if (!psmSet.loadFromFile(argv[1])) {
    ERROR_MSG("Loading PSM file [" << argv[1] << "]");
    return -1;
  }

  SpectrumPairSet pairs;
  SpecSet starSpectra;

  string pairsFileName = argv[2];
  DEBUG_MSG("Loading pairs file [" << pairsFileName << "]...");
  if (!pairs.loadFromBinaryFile(pairsFileName.c_str())) {
    ERROR_MSG("Could not load " << pairsFileName);
    return false;
  }

  string starsFileName = argv[3];
  DEBUG_MSG("Loading stars file [" << starsFileName << "]...");
  if (!starSpectra.loadPklBin(starsFileName.c_str())) {
    ERROR_MSG("Could not load " << starsFileName);
    return false;
  }

  bool useOrig = commandLineParams.exists("USE_ORIGINAL_ANNO");
  DEBUG_VAR(useOrig);

  map<int, float> mapSpecToMass;
  for (int iStar = 0; iStar < starSpectra.size(); iStar++) {
    mapSpecToMass[starSpectra[iStar].scan] = starSpectra[iStar].parentMass;
  }

  //--------------------------------------------------------------------------
  // Find the mods that are supported by an annotation
  //--------------------------------------------------------------------------
  map<int, string> mapSpecToPvrNum;
  map<int,int> mapScanToIndex;
  map<int, set<int> > annotatedMods;
  map<int,float> mapScanToPvalue;
  for (int iPsm = 0; iPsm < psmSet.size(); iPsm++) {
    int scan = psmSet[iPsm]->m_scanNum;
    mapScanToIndex[scan] = iPsm;
    mapSpecToPvrNum[psmSet[iPsm]->m_scanNum] = psmSet[iPsm]->m_peptideRegionGroup;
    mapScanToPvalue[psmSet[iPsm]->m_scanNum] = psmSet[iPsm]->m_pValue;

    vector<float> modifications;
    vector<unsigned int> positions;
    vector<unsigned int> lengths;
    psmSet[iPsm]->getModificationsAndPositions(modifications, positions, lengths, useOrig);
    bool verified = false;
    int totalMods = 0;
    for (int iMod = 0; iMod < modifications.size(); iMod++) {
      float intModValue = roundedMass(modifications[iMod]);
      annotatedMods[scan].insert(intModValue);
    }
  }


  //--------------------------------------------------------------------------
  // Sum all the scores for all the mods for all spectra
  //    Scans must be in the same PVR
  //    Mod must be supported by the annotation for the PSM
  //--------------------------------------------------------------------------
  map<int, map<int, float> > mapSpecToMassToScore;
  map<int, float> mapMassToBestScore;
  map<int, pair<float, int> > mapMassToBestPair;
  DEBUG_VAR(pairs.size()); 
  for (int iPair = 0; iPair < pairs.size(); iPair++) {
    int spec1 = pairs[iPair].spec1 + 1;  // These are 0 based
    int spec2 = pairs[iPair].spec2 + 1;
    
    // First check to see if the two are in the same PVR
    string pvr1 = mapSpecToPvrNum[spec1];
    string pvr2 = mapSpecToPvrNum[spec2];
    if (pvr1 != pvr2) {
      continue;
    }

    float rawDiff1 = mapSpecToMass[spec1] - mapSpecToMass[spec2];
    if (DEBUG_) DEBUG_VAR(rawDiff1);

    // not interested in 0 mass diffs
    if (abs(rawDiff1) < 0.5) {
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
    
    if (pairs[iPair].score1 > mapMassToBestScore[massDiff1]) {
      mapMassToBestScore[massDiff1] = pairs[iPair].score1;
      mapMassToBestPair[massDiff1] = make_pair<int, int>(spec1, spec2);
    }

    if (pairs[iPair].score2 > mapMassToBestScore[massDiff2]) {
      mapMassToBestScore[massDiff2] = pairs[iPair].score2;
      mapMassToBestPair[massDiff2] = make_pair<int, int>(spec2, spec1);
    }

    // Only consider those mass diffs which are supported by annotation
    if (annotatedMods[spec1].find(massDiff1) != annotatedMods[spec1].end()) {
      if (massDiff1 != 0) {
        mapSpecToMassToScore[spec1][massDiff1] += pairs[iPair].score1;
        if (spec1 == DEBUG_SCAN || massDiff1 == DEBUG_MASS) {
          DEBUG_MSG(spec1 << "  " << massDiff1 << "  " << mapSpecToMassToScore[spec1][massDiff1]);
        }
      }
    }
    // Only consider those mass diffs which are supported by annotation
    if (annotatedMods[spec2].find(massDiff2) != annotatedMods[spec2].end()) {
      if (massDiff2 != 0) {
        mapSpecToMassToScore[spec2][massDiff2] += pairs[iPair].score2;
        if (spec2 == DEBUG_SCAN || massDiff2 == DEBUG_MASS) {
          DEBUG_MSG(spec2 << "  " << massDiff2 << "  " << mapSpecToMassToScore[spec2][massDiff2]);
        }
      }
    }
  }

  // Rearrange for output
  map<int, set<pair<float,int> > > mapMassToBestScan;  // mass to score/scan
  map<int, map<int, float> >::iterator itr1 = mapSpecToMassToScore.begin();
  map<int, map<int, float> >::iterator itrEnd1 = mapSpecToMassToScore.end();
  for (; itr1 != itrEnd1; itr1++) {
    int scan = itr1->first;

    map<int, float>::iterator itr2 = itr1->second.begin();
    map<int, float>::iterator itrEnd2 = itr1->second.end();
    for (; itr2 != itrEnd2; itr2++) {
      int mass = itr2->first;
      float score = itr2->second;

      if (mass == DEBUG_MASS) {
        DEBUG_MSG(mass << "  " << scan << "  " << score);
      }
    
      mapMassToBestScan[mass].insert(make_pair(score, scan));
    }  
  }
  
  string outputFile;
  if (commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    outputFile += commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outputFile += "/";
  }
  outputFile += argv[4];

#if 0
  DEBUG_TRACE;
  DEBUG_MSG("Saving best mods to file [" << outputFile << "]");
  psmSetBest.saveToFile(outputFile.c_str(), true, true);

#else
  DEBUG_TRACE;
  ofstream ofs(outputFile.c_str());
  if (!ofs) {
    ERROR_MSG("Could not write to " << outputFile);
  }
  ofs << "mass\trank\tscan\tannotation\tscore\tpvalue\tbest_score\tbest_pair_1\tbest_pair_2" << endl;

  map<int, set<pair<float,int> > >::iterator itrBest = mapMassToBestScan.begin();
  map<int, set<pair<float,int> > >::iterator itrBestEnd = mapMassToBestScan.end();
  for (; itrBest != itrBestEnd; itrBest++) {
    int mass = itrBest->first;
    set<pair<float,int> >::iterator itrSet = itrBest->second.begin();
    set<pair<float,int> >::iterator itrSetEnd = itrBest->second.end();
    int setSize = itrBest->second.size();
    for (int rank = setSize; itrSet != itrSetEnd; itrSet++, rank--) {
      int scan = (*itrSet).second;
      float score = (*itrSet).first;
      int index = mapScanToIndex[scan];
      ofs << mass << "\t"  << rank << "\t" << scan << "\t";
      if (useOrig) {
        ofs << psmSet[index]->m_origAnnotation << "\t";
      } else {
        ofs << psmSet[index]->m_annotation << "\t";
      }
      ofs << score << "\t"
            << mapScanToPvalue[scan] << "\t"
            << mapMassToBestScore[mass] << "\t" 
            << mapMassToBestPair[mass].first << "\t" 
            << mapMassToBestPair[mass].second << endl;
    }
  }
#endif

  return 0;
}

