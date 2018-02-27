//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "abruijn.h"
#include "ClusterData.h"
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ReportTableClusterConsensus.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include "tuple.h"

#include <stdlib.h>

#define DEBUG_ 0
#define DEBUG_PERCENT_MATCH 0
#define DEBUG_TAG_MATCH 0
#define DEBUG_CONTIG_INDEX 0
#define DEBUG_GAPS 0
#define DEBUG_GAP_NUMBER 3

#define MATCH_THRESHOLD 0.9
#define SIM_MATCH_THRESHOLD 0.8

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;


struct Comparison {
  string cleanAnno1;
  string cleanAnno2;
  string cleanAnno3;
  int   matchClean1v2;
  int   matchClean1v3;
  int   matchClean2v3;
  int   matchMods1v2;
  int   matchMods1v3;
  int   matchMods2v3;
  int   matchLocs1v2;
  int   matchLocs1v3;
  int   matchLocs2v3;

  Comparison() {
    cleanAnno1 = "---";
    cleanAnno2 = "---";
    cleanAnno3 = "---";
    matchClean1v2 = 0;
    matchClean1v3 = 0;
    matchClean2v3 = 0;
    matchMods1v2 = 0;
    matchMods1v3 = 0;
    matchMods2v3 = 0;
    matchLocs1v2 = 0;
    matchLocs1v3 = 0;
    matchLocs2v3 = 0;
  }
};

// -------------------------------------------------------------------------
string convertAnnotation(string annotation)
{
  string convertedAnnotation; 
  for (int iChar = 0; iChar < annotation.length(); iChar++) {
    if (annotation[iChar] == 'I') {
      convertedAnnotation += 'L';
    } else if (annotation[iChar] == 'Q') {
      convertedAnnotation += 'K';
    } else {
      convertedAnnotation += annotation[iChar];
    }
  }
  return convertedAnnotation;
}

// -------------------------------------------------------------------------
float getPercentMatch(string & string1, string & string2)
{
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(string1.size());
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(string2.size());
  int dynaArray[string1.size()+1][string2.size()+1];
  for (int i = 0; i < string1.size() + 1; i++) {
    for (int j = 0; j < string2.size() + 1; j++) {
      dynaArray[i][j] = 0;
    }
  }
  for (int i = 1; i <= string1.size(); i++) {
    for (int j = 1; j <= string2.size(); j++) {
      int up = dynaArray[i-1][j];
      int left = dynaArray[i][j-1];
      int diag = dynaArray[i-1][j-1];
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(i << "  " << j);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(up << "  " << left << "  " << diag);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(string1[i-1] << "  " << string2[j-1]);
      if (string1[i-1] == string2[j-1]) {
        diag++;
      }
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(up << "  " << left << "  " << diag);
      int max1 = max(left, up);
      dynaArray[i][j] = max(max1, diag);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(dynaArray[i][j]);
    }
  }

  if (DEBUG_PERCENT_MATCH) {
    cout << string1 << "  " << string2 << endl;
    for (int i = 0; i <= string1.size(); i++) {
      for (int j = 0; j <= string2.size(); j++) {
        cout << dynaArray[i][j] << "  ";
      }
      cout << endl;
    }
    cout << endl;
  }
  int match = dynaArray[string1.size()][string2.size()];
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(match);
  float percent = (float)match / (float)min(string1.size(), string2.size());
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(percent);
  percent = min(percent, (float)1.0);
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(percent);
  return percent;
}

// -------------------------------------------------------------------------
bool isSameMods(vector<float> & mods1, vector<float> & mods2)
{
  if (mods1.size() != mods2.size()) {
    return false;
  }
  sort(mods1.begin(), mods1.end());
  sort(mods2.begin(), mods2.end());
  for (int i = 0; i < mods1.size(); i++) {
    if ( (int)(mods1[i] + 0.5) != (int)(mods2[i] + 0.5) ) {
      return false;
    }
  }
  return true;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 3 && argc != 4) {
    cerr << "Usage: main_comparepsms psm_file1 psm_file2 [psm_file2]" << endl;
    return -1;
  }
  
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

  PeptideSpectrumMatchSet	psmSet3;
  if (argc == 4 && !psmSet3.loadFromFile(argv[3])) {
    ERROR_MSG("Loading PSM file [" << argv[3] << "]");
    return -4;
  }
  DEBUG_VAR(psmSet3.size());
  
  //----------------------------------------------------------
  // Create the PSM maps
  //----------------------------------------------------------
  map<int, psmPtr> mapPsm1;
  for (int iPSM = 0; iPSM < psmSet1.size(); iPSM++ ) {
    psmPtr psmSpec = psmSet1[iPSM];
    mapPsm1[psmSpec->m_scanNum] = psmSpec;
    //DEBUG_VAR(psmSpec->m_scanNum);
  }
  DEBUG_VAR(mapPsm1.size());
  map<int, psmPtr> mapPsm2;
  for (int iPSM = 0; iPSM < psmSet2.size(); iPSM++ ) {
    psmPtr psmSpec = psmSet2[iPSM];
    mapPsm2[psmSpec->m_scanNum] = psmSpec;
    //DEBUG_VAR(psmSpec->m_scanNum);
  }
  DEBUG_VAR(mapPsm2.size());
  map<int, psmPtr> mapPsm3;
  for (int iPSM = 0; iPSM < psmSet3.size(); iPSM++ ) {
    psmPtr psmSpec = psmSet3[iPSM];
    mapPsm3[psmSpec->m_scanNum] = psmSpec;
    //DEBUG_VAR(psmSpec->m_scanNum);
  }
  DEBUG_VAR(mapPsm3.size());

  //----------------------------------------------------------
  // Compare the PSMs
  //----------------------------------------------------------
  map<int, Comparison>  mapScanCompare;

  map<int, psmPtr>::iterator itrp = mapPsm1.begin();
  map<int, psmPtr>::iterator itrpEnd = mapPsm1.end();
  for ( ; itrp != itrpEnd; itrp++) {
    psmPtr psmSpec1 = itrp->second;
    int scanNum = psmSpec1->m_scanNum;

    string convertAnno1 = convertAnnotation(psmSpec1->m_annotation);
    string cleanAnno1;
    PeptideSpectrumMatch::getUnmodifiedPeptide(convertAnno1, cleanAnno1);
    mapScanCompare[scanNum].cleanAnno1 = cleanAnno1;
    //DEBUG_VAR(cleanAnno1);

    if (mapPsm2.find(scanNum) == mapPsm2.end()) {
      continue;
    }
    psmPtr psmSpec2 = mapPsm2[scanNum];
    string convertAnno2 = convertAnnotation(psmSpec2->m_annotation);
    string cleanAnno2;
    PeptideSpectrumMatch::getUnmodifiedPeptide(convertAnno2, cleanAnno2);
    mapScanCompare[scanNum].cleanAnno2 = cleanAnno2;
    //DEBUG_VAR(cleanAnno2);

    if (cleanAnno1 == cleanAnno2) {
      mapScanCompare[scanNum].matchClean1v2 = 1;
    }

    vector<float> modifications1;
    vector<unsigned int> positions1;
    vector<unsigned int> lengths1;
    psmSpec1->getModificationsAndPositions(modifications1, positions1, lengths1);

    vector<float> modifications2;
    vector<unsigned int> positions2;
    vector<unsigned int> lengths2;
    psmSpec2->getModificationsAndPositions(modifications2, positions2, lengths2);

    mapScanCompare[scanNum].matchMods1v2 = isSameMods(modifications1, modifications2);
    mapScanCompare[scanNum].matchLocs1v2 = (mapScanCompare[scanNum].matchMods1v2 && (positions1 == positions2));
  }
  DEBUG_VAR(mapScanCompare.size());

  itrp = mapPsm2.begin();
  itrpEnd = mapPsm2.end();
  for ( ; itrp != itrpEnd; itrp++) {
    psmPtr psmSpec2 = itrp->second;
    int scanNum = psmSpec2->m_scanNum;

    // If there is a PSM in set 1 for this scan.. we already did it above
    if (mapPsm1.find(scanNum) != mapPsm1.end()) {
      continue;
    }

    string convertAnno2 = convertAnnotation(psmSpec2->m_annotation);
    string cleanAnno2;
    PeptideSpectrumMatch::getUnmodifiedPeptide(convertAnno2, cleanAnno2);
    mapScanCompare[scanNum].cleanAnno2 = cleanAnno2;
    DEBUG_VAR(cleanAnno2);
  }
  DEBUG_VAR(mapScanCompare.size());


  itrp = mapPsm3.begin();
  itrpEnd = mapPsm3.end();
  for ( ; itrp != itrpEnd; itrp++) {
    psmPtr psmSpec3 = itrp->second;
    int scanNum = psmSpec3->m_scanNum;

    string convertAnno3 = convertAnnotation(psmSpec3->m_annotation);
    string cleanAnno3;
    PeptideSpectrumMatch::getUnmodifiedPeptide(convertAnno3, cleanAnno3);
    mapScanCompare[scanNum].cleanAnno3 = cleanAnno3;

    vector<float> modifications3;
    vector<unsigned int> positions3;
    vector<unsigned int> lengths3;
    psmSpec3->getModificationsAndPositions(modifications3, positions3, lengths3);

    if (mapPsm1.find(scanNum) != mapPsm1.end()) {
      if (mapScanCompare[scanNum].cleanAnno1 == cleanAnno3) {
        mapScanCompare[scanNum].matchClean1v3 = 1;
      }
      psmPtr psmSpec1 = mapPsm1[scanNum];
      vector<float> modifications1;
      vector<unsigned int> positions1;
      vector<unsigned int> lengths1;
      psmSpec1->getModificationsAndPositions(modifications1, positions1, lengths1);
      mapScanCompare[scanNum].matchMods1v3 = isSameMods(modifications1, modifications3);
      mapScanCompare[scanNum].matchLocs1v3 = (mapScanCompare[scanNum].matchMods1v3 && (positions1 == positions3));
    }

    if (mapPsm2.find(scanNum) != mapPsm2.end()) {
      if (mapScanCompare[scanNum].cleanAnno2 == cleanAnno3) {
        mapScanCompare[scanNum].matchClean2v3 = 1;
      }
      psmPtr psmSpec2 = mapPsm2[scanNum];
      vector<float> modifications2;
      vector<unsigned int> positions2;
      vector<unsigned int> lengths2;
      psmSpec2->getModificationsAndPositions(modifications2, positions2, lengths2);
      mapScanCompare[scanNum].matchMods2v3 = isSameMods(modifications2, modifications3);
      mapScanCompare[scanNum].matchLocs2v3 = (mapScanCompare[scanNum].matchMods2v3 && (positions2 == positions3));
    }
  }
  DEBUG_VAR(mapScanCompare.size());
  
  int totalClean1v2 = 0;
  int totalClean1v3 = 0;
  int totalClean2v3 = 0;
  int totalMods1v2 = 0;
  int totalMods1v3 = 0;
  int totalMods2v3 = 0;
  int totalLocs1v2 = 0;
  int totalLocs1v3 = 0;
  int totalLocs2v3 = 0;

  //----------------------------------------------------------
  // Output results
  //----------------------------------------------------------
  if (argc == 4) {
    cout << "scan#\tanno1\tanno2\tanno3\tscore1\tpvalue1\tscore2\tpvalue2\tscore3\tpvalue3\tclean1v2\tclean1v3\tclean2v3\tmods1v2\tmods1v3\tmods2v3\tlocs1v2\tlocs1v3\tlocs2v3" << endl;
  } else {
    cout << "scan#\tanno1\tanno2\tscore1\tpvalue1\tscore2\tpvalue2\tclean1v2\tmods1v2\tlocs1v2" << endl;
  }

  map<int, Comparison>::iterator itrm = mapScanCompare.begin();
  map<int, Comparison>::iterator itrmEnd = mapScanCompare.end();
  for ( ; itrm != itrmEnd; itrm++) {
    Comparison & comparison = itrm->second;
    int scanNum = itrm->first;
    cout << scanNum << "\t";
    
    if (mapPsm1.find(scanNum) != mapPsm1.end()) {
      cout << mapPsm1[scanNum]->m_annotation << "\t";
    } else {
      cout << "---" << "\t";
    }
    if (mapPsm2.find(scanNum) != mapPsm2.end()) {
      cout << mapPsm2[scanNum]->m_annotation << "\t";
    } else {
      cout << "---" << "\t";
    }

    if (argc == 4) {
      if (mapPsm3.find(scanNum) != mapPsm3.end()) {
        cout << mapPsm3[scanNum]->m_annotation << "\t";
      } else {
        cout << "---" << "\t";
      }
    }

    if (mapPsm1.find(scanNum) != mapPsm1.end()) {
      cout << mapPsm1[scanNum]->m_score << "\t";
      cout << mapPsm1[scanNum]->m_pValue << "\t";
    } else {
      cout << "-1.0" << "\t";
      cout << "-1.0" << "\t";
    }

    if (mapPsm2.find(scanNum) != mapPsm2.end()) {
      cout << mapPsm2[scanNum]->m_score << "\t";
      cout << mapPsm2[scanNum]->m_pValue << "\t";
    } else {
      cout << "-1.0" << "\t";
      cout << "-1.0" << "\t";
    }

    if (argc == 4) {
      if (mapPsm3.find(scanNum) != mapPsm3.end()) {
        cout << mapPsm3[scanNum]->m_score << "\t";
        cout << mapPsm3[scanNum]->m_pValue << "\t";
      } else {
        cout << "-1.0" << "\t";
        cout << "-1.0" << "\t";
      }
    }

    cout << comparison.matchClean1v2 << "\t";
    if (argc == 4) {
      cout << comparison.matchClean1v3 << "\t";
      cout << comparison.matchClean2v3 << "\t";
    }
    cout << comparison.matchMods1v2 << "\t";
    if (argc == 4) {
      cout << comparison.matchMods1v3 << "\t";
      cout << comparison.matchMods2v3 << "\t";
    }
    cout << comparison.matchLocs1v2 << "\t";
    if (argc == 4) {
      cout << comparison.matchLocs1v3 << "\t";
      cout << comparison.matchLocs2v3 << "\t";
    }

    if (comparison.matchClean1v2) totalClean1v2++;
    if (comparison.matchClean1v3) totalClean1v3 ++;
    if (comparison.matchClean2v3) totalClean2v3 ++;
    if (comparison.matchMods1v2) totalMods1v2 ++;
    if (comparison.matchMods1v3) totalMods1v3 ++;
    if (comparison.matchMods2v3) totalMods2v3 ++;
    if (comparison.matchLocs1v2) totalLocs1v2 ++;
    if (comparison.matchLocs1v3) totalLocs1v3 ++;
    if (comparison.matchLocs2v3) totalLocs2v3 ++;

    cout << endl;
  }

  cout << endl;
  cout << "Total 1:\t" << psmSet1.size() << endl;
  cout << "Total 2:\t" << psmSet2.size() << endl;
  if (argc == 4) {
    cout << "Total 3:\t" << psmSet3.size() << endl;
  }
  cout << "Match Clean 1v2:\t" << totalClean1v2 << endl;
  if (argc == 4) {
    cout << "Match Clean 1v3:\t" << totalClean1v3 << endl;
    cout << "Match Clean 2v3:\t" << totalClean2v3 << endl;
  }
  cout << "Match Mods 1v2:\t" << totalMods1v2 << endl;
  if (argc == 4) {
    cout << "Match Mods 1v3:\t" << totalMods1v3 << endl;
    cout << "Match Mods 2v3:\t" << totalMods2v3 << endl;
  }
  cout << "Match Locs 1v2:\t" << totalLocs1v2 << endl;
  if (argc == 4) {
    cout << "Match Locs 1v3:\t" << totalLocs1v3 << endl;
    cout << "Match Locs 2v3:\t" << totalLocs2v3 << endl;
  }
  
  return 0;
}
