//
//  main_specdump - stand alone executable for dumping spectrum data
//
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

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;


// -------------------------------------------------------------------------
bool isSubPeptide(string pep1, string pep2)
{
  for (int i = 0; i < pep1.length(); i++) {
    if (pep1[i] == 'I') pep1[i] = 'L';
    if (pep1[i] == 'Q') pep1[i] = 'K';
  }
  for (int j = 0; j < pep2.length(); j++) {
    if (pep2[j] == 'I') pep2[j] = 'L';
    if (pep2[j] == 'Q') pep2[j] = 'K';
  }
  size_t pos = pep1.find(pep2);
  if (pos == string::npos) {
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------
void findModsInString(string & annotation, map<string,map<int,int> > & mapModCount, set<int> & setAllMods)
{
  //cout << annotation << endl;
  string modStringAA;
  string modString; 
  float modValue; 
  static string aminoAcids("(,)ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  bool startModAA = false;
  bool startModValue = false;
  for (int iChar = 0; iChar < annotation.length(); iChar++) {
    if (annotation[iChar] =='(') {
      startModAA = true;
    }
    else if (startModAA  && annotation[iChar] != ',') {
    	modStringAA += annotation[iChar];
    }
    else if (startModValue && aminoAcids.find_first_of(annotation[iChar]) == string::npos) {
    	modString += annotation[iChar];
    }
    else if (annotation[iChar] ==')') {
      startModValue = false;
      if (modString.empty()) {
      }
      //cout << modString << endl;
      sscanf(modString.c_str(), "%f", &modValue);
      int intModValue = 0;
      if (modValue < 0) {
        intModValue = int(modValue - 0.5);
      } else {
        intModValue = int(modValue + 0.5);
      }
      
      mapModCount[modStringAA][intModValue] = mapModCount[modStringAA][intModValue] + 1;
      setAllMods.insert(intModValue);
      //cout << modStringAA << "  " << intModValue << endl;
      modString.clear();
      modStringAA.clear();
    }
    else if (annotation[iChar] ==',') {
      startModAA = false;
      startModValue = true;
    }
  }
  return;
}

// -------------------------------------------------------------------------
void outputModMatrix(map<string,map<int,int> > & mapModCount, set<int> & setAllMods)
{
  cout << "   \t";
  map<string,map<int,int> >::iterator itrc1 = mapModCount.begin();
  map<string,map<int,int> >::iterator itrc_end1 = mapModCount.end();
  for ( ; itrc1 != itrc_end1; itrc1++) {
    map<int,int>::iterator itrc2 = itrc1->second.begin();
    map<int,int>::iterator itrc_end2 = itrc1->second.end();
    cout << itrc1->first << "\t";
  }
  cout << "Total" << endl;
  
  set<int>::iterator itrs = setAllMods.begin();
  set<int>::iterator itrs_end = setAllMods.end();
  for ( ; itrs != itrs_end; itrs++) {
    cout << *itrs;
    int rowTotal = 0;
    map<string,map<int,int> >::iterator itrc1 = mapModCount.begin();
    map<string,map<int,int> >::iterator itrc_end1 = mapModCount.end();
    for ( ; itrc1 != itrc_end1; itrc1++) {
      cout << "\t" << mapModCount[itrc1->first][*itrs];
      rowTotal += mapModCount[itrc1->first][*itrs];
    }
    cout << "\t" << rowTotal << endl;
  }
  return;
}

// -------------------------------------------------------------------------
string stripAnnotation(string & annotation)
{
  string cleanAnnotation; 
  static string aminoAcids("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  for (int iChar = 0; iChar < annotation.length(); iChar++)
  {
    if (aminoAcids.find_first_of(annotation[iChar]) != string::npos)
    {
      cleanAnnotation += annotation[iChar];
    }
  }
  return cleanAnnotation;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 6) {
    cerr << "Usage: main_scorecontig msgf_file db_file contig_psm cluster_table cluster_bin" << endl;
    return -1;
  }

  PeptideSpectrumMatchSet	psmSetMsgf;
  if (!psmSetMsgf.loadMSGFDBResultsFile(argv[1]), false) {
    ERROR_MSG("Loading MSGF file [" << argv[1] << "]");
    return -2;
  }

  //------------------------------------------------------------------------------------
  // Adrians scan number re-clustering code
  //------------------------------------------------------------------------------------
  ClusterData clusterInfo;

  // Specify path to project directory as INPUT_CLUSTERS_DIR (ie "svn/sps/qc/test_sps_silac")
  DEBUG_MSG("Loading clusters from " << argv[5]);
  if (!clusterInfo.loadData(argv[5])) {
    ERROR_MSG("Failed to load clusters from " << argv[5]);
    return false;
  }
  map<int, list<int> > clusterInfoMap;
  for (map<unsigned, list<pair<unsigned, unsigned> > >::iterator clustIt = clusterInfo.data.begin(); 
       clustIt != clusterInfo.data.end(); 
       clustIt++) {
    list<int> scans;

    // scans are always 1-based, these indices are 0-based
    for (list<pair<unsigned, unsigned> >::iterator scanIt = clustIt->second.begin(); 
         scanIt != clustIt->second.end(); 
         scanIt++) {
      scans.push_back(((int)scanIt->second) + 1);
    }
    clusterInfoMap[((int)clustIt->first) + 1] = scans;
  }
  // this method will assign each cluster the "consensus" PSM if one of its un-clustered spectra was identified
  psmSetMsgf.cluster(clusterInfoMap);

  //------------------------------------------------------------------------------------
  // now scan#s in m_specIDs will match 1-based indices of star spectra
  //------------------------------------------------------------------------------------


  // Save mostly for debug
  psmSetMsgf.saveToFile("msgf.psm");

  DB_fasta db;
  if (!db.Load(argv[2])) {
    ERROR_MSG("Loading DB file [" << argv[2] << "]");
    return -3;
  }

  PeptideSpectrumMatchSet	psmSetContig;
  if (!psmSetContig.loadFromFile(argv[3])) {
    ERROR_MSG("Loading contig PSM file [" << argv[3] << "]");
    return -4;
  }


  vector<vector<string> > linesCluster;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(argv[4],
                                               ";",
                                                "",
                                                linesCluster)) {
    ERROR_MSG("Loading table report file [" << argv[4] << "]");
    return -5;
  }

  map<int, int> mapContigToDbIndex;
  map<string,map<int,int> > mapModCount;
  set<int> setAllMods;
  for (int iPSM = 0; iPSM < psmSetContig.size(); iPSM++ ) {
    const psmPtr & psm = psmSetContig[iPSM];
    mapContigToDbIndex[psm->m_scanNum] = psm->m_dbIndex;
    findModsInString(psm->m_origAnnotation, mapModCount, setAllMods);
    //cout << psm->m_scanNum << "  " << psm->m_dbIndex << endl;
  }



  vector<float> match(8);

  map<int, tuple<int, int, int> > mapContigStats;  // match, mismatch, unknown (spectra)
  map<int, tuple<int, int, int> > mapReportScanToDb;
  map<int, int> mapScanCount;
  map<int, string> mapReportScanToAnnotation;
  map<int, int> mapScanToContig;
  
  static const int SCAN_INDEX = 0;
  static const int CONTIG_INDEX = 1;
  static const int ANNOTATION_INDEX = 4;
  cout << "linesCluster.size() = " << linesCluster.size() << endl;
  // First line is header but without '#'
  for (int iLine = 1; iLine < linesCluster.size(); iLine++ ) {
    int scanNumber;
    sscanf(linesCluster[iLine][SCAN_INDEX].c_str(), "%d", &scanNumber);
    int contigIndex;
    sscanf(linesCluster[iLine][CONTIG_INDEX].c_str(), "%d", &contigIndex);
    string & annotation = linesCluster[iLine][ANNOTATION_INDEX];

    mapContigStats[contigIndex] = make_tuple<int, int, int>(0, 0, 0);
    
    // Save for later viewing (debug) of annotation 
    mapReportScanToAnnotation[scanNumber]= annotation;
    mapScanToContig[scanNumber] = contigIndex; 
    //DEBUG_VAR(contigIndex);    
    mapScanCount[contigIndex]= mapScanCount[contigIndex] + 1;
                    
    if (!annotation.empty()) {
      int dbIndex = mapContigToDbIndex[contigIndex];
      string cleanAnnotation = stripAnnotation(annotation);
      if (cleanAnnotation.empty()) {
        match[0]++;
        continue;
      }
      string protein(db[dbIndex]);
      size_t posProtein = protein.find(cleanAnnotation);
      //cout << scanNumber << "  " << contigIndex <<  "  " << cleanAnnotation << endl;
      if (posProtein != string::npos) {
        //cout << scanNumber << "  " << cleanAnnotation << "  " << dbIndex <<  "  " << posProtein << endl;
        mapReportScanToDb[scanNumber] = 
                make_tuple<int, int, int>(dbIndex, posProtein, posProtein + cleanAnnotation.length() - 1);
      } else {
        cout << "Unable to find " << scanNumber << "  [" << cleanAnnotation << "] in [" << protein <<  "]" << endl;
      }
    } else {
      match[1]++;
    }
    
  } // for (int iLine = 0; iLine < linesCluster.size(); iLine++ )

  cout << "mapReportScanToDb.size() = " << mapReportScanToDb.size() << endl;
  cout << "match[0] = " << match[0] << endl;
  cout << "match[2] = " << match[2] << endl;

  map<int, tuple<int, int, int> > mapMsgfScanToDb;

  cout << "psmSetMsgf.size() = " << psmSetMsgf.size() << endl;
  for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ ) {

    const psmPtr & psm = psmSetMsgf[iPSM];
    string annotationMsgf = psm->m_annotation;
    string cleanAnnotationMsgf = stripAnnotation(annotationMsgf);

    int scanNumberMsgf = psm->m_scanNum;
    
    //cout << annotationMsgf << "  " << cleanAnnotationMsgf << endl;
    list<sps::tuple<int,float,string> > matches;
    int findResult = db.find(cleanAnnotationMsgf.c_str(), matches);
    //cout << findResult << "  " << matches.size() << endl;

    map<int, tuple<int, int, int> >::iterator itrFind = mapReportScanToDb.find(scanNumberMsgf);

    // Can we find this scanNumber in the report?
    if (itrFind == mapReportScanToDb.end()) {
      match[2]++;
      continue;
    }       

    list<tuple<int,float,string> >::iterator itr = matches.begin();
    list<tuple<int,float,string> >::iterator itr_end = matches.end();
    int iMatch = 0;
    for ( ; itr != itr_end; itr++) {
      iMatch++;
      //cout << itr->m0 << "  " << itr->m1 << "  " << itr->m2 << endl;
      //cout << db[itr->m0] << endl;
      int dbIndexMsgf = itr->m0;
      string proteinMsgf(db[dbIndexMsgf]);
      size_t posProteinMsgf = proteinMsgf.find(cleanAnnotationMsgf);
      
      if (posProteinMsgf == string::npos) {
        cout << "Unable to find " << psm->m_scanNum << "  [" << cleanAnnotationMsgf << "] in [" << proteinMsgf <<  "]" << endl;
        continue;
      }

      int startPosMsgf = posProteinMsgf;
      int endPosMsgf = posProteinMsgf + cleanAnnotationMsgf.length() - 1;
      
      //cout << psm->m_scanNum << "  " << cleanAnnotationMsgf << "  " << itr->m0 <<  "  " << posProtein << endl;
      mapMsgfScanToDb[psm->m_scanNum] = make_tuple<int, int, int>(itr->m0, startPosMsgf, endPosMsgf);
      
      int scanNumberReport = itrFind->first;
      int dbIndexReport = itrFind->second.m0;
      int startPosReport = itrFind->second.m1;
      int endPosReport = itrFind->second.m2;

      string cleanAnnotationReport = stripAnnotation(mapReportScanToAnnotation[scanNumberReport]);
      bool isSubString = isSubPeptide(cleanAnnotationMsgf, cleanAnnotationReport);
		
      if (dbIndexMsgf != dbIndexReport) {
        if (iMatch == matches.size()) {
          // This is last match so apparently none of them match db indexes
          // So lets set the DB Index to same as Report and look there
          dbIndexMsgf = itrFind->second.m0;
          string proteinMsgf(db[dbIndexMsgf]);
          // Now see if we can find the peptide there
          size_t posProteinMsgf = proteinMsgf.find(cleanAnnotationMsgf);
          if (posProteinMsgf == string::npos) {
            if (isSubString) {
              //DEBUG_VAR(scanNumberReport);    
              //DEBUG_VAR(mapScanToContig[scanNumberReport]);    
              mapContigStats[mapScanToContig[scanNumberReport]].m0++;
              match[4]++;
            } else {
              cout << mapScanToContig[scanNumberReport] << "  " << scanNumberReport << "  " 
            	   << cleanAnnotationMsgf << "  " << cleanAnnotationReport << endl;

              //DEBUG_VAR(scanNumberReport);    
              //DEBUG_VAR(mapScanToContig[scanNumberReport]);    
              mapContigStats[mapScanToContig[scanNumberReport]].m1++;
              match[3]++;
            }
            break;
          }
        } else {
          continue;
        }	
      }
      
      //cout << scanNumber << "  " << dbIndex << "  " << startPos <<  "  " << endPos  << "  " << itrFind->second.m1 << "  " << itrFind->second.m2 << endl;
      if (startPosReport >= startPosMsgf && endPosReport <= endPosMsgf) {
        mapContigStats[mapScanToContig[scanNumberReport]].m0++;
        match[7]++;
        break;
      }
      else if (endPosReport < startPosMsgf || startPosReport > endPosMsgf) {
        mapContigStats[mapScanToContig[scanNumberReport]].m0++;
        match[6]++;
        break;
      } else {
        mapContigStats[mapScanToContig[scanNumberReport]].m0++;
        match[5]++;
        break;
      }

    } // for ( ; itr != itr_end; itr++) matches

  } // for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ )

  cout << "mapMsgfScanToDb.size() = " << mapMsgfScanToDb.size() << endl;
  cout << endl;
  cout << "No anno in report      : " << match[0] << endl;
  cout << "No AA's in annotation  : " << match[1] << endl;
  cout << "Scan not in report     : " << match[2] << endl;
  cout << "Peptide not in protein : " << match[3] << endl;
  cout << "Peptide is sub of msgf : " << match[4] << endl;
  cout << "Same Peptide Overlap   : " << match[5] << endl;
  cout << "Same Peptide Diff Loc  : " << match[6] << endl;
  cout << "Same Peptide Same Loc  : " << match[7] << endl;
  cout << endl;


  map<int, tuple<int,int,int> >::iterator itr2 = mapReportScanToDb.begin();
  map<int, tuple<int,int,int> >::iterator itr2_end = mapReportScanToDb.end();
  for ( ; itr2 != itr2_end; itr2++) {
	map<int, tuple<int,int,int> >::iterator itr2_find = mapMsgfScanToDb.find(itr2->first); 
	if (itr2_find == mapMsgfScanToDb.end()) {
	  mapContigStats[mapScanToContig[itr2->first]].m2++;
	}
  }

  map<int, tuple<int,int,int> >::iterator itr = mapContigStats.begin();
  map<int, tuple<int,int,int> >::iterator itr_end = mapContigStats.end();
  int iMatch = 0;
  cout << "Contig\t#spec\tmatch\tmis\tunk\t%match\t%mis\t%unk" << endl;
  for ( ; itr != itr_end; itr++) {
    cout << itr->first << "\t" << mapScanCount[itr->first] << "\t";
    cout << itr->second.m0 << "\t" << itr->second.m1 << "\t" << itr->second.m2 << "\t";
    float total = itr->second.m0 + itr->second.m1 + itr->second.m2;
    float percentMatch = 0.0;
    float percentMismatch = 0.0;
    float percentUnknown = 0.0;
    if (total != 0) { 
      percentMatch = ((float)itr->second.m0 / total) * 100.0;
      percentMismatch = ((float)itr->second.m1 / total) * 100.0;
      percentUnknown = ((float)itr->second.m2 / total) * 100.0;
    } 
    cout << percentMatch << "\t" << percentMismatch << "\t" << percentUnknown << "\t" << endl;
    float denom = percentMatch + percentMismatch;
  }
  
  cout << endl;

  outputModMatrix(mapModCount, setAllMods);

  
  return 0;
}
