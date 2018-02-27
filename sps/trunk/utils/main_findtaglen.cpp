//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "abruijn.h"
#include "clusters.h"
#include "ClusterData.h"
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ReportTableClusterConsensus.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include "tags.h"
#include "tuple.h"

#include <stdlib.h>

using namespace std;
using namespace sps;
using namespace specnets;


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

  if (argc != 5 && argc != 8) {
    cerr << "Usage: main_scorecontig contig_psm_file contig_score_file cluster_file contig_spec_file taglength taggap numtags" << endl;
    return -1;
  }

  PeptideSpectrumMatchSet	psmSetContig;
  if (!psmSetContig.loadFromFile(argv[1])) {
    ERROR_MSG("Loading MSGF PSM file [" << argv[1] << "]");
    return -2;
  }

  map<int, string> mapContigAnno;
  for (int i = 0; i < psmSetContig.size(); i++) {
    mapContigAnno[psmSetContig[i]->m_scanNum] = psmSetContig[i]->m_origAnnotation;
  }

  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(argv[2], "\t", "", lines)) {
      ERROR_MSG("Unable to open contig score file! " << argv[2]);
      return -3;
  }

  map<int, int> mapContigCorrect;
  for (int i = 1; i < lines.size(); i++) {  // skip header line
    int correct;
    sscanf(lines[i][0].c_str(), "%d", &correct);
    int contigNum;
    sscanf(lines[i][1].c_str(), "%d", &contigNum);
    mapContigCorrect[contigNum] = correct;
  }
  //DEBUG_TRACE;

//  if (contigs.Load("assembly/path_spectra_as_cluster.txt",
//                   "assembly/sps_seqs.pklbin") <= 0) {

  Clusters contigs;
  if (contigs.Load(argv[3], argv[4]) <= 0) {
    ERROR_MSG("Problem loading contig information");
    exit(-4);
  }

  unsigned int maxNumTags = 100;
  float peakTol = 0.5;
  list<Tag> tags;
  AAJumps jumps(1);

  map<int, int>  mapTagDepthCount;

if (argc == 8) {
  int tagLen = 7;
  int numJumps = 0;
  sscanf(argv[5], "%d", &tagLen);
  sscanf(argv[6], "%d", &numJumps);
  sscanf(argv[7], "%d", &maxNumTags);

  for (int i = 0; i < contigs.consensus.size(); i++) {
    int scan = i+1;
    int correct = mapContigCorrect[scan];
    string anno = mapContigAnno[scan];
    for (int j = 0; j < anno.length(); j++) {
      if (anno[j] == 'L') anno[j] = 'I';
    }

    ExtractTags(contigs.consensus[i], tags, peakTol, tagLen, numJumps, maxNumTags);
    if (correct == 1) {
      int tagNum = 1;
      list<Tag>::iterator itr = tags.begin();
      list<Tag>::iterator itr_end = tags.end();
      for (; itr != itr_end; itr++) {
        string tagSeq;
        for (int s = 0; s < itr->sequence.size(); s++) {
          int index = itr->sequence[s];
          tagSeq += jumps.aaLetters[(char)index];
        }
        size_t pos = anno.find(tagSeq);
        if (pos != string::npos) {
          mapTagDepthCount[tagNum]++;
          DEBUG_MSG(scan << "\t" << tags.size() << "\t" << tagNum << "\t" << tagSeq << "\t" <<  anno); 
          break;
        }
        tagNum++;
      }
      //break;
    }
  }

  int totalTags = 0;
  map<int, int>::iterator itrMap = mapTagDepthCount.begin();
  map<int, int>::iterator itrMapEnd = mapTagDepthCount.end();
  for (; itrMap != itrMapEnd; itrMap++) {
    DEBUG_MSG(itrMap->first << "\t" << itrMap->second); 
    totalTags += itrMap->second;
  }
  DEBUG_VAR(totalTags); 

} else {

  map<int, float>  mapContigHow;

 for (int x = 0; x < 3; x++) {

  int len;
  int jump;

  switch (x) {
    case 0:
      len = 6;
      jump = 2;
      maxNumTags = 100;
      break;
    case 1:
      len = 5;
      jump = 0;
      maxNumTags = 10;
      break;
    case 2:
      len = 5;
      jump = 1;
      maxNumTags = 5;
      break;
    case 3:
      len = 6;
      jump = 1;
      maxNumTags = 50;
      break;
    case 4:
      len = 4;
      jump = 1;
      maxNumTags = 1;
      break;
    case 5:
      len = 3;
      jump = 1;
      maxNumTags = 1;
      break;
    case 6:
      len = 6;
      jump = 0;
      maxNumTags = 50;
      break;
    case 7:
      len = 5;
      jump = 0;
      maxNumTags = 1;
      break;
  }

  int nNew = 0;

  for (int i = 0; i < contigs.consensus.size(); i++) {
    int scan = i+1;
    int correct = mapContigCorrect[scan];
    string anno = mapContigAnno[scan];
    for (int j = 0; j < anno.length(); j++) {
      if (anno[j] == 'L') anno[j] = 'I';
    }

    ExtractTags(contigs.consensus[i], tags, peakTol, len, jump, maxNumTags);
    if (correct == 1) {
      int tagNum = 1;
      list<Tag>::iterator itr = tags.begin();
      list<Tag>::iterator itr_end = tags.end();
      for (; itr != itr_end; itr++) {
        string tagSeq;
        for (int s = 0; s < itr->sequence.size(); s++) {
          int index = itr->sequence[s];
          tagSeq += jumps.aaLetters[(char)index];
        }
        size_t pos = anno.find(tagSeq);
        if (pos != string::npos) {

          map<int, float>::iterator itrFind = mapContigHow.find(scan);
          if (itrFind == mapContigHow.end()) {
            mapContigHow[scan] = (float)len + (float)jump / 10.0;
            nNew++;
          }
          break;
        }
        tagNum++;
      }
    }
  } // for (int i = 0; i < contigs.consensus.size(); i++) {

  DEBUG_MSG(len << "\t" << jump << "\t" << maxNumTags << "\t" << nNew); 

 } // for (int x = 0; x < 8; x++) {

} // if (argc == 7) {

  return 0;
}
