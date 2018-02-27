//
//  main_specdump - stand alone executable for determining where we lost IDs 
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

#define DEBUG_THIS 1

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;


// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 7 && argc != 2) {
    cerr << "Usage: main_scorecontig msgf_psm_file [tag_psm_file abinfo_file contig_psm_file spectrum_psm_file spectrum_fdr_file]" << endl;
    return -1;
  }
  string tagfile("assembly/tagsearchpsm.txt");
  string compfile("assembly/component_info.bin");
  string contigfile("homology/contigs_psm.txt");
  string specfile("homology/spectrum_psm.txt");
  string specfdrfile("homology/spectrum_psm_fdr.txt");

  if (argc != 2) {
    tagfile = argv[2];
    compfile = argv[3];
    contigfile = argv[4];
    specfile = argv[5];
    specfdrfile= argv[6];
  }
  DEBUG_VAR(tagfile);
  DEBUG_VAR(compfile);
  DEBUG_VAR(contigfile);
  DEBUG_VAR(specfile);
  DEBUG_VAR(specfdrfile);

  PeptideSpectrumMatchSet	psmSetMsgf;
  if (!psmSetMsgf.loadFromFile(argv[1])) {
    ERROR_MSG("Loading MSGF PSM file [" << argv[1] << "]");
    return -2;
  }

  PeptideSpectrumMatchSet	psmSetTag;
  if (!psmSetTag.loadFromFile(tagfile.c_str())) {
    ERROR_MSG("Loading Tag PSM file [" << argv[2] << "]");
    return -3;
  }

  PeptideSpectrumMatchSet	psmSetContig;
  if (!psmSetContig.loadFromFile(contigfile.c_str())) {
    ERROR_MSG("Loading Contig PSM file [" << argv[4] << "]");
    return -3;
  }

  PeptideSpectrumMatchSet	psmSetSpec;
  if (!psmSetSpec.loadFromFile(specfile.c_str())) {
    ERROR_MSG("Loading Spectrum PSM file [" << argv[5] << "]");
    return -3;
  }

  PeptideSpectrumMatchSet	psmSetSpecFdr;
  if (!psmSetSpecFdr.loadFromFile(specfdrfile.c_str())) {
    ERROR_MSG("Loading Spectrum FDR PSM file [" << argv[6] << "]");
    return -3;
  }

  map<int, set<string> > mapTagProteins;
  for (int iPSM = 0; iPSM < psmSetTag.size(); iPSM++ ) {
    psmPtr psmTag = psmSetTag[iPSM];
    mapTagProteins[psmTag->m_scanNum].insert(psmTag->m_protein);
  } // for (int iPSM = 0; iPSM < psmSetTag.size(); iPSM++ )

  abinfo_t contigAbinfo;
  if (!Load_abinfo(compfile.c_str(), contigAbinfo)) {
    ERROR_MSG("Loading Abinfo file [" << compfile.c_str() << "]");
    return -6;
  }
  //DEBUG_VAR(contigAbinfo.size());

  map<int, int> specToContig;  
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > >::iterator itr = contigAbinfo.begin();
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > >::iterator itr_end = contigAbinfo.end();

  for ( ; itr != itr_end; itr++) {
    int contigIndex = itr->first + 1;  // We prefer 1-based contigs
    //DEBUG_MSG("C = " << contigIndex);
    vector<int> specs = itr->second.first.first;
    for (int i = 0; i < specs.size(); i++) {
      int specIndex = specs[i] + 1;  // We prefer 1-based spectra
      //DEBUG_MSG("  S = " << specIndex);
      specToContig[specIndex] = contigIndex;
    }
  }

  map<int, psmPtr> mapScanMsgf;
  for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ ) {
    psmPtr psmSpec = psmSetMsgf[iPSM];
    mapScanMsgf[psmSpec->m_scanNum] = psmSpec;
  }

  map<int, psmPtr> mapScanContig;
  for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ ) {
    psmPtr psmSpec = psmSetMsgf[iPSM];
    mapScanMsgf[psmSpec->m_scanNum] = psmSpec;
  }

  map<int, psmPtr> mapScanSpec;
  for (int iPSM = 0; iPSM < psmSetSpec.size(); iPSM++ ) {
    psmPtr psmSpec = psmSetSpec[iPSM];
    mapScanSpec[psmSpec->m_scanNum] = psmSpec;
  }

  map<int, psmPtr> mapScanSpecFdr;
  for (int iPSM = 0; iPSM < psmSetSpecFdr.size(); iPSM++ ) {
    psmPtr psmSpec = psmSetSpecFdr[iPSM];
    mapScanSpecFdr[psmSpec->m_scanNum] = psmSpec;
  }

  map<int, bool> mapFoundTag;
  map<int, bool> mapFoundTagContig;
  map<int, bool> mapFoundContig;
  map<int, bool> mapFoundSpec;
  map<int, bool> mapFoundSpecFdr;

  //---------------------------------------------------------------------
  // Determine if we have a tag for each spectra and for each contig
  //---------------------------------------------------------------------
  map<int, psmPtr>::iterator itrMsgf = mapScanMsgf.begin();
  map<int, psmPtr>::iterator itrMsgfEnd = mapScanMsgf.end();
  for (; itrMsgf != itrMsgfEnd; itrMsgf++ ) {

    int specIndex = itrMsgf->second->m_scanNum;

    map<int, int>::iterator itrContig = specToContig.find(specIndex);
    if (itrContig == specToContig.end()) {
      continue;
    }
    mapFoundTag[specIndex] = false;

    int contigIndex = itrContig->second;
    //DEBUG_MSG(contigIndex << "  " << specIndex);

    // See if there were any tags for this contig
    map<int, set<string>  >::iterator itrTag = mapTagProteins.find(contigIndex);
    if (itrTag == mapTagProteins.end()) {
      //DEBUG_MSG("  " << "no tags for contig [" << contigIndex << "]");
      continue;
    }

    mapFoundTag[specIndex] = false;
    mapFoundTagContig[specIndex] = false;
    mapFoundContig[specIndex] = false;
    mapFoundSpec[specIndex] = false;
    mapFoundSpecFdr[specIndex] = false;

    set<string>::iterator itrProtein = itrTag->second.begin();
    set<string>::iterator itrProteinEnd = itrTag->second.end();
    for ( ; itrProtein != itrProteinEnd; itrProtein++) {
      if (*itrProtein == itrMsgf->second->m_protein) {
        //DEBUG_MSG("  " << "tag found contig [" << contigIndex << "] protein [" << *itrProtein << "]");
        //DEBUG_MSG(contigIndex << "  " << specIndex << "  " << *itrProtein);
        mapFoundTag[specIndex] = true;
        mapFoundTagContig[contigIndex] = true;
        break;
      }
    }

    if (mapScanMsgf.find(specIndex) != mapScanMsgf.end()) {
      mapFoundContig[specIndex] = true;
    }
    if (mapScanSpec.find(specIndex) != mapScanSpec.end()) {
      mapFoundSpec[specIndex] = true;
    }
    if (mapScanSpecFdr.find(specIndex) != mapScanSpecFdr.end()) {
      mapFoundSpecFdr[specIndex] = true;
    }

  } // for (int iPSM = 0; iPSM < mapScanMsgf.size(); iPSM++ )

#if 0
  DEBUG_VAR(mapScanMsgf.size());
  DEBUG_VAR(specToContig.size());
  DEBUG_VAR(mapFoundTag.size());
  DEBUG_VAR(mapFoundTagContig.size());
  DEBUG_VAR(mapFoundContig.size());
  DEBUG_VAR(mapFoundSpec.size());
  DEBUG_VAR(mapFoundSpecFdr.size());

  map<int, bool>::iterator itrContig = mapFoundContig.begin();
  map<int, bool>::iterator itrContigEnd = mapFoundContig.end();
  for (; itrContig  != itrContigEnd; itrContig++ ) {
    DEBUG_MSG(itrContig->first);
  }
#endif

  cout << "Contig#\tSpec#\tTag\tContig\tSpec\tSpecFDR\tScore\tpValue\tAnno\tMSGFAnno" << endl;
  map<int, bool>::iterator itrTag = mapFoundTag.begin();
  map<int, bool>::iterator itrTagEnd = mapFoundTag.end();
  for (; itrTag != itrTagEnd; itrTag++ ) {
    float score = 0.0;
    float pValue = 0.0;
    string anno;
    if (mapFoundSpec[itrTag->first]) {
      pValue = mapScanSpec[itrTag->first]->m_pValue;
      anno = mapScanSpec[itrTag->first]->m_annotation;
      score = mapScanSpec[itrTag->first]->m_score;
    }

    string annoMsgf;
    if (mapScanMsgf.find(itrTag->first) != mapScanMsgf.end()) {
      annoMsgf = mapScanMsgf[itrTag->first]->m_annotation;
    }

    cout << specToContig[itrTag->first] << "\t" << itrTag->first << "\t" << 
              mapFoundTag[itrTag->first] << "\t" << mapFoundContig[itrTag->first] << "\t" << 
              mapFoundSpec[itrTag->first] << "\t" << mapFoundSpecFdr[itrTag->first] << "\t" << 
              score << "\t" << pValue << "\t" << anno << "\t" << annoMsgf << endl;
  }

  return 0;
}
