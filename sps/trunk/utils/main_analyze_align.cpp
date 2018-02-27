//
//  main_analyze_align - stand alone executable for dumping spectrum data
//

#include "abruijn.h"
#include "denovo.h"
#include "ClusterData.h"
#include "CommandLineParser.h"
#include "ClusterSet.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ReportTableClusterConsensus.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include "tuple.h"

#include <stdlib.h>

#define USING_CANCER 0

#define DEBUG_ 0
#define DEBUG_PERCENT_MATCH 0
#define DEBUG_TAG_MATCH 0
#define DEBUG_CONTIG_INDEX -1
#define DEBUG_GAPS 0
#define DEBUG_GAP_NUMBER 0
#define DEBUG_CLUSTERS 0

#define MATCH_THRESHOLD 0.9
#define SIM_MATCH_THRESHOLD 0.8

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;

namespace sps
{
  bool operator<(const tuple<string, string, float> t1,
  const tuple<string, string, float> t2)
  {
    if (t1.m0 <= t2.m0)
      return true;
    if (t1.m1 <= t2.m1)
      return true;
    if (t1.m2 <= t2.m2)
      return true;
    return false;
  }
}

struct ContigResult
{
  int contig;
  int numSpectra;
  int contigPeaks;
  bool tagMatch;
  bool tagMatchGood;
  bool contigPsmGood;
  float alignScore;
  float bestScore;
  float actualScore;
  bool passFdr;
  bool oneSpecPassFdr;
  int contigMatch;
  int maxContigGap;
  float specProb;
  float fdr;
  bool orient;
  string tag;
  string annotation;
  string denovoAnno;
  string protein;

  ContigResult()
  {
    contig = -1;
    numSpectra = 0;
    contigPeaks = -1;
    tagMatch = false;
    tagMatchGood = false;
    contigPsmGood = false;
    alignScore = 0.0;
    bestScore = 0.0;
    actualScore = 0.0;
    passFdr = false;
    oneSpecPassFdr = false;
    contigMatch = -1;
    maxContigGap = -1;
    specProb = -1;
    fdr = -1;
    orient = 0;
    tag = "---";
    annotation = "---";
    denovoAnno = "---";
    protein = "---";
  }
};

struct SpectrumResult
{
  int contig;
  int spectrum;
  int specPeaks;
  bool passFdr;
  int specMatch;
#if USING_CANCER
  int numOrigSpectra;
  bool normal;
  bool cancer;
#endif
  string specAnno;
  string specAnnoOrig;
  string specAnnoClean;
  int specAnnoLength;
  float specTotalMass;
  float starParentMass;
  int numMods;
  int charge;
  int maxSpecGap;
  int suspiciousMods;
  int maxDbGap;
  string msgfAnno;
  string msgfAnnoClean;
  float msgfTotalMass;
  string specProtein;
  bool isDecoy;
  string msgfProtein;
  float alignScore;
  float specProb;
  float msgfProb;
  float fdr;
  float spChange;
  float matchPercent;
  bool orient;

  SpectrumResult()
  {
    contig = -1;
    spectrum = -1;
    specPeaks = -1;
    passFdr = false;
    specMatch = -1;
#if USING_CANCER
    numOrigSpectra = 0;
    normal = false;
    cancer = false;
#endif
    specAnno = "---";
    specAnnoOrig = "---";
    specAnnoClean = "---";
    specAnnoLength = -1;
    specTotalMass = 0.0;
    starParentMass = 0.0;
    numMods = -1;
    maxSpecGap = -1;
    suspiciousMods = 0;
    maxDbGap = -1;
    msgfAnno = "---";
    msgfAnnoClean = "---";
    msgfTotalMass = 0.0;
    specProtein = "---";
    isDecoy = false;
    msgfProtein = "---";
    specProb = 10.0;
    alignScore = -1;
    fdr = -1;
    msgfProb = -1;
    spChange = 0;
    matchPercent = 0;
    orient = 0;
  }
};

PeptideSpectrumMatchSet psmSetSpecFdr;
abinfo_t contigAbinfo;

map<int, SpectrumResult> specResults;
map<int, ContigResult> contigResults;
map<int, int> specToContig;
map<int, set<int> > contigToSpecSet;
map<int, int> mapContigType;
map<int, set<int> > mapContigVarSet;
map<int, set<int> > mapContigVarSetMsgf;
map<int, string> mapContigToMsgfProtein;
map<int, int> mapContigMsgfHybrid;
map<int, int> mapSpecToMatchType;
map<int, int> mapSpectrumToGoodMsgfID;
int stars_in_contigs_msgf_id = 0;
map<int, psmPtr> mapScanMsgf;
set<int> specPSMSet;
set<int> setNoMsgf;
set<int> setExactMatch;
set<int> setInexactMatch;
set<int> setAssociatedPeptide;
set<int> setAssociatedBadPeptide;
set<int> setDiffPeptide;

int nSimGood = 0;
int nSimBad = 0;
int nGoodContig = 0;
int nBadContig = 0;
int nHybridContig = 0;

int totalContigsWithGoodTags = 0;
int totalSpectraWithGoodTags = 0;
int totalSpectraInContigsWithGoodFdr = 0;
int totalContigsWithGoodFdrSpecs = 0;

int stars_id_msgf_noid_nocontigpsm = 0;
int stars_id_msgf_noid_belowfdr = 0;
int stars_id_msgf_noid_nogoodtag = 0;
int stars_id_msgf_noid_hybridcontig = 0;
int stars_id_msgf_noid_badcontig = 0;

int totalSpecsInContigs = 0;

// -------------------------------------------------------------------------
class PsmSetCustom : public PeptideSpectrumMatchSet
{
public:
  bool load(const char * fileName,
            map<int, set<tuple<string, string, float> > > & mapTuple
            ,
            map<int, ContigResult> & contigResults)
  {
    vector<int> mapIndex;
    ifstream ifsTagFile;
    readHeader(fileName, ifsTagFile, "\t", "", mapIndex);
    psmPtr nextPsm;
    bool more = true;
    while (more) {
      psmPtr psmSpec(new PeptideSpectrumMatch);
      more = readNextPsm(ifsTagFile, mapIndex, "\t", "", false, false, psmSpec);
      if (more) {
        mapTuple[psmSpec->m_scanNum].insert(make_tuple<string, string, float>(psmSpec->m_annotation,
                                                                              psmSpec->m_protein,
                                                                              psmSpec->m_score));
        if (contigResults[psmSpec->m_scanNum].bestScore < psmSpec->m_score) {
          contigResults[psmSpec->m_scanNum].bestScore = psmSpec->m_score;
        }
      }
    }
    ifsTagFile.close();

    return true;
  }

};

// -------------------------------------------------------------------------
void addFlankingAas(psmPtr psm, DB_fasta & dbFasta)
{
  string cleanAnno;
  PeptideSpectrumMatch::getUnmodifiedPeptide(psm->m_annotation, cleanAnno);
  string dbSequence = dbFasta.getSequence(psm->m_dbIndex);
  string prevAA = "_";
  string postAA = "_";

  int start = dbSequence.find(cleanAnno);
  if (start != string::npos) {
    if (start != 0) {
      prevAA = dbSequence[start - 1];
    }
    int end = start + cleanAnno.length() - 1;
    if (end < dbSequence.length() - 1) {
      postAA = dbSequence[end + 1];
    }
  }

  string newAnnotation = prevAA + "." + psm->m_annotation + "." + postAA;
  //DEBUG_VAR(newAnnotation);
  psm->m_annotation = newAnnotation;
  //DEBUG_VAR(psm->m_annotation);
}

// -------------------------------------------------------------------------
string convertAnnotation(string annotation)
{
  string convertedAnnotation;
  for (int iChar = 0; iChar < annotation.length(); iChar++) {
    if (annotation[iChar] == 'I') {
      convertedAnnotation += 'L';
    }
    else if (annotation[iChar] == 'Q') {
      convertedAnnotation += 'K';
    }
    else {
      convertedAnnotation += annotation[iChar];
    }
  }
  return convertedAnnotation;
}

// -------------------------------------------------------------------------
float getPercentMatch(string & string1, string & string2)
{
  if (DEBUG_PERCENT_MATCH)
    DEBUG_VAR(string1.size());
  if (DEBUG_PERCENT_MATCH)
    DEBUG_VAR(string2.size());
  int dynaArray[string1.size() + 1][string2.size() + 1];
  for (int i = 0; i < string1.size() + 1; i++) {
    for (int j = 0; j < string2.size() + 1; j++) {
      dynaArray[i][j] = 0;
    }
  }
  for (int i = 1; i <= string1.size(); i++) {
    for (int j = 1; j <= string2.size(); j++) {
      int up = dynaArray[i - 1][j];
      int left = dynaArray[i][j - 1];
      int diag = dynaArray[i - 1][j - 1];
      if (DEBUG_PERCENT_MATCH)
        DEBUG_MSG(i << "  " << j);
      if (DEBUG_PERCENT_MATCH)
        DEBUG_MSG(up << "  " << left << "  " << diag);
      if (DEBUG_PERCENT_MATCH)
        DEBUG_MSG(string1[i - 1] << "  " << string2[j - 1]);
      if (string1[i - 1] == string2[j - 1]) {
        diag++;
      }
      if (DEBUG_PERCENT_MATCH)
        DEBUG_MSG(up << "  " << left << "  " << diag);
      int max1 = max(left, up);
      dynaArray[i][j] = max(max1, diag);
      if (DEBUG_PERCENT_MATCH)
        DEBUG_MSG(dynaArray[i][j]);
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
  if (DEBUG_PERCENT_MATCH)
    DEBUG_VAR(match);
  float percent = (float)match / (float)min(string1.size(), string2.size());
  if (DEBUG_PERCENT_MATCH)
    DEBUG_VAR(percent);
  percent = min(percent, (float)1.0);
  if (DEBUG_PERCENT_MATCH)
    DEBUG_VAR(percent);
  return percent;
}

//---------------------------------------------------------------------
// Output the spectrum results table
//---------------------------------------------------------------------
void outputSpectrumLevelResults(ostream & ofs)
{
#if USING_CANCER
  ofs << "contig#\tscan#\ttag_exist\ttag_good\tcontig_good\tcontig_1fdr\tpass_fdr\tdecoy\tcontig_match\tmsgf_hybrid\tspec_match\ttag\tnum_orig_spectra\tnormal\tcancer\tspec_anno_orig\tspec_anno\tanno_clean\tlength\tparent_mass\tspec_mass\tmsgf_mass\tmods\tcharge\tmsgf_anno\tmsgf_anno_clean\tspec_prot\tmsgf_prot\tspec_prob\talign_score\tfdr\tmsgf_prob\tSP_change\tmatch%\torient" << endl;
#else
  ofs
      << "contig#\tscan#\ttag_exist\ttag_good\tcontig_good\tcontig_1fdr\tpass_fdr\tdecoy\tcontig_match\tmsgf_hybrid\tspec_match\ttag\tspec_anno_orig\tspec_anno\tanno_clean\tlength\tparent_mass\tspec_mass\tmsgf_mass\tmods\tcharge\tmsgf_anno\tmsgf_anno_clean\tspec_prot\tmsgf_prot\tspec_prob\talign_score\tfdr\tmsgf_prob\tSP_change\tmatch%\torient"
      << endl;
#endif

  map<int, int>::iterator itrs = specToContig.begin();
  map<int, int>::iterator itrsEnd = specToContig.end();
  for (; itrs != itrsEnd; itrs++) {

    int specIndex = itrs->first;
    int contigIndex = specResults[specIndex].contig;

    ofs << contigIndex << "\t";
    ofs << specIndex << "\t";

    //ofs << contigResults[contigIndex].contigPeaks << "\t";
    //ofs << specResults[specIndex].specPeaks << "\t";

    ofs << contigResults[contigIndex].tagMatch << "\t";
    ofs << contigResults[contigIndex].tagMatchGood << "\t";
    ofs << contigResults[contigIndex].contigPsmGood << "\t";
    ofs << contigResults[contigIndex].passFdr << "\t";
    ofs << specResults[specIndex].passFdr << "\t";
    ofs << specResults[specIndex].isDecoy << "\t";

    switch (mapContigType[contigIndex]) {
    case 1:
      ofs << "good\t";
      break;
    case 3:
      ofs << "bad\t";
      break;
    case 8:
      ofs << "hybrid\t";
      break;
    default:
      ofs << "unk\t";
      break;
    }

    switch (mapContigMsgfHybrid[contigIndex]) {
    case 0:
      ofs << "homog\t";
      break;
    case 1:
      ofs << "hybrid\t";
      break;
    case -1:
    default:
      ofs << "unk\t";
      break;
    }

    //DEBUG_MSG(specIndex << "  " << mapSpecToMatchType[specIndex]);
    switch (mapSpecToMatchType[specIndex]) {
    case 1:
    case 2:
      ofs << "match\t";
      break;
    case 3:
      ofs << "nonmatch\t";
      break;
    case 4:
      ofs << "simgood\t";
      break;
    case 5:
      ofs << "simbad\t";
      break;
    case 6:
      ofs << "msgfbad\t";
      break;
    default:
      ofs << "unk\t";
      break;
    }

    ofs << contigResults[contigIndex].tag << "\t";
#if USING_CANCER
    ofs << specResults[specIndex].numOrigSpectra << "\t";
    ofs << specResults[specIndex].normal << "\t";
    ofs << specResults[specIndex].cancer << "\t";
#endif
    ofs << specResults[specIndex].specAnnoOrig << "\t";
    ofs << specResults[specIndex].specAnno << "\t";
    ofs << specResults[specIndex].specAnnoClean << "\t";
    ofs << specResults[specIndex].specAnnoLength << "\t";

    ofs << specResults[specIndex].starParentMass << "\t";
    ofs << specResults[specIndex].specTotalMass << "\t";
    ofs << specResults[specIndex].msgfTotalMass << "\t";
    ofs << specResults[specIndex].numMods << "\t";
    ofs << specResults[specIndex].charge << "\t";

#if 0
    ofs << contigResults[contigIndex].maxContigGap << "\t";
    ofs << specResults[specIndex].maxSpecGap << "\t";
    ofs << specResults[specIndex].suspiciousMods << "\t";
#endif

    ofs << specResults[specIndex].msgfAnno << "\t";
    ofs << specResults[specIndex].msgfAnnoClean << "\t";
    ofs << specResults[specIndex].specProtein << "\t";
    ofs << specResults[specIndex].msgfProtein << "\t";
    ofs << specResults[specIndex].specProb << "\t";
    ofs << specResults[specIndex].alignScore << "\t";
    ofs << specResults[specIndex].fdr << "\t";
    ofs << specResults[specIndex].msgfProb << "\t";
    ofs << specResults[specIndex].spChange << "\t";
    ofs << specResults[specIndex].matchPercent << "\t";
    ofs << specResults[specIndex].orient << endl;

    if (contigResults[contigIndex].tagMatchGood) {
      totalSpectraWithGoodTags++;
    }

    if (contigResults[contigIndex].oneSpecPassFdr) {
      totalSpectraInContigsWithGoodFdr++;
    }

    if (!contigResults[contigIndex].contigPsmGood
        && !specResults[specIndex].passFdr
        && (specResults[specIndex].msgfProb > 0)) {
      stars_id_msgf_noid_nocontigpsm++;
    }

    if (!contigResults[contigIndex].contigPsmGood
        && !specResults[specIndex].passFdr
        && (specResults[specIndex].msgfProb > 0)
        && !contigResults[contigIndex].tagMatchGood) {
      stars_id_msgf_noid_nogoodtag++;
    }

    if (!contigResults[contigIndex].contigPsmGood
        && !specResults[specIndex].passFdr
        && (specResults[specIndex].msgfProb > 0)
        && contigResults[contigIndex].tagMatchGood
        && (mapContigType[contigIndex] == 8)) {
      stars_id_msgf_noid_hybridcontig++;
    }

    if (!contigResults[contigIndex].contigPsmGood
        && !specResults[specIndex].passFdr
        && (specResults[specIndex].msgfProb > 0)
        && contigResults[contigIndex].tagMatchGood
        && (mapContigType[contigIndex] == 3)) {
      stars_id_msgf_noid_badcontig++;
    }

    if (contigResults[contigIndex].contigPsmGood
        && !specResults[specIndex].passFdr
        && (specResults[specIndex].msgfProb > 0)) {
      stars_id_msgf_noid_belowfdr++;
    }

    if (specResults[specIndex].passFdr) {
      mapContigVarSet[contigIndex].insert((int)specResults[specIndex].specTotalMass);
    }
    if (!specResults[specIndex].msgfAnno.empty()) {
      mapContigVarSetMsgf[contigIndex].insert((int)specResults[specIndex].msgfTotalMass);
    }
  }

  return;
}

//---------------------------------------------------------------------
// Output the contig results table
//---------------------------------------------------------------------
void outputContigLevelResults(ostream & ofs)
{
  ofs
      << "contig#\tnum_spectra\tnum_var\tnum_var_msgf\tcontig_peaks\ttag_exist\ttag_good\tcontig_good\tpass_fdr\tone_spec_fdr\ttag\tdenovo\tannotation\tmsgf_contig\tcontig_type\tactual_score\tbest_score\tfdr\tmax_contig_gap\torient"
      << endl;

  map<int, set<int> >::iterator itrc = contigToSpecSet.begin();
  map<int, set<int> >::iterator itrcEnd = contigToSpecSet.end();
  for (; itrc != itrcEnd; itrc++) {

    int contigIndex = itrc->first;

    ofs << contigIndex << "\t";
    ofs << contigResults[contigIndex].numSpectra << "\t";

    ofs << mapContigVarSet[contigIndex].size() << "\t";
    ofs << mapContigVarSetMsgf[contigIndex].size() << "\t";

    ofs << contigResults[contigIndex].contigPeaks << "\t";
    ofs << contigResults[contigIndex].tagMatch << "\t";
    ofs << contigResults[contigIndex].tagMatchGood << "\t";
    ofs << contigResults[contigIndex].contigPsmGood << "\t";
    ofs << contigResults[contigIndex].passFdr << "\t";
    ofs << contigResults[contigIndex].oneSpecPassFdr << "\t";
    ofs << contigResults[contigIndex].tag << "\t";
    ofs << contigResults[contigIndex].denovoAnno << "\t";
    ofs << contigResults[contigIndex].annotation << "\t";

    switch (mapContigMsgfHybrid[contigIndex]) {
    case 0:
      ofs << "homog\t";
      break;
    case 1:
      ofs << "hybrid\t";
      break;
    case -1:
    default:
      ofs << "unk\t";
      break;
    }

    switch (mapContigType[contigIndex]) {
    case 1:
      ofs << "good\t";
      break;
    case 3:
      ofs << "bad\t";
      break;
    case 8:
      ofs << "hybrid\t";
      break;
    default:
      ofs << "unk\t";
      break;
    }

    ofs << contigResults[contigIndex].actualScore << "\t";
    ofs << contigResults[contigIndex].bestScore << "\t";
    ofs << contigResults[contigIndex].fdr << "\t";
    ofs << contigResults[contigIndex].maxContigGap << "\t";
    ofs << contigResults[contigIndex].orient << endl;
  }

  return;
}

//---------------------------------------------------------------------
// Output the summary results table
//---------------------------------------------------------------------
void outputSummaryResults(ofstream & ofs)
{
  ofs << "<html>" << endl;
  ofs << "  <head>" << endl;
  ofs << "    <title>Summary Report</title>" << endl;
  ofs << "  </head>" << endl;

  ofs << "  <body>" << endl;
  ofs << "    <h1>Summary Report</h1>" << endl;

  ofs << "    <h2>Spectrum Level Summary</h2>" << endl;

  ofs << "    <table border=\"1\">" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in All contigs</td>" << endl;
  ofs << "        <td>" << totalSpecsInContigs << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with Good Tag</td>" << endl;
  ofs << "        <td>" << totalSpectraWithGoodTags << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with 1 PSM That Passed FDR</td>" << endl;
  ofs << "        <td>" << totalSpectraInContigsWithGoodFdr << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs Passed FDR</td>" << endl;
  ofs << "        <td>" << psmSetSpecFdr.size() << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "    </table>" << endl;

  ofs << "<br><br>" << endl;

  ofs << "    <table border=\"1\">" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs Identified by MSGFDB</td>" << endl;
  ofs << "        <td>" << stars_in_contigs_msgf_id << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with Matching MSGF and BSAAP Ids</td>"
      << endl;
  ofs << "        <td>" << setExactMatch.size() << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs
      << "        <td>Stars in Contigs with Non-matching MSGF and BSAAP Ids</td>"
      << endl;
  ofs << "        <td>" << setDiffPeptide.size() << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with MSGFDB Id but no BSAAP Id</td>"
      << endl;
  ofs << "        <td>"
      << stars_in_contigs_msgf_id - setExactMatch.size() - setDiffPeptide.size()
      << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "    </table>" << endl;

  ofs << "<br><br>" << endl;

  ofs << "    <table border=\"1\">" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with MSGFDB Id but no PSM in Contig</td>"
      << endl;
  ofs << "        <td>" << stars_id_msgf_noid_nocontigpsm << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs
      << "        <td>Stars in Contigs with MSGFDB Id but no Good Contig Tag</td>"
      << endl;
  ofs << "        <td>" << stars_id_msgf_noid_nogoodtag << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with MSGFDB Id but Hybrid Contig</td>"
      << endl;
  ofs << "        <td>" << stars_id_msgf_noid_hybridcontig << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with MSGFDB Id but Bad Contig</td>"
      << endl;
  ofs << "        <td>" << stars_id_msgf_noid_badcontig << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs
      << "        <td>Stars in Contigs with MSGFDB Id but BSAAP Id is Below FDR</td>"
      << endl;
  ofs << "        <td>" << stars_id_msgf_noid_belowfdr << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with no MSGFDB Id</td>" << endl;
  ofs << "        <td>" << setNoMsgf.size() << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with Id Similar to Good Id</td>" << endl;
  ofs << "        <td>" << nSimGood << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Stars in Contigs with Id Similar to Bad Id</td>" << endl;
  ofs << "        <td>" << nSimBad << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "    </table>" << endl;

  ofs << "    <h2>Contig Level Summary</h2>" << endl;

  ofs << "    <table border=\"1\">" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Contigs Total</td>" << endl;
  ofs << "        <td>" << contigAbinfo.size() << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Contigs with Correct Tag</td>" << endl;
  ofs << "        <td>" << totalContigsWithGoodTags << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Contigs with at Least 1 Id that Passed FDR</td>" << endl;
  ofs << "        <td>" << totalContigsWithGoodFdrSpecs << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Contigs with Correct Tag</td>" << endl;
  ofs << "        <td>" << nGoodContig << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Contigs with Incorrect Tag</td>" << endl;
  ofs << "        <td>" << nBadContig << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Contigs Hybrid</td>" << endl;
  ofs << "        <td>" << nHybridContig << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "      <tr>" << endl;
  ofs << "        <td>Contigs Unknown</td>" << endl;
  ofs << "        <td>"
      << contigAbinfo.size() - nGoodContig - nBadContig - nHybridContig
      << "</td>" << endl;
  ofs << "      </tr>" << endl;

  ofs << "    </table>" << endl;

  ofs << "  </body>" << endl;
  ofs << "</html>" << endl;

  return;
}

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));
  LoggerCleaner loggerCleaner; // Clears loggers on end of scope

  if (argc < 3) {
    cerr << "Usage: main_analyze_align <msgf_psm_file> <database_file> <out_file>"
        << endl;
    return -1;
  }

  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("notags", "NO_TAGS", 0));
  listOptions.push_back(CommandLineParser::Option("tagfile", "TAG_FILE", 1));
  listOptions.push_back(CommandLineParser::Option("peaktol",
                                                  "PEAK_TOLERANCE",
                                                  1));

  listOptions.push_back(CommandLineParser::Option("aafile",
                                                  "AMINO_ACID_FILE",
                                                  1));
  listOptions.push_back(CommandLineParser::Option("cmpfile",
                                                  "COMPONENT_FILE",
                                                  1));
  listOptions.push_back(CommandLineParser::Option("cspecfile",
                                                  "CONTIG_SPECS_FILE",
                                                  1));

  listOptions.push_back(CommandLineParser::Option("callpsm",
                                                  "CONTIG_ALL_PSM_FILE",
                                                  1));
  listOptions.push_back(CommandLineParser::Option("callpsmd",
                                                  "CONTIG_ALL_PSM_DEBUG_FILE",
                                                  1));
  listOptions.push_back(CommandLineParser::Option("cpeaks",
                                                  "CONTIG_PEAKS_FILE",
                                                  1));

  listOptions.push_back(CommandLineParser::Option("sfulltgt",
                                                  "SPECTRUM_FULL_TGT_FILE",
                                                  1));
  listOptions.push_back(CommandLineParser::Option("sfulldec",
                                                  "SPECTRUM_FULL_DEC_FILE",
                                                  1));
  listOptions.push_back(CommandLineParser::Option("starpeak",
                                                  "STAR_PEAKS_FILE",
                                                  1));

  listOptions.push_back(CommandLineParser::Option("starspec",
                                                  "STAR_SPECS_FILE",
                                                  1));

  listOptions.push_back(CommandLineParser::Option("contigfdrpsm",
                                                  "CONTIG_FDR_PSM_FILE",
                                                  1));
  listOptions.push_back(CommandLineParser::Option("spectrafdrpsm",
                                                  "SPECTRUM_FDR_PSM_FILE",
                                                  1));

  listOptions.push_back(CommandLineParser::Option("clustmap",
                                                  "CLUSTER_MAPPING_FILE",
                                                  1));

  listOptions.push_back(CommandLineParser::Option("outputdir",
                                                  "RESULTS_OUTPUT_DIR",
                                                  1));

  CommandLineParser clp(argc, argv, 3, listOptions);
  string parserError = "";
  if (!clp.validate(parserError)) {
    cerr << "Usage: main_analyze_align msgf_psm_file database_file out_file"
        << endl;
    cerr << "  " << parserError << endl;
    return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  string tagfile("assembly/tagsearchpsm.txt");
  if (commandLineParams.exists("TAG_FILE")) {
    tagfile = commandLineParams.getValue("TAG_FILE");
  }
  if (commandLineParams.exists("NO_TAGS")) {
    tagfile = "";
  }

  string compfile("assembly/component_info.bin");
  if (commandLineParams.exists("COMPONENT_FILE")) {
    compfile = commandLineParams.getValue("COMPONENT_FILE");
  }
  string contigspecfile("assembly/sps_seqs.pklbin");
  if (commandLineParams.exists("CONTIG_SPECS_FILE")) {
    contigspecfile = commandLineParams.getValue("CONTIG_SPECS_FILE");
  }

  string contigsallpsmfile("homology/contigs_psm_tgt.txt");
  if (commandLineParams.exists("CONTIG_ALL_PSM_FILE")) {
    contigsallpsmfile = commandLineParams.getValue("CONTIG_ALL_PSM_FILE");
  }
  string contigsallpsmfiledebug("homology/contigs_psm_tgt.txt.debug");
  if (commandLineParams.exists("CONTIG_ALL_PSM_DEBUG_FILE")) {
    contigsallpsmfiledebug =
        commandLineParams.getValue("CONTIG_ALL_PSM_DEBUG_FILE");
  }
  string contigpeaksfile("homology/contigs_midx_tgt.pklbin");
  if (commandLineParams.exists("CONTIG_PEAKS_FILE")) {
    contigpeaksfile = commandLineParams.getValue("CONTIG_PEAKS_FILE");
  }

  string specfulltgtfile("homology/spectrum_psm_tgt.txt");
  if (commandLineParams.exists("SPECTRUM_FULL_TGT_FILE")) {
    specfulltgtfile = commandLineParams.getValue("SPECTRUM_FULL_TGT_FILE");
  }
  string specfulldecfile("homology/spectrum_psm_dec.txt");
  if (commandLineParams.exists("SPECTRUM_FULL_DEC_FILE")) {
    specfulldecfile = commandLineParams.getValue("SPECTRUM_FULL_DEC_FILE");
  }
  string starpeaksfile("homology/spectrum_midx_tgt.pklbin");
  if (commandLineParams.exists("STAR_PEAKS_FILE")) {
    starpeaksfile = commandLineParams.getValue("STAR_PEAKS_FILE");
  }

  string aafile("homology/specprotalign_aa_masses.txt");
  if (commandLineParams.exists("AMINO_ACID_FILE")) {
    aafile = commandLineParams.getValue("AMINO_ACID_FILE");
  }

  string starspecfile("spectra/stars.pklbin");
  if (commandLineParams.exists("STAR_SPECS_FILE")) {
    starspecfile = commandLineParams.getValue("STAR_SPECS_FILE");
  }

  string contigsfdrcleanfile("");
  if (commandLineParams.exists("CONTIG_FDR_PSM_FILE")) {
    contigsfdrcleanfile = commandLineParams.getValue("CONTIG_FDR_PSM_FILE");
  }
  string specfdrfile("homology/spectrum_psm_fdr.txt");
  if (commandLineParams.exists("SPECTRUM_FDR_PSM_FILE")) {
    specfdrfile = commandLineParams.getValue("SPECTRUM_FDR_PSM_FILE");
  }

  string clustermappingfile("spectra/specs_ms.clust");
  if (commandLineParams.exists("CLUSTER_MAPPING_FILE")) {
    clustermappingfile = commandLineParams.getValue("CLUSTER_MAPPING_FILE");
  }

  DEBUG_MSG("Loading MSGF PSM file ...");
  PeptideSpectrumMatchSet psmSetMsgf;
  if (!psmSetMsgf.loadFromFile(argv[1])) {
    ERROR_MSG("Loading MSGF PSM file [" << argv[1] << "]");
    return -2;
  }
  //psmSetMsgf.saveToFile("psm_msgf.psm");

  DEBUG_MSG("Loading Database ...");
  DB_fasta dbFasta;
  if (dbFasta.Load(argv[2]) <= 0) {
    ERROR_MSG("Loading FASTA file [" << argv[2] << "]");
    return -3;
  }

  PeptideSpectrumMatchSet psmSetContigTag;
  if (!tagfile.empty()) {
    DEBUG_MSG("Loading tag file [" << tagfile << "]");
    if (!psmSetContigTag.loadFromFile(tagfile.c_str())) {
      ERROR_MSG("Loading Tag PSM file [" << tagfile << "]");
      return -4;
    }
  }
  //psmSetContigTag.saveToFile("psm_contig_tag.psm");

  DEBUG_MSG("Mapping Contigs to Protein Sets ...");
  map<int, set<string> > mapContigToProteinSet;
  for (int iPSM = 0; iPSM < psmSetContigTag.size(); iPSM++) {
    psmPtr psmTag = psmSetContigTag[iPSM];
    mapContigToProteinSet[psmTag->m_scanNum].insert(psmTag->m_protein);
  } // for (int iPSM = 0; iPSM < psmSetContigTag.size(); iPSM++ )

  PeptideSpectrumMatchSet psmSetContigFdr;
  if (!contigsfdrcleanfile.empty()) {
    DEBUG_MSG("Loading Contig PSM file [" << contigsfdrcleanfile << "]");
    if (!psmSetContigFdr.loadFromFile(contigsfdrcleanfile.c_str())) {
      ERROR_MSG("Loading Contig PSM file [" << contigsfdrcleanfile << "]");
      return -5;
    }
    //psmSetContigFdr.saveToFile("psm_contig.psm");
  }

  DEBUG_MSG("Loading Spectra PSM file [" << specfdrfile << "]");
  if (!psmSetSpecFdr.loadFromFile(specfdrfile.c_str())) {
    ERROR_MSG("Loading Spectra PSM file [" << specfdrfile << "]");
    return -6;
  }
  //psmSetSpecFdr.saveToFile("psm_spec.psm");

  DEBUG_MSG("Loading Spectra Full PSM file [" << specfulltgtfile << "]");
  PeptideSpectrumMatchSet psmSetSpecFullTgt;
  if (!psmSetSpecFullTgt.loadFromFile(specfulltgtfile.c_str())) {
    ERROR_MSG("Loading Spectra Full PSM file [" << specfulltgtfile << "]");
    return -7;
  }
  //psmSetSpecFullTgt.saveToFile("psm_full_spec.psm");

  DEBUG_MSG("Loading Spectra Full PSM file [" << specfulldecfile << "]");
  PeptideSpectrumMatchSet psmSetSpecFullDec;
  if (!psmSetSpecFullDec.loadFromFile(specfulldecfile.c_str())) {
    ERROR_MSG("Loading Spectra Full PSM file [" << specfulldecfile << "]");
    return -8;
  }
  //psmSetSpecFullDec.saveToFile("psm_full_spec.psm");

  DEBUG_MSG("Loading Abinfo file [" << compfile.c_str() << "]");
  if (!Load_abinfo(compfile.c_str(), contigAbinfo)) {
    ERROR_MSG("Loading Abinfo file [" << compfile.c_str() << "]");
    return -9;
  }
  //DEBUG_VAR(contigAbinfo.size());

  DEBUG_MSG("Loading Contig Spectra file [" << contigspecfile << "]");
  SpecSet contigSpectra;
  if (contigSpectra.loadPklBin(contigspecfile.c_str()) <= 0) {
    ERROR_MSG("Loading Contig Spectra file [" << contigspecfile << "]");
    return -10;
  }
  //DEBUG_VAR(contigSpectra.size());

  DEBUG_MSG("Loading Star Spectra file [" << starspecfile << "]");
  SpecSet starSpectra;
  if (starSpectra.loadPklBin(starspecfile.c_str()) <= 0) {
    ERROR_MSG("Loading Star Spectra file [" << starspecfile << "]");
    return -11;
  }
  //DEBUG_VAR(starSpectra.size());

  DEBUG_MSG("Loading contig peaks file [" << contigpeaksfile << "]");
  SpecSet contigPeakMatches;
  if (contigPeakMatches.loadPklBin(contigpeaksfile.c_str()) <= 0) {
    ERROR_MSG("Loading contig peaks file [" << contigpeaksfile << "]");
    return -12;
  }
  //DEBUG_VAR(contigPeaks.size());

  DEBUG_MSG("Loading star peaks file [" << starpeaksfile << "]");
  SpecSet starPeaks;
  if (starPeaks.loadPklBin(starpeaksfile.c_str()) <= 0) {
    ERROR_MSG("Loading star peaks file [" << starpeaksfile << "]");
    return -13;
  }
  //DEBUG_VAR(starPeaks.size());

  DEBUG_MSG("Loading Contigs PSM file [" << contigsallpsmfile << "]");
  PeptideSpectrumMatchSet psmSetContigsAll;
  if (!psmSetContigsAll.loadFromFile(contigsallpsmfile.c_str())) {
    ERROR_MSG("Loading Contigs PSM file [" << contigsallpsmfile << "]");
    return -14;
  }
  //psmSetContigsAll.saveToFile("psm_msgf.psm");

  DEBUG_MSG("Iterating through all contigs/spectra...");

  //----------------------------------------------------------
  // Iterate through the abinfo and get all contigs and spectra
  //----------------------------------------------------------
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int>, vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int>, vector<double> > // Spectrum index, peak mass
          > > >::iterator itr = contigAbinfo.begin();
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int>, vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int>, vector<double> > // Spectrum index, peak mass
          > > >::iterator itr_end = contigAbinfo.end();
  for (; itr != itr_end; itr++) {
    int contigIndex = itr->first + 1; // We prefer 1-based contigs

    ContigResult newContigResult;
    contigResults[contigIndex] = newContigResult;

    //DEBUG_MSG("C = " << contigIndex);
    vector<int> specs = itr->second.first.first;
//    totalSpecsInContigs += specs.size();
    for (int i = 0; i < specs.size(); i++) {

      if (specs[i] >= starSpectra.size()) {
        static bool warningFlag = true;
        if (warningFlag) {
          WARN_MSG(specs[i] << " out of bounds!");
          WARN_MSG("  Further warnings of this type have been supressed");
          warningFlag = false;
        }
        continue;
      }

      int specIndex;
#if 0
      specIndex = specs[i] + 1; // We prefer 1-based spectra
#else
      int verIndex = specs[i];
      //DEBUG_MSG("  V = " << verIndex);
      specIndex = starSpectra[verIndex].scan;
#endif
      //DEBUG_MSG("  S = " << specIndex);
      if (specToContig.find(specIndex) != specToContig.end()) {
        //WARN_MSG(specIndex << " was found twice!"<< "  " << contigIndex << " and " << specToContig[specIndex ]);
      }
      else {
        specToContig[specIndex] = contigIndex;
        contigToSpecSet[contigIndex].insert(specIndex);
        totalSpecsInContigs++;
        contigResults[contigIndex].numSpectra++;
      }
    }
  }
  DEBUG_VAR(specToContig.size());
  DEBUG_VAR(contigToSpecSet.size());

  DEBUG_MSG("Loading Contigs PSM file DEBUG ...");
  map<int, set<tuple<string, string, float> > > mapContigToTupleAnnoProtScore2;
  PsmSetCustom psmSetContigsAllDebug;
  if (!psmSetContigsAllDebug.load(contigsallpsmfiledebug.c_str(),
                                  mapContigToTupleAnnoProtScore2,
                                  contigResults)) {
    ERROR_MSG("Loading Contigs PSM file DEBUG [" << contigsallpsmfiledebug
        << "]");
    return -14;
  }
  //psmSetContigsAll.saveToFile("psm_msgf.psm");

  DEBUG_MSG("Setting up contig results map...");

  //----------------------------------------------------------
  // Setup the contig results vector
  //----------------------------------------------------------
  map<int, set<int> >::iterator itrc = contigToSpecSet.begin();
  map<int, set<int> >::iterator itrcEnd = contigToSpecSet.end();
  for (; itrc != itrcEnd; itrc++) {
    if (DEBUG_GAPS)
      DEBUG_VAR(itrc->first);
    contigResults[itrc->first].contig = itrc->first;
    contigResults[itrc->first].tagMatch =
        (mapContigToProteinSet.find(itrc->first) != mapContigToProteinSet.end());
    int contigIndex0 = itrc->first - 1;
    contigResults[itrc->first].contigPeaks = contigSpectra[contigIndex0].size();
#if 0
    if (DEBUG_GAPS) DEBUG_VAR(contigResults[itrc->first].contigPeaks);
    contigResults[itrc->first].maxContigGap = -2;
    if (DEBUG_GAPS) DEBUG_VAR(contigPeakMatches[contigIndex0].size());
    if (contigPeakMatches[contigIndex0].size() > 0) {
      for (int i = 0; i < contigPeakMatches[contigIndex0].size() - 1; i++) {
        if (DEBUG_GAPS && DEBUG_GAP_NUMBER == i) DEBUG_VAR(i+1);
        if (DEBUG_GAPS && DEBUG_GAP_NUMBER == i) DEBUG_VAR(contigPeakMatches[contigIndex0][i+1][0]);
        if (DEBUG_GAPS && DEBUG_GAP_NUMBER == i) DEBUG_VAR(contigPeakMatches[contigIndex0][i][0]);
        float peakGapContig = contigPeakMatches[contigIndex0][i+1][0] - contigPeakMatches[contigIndex0][i][0];
        if (DEBUG_GAPS && DEBUG_GAP_NUMBER == i) DEBUG_VAR(peakGapContig);
        if (DEBUG_GAPS && DEBUG_GAP_NUMBER == i) DEBUG_VAR(contigResults[itrc->first].maxContigGap);
        if (peakGapContig > contigResults[itrc->first].maxContigGap) {
          contigResults[itrc->first].maxContigGap = peakGapContig;
        }
      }
    }
#endif
  }

  DEBUG_MSG("Setting up spectrum results map...");

  //----------------------------------------------------------
  // Setup the spectrum results vector
  //----------------------------------------------------------
  map<int, int>::iterator itrx = specToContig.begin();
  map<int, int>::iterator itrxEnd = specToContig.end();
  for (; itrx != itrxEnd; itrx++) {
    //DEBUG_VAR(itrx->first);
    SpectrumResult newSpecResult;
    specResults[itrx->first] = newSpecResult;
    specResults[itrx->first].contig = itrx->second;
    specResults[itrx->first].spectrum = itrx->first;

    int specIndex0 = itrx->first - 1;
    //DEBUG_VAR(specIndex0);
    if (specIndex0 < 0 || specIndex0 >= starSpectra.size()) {
      continue;
    }

    specResults[itrx->first].specPeaks = starSpectra[specIndex0].size();
    specResults[itrx->first].starParentMass =
        starSpectra[specIndex0].parentMass;
    specResults[itrx->first].charge = starSpectra[specIndex0].parentCharge;
  }

#if USING_CANCER
  //----------------------------------------------------------
  // Find the original files and scans for the spectra
  //----------------------------------------------------------
  DEBUG_MSG("Getting cluster information...");
  ClusterSet clusterSet;
  clusterSet.loadBinaryFile(clustermappingfile);
  map<int, list<pair<int, string> > > clusterInfo;
  clusterSet.getScanToFileMapping(clusterInfo);
  DEBUG_VAR(clusterInfo.size());

  map<int, list<pair<int, string> > >::iterator itrClust = clusterInfo.begin();
  map<int, list<pair<int, string> > >::iterator itrClustEnd = clusterInfo.end();
  for (; itrClust != itrClustEnd; itrClust++) {
    if (DEBUG_CLUSTERS) DEBUG_VAR(itrClust->first);
    if (specResults.find(itrClust->first) != specResults.end()) {
      specResults[itrClust->first].numOrigSpectra = itrClust->second.size();
      if (DEBUG_CLUSTERS) DEBUG_MSG("   Y " << itrClust->first);
      list<pair<int, string> >::iterator itrList = itrClust->second.begin();
      list<pair<int, string> >::iterator itrListEnd = itrClust->second.end();
      for (; itrList != itrListEnd; itrList++) {
        if (DEBUG_CLUSTERS) DEBUG_MSG("       " << itrList->second);
        if (itrList->second.find("TCGA") != string::npos) {
          specResults[itrClust->first].cancer = true;
        }
        else {
          specResults[itrClust->first].normal = true;
        }
      }
      if (DEBUG_CLUSTERS) DEBUG_MSG("XX " << itrClust->first << "  "
          << specResults[itrClust->first].normal << " "
          << specResults[itrClust->first].cancer);
    }
    else {
      if (DEBUG_CLUSTERS) DEBUG_MSG("   N " << itrClust->first);
    }
  }
#endif

  DEBUG_MSG("Creating map of contigs to annos...");

  //----------------------------------------------------------
  // Create the map of contigs to their annotations
  //----------------------------------------------------------
  map<int, set<tuple<string, string, float> > > mapContigToTupleAnnoProtScore;
  for (int iPSM = 0; iPSM < psmSetContigsAll.size(); iPSM++) {
    psmPtr psmSpec = psmSetContigsAll[iPSM];
    mapContigToTupleAnnoProtScore[psmSpec->m_scanNum].insert(make_tuple<string,
        string, float>(psmSpec->m_annotation,
                       psmSpec->m_protein,
                       psmSpec->m_score));
  }

  DEBUG_MSG("Populating contig fields...");

  //----------------------------------------------------------
  // Populate some contig fields
  // We use FDR contigs here because they are singular
  //----------------------------------------------------------
  for (int iPSM = 0; iPSM < psmSetContigFdr.size(); iPSM++) {
    psmPtr psmSpec = psmSetContigFdr[iPSM];
    // Ignore decoys
    if (psmSpec->m_isDecoy || psmSpec->m_protein.find("XXX") != string::npos) {
      continue;
    }
    contigResults[psmSpec->m_scanNum].passFdr = true;
    contigResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
    contigResults[psmSpec->m_scanNum].actualScore = psmSpec->m_score;
    contigResults[psmSpec->m_scanNum].fdr = psmSpec->m_fdr;
    contigResults[psmSpec->m_scanNum].orient = psmSpec->m_matchOrientation;
    contigResults[psmSpec->m_scanNum].annotation = psmSpec->m_annotation;
    contigResults[psmSpec->m_scanNum].protein = psmSpec->m_protein;
  }

  DEBUG_MSG("Populating more contig fields...");

  //----------------------------------------------------------
  // Populate more contig fields
  //   Now we use all contigs and pick highest scoring
  //   And don't replace anything from the FDR contigs
  //----------------------------------------------------------
  for (int iPSM = 0; iPSM < psmSetContigsAll.size(); iPSM++) {
    psmPtr psmSpec = psmSetContigsAll[iPSM];
    if (contigResults.find(psmSpec->m_scanNum) == contigResults.end()) {
      contigResults[psmSpec->m_scanNum].passFdr = false;
      contigResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
      contigResults[psmSpec->m_scanNum].orient = psmSpec->m_matchOrientation;
      contigResults[psmSpec->m_scanNum].annotation = psmSpec->m_annotation;
      contigResults[psmSpec->m_scanNum].protein = psmSpec->m_protein;
      contigResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
    }
    else if (psmSpec->m_score > contigResults[psmSpec->m_scanNum].alignScore) {
      contigResults[psmSpec->m_scanNum].passFdr = false;
      contigResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
      contigResults[psmSpec->m_scanNum].orient = psmSpec->m_matchOrientation;
      contigResults[psmSpec->m_scanNum].annotation = psmSpec->m_annotation;
      contigResults[psmSpec->m_scanNum].protein = psmSpec->m_protein;
      contigResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
    }
  }

  // Load amino acid masses
  AAJumps jumps(2);
  jumps.loadJumps(aafile.c_str());

  float peakTol = commandLineParams.getValueFloat("PEAK_TOLERANCE", 0.4);
  DEBUG_VAR(peakTol);

#if 0
  DEBUG_MSG("Finding denovo annotations for contig spectra...");

  //----------------------------------------------------------
  // Find the denovo annotation for the contig spectra
  //----------------------------------------------------------
  for (int iSpec = 0; iSpec < contigSpectra.size(); iSpec++ ) {

    vector<vector<int> > neighs;
    getNeighborsL(contigSpectra[iSpec], jumps, peakTol, neighs);
    DEBUG_VAR(neighs.size());

    list<Spectrum> seqs;
    list<float> seqScores;
    float minPercMaxScore = 0.99;

    denovo_LtoR(contigSpectra[iSpec], neighs, seqs, seqScores, minPercMaxScore);
    DEBUG_VAR(seqs.size());
    DEBUG_VAR(seqScores.size());

    if (seqs.size() == 0) {
      continue;
    }
    Spectrum & spec = *seqs.begin();

    string sequence;
    jumps.getPeptideFromSpectrum(spec, sequence, peakTol);
    DEBUG_VAR(sequence);
    contigResults[contigSpectra[iSpec].scan].denovoAnno = sequence;
  }
#endif  

  DEBUG_MSG("Create the full PSM map...");

  //----------------------------------------------------------
  // Create the full PSM map
  //----------------------------------------------------------
  map<int, psmPtr> mapScanFullPsm;
  for (int iPSM = 0; iPSM < psmSetSpecFullTgt.size(); iPSM++) {
    psmPtr psmSpec = psmSetSpecFullTgt[iPSM];

    psmSpec->changeGapAnnosToSingle();
    psmSpec->addFixedCysteineMods();

    if (mapScanFullPsm.find(psmSpec->m_scanNum) == mapScanFullPsm.end()) {
      mapScanFullPsm[psmSpec->m_scanNum] = psmSpec;
      specResults[psmSpec->m_scanNum].specAnno = psmSpec->m_annotation;
      specResults[psmSpec->m_scanNum].specAnnoOrig = psmSpec->m_origAnnotation;
      specResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
      specResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
    }
    else {
      if (psmSpec->m_score > mapScanFullPsm[psmSpec->m_scanNum]->m_score) {
        mapScanFullPsm[psmSpec->m_scanNum] = psmSpec;
        specResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
        specResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
      }
    }

    if (specResults[psmSpec->m_scanNum].specProb < 0.0) {
      specResults[psmSpec->m_scanNum].specProb = 1.0;
    }

  }
  for (int iPSM = 0; iPSM < psmSetSpecFullDec.size(); iPSM++) {
    psmPtr psmSpec = psmSetSpecFullDec[iPSM];

    psmSpec->changeGapAnnosToSingle();
    psmSpec->addFixedCysteineMods();

    if (mapScanFullPsm.find(psmSpec->m_scanNum) == mapScanFullPsm.end()) {
      mapScanFullPsm[psmSpec->m_scanNum] = psmSpec;
      specResults[psmSpec->m_scanNum].specAnno = psmSpec->m_annotation;
      specResults[psmSpec->m_scanNum].specAnnoOrig = psmSpec->m_origAnnotation;
      specResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
      specResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
      specResults[psmSpec->m_scanNum].isDecoy = true;
    }
    else {
      if (psmSpec->m_score > mapScanFullPsm[psmSpec->m_scanNum]->m_score) {
        mapScanFullPsm[psmSpec->m_scanNum] = psmSpec;
        specResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
        specResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
        specResults[psmSpec->m_scanNum].isDecoy = true;
      }
    }

    if (specResults[psmSpec->m_scanNum].specProb < 0.0) {
      specResults[psmSpec->m_scanNum].specProb = 1.0;
    }

  }

  DEBUG_MSG("Populating one spectrum pass FDR flag...");

  //----------------------------------------------------------
  // Set one spectrum pass FDR contig flag
  //----------------------------------------------------------
  for (int iPSM = 0; iPSM < psmSetSpecFdr.size(); iPSM++) {
    psmPtr psmFdr = psmSetSpecFdr[iPSM];
    // Ignore decoys
    if (psmFdr->m_isDecoy || psmFdr->m_protein.find("XXX") != string::npos) {
      continue;
    }
    psmFdr->changeGapAnnosToSingle();
    //psmFdr->addFixedCysteineMods();
    specResults[psmFdr->m_scanNum].specAnno = psmFdr->m_annotation;
    //DEBUG_VAR(specResults[psmFdr->m_scanNum].specAnno);
    specResults[psmFdr->m_scanNum].specAnnoOrig = psmFdr->m_origAnnotation;
    specResults[psmFdr->m_scanNum].passFdr = true;
    specResults[psmFdr->m_scanNum].fdr = psmFdr->m_fdr;
    specResults[psmFdr->m_scanNum].alignScore = psmFdr->m_score;
    specResults[psmFdr->m_scanNum].specProb = psmFdr->m_pValue;
    int contigIndex = specToContig[psmFdr->m_scanNum];
    if (contigResults[contigIndex].oneSpecPassFdr == false) {
      contigResults[contigIndex].oneSpecPassFdr = true;
      totalContigsWithGoodFdrSpecs++;
    }

    if (specResults[psmFdr->m_scanNum].specProb < 0.0 ||
        specResults[psmFdr->m_scanNum].specProb > 1.0) {
      specResults[psmFdr->m_scanNum].specProb = 1.0;
    }

  }

  DEBUG_MSG("Populate some spectrum fields...");

  //----------------------------------------------------------
  // Populate some spectrum fields
  //----------------------------------------------------------
  map<int, psmPtr>::iterator itrp = mapScanFullPsm.begin();
  map<int, psmPtr>::iterator itrpEnd = mapScanFullPsm.end();
  for (; itrp != itrpEnd; itrp++) {
    psmPtr psmSpec = itrp->second;

//    if (!psmSpec->m_isDecoy) {
//      addFlankingAas(psmSpec, dbFasta);
//    }

    psmSpec->changeGapAnnosToSingle();
    psmSpec->addFixedCysteineMods();

    if (specResults[psmSpec->m_scanNum].specAnno.empty()) {
      specResults[psmSpec->m_scanNum].specAnno = psmSpec->m_annotation;
    }
    if (specResults[psmSpec->m_scanNum].specAnnoOrig.empty()) {
      specResults[psmSpec->m_scanNum].specAnnoOrig = psmSpec->m_origAnnotation;
    }
    PeptideSpectrumMatch::getUnmodifiedPeptide(specResults[psmSpec->m_scanNum].specAnno,
                                               specResults[psmSpec->m_scanNum].specAnnoClean);
    specResults[psmSpec->m_scanNum].specAnnoLength =
        specResults[psmSpec->m_scanNum].specAnnoClean.length();
    specResults[psmSpec->m_scanNum].specProtein = psmSpec->m_protein;
    specResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
    specResults[psmSpec->m_scanNum].orient = psmSpec->m_matchOrientation;

    specResults[psmSpec->m_scanNum].specTotalMass =
        jumps.getPeptideMass(psmSpec->m_annotation);

    if (specResults[psmSpec->m_scanNum].specProb < 0.0) {
      specResults[psmSpec->m_scanNum].specProb = 1.0;
    }

    vector<float> modifications;
    vector<unsigned int> positions;
    vector<unsigned int> lengths;
    // get mods from ORIGINAL (true = useOriginal)
    psmSpec->getModificationsAndPositions(modifications,
                                          positions,
                                          lengths,
                                          true);
    int numMods = modifications.size();
    if (positions.size() > 0) {
      if (positions[0] == 0 && fabs(modifications[0]) == 1.0) {
        numMods--;
      }
      if (positions.size() > 1 && lengths[lengths.size() - 1] == 0
          && fabs(modifications[positions.size() - 1]) == 1.0) {
        numMods--;
      }

      for (int i = 0; i < (int)positions.size() - 1; i++) {
        if ((positions[i] >= positions[i + 1] - lengths[i + 1])
            && (modifications[i] + modifications[i + 1] < 2.0)
            && (modifications[i] + modifications[i + 1] > -2.0)) {
          //DEBUG_MSG(psmSpec->m_scanNum << "  " << positions[i] << "  " << positions[i+1] << "  " << lengths[i] << "  " << lengths[i+1])
          specResults[psmSpec->m_scanNum].maxSpecGap = max(lengths[i],
                                                           lengths[i + 1]);
          if (lengths[i] + lengths[i + 1] >= 8 && lengths[i] != 1
              && lengths[i + 1] != 1) {
            specResults[psmSpec->m_scanNum].suspiciousMods = 1;
          }
        }
      }

    }
    specResults[psmSpec->m_scanNum].numMods = numMods;
  }

  DEBUG_MSG("Create MSGF ID map...");

  //----------------------------------------------------------
  // Create MSGF ID map
  //----------------------------------------------------------
  for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++) {

    psmPtr psmMsgf = psmSetMsgf[iPSM];
    mapScanMsgf[psmMsgf->m_scanNum] = psmMsgf;

    if (specToContig.find(psmMsgf->m_scanNum) != specToContig.end()) {
      stars_in_contigs_msgf_id++;
    }

    specResults[psmMsgf->m_scanNum].msgfAnno = psmMsgf->m_annotation;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmMsgf->m_annotation,
                                               specResults[psmMsgf->m_scanNum].msgfAnnoClean);
    specResults[psmMsgf->m_scanNum].msgfProtein = psmMsgf->m_protein;

    specResults[psmMsgf->m_scanNum].msgfProb = psmMsgf->m_pValue;
    /*
     if (psmMsgf->m_score > 0) {
     specResults[psmMsgf->m_scanNum].msgfProb = -log10(psmMsgf->m_score);
     } else if (psmMsgf->m_pValue > 0) {
     specResults[psmMsgf->m_scanNum].msgfProb = psmMsgf->m_pValue;
     }
     */
    specResults[psmMsgf->m_scanNum].msgfTotalMass =
        jumps.getPeptideMass(psmMsgf->m_annotation);

    if (specResults[psmMsgf->m_scanNum].alignScore == 0.0
        || specResults[psmMsgf->m_scanNum].msgfProb == 0.0) {
      specResults[psmMsgf->m_scanNum].spChange = 0.0;
    }
    else {
      specResults[psmMsgf->m_scanNum].spChange =
          specResults[psmMsgf->m_scanNum].msgfProb
              - specResults[psmMsgf->m_scanNum].specProb;
    }

  } // for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ )

  int nUnkonwnContigs = contigAbinfo.size();

  DEBUG_MSG("Finding contig/spectra tags...");
  //----------------------------------------------------------
  // Figure out if each contig had a tag that matched an MSGF 
  //          annotation for a spectra in this contig
  //----------------------------------------------------------
  for (int iTagPSM = 0; iTagPSM < psmSetContigTag.size(); iTagPSM++) {
    psmPtr psmTagContig = psmSetContigTag[iTagPSM];
    int contigIndex = psmTagContig->m_scanNum;
    if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
      DEBUG_VAR(contigIndex);

    string cleanTag = convertAnnotation(psmTagContig->m_origAnnotation);
    if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
      DEBUG_MSG("  " << cleanTag);

    string cleanContig;
    PeptideSpectrumMatch::getUnmodifiedPeptide(contigResults[contigIndex].annotation,
                                               cleanContig);
    if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
      DEBUG_MSG("  " << cleanContig);

    //DEBUG_MSG(psmTagContig->m_scanNum << "  " << psmTagContig->m_protein);
    set<int>::iterator itrs = contigToSpecSet[contigIndex].begin();
    set<int>::iterator itrsEnd = contigToSpecSet[contigIndex].end();
    for (; itrs != itrsEnd; itrs++) {

      int specIndex = *itrs;
      if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
        DEBUG_MSG("  " << specIndex);
      if (mapScanMsgf.find(specIndex) == mapScanMsgf.end()) {
        continue;
      }

      string msgfAnnoCleanConverted =
          convertAnnotation(specResults[specIndex].msgfAnnoClean);
      if (msgfAnnoCleanConverted.find(cleanTag) == string::npos) {
        continue;
      }

      if (!contigResults[contigIndex].tagMatchGood) {
        if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
          DEBUG_MSG("    YES");
        contigResults[contigIndex].tagMatchGood = true;
        contigResults[contigIndex].tag = cleanTag;
        totalContigsWithGoodTags++;
      }

      //----------------------------------------------------------
      // Figure out if each contig had a alignment that matched an MSGF 
      //          annotation for a spectra in this contig
      //----------------------------------------------------------
      set<tuple<string, string, float> >::iterator itrs2 =
          mapContigToTupleAnnoProtScore[contigIndex].begin();
      set<tuple<string, string, float> >::iterator itrs2End =
          mapContigToTupleAnnoProtScore[contigIndex].end();
      for (; itrs2 != itrs2End; itrs2++) {

        string contigProt = itrs2->m1;
        if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
          DEBUG_MSG("  " << contigProt);

        // Proteins must match also
        if (contigProt != specResults[specIndex].msgfProtein) {
          continue;
        }

        string contigAnno;
        PeptideSpectrumMatch::getUnmodifiedPeptide(string(itrs2->m0),
                                                   contigAnno);
        if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
          DEBUG_MSG("  " << contigAnno);
        string cleanAnno = convertAnnotation(contigAnno);
        if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
          DEBUG_MSG("  " << cleanAnno);

        // Locate the msgf annotation in the daatbase
        string dbProtein =
            convertAnnotation(string(dbFasta[(char *)specResults[specIndex].msgfProtein.c_str()]));
        if (dbProtein.empty()) {
          ERROR_MSG(specResults[specIndex].msgfProtein << "not found");
          continue;
        }

        // Locate the msgf annotation in the database protein
        unsigned int msgfAnnoLoc = dbProtein.find(msgfAnnoCleanConverted);
        if (msgfAnnoLoc == string::npos) {
          continue;
        }
        if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
          DEBUG_VAR(msgfAnnoLoc);
        unsigned int msgfAnnoLocEnd = msgfAnnoLoc
            + specResults[specIndex].msgfAnnoClean.size() - 1;
        if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
          DEBUG_VAR(msgfAnnoLocEnd);

        // Locate our annotation in the database protein
        unsigned int ourAnnoLoc = dbProtein.find(cleanAnno);
        if (ourAnnoLoc == string::npos) {
          continue;
        }
        if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
          DEBUG_VAR(ourAnnoLoc);
        unsigned int ourAnnoLocEnd = ourAnnoLoc + cleanContig.size();
        if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
          DEBUG_VAR(ourAnnoLocEnd);

        if ((msgfAnnoLoc >= ourAnnoLoc && msgfAnnoLoc <= ourAnnoLocEnd)
            || (msgfAnnoLocEnd >= ourAnnoLoc && msgfAnnoLocEnd <= ourAnnoLocEnd)
            || (msgfAnnoLoc < ourAnnoLoc && msgfAnnoLocEnd > ourAnnoLocEnd)) {
          contigResults[contigIndex].contigPsmGood = true;
          contigResults[contigIndex].actualScore = itrs2->m2;
          if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
            DEBUG_MSG("    YES2");
          break; // Only need one
        }

        if (msgfAnnoCleanConverted.find(cleanAnno) != string::npos) {
          contigResults[contigIndex].contigPsmGood = true;
          //DEBUG_MSG("    YES2");
          break;// Only need one
        }

      } // for ( ; itrs2 != itrs2End ; itrs2++) {

      if (contigResults[contigIndex].actualScore == 0) {

        //----------------------------------------------------------
        // Figure out if each contig had a alignment that matched an MSGF 
        //          annotation for a spectra in this contig
        //----------------------------------------------------------
        set<tuple<string, string, float> >::iterator itrs2 =
            mapContigToTupleAnnoProtScore2[contigIndex].begin();
        set<tuple<string, string, float> >::iterator itrs2End =
            mapContigToTupleAnnoProtScore2[contigIndex].end();
        for (; itrs2 != itrs2End; itrs2++) {

          string contigProt = itrs2->m1;
          if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
            DEBUG_MSG("  " << contigProt);

          // Proteins must match also
          if (contigProt != specResults[specIndex].msgfProtein) {
            continue;
          }

          string contigAnno;
          PeptideSpectrumMatch::getUnmodifiedPeptide(string(itrs2->m0),
                                                     contigAnno);
          if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
            DEBUG_MSG("  " << contigAnno);
          string cleanAnno = convertAnnotation(contigAnno);
          if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
            DEBUG_MSG("  " << cleanAnno);

          // Locate the msgf annotation in the daatbase
          string dbProtein =
              convertAnnotation(string(dbFasta[(char *)specResults[specIndex].msgfProtein.c_str()]));
          if (dbProtein.empty()) {
            ERROR_MSG(specResults[specIndex].msgfProtein << "not found");
            continue;
          }

          // Locate the msgf annotation in the database protein
          unsigned int msgfAnnoLoc = dbProtein.find(msgfAnnoCleanConverted);
          if (msgfAnnoLoc == string::npos) {
            continue;
          }
          if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
            DEBUG_VAR(msgfAnnoLoc);
          unsigned int msgfAnnoLocEnd = msgfAnnoLoc
              + specResults[specIndex].msgfAnnoClean.size() - 1;
          if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
            DEBUG_VAR(msgfAnnoLocEnd);

          // Locate our annotation in the database protein
          unsigned int ourAnnoLoc = dbProtein.find(cleanAnno);
          if (ourAnnoLoc == string::npos) {
            continue;
          }
          if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
            DEBUG_VAR(ourAnnoLoc);
          unsigned int ourAnnoLocEnd = ourAnnoLoc + cleanContig.size();
          if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
            DEBUG_VAR(ourAnnoLocEnd);

          if ((msgfAnnoLoc >= ourAnnoLoc && msgfAnnoLoc <= ourAnnoLocEnd)
              || (msgfAnnoLocEnd >= ourAnnoLoc
                  && msgfAnnoLocEnd <= ourAnnoLocEnd)
              || (msgfAnnoLoc < ourAnnoLoc && msgfAnnoLocEnd > ourAnnoLocEnd)) {
            contigResults[contigIndex].actualScore = itrs2->m2;
            if (DEBUG_TAG_MATCH && DEBUG_CONTIG_INDEX == contigIndex)
              DEBUG_MSG("    YES2");
            break;
          }

        } // for ( ; itrs2 != itrs2End ; itrs2++) {

      }

    }

  } // for (int iTagPSM = 0; iTagPSM < psmSetContigTag.size(); iTagPSM++ ) {

  DEBUG_VAR(totalContigsWithGoodTags);

  int nNoMatch = 0;
  int nDiffProtein = 0;
  int nExactMatch = 0;
  int nMatchG80 = 0;
  int nDiffPeptide = 0;

  DEBUG_MSG("Matching MSGF results...");
  //---------------------------------------------------------------------
  // Determine if the spectra from specnets match MSGFDB results
  //---------------------------------------------------------------------
  map<int, psmPtr>::iterator itrFull = mapScanFullPsm.begin();
  map<int, psmPtr>::iterator itrFullEnd = mapScanFullPsm.end();
  for (; itrFull != itrFullEnd; itrFull++) {

    const psmPtr & psmSpec = itrFull->second;
    int specIndex = psmSpec->m_scanNum;
    //DEBUG_VAR(specIndex);

    specPSMSet.insert(specIndex);
    int contigIndex = specToContig[specIndex];
    //DEBUG_VAR(contigIndex);

    string temp;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSpec->m_annotation, temp);
    string cleanAnnotationSpec = convertAnnotation(temp);
    //DEBUG_VAR(cleanAnnotationSpec);

    map<int, psmPtr>::iterator itrFind = mapScanMsgf.find(specIndex);
    if (itrFind == mapScanMsgf.end()) {
      if (specResults[specIndex].passFdr)
        setNoMsgf.insert(specIndex);
      continue;
    }

    psmPtr psmMsgf = itrFind->second;
    string temp2;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmMsgf->m_annotation, temp2);
    string cleanAnnotationMsgf = convertAnnotation(temp2);
    //DEBUG_VAR(cleanAnnotationMsgf);

    mapSpectrumToGoodMsgfID[specIndex] = 0;

    map<int, set<string> >::iterator itrTag =
        mapContigToProteinSet.find(contigIndex);
    if (itrTag != mapContigToProteinSet.end()) {
      set<string>::iterator itrProtein = itrTag->second.begin();
      set<string>::iterator itrProteinEnd = itrTag->second.end();
      for (; itrProtein != itrProteinEnd; itrProtein++) {
        if (itrProtein->find(psmMsgf->m_protein) != string::npos) {
          //DEBUG_MSG("  " << "found protein [" << contigIndex << "] protein [" << *itrProtein << "]");
          mapSpectrumToGoodMsgfID[specIndex] = 1;
          break;
        }
      }
    }

    float percentMatch = getPercentMatch(cleanAnnotationSpec,
                                         cleanAnnotationMsgf);
    specResults[specIndex].matchPercent = percentMatch;

    if (percentMatch >= MATCH_THRESHOLD) {
      //DEBUG_MSG("MATCH");
      nExactMatch++;
      if (specResults[specIndex].passFdr) {
        setExactMatch.insert(specIndex);
      }
      mapSpecToMatchType[specIndex] = 1;
      mapSpectrumToGoodMsgfID[specIndex] = 1;
      continue;
    }

    if (specResults[specIndex].passFdr) {
      //DEBUG_MSG("NON MATCH");
      setDiffPeptide.insert(specIndex);
    }
    mapSpecToMatchType[specIndex] = 3;
    //DEBUG_MSG(specIndex << "  " << mapSpecToMatchType[specIndex]);

  } // for (int iPSM = 0; iPSM < psmSetSpecFull.size(); iPSM++ )

  DEBUG_MSG("Finding similar spectra...");
  //---------------------------------------------------------------------
  // Find similar spectra in contig to figure out if unknowns are good or bad
  //    and if MSGFDB got the wrong answer
  //---------------------------------------------------------------------
  // Find spectra that can be matched to others in the contig
  map<int, psmPtr>::iterator itrm = mapScanFullPsm.begin();
  map<int, psmPtr>::iterator itrmEnd = mapScanFullPsm.end();
  for (; itrm != itrmEnd; itrm++) {
    const psmPtr & psmSpec = itrm->second;

//  for (int iPSM = 0; iPSM < psmSetSpecFull.size(); iPSM++ ) {
//    const psmPtr & psmSpec = psmSetSpecFull[iPSM];

    int specIndex = psmSpec->m_scanNum;
    int contigIndex = specToContig[specIndex];
    int matchType = mapSpecToMatchType[specIndex];

    string temp;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psmSpec->m_annotation, temp);
    string cleanAnnotationSpec1 = convertAnnotation(temp);
    // Only compare to previously uncategorized spectra (either match or non-match)
    if (matchType != 1 && matchType != 3) {
      continue;
    }

    map<int, psmPtr>::iterator itrm2 = mapScanFullPsm.begin();
    map<int, psmPtr>::iterator itrm2End = mapScanFullPsm.end();
    for (; itrm2 != itrm2End; itrm2++) {
      const psmPtr & psmSpec2 = itrm2->second;

//    for (int iPSM2 = 0; iPSM2 < psmSetSpecFull.size(); iPSM2++ ) {
//      const psmPtr & psmSpec2 = psmSetSpecFull[iPSM2];

      int specIndex2 = psmSpec2->m_scanNum;
      if (specIndex == specIndex2) {
        continue;
      }
      int contigIndex2 = specToContig[specIndex2];
      if (contigIndex != contigIndex2) {
        continue;
      }

      // Only set previously uncategorized spectra
      int matchType2 = mapSpecToMatchType[specIndex2];
      if (matchType2 == 1 || matchType2 == 3 || matchType2 == 4) {
        continue;
      }

      string temp;
      PeptideSpectrumMatch::getUnmodifiedPeptide(psmSpec2->m_annotation, temp);
      string cleanAnnotationSpec2 = convertAnnotation(temp);
      float percentMatch = getPercentMatch(cleanAnnotationSpec1,
                                           cleanAnnotationSpec2);
      //DEBUG_MSG(contigIndex << "  " << contigIndex2 << "  " << specIndex << "  " << specIndex2);
      //DEBUG_MSG("  " << cleanAnnotationSpec1 << "  " << cleanAnnotationSpec2 << "  " << percentMatch << "  " << matchType << "  " << matchType2);
      specResults[specIndex2].matchPercent = -percentMatch;
      if (percentMatch > SIM_MATCH_THRESHOLD) {
        if (matchType == 1) {
          if (specResults[specIndex2].passFdr)
            nSimGood++;
          mapSpecToMatchType[specIndex2] = 4;
          mapSpectrumToGoodMsgfID[specIndex2] = 1;
        }
        else if (matchType == 3) {
          if (specResults[specIndex2].passFdr)
            nSimBad++;
          mapSpecToMatchType[specIndex2] = 5;
        }
      }
    }
  }

  DEBUG_MSG("Finding contig hybrids...");
  //---------------------------------------------------------------------
  // Figure out if the contigs are good, bad or hybrid
  //---------------------------------------------------------------------

  map<int, int>::iterator itrs = specToContig.begin();
  map<int, int>::iterator itrsEnd = specToContig.end();
  for (; itrs != itrsEnd; itrs++) {

    int specIndex = itrs->first;
    int contigIndex = itrs->second;

    map<int, psmPtr>::iterator itrFindMsgf = mapScanMsgf.find(specIndex);
    if (itrFindMsgf == mapScanMsgf.end()) {
      map<int, string>::iterator itrFindProt =
          mapContigToMsgfProtein.find(contigIndex);
      if (itrFindProt == mapContigToMsgfProtein.end()) {
        mapContigMsgfHybrid[contigIndex] = -1;
      }
      continue;
    }

    const psmPtr & psmSpec = itrFindMsgf->second;

    map<int, string>::iterator itrFindProt =
        mapContigToMsgfProtein.find(contigIndex);
    if (itrFindProt == mapContigToMsgfProtein.end()) {
      mapContigToMsgfProtein[contigIndex] = psmSpec->m_protein;
      mapContigMsgfHybrid[contigIndex] = 0;
    }
    else {
      string prevProtein = itrFindProt->second;
      if (itrFindMsgf->second->m_protein != prevProtein) {
        mapContigMsgfHybrid[contigIndex] = 1;
      }
    }

  }

  //---------------------------------------------------------------------
  // Figure out if the contigs are good, bad or hybrid
  //---------------------------------------------------------------------
  map<int, vector<int> > mapContigVecType;
  map<int, int> mapContigHasTag;

  itrm = mapScanFullPsm.begin();
  for (; itrm != itrmEnd; itrm++) {
    const psmPtr & psmSpec = itrm->second;

//  for (int iPSM = 0; iPSM < psmSetSpecFull.size(); iPSM++ ) {
//    const psmPtr & psmSpec = psmSetSpecFull[iPSM];

    int specIndex = psmSpec->m_scanNum;
    int contigIndex = specToContig[specIndex];
    int matchType = mapSpecToMatchType[specIndex];
    int contigType = 0;
    if (matchType == 1 || matchType == 2 || matchType == 4 || matchType == 6)
      contigType = 1;
    if (matchType == 3 || matchType == 5)
      contigType = 3;
    mapContigVecType[contigIndex].push_back(contigType);
  }

  map<int, vector<int> >::iterator itrv = mapContigVecType.begin();
  map<int, vector<int> >::iterator itrvEnd = mapContigVecType.end();
  for (; itrv != itrvEnd; itrv++) {
    int contigIndex = itrv->first;
    vector<int> & vecType = itrv->second;
    for (int i = 0; i < vecType.size(); i++) {
      if (vecType[i] != 0) {
        if (mapContigType.find(contigIndex) == mapContigType.end()) {
          mapContigType[contigIndex] = vecType[i];
        }
        else if (mapContigType[contigIndex] != vecType[i]) {
          mapContigType[contigIndex] = 8;
        }
      }
    }
    if (mapContigType[contigIndex] == 1)
      nGoodContig++;
    if (mapContigType[contigIndex] == 3)
      nBadContig++;
    if (mapContigType[contigIndex] == 8)
      nHybridContig++;
  }

  map<int, set<int> > mapContigVarSet;
  map<int, set<int> > mapContigVarSetMsgf;

  //---------------------------------------------------------------------
  // Output result tables
  //---------------------------------------------------------------------
  DEBUG_MSG("Outputting results tables...");

  if (!commandLineParams.exists("RESULTS_OUTPUT_DIR")) {
    ofstream ofs(argv[3]);
    if (!ofs || !ofs.good()) {
      ERROR_MSG("Unable to open file [" << argv[3] << "]");
      return false;
    }
    outputSpectrumLevelResults(ofs);
    ofs << endl;
    outputContigLevelResults(ofs);
    ofs << endl;
    outputSummaryResults(ofs);
    ofs.close();
  }
  else {

    string outDir = commandLineParams.getValue("RESULTS_OUTPUT_DIR");
    outDir += "/";

    string filename = outDir + "spectrum_level_results_table.txt";
    DEBUG_MSG("Saving contig results table to [" << filename << "]");
    ofstream specofs(filename.c_str());
    if (!specofs || !specofs.good()) {
      ERROR_MSG("Unable to open file [" << filename << "]");
      return false;
    }
    outputSpectrumLevelResults(specofs);
    specofs.close();

    filename = outDir + "contig_level_results_table.txt";
    DEBUG_MSG("Saving contig results table to [" << filename << "]");
    ofstream contigofs(filename.c_str());
    if (!contigofs || !contigofs.good()) {
      ERROR_MSG("Unable to open file [" << filename << "]");
      return false;
    }
    outputContigLevelResults(contigofs);
    contigofs.close();

    filename = outDir + "summary_results_table.txt";
    DEBUG_MSG("Saving contig results table to [" << filename << "]");
    ofstream summaryofs(filename.c_str());
    if (!summaryofs || !summaryofs.good()) {
      ERROR_MSG("Unable to open file [" << filename << "]");
      return false;
    }
    outputSummaryResults(summaryofs);
    summaryofs.close();
  }

  DEBUG_MSG("Done.");
  return 0;
}

