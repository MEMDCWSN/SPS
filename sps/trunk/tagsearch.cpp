#include "spectrum.h"
#include "inputParams.h"
#include "batch.h"
#include "alignment_scoring.h"
#include "db_fasta.h"
#include "tags.h"
#include "tuple.h"
#include <ctime>

using namespace specnets;

struct ProteinMatch
{
  int protIdx;
  float startMass;
  float score, percScore;
  string peptide;

  ProteinMatch()
  {
    protIdx = 0;
    startMass = 0.0;
    score = 0.0;
    percScore = 0.0;
    peptide = "";
  }
  void set(int inProtIdx,
           float inStartMass,
           float inScore,
           float inPercScore,
           string & inPeptide)
  {
    protIdx = inProtIdx;
    startMass = inStartMass;
    score = inScore;
    percScore = inPercScore;
    peptide = inPeptide;
  }
  ProteinMatch &operator=(const ProteinMatch &other)
  {
    protIdx = other.protIdx;
    startMass = other.startMass;
    score = other.score;
    percScore = other.percScore;
    peptide = other.peptide;
    return *this;
  }
};

int main(int argc, char **argv)
{
  InputParams params;
  vector<const char *> paramStrs;
  paramStrs.resize(2);
  paramStrs[0] = "INPUT_FASTA";
  paramStrs[1] = "OUTPUT_PEPTIDES";

  if (argc == 1)
    params.readParams("tagsearch.params");
  else
    params.readParams(argv[1]);
  if (!params.confirmParams(paramStrs)) {
    cerr << "ERROR: Parameters file ";
    if (argc == 1)
      cerr << "tagsearch.params";
    else
      cerr << argv[1];
    cerr
        << " is incomplete. One of the following is missing: INPUT_FASTA, OUTPUT_PEPTIDES\n";
    return -1;
  }

  const char *dbFN = params.getValue(paramStrs[0]);
  const char *pepsFN = params.getValue(paramStrs[1]);

  int tagLen =
      params.paramPresent("TAG_LEN") ? (int) params.getValueInt("TAG_LEN") : 6;
  int numJumps =
          params.paramPresent("DOUBLE_AA_JUMPS") ? (int) params.getValueInt("DOUBLE_AA_JUMPS") : 1;
  short minMatchFlanking =
          params.paramPresent("MATCH_TAG_FLANKING_MASSES") ? (short) params.getValueInt("MATCH_TAG_FLANKING_MASSES") : 1;
  unsigned int maxNumTags =
          params.paramPresent("MAX_NUM_TAGS") ? (unsigned int) params.getValueInt("MAX_NUM_TAGS") : 0;
  float peakTol =
          params.paramPresent("TOLERANCE_PEAK") ? (float) params.getValueDouble("TOLERANCE_PEAK") : 0.5;
  float pmTol =
          params.paramPresent("TOLERANCE_PM") ? (float) params.getValueDouble("TOLERANCE_PM") : 1;
  int specType =
          params.paramPresent("SPEC_TYPE_MSMS") ? ((int) params.getValueInt("SPEC_TYPE_MSMS") ? 1 : 0) : 0;
  bool maxParsimony =
          (bool) params.paramPresent("MAX_PARSIMONY") ? ((int) params.getValueInt("MAX_PARSIMONY") ? 1 : 0) : 0;
  if (params.paramPresent("AMINO_ACID_MASSES")) {
    AAJumps jumps(-1);
    jumps.loadJumps(params.getValue("AMINO_ACID_MASSES"), true);
  }
  unsigned short minPercEI =
      params.paramPresent("MIN_PERC_EXPINT") ? (unsigned short) (10000
          * params.getValueDouble("MIN_PERC_EXPINT")) : 0;
  unsigned short minPercTP =
      params.paramPresent("MIN_PERC_TP") ? (unsigned short) (10000
          * params.getValueDouble("MIN_PERC_TP")) : 0;
  unsigned int specIdx, peakIdx, searchIdx;

  cerr << "Got minPercEI = " << minPercEI << ", minPercTP = " << minPercTP
      << endl;

  SpecSet specSet;
  short loadOk = 0;
  if (params.paramPresent("INPUT_SPECS"))
    loadOk = specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
  else if (params.paramPresent("INPUT_SPECS_PKLBIN"))
    loadOk = specSet.loadPklBin(params.getValue("INPUT_SPECS_PKLBIN"));
  if (loadOk <= 0 or specSet.size() == 0) {
    cerr << "Error loading spectra!\n";
    return -1;
  }

  vector<unsigned int> specsToSearch;
  if (params.paramPresent("SUBSET_IDX"))
    Load_binArray(params.getValue("SUBSET_IDX"), specsToSearch);
  else {
    specsToSearch.resize(specSet.size());
    for (searchIdx = 0; searchIdx < specSet.size(); searchIdx++)
      specsToSearch[searchIdx] = searchIdx;
  }

  // Variables for evaluating tagging efficiency
  SpecSet correctPepsAsSpecs;
  vector < list<Tag> > correctTags;
  vector<unsigned int> numCorrectTags;
  if (params.paramPresent("CORRECT_PEPTIDE_SPECS")) {
    correctPepsAsSpecs.loadPklBin(params.getValue("CORRECT_PEPTIDE_SPECS"));
    if (correctPepsAsSpecs.size() != specSet.size())
      correctPepsAsSpecs.resize(0);
    else {
      correctTags.resize(correctPepsAsSpecs.size());
      numCorrectTags.resize(correctPepsAsSpecs.size());
    }
  }

  DB_fasta db;
  if (db.Load(dbFN) <= 0)
    return -1;
  db.replaceAA('L', 'I');
  db.replaceAA('K', 'Q');
  DB_fasta dbOri;
  if (dbOri.Load(dbFN) <= 0)
    return -1;

  BufferedLineReader prevPeps;
  if (params.paramPresent("INPUT_PEPTIDES"))
    prevPeps.Load(params.getValue("INPUT_PEPTIDES"));

  DB_index index(db, 1024, 1024, tagLen);
  FILE *pepsFID = (FILE *) fopen(pepsFN, "wb");
  fprintf(pepsFID,
          "SpecIdx\tPeptide\tTag score\tPerc peptide score\tProtIdx\tProtein\n");

  unsigned int debugIdx = 149; // Set to specSet.size() to skip debug output

  list < Tag > tags;
  AAJumps jumps(1);
  int matchedTags = 0, matchedRevTags = 0;
  char *tagText = (char *) malloc(tagLen + 1);
  tagText[tagLen] = (char) 0;
  list < sps::tuple<int, float, string> > matches; // List of matched (protein index, peptide start, peptide string)
  char buffer2k[2048]; // Buffer for manipulation of peptide strings
  vector < list<ProteinMatch> > proteinHits(specSet.size());
  ProteinMatch curMatch;
  vector<bool> matchedProts(db.size());
  for (unsigned int pIdx = 0; pIdx < matchedProts.size(); pIdx++)
    matchedProts[pIdx] = false;
  vector<float> pepMasses; // Peptide masses derived from peptide sequence
  Spectrum pepSpec; // Peptide spectrum derived from peptide sequence
  for (searchIdx = 0; searchIdx < specsToSearch.size(); searchIdx++) {
    specIdx = specsToSearch[searchIdx];
    cout << "Processing spectrum " << specIdx << "... ";
    cout.flush();

    float specScore = 0, matchTagScore = 0;
    for (peakIdx = 0; peakIdx < specSet[specIdx].size(); peakIdx++)
      specScore += specSet[specIdx][peakIdx][1];
    specSet[specIdx].addZPMpeaks(peakTol, ((float) specType)
        * AAJumps::massHion, true);

    if (specIdx == debugIdx) {
      cerr << "Spectrum " << debugIdx << " has peptide ID "
          << prevPeps.getline(specIdx) << "|\n";
    }

    if (specIdx < prevPeps.size() and strlen(prevPeps.getline(specIdx)) > 0) {
      // Matching previously assigned peptide
      if (strlen(prevPeps.getline(specIdx)) < 2048) {
        const char* constLine = prevPeps.getline(specIdx);
        getPepSeq(constLine, buffer2k);
        dbOri.find(buffer2k, matches); // Search unmodified peptide string
        string oriPep(prevPeps.getline(specIdx));
        for (list<sps::tuple<int, float, string> >::iterator iter =
            matches.begin(); iter != matches.end(); iter++) {
          iter->m2 = oriPep;
          if (specIdx == debugIdx)
            cerr << " -- match (" << iter->m0 << "," << iter->m2 << ")\n";
        }
        if (matches.empty()) {
          matches.push_back(sps::make_tuple<int, float, string>(-1, -1, oriPep));
          if (specIdx == debugIdx)
            cerr << " -- No protein match: (" << matches.front().m0 << ","
                << matches.front().m2 << ")\n";
        }
        matchTagScore = MatchSpecToPeptide(specSet[specIdx],
                                           prevPeps.getline(specIdx),
                                           peakTol,
                                           0,
                                           false);
        if (specIdx == debugIdx)
          cerr << " -- matchTagScore = " << matchTagScore << "\n";
      }
    }
    else {
      // Tag-based matches
      tags.clear();
      proteinHits[specIdx].clear();
      ExtractTags(specSet[specIdx], tags, peakTol, tagLen, numJumps, maxNumTags);
      cout << tags.size() << " tags...";
      cout.flush();

      if (correctPepsAsSpecs.size() > 0) {
        SpecSet curSpec(1);
        curSpec[0] = correctPepsAsSpecs[specIdx];
        vector < list<unsigned int> > matchedTagsIdx;
        vector < list<float> > matchedTagsScores;
        vector < list<unsigned int> > matchedTagsPos;

        MatchTagsToSpectra(curSpec,
                           tags,
                           peakTol,
                           2,
                           0,
                           0,
                           1,
                           matchedTagsIdx,
                           matchedTagsScores,
                           matchedTagsPos);
        numCorrectTags[specIdx] = matchedTagsIdx.size();

        list<Tag>::iterator tagIter = tags.begin();
        list<unsigned int>::iterator matchIter = matchedTagsIdx[0].begin();
        for (unsigned int tagCount = 0; tagCount < tags.size(); tagCount++, tagIter++)
          if (*matchIter == tagCount) {
            correctTags[specIdx].push_back(*tagIter);
            matchIter++;
          }
      }

      matches.clear();
      matchTagScore = 0;
      for (list<Tag>::iterator tagIter = tags.begin(); tagIter != tags.end(); tagIter++) {
        if (tagIter->score < matchTagScore - 0.0001)
          break;
        for (int cIdx = 0; cIdx < tagLen; cIdx++)
          tagText[cIdx] = jumps.aaLetters[tagIter->sequence[cIdx]];
        if (index.find(tagText,
                       db,
                       matches,
                       minMatchFlanking,
                       tagIter->flankingPrefix,
                       peakTol,
                       tagIter->flankingSuffix,
                       pmTol)) {
          matchedTags++;
          matchTagScore = tags.front().score;
        }

        for (int cIdx = 0; cIdx < tagLen; cIdx++)
          tagText[cIdx] = jumps.aaLetters[tagIter->sequence[tagLen - 1 - cIdx]];
        if (index.find(tagText,
                       db,
                       matches,
                       minMatchFlanking,
                       tagIter->flankingSuffix + 18,
                       pmTol,
                       tagIter->flankingPrefix - 18,
                       peakTol)) {
          matchedRevTags++;
          matchTagScore = tags.front().score;
        }
      }
      cout << matches.size() << " matches...";
      cout.flush();

      // Convert Q/I-only matches to Q/K,I/L matches (dbOri)
      for (list<sps::tuple<int, float, string> >::iterator match =
          matches.begin(); match != matches.end(); match++) {
        string prot(db.getSequence(match->m0));
        string::size_type loc = prot.find(match->m2);
        if (loc != string::npos) {
          char *protOri = dbOri.getSequence(match->m0);
          strncpy(buffer2k, &protOri[loc], match->m2.size());
          buffer2k[match->m2.size()] = (char) 0;
          match->m2 = string(buffer2k);
        }
      }
    }

    for (list<sps::tuple<int, float, string> >::iterator match =
        matches.begin(); match != matches.end(); match++) {
      if (maxParsimony or params.paramPresent("OUTPUT_MATCHED_PROTS")) {
        curMatch.set(match->m0, match->m1, matchTagScore, 100 * matchTagScore
            / specScore, match->m2);
        proteinHits[specIdx].push_back(curMatch);
        if (specIdx == debugIdx)
          cerr << " -- proteinHits[specIdx] = " << proteinHits[specIdx].size()
              << " entries\n";
      }
      if (not maxParsimony) {
        fprintf(pepsFID,
                "%d\t%s\t%.1f\t%.1f\t%d",
                specIdx + 1,
                match->m2.c_str(),
                matchTagScore,
                100 * matchTagScore / specScore,
                match->m0);
        fprintf(pepsFID,
                "\t%s %s\n",
                db.getID(match->m0),
                db.getDesc(match->m0));
        matchedProts[match->m0] = true;
      }
    }
    cout << "done\n";
    cout.flush();
  }

  vector < list<pair<int, float> > > simpleProteinHits(specSet.size());
  if (maxParsimony or params.paramPresent("OUTPUT_MATCHED_PROTS")) {
    for (specIdx = 0; specIdx < specSet.size(); specIdx++) {
      simpleProteinHits[specIdx].clear();
      for (list<ProteinMatch>::iterator iter = proteinHits[specIdx].begin(); iter
          != proteinHits[specIdx].end(); iter++)
        if (iter->protIdx >= 0) {
          simpleProteinHits[specIdx].push_back(make_pair<int, float> (iter->protIdx,
                                                                    iter->startMass));
        }
    }
  }
  if (maxParsimony) {
    //cerr<<" -->> simpleProteinHits["<<debugIdx<<"] = "<<simpleProteinHits[debugIdx].size()<<" entries\n"; cerr.flush();
    //cerr<<" -->> proteinHits["<<debugIdx<<"] = "<<proteinHits[debugIdx].size()<<" entries: "<<proteinHits[debugIdx].front().protIdx<<","<<proteinHits[debugIdx].front().peptide<<"\n"; cerr.flush();
    MaximumParsimony( simpleProteinHits);
    //cerr<<" -->> proteinHits["<<debugIdx<<"] = "<<proteinHits[debugIdx].size()<<" entries\n"; cerr.flush();

    for (searchIdx = 0; searchIdx < specsToSearch.size(); searchIdx++) {
      specIdx = specsToSearch[searchIdx];
      if (simpleProteinHits[specIdx].empty())
        continue;
      for (list<ProteinMatch>::iterator iter = proteinHits[specIdx].begin(); iter
          != proteinHits[specIdx].end(); iter++)
        if (simpleProteinHits[specIdx].front().first == iter->protIdx) {
          curMatch = *iter;
          break;
        }
      proteinHits[specIdx].clear();
      proteinHits[specIdx].push_back(curMatch);

      fprintf(pepsFID,
              "%d\t%s\t%.1f\t%.1f\t%d",
              specIdx + 1,
              curMatch.peptide.c_str(),
              curMatch.score,
              curMatch.percScore,
              curMatch.protIdx);
      fprintf(pepsFID,
              "\t%s %s\n",
              db.getID(curMatch.protIdx),
              db.getDesc(curMatch.protIdx));
      matchedProts[curMatch.protIdx] = true;
    }
  }

  // Output sets of matched spectrum/protein indices
  if (params.paramPresent("OUTPUT_MIDX")) {
    SpecSet midx;
    midx.resize(specSet.size());
    vector < vector<int> > proteinMatches(specSet.size()); // Best match[i]: protein / #mods / b/y-match (col 0/1/2)
    vector < vector<unsigned short> > eitp(specSet.size()); // Best match[i]: 1000 * %explained_score / 1000 %matched_ions (col 0/1)

    for (specIdx = 0; specIdx < specSet.size(); specIdx++) {
      proteinMatches[specIdx].resize(3);
      proteinMatches[specIdx][0] = -1; // Assume not matched to any protein
      proteinMatches[specIdx][1] = -1; // Number of modifications is not tracked here
      proteinMatches[specIdx][2] = 0; // Always a b/prefix match (spectrum is reversed otherwise)
      eitp[specIdx].resize(2);
      eitp[specIdx][0] = 0;
      eitp[specIdx][1] = 0;

      if (specIdx == debugIdx)
        cerr << " -->> proteinHits[specIdx] = " << proteinHits[specIdx].size()
            << " entries\n";

      if (not proteinHits[specIdx].empty()) {
        curMatch = proteinHits[specIdx].front();
        proteinMatches[specIdx][0] = curMatch.protIdx;

        vector<int> idxPep;
        float specScore = 0, matchScore = 0;
        for (peakIdx = 0; peakIdx < specSet[specIdx].size(); peakIdx++)
          specScore += specSet[specIdx][peakIdx][1];
        MatchSpecToPeptide(specSet[specIdx],
                           curMatch.peptide.c_str(),
                           peakTol,
                           0,
                           true,
                           &idxPep);
        eitp[specIdx][1]
            = (unsigned short) round(10000.0 * ((float) idxPep.size())
                / ((float) curMatch.peptide.length() + 1));

        if (specIdx == debugIdx)
          cerr << " -->> eitp[specIdx][1] = " << eitp[specIdx][1]
              << ", idxPep.size() = " << idxPep.size() << "\n";

        if (curMatch.protIdx < 0) {
          midx[specIdx].resize(0); // Not matched to a protein sequence
          for (peakIdx = 0; peakIdx < idxPep.size(); peakIdx++)
            matchScore += specSet[specIdx][peakIdx][1];
        }
        else {
          // Get indices for the matched protein masses
          string prot(dbOri.getSequence(curMatch.protIdx));
          string::size_type loc = prot.find(curMatch.peptide);
          if (specIdx == debugIdx) {
            cerr << "Peptide: " << curMatch.peptide << " matches with loc = "
                << loc << "\n";
            cerr << "Protein: " << prot << "\n";
          }
          if (loc != string::npos) {
            midx[specIdx].resize(idxPep.size());
            for (peakIdx = 0; peakIdx < idxPep.size(); peakIdx++) {
              midx[specIdx][peakIdx].set(peakIdx, loc + idxPep[peakIdx]);
              matchScore += specSet[specIdx][peakIdx][1];
            }
          }
        }
        eitp[specIdx][0] = (unsigned short) round(10000.0 * matchScore
            / specScore);
        if (specIdx == debugIdx)
          cerr << " -->> eitp[specIdx][0] = " << eitp[specIdx][0]
              << ", matchScore= " << matchScore << ", specScore= " << specScore
              << "\n";
        if (eitp[specIdx][0] < minPercEI or eitp[specIdx][1] < minPercTP) {
          eitp[specIdx][0] = 0;
          eitp[specIdx][1] = 0;
          specSet[specIdx].resize(0);
          midx[specIdx].resize(0);
          proteinMatches[specIdx][0] = -1;
        }
        if (specIdx == debugIdx)
          cerr << " -->> eitp[specIdx][0] = " << eitp[specIdx][0]
              << ", matchScore= " << matchScore << "\n";
      }
    }

    midx.SaveSpecSet_pklbin(params.getValue("OUTPUT_MIDX"));
    if (params.paramPresent("OUTPUT_SPECS"))
      specSet.SaveSpecSet_pklbin(params.getValue("OUTPUT_SPECS"));
    if (params.paramPresent("OUTPUT_MATCHED_PROTS"))
      Save_binArray(params.getValue("OUTPUT_MATCHED_PROTS"), proteinMatches);
    if (params.paramPresent("OUTPUT_EITP"))
      Save_binArray(params.getValue("OUTPUT_EITP"), eitp);
  }

  if (params.paramPresent("OUTPUT_MATCHED_PROTS")
      and not params.paramPresent("OUTPUT_MIDX"))
    Save_binListArray<int, float, list<pair<int, float> > , list<pair<int,
        float> >::iterator> (params.getValue("OUTPUT_MATCHED_PROTS"),
                             simpleProteinHits);

  free(tagText);
  fclose(pepsFID);

  if (correctPepsAsSpecs.size() > 0) {
    Save_binArray("tags_countCorrect.bla", numCorrectTags);
    FILE *tagsFID = (FILE *) fopen("tags_correct.csv", "wb");
    if (tagsFID) {
      string buffer;
      for (specIdx = 0; specIdx < correctTags.size(); specIdx++)
        for (list<Tag>::iterator tagIter = correctTags[specIdx].begin(); tagIter
            != correctTags[specIdx].end(); tagIter++)
          fprintf(tagsFID,
                  "%.1f\t%s\t%.1f\t%.1f\n",
                  tagIter->flankingPrefix,
                  tagIter->asString(buffer).c_str(),
                  tagIter->flankingSuffix,
                  tagIter->score);
      fclose(tagsFID);
    }
  }

  if (params.paramPresent("OUTPUT_FASTA")) {
    pepsFID = (FILE *) fopen(params.getValue("OUTPUT_FASTA"), "wb");
    for (unsigned int pIdx = 0; pIdx < matchedProts.size(); pIdx++)
      if (matchedProts[pIdx]) {
        fprintf(pepsFID,
                ">%s %s\n%s\n",
                db.getID(pIdx),
                db.getDesc(pIdx),
                db[pIdx]);
      }
    fclose(pepsFID);
  }

  return 0;
}

