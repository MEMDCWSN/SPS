// Header Includes
#include "ExecTagSearch.h"

// Module Includes
#include "Logger.h"
//#include "FileUtils.h"

// SpecNets Includes
#include "tuple.h"

using namespace std;
using namespace specnets;

static bool DEBUG_TAGSEARCH_CORRECT = false;

namespace specnets
{
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
    void set(int   inProtIdx,
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

  ExecTagSearch::ExecTagSearch(void) :
    m_spectra(0x0), m_db(0x0), m_specsToSearch(0x0), ownInput(true),
        ownOutput(true)
  {
    m_name = "ExecTagSearch";
    m_type = "ExecTagSearch";
  }

  // -------------------------------------------------------------------------

  ExecTagSearch::ExecTagSearch(const ParameterList & inputParams) :
    ExecBase(inputParams), m_spectra(0x0), m_db(0x0), m_specsToSearch(0x0),
        ownInput(true), m_psmSet(0x0), ownOutput(true)
  {
    m_name = "ExecTagSearch";
    m_type = "ExecTagSearch";
  }

  // -------------------------------------------------------------------------

  ExecTagSearch::ExecTagSearch(const ParameterList & inputParams,
                               SpecSet * spectra,
                               DB_fasta * db,
                               vector<unsigned int> * specsToSearch,
                               PeptideSpectrumMatchSet * psmSet) :
    ExecBase(inputParams), m_spectra(spectra), m_db(db),
        m_specsToSearch(specsToSearch), ownInput(false), m_psmSet(psmSet), ownOutput(false)
  {
    m_name = "ExecTagSearch";
    m_type = "ExecTagSearch";
  }

  // -------------------------------------------------------------------------

  ExecTagSearch::~ExecTagSearch(void)
  {
    if (ownInput) {
      if (m_spectra)
        delete m_spectra;
      if (m_db)
        delete m_db;
      if (m_specsToSearch)
        delete m_specsToSearch;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecTagSearch::clone(const ParameterList & inputParams) const
  {
    return new ExecTagSearch(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecTagSearch::invoke(void)
  {
    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");
    float pmTol = (float) m_params.getValueDouble("TOLERANCE_PM");

    int numJumps =
            m_params.exists("DOUBLE_AA_JUMPS") ? (int) m_params.getValueInt("DOUBLE_AA_JUMPS") : 1;
    DEBUG_VAR(numJumps);
    int tagLen =
        m_params.exists("TAG_LEN") ? (int) m_params.getValueInt("TAG_LEN") : 6;
    DEBUG_VAR(tagLen);
    short minMatchFlanking =
            m_params.exists("MATCH_TAG_FLANKING_MASSES") ? (short) m_params.getValueInt("MATCH_TAG_FLANKING_MASSES") : 1;
    bool topScoringOnly =
            m_params.exists("TAG_MATCH_TOP_SCORING_ONLY") ? m_params.getValueBool("TAG_MATCH_TOP_SCORING_ONLY") : true;
    DEBUG_VAR(topScoringOnly);
    int maxNumTags =
            m_params.exists("MAX_NUM_TAGS") ? m_params.getValueInt("MAX_NUM_TAGS") : 0;
    DEBUG_VAR(maxNumTags);
    int specType =
            m_params.exists("SPEC_TYPE_MSMS") ? ((int) m_params.getValueInt("SPEC_TYPE_MSMS") ? 1 : 0) : 0;
    bool maxParsimony =
            m_params.exists("MAX_PARSIMONY") ? m_params.getValueBool("MAX_PARSIMONY", false) : false;
    DEBUG_VAR(maxParsimony);
    unsigned short minPercEI =
        m_params.exists("MIN_PERC_EXPINT") ? (unsigned short) (10000 * m_params.getValueDouble("MIN_PERC_EXPINT")) : 0;
    unsigned short minPercTP =
        m_params.exists("MIN_PERC_TP") ? (unsigned short) (10000 * m_params.getValueDouble("MIN_PERC_TP")) : 0;

    float peakSkipPenalty = (float)m_params.getValueDouble("TAG_PEAK_SKIP_PENALTY");
    DEBUG_VAR(peakSkipPenalty);
    unsigned int specIdx, peakIdx, searchIdx; searchIdx;

    unsigned int scanFirst = m_params.getValueInt("SCAN_FIRST", 0);
    unsigned int scanLast = m_params.getValueInt("SCAN_LAST", 
					(*m_spectra)[m_spectra->size()-1].scan);
    DEBUG_VAR(scanFirst);
    DEBUG_VAR(scanLast);
    
    unsigned int debugIdx = m_params.getValueInt("DEBUG_TAGSEARCH_SINGLE", m_spectra->size());
    DEBUG_VAR(debugIdx);

    if (!m_spectra or !m_db)
      return false;

    if (!m_spectra or m_spectra->size() == 0) {
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

    if (!m_db or m_db->size() == 0) {
      ERROR_MSG("ERROR: empty database");
      return false;
    }

    AAJumps jumps(1);

    vector<unsigned int> specsToSearch;
    if (m_specsToSearch == (vector<unsigned int> *) 0
        or m_specsToSearch->empty()) {
      if (m_specsToSearch == (vector<unsigned int> *) 0)
        m_specsToSearch = &specsToSearch;
      m_specsToSearch->resize(m_spectra->size());
      for (searchIdx = 0; searchIdx < m_spectra->size(); searchIdx++)
        (*m_specsToSearch)[searchIdx] = searchIdx;
    }

    // Variables for evaluating tagging efficiency
    if (m_params.exists("CORRECT_PEPTIDE_SPECS")) {
      DEBUG_VAR(m_params.exists("CORRECT_PEPTIDE_SPECS"));
      m_correctPepsAsSpecs.loadPklBin(m_params.getValue("CORRECT_PEPTIDE_SPECS").c_str());
      if (m_correctPepsAsSpecs.size() != m_spectra->size())
        m_correctPepsAsSpecs.resize(0);
      else {
        m_correctTags.resize(m_correctPepsAsSpecs.size());
        m_numCorrectTags.resize(m_correctPepsAsSpecs.size());
      }
    } else if (m_params.exists("CORRECT_PEPTIDE_PSMS")) {
      DEBUG_VAR(m_params.exists("CORRECT_PEPTIDE_PSMS"));
      PeptideSpectrumMatchSet  psmCorrect;
      if (!psmCorrect.loadFromFile(m_params.getValue("CORRECT_PEPTIDE_PSMS").c_str())) {
        ERROR_MSG("Loading Correct PSM file [" << m_params.getValue("CORRECT_PEPTIDE_PSMS") << "]");
        m_correctPepsAsSpecs.resize(0);
        m_correctTags.resize(0);
        m_numCorrectTags.resize(0);
      } else {
        for (int iPSM = 0; iPSM < psmCorrect.size(); iPSM++ ) {
          psmPtr psmc = psmCorrect[iPSM];
          vector<float> masses;
          if (specIdx == debugIdx) DEBUG_VAR(psmc->m_annotation);
          jumps.getPRMMasses(psmc->m_annotation,
                             masses,
                             0.0,  // Offset for masses
                             0x0,  // token vector pointer (we don't care)
                             true);  // Add zero mass (why not?)
          Spectrum specCorrect;
          for (int iMass = 0; iMass < masses.size(); iMass++) {
            if (specIdx == debugIdx) DEBUG_VAR(psmc->m_annotation);
            specCorrect.insertPeak(masses[iMass], 1.0, 1.0);
            specCorrect.scan = psmc->m_scanNum; 
          }
          m_correctPepsAsSpecs.push_back(specCorrect);
        } // for (int iPSM = 0; iPSM < psmSetContigTag.size(); iPSM++ )

        m_correctTags.resize(m_correctPepsAsSpecs.size());
        m_numCorrectTags.resize(m_correctPepsAsSpecs.size());
      }
    }

    //    DB_fasta db; if(db.Load(dbFN)<=0) return -1;
    //DB_fasta dbOri; // Same database as m_db but without ignoring I/L, Q/K distinctions
    //    dbOri.Load(m_params.getValue("INPUT_FASTA").c_str());
    m_dbOri = *m_db;
    DEBUG_VAR(m_db->size());
    m_db->replaceAA('L', 'I');
    m_db->replaceAA('K', 'Q');

    BufferedLineReader prevPeps;
    if (m_params.exists("INPUT_PEPTIDES"))
      prevPeps.Load(m_params.getValue("INPUT_PEPTIDES").c_str());

    DB_index index(*m_db, 1024, 1024, tagLen);

    FILE *pepsFID = (FILE *) 0;
    if (m_params.exists("OUTPUT_PEPTIDES")) {
      pepsFID = fopen(m_params.getValue("OUTPUT_PEPTIDES").c_str(), "wb");
      if (pepsFID)
        fprintf(pepsFID,
                "SpecIdx\tPeptide\tTag score\tPerc peptide score\tProtIdx\tProtein\n");
    }

    list<Tag> tags;
    int matchedTags = 0, matchedRevTags = 0;
    char *tagText = (char *) malloc(tagLen + 1);
    tagText[tagLen] = (char) 0;
    list<sps::tuple<int, float, string> > matches; // List of matched (protein index, peptide start, peptide string)
    list<pair<string, bool> > tagList; // List of which tages matched and their orientation
    char buffer2k[2048]; // Buffer for manipulation of peptide strings

    vector<float> pepMasses; // Peptide masses derived from peptide sequence
    Spectrum pepSpec; // Peptide spectrum derived from peptide sequence

    for (searchIdx = 0; searchIdx < m_specsToSearch->size(); searchIdx++) {
      specIdx = (*m_specsToSearch)[searchIdx];
      int scanNum = (*m_spectra)[specIdx].scan;
      DEBUG_VAR(scanNum);
      if (scanNum < scanFirst || scanNum > scanLast) {
        continue;
      }
      
      DEBUG_MSG("Processing spectrum [" << specIdx <<
                "] scan [" << scanNum << "]");

      if (debugIdx != m_spectra->size() && scanNum != debugIdx) {
        continue;
      }

      float specScore = 0, matchTagScore = 0;
      for (peakIdx = 0; peakIdx < (*m_spectra)[specIdx].size(); peakIdx++) {
        specScore += (*m_spectra)[specIdx][peakIdx][1];
      }
      (*m_spectra)[specIdx].addZPMpeaks(peakTol, ((float) specType) * AAJumps::massHion, true);

      if (specIdx == debugIdx) {
        if (prevPeps.size() > specIdx)
          DEBUG_MSG("Spectrum " << debugIdx << " has peptide ID " << prevPeps.getline(specIdx) << "|");
      }

      if ( not (*m_spectra)[specIdx].psmList.empty() ) {  // Matching previously assigned peptide
        DEBUG_MSG("Matching previously assigned peptide...");
        string curPeptide = (*m_spectra)[specIdx].psmList.front()->m_annotation;
        if ( curPeptide.size() < 2048) {
          getPepSeq( curPeptide.c_str() , buffer2k);
          m_dbOri.find(buffer2k, matches); // Search unmodified peptide string
          for (list<sps::tuple<int, float, string> >::iterator iter = matches.begin();
               iter != matches.end();
               iter++) {
            iter->m2 = curPeptide;
            if (specIdx == debugIdx)
              DEBUG_MSG(" -- match (" << iter->m0 << "," << iter->m2 << ")");
          }
          if (matches.empty()) {
            matches.push_back(sps::make_tuple<int, float, string>(-1, -1.0, curPeptide));
            tagList.push_back(make_pair<string,int>("",0));
            if (specIdx == debugIdx)
              DEBUG_MSG(" -- No protein match: (" << matches.front().m0 << "," << matches.front().m2 << ")");
          }
          matchTagScore = MatchSpecToPeptide((*m_spectra)[specIdx],
                                             (*m_spectra)[specIdx].psmList.front()->m_annotation.c_str(),
                                             peakTol,
                                             0,
                                             false);
          if (specIdx == debugIdx)
            DEBUG_MSG(" -- matchTagScore = " << matchTagScore);
        }
      }
      else
      {
        DEBUG_MSG("Tag based matching...");
        // Tag-based matches
        tags.clear();

#if 1
        Spectrum & specCurr = (*m_spectra)[specIdx];
//        specCurr.output(cerr);
        // Clean up any pre-existing endpoints so we don't add too many
        list<int> peaksToRemove;
        for (int iPeak = 0; iPeak < specCurr.size(); iPeak++) {
          if (specCurr[iPeak][0] < 57.0 || specCurr[iPeak][0] > specCurr.parentMass - 57.0) {
            peaksToRemove.push_back(iPeak);
          }
        }
        if (peaksToRemove.size() != 0) {
          specCurr.removePeaks(peaksToRemove);
        }
        bool specTypeMSMS = m_params.getValueBool("SPEC_TYPE_MSMS", false);
        float ionOffset = specTypeMSMS ? AAJumps::massHion : 0;
        specCurr.addZPMpeaks(peakTol, ionOffset, false);
//        specCurr.output(cerr);
#endif

        if(maxNumTags>=0)
//        	ExtractTags((*m_spectra)[specIdx],
        	ExtractTagsAllJumps((*m_spectra)[specIdx],
        			tags,
        			peakTol,
        			tagLen,
        			numJumps,
        			(unsigned int)maxNumTags,
                            peakSkipPenalty);
        DEBUG_MSG(tags.size() << " tags...");
        list<Tag>::iterator itr = tags.begin();
        list<Tag>::iterator itr_end = tags.end();
        for (; itr != itr_end; itr++) {
          string seq;
          for (int c = 0; c < itr->sequence.size(); c++) {
            seq += jumps.aaLetters[itr->sequence[c]];
          }
          if (specIdx == debugIdx)
            DEBUG_MSG("[" << itr->flankingPrefix << "]" << seq << "[" << itr->flankingSuffix << "]" << "\t" << itr->score);
        }

        //-----------------------------------------------
        // MATCH THE TAGS TO "CORRECT" SPECTRUM AS A TEST
        //-----------------------------------------------
        DEBUG_VAR(m_correctPepsAsSpecs.size());
        if (m_correctPepsAsSpecs.size() > 0) {
          SpecSet curSpec(1);
          curSpec[0] = m_correctPepsAsSpecs[specIdx];
          vector<list<unsigned int> > matchedTagsIdx;
          vector<list<float> > matchedTagsScores;
          vector<list<unsigned int> > matchedTagsPos;

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
          m_numCorrectTags[specIdx] = matchedTagsIdx.size();

          list<Tag>::iterator tagIter = tags.begin();
          list<unsigned int>::iterator matchIter = matchedTagsIdx[0].begin();
          for (unsigned int tagCount = 0; tagCount < tags.size(); tagCount++, tagIter++)
            if (*matchIter == tagCount) {
              m_correctTags[specIdx].push_back(*tagIter);
              matchIter++;
            }
        }
#if 0
        if (m_correctPepsAsSpecs.size() > 0) {
          SpecSet curSpec(1);
          curSpec[0] = m_correctPepsAsSpecs[specIdx];
          vector<list<unsigned int> > matchedTagsIdx;
          vector<list<float> > matchedTagsScores;
          vector<list<unsigned int> > matchedTagsPos;

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
          if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(tags.size());
          if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(matchedTagsIdx[0].size());
          if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(matchedTagsScores[0].size());
          if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(matchedTagsPos[0].size());
          m_numCorrectTags[specIdx] = tags.size();
          if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(m_numCorrectTags[specIdx]);

          list<unsigned int>::iterator matchIter = matchedTagsIdx[0].begin();
          list<unsigned int>::iterator matchIterEnd = matchedTagsIdx[0].end();
          list<unsigned int>::iterator posIter = matchedTagsPos[0].begin();
          list<float>::iterator scoreIter = matchedTagsScores[0].begin();
          for (; matchIter != matchIterEnd; matchIter++, posIter++, scoreIter++) {

            if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(*matchIter );
            if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(*posIter);
            if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(*scoreIter);

            list<Tag>::iterator tagIter = tags.begin();
            for (int iTag = 0; iTag != !*matchIter; iTag++, tagIter++);

            psmPtr p(new PeptideSpectrumMatch);
            p->m_scanNum = curSpec[0].scan;
            if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(p->m_scanNum);

              string seq;
              for (int c = 0; c < tagIter->sequence.size(); c++) {
                seq += jumps.aaLetters[tagIter->sequence[c]];
              }
              p->m_origAnnotation = seq;
              if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(p->m_origAnnotation);

              char anno[256];
              sprintf(anno,
                    "[%.1f]%s[%.1f]",
                    tagIter->flankingPrefix,
                    seq.c_str(),
                    tagIter->flankingSuffix,
                    tagIter->score);

              p->m_annotation = string(anno);
              if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(p->m_annotation);
              p->m_startMass = *posIter;
              if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(p->m_startMass);
              p->m_score = tagIter->score;
              if (DEBUG_TAGSEARCH_CORRECT) DEBUG_VAR(p->m_score);
              //p->m_matchOrientation = ???;
              m_psmCorrect.push_back(p);
          }
        }
        //-----------------------------------------------
#endif
        if (specIdx == debugIdx)
          DEBUG_MSG(" -1- (*m_spectra)[" << debugIdx << "]size() = " << (*m_spectra)[debugIdx].size());

        matches.clear();
        tagList.clear();
        matchTagScore = 0;
        for (list<Tag>::iterator tagIter = tags.begin(); tagIter != tags.end(); tagIter++) {

          if ((maxParsimony || topScoringOnly) && tagIter->score < matchTagScore - 0.0001)
            break;

	  // Match tags in forward orientation
          for (int cIdx = 0; cIdx < tagLen; cIdx++)
            tagText[cIdx] = jumps.aaLetters[tagIter->sequence[cIdx]];

          if (specIdx == debugIdx) DEBUG_MSG(tagText << "  " << tagIter->score);

          int preSize = matches.size();
          if (specIdx == debugIdx) DEBUG_VAR(preSize);
          if (index.findAll(tagText,
                         (*m_db),
                         matches,
                         minMatchFlanking,
                         tagIter->flankingPrefix,
                         peakTol,
                         tagIter->flankingSuffix,
                         pmTol)) {
            matchedTags++;
            matchTagScore = tagIter->score;
          }
          int postSize = matches.size();
          if (specIdx == debugIdx) DEBUG_VAR(postSize);
          for (int i = 0; i < postSize - preSize; i++) {
            tagList.push_back(make_pair<string,int>(tagText,0));
          }

	  // Match tags in reverse orientation
          for (int cIdx = 0; cIdx < tagLen; cIdx++)
            tagText[cIdx] = jumps.aaLetters[tagIter->sequence[tagLen - 1 - cIdx]];

          if (specIdx == debugIdx) DEBUG_MSG(tagText << "  " << tagIter->score);

          preSize = matches.size();
          if (specIdx == debugIdx) DEBUG_VAR(preSize);
          if (index.findAll(tagText,
                         (*m_db),
                         matches,
                         minMatchFlanking,
                         tagIter->flankingSuffix + 18,
                         pmTol,
                         tagIter->flankingPrefix - 18,
                         peakTol)) {
            matchedRevTags++;
            matchTagScore = tagIter->score;
          }
          postSize = matches.size();
          if (specIdx == debugIdx) DEBUG_VAR(postSize);
          for (int i = 0; i < postSize - preSize; i++) {
            tagList.push_back(make_pair<string,int>(tagText,1));
          }
        }

        if (specIdx == debugIdx)
          DEBUG_MSG(" -2- (*m_spectra)[" << debugIdx << "]size() = " << (*m_spectra)[debugIdx].size());

        // Convert Q/I-only matches to Q/K,I/L matches (m_dbOri)
        for (list<sps::tuple<int,float,string> >::iterator match = matches.begin();
            match != matches.end();
            match++) {
          string prot(m_db->getSequence(match->m0));
          string::size_type loc = prot.find(match->m2);
          if (loc != string::npos) {
            char *protOri = m_dbOri.getSequence(match->m0);
            strncpy(buffer2k, &protOri[loc], match->m2.size());
            buffer2k[match->m2.size()] = (char) 0;
            match->m2 = string(buffer2k);
          }
        }
      }

      if (specIdx == debugIdx)
        DEBUG_MSG(" -3- (*m_spectra)[" << debugIdx << "]size() = " << (*m_spectra)[debugIdx].size());

      DEBUG_VAR(matches.size());
      DEBUG_VAR(tagList.size());
      list<pair<string, bool> >::iterator tliter = tagList.begin();
      for (list<sps::tuple<int, float, string> >::iterator match = matches.begin();
           match != matches.end();
           match++,tliter++) {

        //DEBUG_VAR(match->m0);
        if (match->m0 >= 0) {
          psmPtr p(new PeptideSpectrumMatch);
          //p->m_spectrumFile = ;
          p->m_scanNum = (*m_spectra)[specIdx].scan;
          //DEBUG_VAR(p->m_scanNum);
          // Removed description for now because of an apparent parsing problem
          //p->m_protein = string(m_dbOri.getID(match->m0)) + " " + string(m_dbOri.getDesc(match->m0));
          p->m_protein = string(m_dbOri.getID(match->m0));
          //DEBUG_VAR(p->m_protein);
          p->m_annotation = match->m2;
          //DEBUG_VAR(p->m_annotation);
          p->m_origAnnotation = tliter->first;
          p->m_dbIndex = match->m0;
          p->m_startMass = match->m1;
          p->m_score = matchTagScore;
          p->m_matchOrientation = tliter->second;
          m_psmSet->push_back(p);
          //(*m_spectra)[specIdx].psmList.push_back(p);
        }
        if (pepsFID)
          fprintf(pepsFID,
                  "%d\t%s\t%.1f\t%.1f\t%d",
                  specIdx + 1,
                  match->m2.c_str(),
                  matchTagScore,
                  100 * matchTagScore / specScore,
                  match->m0);
        if (pepsFID)
          fprintf(pepsFID,
                  "\t%s %s\n",
                  m_db->getID(match->m0),
                  m_db->getDesc(match->m0));
      }
      DEBUG_MSG("done");
    }
    DEBUG_VAR(m_psmSet->size());

    if (maxParsimony) {
      m_spectra->maximumParsimony();
    }

    free(tagText);
    if (pepsFID)
      fclose(pepsFID);

    if (m_specsToSearch == &specsToSearch)
      m_specsToSearch = (vector<unsigned int> *) 0;

    // Put the DB back the way it was
    *m_db = m_dbOri;

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecTagSearch::loadInputData(void)
  {
    if (ownInput) {
      if (!m_spectra)
        m_spectra = new SpecSet;
      if (!m_db)
        m_db = new DB_fasta;
      if (!m_specsToSearch)
        m_specsToSearch = new vector<unsigned int> ;
    }
    if (ownOutput) {
      m_psmSet = new PeptideSpectrumMatchSet;
    }
    m_spectra->resize(0);
    m_specsToSearch->resize(0);

    if (m_params.exists("AMINO_ACID_MASSES")) {
      AAJumps tmpJumps(-1);
      tmpJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true); // Set global defaults for amino acid masses
    }

    if (!m_params.exists("INPUT_SPECS_PKLBIN")) {
      ERROR_MSG("Parameters are incomplete. INPUT_SPECS_PKLBIN is missing.");
      return false;
    }
    else if (m_spectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str())
        <= 0 or m_spectra->size() == 0) {
      ERROR_MSG("Error reading input spectra from " <<  m_params.getValue("INPUT_SPECS_PKLBIN"));
      return false;
    }

    if (!m_params.exists("INPUT_FASTA")) {
      ERROR_MSG("Parameters are incomplete. INPUT_FASTA is missing.");
      return false;
    }
    else if (m_db->Load(m_params.getValue("INPUT_FASTA").c_str()) <= 0) {
      ERROR_MSG("Error reading database sequences from " << m_params.getValue("INPUT_FASTA"));
      return false;
    }

    if (m_params.exists("SUBSET_IDX")) {
      if (Load_binArray(m_params.getValue("SUBSET_IDX").c_str(),
                        *m_specsToSearch) <= 0) {
        ERROR_MSG("Error reading spectrum indices from " << m_params.getValue("SUBSET_IDX"));
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecTagSearch::saveOutputData(void)
  {
    if (m_correctPepsAsSpecs.size() > 0) {
      AAJumps jumps(1);

      //Save_binArray("tags_countCorrect.bla", m_numCorrectTags);

      FILE *countsFID = (FILE *) fopen("tags_correct.txt", "w");
      if (countsFID) {
        for (int i = 0; i < m_numCorrectTags.size(); i++) {
          fprintf(countsFID, "%u\n", m_numCorrectTags[i]);
        }
        fclose(countsFID);
      }

      FILE *tagsFID = (FILE *) fopen("tags_correct.csv", "w");
      if (tagsFID) {
        for (int specIdx = 0; specIdx < m_correctTags.size(); specIdx++) {
          for (list<Tag>::iterator tagIter = m_correctTags[specIdx].begin();
               tagIter != m_correctTags[specIdx].end();
               tagIter++) {
            string seq;
            for (int c = 0; c < tagIter->sequence.size(); c++) {
              seq += jumps.aaLetters[tagIter->sequence[c]];
            }
            fprintf(tagsFID,
                    "%.1f\t%s\t%.1f\t%.1f\n",
                    tagIter->flankingPrefix,
                    seq.c_str(),
                    tagIter->flankingSuffix,
                    tagIter->score);
          }
        }
        fclose(tagsFID);
      }

      m_psmCorrect.saveToFile("tag_correct.psm");

    }

    if (m_spectra and m_params.exists("OUTPUT_SPECS")) {
      m_spectra->savePklBin(m_params.getValue("OUTPUT_SPECS").c_str());
    }

    if(m_params.exists("OUTPUT_PSM")) {
      m_psmSet->saveToFile(m_params.getValue("OUTPUT_PSM").c_str(), true);
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecTagSearch::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecTagSearch::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecTagSearch::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecTagSearch::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecTagSearch::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");

    m_isValid = true;
    return true;
  }

} // namespace specnets

