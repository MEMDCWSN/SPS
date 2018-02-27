#include "ExecSpecProtAlign.h"

#include "alignment_modmut.h"
#include "AlignmentPenaltyBased.h"
#include "ExecTagSearch.h"
#include "FdrPeptide.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"

// SpecNets Includes
#include "tags.h"
#include <limits.h>
#include <time.h>

static bool DEBUG_SPECPROTALIGN = false;
static bool DEBUG_SPECPROTALIGN_SPECTRA = false;
static bool DEBUG_SPECPROTALIGN_TIME = false;
static bool DEBUG_SPECPROTALIGN_SPRINKLE = false;
static bool DEBUG_SPECPROTALIGN_ANNO = false;
static bool DEBUG_SPECPROTALIGN_MERGE = false;
static bool DEBUG_SPECPROTALIGN_SPECPROB = false;
static bool DEBUG_SPECPROTALIGN_SPLIT = false;
static bool DEBUG_SPECPROTALIGN_RANGE = false;
static bool DEBUG_SPECPROTALIGN_COMPUTE_GAP_ANNOS = false;
static bool DEBUG_SPECPROTALIGN_TAG_SEEDING = false;
static bool DEBUG_SPECPROTALIGN_TAG_SEEDING2 = false;

const float AVERAGE_AA_MASS = 110.0;
const float MAX_SPEC_PROB_PEAK = 100.0;

using namespace std;
using namespace specnets;

namespace specnets
{
  class LessThanPredicate
  {
  public:
    LessThanPredicate(float value) :
      m_value(value)
    {
    }
    bool operator()(const psmPtr & value)
    {
      return value->m_score < m_value;
    }
  private:
    float m_value;
  };

  // -------------------------------------------------------------------------
  bool PsmDbIndexAnnoSort(psmPtr left, psmPtr right)
  {
    if (left->m_dbIndex == right->m_dbIndex) {
      return left->m_annotation < right->m_annotation;
    } else {
      return left->m_dbIndex < right->m_dbIndex;
    }
  }

  // -------------------------------------------------------------------------
  bool PsmDbIndexAnnoUnique(psmPtr left, psmPtr right)
  {
    return (left->m_dbIndex == right->m_dbIndex) && (left->m_annotation == right->m_annotation);
  }

  // -------------------------------------------------------------------------
  void DebugPsm(psmPtr & psm)
  {
    DEBUG_MSG(psm->m_scanNum << "  " << 
              psm->m_matchOrientation << "  " << 
              psm->m_annotation << "  " << 
              psm->m_score << "  " << 
              psm->m_pValue << "  " << 
              psm->m_dbIndex << "  " << 
              psm->m_protein);
    return;
  } 

  // -------------------------------------------------------------------------
  float averageTopPeaks(Spectrum & spec, int topN)
  {
    float average = 0.0;
    vector<float> vecPeaks;
    for (int i = 0; i < spec.size(); i++) {
      vecPeaks.push_back(spec[i][1]);
    }
    sort(vecPeaks.begin(), vecPeaks.end(), std::greater<int>());

    float peaksCounted = 0.0;
    for (int i = 0; i < vecPeaks.size() && i < topN; i++) {
      average += vecPeaks[i];
      peaksCounted += 1.0;
    }
    average /= peaksCounted;
    return average;
  }

  //-----------------------------------------------------------------------------
  void ExecSpecProtAlign::loadScanSpecificPenalties(
                                  string & filename,
                                  map<int, map<int, float> > & penalties)
  {
    vector<string> header;
    vector<vector<string> > lines;
    vector<string> requiredHeader;
    vector<int> requiredHeaderIndex;
    DelimitedTextReader::loadDelimitedFile(filename.c_str(),
                              "\t",
                              "",
                              header,
                              lines,
                              requiredHeader,
                              requiredHeaderIndex);
    //DEBUG_VAR(lines.size());
    for (int i = 0; i < lines.size(); i++) {
      //DEBUG_VAR(lines[i][0]);
      string scanString = lines[i][0];
      string massString = lines[i][1];
      string scoreString = lines[i][2];
      int scan = atoi(scanString.c_str());
      int mass = atoi(massString.c_str());
      float score = atof(scoreString.c_str());
      //DEBUG_MSG(scan << "  " << mass  << "  " << score);
      penalties[scan][mass] = score;
    }
    return;
  }

  //-----------------------------------------------------------------------------
  float roundMass(float mass)
  {
    int sign = mass < 0 ? -1 : 1;
    return (int)(fabs(mass) * 0.9995 + 0.5) * sign;
    //return abs((float)((int)(mass / m_resolution)) * m_resolution);
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlign::ExecSpecProtAlign(void) :
    m_inputSpectra(0x0), m_prmSpectra(0x0), m_db(0x0), 
    m_penaltyMatrixBlosum(0x0), m_penaltyMatrixMods(0x0),
    m_scanSpecificPenalties(0x0),
    m_contigAbinfo(0x0), m_filterPsmSet(0x0), 
    ownInput(true), 
    m_matchedSpectraAll(0x0),
    ownOutput(true)
  {
    ownSpectra = false;
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       SpecSet * prmSpectra,
                                       DB_fasta * db,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods,
                                       map<int, map<int, float> > * scanSpecificPenalties,
                                       abinfo_t * contigAbinfo,
                                       PeptideSpectrumMatchSet * filterPsmSet) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), m_prmSpectra(prmSpectra), m_db(db),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum),
        m_penaltyMatrixMods(penaltyMatrixMods), 
        m_scanSpecificPenalties(scanSpecificPenalties),
        m_contigAbinfo(contigAbinfo), m_filterPsmSet(filterPsmSet), 
        ownInput(false),
        m_matchedSpectraAll(0x0), ownOutput(true)
  {
    ownSpectra = false;
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputSpectra(0x0), m_prmSpectra(0x0), m_db(0x0), 
       m_penaltyMatrixBlosum(0x0), m_penaltyMatrixMods(0x0), 
       m_scanSpecificPenalties(0x0),
       m_contigAbinfo(0x0), m_filterPsmSet(0x0), 
       ownInput(true),
       m_matchedSpectraAll(0x0), ownOutput(true)
  {
    ownSpectra = true;
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       SpecSet * prmSpectra,
                                       DB_fasta * db,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods,
                                       map<int, map<int, float> > * scanSpecificPenalties,
                                       abinfo_t * contigAbinfo,
                                       PeptideSpectrumMatchSet * filterPsmSet,
                                       SpecSet * matchedSpectraAll) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), m_prmSpectra(prmSpectra), m_db(db),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum),
        m_penaltyMatrixMods(penaltyMatrixMods), 
        m_scanSpecificPenalties(scanSpecificPenalties),
        m_contigAbinfo(contigAbinfo), m_filterPsmSet(filterPsmSet), 
        ownInput(false),
        m_matchedSpectraAll(matchedSpectraAll), ownOutput(false)
  {
    ownSpectra = false;
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlign::~ExecSpecProtAlign(void)
  {
    DEBUG_TRACE;

    if (ownInput && ownSpectra)
    {
      if (m_inputSpectra) {
        delete m_inputSpectra;
      }
    }
    if (ownInput)
    {
      if (m_prmSpectra) {
        delete m_prmSpectra;
      }
      if (m_db) {
        delete m_db;
      }
      if (m_penaltyMatrixBlosum) {
        delete m_penaltyMatrixBlosum;
      }
      if (m_penaltyMatrixMods) {
        delete m_penaltyMatrixMods;
      }
      if (m_scanSpecificPenalties) {
        delete m_scanSpecificPenalties;
      }
    }
    if (ownOutput)
    {
      if (m_matchedSpectraAll) {
        delete m_matchedSpectraAll;
      }
    }

  }

  // -------------------------------------------------------------------------
  ExecBase * ExecSpecProtAlign::clone(const ParameterList & inputParams) const
  {
    return new ExecSpecProtAlign(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecSpecProtAlign::invoke(void)
  {
    if (!m_inputSpectra || m_inputSpectra->size() == 0) {
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

    // More readable code if we use a reference instead of pointer
    SpecSet & inputSpectra = *m_inputSpectra;

    if (!m_db or m_db->size() == 0) {
      ERROR_MSG("ERROR: empty database");
      return false;
    }

    if (!m_matchedSpectraAll) {
      m_matchedSpectraAll = new SpecSet;
      ownOutput = true;
    }

    DEBUG_SPECPROTALIGN = m_params.getValueBool("DEBUG_SPECPROTALIGN");
    DEBUG_VAR(DEBUG_SPECPROTALIGN);
    DEBUG_SPECPROTALIGN_SPECTRA = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPECTRA");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPECTRA);
    DEBUG_SPECPROTALIGN_TIME = m_params.getValueBool("DEBUG_SPECPROTALIGN_TIME");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_TIME);
    DEBUG_SPECPROTALIGN_SPRINKLE = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPRINKLE");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPRINKLE);
    DEBUG_SPECPROTALIGN_ANNO = m_params.getValueBool("DEBUG_SPECPROTALIGN_ANNO");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_ANNO);
    DEBUG_SPECPROTALIGN_MERGE = m_params.getValueBool("DEBUG_SPECPROTALIGN_MERGE");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_MERGE);
    DEBUG_SPECPROTALIGN_RANGE= m_params.getValueBool("DEBUG_SPECPROTALIGN_RANGE");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_RANGE);
    DEBUG_SPECPROTALIGN_COMPUTE_GAP_ANNOS= m_params.getValueBool("DEBUG_SPECPROTALIGN_COMPUTE_GAP_ANNOS");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_COMPUTE_GAP_ANNOS);
    DEBUG_SPECPROTALIGN_TAG_SEEDING= m_params.getValueBool("DEBUG_SPECPROTALIGN_TAG_SEEDING");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_TAG_SEEDING);

    DEBUG_VAR(m_params.getValue("DEBUG_SPECPROTALIGN_SINGLESPECTRUM"));
    DEBUG_VAR(m_params.getValue("DEBUG_SPECPROTALIGN_MODSPECTRUM"));
    DEBUG_VAR(m_params.getValue("DEBUG_SPECPROTALIGN_MODCONTIG"));
    DEBUG_VAR(m_params.getValue("DEBUG_SPECPROTALIGN_SINGLEPROTID"));
    DEBUG_VAR(m_params.getValue("DEBUG_SPECPROTALIGN_SKIPTAGCREATION"));

    // -----------------------------------------------------
    // Get all the parameters for alignment
    // -----------------------------------------------------
    unsigned int startIdx = m_params.getValueInt("IDX_START", 0);
    unsigned int endIdx = m_params.exists("IDX_END") ? min(inputSpectra.size()
        - 1, (unsigned int)m_params.getValueInt("IDX_END"))
        : inputSpectra.size() - 1;

    unsigned int scanFirst = m_params.getValueInt("SCAN_FIRST", 0);
    unsigned int scanLast = m_params.getValueInt("SCAN_LAST", 
					inputSpectra[inputSpectra.size()-1].scan);
    DEBUG_VAR(startIdx);
    DEBUG_VAR(endIdx);
    DEBUG_VAR(scanFirst);
    DEBUG_VAR(scanLast);
    
    float pmTol = (float)m_params.getValueDouble("TOLERANCE_PM", 0.4);
    float peakTol = (float)m_params.getValueDouble("TOLERANCE_PEAK", 1.5);
    bool enforceEndpeaks = m_params.getValueInt("ENFORCE_ENDPEAKS");
    float thresholdPercent = m_params.getValueDouble("ALIGNMENT_SCORE_THRESHOLD", 0.75);

    DEBUG_VAR(pmTol);
    DEBUG_VAR(peakTol);
    DEBUG_VAR(enforceEndpeaks);
    DEBUG_VAR(thresholdPercent);

    float maxModMass = m_params.getValueDouble("MAX_MOD_MASS", 100);
    float minModMass = m_params.getValueDouble("MIN_MOD_MASS", -100);
    int maxNumMods = m_params.getValueInt("MAX_NUM_MODS", 2);
    int minNumMatchPeaks = m_params.getValueInt("MIN_MATCHED_PEAKS_DB", 6);
    bool tagParsimony = m_params.getValueBool("MAX_PARSIMONY", 1);
    float resolution = m_params.getValueDouble("ALIGNMENT_RESOLUTION", 1.0);
    
    DEBUG_VAR(maxModMass);
    DEBUG_VAR(minModMass);
    DEBUG_VAR(maxNumMods);
    DEBUG_VAR(minNumMatchPeaks);
    DEBUG_VAR(resolution);

    bool penaltyAlign = m_params.getValueBool("PENALTY_ALIGNMENT", false);
    // Make sure tag parsimony is false if we are using penalty alignment
    if (penaltyAlign)  {
      tagParsimony = false;
    }
    DEBUG_VAR(tagParsimony);
    float penaltyAlpha = m_params.getValueFloat("PENALTY_ALIGNMENT_ALPHA", 1.0);
    float penaltyBeta = m_params.getValueFloat("PENALTY_ALIGNMENT_BETA", 1000.0);
    bool alignall = m_params.getValueBool("ALIGN_ALL", false);
    int maxDbGapAas = m_params.getValueInt("MAX_ALIGN_DB_GAP_AAS", 8);
    int maxSpectrumGapDaltons = m_params.getValueInt("MAX_ALIGN_SPECTRUM_GAP_DALTONS", 1500);

    DEBUG_VAR(penaltyAlign);
    DEBUG_VAR(penaltyAlpha);
    DEBUG_VAR(penaltyBeta);
    DEBUG_VAR(alignall);
    DEBUG_VAR(maxDbGapAas);
    DEBUG_VAR(maxSpectrumGapDaltons);

    bool roundAnno = m_params.exists("SPECPROTALIGN_ROUND_ANNOTATION_MAX");
    float roundAnnoMax = m_params.getValueFloat("SPECPROTALIGN_ROUND_ANNOTATION_MAX", 1.0);
    DEBUG_VAR(roundAnno);
    DEBUG_VAR(roundAnnoMax);
    
    string clusteredSpecFileName = m_params.getValue("CLUSTERED_SPECTRUM_FILENAME");

    int modIdx, specIdx, protIdx;
    Spectrum tmpSpec;
    tmpSpec.reserve(1024);
    Spectrum cSpec;
    cSpec.reserve(1024);
    Spectrum cSpecRev;
    cSpecRev.reserve(1024);
    AMM_match *bestMatch;

#if 0
    // Create decoy database on the fly (in place)
    //   If we create it on the fly there is no way to see it (or use it later)
    bool usingDecoyDb = m_params.getValueInt("USE_DECOY_DATABASE", 0);
    DEBUG_VAR(usingDecoyDb);
    if (usingDecoyDb) {
      m_db->replaceDecoyShuffled();
    }
#endif

    // -----------------------------------------------------
    // For contig alignment - do tag matching inside here
    // -----------------------------------------------------
    bool skipTagCreation = m_params.getValueBool("DEBUG_SPECPROTALIGN_SKIPTAGCREATION");
    if (!enforceEndpeaks && !skipTagCreation) {
      createTags();

      bool useSeeding = m_params.getValueBool("SPECPROTALIGN_USE_TAG_SEEDING", false);
      DEBUG_VAR(useSeeding);

      // We do not (can not) use seed tags on decoys
      //   usingDecoyDb may not be set (database may exist as decoys before start)
      string protein(m_db->getID(0));
      if (useSeeding && protein.find("XXX") == string::npos) {
        createSeedTags();
      }
      
    } else if (enforceEndpeaks) {
      if (m_params.exists("INPUT_PSM")) {
        int scanFirst = m_params.getValueInt("SCAN_FIRST", -1);
        int scanLast = m_params.getValueInt("SCAN_LAST", -1);
        DEBUG_VAR(scanFirst);
        DEBUG_VAR(scanLast);
        PeptideSpectrumMatchSet psmSet;
        string fileName = m_params.getValue("INPUT_PSM");
        if (!psmSet.loadFromFile(fileName.c_str(), scanFirst, scanLast)) {
          ERROR_MSG("Error reading input PSM file.");
          return false;
        }
        psmSet.addSpectra(m_inputSpectra);
      } 
    }
    
    // -----------------------------------------------------
    if (!enforceEndpeaks) { // This is a proxy (for now) to check for contig alignment only
    // Load amino acid masses
      AAJumps jumps(1);
      if (m_params.exists("AMINO_ACID_MASSES")) {
        string fileName = m_params.getValue("AMINO_ACID_MASSES");
        jumps.loadJumps(fileName.c_str(), true);
      }
      float knownModPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_KNOWN_PENALTY", 0.01);
      float unknownPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_PENALTY", 1.0);
      float unknownMultiplier = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER", 1.5);
      float minModMass = m_params.getValueDouble("MIN_MOD_MASS", -100.0);
      float maxModMass = m_params.getValueDouble("MAX_MOD_MASS", 100.0);

      if (ownInput && m_penaltyMatrixMods != 0) {
        delete m_penaltyMatrixMods;
      }
      m_penaltyMatrixMods = new PenaltyMatrix(jumps,
                                    resolution,
                                    knownModPenalty, 
                                    unknownPenalty,
                                    unknownMultiplier,
                                    minModMass,
                                    maxModMass);
      if (ownInput && m_penaltyMatrixBlosum != 0) {
        delete m_penaltyMatrixBlosum;
      }
      m_penaltyMatrixBlosum = new PenaltyMatrix(jumps,
                                    resolution,
                                    knownModPenalty, 
                                    unknownPenalty,
                                    unknownMultiplier,
                                    minModMass,
                                    maxModMass);
    }

    // -----------------------------------------------------
    // Construct the penalty based alignment class
    //   this will initialize all the scoring vectors
    // -----------------------------------------------------
    DEBUG_MSG("Constructing penalty based alignment object...");
    AlignmentPenaltyBased * apb = 0x0;
    if (penaltyAlign) {

      clock_t startTime = clock();
      apb = new AlignmentPenaltyBased(maxSpectrumGapDaltons,
                              m_penaltyMatrixMods,
                              m_penaltyMatrixBlosum,
                              penaltyAlpha,
                              penaltyBeta,
                              maxModMass,
                              minModMass);

       apb->setDebugFlags(
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_RANGE", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_RANGE2", false),

           m_params.getValueBool("PENALTY_ALIGN_DEBUG_AAS", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_SPECS", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN1", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN2", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN3", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN4", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN_EXACT", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN_NTERM", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN_CTERM", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN_GAP", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN_GAP_ANNO", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN_GAP_CACHE", false),
           m_params.getValueBool("PENALTY_ALIGN_DEBUG_ALIGN_LARGE_GAP", false));

      float elapsed = ((float)(clock() - startTime)) / CLOCKS_PER_SEC;
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsed << " seconds elapsed for APB ctor");
    }

    clock_t startTimeTotal = clock();
    
    DEBUG_VAR(inputSpectra.size());
    for (specIdx = 0; specIdx < inputSpectra.size(); specIdx++) {

      int scanNum = inputSpectra[specIdx].scan;
      //DEBUG_MSG("Index [" << specIdx << "]  Scan ["  << scanNum << "]");
      if (scanNum < scanFirst || scanNum > scanLast) {
        continue;
      }

      if ( m_params.exists("DEBUG_SPECPROTALIGN_SINGLESPECTRUM")) {
	      int scanDebug = m_params.getValueInt("DEBUG_SPECPROTALIGN_SINGLESPECTRUM");
        if (scanDebug < scanNum) {
          break;
        }
        if (scanDebug != scanNum) {
          m_matchedSpectraAll->push_back(inputSpectra[specIdx]);
          m_matchedSpectraAllSprinkled.push_back(inputSpectra[specIdx]);
          (*m_matchedSpectraAll)[m_matchedSpectraAll->size()-1].psmList.clear();
          continue;
        }

        map<int, float>::iterator itrm2 = (*m_scanSpecificPenalties)[scanNum].begin();
        map<int, float>::iterator itrmEnd2 = (*m_scanSpecificPenalties)[scanNum].end();
        for ( ; itrm2 != itrmEnd2; itrm2++) {
          DEBUG_MSG(scanNum << "  " << itrm2->first  << "  " << itrm2->second);
        }
      }

      // Debug parameter for skipping every Nth contig
      if (!enforceEndpeaks && m_params.exists("DEBUG_SPECPROTALIGN_MODCONTIG")) {
        int scanModDebug = m_params.getValueInt("DEBUG_SPECPROTALIGN_MODCONTIG");
        if (scanNum % scanModDebug != 0) {
          m_matchedSpectraAll->push_back(inputSpectra[specIdx]);
          m_matchedSpectraAllSprinkled.push_back(inputSpectra[specIdx]);
          (*m_matchedSpectraAll)[m_matchedSpectraAll->size()-1].psmList.clear();
          continue;
        }
      }

      // Debug parameter for skipping every Nth spectrum
      if (enforceEndpeaks && m_params.exists("DEBUG_SPECPROTALIGN_MODSPECTRUM")) {
        int scanModDebug = m_params.getValueInt("DEBUG_SPECPROTALIGN_MODSPECTRUM");
        if (scanNum % scanModDebug != 0) {
          m_matchedSpectraAll->push_back(inputSpectra[specIdx]);
          m_matchedSpectraAllSprinkled.push_back(inputSpectra[specIdx]);
          (*m_matchedSpectraAll)[m_matchedSpectraAll->size()-1].psmList.clear();
          continue;
        }
      }

      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(inputSpectra[specIdx].parentMass);
      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(inputSpectra[specIdx].size());

      // Apparently it is possible to have a contig that is empty
      if (inputSpectra[specIdx].size() == 0) {
        m_matchedSpectraAll->push_back(inputSpectra[specIdx]);
        (*m_matchedSpectraAll)[m_matchedSpectraAll->size()-1].psmList.clear();
        m_matchedSpectraAllSprinkled.push_back(inputSpectra[specIdx]);
        continue;
      }

      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[specIdx].psmList.size());

      // No need to do anything if there are no PSM (and not aligning to all positions)
      if (!alignall && inputSpectra[specIdx].psmList.size() == 0) {
        m_matchedSpectraAll->push_back(inputSpectra[specIdx]);
        m_matchedSpectraAllSprinkled.push_back(inputSpectra[specIdx]);
        continue;
      }

      cSpec = inputSpectra[specIdx];
      cSpec.psmList.clear();

      float avgPeakIntensity = 0;
      // Reverse will be created here
//      prepareSpectrum(cSpec, cSpecRev, avgPeakIntensity, false, false, false);
      // If prepareSpectrum returns false.. there is something wrong with spectrum
      if (!prepareSpectrum(cSpec, cSpecRev, avgPeakIntensity, true, enforceEndpeaks, false)) {
        continue;
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeakIntensity)

      m_matchedSpectraAllSprinkled.push_back(cSpec);

      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[specIdx].psmList.size())

      // -----------------------------------------------------
      // Coalesce the starting positions on all proteins
      // -----------------------------------------------------
      map<int, set<float> > matchedProtMapForward;
      map<int, set<float> > matchedProtMapReverse;
      if (alignall) {
        DEBUG_MSG("Inserting dummy start mass for every DB index")
        // If ALIGN_ALL is specified insert dummy start mass for every db idx
        // This will cause alignment to every protein in database
        for (int idx = 0; idx < m_db->size(); idx++) {
          matchedProtMapForward[idx].insert(0.0);
          matchedProtMapReverse[idx].insert(0.0);
        }
      } else if (inputSpectra[specIdx].psmList.size() != 0) {
        list<psmPtr>::iterator itr = inputSpectra[specIdx].psmList.begin();
        list<psmPtr>::iterator itrEnd = inputSpectra[specIdx].psmList.end();
        for (; itr != itrEnd; itr++) {
          int idx = (*itr)->m_dbIndex;
          // Only use the forward orientation tags(if tag parsimony was not done)
          if (tagParsimony || (*itr)->m_matchOrientation == 0) {
            if (matchedProtMapForward.find(idx) == matchedProtMapForward.end()) {
              set<float> newList;
              matchedProtMapForward[idx] = newList;
            }
            matchedProtMapForward[idx].insert((*itr)->m_startMass * AA_ROUNDING);
            if (DEBUG_SPECPROTALIGN_RANGE) DEBUG_MSG(idx << "  " << m_db->getID(idx) << "  " << (*itr)->m_startMass)
          }
          if (tagParsimony || (*itr)->m_matchOrientation == 1) {
            if (matchedProtMapReverse.find(idx) == matchedProtMapReverse.end()) {
              set<float> newList;
              matchedProtMapReverse[idx] = newList;
            }
            matchedProtMapReverse[idx].insert((*itr)->m_startMass * AA_ROUNDING);
            if (DEBUG_SPECPROTALIGN_RANGE) DEBUG_MSG(idx << "  " << m_db->getID(idx) << "  " << (*itr)->m_startMass)
          }
        } // for (; itr != itrEnd; itr++)
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(matchedProtMapForward.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(matchedProtMapReverse.size());

      DEBUG_MSG("Matching as b...");
      clock_t startTime = clock();
      // -----------------------------------------------------
      // Forward matching
      // -----------------------------------------------------
      DEBUG_VAR(matchedProtMapReverse.size());
      map<int, set<float> >::iterator itrMap = matchedProtMapForward.begin();
      map<int, set<float> >::iterator itrMapEnd = matchedProtMapForward.end();
      for (; itrMap != itrMapEnd; itrMap++) {

        int protIdx = itrMap->first;
        if ( m_params.exists("DEBUG_SPECPROTALIGN_SINGLEPROTID")) {
          int protIdDebug = m_params.getValueInt("DEBUG_SPECPROTALIGN_SINGLEPROTID");
          if (protIdDebug != protIdx) {
            continue;
          }
        }

        //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(protIdx);
        // If parsimony was done in TagSearch then don't use start positions
        if (tagParsimony) {
          itrMap->second.clear();
        }
        //DEBUG_VAR(protIdx);

        Spectrum dbSpec;
        dbSpec.reserve(1024);
        dbSpec = m_db->getMassesSpec(protIdx);
        // Round all the db spectra to "center them at 0.0"
        dbSpec.roundPeaks(AA_ROUNDING, false, false);

        if (!penaltyAlign) {
          scoreOverlapAMM(cSpec,
                          dbSpec,
                          protIdx,
                          0,
                          itrMap->second,
                          maxNumMods,
                          minNumMatchPeaks,
                          pmTol,
                          peakTol,
                          57,
                          maxModMass,
                          minModMass,
                          enforceEndpeaks);
        } else {
          string dbString = m_db->getSequence(protIdx);
          apb->computeAlignment(cSpec,
                             dbSpec,
                             m_db->getSequence(protIdx),
                             protIdx,
                             0,
                             itrMap->second,
                             (*m_scanSpecificPenalties)[scanNum],
                             minNumMatchPeaks,
                             maxDbGapAas,
                             pmTol,
                             peakTol,
                             avgPeakIntensity,
                             enforceEndpeaks);
        }
        //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());

      } // for (; itrMap !=  itrMapEnd; itrMap++) {

      float elapsed = ((float)(clock() - startTime)) / CLOCKS_PER_SEC;
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsed << " seconds elapsed");
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_VAR(matchedProtMapForward.size());
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG((matchedProtMapForward.size() == 0 ? 0 : elapsed / (float)matchedProtMapForward.size()) << " per protein");

      list<psmPtr>::iterator litr = cSpec.psmList.begin();
      list<psmPtr>::iterator litrEnd = cSpec.psmList.end();
      for (; litr != litrEnd; litr++) {
        (*litr)->m_scanNum = cSpec.scan;
        (*litr)->m_protein = m_db->getID((*litr)->m_dbIndex);
        m_allPsms.push_back(*litr);
      }

      DEBUG_MSG("Matching as y...");
      startTime = clock();
      // -----------------------------------------------------
      // Reverse matching
      // -----------------------------------------------------
      DEBUG_VAR(matchedProtMapReverse.size());
      itrMap = matchedProtMapReverse.begin();
      itrMapEnd = matchedProtMapReverse.end();
      for (; itrMap != itrMapEnd; itrMap++) {

        int protIdx = itrMap->first;
        if ( m_params.exists("DEBUG_SPECPROTALIGN_SINGLEPROTID")) {
          int protIdDebug = m_params.getValueInt("DEBUG_SPECPROTALIGN_SINGLEPROTID");
          if (protIdDebug != protIdx) {
            continue;
          }
        }

        //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(protIdx);
        // If parsimony was done in TagSearch then don't use start positions
        if (tagParsimony) {
          itrMap->second.clear();
        }
        //DEBUG_VAR(protIdx);

        Spectrum dbSpec;
        dbSpec.reserve(1024);
        dbSpec = m_db->getMassesSpec(protIdx);
        // Round all the db spectra to "center them at 0.0"
        dbSpec.roundPeaks(AA_ROUNDING, false, false);
        if (!penaltyAlign) {
          scoreOverlapAMM(cSpecRev,
                          dbSpec,
                          protIdx,
                          1,
                          itrMap->second,
                          maxNumMods,
                          minNumMatchPeaks,
                          pmTol,
                          peakTol,
                          57,
                          maxModMass,
                          minModMass,
                          enforceEndpeaks);
        } else {
          apb->computeAlignment(cSpecRev,
                               dbSpec,
                               m_db->getSequence(protIdx),
                               protIdx,
                               1,
                               itrMap->second,
                               (*m_scanSpecificPenalties)[scanNum],
                               minNumMatchPeaks,
                               maxDbGapAas,
                               pmTol,
                               peakTol,
                               avgPeakIntensity,
                               enforceEndpeaks);
        }
        //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());

      } // for (; itr != itrEnd; itr++) {

      elapsed = ((float)(clock() - startTime)) / CLOCKS_PER_SEC;
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsed << " seconds elapsed");
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_VAR(matchedProtMapReverse.size());
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG((matchedProtMapReverse.size() == 0 ? 0 : elapsed / (float)matchedProtMapReverse.size()) << " per protein");

      litr = cSpecRev.psmList.begin();
      litrEnd = cSpecRev.psmList.end();
      for (; litr != litrEnd; litr++) {
        (*litr)->m_scanNum = cSpecRev.scan;
        (*litr)->m_protein = m_db->getID((*litr)->m_dbIndex);
        m_allPsms.push_back(*litr);
      }

      DEBUG_MSG("Matching Complete");

      if (m_allPsms.size() > 0 ) {
        sort(m_allPsms.m_psmSet.begin(), 
             m_allPsms.m_psmSet.end(), 
             PsmDbIndexAnnoSort);
      }
      
      // -----------------------------------------------------
      // Determine the top PSMs to pass on
      // -----------------------------------------------------
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());

      // Find the top score
      float bestScore = -(float)INT_MAX;
      litr = cSpec.psmList.begin();
      litrEnd = cSpec.psmList.end();
      for (; litr != litrEnd; litr++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(*litr);
        if ((*litr)->m_score > bestScore) {
          bestScore = (*litr)->m_score;
        }
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(bestScore);

      bool bestScoreForward = true;
      litr = cSpecRev.psmList.begin();
      litrEnd = cSpecRev.psmList.end();
      for (; litr != litrEnd; litr++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(*litr);
        if ((*litr)->m_score > bestScore) {
          bestScore = (*litr)->m_score;
          bestScoreForward = false;
        }
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(bestScore);

      float thresholdScore = bestScore > 0 ? bestScore * thresholdPercent
          : bestScore / thresholdPercent;
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(thresholdScore);
      cSpec.psmList.remove_if(LessThanPredicate(thresholdScore));
      cSpecRev.psmList.remove_if(LessThanPredicate(thresholdScore));
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(bestScoreForward);
      if (bestScoreForward) {
        // "Un-round" all the spectra to put them back like they were before
        //cSpec.roundPeaks(1.0 / AA_ROUNDING, false, false);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
        // Add all the reverse PSMs
        for (list<psmPtr>::iterator iter = cSpecRev.psmList.begin(); iter
            != cSpecRev.psmList.end(); iter++) {
          cSpec.psmList.push_back(*iter);
        }
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
        m_matchedSpectraAll->push_back(cSpec);
      }
      else
      {
        // "Un-round" all the spectra to put them back like they were before
        //cSpecRev.roundPeaks(1.0 / AA_ROUNDING, false, false);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());
        // Add all the forward PSMs
        for (list<psmPtr>::iterator iter = cSpec.psmList.begin(); iter
            != cSpec.psmList.end(); iter++) {
          cSpecRev.psmList.push_back(*iter);
        }
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());
        m_matchedSpectraAll->push_back(cSpecRev);
      }

    } // for (specIdx = 0; specIdx < inputSpectra.size(); specIdx++) {

    float elapsedTotal = ((float)(clock() - startTimeTotal)) / CLOCKS_PER_SEC;
    if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsedTotal << " total seconds elapsed");
    if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsedTotal / (float)inputSpectra.size() << " average per spectrum");
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectraAll->size());

    // -----------------------------------------------------
    // We could have multiple hits on one protein
    // So we perform a sort() and unique() to eliminate them
    // -----------------------------------------------------
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      //if (DEBUG_SPECPROTALIGN)DEBUG_VAR(i);
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
      (*m_matchedSpectraAll)[i].psmList.sort(PsmDbIndexAnnoSort);
      (*m_matchedSpectraAll)[i].psmList.unique(PsmDbIndexAnnoUnique);
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
    }

    // -----------------------------------------------------
    // Set backwards pointers/scans for PSM--->Spectrum link
    // -----------------------------------------------------
    list<psmPtr>::iterator iter;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      //(*m_matchedSpectraAll)[i].scan = startIdx + i + 1;
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].scan);
      for (iter = (*m_matchedSpectraAll)[i].psmList.begin(); iter
          != (*m_matchedSpectraAll)[i].psmList.end(); iter++) {
        (*iter)->m_spectrum = &(*m_matchedSpectraAll)[i];
        //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*iter)->m_matchedPeaks.size());

        string dbSeqString = m_db->getSequence((*iter)->m_dbIndex);
        if (tagParsimony) {
          string annotation;
          (*iter)->getAnnotationFromMatchedPeaks(m_db->getMassesSpec((*iter)->m_dbIndex),
                                                 dbSeqString,
                                                 annotation);
          (*iter)->m_annotation = (*iter)->m_origAnnotation;
          (*iter)->m_annotation = annotation;
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*iter)->m_annotation);
        }

        // For now copy annotation to original annotation(we'll fix this later)
        (*iter)->m_origAnnotation = (*iter)->m_annotation;
        //        (*iter)->m_annotation = stringAnnotation;

        vector<float> modifications;
        (*iter)->getModifications(modifications);
        (*iter)->m_numMods = modifications.size();
        (*iter)->m_charge = cSpec.parentCharge;
        if (clusteredSpecFileName.empty()) {
          (*iter)->m_spectrumFile = cSpec.fileName;
        } else {
          (*iter)->m_spectrumFile = clusteredSpecFileName;
        }
      }
    }

    if (enforceEndpeaks) {   // This is a proxy (for now) to check for star alignment only
      PeptideSpectrumMatchSet psmSetTemp;
      psmSetTemp.getPSMSet(m_matchedSpectraAll);
      for (int i = 0; i < psmSetTemp.size(); i++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(psmSetTemp[i]);
      }

      // Get 1 best scoring PSM for each scan
      map<int, psmPtr> psmMap;      
      for (int i = 0; i < psmSetTemp.size(); i++) {
        int scanNum = psmSetTemp[i]->m_scanNum;
        map<int, psmPtr>::iterator itrFind = psmMap.find(scanNum);
        if (itrFind == psmMap.end() || psmSetTemp[i]->m_score > itrFind->second->m_score) {
          psmMap[scanNum] = psmSetTemp[i];
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR(scanNum);
          if (DEBUG_SPECPROTALIGN) DebugPsm(psmMap[scanNum]);
        }
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(psmMap.size());

      // Convert the map of best PSMs back to a PSM set
      PeptideSpectrumMatchSet psmSetTemp2;
      map<int, psmPtr>::iterator itrMap = psmMap.begin();
      map<int, psmPtr>::iterator itrMapEnd = psmMap.end();
      for (; itrMap != itrMapEnd; itrMap++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(itrMap->second);
        psmSetTemp2.push_back(itrMap->second);
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(psmSetTemp2.size());

      // Now add back in all the PSMs that were equally as good (and yet different)
      //   These might have better spectral probs
      for (int i = 0; i < psmSetTemp.size(); i++) {
        int scanNum = psmSetTemp[i]->m_scanNum;
        if (psmMap[scanNum]->m_score == psmSetTemp[i]->m_score &&
            psmMap[scanNum]->m_annotation != psmSetTemp[i]->m_annotation) {
          if (DEBUG_SPECPROTALIGN) DebugPsm(psmSetTemp[i]);
          psmSetTemp2.push_back(psmSetTemp[i]);
        }
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(psmSetTemp2.size());

      for (int i = 0; i < psmSetTemp2.size(); i++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(psmSetTemp2[i]);
      }

      bool useRawScoreForProb = m_params.getValueBool("USE_ALIGNMENT_SCORE_FOR_PROBABILITY", true);
      DEBUG_VAR(useRawScoreForProb);

#if 0
      //----------------------------------------------
      // Now recompute using sprinkled spectrum
      //----------------------------------------------
      if (penaltyAlign) {
        recomputeWithSprinkled(psmSetTemp2, specIdx, apb);
      }

      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(psmSetTemp2.size());
      for (int i = 0; i < psmSetTemp2.size(); i++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(psmSetTemp2[i]);
      }
#endif

      //------------------------------------------
      // Compute the spectral probabilities
      //------------------------------------------
      computeSpectralProbabilities(apb, &psmSetTemp2);
      for (int i = 0; i < psmSetTemp2.size(); i++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(psmSetTemp2[i]);
      }

      //------------------------------------------
      // Find the PSMs with the best p-value now
      //------------------------------------------
      psmMap.clear();
      for (int i = 0; i < psmSetTemp2.size(); i++) {
        int scanNum = psmSetTemp2[i]->m_scanNum;
        map<int, psmPtr>::iterator itrFind = psmMap.find(scanNum);
        if (itrFind == psmMap.end() || psmSetTemp2[i]->m_pValue > itrFind->second->m_pValue) {
          psmMap[scanNum] = psmSetTemp2[i];
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR(scanNum);
          if (DEBUG_SPECPROTALIGN) DebugPsm(psmMap[scanNum]);
        }
      }

      // Convert the map of best PSMs back to a PSM set
      psmSetTemp.m_psmSet.clear();
      for (itrMap = psmMap.begin(); itrMap != itrMapEnd; itrMap++) {
      
        string annotationMergedMods;
        mergeModifications(itrMap->second, annotationMergedMods);
        itrMap->second->m_annotation = annotationMergedMods;
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(itrMap->second->m_annotation);

        psmSetTemp.push_back(itrMap->second);
      }

      m_matchedSpectraAll->clearPsms(); // Clear all PSMs
      psmSetTemp.addSpectra(m_matchedSpectraAll); // Associate the new PSMs (with pvalue)
    }

    // Free up the penalty alignment object 
    if (penaltyAlign) {
      delete apb;
    }
      
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecProtAlign::loadInputData(void)
  {
    if (ownInput) {
      if (!m_inputSpectra)
        m_inputSpectra = new SpecSet;
      if (!m_prmSpectra)
        m_prmSpectra = new SpecSet;
      if (!m_db)
        m_db = new DB_fasta;
      if (!m_filterPsmSet)
        m_filterPsmSet = new PeptideSpectrumMatchSet;
      if (!m_contigAbinfo) {
        m_contigAbinfo = new abinfo_t;
      }
    }
    m_inputSpectra->resize(0);

    if (ownOutput) {
      if (!m_matchedSpectraAll)
        m_matchedSpectraAll = new SpecSet;
    }
    m_matchedSpectraAll->resize(0);

    // Load amino acid masses
    AAJumps jumps(1);
    if (m_params.exists("AMINO_ACID_MASSES")) {
      string filename = m_params.getValue("AMINO_ACID_MASSES");
      DEBUG_MSG("Loading: " << filename);
      jumps.loadJumps(filename.c_str(), true);
    }

    if (m_params.exists("INPUT_SPECS_PKLBIN")) {
      string filename = m_params.getValue("INPUT_SPECS_PKLBIN");
      DEBUG_MSG("Loading: " << filename);
      if (m_inputSpectra->loadPklBin(filename.c_str()) < 0) {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
      string psmFilename = m_params.getValue("INPUT_PSM");
      if (!psmFilename.empty()) {
        int scanFirst = m_params.getValueInt("SCAN_FIRST", -1);
        int scanLast = m_params.getValueInt("SCAN_LAST", -1);
        DEBUG_VAR(scanFirst);
        DEBUG_VAR(scanLast);
        PeptideSpectrumMatchSet psmSet;
        string filename = m_params.getValue("INPUT_PSM");
        DEBUG_MSG("Loading: " << filename);
        if (!psmSet.loadFromFile(filename.c_str(), scanFirst, scanLast)) {
          ERROR_MSG("Error reading input PSM file.");
          return false;
        }
        psmSet.addSpectra(m_inputSpectra);
      } 
    }

    if (m_params.exists("INPUT_PRM_PKLBIN")) {
      DEBUG_MSG("Loading: " << m_params.getValue("INPUT_PRM_PKLBIN"));
      if (m_prmSpectra->loadPklBin(m_params.getValue("INPUT_PRM_PKLBIN").c_str()) < 0) {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    // Create decoy database on the fly (in place)
    bool usingDecoyDb = m_params.getValueInt("USE_DECOY_DATABASE", 0);
    DEBUG_VAR(usingDecoyDb);
    if (!usingDecoyDb) {
      string filename = m_params.getValue("FASTA_DATABASE");
      DEBUG_MSG("Loading: " << filename);
      if (!m_params.exists("FASTA_DATABASE")) {
        ERROR_MSG("Parameters are incomplete. FASTA_DATABASE is missing.");
        return false;
      } else if (m_db->Load(filename.c_str()) <= 0) {
        ERROR_MSG("Error reading database sequences from "
            << m_params.getValue("FASTA_DATABASE"));
        return false;
      }
    } else {
      string filename = m_params.getValue("FASTA_DATABASE_DECOY");
      DEBUG_MSG("Loading: " << filename);
      if (!m_params.exists("FASTA_DATABASE_DECOY")) {
        ERROR_MSG("Parameters are incomplete. FASTA_DATABASE_DECOY is missing.");
        return false;
      } else if (m_db->Load(filename.c_str()) <= 0) {
        ERROR_MSG("Error reading database sequences from "
            << m_params.getValue("FASTA_DATABASE_DECOY"));
        return false;
      }
    }
    
    if (m_params.exists("INPUT_ABINFO")) {
      string filename = m_params.getValue("INPUT_ABINFO");
      DEBUG_MSG("Loading Abinfo.. [" << filename << "]");
      if (!Load_abinfo(filename.c_str(), *m_contigAbinfo)) {
        ERROR_MSG("Error reading contig abinfo file [" << filename.c_str() << "]");
        return false;
      }
    }

    string psmFilename = m_params.getValue("INPUT_FILTER_PSMS");
    if (!psmFilename.empty()) {
      string filename = m_params.getValue("INPUT_FILTER_PSMS");
      DEBUG_MSG("Loading PSM file [" << filename << "]...");
      if (!m_filterPsmSet->loadFromFile(filename.c_str())) {
        ERROR_MSG("Error reading input filter PSM file.");
        return false;
      }
    } 
    
    psmFilename = m_params.getValue("INPUT_FILTER_MSGFDB_PSMS");
    if (!psmFilename.empty()) {
      string filename = m_params.getValue("INPUT_FILTER_MSGFDB_PSMS");
      DEBUG_MSG("Loading MSGFDB PSM file [" << filename << "]...");
      if (!m_filterPsmSet->loadMSGFDBResultsFile(filename.c_str())) {
        ERROR_MSG("Error reading input filter MSGFDB PSM file.");
        return false;
      }
    } 

    //---------------------------------------------------------------------------
    // Load penalty matrices for new alignment
    //---------------------------------------------------------------------------

    bool penaltyAlign = m_params.exists("PENALTY_ALIGNMENT")
        ? m_params.getValueBool("PENALTY_ALIGNMENT", false) : false;
    DEBUG_VAR(penaltyAlign);
    if (penaltyAlign) {
      float resolution = m_params.getValueDouble("ALIGNMENT_RESOLUTION", 1.0);
      float knownModPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_KNOWN_PENALTY", 0.01);
      float unknownPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_PENALTY", 1.0);
      float unknownMultiplier = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER", 1.5);
      float minModMass = m_params.getValueDouble("MIN_MOD_MASS", -100.0);
      float maxModMass = m_params.getValueDouble("MAX_MOD_MASS", 100.0);
      DEBUG_VAR(resolution);
      DEBUG_VAR(knownModPenalty);
      DEBUG_VAR(unknownPenalty);
      DEBUG_VAR(unknownMultiplier);
      DEBUG_VAR(minModMass);
      DEBUG_VAR(maxModMass);
      
      if (!m_penaltyMatrixBlosum) {
        m_penaltyMatrixBlosum = new PenaltyMatrix(jumps,
                                    resolution,
                                    knownModPenalty, 
                                    unknownPenalty,
                                    unknownMultiplier,
                                    minModMass,
                                    maxModMass);
      }
      if (!m_penaltyMatrixMods) {
        m_penaltyMatrixMods = new PenaltyMatrix(jumps,
                                    resolution,
                                    knownModPenalty, 
                                    unknownPenalty,
                                    unknownMultiplier,
                                    minModMass,
                                    maxModMass);
      }

      if (!m_params.exists("BLOSUM_PENALTY_FILE")) {
        WARN_MSG("KNOWN_BLOSUM_PENALTY is missing.");
      } else {
        string blosumFileName = m_params.getValue("BLOSUM_PENALTY_FILE");
        string knowmModsFileName = m_params.getValue("KNOWN_MODS_FILE");
        string cleavagePenaltiesFileName = m_params.getValue("CLEAVAGE_PENALTY_FILE");
        DEBUG_MSG("Loading: " << blosumFileName);
        DEBUG_MSG("Loading: " << knowmModsFileName);
        DEBUG_MSG("Loading: " << cleavagePenaltiesFileName);
        if (!m_penaltyMatrixBlosum->load(blosumFileName,
                                         knowmModsFileName,
                                         cleavagePenaltiesFileName)) {
          ERROR_MSG("Error loading blosum penalties from " << blosumFileName);
          return false;
        }
      }

      if (!m_params.exists("MODS_PENALTY_FILE")) {
        WARN_MSG("MODS_PENALTY_FILE is missing.");
      } else {
        string modFileName = m_params.getValue("MODS_PENALTY_FILE");
        string knowmModsFileName = m_params.getValue("KNOWN_MODS_FILE");
        string cleavagePenaltiesFileName = m_params.getValue("CLEAVAGE_PENALTY_FILE");
        DEBUG_MSG("Loading: " << modFileName);
        DEBUG_MSG("Loading: " << knowmModsFileName);
        DEBUG_MSG("Loading: " << cleavagePenaltiesFileName);
        if (!m_penaltyMatrixMods->load(modFileName,
                                       knowmModsFileName,
                                       cleavagePenaltiesFileName)) {
          ERROR_MSG("Error loading mod penalties from " << modFileName);
          return false;
        }
      }

      if (!m_scanSpecificPenalties)
        m_scanSpecificPenalties = new map<int, map<int, float> >();

      if (m_params.exists("SCAN_SPECIFIC_PENALTIES_FILE")) {
        string filename = m_params.getValue("SCAN_SPECIFIC_PENALTIES_FILE");
        DEBUG_MSG("Loading: " << filename);
        loadScanSpecificPenalties(filename,
                                  *m_scanSpecificPenalties);
        DEBUG_VAR(m_scanSpecificPenalties->size());
      }

    } // if (penaltyAlign)

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecProtAlign::saveOutputData(void)
  {
    string outDir;
    bool isChild = m_params.getValueBool("IS_CHILD"); 
    if (isChild) {
      DEBUG_MSG("SAVING CHILD DATA");
      outDir = m_params.getValue("GRID_DATA_DIR_OUT");
    } else {
      DEBUG_MSG("SAVING PARENT DATA");
      outDir = m_params.getValue("GRID_DATA_DIR");
    }
    if (outDir.empty()) {
      outDir = ".";
    } 
    outDir = outDir + "/";
    DEBUG_VAR(outDir);

    bool penaltyAlign = m_params.getValueBool("PENALTY_ALIGNMENT", false);
    DEBUG_VAR(penaltyAlign);

    if (m_params.exists("OUTPUT_MATCHED_PEAKS_IDX_ALL")) {

      AAJumps aaJumps(1);
      if (m_params.exists("AMINO_ACID_MASSES")) {
        DEBUG_MSG("Loading amino acid masses...");
        string fileName = m_params.getValue("AMINO_ACID_MASSES");
        DEBUG_VAR(fileName);
        if (!aaJumps.loadJumps(fileName.c_str(),
                               false))
        {
          ERROR_MSG("Error reading input amino acid mass file " << m_params.getValue("AMINO_ACID_MASSES"));
          return false;
        }
      }

      bool tagParsimony = m_params.getValueBool("MAX_PARSIMONY", true);
      DEBUG_VAR(tagParsimony);

      DEBUG_MSG("Saving matched peaks all...");
      SpecSet tempMatchedPeaks;
      tempMatchedPeaks.resize(m_matchedSpectraAll->size());
      for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
        //DEBUG_VAR(i);
        //DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
        if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
          list<psmPtr>::iterator litr =
              (*m_matchedSpectraAll)[i].psmList.begin();

          if (tagParsimony) {
            //DEBUG_VAR((*litr)->m_annotation);
            //if (penaltyAlign) DEBUG_VAR((*litr)->m_origAnnotation);
            Spectrum locDbSpec = m_db->getMassesSpec((*litr)->m_dbIndex);
            locDbSpec.roundPeaks(AA_ROUNDING, false, false);

            (*litr)->getMatchedPeaksFromAnnotation(locDbSpec, aaJumps);
            //DEBUG_VAR((*litr)->m_annotation);
          }

          int peakListSize = (*litr)->m_matchedPeaks.size();
          //DEBUG_VAR(peakListSize);
          tempMatchedPeaks[i].resize(peakListSize);
          for (int j = 0; j < peakListSize; j++) {
            tempMatchedPeaks[i][j].set((*litr)->m_matchedPeaks[j][0],
                                       (*litr)->m_matchedPeaks[j][1]);
            //DEBUG_MSG((*litr)->m_matchedPeaks[j][0] << "  " << (*litr)->m_matchedPeaks[j][1]);
          }
        }
      }
      string saveName = outDir + "/" + m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str();
      tempMatchedPeaks.savePklBin(saveName.c_str());
    }

    DEBUG_VAR(m_matchedSpectraAll->size());
    string specsName = outDir + m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str();
    string psmName = outDir + m_params.getValue("OUTPUT_PSM_ALL").c_str();
    string peaksName = outDir + m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str();
    DEBUG_VAR(specsName);
    DEBUG_VAR(psmName);
    DEBUG_VAR(peaksName);

    if (m_matchedSpectraAll and m_params.exists("OUTPUT_MATCHED_SPECS_ALL")) {
      if (m_params.exists("OUTPUT_PSM_ALL")
          && m_params.exists("OUTPUT_MATCHED_PEAKS_IDX_ALL")) {
        DEBUG_MSG("Saving matched specs all (3 files)...");
        m_matchedSpectraAll->savePklBin(specsName.c_str(),
                                        psmName.c_str(),
                                        peaksName.c_str());

      } else if (m_params.exists("OUTPUT_PSM_ALL")) {
        DEBUG_MSG("Saving matched specs all (2 files)...");
        DEBUG_VAR(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str());
        DEBUG_VAR(m_params.getValue("OUTPUT_PSM_ALL").c_str());
        m_matchedSpectraAll->savePklBin(specsName.c_str(),
                                        psmName.c_str());

      } else {
        DEBUG_MSG("Saving matched specs all (1 file)...");
        m_matchedSpectraAll->savePklBin(specsName.c_str());
      }
    }

    if (m_params.exists("OUTPUT_PSM_ALL")) {
      m_allPsms.saveToFile((psmName + ".debug").c_str());
    }
    DEBUG_VAR(m_params.getValue("OUTPUT_SPECTRA_SPRINKLED"));
    if (m_matchedSpectraAll && m_params.exists("OUTPUT_SPECTRA_SPRINKLED")) {
      string sprinkName = outDir + m_params.getValue("OUTPUT_SPECTRA_SPRINKLED").c_str();
      m_matchedSpectraAllSprinkled.savePklBin(sprinkName.c_str());
      DEBUG_MSG("Saving matched sprinkled specs all (1 file)...");
    }
    
    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecProtAlign::saveInputData(std::vector<std::string> & filenames)
  {
    string inDir = m_params.getValue("GRID_DATA_DIR_IN");
    if (inDir.empty()) {
      inDir = ".";
    }
    string baseDirectory = inDir + "/";
    DEBUG_VAR(baseDirectory);
    string baseFilename = baseDirectory + getName();

    string specFilename;
    bool enforceEndpeaks = m_params.getValueInt("ENFORCE_ENDPEAKS");
    if (!enforceEndpeaks) {
      specFilename = baseFilename + "_stars.pklbin";  // Individual files for contigs
    } else {
      specFilename = baseDirectory + "stars.pklbin";  // Just one for all nodes for spectra
    }

    string psmFilename = baseFilename + "_tags.txt";

    if (!fileExists(specFilename)) {
      DEBUG_MSG("Saving " << specFilename);
      m_inputSpectra->savePklBin(specFilename.c_str());
    } else {
      DEBUG_MSG("Not Saving " << specFilename << " (already exists)");
    }
    m_params.setValue("INPUT_SPECS_PKLBIN", specFilename);

    string dbDecoyFilename = baseDirectory + "decoy.fasta";
    if (!fileExists(dbDecoyFilename)) {
      DB_fasta dbDecoy;
      dbDecoy = *m_db;
      dbDecoy.replaceDecoyShuffled();
      dbDecoy.Save(dbDecoyFilename.c_str());
    } else {
      DEBUG_MSG("Not Saving " << dbDecoyFilename << " (already exists)");
    }
    m_params.setValue("FASTA_DATABASE_DECOY", dbDecoyFilename);

    string prmFilename = baseDirectory + "prm_spectra.pklbin";
    if (enforceEndpeaks) {
      if (m_prmSpectra != 0x0) {
        m_params.setValue("INPUT_PRM_PKLBIN", prmFilename );
        if (!fileExists(prmFilename))
        {
          DEBUG_MSG("Saving " << prmFilename );
          m_prmSpectra->savePklBin(prmFilename.c_str());
        } else {
          DEBUG_MSG("Not Saving " << prmFilename  << " (already exists)");
        }
      }
    }


    m_params.removeParam("INPUT_FILTER_PSMS"); 
    m_params.removeParam("INPUT_FILTER_MSGFDB_PSMS"); 
    
    if (!enforceEndpeaks) {
      m_params.setValue("OUTPUT_TAG_PSM", getName() + "_tags.txt");
    }

    string paramDir = m_params.getValue("GRID_DATA_DIR_PARAMS");
    if (paramDir.empty()) {
      paramDir = ".";
    }
    baseDirectory = paramDir + "/";
    baseFilename = baseDirectory + getName();

    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(specFilename);
    filenames.push_back(psmFilename);
    filenames.push_back(dbDecoyFilename);

    return true;
  }
  

  // -------------------------------------------------------------------------
  bool ExecSpecProtAlign::loadOutputData(void)
  {
    if (m_matchedSpectraAll == 0x0) {
      ownOutput = true;
      m_matchedSpectraAll = new SpecSet;
    }

    string outDir = m_params.getValue("GRID_DATA_DIR_INTERMEDIATE");
    if (outDir.empty()) {
      outDir = ".";
    }
    outDir = outDir + "/";
    DEBUG_VAR(outDir);

    DEBUG_VAR(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str());
    DEBUG_VAR(m_matchedSpectraAll->size());
    string specsName = outDir + m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str();
    string psmName = outDir + m_params.getValue("OUTPUT_PSM_ALL").c_str();
    string peaksName = outDir + m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str();
    
    if (m_params.exists("OUTPUT_MATCHED_SPECS_ALL")) {
      m_matchedSpectraAll->loadPklBin(specsName.c_str(),
                                      psmName.c_str(),
                                      peaksName.c_str());
      DEBUG_VAR(m_matchedSpectraAll->size());				      
    }

    string sprinkledName = outDir + m_params.getValue("OUTPUT_SPECTRA_SPRINKLED").c_str();
    
    DEBUG_VAR(sprinkledName);
    if (m_params.exists("OUTPUT_SPECTRA_SPRINKLED")) {
      m_matchedSpectraAllSprinkled.loadPklBin(sprinkledName.c_str());
      DEBUG_VAR(m_matchedSpectraAllSprinkled.size());				      
    }

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::splitContigs(int numSplit,
                                       int startBaseIdx, 
                                       int endBaseIdx)
  {
    int spectraSize = m_inputSpectra->size();
    int numSpectraPerSplit = (endBaseIdx - startBaseIdx + 1) / numSplit;
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(numSpectraPerSplit);
    int extraSpectra = (endBaseIdx - startBaseIdx + 1) % numSplit;
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(extraSpectra);

    int indexStart = startBaseIdx;
    int indexEnd = startBaseIdx + numSpectraPerSplit - 1;
    if (extraSpectra > 0) {
      indexEnd++;
      extraSpectra--;
    }

    for (int i = 0; i < numSplit; i++) {
      if (startBaseIdx >= spectraSize) {
        break;
      }

      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(i);

      // Copy the parameters
      ParameterList childParams(m_params);
      childParams.removeParam("GRID_EXECUTION"); // necessary for Proteosafe
      // The child params must know only about the input dir
      
      // Set the start and end indices
      char buf[128];
      sprintf(buf, "%d", indexStart);
      childParams.setValue("IDX_START", buf);
      sprintf(buf, "%d", indexEnd);
      childParams.setValue("IDX_END", buf);

      SpecSet * starSpecSet = new SpecSet;
      for (int iSpec = indexStart; iSpec <= indexEnd; iSpec++) {
        starSpecSet->push_back((*m_inputSpectra)[iSpec]);
      }

      sprintf(buf, "%d", (*m_inputSpectra)[indexStart].scan );
      childParams.setValue("SCAN_FIRST", buf);
      sprintf(buf, "%d", (*m_inputSpectra)[indexEnd].scan);
      childParams.setValue("SCAN_LAST", buf);

      // Make a clone of this module
      ExecSpecProtAlign * theClone = new ExecSpecProtAlign(childParams,
                                                  starSpecSet,
                                                  m_prmSpectra,
                                                  m_db,
                                                  m_penaltyMatrixBlosum,
                                                  m_penaltyMatrixMods,
                                                  m_scanSpecificPenalties,
                                                  m_contigAbinfo, 
                                                  m_filterPsmSet);
      theClone->m_params.setValue("IS_CHILD", "1"); 
      
      // Special flag because we are NOT sharing the input spectra with our parent
      theClone->ownSpectra = true;

      // Give it a new name based on the split
      theClone->setName(makeName(m_name, i));

      // Have to set up the output files also so the params will be correct on reload
      string baseName = theClone->getName();
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(baseName);

      theClone->m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", 
                                  baseName + "_contigs_all.pklbin");
      theClone->m_params.setValue("OUTPUT_PSM_ALL",
                                  baseName + "_psm_all.txt");
      theClone->m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", 
                                  baseName + "_midx_all.pklbin");
      //theClone->m_params.removeParam("OUTPUT_MATCHED_PEAKS_IDX_ALL");
      theClone->m_params.setValue("OUTPUT_SPECTRA_SPRINKLED", 
                                  baseName + "_sprinkled.pklbin");

      std::string suffix("");
      char bufSplit[128];
      sprintf(bufSplit, "%d", i + 1);
      theClone->m_params.setValue("NUM_SPLIT", bufSplit);

      m_subModules.push_back(theClone);

      indexStart = indexEnd + 1;
      indexEnd = indexStart + numSpectraPerSplit - 1;
      if (extraSpectra > 0) {
        indexEnd++;
        extraSpectra--;
      }

    } // for (int i = 0; i < numSplit; i++)
  }
  
  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::splitStars(int numSplit,
                                     int startBaseIdx, 
                                     int endBaseIdx)
  {
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_TRACE;
    PeptideSpectrumMatchSet psmSet;
    map<int, int> mapPsmCounts;
    int  totalPsms;

    int startScan = (*m_inputSpectra)[startBaseIdx].scan;
    int endScan = (*m_inputSpectra)[endBaseIdx].scan;
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(startScan);
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(endScan);

    string fileName = m_params.getValue("INPUT_TAGS_FOR_SPLITTING");
    DEBUG_VAR(fileName);
    
    psmSet.loadSizesFromFile(fileName.c_str(),
                             mapPsmCounts,
                             totalPsms,
                             startScan,
                             endScan);

    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(totalPsms);

    if (DEBUG_SPECPROTALIGN_SPLIT) {
      int total = 0;
      DEBUG_VAR(mapPsmCounts.size());
      map<int, int>::iterator itrMap = mapPsmCounts.begin();
      map<int, int>::iterator itrMapEnd = mapPsmCounts.end();
      for ( ; itrMap != itrMapEnd; itrMap++) {
        //DEBUG_MSG(itrMap->first << "  " << itrMap->second);
        total += itrMap->second;
      }
      DEBUG_VAR(total);
    }

    int psmsPerSplit = totalPsms / numSplit;
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(psmsPerSplit);

    if (psmsPerSplit == 0) {
      DEBUG_MSG("No PSMs to split");
      return;
    }

    int usedPsms = 0;
    int thisPsms = 0;
    int thisStartIdx = startBaseIdx;
    int thisSplit = 0;

    for (int i = startBaseIdx; i <= endBaseIdx; i++) {
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(i);
      int scanNum = (*m_inputSpectra)[i].scan;
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(scanNum);

      usedPsms += mapPsmCounts[scanNum];
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(usedPsms);
      thisPsms += mapPsmCounts[scanNum];
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(thisPsms);

      // Keep gathering PSMs until we have enough for one split
      if (thisPsms < psmsPerSplit && i != endBaseIdx) {
        continue;
      }

      thisPsms = 0;

      // Copy the parameters
      ParameterList childParams(m_params);
      childParams.removeParam("GRID_EXECUTION"); // necessary for Proteosafe
      // The child params must know only about the input dir

      // Set the start and end indices
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_MSG(thisStartIdx << " to " << i);
      char buf[128];
      sprintf(buf, "%d", thisStartIdx);
      childParams.setValue("IDX_START", buf);
      sprintf(buf, "%d", i);
      childParams.setValue("IDX_END", buf);

      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_MSG((*m_inputSpectra)[thisStartIdx].scan << " to " << (*m_inputSpectra)[i].scan);
      sprintf(buf, "%d", (*m_inputSpectra)[thisStartIdx].scan );
      childParams.setValue("SCAN_FIRST", buf);
      sprintf(buf, "%d", scanNum);
      childParams.setValue("SCAN_LAST", buf);

      // Make a clone of this module
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_MSG("Making clone...");
      ExecBase * theClone = new ExecSpecProtAlign(childParams,
                                                  m_inputSpectra,
                                                  m_prmSpectra,
                                                  m_db,
                                                  m_penaltyMatrixBlosum,
                                                  m_penaltyMatrixMods,
                                                  m_scanSpecificPenalties,
                                                  m_contigAbinfo, 
                                                  m_filterPsmSet);
      theClone->m_params.setValue("IS_CHILD", "1"); 
      
      // Give it a new name based on the split
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(thisSplit);
      theClone->setName(makeName(m_name, thisSplit));

      // Have to set up the output files also so the params will be correct on reload
      string baseName = theClone->getName();
      if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(baseName);

      theClone->m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", 
                                  baseName + "_contigs_all.pklbin");
      theClone->m_params.setValue("OUTPUT_PSM_ALL",
                                  baseName + "_psm_all.txt");
      theClone->m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", 
                                  baseName + "_midx_all.pklbin");
      //theClone->m_params.removeParam("OUTPUT_MATCHED_PEAKS_IDX_ALL");
      theClone->m_params.setValue("OUTPUT_SPECTRA_SPRINKLED", 
                                  baseName + "_sprinkled.pklbin");

      std::string suffix("");
      char bufSplit[128];
      sprintf(bufSplit, "%d", thisSplit);
      theClone->m_params.setValue("NUM_SPLIT", bufSplit);

      thisSplit++;
      thisStartIdx = i + 1; // The next start index

      m_subModules.push_back(theClone);

      // Exit if we have used all the PSMs
      if (usedPsms >= totalPsms) {
        break;
      }

    } // for (int i = startBaseIdx; i <= endBaseIdx; i++) {
  }  

  // -------------------------------------------------------------------------
  vector<ExecBase *> const & ExecSpecProtAlign::split(int numSplit)
  {
    DEBUG_SPECPROTALIGN_SPLIT = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPLIT");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPLIT);

    DEBUG_VAR(numSplit);

    if (numSplit < 2) {
      DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
      return m_subModules;
    }

    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(m_inputSpectra);
    int spectraSize = m_inputSpectra->size();
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(spectraSize);
    if (spectraSize == 0) {
      DEBUG_MSG("Must have at least one spectra");
      return m_subModules;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START")) {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    } else {
      startBaseIdx = 0;
    }
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(startBaseIdx);
    int endBaseIdx;
    if (m_params.exists("IDX_END")) {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
      endBaseIdx = min(endBaseIdx, (int)spectraSize - 1);
    } else {
      endBaseIdx = spectraSize - 1;
    }
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(endBaseIdx);

    bool enforceEndpeaks = m_params.getValueInt("ENFORCE_ENDPEAKS");
    if (DEBUG_SPECPROTALIGN_SPLIT) DEBUG_VAR(enforceEndpeaks);

    if (!enforceEndpeaks) { // This is a proxy to check for contig alignment only
      splitContigs(numSplit, startBaseIdx, endBaseIdx);
    } else {
      splitStars(numSplit, startBaseIdx, endBaseIdx);
    }

    DEBUG_VAR(m_subModules.size());

    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecProtAlign::merge(void)
  {
    if (m_matchedSpectraAll == 0x0) {
      ownOutput = true;
      m_matchedSpectraAll = new SpecSet;
    }

    DEBUG_SPECPROTALIGN_ANNO = m_params.getValueBool("DEBUG_SPECPROTALIGN_ANNO");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_ANNO);
    DEBUG_SPECPROTALIGN_MERGE = m_params.getValueBool("DEBUG_SPECPROTALIGN_MERGE");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_MERGE);

    DEBUG_VAR(m_subModules.size());
    int iSpec = 0;
    for (int i = 0; i < m_subModules.size(); i++) {
      ExecSpecProtAlign * espa = (ExecSpecProtAlign*)m_subModules[i];
      if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR(espa);
      if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR(espa->m_matchedSpectraAll->size());
      for (int j = 0; j < espa->m_matchedSpectraAll->size(); j++) {
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR(iSpec);
        m_matchedSpectraAll->push_back(espa->m_matchedSpectraAll->operator[](j));
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*m_matchedSpectraAll)[iSpec].size());
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*m_matchedSpectraAll)[iSpec].psmList.size());
        iSpec++;
      }
    }

    bool changeGapAnnos = m_params.getValueBool("CHANGE_GAP_ANNOTATIONS_TO_SINGLE");
    DEBUG_VAR(changeGapAnnos);

    // Load amino acid masses
    AAJumps jumps(1);
    if (m_params.exists("AMINO_ACID_MASSES")) {
      string fileName = m_params.getValue("AMINO_ACID_MASSES");
      jumps.loadJumps(fileName.c_str(), true);
    }
    
    // Set backwards pointers/scans for PSM--->Spectrum link
    list<psmPtr>::iterator iter;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*m_matchedSpectraAll)[i].scan);
      //(*m_matchedSpectraAll)[i].scan = i + 1;
      for (iter = (*m_matchedSpectraAll)[i].psmList.begin(); iter
          != (*m_matchedSpectraAll)[i].psmList.end(); iter++) {
        (*iter)->m_spectrum = &(*m_matchedSpectraAll)[i];
        (*iter)->m_scanNum = (*m_matchedSpectraAll)[i].scan;

        if (DEBUG_SPECPROTALIGN_MERGE) DebugPsm(*iter);

        // LARS: This is to handle issues with ProteoSAFe (can't handle gaps)
        if (changeGapAnnos) {
          string annotationOut;
          (*iter)->m_origAnnotation = (*iter)->m_annotation;
          (*iter)->changeGapAnnosToSingle();
          (*iter)->addFixedCysteineMods();
        }
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*iter)->m_annotation);
      }
    }

    DEBUG_VAR(m_subModules.size());
    for (int i = 0; i < m_subModules.size(); i++) {
      ExecSpecProtAlign * espa = (ExecSpecProtAlign*)m_subModules[i];
      if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR(espa->m_matchedSpectraAllSprinkled.size());
      for (int j = 0; j < espa->m_matchedSpectraAllSprinkled.size(); j++) {
        m_matchedSpectraAllSprinkled.push_back(espa->m_matchedSpectraAllSprinkled[j]);
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::computeSpectralProbabilities(AlignmentPenaltyBased * apb, PeptideSpectrumMatchSet * psmSet)
  {
    DEBUG_SPECPROTALIGN_SPECPROB = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPECPROB");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPECPROB);
    DEBUG_SPECPROTALIGN_TIME = m_params.getValueBool("DEBUG_SPECPROTALIGN_TIME");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_TIME);
    DEBUG_SPECPROTALIGN_SPRINKLE = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPRINKLE");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPRINKLE);

    float peakTol = (float)m_params.getValueDouble("TOLERANCE_PEAK");
    float pmTol = (float)m_params.getValueDouble("TOLERANCE_PM");
    //float peakTol = (float)m_params.getValueDouble("ALIGNMENT_RESOLUTION");
    float maxModMass = m_params.getValueDouble("MAX_MOD_MASS");
    float minModMass = m_params.getValueDouble("MIN_MOD_MASS");
    bool penaltyAlign = m_params.getValueBool("PENALTY_ALIGNMENT", false);
    float penaltyAlpha = m_params.getValueFloat("PENALTY_ALIGNMENT_ALPHA", 1.0);
    float penaltyBeta = m_params.getValueFloat("PENALTY_ALIGNMENT_BETA", 1000.0);
    bool alignall = m_params.getValueBool("ALIGN_ALL", false);
    int maxDbGapAas = m_params.getValueInt("MAX_ALIGN_DB_GAP_AAS", 8);
    int maxSpectrumGapDaltons = m_params.getValueInt("MAX_ALIGN_SPECTRUM_GAP_DALTONS", 1500);

    bool useRawScoreForProb = m_params.getValueBool("USE_ALIGNMENT_SCORE_FOR_PROBABILITY", true);
    DEBUG_VAR(useRawScoreForProb);

    DEBUG_VAR(peakTol);
    DEBUG_VAR(maxModMass);
    DEBUG_VAR(minModMass);
    DEBUG_VAR(penaltyAlign);
    DEBUG_VAR(penaltyAlpha);
    DEBUG_VAR(penaltyBeta);
    DEBUG_VAR(alignall);
    DEBUG_VAR(maxDbGapAas);
    DEBUG_VAR(maxSpectrumGapDaltons);
    
    // -----------------------------------------------------
    // Setup mods for Spectral Probability Computation
    // -----------------------------------------------------
    if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(m_penaltyMatrixMods);
    vector<pair<unsigned int, bool> > ntermMods;
    vector<unsigned int> mods;
    m_penaltyMatrixMods->getSpecProbMods(ntermMods, mods);

    AAJumps aminoacids(1); // Used to calculate spectral probabilities
    bool isReversed;

    // -----------------------------------------------------
    // Compute Spectral Probability scores
    // -----------------------------------------------------
    clock_t startTime = clock();

    if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(psmSet->size());

    for (int i = 0; i < psmSet->size(); i++) {
      
      int scanNumber = (*psmSet)[i]->m_scanNum;
      if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(scanNumber);

      if ( m_params.exists("DEBUG_SPECPROTALIGN_SINGLESPECTRUM")) {
        int scanDebug = m_params.getValueInt("DEBUG_SPECPROTALIGN_SINGLESPECTRUM");
        if (scanDebug != scanNumber) {
          continue;
        }
      }

      Spectrum & origSpec = *((*psmSet)[i]->m_spectrum);

      if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_TRACE;
      if (DEBUG_SPECPROTALIGN_SPECPROB) origSpec.outputDebug();
      if (DEBUG_SPECPROTALIGN_SPECPROB) 
          DEBUG_MSG("OR  " << (*psmSet)[i]->m_annotation.c_str() << "  " 
                   << (*psmSet)[i]->m_score << "  " << (*psmSet)[i]->m_pValue);

      if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(origSpec.psmList.size());
      if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(m_prmSpectra);
      if (DEBUG_SPECPROTALIGN_SPECPROB && m_prmSpectra != 0x0) DEBUG_VAR(m_prmSpectra->size());

      if (!penaltyAlign) {

      } else if (m_prmSpectra != 0x0 && m_prmSpectra->size() != 0) {

        if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR((*psmSet)[i]->m_annotation);

        if (DEBUG_SPECPROTALIGN_COMPUTE_GAP_ANNOS) {
          string fullAnnotation;
          apb->computeAllGapAnnotations((*psmSet)[i]->m_annotation, fullAnnotation);
          if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(fullAnnotation);
          (*psmSet)[i]->m_annotation = fullAnnotation;
        }

        if ((*psmSet)[i]->m_score < 0) {
          (*psmSet)[i]->m_pValue = 10.0;
          vector<float> modifications;
          (*psmSet)[i]->getModifications(modifications);
          (*psmSet)[i]->m_numMods = modifications.size();
          if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_MSG("NG  " << (*psmSet)[i]->m_annotation.c_str()
                                                             << "  " << (*psmSet)[i]->m_score
                                                             << "  " << (*psmSet)[i]->m_pValue);
          continue;
        }

        vector<float> modifications;
        vector<unsigned int> positions;
        vector<unsigned int> lengths;
        (*psmSet)[i]->getModificationsAndPositions(modifications, positions, lengths);
        for (int j = 0; j < modifications.size(); j++) {
          if (DEBUG_SPECPROTALIGN_SPECPROB) {
            DEBUG_MSG(j << "  " << modifications[j] << "  " << positions[j] << "  " << lengths[j]);
          }
        }
        float parentMassShift = 0.0;
        if (lengths.size() > 0 &&
            positions[0] == 0 &&
            lengths[0] == 0 &&
            modifications[0] >= -pmTol) {
          DEBUG_MSG(modifications[0] << "  " << positions[0] << "  " << lengths[0]);
          parentMassShift = modifications[0];
        } else if (lengths.size() > 0 &&
            positions[0] == 0 &&
            lengths[0] == 0 &&
            modifications[0] <= pmTol) {
          DEBUG_MSG(modifications[0] << "  " << positions[0] << "  " << lengths[0]);
          parentMassShift = modifications[0];
        }
        if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(parentMassShift);

        string cleanAnnotation;
        PeptideSpectrumMatch::getUnmodifiedPeptide((*psmSet)[i]->m_annotation, cleanAnnotation);
        if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(cleanAnnotation.size());

        float parentMassShift2 = 0.0;
        if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(lengths.size());
        if (lengths.size() > 0 &&
            positions[lengths.size()-1] == cleanAnnotation.size() &&
            lengths[lengths.size()-1] == 0 &&
            modifications[lengths.size()-1] >= -pmTol && modifications[lengths.size()-1] <= 0) {
          DEBUG_MSG(modifications[lengths.size()-1] << "  " << positions[lengths.size()-1] << "  " << lengths[lengths.size()-1]);
          parentMassShift2 = modifications[lengths.size()-1];
        } else if (lengths.size() > 0 &&
            positions[lengths.size()-1] == cleanAnnotation.size() &&
            lengths[lengths.size()-1] == 0 &&
            modifications[lengths.size()-1] >= 0 && modifications[lengths.size()-1] <= pmTol) {
          DEBUG_MSG(modifications[lengths.size()-1] << "  " << positions[lengths.size()-1] << "  " << lengths[lengths.size()-1]);
          parentMassShift2 = modifications[lengths.size()-1];
        }
        if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(parentMassShift2);

        string annotationNoParentMassAdjusts;
        stripParentMassAdjustments((*psmSet)[i], annotationNoParentMassAdjusts);
        if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(annotationNoParentMassAdjusts);
        (*psmSet)[i]->m_annotation = annotationNoParentMassAdjusts;

        Spectrum tempSprinkled = origSpec;

        tempSprinkled.psmList.clear();
        tempSprinkled.psmList.push_back((*psmSet)[i]);

        tempSprinkled.parentMass -= parentMassShift; // reverse the nterm shift
        tempSprinkled.parentMass -= parentMassShift2; // reverse the cterm shift
        if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_VAR(tempSprinkled.parentMass);

        if (parentMassShift == 0.0) {
          tempSprinkled.computeSpectralProbabilities(ntermMods, mods, peakTol, useRawScoreForProb);
          if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_MSG("SP  " << (*psmSet)[i]->m_annotation.c_str()
                                                             << "  " << (*psmSet)[i]->m_score
                                                             << "  " << (*psmSet)[i]->m_pValue);
        } else if (parentMassShift > 0.0) {
          // Shift all the peaks by -1 to counteract the +1
          for (int iPeak = 1; iPeak < tempSprinkled.size(); iPeak++) {
            if (tempSprinkled[iPeak ][0] > 1.0 ) {
              tempSprinkled[iPeak ][0] -= parentMassShift;
            }
            //if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_MSG(i << "  " << tempSprinkled[iPeak ][0]);
          }

          tempSprinkled.computeSpectralProbabilities(ntermMods, mods, peakTol, useRawScoreForProb);
          if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_MSG("M1  " << (*psmSet)[i]->m_annotation.c_str()
                                                             << "  " << (*psmSet)[i]->m_score
                                                             << "  " << (*psmSet)[i]->m_pValue);
        } else if (parentMassShift < 0.0) {
          // Shift all the peaks by 1 to counteract the +1
          for (int iPeak = 1; iPeak < tempSprinkled.size(); iPeak++) {
            if (tempSprinkled[iPeak ][0] > 0.0 ) {
              tempSprinkled[iPeak ][0] += parentMassShift;
            }
          }
          tempSprinkled.computeSpectralProbabilities(ntermMods, mods, peakTol, useRawScoreForProb);
          if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_MSG("P1  " << (*psmSet)[i]->m_annotation.c_str()
                                                             << "  " << (*psmSet)[i]->m_score
                                                             << "  " << (*psmSet)[i]->m_pValue);
        }
        // LARS: Temp fix to "0" pvalues
        if ((*psmSet)[i]->m_pValue == 1e-38) {
          (*psmSet)[i]->m_pValue = 20.0;
        }

      } else {

        if (DEBUG_SPECPROTALIGN_COMPUTE_GAP_ANNOS) {
          string fullAnnotation;
          apb->computeAllGapAnnotations((*psmSet)[i]->m_annotation, fullAnnotation);
          (*psmSet)[i]->m_annotation = fullAnnotation;
        }

        Spectrum temp = origSpec;
        (*psmSet)[i]->m_spectrum = &temp;
        temp.psmList.clear();
        temp.psmList.push_back((*psmSet)[i]);
        temp.computeSpectralProbabilities(ntermMods, mods, peakTol, useRawScoreForProb);

        if (DEBUG_SPECPROTALIGN_SPECPROB) DEBUG_MSG("UN  " << (*psmSet)[i]->m_annotation.c_str() << "  " << (*psmSet)[i]->m_score << "  " << (*psmSet)[i]->m_pValue);
      }

      vector<float> modifications;
      (*psmSet)[i]->getModifications(modifications);
      (*psmSet)[i]->m_numMods = modifications.size();
    }

    float elapsed = ((float)(clock() - startTime)) / CLOCKS_PER_SEC;
    if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsed << " seconds elapsed for spectral probs");

    return;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::createTags(void)
  {
    m_inputSpectra->clearPsms(); // Clear all the target PSMs
    PeptideSpectrumMatchSet psmSetTag;
    ExecTagSearch moduleTagsearch(m_params,
                                  m_inputSpectra,
                                  m_db,
                                  (vector<unsigned int> *)0, //specsToSearch,
                                  &psmSetTag);
    moduleTagsearch.invoke();
    DEBUG_VAR(psmSetTag.size());
    psmSetTag.addSpectra(m_inputSpectra); // Associate PSMs
    if(m_params.exists("OUTPUT_TAG_PSM")) {
      string outDir = m_params.getValue("GRID_DATA_DIR_OUT");
      if (outDir.empty()) {
        outDir = ".";
      }
      string saveName = outDir + "/" + m_params.getValue("OUTPUT_TAG_PSM").c_str() + ".original";
      DEBUG_MSG("Saving tag file [" << saveName << "]");
      psmSetTag.saveToFile(saveName.c_str(), true);
    }
    
    return;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::createSeedTags(void)
  {
    if (!m_filterPsmSet || !m_contigAbinfo) {
      return;
    }
      
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(m_filterPsmSet->size());
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(m_contigAbinfo->size());
    if (m_filterPsmSet->size() == 0) {
      return;
    }

    // Create the Assembled Shifts structure from the contig data and specs    
    vector<vector<sps::tuple<unsigned int, float, bool> > > spectrumLists;
    getAssembledShifts(*m_inputSpectra,
                       *m_contigAbinfo,
                       spectrumLists);

#if 1
    // Create the scan to index number mapping for ease of finding spec
    map<int, int> scanToIndex;               
    for (int i = 0; i < m_inputSpectra->size(); i++) {
      scanToIndex[(*m_inputSpectra)[i].scan] = i;
    }
#endif

    // Create the spec to contig mapping for ease of finding contigs for specs        
    map<int, int> specScanToContigIndex;
    map<int, int> specScanToContigScan;
    map<int, int> scanToOffset;
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(spectrumLists.size());
    for (int i = 0; i < spectrumLists.size(); i++)  {
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_MSG(i << "   " << (*m_inputSpectra)[i].scan);
      for (int j = 0; j < spectrumLists[i].size(); j++)  {
        int spectrumIndex = spectrumLists[i][j].m0;
        int offset = spectrumLists[i][j].m1;
        if (DEBUG_SPECPROTALIGN_TAG_SEEDING && spectrumIndex != 1) // Problem with 1
            DEBUG_MSG("        " << spectrumIndex + 1 << "   " << offset);
        specScanToContigScan[spectrumIndex + 1] = (*m_inputSpectra)[i].scan;
        specScanToContigIndex[spectrumIndex + 1] = i;
        scanToOffset[spectrumIndex + 1] = spectrumLists[i][j].m1;
      }
    }

    // "Normalize" the AA's in the database
    DB_fasta dbCopy;  // Create the DB first
    dbCopy = *m_db;   // Then copy.. if you don't.. you seg fault! (dunno why)
    dbCopy.replaceAA('L', 'I');
    dbCopy.replaceAA('K', 'Q');

    float pmTol = (float)m_params.getValueDouble("TOLERANCE_PM", 0.4);
    float peakTol = (float)m_params.getValueDouble("TOLERANCE_PEAK", 1.5);
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(pmTol);
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(peakTol);

    //------------------------------------------------------------------------
    //  THIS WAS HACK #1 - ELIMINATE THOSE TAGS THE SEEDS REPLACE
    //------------------------------------------------------------------------
#if 0
    // Find the shortest annotation.. that will be our "tag" length
    int shortestAnno = 50;
    for (int i = 0; i < m_filterPsmSet->size(); i++)  {
      int scan = (*m_filterPsmSet)[i]->m_scanNum;
      // Check to see if this scan is part of any contig
      if (specScanToContigScan.find(scan) == specScanToContigScan.end()) {
        continue;
      }
      string anno;
      PeptideSpectrumMatch::getUnmodifiedPeptide((*m_filterPsmSet)[i]->m_annotation, anno);
      //if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) 
      //  DEBUG_MSG(scan << "   " << anno << "   " << specScanToContigScan[scan]);
      if (anno.length() < shortestAnno) {
        shortestAnno = anno.length();
      }
    }
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_VAR(shortestAnno);
    DB_index index(dbCopy, 1024, 1024, shortestAnno);

    //-------------------------------------------------------------------------
    // Find all the dbIndexes according to seed PSMs that contigs belong to
    //-------------------------------------------------------------------------
    map<int, set<int> > mapScanToIndexSet;
    for (int i = 0; i < m_filterPsmSet->size(); i++)  {

      int scan = (*m_filterPsmSet)[i]->m_scanNum;
      // Check to see if this scan is part of any contig
      if (specScanToContigScan.find(scan) == specScanToContigScan.end()) {
        continue;
      }

      string anno;
      PeptideSpectrumMatch::getUnmodifiedPeptide((*m_filterPsmSet)[i]->m_annotation, anno);
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_MSG(anno << "   " << scan << "   c" << specScanToContigScan[scan]);

      // Chop the annotation down to "tag" size
      anno = anno.substr(0, shortestAnno);
      // "Normalize" the AA's in the annotation to match the database
      for (int j = 0; j < shortestAnno; j++) {
        if (anno[j] == 'L') anno[j] = 'I';
        if (anno[j] == 'K') anno[j] = 'Q';
      }
      
      // Find all the tag matches
      list<sps::tuple<int, float, string> > matches;
      index.findAll(anno.c_str(),
                    dbCopy,
                    matches,
                    1,
                    0,
                    pmTol,
                    0,
                    peakTol);
      //if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_VAR(matches.size());

      list<sps::tuple<int, float, string> >::iterator itr = matches.begin();
      list<sps::tuple<int, float, string> >::iterator itrEnd = matches.end();
      for ( ; itr != itrEnd; itr++)  {
        if (DEBUG_SPECPROTALIGN_TAG_SEEDING2)
          DEBUG_MSG(itr->m2 << "  " << itr->m0  << "  " << itr->m1);

        // Add all db indexes of matches to the inclusion set for this scan
        //if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_MSG(specScanToContigScan[scan] << "  " << itr->m0);
        mapScanToIndexSet[specScanToContigScan[scan]].insert(itr->m0);
      }
      
    } // for (int i = 0; i < m_filterPsmSet->size(); i++)  {

    //-------------------------------------------------------------------------
    // Now go through all the PSM lists and remove any that are not
    //   in the inclusion set generated above
    //-------------------------------------------------------------------------
    for (int i = 0; i < m_inputSpectra->size(); i++) {

      int scan = (*m_inputSpectra)[i].scan;
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_VAR(scan);
      // Check to see if there is anything in the set (if not skip it)
      if (mapScanToIndexSet.find(scan) == mapScanToIndexSet.end()) {
        continue;
      }

      if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_VAR((*m_inputSpectra)[i].psmList.size());

      // Make sure we find at least one match to the "right" PSM dbindex
      //   Otherwise we'll end up deleting all the tags
      bool foundOne = false;
      list<psmPtr>::iterator litr = (*m_inputSpectra)[i].psmList.begin();
      list<psmPtr>::iterator litrEnd = (*m_inputSpectra)[i].psmList.end();
      for (; litr != litrEnd; litr++) {
        psmPtr temp = *litr;
        if (mapScanToIndexSet[scan].find(temp->m_dbIndex) != mapScanToIndexSet[scan].end()) {
          foundOne = true;
          break;
        }      
      } // for ( ; itr != itrEnd; itr++)  {

      // Delete any tags that don't have right dbindex according to seed psms
      if (foundOne) {
        for (litr = (*m_inputSpectra)[i].psmList.begin(); litr != litrEnd; ) {
          psmPtr temp = *litr;
          if (mapScanToIndexSet[scan].find(temp->m_dbIndex) == mapScanToIndexSet[scan].end()) {
            litr = (*m_inputSpectra)[i].psmList.erase(litr);
            if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) {
              DEBUG_MSG("ERASING PSM with dbindex [" << temp->m_dbIndex << "]");
            }
          } else {
            ++litr;
            if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) {
              DEBUG_MSG("KEEPING PSM with dbindex [" << temp->m_dbIndex << "]");
            }
          }      
        } // for ( ; itr != itrEnd; itr++)  {
      } // if (foundOne) {
      
    } // for (int i = 0; i < m_inputSpectra.size(); i++) {

    // Overwrite the tag file since we may have changed it
    if(m_params.exists("OUTPUT_TAG_PSM")) {
      PeptideSpectrumMatchSet psmSetTag;
      psmSetTag.getPSMSet(m_inputSpectra);
      string outDir = m_params.getValue("GRID_DATA_DIR_OUT");
      if (outDir.empty()) {
        outDir = ".";
      }
      string saveName = outDir + "/" + m_params.getValue("OUTPUT_TAG_PSM").c_str();
      DEBUG_MSG("Saving tag file [" << saveName << "]");
      psmSetTag.saveToFile(saveName.c_str(), true);
    }
#endif

#if 0
    AAJumps aminoacids(1); // Used to match PSM to spectrum
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_VAR(specIdx);
    bool isReversed;
    float matchScore = MatchSpecToPeptide((*m_inputSpectra[specIdx]),
                                          (*m_filterPsmSet)[i]->m_annotation.c_str(),
                                          peakTol,
                                          0,
                                          false,
                                          &isReversed,
                                          &aminoacids);
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_VAR(matchScore);
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING2) DEBUG_VAR(isReversed);
#endif

    int minSeedTagLength =
        m_params.getValueInt("SPECPROTALIGN_MIN_SEED_TAG_LENGTH", 6);
    if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(minSeedTagLength);
    DB_index indexSeed(dbCopy, 1024, 1024, minSeedTagLength);

    //-------------------------------------------------------------------------
    // Go through the seed PSMs and delete all previous tags from their contigs
    //-------------------------------------------------------------------------
    for (int i = 0; i < m_filterPsmSet->size(); i++)  {

      int scan = (*m_filterPsmSet)[i]->m_scanNum;
      //if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(scan);
      // Check to see if this scan is part of any contig
      if (specScanToContigScan.find(scan) == specScanToContigScan.end()) {
        continue;
      }
      int contigIndex = specScanToContigIndex[scan];
      //if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(contigIndex);
      (*m_inputSpectra)[contigIndex].psmList.clear();
    }

    //------------------------------------------------------------------------
    //  THIS WAS HACK #2 - CREATE TAGS INSIDE THE CONTIGS AT SIZE 8
    //------------------------------------------------------------------------
#if 0
    //-------------------------------------------------------------------------
    // Go through all the spectra and insert tags from the seeds
    //-------------------------------------------------------------------------
    for (int i = 0; i < m_filterPsmSet->size(); i++)  {

      int scan = (*m_filterPsmSet)[i]->m_scanNum;
      //if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(scan);
      // Check to see if this scan is part of any contig
      if (specScanToContigScan.find(scan) == specScanToContigScan.end()) {
        continue;
      }

      string anno;
      PeptideSpectrumMatch::getUnmodifiedPeptide((*m_filterPsmSet)[i]->m_annotation, anno);

      int offset = scanToOffset[scan];
      int contigIndex = specScanToContigIndex[scan];
      float parentMass = (*m_inputSpectra)[contigIndex].parentMass;
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING) {
        DEBUG_MSG(anno << "   " << scan << "   c" <<
                  contigIndex + 1 << "   " << offset << "   " << parentMass);
      }

      int startAA = 0;
      bool foundStart = false;
      int endAA = 0;
      vector<float> peakMasses;
      getMasses((char *)anno.c_str(), peakMasses);
      float totalMass = 0;
      for (int m = 0; m < peakMasses.size(); m++) {
        totalMass += peakMasses[m];
        if (!foundStart && totalMass >= -offset) {
          startAA = m;
          foundStart = true;
        }
        if (totalMass <= parentMass - offset) {
          endAA = m;
        } else {
          break;
        }
      }

      if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_MSG("      " << startAA << "   " << endAA);
      if (endAA - startAA +1 < minSeedTagLength) {
        continue;
      }

      string tag = anno.substr(startAA, endAA - startAA + 1);
      // Chop the tag down to proper size
      tag = tag.substr(0, minSeedTagLength);
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_MSG("      " << tag);

      // "Normalize" the AA's in the annotation to match the database
      for (int j = 0; j < minSeedTagLength; j++) {
        if (tag[j] == 'L') tag[j] = 'I';
        if (tag[j] == 'K') tag[j] = 'Q';
      }

      // Find all the tag matches
      list<sps::tuple<int, float, string> > matches;
      indexSeed.findAll(tag.c_str(),
                    dbCopy,
                    matches,
                    1,
                    0,
                    pmTol,
                    0,
                    peakTol);
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_MSG("      matches = " << matches.size());

      list<sps::tuple<int, float, string> >::iterator itr = matches.begin();
      list<sps::tuple<int, float, string> >::iterator itrEnd = matches.end();
      for ( ; itr != itrEnd; itr++)  {
        if (DEBUG_SPECPROTALIGN_TAG_SEEDING)
          DEBUG_MSG("      " << itr->m2 << "  " << itr->m0  << "  " << itr->m1);

        psmPtr p(new PeptideSpectrumMatch);
        p->m_scanNum = scan;
        p->m_protein = string(m_db->getID(itr->m0));
        p->m_annotation = itr->m2;
        p->m_origAnnotation = tag;
        p->m_dbIndex = itr->m0;
        p->m_startMass = itr->m1;

        (*m_inputSpectra)[contigIndex].psmList.push_back(p);

        // Add all db indexes of matches to the inclusion set for this scan
        //if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_MSG(specScanToContigScan[scan] << "  " << itr->m0);
        //mapScanToIndexSet[specScanToContigScan[scan]].insert(itr->m0);
      }
    }  // for (int i = 0; i < m_filterPsmSet->size(); i++)  {
#else
    //-------------------------------------------------------------------------
    // Go through all the spectra and insert tags from the seeds
    //-------------------------------------------------------------------------
    for (int i = 0; i < m_filterPsmSet->size(); i++)  {

      int scan = (*m_filterPsmSet)[i]->m_scanNum;
      //if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_VAR(scan);
      // Check to see if this scan is part of any contig
      if (specScanToContigScan.find(scan) == specScanToContigScan.end()) {
        continue;
      }

      string anno;
      PeptideSpectrumMatch::getUnmodifiedPeptide((*m_filterPsmSet)[i]->m_annotation, anno);

      int offset = scanToOffset[scan];
      int contigIndex = specScanToContigIndex[scan];
      float parentMass = (*m_inputSpectra)[contigIndex].parentMass;
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING) {
        DEBUG_MSG(anno << "   " << scan << "   c" <<
                  contigIndex + 1 << "   " << offset << "   " << parentMass);
      }
      
      if (anno.size() < minSeedTagLength) {
        continue;
      }

      // Chop the tag down to proper size
      string tag = anno.substr(0, minSeedTagLength);
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_MSG("      " << tag);

      // "Normalize" the AA's in the annotation to match the database
      for (int j = 0; j < minSeedTagLength; j++) {
        if (tag[j] == 'L') tag[j] = 'I';
        if (tag[j] == 'K') tag[j] = 'Q';
      }

      // Find all the tag matches
      list<sps::tuple<int, float, string> > matches;
      indexSeed.findAll(tag.c_str(),
                    dbCopy,
                    matches,
                    1,
                    0,
                    pmTol,
                    0,
                    peakTol);
      if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_MSG("      matches = " << matches.size());

      list<sps::tuple<int, float, string> >::iterator itr = matches.begin();
      list<sps::tuple<int, float, string> >::iterator itrEnd = matches.end();
      for ( ; itr != itrEnd; itr++)  {
        if (DEBUG_SPECPROTALIGN_TAG_SEEDING)
          DEBUG_MSG("      " << itr->m2 << "  " << itr->m0  << "  " << itr->m1);

        psmPtr p(new PeptideSpectrumMatch);
        p->m_scanNum = scan;
        p->m_protein = string(m_db->getID(itr->m0));
        p->m_annotation = itr->m2;
        p->m_origAnnotation = tag;
        p->m_dbIndex = itr->m0;
        // Add in the offset to the tag match location so it lines up with contig start
        p->m_startMass = itr->m1 - offset;

        (*m_inputSpectra)[contigIndex].psmList.push_back(p);

        // Add all db indexes of matches to the inclusion set for this scan
        //if (DEBUG_SPECPROTALIGN_TAG_SEEDING) DEBUG_MSG(specScanToContigScan[scan] << "  " << itr->m0);
        //mapScanToIndexSet[specScanToContigScan[scan]].insert(itr->m0);
      }
    }  // for (int i = 0; i < m_filterPsmSet->size(); i++)  {
#endif

    // Overwrite the tag file since we may have changed it
    if(m_params.exists("OUTPUT_TAG_PSM")) {
      PeptideSpectrumMatchSet psmSetTag;
      psmSetTag.getPSMSet(m_inputSpectra);
      string outDir = m_params.getValue("GRID_DATA_DIR_OUT");
      if (outDir.empty()) {
        outDir = ".";
      }
      string saveName = outDir + "/" + m_params.getValue("OUTPUT_TAG_PSM").c_str();
      DEBUG_MSG("Saving tag file [" << saveName << "]");
      psmSetTag.saveToFile(saveName.c_str(), true);
    } else {
      PeptideSpectrumMatchSet psmSetTag;
      psmSetTag.getPSMSet(m_inputSpectra);
      psmSetTag.saveToFile("aoutput_tag.psm", true);
    }

    return;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecProtAlign::prepareSpectrum(Spectrum & cSpec,
                                          Spectrum & cSpecRev, 
                                          float    & avgPeakIntensity,
                                          bool scaleIntensities, 
                                          bool sprinkle,
                                          bool fullSprinkle)
  {
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(scaleIntensities);
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(sprinkle);

    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Input spectrum");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();

    // Make sure the peak tolerances are set properly
    float peakTol = (float)m_params.getValueDouble("TOLERANCE_PEAK");
    DEBUG_VAR(peakTol);

    // Make sure the peak tolerances is at least 1.0 since that is our "resolution"
    peakTol = max(peakTol, (float)1.0);
    DEBUG_VAR(peakTol);

    cSpec.setPeakTolerance(peakTol);

    // Clean up any pre-existing endpoints so we don't add too many
    list<int> peaksToRemove;
    for (int iPeak = 0; iPeak < cSpec.size(); iPeak++) {
      if (cSpec[iPeak][0] < 57.0 || cSpec[iPeak][0] > cSpec.parentMass - 57.0) {
        peaksToRemove.push_back(iPeak);
      }
    }
    if (peaksToRemove.size() != 0) {
      cSpec.removePeaks(peaksToRemove);
    }
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Endpoints cleaned up");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();

    // If there are no useful peaks then this spectrum is garbage
    if (cSpec.size() == 0) {
      WARN_MSG("Spectrum [" << cSpec.scan << "] contains no usable peaks!");
      return false;
    }

    if (sprinkle) {
      // If we are scaling intensities, then use full sprinkle
      sprinkleSpectrum(cSpec, false, fullSprinkle);
      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Sprinkled");
      if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();
    }

    // Get the reverse spectrum
    cSpec.reverse(0, &cSpecRev);
    cSpecRev.psmList.clear();

    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Reversed");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpecRev.outputDebug();

    //bool enforceEndpeaks = m_params.getValueInt("ENFORCE_ENDPEAKS");
    //DEBUG_VAR(enforceEndpeaks);

    // Round all the spectra to "center them at 0.0"
    cSpec.roundPeaks(AA_ROUNDING, false, false);
    cSpecRev.roundPeaks(AA_ROUNDING, false, false);
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Rounded");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Rounded (reverse)");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpecRev.outputDebug();

    cSpec.setResolution(1, true); // Merge peaks with the same rounded mass
    cSpecRev.setResolution(1, true); // Merge peaks with the same rounded mass

    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Resolution set");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Resolution set (reverse)");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpecRev.outputDebug();

    if (scaleIntensities) {
      scalePeakIntensities(cSpec, avgPeakIntensity);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeakIntensity)
      float dummy;
      scalePeakIntensities(cSpecRev, dummy);
    }
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeakIntensity)

    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Scaled");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Scaled (reverse)");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpecRev.outputDebug();

    string avgType = m_params.getValue("SPECPROTALIGN_AVERAGE_TYPE");
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgType);

    if (avgType == "topk") {
      //int expectNumberOfPeaks = (int)(cSpec.parentMass / AVERAGE_AA_MASS * 2.0);
      avgPeakIntensity = averageTopPeaks(cSpec, 20);
    } else if (avgType == "avgaboveavg") {
      float avgPeak = 0.0;
      for (int i = 0; i < cSpec.size(); i++) {
        avgPeak += cSpec[i][1];
      }
      avgPeak /= (float)cSpec.size();
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeak);
      int peaksAboveAvg = 0;
      for (int i = 0; i < cSpec.size(); i++) {
        if (cSpec[i][1] > avgPeak) {
          peaksAboveAvg++;
          avgPeakIntensity += cSpec[i][1];
        }
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(peaksAboveAvg);
      avgPeakIntensity /= (float)peaksAboveAvg;
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeakIntensity);
    } else {
      for (int i = 0; i < cSpec.size(); i++) {
        avgPeakIntensity += cSpec[i][1];
      }
      avgPeakIntensity /= (float)cSpec.size();
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeakIntensity);
    }
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeakIntensity);

    cSpec.addZPMpeaks(peakTol, 0.0, false);
    cSpecRev.addZPMpeaks(peakTol, 0.0, false);

    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("With endpeaks");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("With endpeaks (reverse)");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpecRev.outputDebug();

    return true;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::recomputeWithSprinkled(
                                    PeptideSpectrumMatchSet & psmSetTemp2,
                                    int specIdx,
                                    AlignmentPenaltyBased * apb)
  {
    int minNumMatchPeaks = m_params.getValueInt("MIN_MATCHED_PEAKS_DB");
    int maxDbGapAas = m_params.getValueInt("MAX_ALIGN_DB_GAP_AAS", 8);
    float pmTol = (float)m_params.getValueDouble("TOLERANCE_PM");
    float peakTol = (float)m_params.getValueDouble("TOLERANCE_PEAK");
    bool enforceEndpeaks = m_params.getValueInt("ENFORCE_ENDPEAKS");
    
    SpecSet & inputSpectra = *m_matchedSpectraAll; // for convenience
    for (int i = 0; i < psmSetTemp2.size(); i++) {

      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("-------------------------------------------------");
      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("SCAN NUMBER = " << psmSetTemp2[i]->m_scanNum);
      if (DEBUG_SPECPROTALIGN) DebugPsm(psmSetTemp2[i]);

      Spectrum cSpec;
      int matchingIdx = 0;
      for (specIdx = 0; specIdx < inputSpectra.size(); specIdx++) {
        int scanNum = inputSpectra[specIdx].scan;
        if (inputSpectra[specIdx].scan == psmSetTemp2[i]->m_scanNum) {
          cSpec = inputSpectra[specIdx];
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[specIdx].psmList.size());
          matchingIdx = specIdx;
          break;
        }
      }
      if (cSpec.size() == 0) {
        ERROR_MSG("Couldn't find original spectrum!");
        return;
      }
      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Original Spectrum");
      if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();

      Spectrum cSpecRev;
      float avgPeakIntensity = 0;
      prepareSpectrum(cSpec, cSpecRev, avgPeakIntensity, true, true, true);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeakIntensity)

      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(psmSetTemp2[i]->m_matchOrientation);
      if (psmSetTemp2[i]->m_matchOrientation) {
        cSpec = cSpecRev;
      }

      // Get the starting positions for this protein and the correct orientation
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(matchingIdx);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[matchingIdx].psmList.size());
      set<float> startPositionList;
      list<psmPtr>::iterator itr = inputSpectra[matchingIdx].psmList.begin();
      list<psmPtr>::iterator itrEnd = inputSpectra[matchingIdx].psmList.end();
      for (; itr != itrEnd; itr++) {
        int idx = (*itr)->m_dbIndex;
        if (psmSetTemp2[i]->m_matchOrientation == 0 && (*itr)->m_matchOrientation == 0) {
          startPositionList.insert((*itr)->m_startMass);
          if (DEBUG_SPECPROTALIGN) DEBUG_MSG(idx << "  " << (*itr)->m_startMass)
        } else if (psmSetTemp2[i]->m_matchOrientation == 1 && (*itr)->m_matchOrientation == 1) {
          startPositionList.insert((*itr)->m_startMass);
          if (DEBUG_SPECPROTALIGN) DEBUG_MSG(idx << "  " << (*itr)->m_startMass)
        }
      } // for (; itr != itrEnd; itr++)
      
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(startPositionList.size());

      cSpec.psmList.clear();

      // Replace the spectrum with our new one for spectral prob
      *(psmSetTemp2[i]->m_spectrum) = cSpec;

      int protIdx = psmSetTemp2[i]->m_dbIndex;
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(protIdx);

      Spectrum dbSpec;
      dbSpec.reserve(1024);
      dbSpec = m_db->getMassesSpec(protIdx);
      // Round all the db spectra to "center them at 0.0"
      dbSpec.roundPeaks(AA_ROUNDING, false, false);
      
      if (DEBUG_SPECPROTALIGN_SPECTRA) dbSpec.outputDebug();

      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(psmSetTemp2[i]->m_matchOrientation);
      apb->computeAlignment(cSpec,
                             dbSpec,
                             m_db->getSequence(protIdx),
                             protIdx,
                             psmSetTemp2[i]->m_matchOrientation,
                             startPositionList,
                             (*m_scanSpecificPenalties)[psmSetTemp2[i]->m_scanNum],
                             minNumMatchPeaks,
                             maxDbGapAas,
                             pmTol,
                             peakTol,
                             avgPeakIntensity,
                             enforceEndpeaks);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
      if (cSpec.psmList.size() == 0) {
        ERROR_MSG("No PSM came back from second alignment computation!");
        return;
      }
      // Now replace the PSM from before with our new one
      psmPtr tempPsm = *(cSpec.psmList.begin());
      psmSetTemp2[i]->m_origAnnotation = tempPsm->m_annotation;
      psmSetTemp2[i]->m_annotation = tempPsm->m_annotation;
      psmSetTemp2[i]->m_score = tempPsm->m_score;
      psmSetTemp2[i]->m_startMass = tempPsm->m_startMass;

      // Make sure that the spectrum we are pointing to also has this PSM in its list
      psmSetTemp2[i]->m_spectrum->psmList.push_back(psmSetTemp2[i]);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(psmSetTemp2[i]->m_spectrum->psmList.size());
      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("-------------------------------------------------");
    }
    
    return;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::scalePeakIntensities(Spectrum & cSpec, float & avgPeakIntensity)
  {
    // Need to scale the intensities for compatibility with spectral prob
    float maxScore = 0;
    for (unsigned int i = 0; i < cSpec.size(); i++) {
      maxScore = maxScore < cSpec[i][1] ? cSpec[i][1] : maxScore;
    }
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(maxScore);
    for (unsigned int i = 0; i < cSpec.size(); i++) {
      float scaledPeak = ceil(cSpec[i][1] / maxScore * MAX_SPEC_PROB_PEAK);
      cSpec[i][1] = scaledPeak;
    }

    // Scale the average also
    avgPeakIntensity = avgPeakIntensity / maxScore * MAX_SPEC_PROB_PEAK;
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(avgPeakIntensity)

    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Scaled Intensities");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();

#if 0
    // Remove those peaks whose intensities are now 0
    Spectrum tempSpec = cSpec;
    unsigned int j = 0;
    for (unsigned int i = 0; i < cSpec.size(); i++) {
      if (cSpec[i][1] > 0.00001 && j != i) {
        tempSpec[j++] = cSpec[i];
      }
    }
    if (j == 0) {
      ERROR_MSG("All peaks removed from spectrum!");
      return;
    }
    cSpec = tempSpec;
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(j-1);
    cSpec.resize(j-1);
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(cSpec.size());
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Culled 0 intensities");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();
#endif

    return;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::sprinkleSpectrum(Spectrum & cSpec, bool reversed, bool useFullSprinkle)
  {
    float prmMultiplier = m_params.getValueFloat("SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER", 0.1);
    DEBUG_VAR(prmMultiplier);
    bool useSprinkling = (prmMultiplier != 0.0);
    DEBUG_VAR(useSprinkling);
    if (!useSprinkling) {
      return;
    }

    int rankNumber = m_params.getValueInt("SPECPROTALIGN_PRM_CONTRIBUTION_RANK", 3);
    DEBUG_VAR(rankNumber);
    DEBUG_VAR(useFullSprinkle);

    map<int, int> scanToIndex;
    DEBUG_VAR(m_prmSpectra);
    DEBUG_VAR(m_prmSpectra->size());
    for (int i = 0; i < m_prmSpectra->size(); i++) {
      int scanNum = (*m_prmSpectra)[i].scan;
      scanToIndex[scanNum] = i;
    }

    Spectrum sprinkledSpec;
    int scanNum = cSpec.scan;
    DEBUG_VAR(scanNum - 1);
    Spectrum tempSprinkled = (*m_prmSpectra)[scanToIndex[scanNum]];
    tempSprinkled.setResolution(1, true); // Merge peaks with the same rounded mass

    if (!useFullSprinkle) {
      tempSprinkled.rankFilterPeaks(rankNumber);
    }

    DEBUG_VAR(reversed);
    if (reversed) {
      Spectrum tempSprinkledRev;
      tempSprinkled.reverse(0, &tempSprinkledRev);
      cSpec.addPeaksPartial(tempSprinkledRev, sprinkledSpec, prmMultiplier);
    } else {
      cSpec.addPeaksPartial(tempSprinkled, sprinkledSpec, prmMultiplier);
    }
    cSpec = sprinkledSpec;

    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Sprinkled");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.outputDebug();

    return;

  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::stripParentMassAdjustments(psmPtr & psm, string & annotationOut)
  {
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_TRACE;
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(psm->m_annotation);
    string cleanAnnotation;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psm->m_annotation, cleanAnnotation);
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(cleanAnnotation);

    float pmTol = (float)m_params.getValueDouble("TOLERANCE_PM", 0.4);
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(pmTol);

    float minAdjustment = m_params.getValueInt("MIN_PARENT_MASS_ADUSTMENT", 1.0);
    float pmAdjustment = max(pmTol, minAdjustment);
 
      vector<float> modifications;
      vector<unsigned int> positions;
      vector<unsigned int> lengths;
      psm->getModificationsAndPositions(modifications, positions, lengths);

      if (DEBUG_SPECPROTALIGN_ANNO){
        for (int i = 0; i < modifications.size(); i++) {
          DEBUG_MSG(i << "  " << modifications[i] << "  " << positions[i] << "  " << lengths[i]);
        }
      }

      vector<float> modificationsNew;
      vector<unsigned int> positionsNew;
      vector<unsigned int> lengthsNew;
      for (int i = 0; i < modifications.size(); i++) {
        int modSize = (int)(fabs(modifications[i]) + 0.5);
        if (i == 0 && lengths[i] == 0 && modSize <= pmAdjustment) {
          // Remove the nterm parent mass adjustments
#if 0
          if (modifications.size() > 1 && lengths[i + 1] > 1 && 
              positions[i + 1] - lengths[i + 1] == 0) {
            // merge into the next modification
            modifications[i+1] += modifications[i];
          }
#endif
        } else if (i == modifications.size() - 1 && 
                   lengths[i] == 0 && 
                   modSize <= pmAdjustment) {
          // Remove the cterm parent mass adjustments
#if 0
          if (i >= 1 && positions[i - 1] == cleanAnnotation.size() && lengths[i - 1] > 1) {
            // merge into the previous modification (prev mod already in vector)
            modificationsNew[modificationsNew.size() - 1] += modifications[i]; 
          }
#endif
        } else {
          modificationsNew.push_back(modifications[i]);
          positionsNew.push_back(positions[i]);
          lengthsNew.push_back(lengths[i]);
        }
      }

      PeptideSpectrumMatch::makeAnnotationFromData(cleanAnnotation, 
                             modificationsNew, 
                             positionsNew, 
                             lengthsNew, 
                             annotationOut);
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(annotationOut);
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_TRACE;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::mergeModifications(psmPtr & psm, string & annotationOut)
  {
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_TRACE;
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(psm->m_annotation);
    string cleanAnnotation;
    PeptideSpectrumMatch::getUnmodifiedPeptide(psm->m_annotation, cleanAnnotation);
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(cleanAnnotation);

    vector<float> modifications;
    vector<unsigned int> positions;
    vector<unsigned int> lengths;
    psm->getModificationsAndPositions(modifications, positions, lengths);
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(modifications.size());
    if (modifications.size() < 2) {
      annotationOut = psm->m_annotation;
      return;
    }

    if (DEBUG_SPECPROTALIGN_ANNO) {
      for (int i = 0; i < modifications.size(); i++) {
        DEBUG_MSG(i << "  " << modifications[i] << "  " << positions[i] << "  " << lengths[i]);
      }
    }

    vector<float> modificationsNew;
    vector<unsigned int> positionsNew;
    vector<unsigned int> lengthsNew;

    for (int i = 0; i < modifications.size(); i++) {
      int modSize = (int)(fabs(modifications[i]) + 0.5);

      if (DEBUG_SPECPROTALIGN_ANNO) {
        if (i < modifications.size() - 1) {
          DEBUG_MSG("     " << i << "  " << positions[i] << "  " << positions[i + 1] << "  "
                            << modifications[i] << "  " << modifications[i + 1] << "  "
                            << fabs(modifications[i] + modifications[i+1]));
        }
      }

      modifications[i] = roundMass(modifications[i]);

      if (i < modifications.size() - 1) {
        string strAA = cleanAnnotation.substr(positions[i] - lengths[i], lengths[i]);
        string strAA2 = cleanAnnotation.substr(positions[i + 1] - lengths[i + 1], lengths[i + 1]);
        if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(strAA);
        if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(strAA2);
#if 0
        // Removed this merge because we now only allow nterm + known
        if (positions[i] == 0 &&         // Nterm gap
            lengths[i + 1] > 1 &&        // Must be gap
            positions[i + 1] - lengths[i + 1] == 0) {
          // Merge Nterm gap into following gap mod since we can't tell
          modifications[i+1] += modifications[i];
          if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_MSG("Merged Nterm");
        } else
#endif
          if (!m_penaltyMatrixMods->isKnown(strAA, modifications[i]) &&
                   !m_penaltyMatrixMods->isKnown(strAA2, modifications[i+1]) && 
                   positions[i] == positions[i + 1] - lengths[i + 1] &&
                   fabs(modifications[i] + modifications[i+1]) < 1.5) {
          // The mods must both be > 1 in length
          // The mods must be next to each other
          // The mods must add up to < 1.5
          // merge into the next modification
          modifications[i+1] += modifications[i];
          lengths[i+1] += lengths[i];
          if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_MSG("Merged");
        } else if (fabs(modifications[i]) > 0.5) {
          modificationsNew.push_back(modifications[i]);
          positionsNew.push_back(positions[i]);
          lengthsNew.push_back(lengths[i]);
        }
      } else if (fabs(modifications[i]) > 0.5) {
          modificationsNew.push_back(modifications[i]);
          positionsNew.push_back(positions[i]);
          lengthsNew.push_back(lengths[i]);
      }
    }

    PeptideSpectrumMatch::makeAnnotationFromData(cleanAnnotation,
                           modificationsNew,
                           positionsNew,
                           lengthsNew,
                           annotationOut);
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(annotationOut);
    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_TRACE;
  }

  // -------------------------------------------------------------------------
  bool ExecSpecProtAlign::validateParams(std::string & error)
  {
    m_isValid = false;

#if 0
    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    VALIDATE_PARAM_EXIST("ENFORCE_ENDPEAKS");
    VALIDATE_PARAM_EXIST("MAX_MOD_MASS");
    VALIDATE_PARAM_EXIST("MIN_MOD_MASS");
    VALIDATE_PARAM_EXIST("MAX_ALIGN_DB_GAP_AAS");
    VALIDATE_PARAM_EXIST("MAX_ALIGN_SPECTRUM_GAP_DALTONS");
    VALIDATE_PARAM_EXIST("MAX_NUM_MODS");
    VALIDATE_PARAM_EXIST("MIN_MATCHED_PEAKS_DB");
    VALIDATE_PARAM_EXIST("MAX_PARSIMONY");
#endif

    m_isValid = true;
    return true;
  }

}
