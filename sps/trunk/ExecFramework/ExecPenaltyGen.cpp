#include "ExecPenaltyGen.h"

#include "DelimitedTextReader.h"
#include "Logger.h"
#include "ParameterList.h"
#include "PenaltyMatrix.h"
#include "PeptideSpectrumMatchSet.h"
#include "SpectrumPairSet.h"
#include "spectrum.h"

static bool DEBUG_SCAN_SPECIFIC = false;

using namespace std;
using namespace specnets;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecPenaltyGen::ExecPenaltyGen(void) :
    m_filteredPairs(0x0), 
    m_prmSpectra(0x0), 
    m_penaltyMatrixMods(0x0),
    m_scanSpecificPenalties(0x0),
    ownInput(true), ownOutput(true)
  {
    DEBUG_TRACE
    m_name = "ExecPenaltyGen";
    m_type = "ExecPenaltyGen";
  }

  // -------------------------------------------------------------------------
  ExecPenaltyGen::ExecPenaltyGen(const ParameterList & inputParams) :
    ExecBase(inputParams), 
        m_filteredPairs(0x0), 
        m_prmSpectra(0x0), 
        m_penaltyMatrixMods(0x0),
        m_scanSpecificPenalties(0x0),
    ownInput(true), ownOutput(true)
  {
    DEBUG_TRACE
    m_name = "ExecPenaltyGen";
    m_type = "ExecPenaltyGen";
  }

  // -------------------------------------------------------------------------
  ExecPenaltyGen::ExecPenaltyGen(const ParameterList & inputParams,
                                 SpectrumPairSet *     filteredPairs,
                                 SpecSet *             prmSpectra,
                                 PenaltyMatrix *       penaltyMatrixMods,
                                 map<int, map<int, float> > * scanSpecificPenalties) :
    ExecBase(inputParams), 
        m_filteredPairs(filteredPairs), 
        m_prmSpectra(prmSpectra), 
        m_penaltyMatrixMods(penaltyMatrixMods),
        m_scanSpecificPenalties(scanSpecificPenalties),
        ownInput(false), ownOutput(false)
  {
    DEBUG_TRACE
    m_name = "ExecPenaltyGen";
    m_type = "ExecPenaltyGen";
  }

  // -------------------------------------------------------------------------
  ExecPenaltyGen::~ExecPenaltyGen(void)
  {
    if (ownInput) {
      if (m_filteredPairs) {
        delete m_filteredPairs;
      }
      if (m_prmSpectra) {
        delete m_prmSpectra;
      }
    }
    if (ownOutput) {
      if (m_penaltyMatrixMods) {
        delete m_penaltyMatrixMods;
      }
      if (m_scanSpecificPenalties) {
        delete m_scanSpecificPenalties;
      }
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecPenaltyGen::clone(const ParameterList & inputParams) const
  {
    return new ExecPenaltyGen(inputParams);
  }
  // -------------------------------------------------------------------------
  bool ExecPenaltyGen::invoke(void)
  {
    DEBUG_TRACE;
    float resolution = m_params.getValueDouble("ALIGNMENT_RESOLUTION", 1.0);
    DEBUG_VAR(resolution);
    float maxPeakEquivalents = m_params.getValueDouble("MAX_PEAK_EQUIVALENTS", 1.5);
    DEBUG_VAR(maxPeakEquivalents);
    float minPeakEquivalents = m_params.getValueDouble("MIN_PEAK_EQUIVALENTS", 1.0);
    DEBUG_VAR(minPeakEquivalents);

    float minModMass = m_params.getValueDouble("MIN_MOD_MASS", -100.0);
    DEBUG_VAR(minModMass);
    float maxModMass = m_params.getValueDouble("MAX_MOD_MASS", 100.0);
    DEBUG_VAR(maxModMass);
    float maxDiffMass = max(-minModMass, maxModMass);
    DEBUG_VAR(maxDiffMass);

    map<float, float> modFreqs;
    m_filteredPairs->getModificationFrequencies(resolution, maxDiffMass, modFreqs);

    float minFrequency = m_params.getValueDouble("MIN_PENALTY_FREQUENCY", 0.005);
    DEBUG_VAR(minFrequency);

    float avgIntensity = m_prmSpectra->averageIntensity();
    DEBUG_VAR(avgIntensity);

    m_penaltyMatrixMods->createFromModificationFreqs(modFreqs,
                                                     minPeakEquivalents,
                                                     maxPeakEquivalents,
                                                     minFrequency,
                                                     avgIntensity);

    string modFileName = m_params.getValue("OUTPUT_MOD_PENALTIES");
    DEBUG_VAR(modFileName);
    m_penaltyMatrixMods->saveMatrix(modFileName);

    string knowmModsFileName = m_params.getValue("OUTPUT_KNOWN_MODS");
    DEBUG_VAR(knowmModsFileName);
    if (!knowmModsFileName.empty()) {
      m_penaltyMatrixMods->saveKnownMods(knowmModsFileName);
    }
    string cleavagePenaltiesFileName = m_params.getValue("OUTPUT_CLEAVAGE_PEN");
    DEBUG_VAR(cleavagePenaltiesFileName);
    if (!cleavagePenaltiesFileName.empty()) {
      m_penaltyMatrixMods->saveCleavagePenalties(cleavagePenaltiesFileName);
    }

    string scanSpecificFileName = m_params.getValue("OUTPUT_SCAN_SPECIFIC_PENALTIES");
    DEBUG_VAR(scanSpecificFileName);
    if (scanSpecificFileName.empty()) {
      return true;
    }

    DEBUG_SCAN_SPECIFIC = m_params.getValueBool("DEBUG_SCAN_SPECIFIC", false);
    DEBUG_VAR(DEBUG_SCAN_SPECIFIC);

    // Scan specific penalties
    map<int, float> mapSpecToMass;
    for (int iStar = 0; iStar < m_prmSpectra->size(); iStar++) {
      mapSpecToMass[(*m_prmSpectra)[iStar].scan] = (*m_prmSpectra)[iStar].parentMass;
    }

    // This will make access easier
    map<int, map<int, float> > & mapSpecToMassToScore = *m_scanSpecificPenalties;

    float minScore = 100.0;
    float maxScore = 0.0;
    vector<float> pairScores;
    DEBUG_VAR(m_filteredPairs->size());
    for (int iPair = 0; iPair < m_filteredPairs->size(); iPair++) {
      int spec1 = (*m_filteredPairs)[iPair].spec1 + 1;  // These are 0 based
      int spec2 = (*m_filteredPairs)[iPair].spec2 + 1;
      if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(spec1);
      if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(spec2);

      if (mapSpecToMass[spec1] == 0.0 || mapSpecToMass[spec2] == 0.0) {
        continue;
      }

      float rawDiff1 = mapSpecToMass[spec1] - mapSpecToMass[spec2];
      if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(rawDiff1);

      // not interested in very small mass diffs
      if (abs(rawDiff1) < resolution * 2.0) {
        continue;
      }

      int   sign1    = rawDiff1 < 0 ? -1 : 1;
      if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(sign1);
      float rawDiff2 = mapSpecToMass[spec2] - mapSpecToMass[spec1];
      if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(rawDiff2);
      int   sign2    = rawDiff2 < 0 ? -1 : 1;
      if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(sign2);

      int massDiff1 = (int)(abs(rawDiff1) + 0.5) * sign1;
      if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(massDiff1);
      int massDiff2 = (int)(abs(rawDiff2) + 0.5) * sign2;
      if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(massDiff2);

      if (massDiff1 >= minModMass &&
          massDiff1 <= maxModMass &&
          (*m_filteredPairs)[iPair].score1 > mapSpecToMassToScore[spec1][massDiff1]) {
        mapSpecToMassToScore[spec1][massDiff1] = (*m_filteredPairs)[iPair].score1;
      }
      if (massDiff2 >= minModMass &&
          massDiff2 <= maxModMass &&
          (*m_filteredPairs)[iPair].score2 > mapSpecToMassToScore[spec2][massDiff2]) {
        mapSpecToMassToScore[spec2][massDiff2] = (*m_filteredPairs)[iPair].score2;
      }

      minScore = (*m_filteredPairs)[iPair].score1 < minScore ? (*m_filteredPairs)[iPair].score1 : minScore;
      maxScore = (*m_filteredPairs)[iPair].score1 > maxScore ? (*m_filteredPairs)[iPair].score1 : maxScore;
      minScore = (*m_filteredPairs)[iPair].score2 < minScore ? (*m_filteredPairs)[iPair].score2 : minScore;
      maxScore = (*m_filteredPairs)[iPair].score2 > maxScore ? (*m_filteredPairs)[iPair].score2 : maxScore;

      pairScores.push_back((*m_filteredPairs)[iPair].score1);
      pairScores.push_back((*m_filteredPairs)[iPair].score2);
    }
    DEBUG_VAR(minScore);
    DEBUG_VAR(maxScore);

    sort(pairScores.begin(), pairScores.end());
    vector<float> scaleVector(100);
    for (int i = 0; i < 100; i++) {
      scaleVector[i] = pairScores[pairScores.size() - 1 - pairScores.size() / 100 * i];
      if (DEBUG_SCAN_SPECIFIC) DEBUG_MSG(i << "  " << scaleVector[i]);
    }

    // We will scale the AlignGF scores to 1/2 min peak equivs
    minPeakEquivalents /= 2.0;
    if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(maxPeakEquivalents);
    if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(minPeakEquivalents);
    float penaltyRange = maxPeakEquivalents - minPeakEquivalents;
    if (DEBUG_SCAN_SPECIFIC) DEBUG_VAR(penaltyRange);

    // Scale the penalties from 1/2 min to max peak equivalents
    //    Just like "found" penalties
    map<int, map<int, float> >::iterator itr1 = mapSpecToMassToScore.begin();
    map<int, map<int, float> >::iterator itrEnd1 = mapSpecToMassToScore.end();
    for (; itr1 != itrEnd1; itr1++) {
      int scan = itr1->first;
      // Scale all the masses
      map<int, float>::iterator itr2 = itr1->second.begin();
      map<int, float>::iterator itrEnd2 = itr1->second.end();
      for (itr2 = itr1->second.begin(); itr2 != itrEnd2; itr2++) {
        int mass = itr2->first;
        float score = itr2->second;

        int rank = 100;
        for (int i = 0; i < 100; i++) {
          if (itr2->second >= scaleVector[i]) {
            rank = i;
            break;
          }
        }
        itr2->second = minPeakEquivalents + rank * 0.01 * penaltyRange;
        if (DEBUG_SCAN_SPECIFIC) DEBUG_MSG(scan << "  " << mass << "  " << score << "  " << rank << "  " << itr2->second);
      }
    }

    ofstream ofs(scanSpecificFileName.c_str());
    if (!ofs) {
      ERROR_MSG("Could not write to " << scanSpecificFileName);
    }
    ofs << "scan\tmass\tscore" << endl;

    for (itr1 = mapSpecToMassToScore.begin(); itr1 != itrEnd1; itr1++) {
      int scan = itr1->first;

      map<int, float>::iterator itr2 = itr1->second.begin();
      map<int, float>::iterator itrEnd2 = itr1->second.end();
      for (; itr2 != itrEnd2; itr2++) {
        int mass = itr2->first;
        float score = itr2->second;
        ofs << scan << "\t"  << mass << "\t" << -score << endl;
      }
    }
    
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPenaltyGen::loadInputData(void)
  {
    if (ownInput)
    {
      if (!m_filteredPairs) {
        m_filteredPairs = new SpectrumPairSet;
      }
      if (!m_prmSpectra) {
        m_prmSpectra = new SpecSet;
      }
      if (!m_scanSpecificPenalties) {
        m_scanSpecificPenalties = new map<int, map<int, float> >();
      }

    }
    if (ownOutput)
    {
      // Load amino acid masses
      AAJumps jumps(1);
      if (m_params.exists("AMINO_ACID_MASSES")) {
        DEBUG_MSG("Loading: " << m_params.getValue("AMINO_ACID_MASSES"));
        jumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true);
      }

      float resolution = m_params.getValueDouble("ALIGNMENT_RESOLUTION", 1.0);
      float minFrequency = m_params.getValueDouble("MIN_PENALTY_FREQUENCY", 0.005);
      float unknownPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_PENALTY",
                                               1.0);
      float knownModPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_KNOWN_PENALTY",
                                                1.0);
      float minModMass = m_params.getValueDouble("MIN_MOD_MASS", -100.0);
      float maxModMass = m_params.getValueDouble("MAX_MOD_MASS", 100.0);
      float unknownMultiplier =
          m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER", 2.0);
      DEBUG_VAR(resolution);
      m_penaltyMatrixMods= new PenaltyMatrix(jumps,
                                        resolution,
                                        knownModPenalty,
                                        unknownPenalty,
                                        unknownMultiplier,
                                        minModMass,
                                        maxModMass);
    }
    
    string filteredPairsFilename = m_params.getValue("INPUT_FILTERED_PAIRS");
    DEBUG_VAR(filteredPairsFilename);
    if (!filteredPairsFilename.empty()) {
      m_filteredPairs->loadFromBinaryFile(filteredPairsFilename);
      DEBUG_VAR(m_filteredPairs->size());
    }
        
    string prmPklbinFilename = m_params.getValue("INPUT_PRM_PKLBIN");
    if (!prmPklbinFilename.empty()) {
      DEBUG_MSG("Loading: " << prmPklbinFilename);
      if (m_prmSpectra->loadPklBin(prmPklbinFilename.c_str()) < 0) {
        ERROR_MSG("Error reading input PRM spectra.");
        return false;
      }
    }

    string knowmModsFileName = m_params.getValue("INPUT_KNOWN_MODS");
    DEBUG_VAR(knowmModsFileName);
    if (!knowmModsFileName.empty()) {
      m_penaltyMatrixMods->loadKnownModifications(knowmModsFileName);
    }
    
    string cleavagePenaltiesFileName = m_params.getValue("INPUT_CLEAVAGE_PENALTIES");
    DEBUG_VAR(cleavagePenaltiesFileName);
    if (!cleavagePenaltiesFileName.empty()) {
      m_penaltyMatrixMods->loadCleavagePenalties(cleavagePenaltiesFileName);
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPenaltyGen::saveOutputData(void)
  {
    DEBUG_TRACE;

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPenaltyGen::saveInputData(std::vector<std::string> & filenames)
  {
    return true;
  }
  

  // -------------------------------------------------------------------------
  bool ExecPenaltyGen::loadOutputData(void)
  {
    return true;
  }
 
  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecPenaltyGen::split(int numSplit)
  {
    //EMPTY
  }

  // -------------------------------------------------------------------------
  bool ExecPenaltyGen::merge(void)
  {
    //EMPTY
  }
 
  // -------------------------------------------------------------------------
  bool ExecPenaltyGen::validateParams(std::string & error)
  {
    m_isValid = false;

    //VALIDATE_PARAM_EXIST("INPUT_KNOWN_MODS");
    //VALIDATE_PARAM_EXIST("INPUT_CLEAVAGE_PENALTIES");
    VALIDATE_PARAM_EXIST("OUTPUT_MOD_PENALTIES");
    //VALIDATE_PARAM_EXIST("OUTPUT_MOD_KNOWN");
    //VALIDATE_PARAM_EXIST("OUTPUT_MOD_CLEAVAGE");

    m_isValid = true;
    return true;
  }

}


