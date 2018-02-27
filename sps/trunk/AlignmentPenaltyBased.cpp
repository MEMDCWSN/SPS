#include "AlignmentPenaltyBased.h"
#include "Logger.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <limits.h>
#include <stdio.h>

static bool DEBUG_RANGE = false;
static bool DEBUG_RANGE2 = false;
static bool DEBUG_AAS = false;
static bool DEBUG_SPECS = false;
static bool DEBUG_ALIGN = false;
static bool DEBUG_ALIGN1 = false;
static bool DEBUG_ALIGN2 = false;
static bool DEBUG_ALIGN3 = false;
static bool DEBUG_ALIGN4 = false;
static bool DEBUG_ALIGN_EXACT = false;
static bool DEBUG_ALIGN_NTERM = false;
static bool DEBUG_ALIGN_CTERM = false;
static bool DEBUG_ALIGN_TERM = false;
static bool DEBUG_ALIGN_GAP = false;
static bool DEBUG_GAP_ANNO = false;
static bool DEBUG_CACHE = false;
static bool DEBUG_ALIGN_LARGE = false;
static bool DEBUG_GAP_VECTOR = false;

#define DEBUG_MATCH {                                               \
  DEBUG_MSG(strAA << "  " << deltaDbMass <<                         \
                     "  " << deltaSpecMass << "  " << deltaMasses); \
  DEBUG_VAR(predDbSpecIdx);                                         \
  DEBUG_VAR(predSpecIdx);                                           \
  DEBUG_VAR(dbSpecIdx);                                             \
  DEBUG_VAR(specIdx);                                               \
  DEBUG_VAR(spec[specIdx][1]);                                      \
  DEBUG_VAR(penalty);                                               \
  DEBUG_VAR(score);                                                 \
}

#define DEBUG_MATCH2 {                                              \
  DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);               \
  DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);                       \
  DEBUG_VAR(matchType[dbSpecIdx][specIdx]);                         \
  DEBUG_TRACE;                                                      \
}

// Currently the R parameter (for penalty reversal is not used)
#define REPLACE_SCORE_IF_BETTER_WITH_CLEAVPEN(X, Y, Z, R) {         \
  float penalty = X;                                                \
  penalty += computeCleavagePenalty(dbSeq,                          \
                         predDbSpecIdx,                             \
                         dbSpecIdx,                                 \
                         predSpecIdx == 0,                          \
                         specIdx == spec.size() - 1,                \
                         R,                                         \
                         Z) * avgPeakIntensity;                     \
  float score = matchMatrix[predDbSpecIdx][predSpecIdx] +           \
                         spec[specIdx][1] + penalty;                \
  if (Z) DEBUG_MATCH;                                               \
  if (score > matchMatrix[dbSpecIdx][specIdx]) {                    \
    matchMatrix[dbSpecIdx][specIdx] = score;                        \
    matchPtr[dbSpecIdx][specIdx] =                                  \
             make_pair<int,int>(predDbSpecIdx, predSpecIdx);        \
    matchType[dbSpecIdx][specIdx] = Y;                              \
    if (Z) DEBUG_MATCH2;                                            \
  }                                                                 \
}

#define REPLACE_SCORE_IF_BETTER(X, Y, Z) {                          \
  float penalty = X;                                                \
  float score = matchMatrix[predDbSpecIdx][predSpecIdx] +           \
                         spec[specIdx][1] + penalty;                \
  if (Z) DEBUG_MATCH;                                               \
  if (score > matchMatrix[dbSpecIdx][specIdx]) {                    \
    matchMatrix[dbSpecIdx][specIdx] = score;                        \
    matchPtr[dbSpecIdx][specIdx] =                                  \
             make_pair<int,int>(predDbSpecIdx, predSpecIdx);        \
    matchType[dbSpecIdx][specIdx] = Y;                              \
    if (Z) DEBUG_MATCH2;                                            \
  }                                                                 \
}


const int KMER_LENGTH = 4;
const float MINIMUM_GAP_MASS_TO_SPLIT = 2.05;  // A little more than 2 (be more scientific later)

// Assuming vectors of size ~1000 then the cache size will end up in the MB region (not GB)
#define USE_CACHE 0
const float CACHE_LIMIT = 10000.0;

// Used in the conversion from parent mass to penalty multiplier
const float PENALTY_PM_DIVIDER = 1000.0;


namespace specnets
{
  // -------------------------------------------------------------------------
  void debugVector(vector<float> & scores)
  {
    for (int i = 0; i < scores.size(); i++) {
      cout << scores[i] << "\t";
      if (i != 0 && i % 50 == 0) {
        cout << endl;
      }
    }
    cout << endl;
    return;
  }

  //--------------------------------------------------------------------------
  string getCleanAnnotation(string & annotationIn)
  {
    string cleanAnnotation; 
    static string aminoAcids("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    for (int iChar = 0; iChar < annotationIn.length(); iChar++) {
      if (aminoAcids.find_first_of(annotationIn[iChar]) != string::npos) {
        cleanAnnotation += annotationIn[iChar];
      }
    }
    return cleanAnnotation;
  }
    
  //--------------------------------------------------------------------------
  float getScanSpecificPenalty(map<int, float> & scanSpecificPenalties,
                               float deltaMasses,
                               float avgPeakIntensity)
  {
    int sign        = deltaMasses < 0 ? -1 : 1;
    int roundedMass = (int)(abs(deltaMasses) + 0.5) * sign;
    if (scanSpecificPenalties.find(roundedMass) == scanSpecificPenalties.end()) {
      //DEBUG_MSG(deltaMasses << "  Max penalty");
      return -1000000.0; // A very large penalty
    }
    //DEBUG_MSG(deltaMasses << "  " << scanSpecificPenalties[roundedMass]);
    return scanSpecificPenalties[roundedMass] * avgPeakIntensity;
  }

  //--------------------------------------------------------------------------
  void debugAAs(PenaltyMatrix * modPenaltyMatrix,
		PenaltyMatrix * blossumPenaltyMatrix)
  {
    for (int i = 0; i < 26; i++) {
      string strAA;  // Create a dummy string
      strAA.append(1, char('A' + i)); // Append the DB AA char
      DEBUG_MSG(strAA << " = " << modPenaltyMatrix->getMass(strAA));
    }
    for (int i = 0; i < 26; i++) {
      string strAA;  // Create a dummy string
      strAA.append(1, char('A' + i)); // Append the DB AA char
      DEBUG_MSG(strAA << " = " << blossumPenaltyMatrix->getMass(strAA));
    }
    return;
  }

  //--------------------------------------------------------------------------
  void debugDBSpec(char * dbSeq, Spectrum & dbSpec)
  {
    DEBUG_MSG(0 << ", " << dbSpec[0][0] << ", " << dbSpec[0][0]);
    for (int dbSpecIdx = 1; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      DEBUG_MSG(dbSeq[dbSpecIdx - 1] << ", " << dbSpecIdx << ", " << dbSpec[dbSpecIdx][0] << ", " << dbSpec[dbSpecIdx][0] - dbSpec[dbSpecIdx-1][0]);
    }
    return;
  }
  
  //--------------------------------------------------------------------------
  void debugSpec(Spectrum & spec)
  {
    DEBUG_MSG(0 << ", " << spec[0][0] << ", " << 0 << ", " << spec[0][1]);
    for (int specIdx = 1; specIdx < spec.size(); specIdx++) {
      DEBUG_MSG(specIdx << ", " << spec[specIdx][0] << ", " << spec[specIdx][0] - spec[specIdx-1][0] << ", " << spec[specIdx][1]);
    }
    return;
  }

  //--------------------------------------------------------------------------
  void debugMatrices(Spectrum & 	               spec,
                     char *                      dbSeq,
                     vector<vector< float > > &  matchMatrix,
                     vector<vector< string > > & matchType)
  {
    DEBUG_VAR(spec.size());
    DEBUG_VAR(matchMatrix.size());
    DEBUG_VAR(strlen(dbSeq));
    for (int dbSpecIdx = 0; dbSpecIdx < matchMatrix.size(); dbSpecIdx++) {
      bool outputRow = false;
      // Only output "interesting" rows
      for (int specIdx = 1; specIdx < spec.size(); specIdx++) {
        if (matchMatrix[dbSpecIdx][specIdx] != -(float)INT_MAX) {
          outputRow = true;
          break;
        }
      }
      if (outputRow) {
        cout << dbSpecIdx << "\t" << dbSeq[dbSpecIdx] << "\t";
        for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
          cout << matchMatrix[dbSpecIdx][specIdx] << "\t";
        }
        cout << endl;
      }
    }

    for (int dbSpecIdx = 0; dbSpecIdx < matchMatrix.size(); dbSpecIdx++) {
      bool outputRow = false;
      // Only output "interesting" rows
      for (int specIdx = 1; specIdx < spec.size(); specIdx++) {
        if (matchMatrix[dbSpecIdx][specIdx] != -(float)INT_MAX) {
          outputRow = true;
          break;
        }
      }
      if (outputRow) {
        cout << dbSpecIdx << "\t";
        for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
          cout << matchType[dbSpecIdx][specIdx] << "\t";
        }
        cout << endl;
      }
    }
    return;
  }

  //--------------------------------------------------------------------------
  AlignmentPenaltyBased::AlignmentPenaltyBased(int             maxSpecGap,
                                               PenaltyMatrix * modPenaltyMatrix,
                                               PenaltyMatrix * blossumPenaltyMatrix,
                                               float           penaltyAlpha,
                                               float           penaltyBeta,
                                               float           maxMod,
                                               float           minMod)
  : m_maxSpecGap(maxSpecGap + 1),
    m_modPenaltyMatrix(modPenaltyMatrix),
    m_blossumPenaltyMatrix(blossumPenaltyMatrix),
    m_penaltyAlpha(penaltyAlpha),
    m_penaltyBeta(penaltyBeta),
    m_maxMod(maxMod),
    m_minMod(minMod)
  {
    if (m_modPenaltyMatrix != 0x0 || m_modPenaltyMatrix != 0x0) {
      m_mapGapVectors.clear();
      createGapVectors();
    }
    m_cacheHitsTotal = 0;
    m_cacheMissesTotal = 0;
  }

  //--------------------------------------------------------------------------
  AlignmentPenaltyBased::~AlignmentPenaltyBased()
  {
    // Nothing here
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::clearCache(void)
  {
    m_mapCachedGapVectors.clear();
  }

  //--------------------------------------------------------------------------
  bool AlignmentPenaltyBased::isCached(string & aaString)
  {
    return m_mapCachedGapVectors.find(aaString) != m_mapCachedGapVectors.end();
  }
  
  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::sortStringLetters(string & toSort)
  {
    int startIndex = 0;
    int endIndex = toSort.length() - 1;
    if (toSort[0] == m_modPenaltyMatrix->ntermChar() ||
        toSort[0] == '=') {
      startIndex += 2;
    }
    if (toSort[toSort.length() - 1] == m_modPenaltyMatrix->ctermChar()) endIndex -= 2;
    // Initialize the vector
    for (int i = startIndex; i <= endIndex; i++) {
      for (int j = i + 1; j <= endIndex; j++) {
        if (toSort[i] > toSort[j]) {
          char temp = toSort[i];
          toSort[i] = toSort[j];
          toSort[j] = temp;
        }
      }
    }
  }
  
  //---------------------------------------------------------------------------------------
  // This method adds to the cache all the maxLength kMers that make up the protein string
  //---------------------------------------------------------------------------------------
  void AlignmentPenaltyBased::cacheProteinStrings(string & proteinString, int maxLength)
  {
    if (DEBUG_CACHE) DEBUG_VAR(proteinString);
    m_mapCachedGapVectors.clear();
    vector<float> scoresCombined;
    for (int iLength = maxLength - m_kMer + 1; iLength <= maxLength; iLength++) {
      if (DEBUG_CACHE) DEBUG_VAR(iLength);
      for (int iStart = 0; iStart < proteinString.length() - iLength; iStart++) {
	
        if (DEBUG_CACHE) DEBUG_VAR(iStart);
        string full = proteinString.substr(iStart, iLength);
        sortStringLetters(full);
        if (DEBUG_CACHE) DEBUG_VAR(full);
        if (isCached(full)) {
          continue;
        }
	
        string left = full.substr(0, m_kMer);
        sortStringLetters(left);
        if (DEBUG_CACHE) DEBUG_VAR(left);
        vector<float> scoresLeft = m_mapGapVectors[left];

        string right = full.substr(m_kMer);
        sortStringLetters(right);
        if (DEBUG_CACHE) DEBUG_VAR(right);
        vector<float> scoresRight = m_mapGapVectors[left];
	
        initializeGapVector(scoresCombined);
        combineVectors(scoresLeft, scoresRight, scoresCombined);
	 m_mapCachedGapVectors[full] = scoresCombined;
      }
    }
    if (DEBUG_CACHE) DEBUG_MSG("Done cacheing protein strings");
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::computeAlignment(
            Spectrum & 	spec,
            Spectrum & 	dbSpec,
            char * 		  dbSeq,
            int 		    dbIndex,
            int 		    matchOrientation,
            set<float> &    startRange,
            map<int, float> & scanSpecificPenalties,
            int 		    minMatchedPeaks,
            int 		    maxGapSize,
            float 		  pmTolerance,
            float 		  tolerance,
            float       avgPeakIntensity,
            bool 		    enforceEndpeaks)
  {
    if (DEBUG_CACHE) DEBUG_MSG("START PENALTY ALIGNMENT");
    if (m_mapCachedGapVectors.size() > CACHE_LIMIT * 0.8) {
      clearCache();
    }
    m_cacheHits = 0;
    m_cacheMisses = 0;
    
    if (spec.size() == 0) {
      WARN_MSG("Spectrum size is 0");
      return;
    }

    if (avgPeakIntensity == 0) {
      WARN_MSG("Average peak intensity is 0");
      return;
    }

    // Cap the length modifier at 2000 daltons
    //   Empirically we note that the peak intensities don't continue to rise
    //   indefinitely, instead they taper off after ~2000 daltons
    float lengthModifierMass = min(spec.parentMass, (float)2000.0);

    // We modify the avgPeakIntensity by a length modifier
    //   This is equivalent to modifying the peak equivalents for penalties
    float penaltyLengthModifier = lengthModifierMass / PENALTY_PM_DIVIDER;
    avgPeakIntensity *= penaltyLengthModifier;
    DEBUG_VAR(penaltyLengthModifier);

    if (DEBUG_ALIGN) {
      DEBUG_VAR(DEBUG_AAS);
      DEBUG_VAR(DEBUG_RANGE);
      DEBUG_VAR(DEBUG_RANGE2);
      DEBUG_VAR(DEBUG_SPECS);
      DEBUG_VAR(DEBUG_ALIGN);
      DEBUG_VAR(DEBUG_ALIGN1);
      DEBUG_VAR(DEBUG_ALIGN2);
      DEBUG_VAR(DEBUG_ALIGN3);
      DEBUG_VAR(DEBUG_ALIGN4);
      DEBUG_VAR(DEBUG_ALIGN_EXACT);
      DEBUG_VAR(DEBUG_ALIGN_NTERM);
      DEBUG_VAR(DEBUG_ALIGN_CTERM);
      DEBUG_VAR(DEBUG_GAP_ANNO);
      DEBUG_VAR(DEBUG_CACHE);
      DEBUG_VAR(DEBUG_ALIGN_LARGE);
      DEBUG_VAR(DEBUG_GAP_VECTOR);

      DEBUG_VAR(enforceEndpeaks);
      DEBUG_VAR(maxGapSize);
      DEBUG_VAR(startRange.size());
      DEBUG_VAR(pmTolerance);
      DEBUG_VAR(tolerance);
      DEBUG_VAR(spec.parentMass);
      DEBUG_VAR(avgPeakIntensity);
    }

    if (DEBUG_AAS) debugAAs(m_modPenaltyMatrix, m_blossumPenaltyMatrix);
    if (DEBUG_SPECS) debugDBSpec(dbSeq, dbSpec);
    if (DEBUG_SPECS) debugSpec(spec);

    vector<vector< float > >         matchMatrix;
    vector<vector< pair<int,int> > > matchPtr;
    vector<vector< string > >        matchType;
    vector<vector<char> >            startFlags(spec.size());

    // Initialize the starting flags
    int firstValidDbIndex;
    setStartingFlagArray(spec, dbSpec, startRange, startFlags, firstValidDbIndex, tolerance);
    if (DEBUG_ALIGN) DEBUG_VAR(firstValidDbIndex);

    matchMatrix.resize(dbSpec.size());
    matchPtr.resize(dbSpec.size());
    matchType.resize(dbSpec.size());

    char dbGapString[strlen(dbSeq) + 1];
    
    //---------------------------------
    //  Initialize the matrices
    //---------------------------------
    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      matchMatrix[dbSpecIdx].resize(spec.size());
      matchPtr[dbSpecIdx].resize(spec.size());
      matchType[dbSpecIdx].resize(spec.size());
      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
        matchMatrix[dbSpecIdx][specIdx] = -(float)INT_MAX;
        matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(-1,-1);
        matchType[dbSpecIdx][specIdx] = "?";
      }
    }

    //---------------------------------
    // Compute the match between the two spectra
    //---------------------------------
    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {

      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
    
        if (DEBUG_RANGE2) DEBUG_MSG(dbSpecIdx << ", " << specIdx << ", " << (int)startFlags[specIdx][dbSpecIdx]);

        if (!enforceEndpeaks or specIdx == 0) {
          // Matching can start on any mass, assume no predecessor 
          matchMatrix[dbSpecIdx][specIdx] = spec[specIdx][1];
        }

        if (DEBUG_RANGE2) DEBUG_MSG(matchMatrix[dbSpecIdx][specIdx]);

        for (int predDbSpecIdx = dbSpecIdx - 1;
             (predDbSpecIdx >= 0) && (predDbSpecIdx >= firstValidDbIndex) && (predDbSpecIdx >= dbSpecIdx - maxGapSize * 2);
              predDbSpecIdx--) {

          for (int predSpecIdx = specIdx - 1;
        	   (predSpecIdx >= 0);
        	   predSpecIdx--) {

            if (!startFlags[predSpecIdx][predDbSpecIdx]) {
              // if location is not a valid start, skip scoring
              if (DEBUG_RANGE2) DEBUG_MSG("Skipping " << predDbSpecIdx << ", " << predSpecIdx);
              continue;
            }
            if (matchMatrix[predDbSpecIdx][predSpecIdx] == -(float)INT_MAX) {
              // Must build on a real score (can make anything out of -infinity)
              if (DEBUG_RANGE2) DEBUG_MSG("Skipping " << predDbSpecIdx << ", " << predSpecIdx);
              continue;
            }
        
            //makeDbString(dbGapString, dbSeq, predDbSpecIdx, dbSpecIdx);
            int dbGapLength = dbSpecIdx - predDbSpecIdx;//strlen(dbGapString);
            //int dbGapLength = dbSpecIdx - predDbSpecIdx;
            float dbGapMass = dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0];
            float dbAAMinMass = (dbSpecIdx - predDbSpecIdx) * 57.0;
            float specGapMass = spec[specIdx][0] - spec[predSpecIdx][0];

            if (DEBUG_ALIGN1) DEBUG_MSG(predDbSpecIdx << "," << dbSpecIdx << "  " << predSpecIdx << "," << specIdx);
            if (DEBUG_ALIGN1) DEBUG_MSG(dbGapLength << "  " << dbGapMass << "  " << specGapMass << "  " << dbAAMinMass);
            
            float deltaDbMass = dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0];
            float deltaSpecMass = spec[specIdx][0] - spec[predSpecIdx][0];
            float deltaMasses = deltaSpecMass - deltaDbMass;

            makeDbString(dbGapString, dbSeq, predDbSpecIdx, dbSpecIdx);
            string strAA(dbGapString);
            if (DEBUG_ALIGN1) DEBUG_MSG(strAA << "  " << deltaDbMass << "  " << deltaSpecMass << "  " << deltaMasses);


            if (abs(deltaMasses) < tolerance ||
                abs(deltaMasses) < pmTolerance && specIdx == spec.size() - 1) {
              // An exact match regardless of number of AAs
              // The second case takes care of cterm matches within parent mass tolerance
              if (DEBUG_ALIGN_EXACT) DEBUG_MSG("Exact Match");
              
              string pmGapString("");
              // Although this is an "exact" match.. there may be a gap 
              //    due to a parent mass error.. if so we need to add anno for it
              //    if we do not, later processing (spec prob) may not work
              int pmGapMass = int(abs(deltaMasses + 0.5)) * (deltaMasses > 0 ? 1 : -1);
              if (pmGapMass > 0) {
                pmGapString = "[" + massToString(pmGapMass) + "]" ;
              }

              REPLACE_SCORE_IF_BETTER_WITH_CLEAVPEN(0.0,
                                      strAA + pmGapString,
                                      DEBUG_ALIGN_EXACT,
                                      false);

            } else if (!enforceEndpeaks) {
              continue;

            } else if (dbGapLength > maxGapSize &&
                       deltaMasses > m_minMod &&
                       deltaMasses < m_maxMod) {

              // A very large gap with a single mod
              if (DEBUG_ALIGN_LARGE) DEBUG_MSG("Large Gap Match");
              // Loop over all AAs in the gap and find the minimum penalty for
              //   a mod of (AA,m) for all AA in the gap
              float minPenalty = -1000000.0;
              for (int iaa = 0; iaa < strAA.length(); iaa++) {
                float modPenalty = (*m_modPenaltyMatrix)(strAA[iaa], deltaMasses, avgPeakIntensity);
                if (DEBUG_ALIGN_LARGE) DEBUG_MSG(strAA[iaa] << "  " <<
                                                 deltaMasses << "  " <<
                                                 avgPeakIntensity << "  " << modPenalty);
                if (modPenalty > minPenalty) {
                  minPenalty = modPenalty;
                }
              }
              // Check for scan specific mass penalty
              float scanSpecificPenalty = getScanSpecificPenalty(scanSpecificPenalties,
                                                                 deltaMasses,
                                                                 avgPeakIntensity);
              if (scanSpecificPenalty > minPenalty) {
                minPenalty = scanSpecificPenalty;
              }

              REPLACE_SCORE_IF_BETTER(minPenalty,
                                      "(" + strAA + "," + massToString(deltaMasses) + ")",
                                      DEBUG_ALIGN_LARGE);

            } else if (dbGapLength > maxGapSize) {
              // A very large gap that we won't handle
              continue;

            } else if (dbSpecIdx - predDbSpecIdx == 1) {

              if (deltaDbMass > 400.0) continue;
              if (deltaSpecMass > 400.0) continue;
            
              if (deltaMasses < m_minMod || deltaMasses > m_maxMod) {
                continue;
              }

              if (predSpecIdx == 0 &&
                  m_modPenaltyMatrix->isNterm(strAA, deltaMasses)) {
              
                if (DEBUG_ALIGN_NTERM) DEBUG_MSG("Nterm AA Modification Match");
                REPLACE_SCORE_IF_BETTER(m_modPenaltyMatrix->getKnownPenalty(avgPeakIntensity),
                                        "(" + strAA + "," + massToString(deltaMasses) + ")",
                                        DEBUG_ALIGN_NTERM);

              } else if (specIdx == spec.size() - 1 &&
                         m_modPenaltyMatrix->isCterm(strAA, deltaMasses)) {
                if (DEBUG_ALIGN_NTERM) DEBUG_MSG("Cterm AA Modification Match");
                REPLACE_SCORE_IF_BETTER_WITH_CLEAVPEN(m_modPenaltyMatrix->getKnownPenalty(avgPeakIntensity),
                                        "(" + strAA + "," + massToString(deltaMasses) + ")",
                                        DEBUG_ALIGN_NTERM,
                                        true);

              } else if (m_modPenaltyMatrix->isInMatrix(strAA, deltaMasses)) {

                if (DEBUG_ALIGN2) DEBUG_MSG("Known/Putative Modification Match");
                REPLACE_SCORE_IF_BETTER((*m_modPenaltyMatrix)(strAA, deltaMasses, avgPeakIntensity),
                                       "(" + strAA + "," + massToString(deltaMasses) + ")",
                                       DEBUG_ALIGN2);

                // Check for scan specific mass penalty
                float scanSpecificPenalty = getScanSpecificPenalty(scanSpecificPenalties,
                                                                   deltaMasses,
                                                                   avgPeakIntensity);
                REPLACE_SCORE_IF_BETTER(scanSpecificPenalty,
                                       "(" + strAA + "," + massToString(deltaMasses) + ")",
                                       DEBUG_ALIGN2);


              } else if (predSpecIdx == 0 &&
                         m_modPenaltyMatrix->isNterm(deltaMasses)) {
              
                if (DEBUG_ALIGN_NTERM) DEBUG_MSG("Nterm Modification Match");
                REPLACE_SCORE_IF_BETTER(m_modPenaltyMatrix->getKnownPenalty(avgPeakIntensity),
                                        "[" + massToString(deltaMasses) + "]" + strAA,
                                        DEBUG_ALIGN_NTERM);
                
              } else if (specIdx == spec.size() - 1 &&
                         m_modPenaltyMatrix->isCterm(deltaMasses)) {
              
                if (DEBUG_ALIGN_NTERM) DEBUG_MSG("Cterm Modification Match");
                REPLACE_SCORE_IF_BETTER_WITH_CLEAVPEN(m_modPenaltyMatrix->getKnownPenalty(avgPeakIntensity),
                                        strAA + "[" + massToString(deltaMasses) + "]",
                                        DEBUG_ALIGN_NTERM,
                                        true);

              } else if (!m_modPenaltyMatrix->isInMatrix(strAA, deltaMasses) &&
                         m_modPenaltyMatrix->getMass(strAA) + deltaMasses >= 57.0) {

                if (DEBUG_ALIGN2) DEBUG_MSG("Unknown Match");
                if (specIdx == spec.size() - 1) {
                  REPLACE_SCORE_IF_BETTER_WITH_CLEAVPEN(m_modPenaltyMatrix->getUnknownPenalty(avgPeakIntensity),
                                          "(" + strAA + "," + massToString(deltaMasses) + ")",
                                          DEBUG_ALIGN2,
                                          true);
                  // Check for scan specific mass penalty
                  float scanSpecificPenalty = getScanSpecificPenalty(scanSpecificPenalties,
                                                                     deltaMasses,
                                                                     avgPeakIntensity);
                  REPLACE_SCORE_IF_BETTER_WITH_CLEAVPEN(scanSpecificPenalty,
                                         "(" + strAA + "," + massToString(deltaMasses) + ")",
                                         DEBUG_ALIGN2,
                                         true);
                } else {
                  REPLACE_SCORE_IF_BETTER_WITH_CLEAVPEN(m_modPenaltyMatrix->getUnknownPenalty(avgPeakIntensity),
                                          "(" + strAA + "," + massToString(deltaMasses) + ")",
                                          DEBUG_ALIGN2,
                                          false);
                  // Check for scan specific mass penalty
                  float scanSpecificPenalty = getScanSpecificPenalty(scanSpecificPenalties,
                                                                     deltaMasses,
                                                                     avgPeakIntensity);
                  REPLACE_SCORE_IF_BETTER_WITH_CLEAVPEN(scanSpecificPenalty,
                                         "(" + strAA + "," + massToString(deltaMasses) + ")",
                                         DEBUG_ALIGN2,
                                         false);
                }

              }

            } else if (dbAAMinMass < specGapMass && 
                       (dbSpecIdx - predDbSpecIdx <= maxGapSize) && 
                       (spec[specIdx][0] - spec[predSpecIdx][0] < m_maxSpecGap - 2) &&
                       specGapMass > dbGapMass + m_minMod &&
                       specGapMass < dbGapMass + m_maxMod &&
                       matchMatrix[predDbSpecIdx][predSpecIdx] + spec[specIdx][1] > matchMatrix[dbSpecIdx][specIdx]) {
 
              if (DEBUG_ALIGN_GAP) DEBUG_MSG("Gap Match");
              if (DEBUG_ALIGN_GAP) DEBUG_MSG(dbGapString << "  " << deltaMasses);
              string matchString;
              float penalty = -(float)INT_MAX;
              int specGapLength = (int)round(specGapMass);
              if (DEBUG_ALIGN_GAP) DEBUG_VAR(specGapLength);
              makeDbString(dbGapString, dbSeq, predDbSpecIdx, dbSpecIdx);
              int ntermOrCterm = predSpecIdx == 0 ? 1 : 0;  // 1 if nterm, 0 if not
              ntermOrCterm = specIdx == spec.size() - 1 ? 2 : ntermOrCterm;  // 2 if cterm
              if (DEBUG_ALIGN_GAP) DEBUG_VAR(ntermOrCterm);
              penalty = getGapPenalty(specGapLength, dbGapString, avgPeakIntensity, ntermOrCterm, false);
              if (DEBUG_ALIGN_GAP) DEBUG_VAR(penalty);
              float scanSpecificPenalty = getScanSpecificPenalty(scanSpecificPenalties,
                                                                 deltaMasses,
                                                                 avgPeakIntensity);
              if (DEBUG_ALIGN_GAP) DEBUG_VAR(scanSpecificPenalty);
              penalty = max(penalty, scanSpecificPenalty); // remember penalties are neg
              if (DEBUG_ALIGN_GAP) DEBUG_VAR(penalty);

              if (penalty == 0.0) {
                matchString = dbGapString;
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(matchString);
              } else {
                char buffer[20];
                sprintf(buffer, "%0.2f", deltaMasses);
                // Gap annotation (place holder because we don't compute exact)
                matchString = "(" + string(dbGapString) + "," + string(buffer) + ")";
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(matchString);
              }

              float bestGapPenalty = penalty;
              if (DEBUG_ALIGN_NTERM) DEBUG_VAR(bestGapPenalty);

              // If we are at nterm then loop over all the nterm mods and try them with gap
              if (predSpecIdx == 0) {
                const set<float> & setNtermMods = m_modPenaltyMatrix->getNtermMods();
                if (DEBUG_ALIGN_NTERM) DEBUG_VAR(setNtermMods.size());

                if (DEBUG_ALIGN_NTERM) DEBUG_MSG("Try Nterm Match");

                set<float>::const_iterator itrs = setNtermMods.begin();
                set<float>::const_iterator itrsEnd = setNtermMods.end();
                for ( ; itrs != itrsEnd; itrs++) {

                  float ntermMod = *itrs;
                  float newSpecGapMass = specGapMass - ntermMod;
                  float newDeltaMasses = deltaMasses - ntermMod;

                  if (DEBUG_ALIGN_NTERM) DEBUG_TRACE;
                  if (DEBUG_ALIGN_NTERM) DEBUG_VAR(ntermMod);
                  if (DEBUG_ALIGN_NTERM) DEBUG_VAR(newSpecGapMass);

                  if (dbAAMinMass < newSpecGapMass &&
                      newSpecGapMass > dbGapMass + m_minMod &&
                      newSpecGapMass < dbGapMass + m_maxMod) {

                    int specGapLength = (int)round(newSpecGapMass);
                    int ntermOrCterm = 0;  // Can't allow double nterm mods
                    if (fabs(ntermMod) <= pmTolerance && deltaMasses < 57.0) {
                      ntermOrCterm = 1;  // However the 1.0 masses are not TRUE nterm mods
                                         // They are parent mass adjustments
                    }
                    float ntermPenalty = getGapPenalty(specGapLength, 
                                                       dbGapString, 
                                                       avgPeakIntensity, 
                                                       ntermOrCterm, 
                                                       true);

                    float scanSpecificPenalty = getScanSpecificPenalty(scanSpecificPenalties,
                                                                       newDeltaMasses,
                                                                       avgPeakIntensity);
                    if (DEBUG_ALIGN_GAP) DEBUG_VAR(scanSpecificPenalty);
                    ntermPenalty = max(ntermPenalty, scanSpecificPenalty); // remember penalties are neg
                    if (DEBUG_ALIGN_GAP) DEBUG_VAR(ntermPenalty);

                    if (DEBUG_ALIGN_NTERM) DEBUG_VAR(specGapLength);
                    if (DEBUG_ALIGN_NTERM) DEBUG_VAR(ntermOrCterm);
                    if (DEBUG_ALIGN_NTERM) DEBUG_VAR(ntermPenalty);

                    if (ntermPenalty > bestGapPenalty) {
                      if (ntermPenalty == 0.0) {
                        matchString = "[" + massToString(ntermMod)+ "]" + strAA;
                      } else {
                        matchString = "[" + massToString(ntermMod)+ "]" +
                                      "(" + strAA + "," +
                                      massToString(deltaMasses - ntermMod) + ")";
                      }
                      bestGapPenalty = ntermPenalty;

                      if (DEBUG_ALIGN_NTERM) DEBUG_VAR(matchString);
                      if (DEBUG_ALIGN_NTERM) DEBUG_VAR(bestGapPenalty);
                    } // if (ntermPenalty > bestGapPenalty) {
                  } // if (dbAAMinMass < newSpecGapMass  ...
                } // or ( ; itrs != itrsEnd; itrs++) {
              } // if (predSpecIdx == 0) {

              // If we are at cterm then loop over all the cterm mods and try them with gap
              if (specIdx == spec.size() - 1) {
                const set<float> & setCtermMods = m_modPenaltyMatrix->getCtermMods();
                if (DEBUG_ALIGN_CTERM) DEBUG_VAR(setCtermMods.size());
		
                if (DEBUG_ALIGN_CTERM) DEBUG_MSG("Try Cterm Match");

                set<float>::const_iterator itrs = setCtermMods.begin();
                set<float>::const_iterator itrsEnd = setCtermMods.end();
                for ( ; itrs != itrsEnd; itrs++) {
                  float ctermMod = *itrs;
                  float newSpecGapMass = specGapMass - ctermMod;
                  float newDeltaMasses = deltaMasses - ctermMod;

                  if (DEBUG_ALIGN_CTERM) DEBUG_TRACE;
                  if (DEBUG_ALIGN_CTERM) DEBUG_VAR(ctermMod);
                  if (DEBUG_ALIGN_CTERM) DEBUG_VAR(newSpecGapMass);

                  if (dbAAMinMass < newSpecGapMass &&
                      newSpecGapMass > dbGapMass + m_minMod &&
                      newSpecGapMass < dbGapMass + m_maxMod) {

                    int specGapLength = (int)round(newSpecGapMass);
                    int ntermOrCterm = 0;  // Can't allow double cterm mods
                    if (fabs(ctermMod) <= pmTolerance && deltaMasses < 57.0) {
                      ntermOrCterm = 2;  // However these masses are not TRUE cterm mods
                                         // They are parent mass adjustments
                    }
                    float ctermPenalty = getGapPenalty(specGapLength, 
                                                       dbGapString, 
                                                       avgPeakIntensity, 
                                                       ntermOrCterm,
                                                       false);
                    float scanSpecificPenalty = getScanSpecificPenalty(scanSpecificPenalties,
                                                                       newDeltaMasses,
                                                                       avgPeakIntensity);
                    if (DEBUG_ALIGN_GAP) DEBUG_VAR(scanSpecificPenalty);
                    ctermPenalty = max(ctermPenalty, scanSpecificPenalty); // remember penalties are neg
                    if (DEBUG_ALIGN_GAP) DEBUG_VAR(ctermPenalty);


                    if (DEBUG_ALIGN_CTERM) DEBUG_VAR(specGapLength);
                    if (DEBUG_ALIGN_CTERM) DEBUG_VAR(ntermOrCterm);
                    if (DEBUG_ALIGN_CTERM) DEBUG_VAR(ctermPenalty);

                    if (ctermPenalty > bestGapPenalty) {
                      if (ctermPenalty == 0.0) {
                        matchString = strAA + "[" + massToString(ctermMod)+ "]";
                      } else {
                        matchString = "(" + strAA + "," +
                            massToString(deltaMasses - ctermMod) + ")" +
                            "[" + massToString(ctermMod) + "]";
                      }
                      bestGapPenalty = ctermPenalty;
                      if (DEBUG_ALIGN_CTERM) DEBUG_VAR(matchString);
                      if (DEBUG_ALIGN_CTERM) DEBUG_VAR(bestGapPenalty);
                    } // if (ctermPenalty > bestGapPenalty) {
                  } // if (dbAAMinMass < newSpecGapMass &&
                } // for ( ; itrs != itrsEnd; itrs++) {
              } // if (specIdx == spec.size() - 1) {

              if (DEBUG_ALIGN_GAP) DEBUG_VAR(bestGapPenalty);
              float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                spec[specIdx][1] + bestGapPenalty;

              if (DEBUG_ALIGN_GAP) DEBUG_VAR(score);
              if (bestGapPenalty == 0.0 || specIdx == spec.size() - 1) {
                // If there is no gap penalty, then there was no mod
                //   In this case, compute cleavage penalties
                //   No cleavage penalty if a mod blocks a cleavage
                //bool reverse = (bestGapPenalty != 0.0 && specIdx == spec.size() - 1);
                bool reverse = false;
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(reverse);
                score += computeCleavagePenalty(dbSeq,
                                                 predDbSpecIdx,
                                                 dbSpecIdx,
                                                 predSpecIdx == 0,
                                                 specIdx == spec.size() - 1,
                                                 reverse,
                                                 DEBUG_ALIGN_GAP) * avgPeakIntensity;
              }

              if (DEBUG_ALIGN_GAP) DEBUG_VAR(score);
              if (score > matchMatrix[dbSpecIdx][specIdx]) {
                matchMatrix[dbSpecIdx][specIdx] = score;
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                matchType[dbSpecIdx][specIdx] = matchString;
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
              }
            } // } else if (dbAAMinMass < specGapMass &&

          } // predSpecIdx

        } // predDbSpecIdx

      } // specIdxS

    } // dbSpecIdx

    if (DEBUG_ALIGN3) DEBUG_TRACE;
    if (DEBUG_ALIGN3) debugMatrices(spec, dbSeq, matchMatrix, matchType);
    if (DEBUG_ALIGN3) DEBUG_TRACE;
       
    //---------------------------------
    // Find best scoring match
    //---------------------------------
    float bestScore = -(float)INT_MAX;
    pair<int,int> bestMatchPtr;

    if (DEBUG_ALIGN4) DEBUG_VAR(enforceEndpeaks);
    if (enforceEndpeaks) {
      for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
        if (matchMatrix[dbSpecIdx][spec.size()-1] > bestScore) {
          bestScore = matchMatrix[dbSpecIdx][spec.size()-1];
          bestMatchPtr = make_pair<int,int>(dbSpecIdx,spec.size()-1);
          if (DEBUG_ALIGN4) DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);
        }
      }
    } else {
      for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
        for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
          if (matchMatrix[dbSpecIdx][specIdx] > bestScore) {
            bestScore = matchMatrix[dbSpecIdx][specIdx];
            bestMatchPtr = make_pair<int,int>(dbSpecIdx,specIdx);
            if (DEBUG_ALIGN4) DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);
          }
        }
      }
    }
    
    if (DEBUG_ALIGN4) DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);

    //---------------------------------
    // Create PSM for the match
    //---------------------------------
    if (bestScore != -(float)INT_MAX) {

      // Create a PSM for this match      
      psmPtr p(new PeptideSpectrumMatch);
      p->m_score = bestScore;
      p->m_dbIndex = dbIndex;
      p->m_matchOrientation = matchOrientation;
      p->m_spectrum = &spec;

      int lastDbIndex = bestMatchPtr.first;
      
      if (DEBUG_ALIGN4) {
        DEBUG_VAR(p->m_score);
        DEBUG_VAR(p->m_matchOrientation);
        DEBUG_VAR(bestMatchPtr.first);
        DEBUG_VAR(bestMatchPtr.second);
        DEBUG_VAR(matchType[bestMatchPtr.first][bestMatchPtr.second]);
        float score = matchMatrix[bestMatchPtr.first][bestMatchPtr.second];
        DEBUG_VAR(score);
        float peakMass = spec[bestMatchPtr.second][0];
        DEBUG_VAR(peakMass);
        float peakValue = spec[bestMatchPtr.second][1];
        DEBUG_VAR(peakValue);
      }

      p->m_annotation.clear();
      // Check for gap at end
      if (bestMatchPtr.second != spec.size() - 1) {
        float endMass = spec[spec.size() - 1][0] - spec[bestMatchPtr.second][0];
        p->m_annotation = "[" + massToString(endMass) + "]";
      }
      if (DEBUG_ALIGN4) DEBUG_VAR(p->m_annotation);

      pair<int,int> prevMatchPtr = bestMatchPtr;
      TwoValues<int> prevIndices;
      prevIndices[0] = -1;
      prevIndices[1] = -1;
      bool wasGap = false;
      while (bestMatchPtr.second != -1) {
        TwoValues<int> tmpIndices;
        tmpIndices[0] = bestMatchPtr.second;
        tmpIndices[1] = bestMatchPtr.first;
        p->m_matchedPeaks.push_back(tmpIndices);

        string matchAnnotation = matchType[bestMatchPtr.first][bestMatchPtr.second];

        if (matchAnnotation != "?") {
          p->m_annotation = matchAnnotation + p->m_annotation;
          prevIndices = tmpIndices;
	 }
        if (DEBUG_ALIGN4) DEBUG_VAR(p->m_annotation);
        prevMatchPtr = bestMatchPtr;
        bestMatchPtr = matchPtr[bestMatchPtr.first][bestMatchPtr.second];

        if (DEBUG_ALIGN4) {
          DEBUG_VAR(bestMatchPtr.first);
          DEBUG_VAR(bestMatchPtr.second);
          if (bestMatchPtr.first != -1 || bestMatchPtr.second != -1) {
            DEBUG_VAR(matchType[bestMatchPtr.first][bestMatchPtr.second]);

            float score = matchMatrix[bestMatchPtr.first][bestMatchPtr.second];
            DEBUG_VAR(score);
            float peakMass = spec[bestMatchPtr.second][0];
            DEBUG_VAR(peakMass);
            float peakValue = spec[bestMatchPtr.second][1];
            DEBUG_VAR(peakValue);
          }
        }

      }
      int firstDbIndex = prevMatchPtr.first;
      p->m_startMass = dbSpec[firstDbIndex][0];

      // Check for gap at beginning
      if (prevMatchPtr.second != 0) {
        p->m_annotation = "[" + massToString(spec[prevMatchPtr.second][0]) + "]" + p->m_annotation;
        p->m_startMass -= spec[prevMatchPtr.second][0]; // Adjust start mass backward to account for gap
      }
      if (DEBUG_ALIGN4) DEBUG_VAR(p->m_annotation);

      if (DEBUG_ALIGN4) DEBUG_VAR(firstDbIndex);
      if (DEBUG_ALIGN4) DEBUG_VAR(lastDbIndex);

      if (lastDbIndex - firstDbIndex > 2) {

        if (DEBUG_ALIGN4) DEBUG_VAR(p->m_annotation);
        if (DEBUG_ALIGN4) DEBUG_VAR(p->m_matchedPeaks.size());
        reverse(p->m_matchedPeaks.begin(),p->m_matchedPeaks.end());

        spec.psmList.push_back(p);
      } // if (lastDbIndex - firstDbIndex > 2) {

    }  // if (bestScore != -(float)INT_MAX)

    if (DEBUG_CACHE) DEBUG_VAR(m_cacheHits);
    if (DEBUG_CACHE) DEBUG_VAR(m_cacheMisses);
    m_cacheHitsTotal += m_cacheHits;
    m_cacheMissesTotal += m_cacheMisses;

    if (DEBUG_CACHE) DEBUG_VAR(m_cacheHitsTotal);
    if (DEBUG_CACHE) DEBUG_VAR(m_cacheMissesTotal);

    return;
  }

  //--------------------------------------------------------------------------
  float AlignmentPenaltyBased::getPeptideMass(string & stringPeptide)
  {
    float totalMass = 0.0;
    for (int i = 0; i < stringPeptide.length(); i++) {
      string strAA;
      strAA.append(1, stringPeptide[i]);
      float mass = m_modPenaltyMatrix->getMass(strAA) * 0.9995;
      totalMass += mass;
    }
    return totalMass;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::computeAllGapAnnotations(string & stringAnnotationIn,
                                                       string & stringAnnotationOut)
  {
    if (DEBUG_GAP_ANNO) DEBUG_VAR(stringAnnotationIn);
    psmPtr psmTemp(new PeptideSpectrumMatch);
    psmTemp->m_annotation = stringAnnotationIn;
    
    vector<float> modifications;
    vector<unsigned int> positions;
    vector<unsigned int> lengths;
    psmTemp->getModificationsAndPositions(modifications, positions, lengths);
    
    string cleanAnnotation = getCleanAnnotation(psmTemp->m_annotation);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(cleanAnnotation);
    
    vector<float> newModifications;
    vector<unsigned int> newPositions;
    vector<unsigned int> newLengths;
    psmPtr psmTemp2(new PeptideSpectrumMatch);
    psmTemp2->m_annotation;

    bool hadNterm= false;
    for (int iMod = 0; iMod < modifications.size(); iMod++) {
      if (DEBUG_GAP_ANNO) DEBUG_VAR(positions[iMod]);
      if (DEBUG_GAP_ANNO) DEBUG_VAR(lengths[iMod]);
      // Only concerned with gaps with > minimum mass
      if (lengths[iMod] < 2 || fabs(modifications[iMod]) < MINIMUM_GAP_MASS_TO_SPLIT) {
        if (positions[iMod] == 0 && fabs(modifications[iMod]) != 1.0) {
          hadNterm = true;
        }
        newModifications.push_back(modifications[iMod]);
        newPositions.push_back(positions[iMod]);
        newLengths.push_back(lengths[iMod]);
        continue;
      }
      string gapString = cleanAnnotation.substr(positions[iMod] - lengths[iMod], lengths[iMod]);
      gapString = "(" + gapString + "," + massToString(modifications[iMod]) + ")";
      if (DEBUG_GAP_ANNO) DEBUG_VAR(gapString);
      string cleanGapString = getCleanAnnotation(gapString);
      if (DEBUG_GAP_ANNO) DEBUG_VAR(cleanGapString);
      
      bool nterm = !hadNterm && positions[iMod] - lengths[iMod] == 0 ? true : false;
      if (DEBUG_GAP_ANNO) DEBUG_VAR(nterm);
      bool knownonly = hadNterm && positions[iMod] - lengths[iMod] == 0 ? true : false;
      if (DEBUG_GAP_ANNO) DEBUG_VAR(knownonly);
      bool cterm = positions[iMod] == cleanAnnotation.length() ? true : false;
      if (DEBUG_GAP_ANNO) DEBUG_VAR(nterm);

      string annotation;
      float aaMass = getPeptideMass(cleanGapString);
      if (DEBUG_GAP_ANNO) DEBUG_VAR(aaMass);
      if (DEBUG_GAP_ANNO) DEBUG_VAR(modifications[iMod]);
      getGapAnnotation(aaMass + modifications[iMod], cleanGapString, annotation, nterm ? 1 : 0, knownonly);
      if (DEBUG_GAP_ANNO) DEBUG_VAR(annotation);
      psmTemp2->m_annotation = annotation;

      vector<float> modifications2;
      vector<unsigned int> positions2;
      vector<unsigned int> lengths2;
      psmTemp2->getModificationsAndPositions(modifications2, positions2, lengths2);
      if (DEBUG_GAP_ANNO) DEBUG_VAR(modifications2.size());
      int modStartUnknown = 1;
      float unusedModMass = 0.0;
      for (int iMod2 = 0; iMod2 < modifications2.size(); iMod2++) {

        int truePosition = positions2[iMod2] + positions[iMod] - lengths[iMod];   // This is a "sub-mod" of the original

        string strAA("X");
        if (DEBUG_GAP_ANNO) DEBUG_VAR(positions2[iMod2]-1);
        strAA[0] = cleanGapString[positions2[iMod2]-1];
        if (DEBUG_GAP_ANNO) DEBUG_VAR(strAA);
        if (m_modPenaltyMatrix->isKnown(strAA, modifications2[iMod2]) ||
            (nterm && m_modPenaltyMatrix->isNterm(strAA, modifications2[iMod2]))) {
          // Push back a mod for all the previous unknown mass
          if (DEBUG_GAP_ANNO) DEBUG_VAR(unusedModMass);
          if (unusedModMass != 0.0) {
            newModifications.push_back(unusedModMass);
            newPositions.push_back(truePosition - 1);
            newLengths.push_back(positions2[iMod2] - modStartUnknown);
            if (DEBUG_GAP_ANNO) DEBUG_VAR(positions2[iMod2] - modStartUnknown);
            unusedModMass = 0.0;
          }
          // Push back the known mod
          newModifications.push_back(modifications2[iMod2]);
          newPositions.push_back(truePosition);
          newLengths.push_back(lengths2[iMod2]);
          if (DEBUG_GAP_ANNO) DEBUG_VAR(modifications2[iMod2]);
          if (DEBUG_GAP_ANNO) DEBUG_VAR(truePosition);
          if (DEBUG_GAP_ANNO) DEBUG_VAR(lengths2[iMod2]);
          modStartUnknown = positions2[iMod2] + 1;
          if (DEBUG_GAP_ANNO) DEBUG_VAR(modStartUnknown);
        } else {
          // If this is the last mod then push back a mod for all the leftovers
          if (iMod2 == modifications2.size() - 1) {
            if (DEBUG_GAP_ANNO) DEBUG_VAR(positions2[iMod2]);
            if (DEBUG_GAP_ANNO) DEBUG_VAR(unusedModMass);
            newModifications.push_back(modifications2[iMod2] + unusedModMass);
            if (DEBUG_GAP_ANNO) DEBUG_VAR(positions[iMod]);
            newPositions.push_back(positions[iMod]);
            if (DEBUG_GAP_ANNO) DEBUG_VAR(modStartUnknown);
            newLengths.push_back(cleanGapString.length() - modStartUnknown + 1);
            if (DEBUG_GAP_ANNO) DEBUG_VAR(cleanGapString.length() - modStartUnknown + 1);
          } else {
            // Collect up the mass
            unusedModMass += modifications2[iMod2];
            if (DEBUG_GAP_ANNO) DEBUG_VAR(unusedModMass);
          }
        }
      }
      if (DEBUG_GAP_ANNO) DEBUG_VAR(modifications2.size());
      for (int iMod2 = 0; iMod2 < modifications2.size(); iMod2++) {
        if (DEBUG_GAP_ANNO) DEBUG_VAR(newModifications[iMod2]);
        if (DEBUG_GAP_ANNO) DEBUG_VAR(newPositions[iMod2]);
        if (DEBUG_GAP_ANNO) DEBUG_VAR(newLengths[iMod2]);
      }
    }
    
    if (newModifications.size() == 0) {
      stringAnnotationOut = stringAnnotationIn;
    } else {
      PeptideSpectrumMatch::makeAnnotationFromData(cleanAnnotation, 
                             newModifications, 
                             newPositions, 
                             newLengths, 
                             stringAnnotationOut);

    }
    if (DEBUG_GAP_ANNO) DEBUG_VAR(stringAnnotationOut);
    
    return;
  }

  //--------------------------------------------------------------------------
  float AlignmentPenaltyBased::computeCleavagePenalty(string dbSeq,
                                                      int firstDbIndex,
                                                      int lastDbIndex,
                                                      bool isPepStart,
                                                      bool isPepEnd,
                                                      bool reverse,
                                                      bool debugFlag)
  {
    if (debugFlag) DEBUG_VAR(firstDbIndex);
    if (debugFlag) DEBUG_VAR(lastDbIndex - 1);
    if (debugFlag) DEBUG_VAR(isPepStart);
    if (debugFlag) DEBUG_VAR(isPepEnd);
    float cleavagePenalty = 0.0;
    if (debugFlag) DEBUG_MSG(dbSeq[firstDbIndex] << "  " <<
                              dbSeq[lastDbIndex - 1] << "  " <<
                              isPepStart << "  " <<
                              isPepEnd);
    if (isPepStart && firstDbIndex != 0) {
      // Check for cleavage penalty BEFORE the first peak match to DB
      char prevAA = dbSeq[firstDbIndex - 1];
      if (debugFlag) DEBUG_VAR(prevAA);
      cleavagePenalty += m_modPenaltyMatrix->getCleavagePenalty(prevAA, 0);
      if (debugFlag) DEBUG_VAR(cleavagePenalty);
    }
    if (debugFlag) DEBUG_VAR(cleavagePenalty);

    // Check for cleavage penalty AT the last peak match to DB
    if (isPepEnd) {
      char endAA = dbSeq[lastDbIndex - 1];
      if (debugFlag) DEBUG_VAR(endAA);
      cleavagePenalty += m_modPenaltyMatrix->getCleavagePenalty(endAA, 2);
      //    The reverse flag is to reverse the "sense" of the last AA penalty
      //    So if you would have had a penalty don't, and vice-versa
      //    This is for and cterm K,42 or similar
      //    So we pretend like it is NOT whatever it was
      if (reverse) {
        if (cleavagePenalty == 0.0) {
          endAA = 'X';
          cleavagePenalty += m_modPenaltyMatrix->getCleavagePenalty(endAA, 2);
        } else {
          cleavagePenalty = 0.0;
        }
      }
      if (debugFlag) DEBUG_VAR(cleavagePenalty);
    }
    if (debugFlag) DEBUG_VAR(cleavagePenalty);

    // Check for all internal cleavage penalties
    for (int c = firstDbIndex + (isPepStart ? 1 : 0);
        c < lastDbIndex - (isPepEnd ? 1 : 0);
        c++) {
      char internalAA = dbSeq[c];
      if (debugFlag) DEBUG_VAR(internalAA);
      if (debugFlag) DEBUG_VAR(internalAA);
      cleavagePenalty += m_modPenaltyMatrix->getCleavagePenalty(internalAA, 1);
      if (debugFlag) DEBUG_VAR(cleavagePenalty);
    }
    if (debugFlag) DEBUG_VAR(cleavagePenalty);

    return cleavagePenalty;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::getGapAnnotation(float    specGapLengthFloat,
                                               string & aaString,
                                               string & stringAnnotation,
                                               int      ntermOrCterm,
                                               bool     knownModsOnly)
  {
    if (DEBUG_GAP_ANNO) DEBUG_VAR(specGapLengthFloat);
    int specGapLength = (int)round(specGapLengthFloat);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(specGapLength);
    if (aaString.empty()) {
      return;
    }

    string strAA;
    strAA.append(1, aaString[0]);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(strAA);
    bool useNterm = (ntermOrCterm == 1) && m_modPenaltyMatrix->existsNterm(strAA);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(useNterm);

    vector<float> scoresSoFar;
    if (useNterm) {
      scoresSoFar = m_mapGapVectors[m_modPenaltyMatrix->ntermString() + strAA];
    } else if (knownModsOnly) {
      scoresSoFar = m_mapGapVectors["=" + strAA];
    } else {
      scoresSoFar = m_mapGapVectors[strAA];
    }

    int aaMass = (int)m_modPenaltyMatrix->getMass(strAA);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(aaMass);

    vector<string> vecAnnotations1(m_maxSpecGap);
    vecAnnotations1[aaMass] = strAA;
    vector<string> vecAnnotations2(m_maxSpecGap);

    // Initialize first annotation vector
    for (int i = 57; i < aaMass; i++) {
      vecAnnotations1[i] = makeAnnotation(strAA, i - aaMass);
    }
    for (int i = aaMass + 1; i < m_maxSpecGap; i++) {
      vecAnnotations1[i] = makeAnnotation(strAA, i - aaMass);
    }

#if 0
    for (int i = 0; i < m_maxSpecGap; i++) {
      DEBUG_MSG(i << "  " << vecAnnotations1[i]);
    }
#endif

    vector<string> & prevAnno = vecAnnotations1;
    vector<string> & newAnno = vecAnnotations2;
    // Add annotations for each additional character
    for (int iChar = 1; iChar < aaString.length(); iChar++) {
      string strNewAA;
      strNewAA.append(1, aaString[iChar]);
      //if (DEBUG_GAP_ANNO) DEBUG_VAR(strNewAA);
      int newAAMass = (int)m_modPenaltyMatrix->getMass(strNewAA);
      //if (DEBUG_GAP_ANNO) DEBUG_VAR(newAAMass);
      vector<float> scoresCombined(m_maxSpecGap);
      initializeGapVector(scoresCombined);
      combineVectorsWithAnnotation(scoresSoFar, m_mapGapVectors[strNewAA], scoresCombined,
                                   vecAnnotations1, vecAnnotations2, strNewAA, newAAMass);
      scoresSoFar = scoresCombined;

      // swap the annotation vectors to avoid copies
      vector<string> & tempAnno = prevAnno;
      prevAnno = newAnno;
      newAnno = tempAnno;
    }

#if 0
    for (int i = 0; i < prevAnno.size(); i++) {
      DEBUG_MSG(i << "  " << prevAnno[i]);
    }
#endif

    if (DEBUG_GAP_ANNO) DEBUG_VAR(specGapLength);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(prevAnno.size());
    if (specGapLength > prevAnno.size()) {
      WARN_MSG("Gap size too large!");
      DEBUG_VAR(specGapLength);
      DEBUG_VAR(prevAnno.size());
      // This should never happen.. but just in case
      stringAnnotation = "(" + aaString +",x.xxx)";
      return;
    }
    stringAnnotation = prevAnno[specGapLength];
    return;
  }

  //--------------------------------------------------------------------------
  float AlignmentPenaltyBased::getGapPenalty(int      specGapLength,
                                             string   aaString,
                                             float    avgPeakIntensity,
                                             int      ntermOrCterm,
                                             bool     knownModsOnly)
  {
    bool DEBUG_GAP_PENALTY = false;
    if (aaString == "SLHV") DEBUG_GAP_PENALTY = true;

    float returnValue = 1.0;  // This gets set in the "fast compute" version
    float returnValue2 = 1.0;  // This gets set in the "fast compute" version
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(specGapLength);
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(aaString);
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(avgPeakIntensity);
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(ntermOrCterm);
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(knownModsOnly);
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(aaString[0]);
    string strAA = "X";
    strAA[0] = aaString[0];
    bool useNterm = (ntermOrCterm == 1) && m_modPenaltyMatrix->existsNterm(strAA);
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(useNterm);
    strAA[0] = aaString[aaString.length() - 1];
    bool useCterm = (ntermOrCterm == 2) && m_modPenaltyMatrix->existsCterm(strAA);
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(useCterm);

    if (aaString.empty()) {
      WARN_MSG("Empty string passed to getGapPenalty!")
      return -(float)INT_MAX;
    }

    // If the length of the string is small enough then we have the answer in map
    string sortedString = aaString;
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(sortedString);
    if (aaString.length() <= m_kMer) {
      
      if (useNterm) {
        sortedString = m_modPenaltyMatrix->ntermString() + sortedString;
      }
      if (!useNterm && !useCterm && knownModsOnly) {
        sortedString = "=" + sortedString;
      }
      if (useCterm) {
        sortedString = sortedString + m_modPenaltyMatrix->ctermString() ;
      }
      sortStringLetters(sortedString);
      if (DEBUG_GAP_PENALTY) DEBUG_VAR(sortedString);

      // See if this string is in the map
      if (m_mapGapVectors.find(sortedString) == m_mapGapVectors.end()) {
        return -(float)INT_MAX;
      }

#if 0
      if (DEBUG_GAP_PENALTY && knownModsOnly) {
        vector<float> temp = m_mapGapVectors[sortedString];
        DEBUG_VAR(left);
        for (int i = 0; i < temp.size(); i++) {
          DEBUG_MSG(i << "  " << temp[i])
        }
      }
#endif
      if (DEBUG_GAP_PENALTY) DEBUG_VAR(m_mapGapVectors[sortedString][specGapLength]);
      return m_mapGapVectors[sortedString][specGapLength] * avgPeakIntensity;
    }

    if (useNterm) {
      sortedString = m_modPenaltyMatrix->ntermString() + sortedString;
    } else if (knownModsOnly) {
      sortedString = "=" + sortedString;
    }
    
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(sortedString);

    if (USE_CACHE) {
      if (isCached(sortedString)) {
        //if (DEBUG_CACHE) DEBUG_MSG("Retrieved answer for [" << sortedString << "] from the cache");
        m_cacheHits++;
        return m_mapCachedGapVectors[sortedString][specGapLength] * avgPeakIntensity;
      }
      //if (DEBUG_CACHE) DEBUG_MSG("Answer for [" << sortedString << "] not in cache");
      m_cacheMisses++;
    }

    string left;
    if (useNterm || knownModsOnly) {
      left = sortedString.substr(0, m_kMer + 1);  // Need to pull of "<" also
    } else {
      left = sortedString.substr(0, m_kMer);
    }
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(left);
    sortStringLetters(left);
    vector<float> scoresSoFar = m_mapGapVectors[left];

#if 0
    if (DEBUG_GAP_PENALTY) {
      DEBUG_VAR(left);
      for (int i = 0; i < scoresSoFar.size(); i++) {
        DEBUG_MSG(i << "  " << scoresSoFar[i])
      }
    }
#endif

    string remaining = sortedString.substr(m_kMer);
    if (useNterm || knownModsOnly) {
      remaining = sortedString.substr(m_kMer + 1);
    }
    if (DEBUG_GAP_PENALTY) DEBUG_VAR(remaining);
    while (!remaining.empty()) {
      string left = remaining.substr(0, m_kMer);
      //DEBUG_VAR(left)
      sortStringLetters(left);
      if (DEBUG_GAP_PENALTY) DEBUG_VAR(left);

      if (!USE_CACHE && (remaining.size() <= m_kMer)) {
#if 0
        if (DEBUG_GAP_PENALTY) {
          DEBUG_VAR(left);
          for (int i = 0; i < m_mapGapVectors[left].size(); i++) {
            DEBUG_MSG(i << "  " << m_mapGapVectors[left][i])
          }
        }
#endif
        if (remaining.size() <= m_kMer && useCterm) {
          left = left + m_modPenaltyMatrix->ctermString();
        }
        if (remaining.size() <= m_kMer && knownModsOnly) {
          left = "=" + left;
        }
        if (DEBUG_GAP_PENALTY) DEBUG_VAR(left);
        
        returnValue = combineVectorsSingleValue(scoresSoFar, m_mapGapVectors[left], specGapLength);
        if (DEBUG_GAP_PENALTY) DEBUG_VAR(returnValue);
      } else {
        if (DEBUG_GAP_PENALTY) DEBUG_TRACE;
        vector<float> scoresCombined(m_maxSpecGap);
        initializeGapVector(scoresCombined);
        combineVectors(scoresSoFar, m_mapGapVectors[left], scoresCombined);
        scoresSoFar = scoresCombined;
      }

      if (remaining.size() <= m_kMer) {
        returnValue2 = scoresSoFar[specGapLength];
        if (DEBUG_GAP_PENALTY) DEBUG_VAR(returnValue2);
      }

      if (remaining.size() > m_kMer) {
        remaining = remaining.substr(m_kMer);
      } else {
        remaining = "";
      }
      if (DEBUG_GAP_PENALTY) DEBUG_VAR(remaining);
    }

    // Cache this string in case user needs it again
    if (USE_CACHE) {
      m_mapCachedGapVectors[sortedString] = scoresSoFar;
    }

    // Return the value
    if (!USE_CACHE) {
      if (DEBUG_GAP_PENALTY) DEBUG_VAR(returnValue);
      return returnValue * avgPeakIntensity;
    }

    return scoresSoFar[specGapLength] * avgPeakIntensity;
  }


  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::combineVectors(vector<float> & scores1,
                                             vector<float> & scores2,
                                             vector<float> & outputScores)
  {
    float unkPenalty = m_modPenaltyMatrix->getUnknownPenalty(m_penaltyAlpha); // multiply by alpha
    vector<float>::iterator itrOut = outputScores.begin();
    vector<float>::iterator itr1 = scores1.begin();
    for (int i = 0; i < scores1.size(); i++, itrOut++, itr1++) {
      vector<float>::iterator itrOut2 = itrOut;
      vector<float>::iterator itr2 = scores2.begin();
      for (int j = 0; j < scores2.size() - i; j++, itrOut2++, itr2++) {
        float sum = *itr1 + *itr2;
        if (*itrOut2 < sum) {
          *itrOut2 = sum;
        }
        // The maximum possible penalty is always the Unknown Penalty
        //     Don't change the unset values
        if (*itrOut2 < unkPenalty && *itrOut2 > -100000.0 ) {
          *itrOut2 = unkPenalty;
        }

      } // for (int j = 0; j < scores2.size(); j++) {
    } // for (int i = 0; i < scores1.size(); i++) {

    return;
  }
  
  //--------------------------------------------------------------------------
  float AlignmentPenaltyBased::combineVectorsSingleValue(vector<float> & scores1,
                                                         vector<float> & scores2,
                                                         int             specGapLength)
  {
    float returnScore = -(float)INT_MAX;
    float unkPenalty = m_modPenaltyMatrix->getUnknownPenalty(m_penaltyAlpha); // multiply by alpha
    vector<float>::iterator itr1 = scores1.begin();
    vector<float>::iterator itr2 = scores2.begin() + specGapLength;
    for (int i = 0; i <= specGapLength; i++, itr1++, itr2--) {
      float sum = *itr1 + *itr2;
      if (returnScore < sum) {
        returnScore = sum;
      }
      // The maximum possible penalty is always the Unknown Penalty
      //     Don't change the unset values
      if (returnScore < unkPenalty) {
        returnScore = unkPenalty;
      }
    } // for (int i = 0; i < scores1.size(); i++) {

    return returnScore;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::combineVectorsWithAnnotation(vector<float>  & scores1,
                                                           vector<float>  & scores2,
                                                           vector<float>  & outputScores,
                                                           vector<string> & annotationsOld,
                                                           vector<string> & annotationsNew,
                                                           string         & newAA,
                                                           int            & newAAMass)
  {
    vector<char> vecFlags(outputScores.size());
    vector<char>::iterator itrFlags = vecFlags.begin();
    vector<string>::iterator itrAnnoOld = annotationsOld.begin();
    vector<string>::iterator itrAnnoNew = annotationsNew.begin();
    vector<float>::iterator itrOut = outputScores.begin();
    vector<float>::iterator itr1 = scores1.begin();
    for (int i = 0; i < scores1.size(); i++, itrOut++, itr1++, itrFlags++, itrAnnoOld++, itrAnnoNew++) {
      vector<float>::iterator itrOut2 = itrOut;
      vector<float>::iterator itr2 = scores2.begin();
      vector<char>::iterator itrFlags2 = itrFlags;
      vector<string>::iterator itrAnnoNew2 = itrAnnoNew;
      for (int j = 0; j < scores2.size() - i; j++, itrOut2++, itr2++, itrFlags2++, itrAnnoNew2++) {
        float sum = *itr1 + *itr2;
        if ((*itrFlags2 == 0) || (*itrOut2 < sum)) {
          *itrFlags2 = 1;   //vecFlags[i + j] = 1;
          *itrOut2 = sum;   //outputScores[i + j] = scores1[i] + scores2[j];
          if (*itr2 == 0) {  // if scores2[j] == 0
            *itrAnnoNew2 = *itrAnnoOld + newAA; //annotationsNew[i + j] = annotationsOld[i] + newAA;
          } else {
            //annotationsNew[i + j] = annotationsOld[i] + makeAnnotation(newAA, j - newAAMass);
            *itrAnnoNew2 = *itrAnnoOld + makeAnnotation(newAA, j - newAAMass);  
          }

          // LARS: Possibly do something about unknowns here (maximum penalty)
        }
      } // for (int j = 0; j < scores2.size(); j++) {

    } // for (int i = 0; i < scores1.size(); i++) {

    return;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::createGapVectors(void)
  {
    // Set this so we know which vectors are pre-created
    //   kMer = 4 means all combinations of 4 AA's (sorted)
    m_kMer = KMER_LENGTH;

    vector<string> aaVec;
    m_modPenaltyMatrix->getAminoAcids(aaVec);
    if (DEBUG_GAP_VECTOR) DEBUG_VAR(aaVec.size());
    sort(aaVec.begin(), aaVec.end());

    //---------------------------------------------------------    
    // Create all the single AA gap vectors
    //---------------------------------------------------------    
    for (int a1 = 0; a1 < aaVec.size(); a1++) {
      string strAA1 = aaVec[a1];
      //DEBUG_VAR(strAA1);
      vector<float> scores1(m_maxSpecGap);
      initializeGapVector(scores1);
      fillVectorSingle(scores1, strAA1, 0, false);
      m_mapGapVectors[strAA1] = scores1;
    }

    if (DEBUG_GAP_VECTOR) DEBUG_TRACE;

    //---------------------------------------------------------    
    // Create all the AA gap vectors
    //---------------------------------------------------------    
    for (int a1 = 0; a1 < aaVec.size(); a1++) {

      string strAA1 = aaVec[a1];
      if (DEBUG_GAP_VECTOR) DEBUG_VAR(strAA1);

      for (int a2 = 0; a2 < aaVec.size() && m_kMer >= 2; a2++) {

        string strAA2 = aaVec[a2];
        string strAACombined2 = strAA1 + strAA2;
        sortStringLetters(strAACombined2);

        // Check to see if we've already done this one
        if (m_mapGapVectors.find(strAACombined2) != m_mapGapVectors.end()) {
          //DEBUG_MSG("Already Done");
          continue;
        }
        if (DEBUG_GAP_VECTOR) DEBUG_VAR(strAACombined2);

        vector<float> scoresCombined2(m_maxSpecGap);
        initializeGapVector(scoresCombined2);

        combineVectors(m_mapGapVectors[strAA1], m_mapGapVectors[strAA2], scoresCombined2);
        m_mapGapVectors[strAACombined2] = scoresCombined2;

        for (int a3 = 0; a3 < aaVec.size() && m_kMer >= 3; a3++) {

          string strAA3 = aaVec[a3];
          string strAACombined3 = strAACombined2 + strAA3;
          sortStringLetters(strAACombined3);

          // Check to see if we've already done this one
          if (m_mapGapVectors.find(strAACombined3) != m_mapGapVectors.end()) {
            //DEBUG_MSG("Already Done");
            continue;
          }
          //DEBUG_VAR(strAACombined3);

          vector<float> scoresCombined3(m_maxSpecGap);
          initializeGapVector(scoresCombined3);
          combineVectors(m_mapGapVectors[strAACombined2], m_mapGapVectors[strAA3], scoresCombined3);
          m_mapGapVectors[strAACombined3] = scoresCombined3;

          for (int a4 = 0; a4 < aaVec.size() && m_kMer >= 4; a4++) {

            string strAA4 = aaVec[a4];
            string strAACombined4 = strAACombined3 + strAA4;
            sortStringLetters(strAACombined4);

            // Check to see if we've already done this one
            if (m_mapGapVectors.find(strAACombined4) != m_mapGapVectors.end()) {
              //DEBUG_MSG("Already Done");
              continue;
            }
            //DEBUG_VAR(strAACombined4);

            vector<float> scoresCombined4(m_maxSpecGap);
            initializeGapVector(scoresCombined4);
            combineVectors(m_mapGapVectors[strAACombined3], m_mapGapVectors[strAA4], scoresCombined4);
            m_mapGapVectors[strAACombined4] = scoresCombined4;
          }

        }  // for (int a3 = 0; a3 < aaVec.size(); a3++) {

      }  // for (int a2 = 0; a2 < aaVec.size(); a2++) {

    } // for (int a1 = 0; a1 < aaVec.size(); a1++) {

    if (DEBUG_GAP_VECTOR) DEBUG_VAR(m_mapGapVectors.size());
    
    //---------------------------------------------------------    
    // Create all the single AA known mod only gap vectors
    //---------------------------------------------------------    
    for (int a1 = 0; a1 < aaVec.size(); a1++) {
      string strAA1 = aaVec[a1];
      //DEBUG_VAR(strAA1);
      vector<float> scores1(m_maxSpecGap);
      initializeGapVector(scores1);
      fillVectorSingle(scores1, strAA1, 0, true);
      // Insert these as "=X"
      m_mapGapVectors["=" + strAA1] = scores1;
    }

    if (DEBUG_GAP_VECTOR) DEBUG_TRACE;

    //---------------------------------------------------------    
    // Create all the known mod only gap vectors
    //---------------------------------------------------------    
    for (int a1 = 0; a1 < aaVec.size(); a1++) {

      string strAA1 = "=" + aaVec[a1];
      if (DEBUG_GAP_VECTOR) DEBUG_VAR(strAA1);

      for (int a2 = 0; a2 < aaVec.size() && m_kMer >= 2; a2++) {

        string strAA2 = "=" + aaVec[a2];
        string strAACombined2 = strAA1 + aaVec[a2];
        sortStringLetters(strAACombined2);

        // Check to see if we've already done this one
        if (m_mapGapVectors.find(strAACombined2) != m_mapGapVectors.end()) {
          //DEBUG_MSG("Already Done");
          continue;
        }
        if (DEBUG_GAP_VECTOR) DEBUG_VAR(strAACombined2);

        vector<float> scoresCombined2(m_maxSpecGap);
        initializeGapVector(scoresCombined2);

        combineVectors(m_mapGapVectors[strAA1], m_mapGapVectors[strAA2], scoresCombined2);
        m_mapGapVectors[strAACombined2] = scoresCombined2;

        for (int a3 = 0; a3 < aaVec.size() && m_kMer >= 3; a3++) {

          string strAA3 = "=" + aaVec[a3];
          string strAACombined3 = strAACombined2 + aaVec[a3];
          sortStringLetters(strAACombined3);

          // Check to see if we've already done this one
          if (m_mapGapVectors.find(strAACombined3) != m_mapGapVectors.end()) {
            //DEBUG_MSG("Already Done");
            continue;
          }
          //DEBUG_VAR(strAACombined3);

          vector<float> scoresCombined3(m_maxSpecGap);
          initializeGapVector(scoresCombined3);
          combineVectors(m_mapGapVectors[strAACombined2], m_mapGapVectors[strAA3], scoresCombined3);
          m_mapGapVectors[strAACombined3] = scoresCombined3;

          for (int a4 = 0; a4 < aaVec.size() && m_kMer >= 4; a4++) {

            string strAA4 = "=" + aaVec[a4];
            string strAACombined4 = strAACombined3 + aaVec[a4];
            sortStringLetters(strAACombined4);

            // Check to see if we've already done this one
            if (m_mapGapVectors.find(strAACombined4) != m_mapGapVectors.end()) {
              //DEBUG_MSG("Already Done");
              continue;
            }
            //DEBUG_VAR(strAACombined4);

            vector<float> scoresCombined4(m_maxSpecGap);
            initializeGapVector(scoresCombined4);
            combineVectors(m_mapGapVectors[strAACombined3], m_mapGapVectors[strAA4], scoresCombined4);
            m_mapGapVectors[strAACombined4] = scoresCombined4;
          }

        }  // for (int a3 = 0; a3 < aaVec.size(); a3++) {

      }  // for (int a2 = 0; a2 < aaVec.size(); a2++) {

    } // for (int a1 = 0; a1 < aaVec.size(); a1++) {

    if (DEBUG_GAP_VECTOR) DEBUG_VAR(m_mapGapVectors.size());
    
    //---------------------------------------------------------    
    // Create all the single Nterm gap vectors
    //---------------------------------------------------------    
    const map<string, set<float> > & mapNtermMods = m_modPenaltyMatrix->getAllNtermMods();
    map<string, set<float> >::const_iterator itrm = mapNtermMods.begin();
    map<string, set<float> >::const_iterator itrmEnd = mapNtermMods.end();
    for(; itrm != itrmEnd; itrm++) {
      string strAA1 = itrm->first;
      // We are not interested in the non-AA nterm mods
      if (strAA1 == m_modPenaltyMatrix->ntermString() || 
          itrm->second.size() == 0) {
        continue;
      }
      vector<float> scores1(m_maxSpecGap);
      initializeGapVector(scores1);
      fillVectorSingle(scores1, strAA1, 1, false);
      // Insert these as "<X"
      m_mapGapVectors[m_modPenaltyMatrix->ntermString() + strAA1] = scores1;
    }

    //---------------------------------------------------------    
    // Create all the Nterm AA gap vectors
    //---------------------------------------------------------    
    itrm = mapNtermMods.begin();
    for(; itrm != itrmEnd; itrm++) {

      // We are not interested in the non-AA nterm mods
      if (itrm->first == m_modPenaltyMatrix->ntermString() || 
          itrm->second.size() == 0) {
        continue;
      }
      string strAA1 = m_modPenaltyMatrix->ntermString() + itrm->first;
      //DEBUG_VAR(strAA1);

      for (int a2 = 0; a2 < aaVec.size(); a2++) {

        string strAA2 = aaVec[a2];
        string strAACombined2 = strAA1 + strAA2;
        sortStringLetters(strAACombined2);

        // Check to see if we've already done this one
        if (m_mapGapVectors.find(strAACombined2) != m_mapGapVectors.end()) {
          //DEBUG_MSG("Already Done");
          continue;
        }
        //DEBUG_VAR(strAACombined2);

        vector<float> scoresCombined2(m_maxSpecGap);
        initializeGapVector(scoresCombined2);

        combineVectors(m_mapGapVectors[strAA1], m_mapGapVectors[strAA2], scoresCombined2);
        m_mapGapVectors[strAACombined2] = scoresCombined2;

        for (int a3 = 0; a3 < aaVec.size(); a3++) {

          string strAA3 = aaVec[a3];
          string strAACombined3 = strAACombined2 + strAA3;
          sortStringLetters(strAACombined3);

          // Check to see if we've already done this one
          if (m_mapGapVectors.find(strAACombined3) != m_mapGapVectors.end()) {
            //DEBUG_MSG("Already Done");
            continue;
          }
          //DEBUG_VAR(strAACombined3);

          vector<float> scoresCombined3(m_maxSpecGap);
          initializeGapVector(scoresCombined3);
          combineVectors(m_mapGapVectors[strAACombined2], m_mapGapVectors[strAA3], scoresCombined3);
          m_mapGapVectors[strAACombined3] = scoresCombined3;

          for (int a4 = 0; a4 < aaVec.size(); a4++) {

            string strAA4 = aaVec[a4];
            string strAACombined4 = strAACombined3 + strAA4;
            sortStringLetters(strAACombined4);

            // Check to see if we've already done this one
            if (m_mapGapVectors.find(strAACombined4) != m_mapGapVectors.end()) {
              //DEBUG_MSG("Already Done");
              continue;
            }
            //DEBUG_VAR(strAACombined4);

            vector<float> scoresCombined4(m_maxSpecGap);
            initializeGapVector(scoresCombined4);
            combineVectors(m_mapGapVectors[strAACombined3], m_mapGapVectors[strAA4], scoresCombined4);
            m_mapGapVectors[strAACombined4] = scoresCombined4;
          }

        }  // for (int a3 = 0; a3 < aaVec.size(); a3++) {

      }  // for (int a2 = 0; a2 < aaVec.size(); a2++) {

    } // for (int a1 = 0; a1 < aaVec.size(); a1++) {

    if (DEBUG_GAP_VECTOR) DEBUG_VAR(m_mapGapVectors.size());

    //---------------------------------------------------------    
    // Create all the single Cterm gap vectors
    //---------------------------------------------------------    
    const map<string, set<float> > & mapCtermMods = m_modPenaltyMatrix->getAllCtermMods();
    itrm = mapCtermMods.begin();
    itrmEnd = mapCtermMods.end();
    for(; itrm != itrmEnd; itrm++) {
      string strAA1 = itrm->first;
      // We are not interested in the non-AA cterm mods
      if (strAA1 == m_modPenaltyMatrix->ctermString() || 
          itrm->second.size() == 0) {
        continue;
      }
      vector<float> scores1(m_maxSpecGap);
      initializeGapVector(scores1);
      fillVectorSingle(scores1, strAA1, 2, false);
      // Insert these as "X>"
      m_mapGapVectors[strAA1 + m_modPenaltyMatrix->ctermString()] = scores1;
    }
    
    //---------------------------------------------------------    
    // Create all the Cterm AA gap vectors
    //---------------------------------------------------------    
    itrm = mapCtermMods.begin();
    for(; itrm != itrmEnd; itrm++) {

      // We are not interested in the non-AA nterm mods
      if (itrm->first == m_modPenaltyMatrix->ctermString() || 
        itrm->second.size() == 0) {
        continue;
      }
      string strAA1 = itrm->first + m_modPenaltyMatrix->ctermString();
      //DEBUG_VAR(strAA1);

      for (int a2 = 0; a2 < aaVec.size(); a2++) {

        string strAA2 = aaVec[a2];
        string strAACombined2 =  strAA2 + strAA1;
        sortStringLetters(strAACombined2);

        // Check to see if we've already done this one
        if (m_mapGapVectors.find(strAACombined2) != m_mapGapVectors.end()) {
          //DEBUG_MSG("Already Done");
          continue;
        }
        //DEBUG_VAR(strAACombined2);

        vector<float> scoresCombined2(m_maxSpecGap);
        initializeGapVector(scoresCombined2);

        combineVectors(m_mapGapVectors[strAA1], m_mapGapVectors[strAA2], scoresCombined2);
        m_mapGapVectors[strAACombined2] = scoresCombined2;

        for (int a3 = 0; a3 < aaVec.size(); a3++) {

          string strAA3 = aaVec[a3];
          string strAACombined3 =  strAA3 + strAACombined2;
          sortStringLetters(strAACombined3);

          // Check to see if we've already done this one
          if (m_mapGapVectors.find(strAACombined3) != m_mapGapVectors.end()) {
            //DEBUG_MSG("Already Done");
            continue;
          }
          //DEBUG_VAR(strAACombined3);

          vector<float> scoresCombined3(m_maxSpecGap);
          initializeGapVector(scoresCombined3);
          combineVectors(m_mapGapVectors[strAACombined2], m_mapGapVectors[strAA3], scoresCombined3);
          m_mapGapVectors[strAACombined3] = scoresCombined3;

          for (int a4 = 0; a4 < aaVec.size(); a4++) {

            string strAA4 = aaVec[a4];
            string strAACombined4 =  strAA4 + strAACombined3;
            sortStringLetters(strAACombined4);

            // Check to see if we've already done this one
            if (m_mapGapVectors.find(strAACombined4) != m_mapGapVectors.end()) {
              //DEBUG_MSG("Already Done");
              continue;
            }
            //DEBUG_VAR(strAACombined4);

            vector<float> scoresCombined4(m_maxSpecGap);
            initializeGapVector(scoresCombined4);
            combineVectors(m_mapGapVectors[strAACombined3], m_mapGapVectors[strAA4], scoresCombined4);
            m_mapGapVectors[strAACombined4] = scoresCombined4;
          }

        }  // for (int a3 = 0; a3 < aaVec.size(); a3++) {

      }  // for (int a2 = 0; a2 < aaVec.size(); a2++) {

    } // for (int a1 = 0; a1 < aaVec.size(); a1++) {

    if (DEBUG_GAP_VECTOR) DEBUG_VAR(m_mapGapVectors.size());

    return;
  }


  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::fillVectorSingle(vector<float> &  scores,
                                               string           strAA,
                                               int              ntermOrCterm,
                                               bool             knownModsOnly)
  {
    //DEBUG_VAR(strAA);
    // No penalty for hitting DB mass exactly
    float aaMass = m_modPenaltyMatrix->getMass(strAA);
    //DEBUG_VAR(aaMass);
    size_t index = (size_t)aaMass;
    scores[index] = 0.0;

    std::map<float, float> penalties;
    m_modPenaltyMatrix->getPenalties(strAA, penalties, 1.0);
    std::map<float, float>::iterator itrm = penalties.begin();
    std::map<float, float>::iterator itrm_end = penalties.end();
    // Loop over all the possible mods (and their penalties)
    for (;itrm != itrm_end; itrm++) {
      float massDiff = itrm->first;
      float penalty = itrm->second;
      if (penalty == 0.0) {
        continue;
      }
      size_t indexPenalty = (size_t)(index + massDiff);
      // Penalties are negative so a "higher" penalty is better
      float newPenalty = penalty * m_penaltyAlpha;
       // Skip anything that is not known if knownModsOnly is set
      if (knownModsOnly && !m_modPenaltyMatrix->isKnown(strAA, massDiff)) {
        continue;
      }
      if (indexPenalty >= 0 &&
          indexPenalty < scores.size()) {
        scores[indexPenalty] = newPenalty;
      }
    }

    if (!knownModsOnly) {
      // Fill in all possible unknown modifications
      int minRealDelta = (int)(m_modPenaltyMatrix->getMass(strAA) - 57.0);
      int startUnk = max((int)(index - minRealDelta), (int)(index + m_minMod)); // min mod is negative
      startUnk = max(0, startUnk);  // Don't go beyond beginning of array
      int endUnk = min((int)(index + m_maxMod), (int)scores.size());
      float unkPenalty = m_modPenaltyMatrix->getUnknownPenalty(m_penaltyAlpha); // multiply by alpha
      for (int iUnk = startUnk; iUnk < endUnk; iUnk++) {
        float newPenalty = unkPenalty;
        if (newPenalty > scores[iUnk]) {
          scores[iUnk] = newPenalty;
        }
      }
    }
    
    if (ntermOrCterm == 1) {
      const set<float> & setNtermMods = m_modPenaltyMatrix->getNtermMods(strAA);
      set<float>::iterator itrs = setNtermMods.begin();
      set<float>::iterator itrs_end = setNtermMods.end();
      for(; itrs != itrs_end; itrs++) {
        size_t indexPenalty = (size_t)(index + *itrs);
        scores[indexPenalty] = m_modPenaltyMatrix->getKnownPenalty(1.0);
      }
    }
    
    if (ntermOrCterm == 2) {
      const set<float> & setCtermMods = m_modPenaltyMatrix->getCtermMods(strAA);
      set<float>::iterator itrs = setCtermMods.begin();
      set<float>::iterator itrs_end = setCtermMods.end();
      for(; itrs != itrs_end; itrs++) {
        size_t indexPenalty = (size_t)(index + *itrs);
        scores[indexPenalty] = m_modPenaltyMatrix->getKnownPenalty(1.0);
      }
    }
    
    return;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::initializeGapVector(vector<float> & scores)
  {
    // Initialize the vector
    for (int i = 0; i < scores.size(); i++) {
      scores[i] = -(float)INT_MAX;
    }
  }

  //--------------------------------------------------------------------------
  string AlignmentPenaltyBased::makeAnnotation(string & strAA, int mass)
  {
    char strMod[100];
    int iMod = (int)(fabs(mass) + 0.5) * (mass < 0 ? -1 : 1);
    sprintf(strMod, "%d", iMod);
    return "(" + strAA + "," + strMod + ")";
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::makeDbString(char * dbGapStringOut, char * dbSeq, int start, int end)
  {
    int length = end - start;
    strncpy(dbGapStringOut, dbSeq + start, length);
    dbGapStringOut[length] = '\0';
    return;
  }

  //--------------------------------------------------------------------------
  string AlignmentPenaltyBased::massToString(float mass)
  {
    char strMod[100];
    int iMod = (int)(fabs(mass) + 0.5) * (mass < 0 ? -1 : 1);
    sprintf(strMod, "%d", iMod);
    return string(strMod);
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::setStartingFlagArray(Spectrum &              spec,
                                                   Spectrum &              dbSpec,
                                                   set<float> &            startRange,
                                                   vector<vector<char> > & startFlags,
                                                   int                   & firstValidDbIndex,
                                                   float                   tolerance)
  {
    //---------------------------------
    //  INITIALIZE THE STARTING FLAGS
    //---------------------------------
    for (int i = 0; i < spec.size(); i++) {
      vector<char> newArray(dbSpec.size());
      startFlags[i] = newArray;
      for (int j = 0; j < dbSpec.size(); j++) {
	    // If the startRange set is empty then allow any start position
	    if (startRange.size() == 0) {
	      startFlags[i][j] = 1;
	    } else {
	      startFlags[i][j] = 0;
	    }
      }
    }

    firstValidDbIndex = dbSpec.size();

    if (DEBUG_RANGE) DEBUG_VAR(m_minMod);
    if (DEBUG_RANGE) DEBUG_VAR(m_maxMod);
    // Find all the valid starting points in the matrix
    set<float>::iterator itr = startRange.begin();
    set<float>::iterator itrEnd = startRange.end();
    for (; itr != itrEnd; itr++) {
      int minIdx1, maxIdx1, minIdx2, maxIdx2;
      if (DEBUG_RANGE) DEBUG_VAR(*itr);
      float minStartMass = *itr + m_minMod * 2;
      float maxStartMass = *itr + m_maxMod * 2;
      if (DEBUG_RANGE) DEBUG_VAR(minStartMass);
      if (DEBUG_RANGE) DEBUG_VAR(maxStartMass);
      float lastMass = 0;
      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
	    float massDiff = spec[specIdx][0] - lastMass;
	    if (DEBUG_RANGE) DEBUG_VAR(massDiff);
	    minStartMass += massDiff;
	    maxStartMass += massDiff;
	    if (DEBUG_RANGE) DEBUG_VAR(minStartMass);
	    if (DEBUG_RANGE) DEBUG_VAR(maxStartMass);
	    lastMass = spec[specIdx][0];
	    // Find peaks close to desired bounds
	    float T = 200.0;

	    list<int> matches1;
	    list<int> matches2;
	    dbSpec.setPeakTolerance(0);
	    dbSpec.findPeaks(minStartMass, T, &matches1);
	    dbSpec.findPeaks(maxStartMass, T, &matches2);
	    if (DEBUG_RANGE) DEBUG_VAR(matches1.size());
	    if (DEBUG_RANGE) DEBUG_VAR(matches1.front());
	    if (DEBUG_RANGE) DEBUG_VAR(matches1.back());
	    if (DEBUG_RANGE) DEBUG_VAR(matches2.size());
	    if (DEBUG_RANGE) DEBUG_VAR(matches2.front());
	    if (DEBUG_RANGE) DEBUG_VAR(matches2.back());
	    minIdx1 = (matches1.size() == 0) ? -1 : matches1.front();
	    maxIdx1 = (matches1.size() == 0) ? -1 : matches1.back();
	    minIdx2 = (matches2.size() == 0) ? -1 : matches2.front();
	    maxIdx2 = (matches2.size() == 0) ? -1 : matches2.back();

	    // If peak is not found, index could be -1.. so make sure at least 0
	    minIdx1 = max<int>(minIdx1,0);
	    minIdx2 = max<int>(minIdx2,0);
	    // Make sure min peak is INSIDE the desired range
	    while (minIdx1 < dbSpec.size() &&
			dbSpec[minIdx1][0] < minStartMass - tolerance - AAJumps::massH2O) minIdx1++;

	    // Make sure max peak is OUTSIDE (above) the desired range
	    while (minIdx2 < dbSpec.size() &&
			dbSpec[minIdx2][0] < maxStartMass + tolerance + AAJumps::massH2O) minIdx2++;

	    // Mark all peaks in range as valid starting points in matrix
	    if (DEBUG_RANGE) DEBUG_VAR(minIdx1);
	    if (DEBUG_RANGE) DEBUG_VAR(minIdx2);
	    for (int i = minIdx1; i <= minIdx2 && i < dbSpec.size(); i++) {
	      startFlags[specIdx][i] = 1;
	    }
           if (minIdx1 < firstValidDbIndex) {
             firstValidDbIndex = minIdx1;
           }
      } // for (int specIdx = 0; specIdx < spec.size(); specIdx++)
    } // for (; itr != itrEnd; itr++) {

    if (DEBUG_RANGE) {
      for (int j = 0; j < dbSpec.size(); j++) {
        cout << j << "\t";
        for (int i = 0; i < spec.size(); i++) {
	  cout << (int)startFlags[i][j] << "\t";
        }
        cout << endl;
      }
    }
    return;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::setDebugFlags(
                       bool debugRange,
                       bool debugRange2,
                       bool debugAas,
                       bool debugSpecs,
                       bool debugAlign,
                       bool debugAlign1,
                       bool debugAlign2,
                       bool debugAlign3,
                       bool debugAlign4,
                       bool debugAlignExact,
                       bool debugAlignNterm,
                       bool debugAlignCterm,
                       bool debugGap,
                       bool debugGapAnno,
                       bool debugCache,
                       bool debugAlignLarge)
  {
    DEBUG_RANGE = debugRange,
    DEBUG_RANGE2 = debugRange2,
    DEBUG_AAS = debugAas,
    DEBUG_SPECS = debugSpecs,
    DEBUG_ALIGN = debugAlign,
    DEBUG_ALIGN1 = debugAlign1,
    DEBUG_ALIGN2 = debugAlign2,
    DEBUG_ALIGN3 = debugAlign3,
    DEBUG_ALIGN4 = debugAlign4,
    DEBUG_ALIGN_EXACT = debugAlignExact,
    DEBUG_ALIGN_NTERM = debugAlignNterm,
    DEBUG_ALIGN_CTERM = debugAlignCterm,
    DEBUG_ALIGN_GAP = debugGap;
    DEBUG_GAP_ANNO = debugGapAnno;
    DEBUG_CACHE = debugCache;
    DEBUG_ALIGN_LARGE = debugAlignLarge;
    return;
  }


} // namespace specnets
