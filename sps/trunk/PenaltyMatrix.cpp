#include "PenaltyMatrix.h"

// Module Includes
#include "DelimitedTextReader.h"

// System Includes
#include <stdio.h>
#include <math.h>

#define DEBUG_PENALTY 0
#define DEBUG_AAS 0
#define DEBUG_SPECPROB_MODS 0

const unsigned int CLEAVAGE_START = 0;
const unsigned int CLEAVAGE_INTERNAL = 1;
const unsigned int CLEAVAGE_END = 2;

using namespace specnets;
using namespace std;

namespace PenaltyMatrix_const {
  std::map<float, float> emptyMap;
  string stringN("N");
  string stringD("D");
  string stringE("E");
  string stringQ("Q");
  string stringK("K");
  string stringI("I");
  string stringL("L");
  string stringNterm("<");
  string stringCterm(">");
  char charNterm = '<';
  char charCterm = '>';
  char charInternal = '|';
}
using namespace PenaltyMatrix_const;

// -------------------------------------------------------------------------
PenaltyMatrix::PenaltyMatrix(AAJumps & jumps, 
                             float resolution /* = 1.0 */, 
                             float knownModPenalty /* = 1.0 */, 
                             float unknownPenalty /* = 1.0 */, 
                             float unknownMultiplier /* = 2.0 */,
                             float minModMass /* = 100.0 */,
                             float maxModMass /* = 100.0 */) :
  m_resolution(resolution), 
  m_allSpectraAveragePeakIntensity(1.0), 
  m_knownModPenalty(-knownModPenalty),           // Penalties are negative always
  m_unknownPenalty(-unknownPenalty),             // Penalties are negative always
  m_unknownMultiplier(unknownMultiplier),
  m_minModMass(minModMass),
  m_maxModMass(maxModMass)
{
  vector<pair<char, float> > vecRefAminoAcids;
  if (jumps.getAllAArefs(vecRefAminoAcids) == 0) {
    WARN_MSG("Initializing with empty amino acids");
  }

#if DEBUG_AAS 
  for (int i = 0; i < vecRefAminoAcids.size(); i++) {
    DEBUG_MSG(vecRefAminoAcids[i].first << " = " << vecRefAminoAcids[i].second);
  }
#endif

  if (DEBUG_PENALTY) DEBUG_VAR(vecRefAminoAcids.size());
  // Create AA sequences (first index) in the penalty matrix
  for(int i = 0; i < vecRefAminoAcids.size(); i++) {
    if (DEBUG_PENALTY) DEBUG_VAR(i);
    char aa = vecRefAminoAcids[i].first;

    if (DEBUG_PENALTY) DEBUG_VAR(aa);
    string strAA("X");  // Create a dummy string
    strAA[0] = aa;      // Set the string to the AA char
    std::map<float, float> newMap;
    if (DEBUG_PENALTY) DEBUG_VAR(strAA);
    penalties[strAA] = newMap; // Create a new map for this AA

    float mass = vecRefAminoAcids[i].second;
    if (DEBUG_PENALTY) DEBUG_VAR(mass);
    // Sanity check.. shouldn't happen
    if (mass == 0.0) {
      continue;
    }
    mapCharMods[strAA] = mass; // map of AA to Mass
    mapModChars[mass] = strAA; // map of Mass to AA
  }

  // Fill in masses for I/L and K/Q if not already in map
  if (mapCharMods.find(stringL) == mapCharMods.end() && mapCharMods.find(stringI) != mapCharMods.end()) {
    if (DEBUG_PENALTY) DEBUG_VAR(stringL);
    mapCharMods[stringL] = mapCharMods[stringI];
    if (DEBUG_PENALTY) DEBUG_VAR(stringL);
    if (DEBUG_PENALTY) DEBUG_VAR(mapCharMods[stringL]);
  }
  if (mapCharMods.find(stringI) == mapCharMods.end() && mapCharMods.find(stringL) != mapCharMods.end()) {
    mapCharMods[stringI] = mapCharMods[stringL];
    if (DEBUG_PENALTY) DEBUG_VAR(stringI);
    if (DEBUG_PENALTY) DEBUG_VAR(mapCharMods[stringI]);
  }

  if (mapCharMods.find(stringK) == mapCharMods.end() && mapCharMods.find(stringQ) != mapCharMods.end()) {
    mapCharMods[stringK] = mapCharMods[stringQ];
    if (DEBUG_PENALTY) DEBUG_VAR(stringK);
    if (DEBUG_PENALTY) DEBUG_VAR(mapCharMods[stringK]);
  }
  if (mapCharMods.find(stringQ) == mapCharMods.end() && mapCharMods.find(stringK) != mapCharMods.end()) {
    mapCharMods[stringQ] = mapCharMods[stringK];
    if (DEBUG_PENALTY) DEBUG_VAR(stringQ);
    if (DEBUG_PENALTY) DEBUG_VAR(mapCharMods[stringQ]);
  }

  m_cleavagePenalties.resize(3);
  
  return;
}

// -------------------------------------------------------------------------
PenaltyMatrix::~PenaltyMatrix(void)
{
}

// -------------------------------------------------------------------------
float PenaltyMatrix::operator()(char & aa, float mass, float averageSpectrumPeakIntensity /* = 0.0*/)
{
  string strAA("X");
  strAA[0] = aa;
  return operator()(strAA, mass, averageSpectrumPeakIntensity);
}

// -------------------------------------------------------------------------
float PenaltyMatrix::operator()(string & seq, float mass, float averageSpectrumPeakIntensity)
{
  if (penalties.find(seq) == penalties.end()) {
    if (averageSpectrumPeakIntensity != 0.0) {
      return m_unknownPenalty * averageSpectrumPeakIntensity;
    }
    return m_unknownPenalty * m_allSpectraAveragePeakIntensity;
  }

  float massRounded = roundMass(mass);

  if (penalties[seq].find(massRounded) == penalties[seq].end()) {
    if (averageSpectrumPeakIntensity != 0.0) {
      return m_unknownPenalty * averageSpectrumPeakIntensity;
    }
    return m_unknownPenalty * m_allSpectraAveragePeakIntensity;
  }

  float basePenalty = penalties[seq][massRounded];
  // Known mods are not adjusted for average peak intensity
  if (m_knownMods[seq].find(massRounded) != m_knownMods[seq].end()) {
    return basePenalty;
  }
  if (averageSpectrumPeakIntensity != 0.0) {
    return basePenalty * averageSpectrumPeakIntensity;
  } else {
    return basePenalty * m_allSpectraAveragePeakIntensity;
  }

  return penalties[seq][massRounded] * averageSpectrumPeakIntensity;
}

// -------------------------------------------------------------------------
char PenaltyMatrix::ntermChar(void)
{
  return charNterm;
}
// -------------------------------------------------------------------------
string PenaltyMatrix::ntermString(void)
{
  return stringNterm;
}
// -------------------------------------------------------------------------
char PenaltyMatrix::ctermChar(void)
{
  return charCterm;
}
// -------------------------------------------------------------------------
string PenaltyMatrix::ctermString(void)
{
  return stringCterm;
}

//-----------------------------------------------------------------------------
bool PenaltyMatrix::existsNterm(string & seq)
{
  if (m_knownNtermMods.find(seq) == m_knownNtermMods.end()) {
    return false;
  }
  return m_knownNtermMods[seq].size() != 0;
}  

//-----------------------------------------------------------------------------
bool PenaltyMatrix::existsCterm(string & seq)
{
  if (m_knownCtermMods.find(seq) == m_knownCtermMods.end()) {
    return false;
  }
  return m_knownCtermMods[seq].size() != 0;
}  

//-----------------------------------------------------------------------------
const map<string, set<float> > & PenaltyMatrix::getKnownMods(void)
{
  return m_knownMods;
}

//-----------------------------------------------------------------------------
const set<float> & PenaltyMatrix::getNtermMods(void)
{
  return m_knownNtermMods[stringNterm];
}

//-----------------------------------------------------------------------------
const set<float> & PenaltyMatrix::getNtermMods(string & strAA)
{
  return m_knownNtermMods[strAA];
}

//-----------------------------------------------------------------------------
const map<string, set<float> > & PenaltyMatrix::getAllNtermMods(void)
{
  return m_knownNtermMods;
}

//-----------------------------------------------------------------------------
const set<float> & PenaltyMatrix::getCtermMods(void)
{
  return m_knownCtermMods[stringCterm];
}

//-----------------------------------------------------------------------------
const set<float> & PenaltyMatrix::getCtermMods(string & strAA)
{
  return m_knownCtermMods[strAA];
}

//-----------------------------------------------------------------------------
const map<string, set<float> > & PenaltyMatrix::getAllCtermMods(void)
{
  return m_knownCtermMods;
}

// -------------------------------------------------------------------------
bool PenaltyMatrix::isInMatrix(string & seq, float mass)
{
  if (penalties.find(seq) == penalties.end()) {
    return false;
  }
  float massRounded = roundMass(mass);
  if (penalties[seq].find(massRounded) == penalties[seq].end()) {
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------
bool PenaltyMatrix::isKnown(string & seq, float mass)
{
  if (m_knownMods.find(seq) == m_knownMods.end()) {
    return false;
  }
  float massRounded = roundMass(mass);
  if (m_knownMods[seq].find(massRounded) == m_knownMods[seq].end()) {
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------
bool PenaltyMatrix::isNterm(float mass)
{
  return isNterm(stringNterm, mass);
}

// -------------------------------------------------------------------------
bool PenaltyMatrix::isNterm(string & seq, float mass)
{
  if (m_knownNtermMods.find(seq) == m_knownNtermMods.end()) {
    return false;
  }
  float massRounded = roundMass(mass);
  if (m_knownNtermMods[seq].find(massRounded) == m_knownNtermMods[seq].end()) {
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------
bool PenaltyMatrix::isCterm(float mass)
{
  return isCterm(stringCterm, mass);
}

// -------------------------------------------------------------------------
bool PenaltyMatrix::isCterm(string & seq, float mass)
{
  if (m_knownCtermMods.find(seq) == m_knownCtermMods.end()) {
    return false;
  }
  float massRounded = roundMass(mass);
  if (m_knownCtermMods[seq].find(massRounded) == m_knownCtermMods[seq].end()) {
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------
void PenaltyMatrix::getPenalties(string & seq, std::map<float, float> & penaltyMap, float averageSpectrumPeakIntensity)
{
  map<string, map<float, float> >::iterator itrMap = penalties.find(seq);

  if (itrMap == penalties.end()) {
    return;
  }

  map<float, float>::iterator itr = itrMap->second.begin();
  map<float, float>::iterator itr_end = itrMap->second.end();
  for (; itr != itr_end; itr++) {

    float basePenalty = itr->second;
    //float basePenalty = -1.0;
    // Known mods are not adjusted for average peak intensity
    if (m_knownMods[itrMap->first].find(itr->first) != m_knownMods[itrMap->first].end()) {
      penaltyMap[itr->first] = basePenalty;
      //float basePenalty = -0.5;
    }
    if (averageSpectrumPeakIntensity != 0.0) {
      penaltyMap[itr->first] = basePenalty * averageSpectrumPeakIntensity;
    } else {
      penaltyMap[itr->first] = basePenalty * m_allSpectraAveragePeakIntensity;
    }
  }

  return;
}

// -------------------------------------------------------------------------
float PenaltyMatrix::getUnknownPenalty(float averageSpectrumPeakIntensity)
{
  if (averageSpectrumPeakIntensity != 0.0) {
    return m_unknownPenalty * averageSpectrumPeakIntensity;
  }
  return m_unknownPenalty * m_allSpectraAveragePeakIntensity;
}

// -------------------------------------------------------------------------
float PenaltyMatrix::getKnownPenalty(float averageSpectrumPeakIntensity)
{
  if (averageSpectrumPeakIntensity != 0.0) {
    return m_knownModPenalty * averageSpectrumPeakIntensity;
  }
  return m_knownModPenalty;
}

// -------------------------------------------------------------------------
void PenaltyMatrix::getAminoAcids(vector<string> & aaVec)
{
  std::map<std::string, float>::iterator itr = mapCharMods.begin();
  std::map<std::string, float>::iterator itr_end = mapCharMods.end();
  for ( ; itr != itr_end; itr++) {
    aaVec.push_back(itr->first);
  }
}

// -------------------------------------------------------------------------
float PenaltyMatrix::getCleavagePenalty(char c, int location)
{
  if (DEBUG_PENALTY) { string XX(" "); XX[0] = c; DEBUG_VAR(XX) }
  if (DEBUG_PENALTY) DEBUG_VAR(location);
  if (m_cleavagePenalties[location].find(c) == m_cleavagePenalties[location].end()) {
    if (DEBUG_PENALTY) DEBUG_TRACE;
    // In the case of start or end cleavages.. absence calls for penalty
    if (location == CLEAVAGE_START) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      if (m_cleavagePenalties[CLEAVAGE_START].size() == 0) {
        return 0.0;
      } else {
        return m_cleavagePenalties[CLEAVAGE_START].begin()->second;
      }
    } else if (location == CLEAVAGE_END) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      if (m_cleavagePenalties[CLEAVAGE_END].size() == 0) {
        return 0.0;
      } else {
        return m_cleavagePenalties[CLEAVAGE_END].begin()->second;
      }
    }
  } else if (location == CLEAVAGE_INTERNAL) {
    if (DEBUG_PENALTY) DEBUG_TRACE;
    // In the case of internal cleavage.. presence calls for penalty
    return m_cleavagePenalties[CLEAVAGE_INTERNAL][c];
  }

  if (DEBUG_PENALTY) DEBUG_TRACE;
  return 0.0;
}

//-----------------------------------------------------------------------------
// Load and process a BLOSUM matrix file
//-----------------------------------------------------------------------------
bool PenaltyMatrix::loadFromBlosum(std::string & filename, float peakEquivalents)
{
  vector<string> aa;

  // Read the tab delimited file
  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(filename.c_str(),
                                                      " \t",
                                                      "",
                                                      lines)) {
    ERROR_MSG("Unable to open file! " << filename);
    return false;
  }

  float maxScore = 0.0;

  int firstLineIdx = -1;

  if (DEBUG_PENALTY) DEBUG_VAR(lines.size());
  for (int i = 0; i < lines.size(); i++) {
    // Discard empty lines
    if (DEBUG_PENALTY) DEBUG_VAR(lines[i].size());
    if (lines[i].size() == 0) {
      continue;
    }
    // Discard comments lines
    if (lines[i][0].compare("#") == 0) {
      continue;
    }

    if (firstLineIdx == -1) {
      firstLineIdx = i;
    }
    if (DEBUG_PENALTY) DEBUG_VAR(firstLineIdx);

    int aaIndex1 = i - firstLineIdx - 1;
    int aaIndex2 = 0;
    if (DEBUG_PENALTY) DEBUG_VAR(aaIndex1);
    if (DEBUG_PENALTY) DEBUG_VAR(aaIndex2);
    for (int j = 0; j < lines[i].size(); j++) {
      // Discard extra spaces and empty strings
      if (lines[i][j].length() == 0 || lines[i][j].compare(" ") == 0) {
        continue;
      }

      if (DEBUG_PENALTY) DEBUG_VAR(lines[i][j]);

      if (aaIndex1 == -1) {
        aa.push_back(lines[i][j]);
      } else {
        if (j == 0) {
          // Do nothing - the first char is the AA
        } else {
          string aa1 = aa[aaIndex1];
          int score;
          sscanf(lines[i][j].c_str(), "%d", &score);
          if (DEBUG_PENALTY) DEBUG_VAR(score);
          float mass1 = getMass(aa[aaIndex1]);
          float mass2 = getMass(aa[aaIndex2]);
          float massDiff = mass2 - mass1;
          // Round to the desired precision
          massDiff = roundMass(massDiff);

          if (massDiff != 0.0) {
            penalties[aa[aaIndex1]][massDiff]  = score;
            // Save the max score so we can subtract it out later so all scores are negative
            if (maxScore < score) {
              maxScore = score;
            }
          }

          aaIndex2++;
        }
      }

    } //for (int j = 0; j < lines[i].size(); j++) {

  } // for (int i = 0; i < lines.size(); i++) {

  if (DEBUG_PENALTY) DEBUG_VAR(maxScore);

  float minPenalty = -1.0;
  map<string, map<float, float> >::iterator itr = penalties.begin();
  map<string, map<float, float> >::iterator itr_end = penalties.end();
  for (; itr != itr_end; itr++) {
    map<float, float>::iterator itr2 = itr->second.begin();
    map<float, float>::iterator itr_end2 = itr->second.end();
    for (; itr2 != itr_end2; itr2++) {
      itr2->second -= maxScore;

      // Make the diagonals (mass diff 0) penalty be 0
      if (itr2->first == 0.0) {
        itr2->second = 0;
      }

      // LARS - DEBUG - EMPTY PENALTIES
      //itr2->second = 0;

      if (itr2->second < minPenalty) {
        minPenalty = itr2->second;
      }
    }
  }

  // Sanity check
  if (minPenalty == 0.0) minPenalty = -1.0;

  // Normalize the penalties to 1.0
  for (map<string, map<float, float> >::iterator itr = penalties.begin(); itr != penalties.end(); itr++) {
    map<float, float>::iterator itr2 = itr->second.begin();
    map<float, float>::iterator itr_end2 = itr->second.end();
    for (; itr2 != itr_end2; itr2++) {
      if (DEBUG_PENALTY) DEBUG_MSG("[" << itr->first << "][" << itr2->first << "] = " << itr2->second);
      itr2->second = itr2->second * peakEquivalents / fabs(minPenalty);
      if (DEBUG_PENALTY) DEBUG_MSG("[" << itr->first << "][" << itr2->first << "] = " << itr2->second);
    }
    if (DEBUG_PENALTY) DEBUG_TRACE;
  }

  m_unknownPenalty = -m_unknownMultiplier * peakEquivalents;
  DEBUG_VAR(m_unknownPenalty);

  return true;
}

//-----------------------------------------------------------------------------
// Create the penalties from modification count data
//-----------------------------------------------------------------------------
bool PenaltyMatrix::createFromModificationFreqs(map<float, float> & modFreq,
                                                float minPeakEquivalents,
                                                float maxPeakEquivalents,
                                                float minFrequency,
                                                float averagePeakIntensity)
{
  m_allSpectraAveragePeakIntensity = averagePeakIntensity;
  if (DEBUG_PENALTY) DEBUG_VAR(m_allSpectraAveragePeakIntensity);

  float maxPenalty = 0.0;
  float minPenalty = 1000.0;

  // Create a set of all the known mass values
  set<float> knownMasses;
  map<string, set<float> >::iterator itrk = m_knownMods.begin();
  map<string, set<float> >::iterator itrk_end = m_knownMods.end();
  for(; itrk != itrk_end; itrk++) {
    set<float>::iterator itrs = itrk->second.begin();
    set<float>::iterator itrs_end = itrk->second.end();
    for(; itrs != itrs_end; itrs++) {
      knownMasses.insert(*itrs);
      knownMasses.insert(-(*itrs));
    }
  }
  itrk = m_knownNtermMods.begin();
  itrk_end = m_knownNtermMods.end();
  for(; itrk != itrk_end; itrk++) {
    set<float>::iterator itrs = itrk->second.begin();
    set<float>::iterator itrs_end = itrk->second.end();
    for(; itrs != itrs_end; itrs++) {
      knownMasses.insert(*itrs);
      knownMasses.insert(-(*itrs));
    }
  }
  itrk = m_knownCtermMods.begin();
  itrk_end = m_knownCtermMods.end();
  for(; itrk != itrk_end; itrk++) {
    set<float>::iterator itrs = itrk->second.begin();
    set<float>::iterator itrs_end = itrk->second.end();
    for(; itrs != itrs_end; itrs++) {
      knownMasses.insert(*itrs);
      knownMasses.insert(-(*itrs));
    }
  }

  // Create a set of all the known mass values
  set<int> setRoundedAAMasses;
  map<string, float>::iterator itrA = mapCharMods.begin();
  map<string, float>::iterator itrA_end = mapCharMods.end();
  for(; itrA != itrA_end; itrA++) {
    setRoundedAAMasses.insert(roundMass(itrA->second));
  }

  // Transmute modification frequencies to penalties
  map<float, float>::iterator itr = modFreq.begin();
  map<float, float>::iterator itr_end = modFreq.end();
  for(; itr != itr_end; itr++) {
    float mass = itr->first;
    float freq = itr->second;

    if (DEBUG_PENALTY) DEBUG_MSG("Mass [" << mass << "] has frequency [" << freq << "]");
    if (mass < m_minModMass || mass > m_maxModMass || freq < minFrequency) {
      if (DEBUG_PENALTY) DEBUG_MSG("Undesirable mod");
      continue;
    }

    // Is this an amino acid mass? If so ignore it.. its not a modification (its an insertion)
    float roundedMass = roundMass(mass);
    if (setRoundedAAMasses.find(roundedMass) != setRoundedAAMasses.end()) {
      if (DEBUG_PENALTY) DEBUG_MSG("Mass is an amino acid insertion [" << mass << "] at resolution [" << m_resolution << "]");
      continue;
    }

    map<string, map<float, float> > ::iterator itrP = penalties.begin();
    map<string, map<float, float> > ::iterator itrP_end = penalties.end();
    for(; itrP != itrP_end; itrP++) {

      if (DEBUG_PENALTY) DEBUG_VAR(itrP->first);
      float aaMass = 0.0;
      if (mapCharMods.find(itrP->first) != mapCharMods.end()) {
        aaMass = mapCharMods[itrP->first];
      }

#if 1
      // We will not consider 1,2 3 dalton mods
      if (roundedMass < 4.0) {
        continue;
      }
#endif

      if (DEBUG_PENALTY) DEBUG_VAR(aaMass);
      float totalRoundedMass = roundMass(aaMass + mass);
      if (DEBUG_PENALTY) DEBUG_VAR(totalRoundedMass);

#if 0
      // Is this an amino acid mass? If so ignore it.. its not a modification (its a mutation)
      if (setRoundedAAMasses.find(totalRoundedMass) != setRoundedAAMasses.end()) {
        if (DEBUG_PENALTY) DEBUG_MSG("Mass is an amino acid mutation [" << totalRoundedMass << "] at resolution [" << m_resolution << "]");
        continue;
      }
#endif

      float modf = modFreq[roundedMass];
      if (DEBUG_PENALTY) DEBUG_VAR(modf);
      if (modf == 0.0) {
        continue;
      }

      // See if this is a known nterm modification
      bool ntermFound = false;
      if (m_knownNtermMods[itrP->first].find(roundedMass) != m_knownNtermMods[itrP->first].end()) {
        ntermFound = true;
      }
      // See if this is a known cterm modification
      bool ctermFound = false;
      if (m_knownCtermMods[itrP->first].find(roundedMass) != m_knownCtermMods[itrP->first].end()) {
        ctermFound = true;
      }
      if (ntermFound || ctermFound) {
        continue;
      }

      // Leaving other amino acids with known mods leads to over-representation of known mods on
      //   unintended amino acids, so we will not allow a known mass (or the negative thereof) on ANY amino acid
      if (m_knownMods[itrP->first].find(roundedMass) != m_knownMods[itrP->first].end()) {
        //penalties[itrP->first][roundedMass] = m_knownModPenalty / -log10(modf);
        //if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has known mod [" << roundedMass
        //                 << "], penalty = " << penalties[itrP->first][roundedMass]);
        continue;
      } else if (knownMasses.find(roundedMass) != knownMasses.end()) {
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has known mass [" << roundedMass << "]");
        // Known Mass will get unknown mod (if not on known AA)
        penalties[itrP->first][roundedMass] = 0.0;
      } else {
        penalties[itrP->first][roundedMass] = modf;
        //penalties[itrP->first][roundedMass] = log10(modf);
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has unknown mass [" << roundedMass
                         << "], penalty = " << penalties[itrP->first][roundedMass]);
      }

      // This amino acid, negative mass combination is too small to be real
      if (aaMass - roundedMass < 57.0) {
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] with mod [" << -roundedMass
                                            << "], not allowed (total mass too small)");
        continue;
      } else if (m_knownMods[itrP->first].find(-roundedMass) != m_knownMods[itrP->first].end()) {
        // This amino acid, negative mass combination is known
        //penalties[itrP->first][-roundedMass] = m_knownModPenalty / -log10(modf);
        //if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has known mod [" << -roundedMass
        //                 << "], penalty = " << penalties[itrP->first][roundedMass]);
        continue;
      } else if (knownMasses.find(-roundedMass) != knownMasses.end()) {
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has known mass [" << -roundedMass << "]");
        // Known Mass will get unknown mod (if not on known AA)
        penalties[itrP->first][-roundedMass] = 0.0;
      } else {
        //penalties[itrP->first][-roundedMass] = log10(modf);
        penalties[itrP->first][-roundedMass] = modf;
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has unknown mod [" << -roundedMass
                         << "], penalty = " << penalties[itrP->first][-roundedMass]);
      }

      if (penalties[itrP->first][roundedMass] > maxPenalty) {
        maxPenalty = penalties[itrP->first][roundedMass];
      }
      if (penalties[itrP->first][-roundedMass] > maxPenalty) {
        maxPenalty = penalties[itrP->first][-roundedMass];
      }

      if (penalties[itrP->first][roundedMass] != 0.0 &&
          penalties[itrP->first][roundedMass] < minPenalty) {
        minPenalty = penalties[itrP->first][roundedMass];
      }
      if (penalties[itrP->first][-roundedMass] != 0.0 &&
          penalties[itrP->first][-roundedMass] < minPenalty) {
        minPenalty = penalties[itrP->first][-roundedMass];
      }

    } // penalties iterator

  } // modFreq iterator

  // Sanity check
  if (maxPenalty == 0.0) maxPenalty = 1.0;
  if (minPenalty == 0.0) minPenalty = 1.0;
  if (DEBUG_PENALTY) DEBUG_VAR(maxPenalty);
  if (DEBUG_PENALTY) DEBUG_VAR(minPenalty);

  float freqSpread = maxPenalty - minPenalty;
  if (freqSpread <= 0.0) freqSpread =1.0;
  if (DEBUG_PENALTY) DEBUG_VAR(freqSpread);

  if (DEBUG_PENALTY) DEBUG_VAR(minPeakEquivalents);
  if (DEBUG_PENALTY) DEBUG_VAR(maxPeakEquivalents);
  float penaltyRange = maxPeakEquivalents - minPeakEquivalents;

  if (DEBUG_PENALTY) DEBUG_VAR(penaltyRange);
  if (DEBUG_PENALTY) DEBUG_VAR(m_knownModPenalty);

  // Normalize the penalties to between minPeakEquivalents and maxPeakEquivalents
  for (map<string, map<float, float> >::iterator itr = penalties.begin(); itr != penalties.end(); itr++) {
    map<float, float>::iterator itr2 = itr->second.begin();
    map<float, float>::iterator itr_end2 = itr->second.end();
    for (; itr2 != itr_end2; itr2++) {
      if (DEBUG_PENALTY) DEBUG_MSG("[" << itr->first << "][" << itr2->first << "] = " << itr2->second << " freq");

      // These were AA's with known masses but not on the known Mod AA
      //    (like K,16 instead of M,16)
      if (itr2->second == 0.0) {
        itr2->second = maxPeakEquivalents;
      } else {
        itr2->second = maxPeakEquivalents - (itr2->second - minPenalty) * penaltyRange / freqSpread;
      }

      if (DEBUG_PENALTY) DEBUG_MSG("[" << itr->first << "][" << itr2->first << "] = " << itr2->second);

      if (itr2->second < minPeakEquivalents) {
        itr2->second = minPeakEquivalents;
      }

      // Penalties are negative
      itr2->second = -itr2->second;

      if (DEBUG_PENALTY) DEBUG_MSG("[" << itr->first << "][" << itr2->first << "] = " << itr2->second);
    }
    if (DEBUG_PENALTY) DEBUG_TRACE;
  }

    if (DEBUG_PENALTY) DEBUG_TRACE;
  // Add in all the known modifications at fixed penalty
  for(itrk = m_knownMods.begin(); itrk != m_knownMods.end(); itrk++) {
    set<float>::iterator itrs = itrk->second.begin();
    set<float>::iterator itrs_end = itrk->second.end();
    for(; itrs != itrs_end; itrs++) {
      penalties[itrk->first][roundMass(*itrs)] = m_knownModPenalty;
    }
  }

  m_unknownPenalty = -m_unknownMultiplier * maxPeakEquivalents;
  if (DEBUG_PENALTY) DEBUG_VAR(m_unknownPenalty);

  return true;
}


//-----------------------------------------------------------------------------
float PenaltyMatrix::getMass(string aa)
{
  if (mapCharMods.size() == 0) {
    WARN_MSG("Internal modification map is empty in getMass() method");
    return 0.0;
  }
  if (mapCharMods.find(aa) == mapCharMods.end()) {
    WARN_MSG("No mass for [" << aa << "]");
  }
  if (mapCharMods[aa] == 0.0) {
    WARN_MSG("Mass of [" << aa << "] is 0.0");
  }

  return mapCharMods[aa];
}

//-----------------------------------------------------------------------------
float PenaltyMatrix::roundMass(float mass)
{
  int sign = mass < 0 ? -1 : 1;
  return (int)(fabs(mass) * 0.9995 + 0.5) * sign;
  //return abs((float)((int)(mass / m_resolution)) * m_resolution);
}

//-----------------------------------------------------------------------------
// Load the cleavage penalties
//-----------------------------------------------------------------------------
bool PenaltyMatrix::loadCleavagePenalties(std::string & filename)
{
  m_cleavagePenalties.resize(3);

  if (DEBUG_PENALTY) DEBUG_TRACE;
  // Read the tab delimited file
  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(filename.c_str(),
                                                      " \t",
                                                      "",
                                                      lines)) {
    ERROR_MSG("Unable to open file! " << filename);
    return false;
  }

  string aa;
  if (DEBUG_PENALTY) DEBUG_VAR(lines.size());
  for (int i = 0; i < lines.size(); i++) {

    // Discard empty lines
    if (lines[i].size() == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }
    // Discard comments lines
    if (lines[i][0].compare("#") == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }

    float penalty = 0.0;
    if (DEBUG_PENALTY) DEBUG_VAR(lines[i].size());
    if (lines[i].size() != 2) {
      WARN_MSG("Malformed cleavage penalty at line " << i+1 << " in file [" << filename << "]");
      continue;
    }

    // First item is the amino acid (plus location indicator)
    string aa = lines[i][0];
    if (DEBUG_PENALTY) DEBUG_VAR(aa);

    // Second item is the penalty
    sscanf(lines[i][1].c_str(), "%f", &penalty);
    if (DEBUG_PENALTY) DEBUG_VAR(penalty);

    if (aa[0] == charNterm) {

      m_cleavagePenalties[CLEAVAGE_START][aa[1]] = penalty;
      if (DEBUG_PENALTY) DEBUG_MSG("Start penalty for amino acid [" << aa[1] << "] is [" << penalty << "]");

    } else if (aa[0] == charCterm) {

      m_cleavagePenalties[CLEAVAGE_END][aa[1]] = penalty;
      if (DEBUG_PENALTY) DEBUG_MSG("End penalty for amino acid [" << aa[1] << "] is [" << penalty << "]");

    } else if (aa[0] == charInternal) {

      m_cleavagePenalties[CLEAVAGE_INTERNAL][aa[1]] = penalty;
      if (DEBUG_PENALTY) DEBUG_MSG("Internal penalty for amino acid [" << aa[1] << "] is [" << penalty << "]");

    } else {
      WARN_MSG("Unknown cleavage location specifier [" << aa[0] << "]");
    }

  } // for (int i = 0; i < lines.size(); i++) {

  return true;
}


//-----------------------------------------------------------------------------
// Load the known modifications
//-----------------------------------------------------------------------------
bool PenaltyMatrix::loadKnownModifications(std::string & filename)
{
  if (DEBUG_PENALTY) DEBUG_TRACE;
  // Read the tab delimited file
  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(filename.c_str(),
                                                      " \t",
                                                      "",
                                                      lines)) {
    ERROR_MSG("Unable to open file! " << filename);
    return false;
  }

  string aa;
  if (DEBUG_PENALTY) DEBUG_VAR(lines.size());
  for (int i = 0; i < lines.size(); i++) {

    // Discard empty lines
    if (lines[i].size() == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }
    // Discard comments lines
    if (lines[i][0].compare("#") == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }

    float mass = 0.0;
    float roundedMass = 0.0;
    if (DEBUG_PENALTY) DEBUG_VAR(lines[i].size());
    for (int j = 0; j < lines[i].size(); j++) {

      if (j == 0) {
        // First item is the amino acid
        aa = lines[i][j];
        if (DEBUG_PENALTY) DEBUG_VAR(aa);
      } else {
        // Rest of the items are masses of the modifications
        sscanf(lines[i][j].c_str(), "%f", &mass);
        if (DEBUG_PENALTY) DEBUG_VAR(mass);
        roundedMass = roundMass(mass);
        if (DEBUG_PENALTY) DEBUG_VAR(roundedMass);

        if (aa[0] == charNterm) {
	  string strAA = stringNterm;
	  if (aa.length() == 2) {
	    strAA[0] = aa[1];
	  }
          m_knownNtermMods[strAA].insert(roundedMass);
          if (DEBUG_PENALTY) DEBUG_MSG("Known nterm modification of amino acid [" << strAA << "] [" << roundedMass << "]");
	} else if (aa[0] == charCterm) {
	  string strAA = stringCterm;
	  if (aa.length() == 2) {
	    strAA[0] = aa[1];
	  }
          m_knownCtermMods[strAA].insert(roundedMass);
          if (DEBUG_PENALTY) DEBUG_MSG("Known cterm modification of amino acid [" << strAA << "] [" << roundedMass << "]");
        } else {
          m_knownMods[aa].insert(roundedMass);
          if (DEBUG_PENALTY) DEBUG_MSG("Known modification of amino acid [" << aa << "] [" << roundedMass << "]");
        }
      }

    } //for (int j = 0; j < lines[i].size(); j++) {

  } // for (int i = 0; i < lines.size(); i++) {

  return true;
}

//-----------------------------------------------------------------------------
// Load a generic matrix
//-----------------------------------------------------------------------------
bool PenaltyMatrix::load(std::string & filename,
                         std::string & knowmModsFileName,
                         std::string & cleavagePenaltiesFileName)
{
  if (DEBUG_PENALTY) DEBUG_TRACE;
  // Read the tab delimited file
  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(filename.c_str(),
                                                      " \t",
                                                      "",
                                                      lines)) {
    ERROR_MSG("Unable to open file! " << filename);
    return false;
  }

  float minPenalty = 0.0;

  if (DEBUG_PENALTY) DEBUG_VAR(lines.size());
  for (int i = 0; i < lines.size(); i++) {

    // Discard empty lines
    if (lines[i].size() == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }
    // Discard comments lines
    if (lines[i][0].compare("#") == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }

    string aa;
    float mass = 0.0;
    float roundedMass = 0.0;
    float penalty = 0.0;
    if (DEBUG_PENALTY) DEBUG_VAR(lines[i].size());
    for (int j = 0; j < lines[i].size(); j++) {

      if (j == 0) {
        // First item is the amino acid
        aa = lines[i][j];
      } else if (j == 1) {
        // Second item is the mass of the modification
        sscanf(lines[i][j].c_str(), "%f", &mass);
        if (DEBUG_PENALTY) DEBUG_VAR(mass);
        roundedMass = roundMass(mass);
        if (DEBUG_PENALTY) DEBUG_VAR(roundedMass);
      } else if (j == 2) {
        // Third item is the penalty
        sscanf(lines[i][j].c_str(), "%f", &penalty);
        if (DEBUG_PENALTY) DEBUG_VAR(penalty);

        if (aa[0] == charNterm) {
          m_knownNtermMods[aa].insert(roundedMass);
	      } else if (aa[0] == charCterm) {
          m_knownCtermMods[aa].insert(roundedMass);
        } else {
          penalties[aa][roundedMass] = penalty;
        }

        if (penalty < minPenalty) {
          minPenalty = penalty;
        }
      }

    } //for (int j = 0; j < lines[i].size(); j++) {

  } // for (int i = 0; i < lines.size(); i++) {

  m_unknownPenalty = minPenalty  * m_unknownMultiplier;
  DEBUG_VAR(m_unknownPenalty);

  if (!knowmModsFileName.empty()) {
    if (!loadKnownModifications(knowmModsFileName)) {
      return false;
    }
  }
  if (!cleavagePenaltiesFileName.empty()) {
    if (!loadCleavagePenalties(cleavagePenaltiesFileName)) {
      return false;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Save the penalty matrix
//-----------------------------------------------------------------------------
bool PenaltyMatrix::saveMatrix(std::string & filename)
{
  ofstream ofs(filename.c_str(), ios_base::out | ios_base::binary);
  if (!ofs) {
    return false;
  }

  std::map<std::string, std::map<float, float> >::iterator itr = penalties.begin();
  std::map<std::string, std::map<float, float> >::iterator itrEnd = penalties.end();
  for ( ; itr != itrEnd; itr++) {
    std::map<float, float>::iterator itrMap = itr->second.begin();
    std::map<float, float>::iterator itrMapEnd = itr->second.end();
    for ( ; itrMap != itrMapEnd; itrMap++) {
      ofs << itr->first << "\t" << itrMap->first <<  "\t" << itrMap->second << endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Save the known mods
//-----------------------------------------------------------------------------
bool PenaltyMatrix::saveKnownMods(std::string & filename)
{
  if (DEBUG_PENALTY) DEBUG_TRACE;
  ofstream ofs(filename.c_str(), ios_base::out | ios_base::binary);
  if (!ofs) {
    return false;
  }

  if (DEBUG_PENALTY) DEBUG_TRACE;
  map<string, set<float> >::iterator itr = m_knownMods.begin();
  map<string, set<float> >::iterator itrEnd = m_knownMods.end();
  for ( ; itr != itrEnd; itr++) {
    if (itr->second.size() == 0) {
      continue;
    }
    ofs << itr->first;
    if (DEBUG_PENALTY) DEBUG_VAR(itr->first);
    set<float>::iterator itrSet = itr->second.begin();
    set<float>::iterator itrSetEnd = itr->second.end();
    for ( ; itrSet != itrSetEnd; itrSet++) {
      ofs << "\t" << *itrSet;
      if (DEBUG_PENALTY) DEBUG_VAR(*itrSet);
    }
    ofs << endl;
  }

  if (DEBUG_PENALTY) DEBUG_TRACE;
  itr = m_knownNtermMods.begin();
  itrEnd = m_knownNtermMods.end();
  for ( ; itr != itrEnd; itr++) {
    if (itr->second.size() == 0) {
      continue;
    }
    if (itr->first != stringNterm) {
      ofs << stringNterm;
    }
    ofs << itr->first;
    if (DEBUG_PENALTY) DEBUG_VAR(itr->first);
    set<float>::iterator itrSet = itr->second.begin();
    set<float>::iterator itrSetEnd = itr->second.end();
    for ( ; itrSet != itrSetEnd; itrSet++) {
      ofs << "\t" << *itrSet;
      if (DEBUG_PENALTY) DEBUG_VAR(*itrSet);
    }
    ofs << endl;
  }

  if (DEBUG_PENALTY) DEBUG_TRACE;
  itr = m_knownCtermMods.begin();
  itrEnd = m_knownCtermMods.end();
  for ( ; itr != itrEnd; itr++) {
    if (itr->second.size() == 0) {
      continue;
    }
    if (itr->first != stringCterm) {
      ofs << stringCterm;
    }
    ofs << itr->first;
    if (DEBUG_PENALTY) DEBUG_VAR(itr->first);
    set<float>::iterator itrSet = itr->second.begin();
    set<float>::iterator itrSetEnd = itr->second.end();
    for ( ; itrSet != itrSetEnd; itrSet++) {
      ofs << "\t" << *itrSet;
      if (DEBUG_PENALTY) DEBUG_VAR(*itrSet);
    }
    ofs << endl;
  }

  if (DEBUG_PENALTY) DEBUG_TRACE;
  ofs.close();

  return true;
}

//-----------------------------------------------------------------------------
// Save the cleavage penalties
//-----------------------------------------------------------------------------
bool PenaltyMatrix::saveCleavagePenalties(std::string & filename)
{
  ofstream ofs(filename.c_str(), ios_base::out | ios_base::binary);
  if (!ofs) {
    return false;
  }

  map<char, float>::iterator itr = m_cleavagePenalties[CLEAVAGE_START].begin();
  map<char, float>::iterator itrEnd = m_cleavagePenalties[CLEAVAGE_START].end();
  for ( ; itr != itrEnd; itr++) {
    ofs << charNterm << itr->first << "\t" << itr->second << endl;
  }

  itr = m_cleavagePenalties[CLEAVAGE_INTERNAL].begin();
  itrEnd = m_cleavagePenalties[CLEAVAGE_INTERNAL].end();
  for ( ; itr != itrEnd; itr++) {
    ofs << charInternal << itr->first << "\t" << itr->second << endl;
  }

  itr = m_cleavagePenalties[CLEAVAGE_END].begin();
  itrEnd = m_cleavagePenalties[CLEAVAGE_END].end();
  for ( ; itr != itrEnd; itr++) {
    ofs << charCterm << itr->first << "\t" << itr->second << endl;
  }

  ofs.close();

  return true;
}


//-----------------------------------------------------------------------------
// Save the amino acid masses
//-----------------------------------------------------------------------------
bool PenaltyMatrix::saveAminoAcids(std::string & filename)
{
  ofstream ofs(filename.c_str(), ios_base::out | ios_base::binary);
  if (!ofs) {
    return false;
  }

  ofs << mapCharMods.size() << endl;
  std::map<std::string, float>::iterator itr =  mapCharMods.begin();
  std::map<std::string, float>::iterator itrEnd = mapCharMods.end();
  for ( ; itr != itrEnd; itr++) {
    ofs << itr ->first << "=" << itr->second << endl;
  }
  ofs.close();

  return true;
}

//-----------------------------------------------------------------------------
void PenaltyMatrix::getSpecProbMods(vector<pair<unsigned int, bool> > & ntermMods,
				    vector<unsigned int> & mods)
{
  map<string, set<float> >::const_iterator itr = m_knownMods.begin();
  map<string, set<float> >::const_iterator itrEnd = m_knownMods.end();
  for (; itr != itrEnd; itr++) {
    if (itr->second.size() == 0) {
      continue;
    }
    string stringAA = itr->first;
    if (DEBUG_SPECPROB_MODS) DEBUG_VAR(stringAA);
    float aaMass = getMass(itr->first);
    if (DEBUG_SPECPROB_MODS) DEBUG_VAR(aaMass);
    if (aaMass == 0)
      continue; // Sanity check
    const set<float> & setMods = itr->second;
    set<float>::const_iterator itrSet = setMods.begin();
    set<float>::const_iterator itrSetEnd = setMods.end();
    for (; itrSet != itrSetEnd; itrSet++) {
      float mod = *itrSet;
      if (DEBUG_SPECPROB_MODS) DEBUG_VAR(mod);
      float totalMass = aaMass + mod;
      if (DEBUG_SPECPROB_MODS) DEBUG_VAR(totalMass);
      mods.push_back(totalMass);
    }
  }
  if (DEBUG_SPECPROB_MODS) DEBUG_VAR(mods.size());

  itr = m_knownNtermMods.begin();
  itrEnd = m_knownNtermMods.end();
  for ( ; itr != itrEnd; itr++) {
    if (itr->second.size() == 0) {
      continue;
    }

    string stringAA = itr->first;
    if (DEBUG_SPECPROB_MODS) DEBUG_VAR(stringAA);
    float aaMass = 0.0;;
    if (itr->first != stringNterm) {
      aaMass = getMass(itr->first);
      if (aaMass == 0)
        continue; // Sanity check
    }
    if (DEBUG_SPECPROB_MODS) DEBUG_VAR(aaMass);

    set<float>::const_iterator itrSet = itr->second.begin();
    set<float>::const_iterator itrSetEnd = itr->second.end();
    for (; itrSet != itrSetEnd; itrSet++) {
      pair<unsigned int, bool> modPair;
      modPair.first = aaMass + *itrSet;
      modPair.second = true;
      if (itr->first == stringNterm) {
        // We won't tell spec prob about these negative parent mass shift mods
        if (aaMass == -1.0) {
         continue;
        }
        modPair.second = false;
      }
      if (DEBUG_SPECPROB_MODS) DEBUG_VAR(modPair.first);
      if (DEBUG_SPECPROB_MODS) DEBUG_VAR(modPair.second);
      ntermMods.push_back(modPair);
    }
  }
  return;
}

    void addKnownModification(string strAA, float mass);
    void addKnownNtermModification(string strAA, float mass);
    void addKnownCtermModification(string strAA, float mass);

    // This method allows manual addition of cleavage penalties
    void addCleavagePenalty(char location, float mass);

  
//-----------------------------------------------------------------------------
void PenaltyMatrix::addKnownModification(string strAA, float mass)
{
  for (int i = 0; i < strAA.length(); i++) {
    string aaChar("X");
    aaChar[0] = strAA[i];
    m_knownMods[aaChar].insert(mass);
  }
  return;
}

//-----------------------------------------------------------------------------
void PenaltyMatrix::addKnownNtermModification(string strAA, float mass)
{
  for (int i = 0; i < strAA.length(); i++) {
    string aaChar("X");
    aaChar[0] = strAA[i];
    m_knownNtermMods[aaChar].insert(mass);
  }
  return;
}

//-----------------------------------------------------------------------------
void PenaltyMatrix::addKnownCtermModification(string strAA, float mass)
{
  for (int i = 0; i < strAA.length(); i++) {
    string aaChar("X");
    aaChar[0] = strAA[i];
    m_knownCtermMods[aaChar].insert(mass);
  }
  return;
}

//-----------------------------------------------------------------------------
void PenaltyMatrix::addCleavagePenalty(int location, char c, float penalty)
{
  if (location != CLEAVAGE_START &&
      location != CLEAVAGE_INTERNAL &&
      location != CLEAVAGE_END) {
    ERROR_MSG("Unknown location for cleavage penalty [" << location << "]");
    return;
  }
  m_cleavagePenalties[location][c] = penalty;
  return;
}


//-----------------------------------------------------------------------------
void PenaltyMatrix::debug(void)
{
  DEBUG_MSG("Known:");
  map<string, set<float> >::iterator itr = m_knownMods.begin();
  map<string, set<float> >::iterator itrEnd = m_knownMods.end();
  for ( ; itr != itrEnd; itr++) {
    DEBUG_MSG("  " << itr->first);
    set<float>::iterator itr2 = itr->second.begin();
    set<float>::iterator itrEnd2 = itr->second.end();
    for ( ; itr2 != itrEnd2; itr2++) {
      DEBUG_MSG("    " << *itr2);
    }
  }

  DEBUG_MSG("Known Nterm:");
  itr = m_knownNtermMods.begin();
  itrEnd = m_knownNtermMods.end();
  for ( ; itr != itrEnd; itr++) {
    DEBUG_MSG("  " << itr->first);
    set<float>::iterator itr2 = itr->second.begin();
    set<float>::iterator itrEnd2 = itr->second.end();
    for ( ; itr2 != itrEnd2; itr2++) {
      DEBUG_MSG("    " << *itr2);
    }
  }

  DEBUG_MSG("Known Cterm:");
  itr = m_knownCtermMods.begin();
  itrEnd = m_knownCtermMods.end();
  for ( ; itr != itrEnd; itr++) {
    DEBUG_MSG("  " << itr->first);
    set<float>::iterator itr2 = itr->second.begin();
    set<float>::iterator itrEnd2 = itr->second.end();
    for ( ; itr2 != itrEnd2; itr2++) {
      DEBUG_MSG("    " << *itr2);
    }
  }

  return;
}
