/*
 * GFTable.cpp
 *
 *  Created on: Jul 19, 2013
 *      Author: aguthals
 */

#include "GFTable.h"
#include "aminoacid.h"
#include "alignment_scoring.h"
#include "Logger.h"

using namespace std;

namespace specnets
{

  bool GFTable::DEBUG = false;

  const double GFTable::PROBABILITY_FACTOR = 1e25;

  const float GFTable::EDGE_SCORE_SCALING_FACTOR = 0.50;

  //vector<vector<vector<char> > > GFTable::m_seenCellCombos;
  //vector<pair<int, int> > GFTable::m_comboOffsets;

  GFTable::GFTable() :
      SpectrumItem(), m_table(0), m_firstPM(0), m_aaMasses(0), m_aaJumps(0x0), m_endMasses(), m_sizes(0), m_useProfile(false), m_nodeScores(0), m_maxScores(0), m_minScores(0), m_edgeScores(0), m_startAlign(0), m_endAlign(0)
  {
  }

  GFTable::GFTable(const GFTable &other) :
      SpectrumItem(), m_table(0), m_firstPM(0), m_aaMasses(0), m_aaJumps(0x0), m_endMasses(), m_sizes(0), m_useProfile(false), m_nodeScores(0), m_maxScores(0), m_minScores(0), m_edgeScores(0), m_startAlign(0), m_endAlign(0)
  {
    this->operator =(other);
  }

  GFTable::GFTable(const Spectrum &spectrum,
                   const AAJumps &jumps,
                   const Spectrum &msSpectrum,
                   const MassErrorModel &massErrorModelTarget,
                   const MassErrorModel &massErrorModelDecoy,
                   const float maxPPMError) :
      SpectrumItem(), m_table(spectrum.size()), m_firstPM(0), m_aaMasses(0), m_aaJumps(0x0), m_endMasses(), m_sizes(0), m_useProfile(false), m_nodeScores(0), m_maxScores(0), m_minScores(0), m_edgeScores(0), m_startAlign(0), m_endAlign(0)
  {
    initialize(spectrum,
               jumps,
               msSpectrum,
               massErrorModelTarget,
               massErrorModelDecoy,
               maxPPMError);
  }

  void GFTable::initialize(const Spectrum &spectrum,
                           const AAJumps &jumps,
                           const Spectrum &msSpectrum,
                           const MassErrorModel &massErrorModelTarget,
                           const MassErrorModel &massErrorModelDecoy,
                           const float maxPPMError)
  {
    if (spectrum.parentMass < AAJumps::massMH)
    {
      WARN_MSG("bad spectrum with parent mass " << spectrum.parentMass << ", skipping GF calculation");
      return;
    }
    if (spectrum.size() == 0)
    {
      clear();
      return;
    }
    //m_binnedSpec = &spectrum;

    const bool ignoreMS2ErrorModel = (massErrorModelTarget.size() == 0
        || massErrorModelDecoy.size() == 0 || msSpectrum.size() == 0);

    m_aaMasses.resize(jumps.size());
    for (int aa = 0; aa < jumps.size(); aa++)
    {
      // Load all rounded AA masses
      m_aaMasses[aa] = floatToInt(jumps[aa]);
    }
    //spectrum.output(cerr);

    m_scan = spectrum.scan;
    m_filename = spectrum.fileName;
    m_index = -1;
    m_aaJumps = &jumps;
    m_startAlign = 0;
    m_endAlign = 0;

    // compute size of node and edge vectors so we can lookup the parent mass at the end each
    // Allow for +/- 1 Da parent mass errors
    int massBinsSize = floatToInt(spectrum.parentMass + spectrum.parentMassTol
        - round(AAJumps::massMH * AA_ROUNDING)) + 1;

    float pepMassUse = spectrum.parentMass - (AAJumps::massMH * AA_ROUNDING);

    m_firstPM = massBinsSize - 1;

    //DEBUG_VAR(massBinsSize);
    //DEBUG_VAR(m_firstPM);

    m_endMasses.clear();
    for (int m = m_firstPM; m < massBinsSize; m++)
    {
      m_endMasses.insert(m);
    }

    m_sizes.resize(massBinsSize);

    m_nodeScores.assign(m_sizes.size(), 0);

    if (massBinsSize > m_table.size())
    {
      m_table.resize(massBinsSize);
    }

    if (massBinsSize > m_edgeScores.size())
    {
      m_edgeScores.resize(massBinsSize);
    }

    vector<float> accMasses(m_nodeScores.size(), -1.0);

    for (int p = 0; p < spectrum.size(); p++)
    {
      int mass = floatToInt(spectrum[p][0]);
      if (mass > 0 && mass < m_firstPM)
      {
        m_nodeScores[mass] += floatToInt(spectrum[p][1]);
      }
      if (mass >= 0 && mass < accMasses.size())
      {
        accMasses[mass] = spectrum[p][0];
      }
    }

    const float absentEdgeScore =
        (ignoreMS2ErrorModel) ? 0 :
            MassErrorModel::getAbsentEdgeScore(massErrorModelTarget,
                                               massErrorModelDecoy);

    //DEBUG_VAR(MassErrorModel::getAbsentEdgeScore(massErrorModelTarget, massErrorModelDecoy));

    //DEBUG_VAR(MassErrorModel::getDangleEdgeScore(massErrorModelTarget, massErrorModelDecoy));

    float minEdgeScore = absentEdgeScore;

    float maxEdgeScore = absentEdgeScore;

    vector<vector<float> > rawEdgeScores(massBinsSize);

    for (int m = 0; m < massBinsSize; m++)
    {
      rawEdgeScores[m].assign(m_aaMasses.size(), absentEdgeScore);

      for (int aa = 0; aa < m_aaMasses.size(); aa++)
      {
        int prevMass = m - m_aaMasses[aa];

        if (prevMass < 0)
        {
          continue;
        }

        rawEdgeScores[m][aa] =
            (ignoreMS2ErrorModel) ? absentEdgeScore :
                MassErrorModel::getEdgeScore(msSpectrum,
                                             spectrum,
                                             accMasses[prevMass],
                                             pepMassUse - accMasses[prevMass],
                                             accMasses[m],
                                             pepMassUse - accMasses[m],
                                             jumps[aa],
                                             maxPPMError,
                                             massErrorModelTarget,
                                             massErrorModelDecoy,
                                             AA_ROUNDING);

        minEdgeScore = min(minEdgeScore, rawEdgeScores[m][aa]);

        maxEdgeScore = max(maxEdgeScore, rawEdgeScores[m][aa]);
      }
    }

    for (int m = 0; m < massBinsSize; m++)
    {
      for (int aa = 0; aa < m_aaMasses.size(); aa++)
      {
        rawEdgeScores[m][aa] -= minEdgeScore;
      }
    }

    //scale edge scores to EDGE_SCORE_SCALING_FACTOR % of max node score path
    vector<float> maxScores(massBinsSize, -1.0);
    maxScores[0] = 0;
    for (int m = 0; m < massBinsSize; m++)
    {
      if (maxScores[m] < 0)
      {
        continue;
      }
      for (int aa = 0; aa < m_aaMasses.size(); aa++)
      {
        int nextMass = m + m_aaMasses[aa];
        if (nextMass >= massBinsSize)
        {
          continue;
        }
        float pathScore = maxScores[m] + (float)m_nodeScores[nextMass];

        maxScores[nextMass] = max(maxScores[nextMass], pathScore);
      }
    }

    if (maxScores[m_firstPM] < 0)
    {
      WARN_MSG("Could not find path in spectrum graph for scan " << spectrum.scan << ", skipping GF calculation");
      clear();
      return;
    }

    const float maxNodeScorePath = maxScores[m_firstPM];
    maxScores.assign(massBinsSize, -1.0);
    maxScores[0] = 0;
    for (int m = 0; m < massBinsSize; m++)
    {
      if (maxScores[m] < 0)
      {
        continue;
      }
      for (int aa = 0; aa < m_aaMasses.size(); aa++)
      {
        int nextMass = m + m_aaMasses[aa];
        if (nextMass >= massBinsSize)
        {
          continue;
        }
        float pathScore = maxScores[m] + (float)rawEdgeScores[nextMass][aa];

        maxScores[nextMass] = max(maxScores[nextMass], pathScore);
      }
    }

    const float maxEdgeScorePath = maxScores[m_firstPM];
    const float edgeScaleFactor =
        (maxNodeScorePath < 0.1 || maxEdgeScorePath < 0.1) ? 3.0 :
            EDGE_SCORE_SCALING_FACTOR * (maxNodeScorePath / maxEdgeScorePath);

    if (DEBUG)
    {
      DEBUG_VAR(maxNodeScorePath);
      DEBUG_VAR(maxEdgeScorePath);
      DEBUG_VAR(edgeScaleFactor);
      DEBUG_VAR(minEdgeScore);
      DEBUG_VAR(maxEdgeScore);
    }

    for (int m = 0; m < massBinsSize; m++)
    {
      m_edgeScores[m].resize(m_aaMasses.size());
      for (int aa = 0; aa < m_aaMasses.size(); aa++)
      {
        m_edgeScores[m][aa] = floatToInt(rawEdgeScores[m][aa]
            * edgeScaleFactor);
      }
    }

    /*
     m_edgeScores.assign(4, 0);

     if (edgeScores != NULL)
     {
     int minEdgeScore = 0;
     for (int e = 0; e < edgeScores->size(); e++)
     {
     minEdgeScore = min(minEdgeScore, (*edgeScores)[e]);
     }

     if (edgeScores->size() != 4)
     {
     ERROR_MSG("Invalid edge score vector size (" << edgeScores->size() << ")");
     }

     for (int e = 0; e < edgeScores->size(); e++)
     {
     m_edgeScores[e] = (*edgeScores)[e] - minEdgeScore;
     }

     //DEBUG_MSG("nn = " << m_edgeScores[0] << ", yn = " << m_edgeScores[1] << ", ny = " << m_edgeScores[2] << ", yy = " << m_edgeScores[3]);
     }
     */

    m_minScores.resize(0);
    m_maxScores.resize(0);
  }

  void GFTable::addPRMScores(const int shift, const Spectrum &spectrum)
  {
    m_startAlign = max(0, shift);
    int end2 = floatToInt(spectrum.parentMass + spectrum.parentMassTol
        - round(AAJumps::massMH * AA_ROUNDING));
    m_endAlign = min(m_firstPM, end2 + shift);
    for (int p = 0; p < spectrum.size(); p++)
    {
      int mass = shift + floatToInt(spectrum[p][0]);
      if (mass > 0 && mass < m_firstPM)
      {
        m_nodeScores[mass] += floatToInt(spectrum[p][1]);
      }
    }
  }

  void GFTable::addMatchingPRMScores(const int shift, const Spectrum &spectrum)
  {
    m_startAlign = max(0, shift);
    int end2 = floatToInt(spectrum.parentMass + spectrum.parentMassTol
        - round(AAJumps::massMH * AA_ROUNDING));
    m_endAlign = min(m_firstPM, end2 + shift);
    for (int p = 0; p < spectrum.size(); p++)
    {
      int mass = shift + floatToInt(spectrum[p][0]);
      if (mass > 0 && mass < m_firstPM && m_nodeScores[mass] > 0)
      {
        m_nodeScores[mass] += floatToInt(spectrum[p][1]);
      }
    }
  }

  void GFTable::computeGeneratingFunction(const int minPathScore,
                                          const int startMass,
                                          const bool useProfile)
  {
    optimizeDimensions(minPathScore, startMass);
    m_useProfile = useProfile;

    if (useProfile)
    {
      if (m_sizes.size() > m_tableProfile.size())
      {
        m_tableProfile.resize(m_sizes.size());
      }
    }
    else
    {
      if (m_sizes.size() > m_table.size())
      {
        m_table.resize(m_sizes.size());
      }
    }

    GFProfileCell defaultNodeProfile;
    GFCell defaultCell;

    for (int m = startMass; m < m_sizes.size(); m++)
    {
      if (useProfile)
      {
        if (m_sizes[m] > m_tableProfile[m].size())
        {
          m_tableProfile[m].resize(m_sizes[m]);
        }
      }
      else
      {
        if (m_sizes[m] > m_table[m].size())
        {
          m_table[m].resize(m_sizes[m]);
        }
      }

      for (int s = 0; s < m_sizes[m]; s++)
      {
        if (useProfile)
        {
          m_tableProfile[m][s] = defaultNodeProfile;
        }
        else
        {
          m_table[m][s] = defaultCell;
        }
      }
    }

    if (useProfile)
    {
      this->m_computeProbabilitiesProfile(startMass);
    }
    else
    {
      this->m_computeProbabilities(startMass);
    }
  }

  void GFTable::clear()
  {
    m_sizes.resize(0);
    m_edgeScores.resize(0);
    m_nodeScores.resize(0);
    m_firstPM = 0;
    m_endMasses.clear();
    m_aaMasses.resize(0);
    m_aaJumps = 0;
  }

  GFTable &GFTable::operator=(const GFTable &other)
  {
    if (this == &other)
    {
      return *this;
    }

    SpectrumItem::operator =((const SpectrumItem &)other);
    m_table = other.m_table;
    m_tableProfile = other.m_tableProfile;
    m_sizes = other.m_sizes;
    m_nodeScores = other.m_nodeScores;
    m_firstPM = other.m_firstPM;
    m_aaMasses = other.m_aaMasses;
    m_aaJumps = other.m_aaJumps;
    m_endMasses = other.m_endMasses;
    return *this;
  }

  int GFTable::getScoreThreshold(const double probThreshold) const
  {
    double totalSpecProb = 0;
    int maxScore = m_sizes[m_sizes.size() - 1] - 1;
    int cutoffScore = 0;

    // Find the maximum score that has a total probability >= to the threshold
    for (int s = maxScore; s >= 0; s--)
    {
      cutoffScore = s;
      for (set<int>::const_iterator mIt = m_endMasses.begin();
          mIt != m_endMasses.end(); mIt++)
      {
        const int m = *mIt;
        if (s >= m_sizes[m])
        {
          continue;
        }

        if (m_useProfile)
        {
          totalSpecProb += getProbability(m_tableProfile[m][s].probability);
        }
        else
        {
          totalSpecProb += getProbability(m_table[m][s].probability);
        }
      }
      if (totalSpecProb >= probThreshold)
      {
        break;
      }
    }
    return cutoffScore + m_minScores[m_sizes.size() - 1];
  }

  int GFTable::encodeDictionary(const int cutoffScore)
  {
    if (DEBUG)
    {
      DEBUG_VAR(cutoffScore);
    }

    if (cutoffScore == 0)
    {
      return cutoffScore;
    }

    m_exploreBackwards(cutoffScore);

    return cutoffScore;
  }

  pair<int, int> GFTable::lookupPeptide(const int startMass,
                                        const int startScore,
                                        const string &peptide) const
  {
    if (m_aaJumps == 0x0)
    {
      ERROR_MSG("Amino acid jumps have not been initialized!");
      abort();
    }
    vector<float> prmMasses;
    vector<string> singleAAJumps;

    if (!m_aaJumps->getRoundedPRMMasses(peptide, prmMasses, 0, &singleAAJumps))
    {
      ERROR_MSG("Cannot parse peptide sequence \'" << peptide << "\'");
      abort();
    }
    if (DEBUG)
    {
      DEBUG_VAR(prmMasses[prmMasses.size() - 1]);
      DEBUG_VAR(m_sizes.size());
    }

    int curMass = startMass;
    int nextMass;
    int curScore = startScore;
    int jumpIdx = 0;
    while (curMass < m_sizes.size() && m_endMasses.count(curMass) == 0
        && jumpIdx < singleAAJumps.size())
    {
      const string &curJumpStr = singleAAJumps[jumpIdx];
      if (DEBUG)
      {
        DEBUG_VAR(curJumpStr);
      }
      if (m_aaJumps->massLookup.count(curJumpStr) == 0)
      {
        ERROR_MSG("Failed to lookup jump \'" << curJumpStr << "\'!");
        abort();
      }
      const int &aa = (int)m_aaJumps->massLookup.at(curJumpStr);

      nextMass = curMass + m_aaMasses[aa];

      if (nextMass >= m_sizes.size())
      {
        break;
      }

      //const int nextMass = mass + m_aaMasses[nextAA];
      //return score + m_nodeScores[nextMass] + m_edgeScores[nextMass][nextAA];

      curScore = getNextScore(curMass, curScore, aa);

      curMass = nextMass;

      if (DEBUG)
      {
        DEBUG_MSG(curMass << " = " << m_nodeScores[curMass] << " ( " << curScore << " ) ");
      }

      jumpIdx++;
    }
    if (DEBUG)
    {
      DEBUG_VAR(curScore);
    }
    return pair<int, int>(curMass, curScore);
  }

  int GFTable::getMatchScore(const string &peptide,
                             int *startMass,
                             int *endMass) const
  {
    if (size() == 0)
    {
      return 0;
    }

    int maxScore = 0;
    if (startMass != 0)
    {
      *startMass = 0;
    }
    if (endMass != 0)
    {
      *endMass = 0;
    }

    pair<int, int> endPt = lookupPeptide(0, m_nodeScores[0], peptide);

    if (m_endMasses.count(endPt.first) > 0 && endPt.second > maxScore)
    {
      maxScore = endPt.second;
      if (startMass != 0)
      {
        *startMass = 0;
      }
      if (endMass != 0)
      {
        *endMass = endPt.first;
      }
    }
    return maxScore;
  }

  double GFTable::getPValue(const int score) const
  {
    double pvalue = 0;
    for (set<int>::const_iterator mIt = m_endMasses.begin();
        mIt != m_endMasses.end(); mIt++)
    {
      const int m = *mIt;
      for (int s = score - m_minScores[m]; s < m_sizes[m]; s++)
      {
        if (DEBUG)
        {
          DEBUG_MSG("m = " << m << ", s = " << s);
          DEBUG_VAR(getProbability(m_table[m][s].probability));
        }

        if (m_useProfile)
        {
          pvalue += getProbability(m_tableProfile[m][s].probability);
        }
        else
        {
          pvalue += getProbability(m_table[m][s].probability);
        }
      }
    }

    return pvalue;
  }

  void GFTable::optimizeDimensions(const int minPathScore, const int startMass)
  {
    if (DEBUG)
    {
      DEBUG_VAR(minPathScore);
      DEBUG_VAR(startMass);
    }
    m_computeMinScores(minPathScore);
    m_computeMaxScores(startMass);

    for (int m = 0; m < m_sizes.size(); m++)
    {
      if (m_maxScores[m] < 0 || m_minScores[m] < 0)
      {
        m_sizes[m] = 0;
        continue;
      }

      // Each 2nd degree vector only needs to be large enough to accommodate all paths ending with minPathScore or higher
      m_sizes[m] = m_maxScores[m] + 1 - m_minScores[m];

      if (DEBUG)
      {
        DEBUG_MSG("Bin " << m << " has scores range [" << m_minScores[m] << " , " << m_maxScores[m] << "]");
      }
    }
  }

  void GFTable::m_exploreBackwards(const int cutoffScore)
  {
    if (!m_useProfile)
    {
      return;
    }

    for (set<int>::const_iterator mIt = m_endMasses.begin();
        mIt != m_endMasses.end(); mIt++)
    {
      const int m = *mIt;
      for (int s = m_sizes[m] - 1; s >= cutoffScore; s--)
      {
        if (m_tableProfile[m][s].prefixes > 0)
        {
          m_tableProfile[m][s].suffixes = 1;
        }
      }
    }

    int nextMass, nextScore;

    for (int m = m_sizes.size() - 1; m >= 0; m--)
    {
      for (int s = m_sizes[m] - 1; s >= 0; s--)
      {

        if (m_tableProfile[m][s].suffixes == 0
            || m_tableProfile[m][s].prefixes == 0)
        {
          m_tableProfile[m][s].prefixes = 0;
          m_tableProfile[m][s].probability = 0;
          continue;
        }

        for (int aa = 0; aa < m_aaMasses.size(); aa++)
        {
          nextMass = m - m_aaMasses[aa];

          if (nextMass < 0)
          {
            continue;
          }

          nextScore = s - m_nodeScores[m];
          if (nextScore < 0 || nextScore >= m_sizes[nextMass]
              || m_tableProfile[nextMass][nextScore].prefixes == 0)
          { // skip nodes that don't exist in the GF tables
            continue;
          }

          m_tableProfile[nextMass][nextScore].suffixes +=
              m_tableProfile[m][s].suffixes;
        }
      }
    }
  }

  void GFTable::m_computeProbabilities(const int startMass)
  {
    int prevScore, currentNodeScore, prevMass, currentPathScore;

    float aaProb = 1.0 / (float)m_aaMasses.size();

    GFCell defaultCell;

    // Compute spectral probabilities
    for (int m = startMass; m < m_sizes.size(); m++)
    {
      currentNodeScore = m_nodeScores[m];

      //DEBUG_VAR(m << " = " << currentNodeScore);

      if (m == startMass)
      {
        currentNodeScore = 0;

        for (int s = 0; s < currentNodeScore; s++)
        {
          m_table[m][s] = defaultCell;
        }

        m_table[m][currentNodeScore] = defaultCell;
        // Initialize base case of the recursion
        setProbability(m_table[m][currentNodeScore].probability, 1.0);

        if (DEBUG)
          DEBUG_MSG("Start (" << m << "," << currentNodeScore << ") = " << getProbability(m_table[m][currentNodeScore].probability));

        for (int s = currentNodeScore + 1; s < m_sizes[m]; s++)
        {
          m_table[m][s] = defaultCell;
        }
        continue;
      }

      for (int s = 0; s < m_sizes[m]; s++)
      {
        m_table[m][s] = defaultCell;
        currentPathScore = s + m_minScores[m];

        if (currentPathScore - currentNodeScore < 0)
        {
          continue;
        }

        // if (DEBUG)
        //  DEBUG_MSG("(" << m << "," << s << ")");

        float incProb = 0;
        int incEdges = 0;

        for (int aa = 0; aa < m_aaMasses.size(); aa++)
        {
          prevMass = m - m_aaMasses[aa];

          // ignore invalid mass/score combinations
          if (prevMass < 0 || (prevMass < m_startAlign && m > m_startAlign)
              || (prevMass < m_endAlign && m > m_endAlign))
          {
            continue;
          }

          prevScore = getPreviousScore(m, currentPathScore, aa)
              - m_minScores[prevMass];

          if (prevScore < 0 || prevScore >= m_sizes[prevMass])
          {
            continue;
          }

          if (m_table[prevMass][prevScore].probability > 0)
          {
            // add up the probability
            incrementProbability(incProb,
                                 m_table[prevMass][prevScore].probability,
                                 aaProb);

            if (DEBUG)
              DEBUG_MSG(" (m= " << prevMass << " ,s= " << prevScore + m_minScores[prevMass] << " ) -" << m_aaJumps->getLabel(aa) << "-> (m= " << m << " ,s= " << s + m_minScores[m] << " )");

            // set the backward edge
            setEdge(incEdges, aa, true);
          }
        }

        if (incProb > 0)
        {
          m_table[m][s].probability = incProb;
          m_table[m][s].rEdges = incEdges;
        }
      }
    }
  }

  void GFTable::m_computeProbabilitiesProfile(const int startMass)
  {
    int prevScore, currentNodeScore, prevMass, currentPathScore;

    float aaProb = 1.0 / (float)m_aaMasses.size();

    // Compute spectral probabilities
    for (int m = startMass; m < m_sizes.size(); m++)
    {
      currentNodeScore = m_nodeScores[m];

      if (m == startMass)
      {
        currentNodeScore = 0;

        for (int s = 0; s < currentNodeScore; s++)
        {
          m_tableProfile[m][s].probability = 0;
          m_tableProfile[m][s].prefixes = 0;
          m_tableProfile[m][s].suffixes = 0;
        }

        // Initialize base case of the recursion
        setProbability(m_tableProfile[m][currentNodeScore].probability, 1.0);
        m_tableProfile[m][currentNodeScore].prefixes = 1;
        m_tableProfile[m][currentNodeScore].suffixes = 0;

        for (int s = currentNodeScore + 1; s < m_sizes[m]; s++)
        {
          m_tableProfile[m][s].probability = 0;
          m_tableProfile[m][s].prefixes = 0;
          m_tableProfile[m][s].suffixes = 0;
        }
        continue;
      }

      for (int s = 0; s < m_sizes[m]; s++)
      {
        m_tableProfile[m][s].probability = 0;
        m_tableProfile[m][s].prefixes = 0;
        m_tableProfile[m][s].suffixes = 0;

        currentPathScore = s + m_minScores[m];

        if (currentPathScore - currentNodeScore < 0)
        {
          continue;
        }

        for (int aa = 0; aa < m_aaMasses.size(); aa++)
        {
          prevMass = m - m_aaMasses[aa];

          // ignore invalid mass/score combinations
          if (prevMass < 0)
          {
            continue;
          }

          prevScore = getPreviousScore(m, currentPathScore, aa)
              - m_minScores[prevMass];

          if (prevScore < 0 || prevScore >= m_sizes[prevMass])
          {
            continue;
          }

          // add up the probability
          incrementProbability(m_tableProfile[m][s].probability,
                               m_tableProfile[prevMass][prevScore].probability,
                               aaProb);

          m_tableProfile[m][s].prefixes +=
              m_tableProfile[prevMass][prevScore].prefixes;
        }
      }
    }
  }

  void GFTable::m_computeMinScores(const int minPathScore)
  {

    if (minPathScore <= 0)
    {
      m_minScores.assign(m_sizes.size(), 0);
      return;
    }
    m_minScores.assign(m_sizes.size(), -1);

    for (int m = m_sizes.size() - 1; m >= 0; m--)
    {
      if (m >= m_firstPM)
      {
        m_minScores[m] = minPathScore;
      }
      else if (m_minScores[m] < 0)
      {
        continue;
      }

      for (int aa = 0; aa < m_aaMasses.size(); aa++)
      {
        int prevScore = max(0, getPreviousScore(m, m_minScores[m], aa));

        int prevMass = m - m_aaMasses[aa];

        if (prevMass < 0)
        {
          continue;
        }

        if (m_minScores[prevMass] < 0)
        {
          m_minScores[prevMass] = prevScore;
        }
        else
        {
          m_minScores[prevMass] = min(m_minScores[prevMass], prevScore);
        }
      }
      if (DEBUG)
      {
        DEBUG_MSG("Min score at " << m << " is " << m_minScores[m]);
      }
    }
  }

  void GFTable::m_computeMaxScores(const int startMass)
  {
    m_maxScores.assign(m_sizes.size(), -1);

    for (int m = startMass; m < m_sizes.size(); m++)
    {
      if (m == startMass)
      {
        m_maxScores[m] = 0;        //m_nodeScores[m];
      }
      else if (m_maxScores[m] < 0 || m_minScores[m] < 0)
      {
        continue;
      }

      if (m_minScores[m] > m_maxScores[m])
      {
        m_maxScores[m] = -1;
        continue;
      }

      for (int aa = 0; aa < m_aaMasses.size(); aa++)
      {
        int nextMass = m + m_aaMasses[aa];
        if (nextMass >= m_sizes.size())
        {
          continue;
        }
        int pathScore = getNextScore(m, m_maxScores[m], aa);

        m_maxScores[nextMass] = max(m_maxScores[nextMass], pathScore);
      }
    }
  }
}
