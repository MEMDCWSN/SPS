/*
 * GFTableIntersection.cpp
 *
 *  Created on: Sep 5, 2013
 *      Author: aguthals
 */

#include "GFTableIntersection.h"

using namespace std;

namespace specnets
{
  bool GFTableIntersection::DEBUG = false;

  const int GFTableIntersection::TABLE_SLICE_BUFFER = 50;

  GFTableIntersection::GFTableIntersection() :
      m_table(0), m_aaJumps(0x0), m_table1(0x0), m_table2(0x0), m_lambdaShift(0), m_begin(0), m_end(0), m_sizes1(0), m_minScores1(0), m_sizes2(0), m_minScores2(0), m_startScore1(0), m_sliceSz(0), m_endDist1(0), m_endDist2(0), m_sliceStartPos(0), m_maxAAMass(0), m_endMass1(0), m_endMass2(0), endScore1Ovlp(0)
  {
  }

  GFTableIntersection::GFTableIntersection(const GFTable &table1,
                                           const GFTable &table2,
                                           const int shift,
                                           const int matchScoreOvlp) :
      m_table(0), m_aaJumps(0x0), m_table1(0x0), m_table2(0x0), m_lambdaShift(0), m_begin(0), m_end(0), m_sizes1(0), m_minScores1(0), m_sizes2(0), m_minScores2(0), m_startScore1(0), m_sliceSz(0), m_endDist1(0), m_endDist2(0), m_sliceStartPos(0), m_maxAAMass(0), m_endMass1(0), m_endMass2(0), endScore1Ovlp(0)
  {
    initialize(table1, table2, shift, matchScoreOvlp);
  }

  void GFTableIntersection::initialize(const GFTable &table1,
                                       const GFTable &table2,
                                       const int shift,
                                       const int matchScoreOvlp)
  {

    if (table1.size() == 0 || table2.size() == 0)
    {
      WARN_MSG("Cannot intersect tables with sizes (" << table1.size() << ", " << table2.size() << ")");
      clear();
      return;
    }

    if (shift < 0)
    {
      ERROR_MSG("Shift must be positive!!");
      clear();
      return;
    }

    m_table1 = &table1;
    m_table2 = &table2;

    m_aaJumps = table1.getAAJumps();

    m_lambdaShift = shift;
    endScore1Ovlp = matchScoreOvlp;

    m_begin = max(0, m_lambdaShift);
    m_end = min(*m_table1->m_endMasses.rbegin(),
                *m_table2->m_endMasses.rbegin() + m_lambdaShift);

    m_startScore1 = 0;

    //DEBUG_VAR(m_begin);
    //DEBUG_VAR(m_end);
    //DEBUG_VAR(m_lambdaShift);

    const int tableSz = max((int)table1.size(),
                            ((int)table2.size()) + m_lambdaShift);

    //DEBUG_VAR(tableSz);

    m_sizes1.assign(tableSz, pair<int, int>(0, 0));
    m_minScores1.assign(tableSz, pair<int, int>(0, 0));

    m_sizes2.assign(tableSz, pair<int, int>(0, 0));
    m_minScores2.assign(tableSz, pair<int, int>(0, 0));

    m_endDist1.resize(0);
    m_endDist2.resize(0);
  }

  void GFTableIntersection::computeGeneratingFunction(const bool useProfile)
  {
    //DEBUG_TRACE;
    m_setDimensions();
    m_initializeSlice();

    m_endDist1.resize(0);
    m_endDist2.resize(0);
    //DEBUG_TRACE;
    const vector<int> &aaMasses = m_table1->m_aaMasses;

    const vector<int> &sizes1 = m_table1->m_sizes;
    const vector<int> &sizes2 = m_table2->m_sizes;

    int scoreDim1, scoreDim2;

    int prevScore1, currentScore1, prevScore2, currentScore2, prevMass1,
        prevMass2, m2, currentPathScore1, currentPathScore2, mOffset,
        prevOffset;

    int end2 = min(*m_table2->m_endMasses.rbegin(),
                   *m_table1->m_endMasses.rbegin() - m_lambdaShift);

    float aaProb = 1.0 / (float)aaMasses.size();

    if (DEBUG)
    {
      DEBUG_VAR(m_sizes1.size());
    }

    m_startScore1 = -1;

    m_endMass1 = *m_table1->m_endMasses.rbegin();

    m_endMass2 = *m_table2->m_endMasses.rbegin() + m_lambdaShift;

    //DEBUG_TRACE;
    /*
     * 2D iteration for beginning of table 1
     */
    for (int m1 = 0; m1 <= m_begin; m1++)
    {
      scoreDim1 = m_sizes1[m1].first;
      scoreDim2 = m_sizes1[m1].second;

      mOffset = m_getSliceOffset(m1);

      if (scoreDim1 == 0 || scoreDim2 == 0)
      {
        continue;
      }

      if (m1 == 0)
      {
        GFTable::setProbability(m_table[mOffset][0][0].probability1, 1.0);
        if (DEBUG)
        {
          DEBUG_MSG("Base case (1) - 0,0,0 == 1");
        }
      }

      if (DEBUG)
      {
        DEBUG_MSG("m1 = " << m1);
      }
      currentScore1 = m_table1->getScore(m1);

      for (int s1 = 0; s1 < scoreDim1; s1++)
      {
        currentPathScore1 = s1 + m_minScores1[m1].first;

        if (currentPathScore1 - currentScore1 < 0)
        {
          continue;
        }

        float incProb = 0;

        for (int aa = 0; aa < aaMasses.size(); aa++)
        {
          prevMass1 = m1 - aaMasses[aa];

          // ignore invalid mass/score combinations
          if (prevMass1 < 0)
          {
            continue;
          }

          prevOffset = m_getSliceOffset(prevMass1);

          prevScore1 = m_table1->getPreviousScore(m1, currentPathScore1, aa)
              - m_minScores1[prevMass1].first;

          if (prevScore1 < 0 || prevScore1 >= m_sizes1[prevMass1].first)
          {
            continue;
          }

          // add up the probability
          GFTable::incrementProbability(incProb,
                                        m_table[prevOffset][prevScore1][0].probability1,
                                        aaProb);
        }

        if (incProb > 0)
        {
          m_table[mOffset][s1][0].probability1 = incProb;
        }
      }

      if (m1 == m_begin)
      {
        m_startScore1 = 0;

        GFTable::setProbability(m_table[mOffset][m_startScore1][0].probability2,
                                1.0);

        if (DEBUG)
        {
          DEBUG_MSG("Base case (2) - " << m1 << "," << m_startScore1 << ",0 == 1");
        }
      }
    }

    //DEBUG_TRACE;
    /**
     * 3D Iteration
     */
    for (int m1 = m_begin + 1; m1 <= m_end; m1++)
    {
      scoreDim1 = m_sizes1[m1].first;
      scoreDim2 = m_sizes1[m1].second;

      currentScore1 = m_table1->getScore(m1);
      m2 = m1 - m_lambdaShift;
      currentScore2 = m_table2->getScore(m2);

      mOffset = m_getSliceOffset(m1);

      if (m_endDist1.size() == 0)
      {
        for (int s1 = 0; s1 < scoreDim1; s1++)
        {
          currentPathScore1 = s1 + m_minScores1[m1].first;

          for (int s2 = 0; s2 < scoreDim2; s2++)
          {

            currentPathScore2 = s2 + m_minScores1[m1].second;

            if (currentPathScore1 - currentScore1 < 0
                || currentPathScore2 - currentScore2 < 0)
            {
              continue;
            }

            float incProb1 = 0;
            for (int aa = 0; aa < aaMasses.size(); aa++)
            {
              prevMass1 = m1 - aaMasses[aa];

              // ignore invalid mass/score combinations
              if (prevMass1 < 0 || prevMass1 < m_begin)
              {
                continue;
              }

              prevMass2 = m2 - aaMasses[aa];

              if (prevMass2 < 0)
              {
                continue;
              }

              prevOffset = m_getSliceOffset(prevMass1);

              prevScore1 = m_table1->getPreviousScore(m1, currentPathScore1, aa)
                  - m_minScores1[prevMass1].first;

              if (prevScore1 < 0 || prevScore1 >= m_sizes1[prevMass1].first)
              {
                continue;
              }

              prevScore2 = m_table2->getPreviousScore(m2, currentPathScore2, aa)
                  - m_minScores1[prevMass1].second;

              if (prevScore2 < 0 || prevScore2 >= m_sizes1[prevMass1].second)
              {
                continue;
              }

              GFTable::incrementProbability(incProb1,
                                            m_table[prevOffset][prevScore1][prevScore2].probability1,
                                            aaProb);
            }
            if (incProb1 > 0)
            {
              m_table[mOffset][s1][s2].probability1 = incProb1;
            }
          }
        }
        if (m1 == m_endMass1)
        {
          m_endDist1 = m_table[mOffset];
        }
      }

      // Need to compute the second probability in separate dimensions since these peptides may begin matching later in spectrum 1
      scoreDim1 = m_sizes2[m1].first;
      scoreDim2 = m_sizes2[m1].second;

      if (m_endDist2.size() == 0)
      {
        for (int s1 = 0; s1 < scoreDim1; s1++)
        {
          currentPathScore1 = s1 + m_minScores2[m1].first;

          for (int s2 = 0; s2 < scoreDim2; s2++)
          {

            currentPathScore2 = s2 + m_minScores2[m1].second;

            if (currentPathScore1 - currentScore1 < 0
                || currentPathScore2 - currentScore2 < 0)
            {
              continue;
            }

            float incProb2 = 0;
            for (int aa = 0; aa < aaMasses.size(); aa++)
            {
              prevMass1 = m1 - aaMasses[aa];

              // ignore invalid mass/score combinations
              if (prevMass1 < 0 || prevMass1 < m_begin)
              {
                continue;
              }

              prevMass2 = m2 - aaMasses[aa];

              if (prevMass2 < 0)
              {
                continue;
              }

              prevOffset = m_getSliceOffset(prevMass1);

              prevScore1 = m_table1->getPreviousScore(m1, currentPathScore1, aa)
                  - m_minScores2[prevMass1].first;

              if (prevScore1 < 0 || prevScore1 >= m_sizes2[prevMass1].first)
              {
                continue;
              }

              prevScore2 = m_table2->getPreviousScore(m2, currentPathScore2, aa)
                  - m_minScores2[prevMass1].second;

              if (prevScore2 < 0 || prevScore2 >= m_sizes2[prevMass1].second)
              {
                continue;
              }

              GFTable::incrementProbability(incProb2,
                                            m_table[prevOffset][prevScore1][prevScore2].probability2,
                                            aaProb);
            }
            if (incProb2 > 0)
            {
              if (DEBUG && m1 == m_end)
              {
                DEBUG_MSG("End ovlp (2) : m1 = " << m1 << ", s1 = " << s1 << ", s2 = " << s2);
              }
              m_table[mOffset][s1][s2].probability2 = incProb2;
            }
          }
        }

        if (m1 == m_endMass2)
        {
          m_endDist2 = m_table[mOffset];
        }
      }
    }
    //DEBUG_TRACE;
    /**
     * Remaining 2D iteration for table 1
     */
    for (int m1 = m_end + 1; m1 < m_sizes1.size() && m1 < sizes1.size(); m1++)
    {
      scoreDim1 = m_sizes1[m1].first;
      scoreDim2 = m_sizes1[m1].second;

      currentScore1 = m_table1->getScore(m1);

      mOffset = m_getSliceOffset(m1);

      for (int s1 = 0; s1 < scoreDim1; s1++)
      {
        currentPathScore1 = s1 + m_minScores1[m1].first;

        if (currentPathScore1 - currentScore1 < 0)
        {
          continue;
        }

        for (int s2 = 0; s2 < scoreDim2; s2++)
        {

          float incProb1 = 0, incProb2 = 0;

          for (int aa = 0; aa < aaMasses.size(); aa++)
          {
            prevMass1 = m1 - aaMasses[aa];

            // ignore invalid mass/score combinations
            if (prevMass1 < 0 || prevMass1 < m_end)
            {
              continue;
            }

            prevScore1 = m_table1->getPreviousScore(m1, currentPathScore1, aa)
                - m_minScores1[prevMass1].first;

            if (prevScore1 < 0 || prevScore1 >= m_sizes1[prevMass1].first)
            {
              continue;
            }

            prevOffset = m_getSliceOffset(prevMass1);

            GFTable::incrementProbability(incProb1,
                                          m_table[prevOffset][prevScore1][s2].probability1,
                                          aaProb);

          }

          if (incProb1 > 0)
          {
            m_table[mOffset][s1][s2].probability1 = incProb1;
          }
        }
      }

      if (m1 == m_endMass1)
      {
        m_endDist1 = m_table[mOffset];
      }
    }
    //DEBUG_TRACE;
    /**
     * Remaining 2D iteration for table 2
     */
    for (int m1 = sizes1.size();
        m1 < m_sizes1.size() && m1 - m_lambdaShift < sizes2.size(); m1++)
    {
      scoreDim1 = m_sizes2[m1].first;
      scoreDim2 = m_sizes2[m1].second;

      m2 = m1 - m_lambdaShift;
      currentScore2 = m_table2->getScore(m2);

      mOffset = m_getSliceOffset(m1);

      for (int s1 = 0; s1 < scoreDim1; s1++)
      {

        for (int s2 = 0; s2 < scoreDim2; s2++)
        {
          currentPathScore2 = s2 + m_minScores2[m1].second;

          if (currentPathScore2 - currentScore2 < 0)
          {
            continue;
          }

          float incProb1 = 0, incProb2 = 0;

          for (int aa = 0; aa < aaMasses.size(); aa++)
          {
            prevMass2 = m2 - aaMasses[aa];

            // ignore invalid mass/score combinations
            if (prevMass2 < 0 || prevMass2 < end2)
            {
              continue;
            }

            prevMass1 = m1 - aaMasses[aa];

            prevOffset = m_getSliceOffset(prevMass1);

            prevScore2 = m_table2->getPreviousScore(m2, currentPathScore2, aa)
                - m_minScores2[prevMass1].second;

            if (prevScore2 < 0 || prevScore2 >= m_sizes2[prevMass1].second)
            {
              continue;
            }

            GFTable::incrementProbability(incProb2,
                                          m_table[prevOffset][s1][prevScore2].probability2,
                                          aaProb);
          }

          if (incProb2 > 0)
          {
            m_table[mOffset][s1][s2].probability2 = incProb2;
          }
        }
      }

      if (m1 == m_endMass2)
      {
        m_endDist2 = m_table[mOffset];
      }
    }
  }

  void GFTableIntersection::clear()
  {
    m_table1 = 0;
    m_table2 = 0;
    m_lambdaShift = 0;
    m_begin = 0;
    m_end = 0;
    m_aaJumps = 0;
  }

  void GFTableIntersection::compress()
  {
    m_table.resize(0);
  }

  double GFTableIntersection::getIntersectingPValue1(const int score1,
                                                     const int score2) const
  {
    if (m_endDist1.size() == 0 || m_endMass1 <= 0)
    {
      ERROR_MSG("Must compute generating function before distribution can be queried");
      abort();
    }
    double pvalue1 = 0;

    for (int s1 = score1 - m_minScores1[m_endMass1].first;
        s1 < m_sizes1[m_endMass1].first; s1++)
    {
      for (int s2 = score2 - m_minScores1[m_endMass1].second;
          s2 < m_sizes1[m_endMass1].second; s2++)
      {
        if (DEBUG)
        {
          //DEBUG_MSG("m = " << m_endMass1 << ", s1 = " << s1 << ", s2 = " << s2);
          //DEBUG_VAR(GFTable::getProbability(m_endDist1[s1][s2].probability1));
        }
        pvalue1 += GFTable::getProbability(m_endDist1[s1][s2].probability1);
      }
    }

    return pvalue1;
  }

  double GFTableIntersection::getIntersectingPValue2(const int score1,
                                                     const int score2) const
  {
    if (m_endDist2.size() == 0 || m_endMass2 <= 0)
    {
      ERROR_MSG("Must compute generating function before distribution can be queried");
      abort();
    }
    double pvalue2 = 0;

    if (DEBUG)
    {
      DEBUG_VAR(m_endMass2);
      DEBUG_VAR(m_minScores2[m_endMass2].first);
      DEBUG_VAR(m_minScores2[m_endMass2].second);
      DEBUG_VAR(m_sizes2[m_endMass2].first);
      DEBUG_VAR(m_sizes2[m_endMass2].second);
      DEBUG_VAR(score1);
      DEBUG_VAR(score2);
    }

    for (int s1 = score1 - m_minScores2[m_endMass2].first;
        s1 < m_sizes2[m_endMass2].first; s1++)
    {
      for (int s2 = score2 - m_minScores2[m_endMass2].second;
          s2 < m_sizes2[m_endMass2].second; s2++)
      {
        if (DEBUG)
        {
          DEBUG_MSG("m = " << m_endMass2 << ", s1 = " << s1 << ", s2 = " << s2);
          DEBUG_VAR(GFTable::getProbability(m_endDist2[s1][s2].probability2));
        }
        pvalue2 += GFTable::getProbability(m_endDist2[s1][s2].probability2);
      }
    }

    return pvalue2;
  }

  pair<double, double> GFTableIntersection::getIntersectingPValue(const string &peptide1,
                                                                  const string &peptide2) const
  {
    int cutoffScore1, cutoffScore2;

    cutoffScore1 = m_table1->lookupPeptide(0, 0, peptide1).second;

    if (DEBUG)
    {
      DEBUG_VAR(cutoffScore1);
    }

    if (peptide2 == peptide1 && m_lambdaShift == 0)
    {
      cutoffScore2 = m_table2->lookupPeptide(0, 0, peptide1).second;

      if (DEBUG)
      {
        DEBUG_VAR(cutoffScore2);
      }

      double pval1 = getIntersectingPValue1(cutoffScore1, cutoffScore2);

      return pair<double, double>(pval1, pval1);
    }

    string peptide2Ovlp = peptide2;
    string peptide1Ovlp = peptide1;

    vector<float> masses;
    m_aaJumps->getRoundedPRMMasses(peptide1, masses, 0);

    const int pm1 = floatToInt(masses[masses.size() - 1]);

    m_aaJumps->getRoundedPRMMasses(peptide2, masses, 0);

    const int pm2 = floatToInt(masses[masses.size() - 1]);

    if (pm1 < pm2 + m_lambdaShift)
    {
      peptide2Ovlp = m_aaJumps->getPeptidePrefix(peptide2,
                                                 (pm2 + m_lambdaShift) - pm1);
    }

    if (pm1 > m_lambdaShift + pm2)
    {
      peptide1Ovlp = m_aaJumps->getPeptidePrefix(peptide1,
                                                 pm1 - (pm2 + m_lambdaShift));
    }

    if (m_lambdaShift > 0)
    {
      peptide1Ovlp = m_aaJumps->getPeptideSuffix(peptide1Ovlp, m_lambdaShift);
    }

    if (DEBUG)
    {
      DEBUG_VAR(peptide1Ovlp);
      DEBUG_VAR(peptide2Ovlp);
      DEBUG_VAR(pm1);
      DEBUG_VAR(pm2);
      DEBUG_VAR(m_lambdaShift);
      DEBUG_VAR(m_startScore1);
    }

    cutoffScore2 = m_table2->lookupPeptide(0, 0, peptide2Ovlp).second;

    if (DEBUG)
    {
      DEBUG_VAR(cutoffScore2);
    }

    double pval1 = getIntersectingPValue1(cutoffScore1, cutoffScore2);

    cutoffScore1 = m_table1->lookupPeptide(m_lambdaShift,
                                           m_startScore1,
                                           peptide1Ovlp).second;

    if (DEBUG)
    {
      DEBUG_VAR(cutoffScore1);
    }

    cutoffScore2 = m_table2->lookupPeptide(0, 0, peptide2).second;

    if (DEBUG)
    {
      DEBUG_VAR(cutoffScore2);
    }

    double pval2 = getIntersectingPValue2(cutoffScore1, cutoffScore2);

    return pair<double, double>(pval1, pval2);
  }

  void GFTableIntersection::m_setDimensions()
  {
    const vector<int> &sizes1 = m_table1->m_sizes;
    const vector<int> &sizes2 = m_table2->m_sizes;

    const vector<int> &minScores1 = m_table1->m_minScores;
    const vector<int> &minScores2 = m_table2->m_minScores;

    int scoreDim1, scoreDim2, minScore1, minScore2;

    for (int m = 0; m < m_sizes1.size(); m++)
    {
      if (m <= m_begin)
      {
        scoreDim1 = sizes1[m];
        minScore1 = minScores1[m];
        scoreDim2 = 1;
        minScore2 = 0;
      }
      else if (m > m_begin && m <= m_end)
      {
        scoreDim1 = sizes1[m];
        minScore1 = minScores1[m];
        scoreDim2 = sizes2[m - m_lambdaShift];
        minScore2 = minScores2[m - m_lambdaShift];
      }
      else if (m > m_end && m < sizes1.size())
      {
        scoreDim1 = sizes1[m];
        minScore1 = minScores1[m];
        scoreDim2 = sizes2[sizes2.size() - 1];
        minScore2 = minScores2[minScores2.size() - 1];
      }
      else if (m >= sizes1.size() && m - m_lambdaShift < sizes2.size())
      {
        scoreDim1 = sizes1[sizes1.size() - 1];
        minScore1 = minScores1[minScores1.size() - 1];

        scoreDim2 = sizes2[m - m_lambdaShift];
        minScore2 = minScores2[m - m_lambdaShift];
      }
      else
      {
        ERROR_MSG("Invalid mass " << m);
        abort();
      }

      if (scoreDim1 == 0 || scoreDim2 == 0)
      {
        scoreDim1 = 0;
        scoreDim2 = 0;
        minScore1 = -1;
        minScore2 = -1;
      }

      m_sizes1[m].first = scoreDim1;
      m_sizes1[m].second = scoreDim2;

      if (DEBUG)
      {
        DEBUG_MSG("m = " << m << " has dims " << scoreDim1 << " x " << scoreDim2 << " (min_scores = " << minScore1 << " , " << minScore2 << ")");
      }

      m_minScores1[m].first = minScore1;
      m_minScores1[m].second = minScore2;
    }

    if (m_lambdaShift == 0)
    {
      m_sizes2 = m_sizes1;
      m_minScores2 = m_minScores1;
    }
    else
    {
      const vector<int> &aaMasses = m_table1->m_aaMasses;

      vector<int> newMaxScores;
      vector<int> newMinScores;
      for (int m = 0; m < m_sizes1.size(); m++)
      {
        if (m <= m_begin)
        {
          scoreDim1 = 1;
          minScore1 = 0;
          scoreDim2 = 1;
          minScore2 = 0;

          if (m == m_begin)
          {
            newMinScores.assign(sizes1.size(), -1);

            // all peptides matching the second spectrum must start here, need to re-compute minimum path
            for (int m1 = sizes1.size() - 1; m1 >= 0; m1--)
            {
              if (m1 >= m_table1->m_firstPM)
              {
                newMinScores[m1] = endScore1Ovlp;
              }
              else if (newMinScores[m1] < 0)
              {
                continue;
              }

              for (int aa = 0; aa < aaMasses.size(); aa++)
              {
                int prevScore = max(0,
                                    m_table1->getPreviousScore(m1,
                                                               newMinScores[m1],
                                                               aa));

                int prevMass = m1 - aaMasses[aa];

                if (prevMass < 0)
                {
                  continue;
                }

                if (newMinScores[prevMass] < 0)
                {
                  newMinScores[prevMass] = prevScore;
                }
                else
                {
                  newMinScores[prevMass] = min(newMinScores[prevMass],
                                               prevScore);
                }
              }
              if (DEBUG)
              {
                DEBUG_MSG("Min score at " << m1 << " is " << newMinScores[m]);
              }
            }

            newMaxScores.assign(sizes1.size(), -1);
            newMaxScores[m] = 0;

            // all peptides matching the second spectrum must start here, need to re-compute optimum path
            for (int m1 = m; m1 < sizes1.size(); m1++)
            {
              int locMinScore1 = newMinScores[m1];

              if (newMaxScores[m1] < 0 || locMinScore1 < 0)
              {
                continue;
              }

              if (locMinScore1 > newMaxScores[m1])
              {
                newMaxScores[m1] = -1;
                continue;
              }

              for (int aa = 0; aa < aaMasses.size(); aa++)
              {
                int nextMass = m1 + aaMasses[aa];
                if (nextMass >= sizes1.size())
                {
                  continue;
                }
                int pathScore = m_table1->getNextScore(m1,
                                                       newMaxScores[m1],
                                                       aa);

                newMaxScores[nextMass] = max(newMaxScores[nextMass], pathScore);
              }
            }
          }
        }
        else if (m > m_begin && m <= m_end)
        {
          int locMinScore1 = newMinScores[m];
          if (newMaxScores[m] < 0 || locMinScore1 < 0)
          {
            scoreDim1 = 0;
            minScore1 = -1;
          }
          else
          {
            scoreDim1 = newMaxScores[m] + 1 - locMinScore1;
            minScore1 = locMinScore1;
          }
          scoreDim2 = sizes2[m - m_lambdaShift];
          minScore2 = minScores2[m - m_lambdaShift];
        }
        else if (m > m_end && m < sizes1.size())
        {
          int locMinScore1 = newMinScores[m];
          if (newMaxScores[m] < 0 || locMinScore1 < 0)
          {
            scoreDim1 = 0;
            minScore1 = -1;
          }
          else
          {
            scoreDim1 = newMaxScores[m] + 1 - locMinScore1;
            minScore1 = locMinScore1;
          }
          scoreDim2 = sizes2[sizes2.size() - 1];
          minScore2 = minScores2[minScores2.size() - 1];
        }
        else if (m >= sizes1.size() && m - m_lambdaShift < sizes2.size())
        {
          int locMinScore1 = newMinScores[minScores1.size() - 1];

          scoreDim1 = newMaxScores[sizes1.size() - 1] + 1 - locMinScore1;
          minScore1 = locMinScore1;

          scoreDim2 = sizes2[m - m_lambdaShift];
          minScore2 = minScores2[m - m_lambdaShift];
        }
        else
        {
          ERROR_MSG("Invalid mass " << m);
          abort();
        }

        if (scoreDim1 == 0 || scoreDim2 == 0)
        {
          scoreDim1 = 0;
          scoreDim2 = 0;
          minScore1 = -1;
          minScore2 = -1;
        }

        m_sizes2[m].first = scoreDim1;
        m_sizes2[m].second = scoreDim2;

        if (DEBUG)
        {
          DEBUG_MSG("m = " << m << " has dims " << scoreDim1 << " x " << scoreDim2 << " (min_scores = " << minScore1 << " , " << minScore2 << ")");
        }

        m_minScores2[m].first = minScore1;
        m_minScores2[m].second = minScore2;
      }
    }
  }

  void GFTableIntersection::m_initializeSlice(const int startMass)
  {
    m_sliceSz = 0;
    m_sliceStartPos = startMass;
    m_maxAAMass = 0;

    for (int aa = 0; aa < m_table1->m_aaMasses.size(); aa++)
    {
      m_maxAAMass = max(m_maxAAMass, m_table1->m_aaMasses[aa]);
    }

    m_sliceSz = m_maxAAMass + TABLE_SLICE_BUFFER + 1;

    m_table.resize(m_sliceSz);

    int scoreDim1, scoreDim2;
    GFIntCell defaultCell;

    for (int m = startMass; m < m_sliceSz && m < m_sizes1.size(); m++)
    {
      scoreDim1 = max(m_sizes1[m].first, m_sizes2[m].first);
      scoreDim2 = max(m_sizes1[m].second, m_sizes2[m].second);

      if (scoreDim1 > m_table[m].size())
      {
        m_table[m].resize(scoreDim1);
      }

      for (int s1 = 0; s1 < scoreDim1; s1++)
      {
        if (scoreDim2 > m_table[m][s1].size())
        {
          m_table[m][s1].resize(scoreDim2);
        }

        for (int s2 = 0; s2 < scoreDim2; s2++)
        {
          m_table[m][s1][s2] = defaultCell;
        }
      }
    }
  }

  void GFTableIntersection::m_advanceSlice()
  {
    for (int m = 0; m <= m_maxAAMass; m++)
    {
      m_table[m].swap(m_table[m + TABLE_SLICE_BUFFER]);
    }

    m_sliceStartPos += TABLE_SLICE_BUFFER;

    int scoreDim1, scoreDim2;
    GFIntCell defaultCell;

    for (int m = m_sliceStartPos + m_maxAAMass + 1;
        (m - m_sliceStartPos) < m_sliceSz && m < m_sizes1.size(); m++)
    {
      int mOffset = m - m_sliceStartPos;

      scoreDim1 = max(m_sizes1[m].first, m_sizes2[m].first);
      scoreDim2 = max(m_sizes1[m].second, m_sizes2[m].second);

      if (scoreDim1 > m_table[mOffset].size())
      {
        m_table[mOffset].resize(scoreDim1);
      }

      for (int s1 = 0; s1 < scoreDim1; s1++)
      {
        if (scoreDim2 > m_table[mOffset][s1].size())
        {
          m_table[mOffset][s1].resize(scoreDim2);
        }

        for (int s2 = 0; s2 < scoreDim2; s2++)
        {
          m_table[mOffset][s1][s2] = defaultCell;
        }
      }
    }
  }
}
