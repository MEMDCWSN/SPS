/*
 * GFTableIntersection.h
 *
 *  Created on: Sep 5, 2013
 *      Author: aguthals
 */

#ifndef GFTABLEINTERSECTION2_H_
#define GFTABLEINTERSECTION2_H_

#include <vector>

#include "GFTable.h"

using namespace std;

namespace specnets
{
  class GFTable;

  struct GFIntCell
  {
  public:

    float probability1;
    float probability2;

    GFIntCell() :
        probability1(0), probability2(0)
    {

    }

    GFIntCell &operator=(const GFIntCell &other)
    {
      if (this == &other)
      {
        return *this;
      }
      probability1 = other.probability1;
      probability2 = other.probability2;
      //zeroProb = other.zeroProb;
      //prefixes = other.prefixes;
      return *this;
    }
  };

  class GFTableIntersection
  {
  protected:
    // ith mass -> jth score -> GFCell
    vector<vector<vector<GFIntCell> > > m_table;

    vector<pair<int, int> > m_sizes1;
    vector<pair<int, int> > m_minScores1;

    vector<pair<int, int> > m_sizes2;
    vector<pair<int, int> > m_minScores2;

    const GFTable *m_table1;
    const GFTable *m_table2;

    int m_lambdaShift;
    int m_begin;
    int m_end;
    int m_startScore1;

    int endScore1Ovlp;

    const AAJumps *m_aaJumps;

    int m_sliceSz;
    int m_maxAAMass;
    int m_sliceStartPos;

    int m_endMass1;
    int m_endMass2;
    vector<vector<GFIntCell> > m_endDist1;
    vector<vector<GFIntCell> > m_endDist2;

  public:

    static bool DEBUG;

    static const int TABLE_SLICE_BUFFER;

    GFTableIntersection();

    GFTableIntersection(const GFTable &table1,
                        const GFTable &table2,
                        const int shift,
                        const int matchScoreOvlp);

    void initialize(const GFTable &table1,
                    const GFTable &table2,
                    const int shift,
                    const int matchScoreOvlp);

    void computeGeneratingFunction(const bool useProfile = false);

    void clear();

    void compress();

    double getIntersectingPValue1(const int score1, const int score2) const;

    double getIntersectingPValue2(const int score1, const int score2) const;

    pair<double, double> getIntersectingPValue(const string &peptide1,
                                               const string &peptide2) const;

  protected:

    void m_setDimensions();

    void m_initializeSlice(const int startMass = 0);

    inline int m_getSliceOffset(const int& m)
    {
      if (m < m_sliceStartPos)
      {
        ERROR_MSG("Mass " << m << " must be greater than offset " << m_sliceStartPos);
        abort();
      }
      if (m - m_sliceStartPos >= m_sliceSz)
      {
        m_advanceSlice();
      }
      return m - m_sliceStartPos;
    }

    void m_advanceSlice();
  };
}

#endif /* GFTABLEINTERSECTION_H_ */
