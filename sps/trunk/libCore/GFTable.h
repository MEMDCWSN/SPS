/*
 * GFTable.h
 *
 *  Created on: Jul 19, 2013
 *      Author: aguthals
 */

#ifndef GFTABLE2_H_
#define GFTABLE2_H_

#include "Cluster.h"
#include "spectrum.h"
#include "PeptideSpectrumMatch.h"
#include "GFTableIntersection.h"
#include "MassErrorModel.h"
//#include "BreakScoreSpectrum.h"

#include <vector>
#include <map>
#include <set>
#include <list>
#include <math.h>

using namespace std;

namespace specnets
{
  class GFTableIntersection;

  struct GFProfileCell
  {
  public:
    // Cummulative probability of all peptides with a prefix at this node
    float probability;
    unsigned long long prefixes;
    unsigned long long suffixes;

    int rEdges;

    GFProfileCell() :
        probability(0), prefixes(0), suffixes(0), rEdges(0)
    {
    }

    GFProfileCell &operator=(const GFProfileCell &other)
    {
      if (this == &other)
      {
        return *this;
      }
      probability = other.probability;
      prefixes = other.prefixes;
      suffixes = other.suffixes;
      rEdges = other.rEdges;
      //zeroProb = other.zeroProb;
      //prefixes = other.prefixes;
      return *this;
    }
  };

  struct GFCell
  {
  public:

    float probability;
    int rEdges;

    GFCell() :
        probability(0), rEdges(0)
    {

    }

    GFCell &operator=(const GFCell &other)
    {
      if (this == &other)
      {
        return *this;
      }
      probability = other.probability;
      rEdges = other.rEdges;
      //zeroProb = other.zeroProb;
      //prefixes = other.prefixes;
      return *this;
    }
  };

  class GFTable : public SpectrumItem
  {

  protected:

    // ith mass -> jth score -> GFCell
    vector<vector<GFCell> > m_table;

    // ith mass -> jth AA -> backwards edge score
    vector<vector<int> > m_edgeScores;

    // ith mass -> jth score -> GFProfileCell
    vector<vector<GFProfileCell> > m_tableProfile;

    vector<int> m_nodeScores;

    // Lower mass range of parent mass
    int m_firstPM;

    int m_startAlign;
    int m_endAlign;

    set<int> m_endMasses;

    // Integer-valued amino acid masses loaded from AAJumps
    vector<int> m_aaMasses;

    vector<int> m_sizes;
    vector<int> m_minScores;
    vector<int> m_maxScores;

    // Pointer to AAJumps that were used to initialize the table
    const AAJumps *m_aaJumps;

    bool m_useProfile;

  public:

    static bool DEBUG;

    static const double PROBABILITY_FACTOR;

    static const float EDGE_SCORE_SCALING_FACTOR;

    inline static void setEdge(int &rEdges, const int &index, const bool &value)
    {
      int mask = 1;
      if (value)
      {
        /*
         DEBUG_VAR(cell.fEdges);
         DEBUG_VAR(index);
         int val2 = mask << index;
         DEBUG_VAR(val2);
         int val3 = cell.fEdges | val2;
         DEBUG_VAR(val3);
         */
        rEdges = rEdges | (mask << index);
      }
      else
      {
        rEdges = rEdges & (~(mask << index));
      }
    }

    inline static bool getEdge(const int &rEdges, const int &index)
    {
      int mask = 1;
      return ((mask << index) & rEdges) != 0;
    }

    inline static double getProbability(const float &cell)
    {
      return ((double)cell) / PROBABILITY_FACTOR;
    }

    inline static void setProbability(float &cell, const double &prob)
    {
      cell = prob * PROBABILITY_FACTOR;
    }

    inline static void incrementProbability(float &celltoInc,
                                            const float &cellRead,
                                            const float& aaProb)
    {
      celltoInc += cellRead * aaProb;
    }

    // Default constructor
    GFTable();

    // Copy constructor
    GFTable(const GFTable &other);

    /** Constructs graph from a spectrum
     * @param spectrum
     * @param jumps set of AA masses to impose edges
     **/
    GFTable(const Spectrum &spectrum,
            const AAJumps &jumps,
            const Spectrum &msSpectrum,
            const MassErrorModel &massErrorModelTarget,
            const MassErrorModel &massErrorModelDecoy,
            const float maxPPMError);

    /** Constructs graph from a spectrum that already has integer-binned masses
     * @param spectrum BreakScoreSpectrum::toIntegerBins should have already been called on this along with some normalization
     * @param jumps set of AA masses to impose edges
     * @param useProfile if true, sacrifice cache performance to compute prefix/suffix counts for every cell (needed for profile construction)
     **/
    void initialize(const Spectrum &spectrum,
                    const AAJumps &jumps,
                    const Spectrum &msSpectrum,
                    const MassErrorModel &massErrorModelTarget,
                    const MassErrorModel &massErrorModelDecoy,
                    const float maxPPMError);

    void addPRMScores(const int shift, const Spectrum &spectrum);

    void addMatchingPRMScores(const int shift, const Spectrum &spectrum);

    void computeGeneratingFunction(const int minPathScore = 0,
                                   const int startMass = 0,
                                   const bool useProfile = false);

    void clear();

    GFTable &operator=(const GFTable &other);

    // Access the vector of GFNodes at a given mass
    inline const vector<GFCell> &operator[](int mass) const
    {
      return m_table[mass];
    }

    // Access the vector of GFNodes at a given mass
    inline vector<GFCell> &operator[](int mass)
    {
      return m_table[mass];
    }

    inline const vector<int> &getAAMasses() const
    {
      return m_aaMasses;
    }

    // Access the score at a given mass
    inline int getScore(const int mass) const
    {
      return m_nodeScores[mass];
    }

    // Get the number of mass bins encoded in the graph
    inline int size() const
    {
      return m_sizes.size();
    }

    inline int getSize(const int m) const
    {
      return m_sizes[m];
    }

    // Get the lower bound of the parent mass
    inline int getPmLowerBound() const
    {
      return m_firstPM;
    }

    inline int getMaxScore() const
    {
      return m_table[m_sizes.size() - 1].size();
    }

    inline const set<int>& getEndMasses() const
    {
      return m_endMasses;
    }

    inline const AAJumps* getAAJumps() const
    {
      return m_aaJumps;
    }

    int getScoreThreshold(const double probThreshold) const;

    /** Removes all nodes/edges (BY MARKING THEM AS INVALID) that do not encode peptides that have probability within the given threshold.
     * @param probThreshold
     * @return minimum score that yields the probability
     **/
    inline int encodeDictionary(const double probThreshold)
    {
      int cutoffScore = getScoreThreshold(probThreshold);
      return encodeDictionary(cutoffScore);
    }

    int encodeDictionary(const int cutoffScore);

    pair<int, int>
    lookupPeptide(const int startMass,
                  const int startScore,
                  const string &peptide) const;

    int getMatchScore(const string &peptide, int *startMass = 0, int *endMass =
                          0) const;

    inline double getPValue(const string &peptide) const
    {
      int score = getMatchScore(peptide);
      if (score == 0)
      {
        return 1.0;
      }
      return getPValue(score);
    }

    inline int getPreviousScore(const int &mass,
                                const int &score,
                                const int &prevAA) const
    {
      const int prevMass = mass - m_aaMasses[prevAA];
      return score - m_nodeScores[mass] - m_edgeScores[mass][prevAA];
    }

    inline int getNextScore(const int &mass,
                            const int &score,
                            const int &nextAA) const
    {
      const int nextMass = mass + m_aaMasses[nextAA];
      return score + m_nodeScores[nextMass] + m_edgeScores[nextMass][nextAA];
    }

    double getPValue(const int score) const;

    inline void computeProbabilities()
    {
      m_computeProbabilities();
    }

    void optimizeDimensions(const int minPathScore, const int startMass = 0);

  protected:

    friend class GFTableIntersection;

    /**
     * Backtracks m_table from a given minimum cutoff and fills in the suffix counts
     * @param cutoffScore
     */
    void m_exploreBackwards(const int cutoffScore);

    void m_computeProbabilities(const int startMass = 0);

    void m_computeProbabilitiesProfile(const int startMass = 0);

    void m_computeMinScores(const int minPathScore);

    void m_computeMaxScores(const int startMass = 0);
  };
}

#endif /* GFTABLE2_H_ */
