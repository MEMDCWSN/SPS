/*
 * SpectrumAligner.h
 *
 *  Created on: Apr 30, 2012
 *      Author: aguthals
 */

#ifndef SPECTRUMALIGNMENT_H_
#define SPECTRUMALIGNMENT_H_

#include "aminoacid.h"
#include "mzrange.h"
#include "spectrum.h"
#include "Logger.h"

#include <vector>

using namespace std;
using namespace specnets;

namespace abruijn
{

  class SpectrumAlignment
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    SpectrumAlignment();

    SpectrumAlignment(const SpectrumAlignment& other);

    SpectrumAlignment &operator=(const SpectrumAlignment &other);

    inline const string& getSpec1ID() const
    {
      return m_spec1ID;
    }

    inline const string& getSpec2ID() const
    {
      return m_spec2ID;
    }

    inline void setSpec1ID(const string& spec1ID)
    {
      m_spec1ID = spec1ID;
    }

    inline void setSpec2ID(const string& spec2ID)
    {
      m_spec2ID = spec2ID;
    }

    inline double getMinScore() const
    {
      return min(m_score1, m_score2);
    }

    inline unsigned int getNumMatchedPeaks() const
    {
      return m_matchedPeaksB.size();
    }

    inline void clearMP()
    {
      m_matchedPeaksB.clear();
      m_matchedPeaksY.clear();
    }

    inline void addMatchedPeaksB(const MZRange& spec1Pk, const MZRange& spec2Pk)
    {
      m_matchedPeaksB.push_back(pair<MZRange, MZRange> (spec1Pk, spec2Pk));
    }

    inline void outputMatchedPeaks(list<pair<MZRange, MZRange> >& outputMatchedPeaksB,
                                   list<pair<MZRange, MZRange> >& outputMatchedPeaksY) const
    {
      outputMatchedPeaksB = m_matchedPeaksB;
      outputMatchedPeaksY = m_matchedPeaksY;
    }

    MZRange getShift(const string& spec1ID, const string& spec2ID) const;

    virtual void addBinaryVersionInfo(map<string, unsigned short>& versions) const
    {
      versions[BIN_VERSION_ID] = BIN_VERSION;
      versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;
    }

    virtual bool saveToBinaryStream(FILE* fp) const;

    virtual bool loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions);

  protected:
    list<pair<MZRange, MZRange> > m_matchedPeaksB;
    list<pair<MZRange, MZRange> > m_matchedPeaksY;

    bool m_matchSuffixMasses;

    string m_spec1ID;
    string m_spec2ID;

    double m_score1;
    double m_score2;

    double m_prob1;
    double m_prob2;

    MZRange m_shiftB;
    MZRange m_shiftY;

    bool m_revSpec2;
  };
}

#endif /* SPECTRUMALIGNMENT_H_ */
