/*
 * AssembledPeak.h
 *
 *  Created on: Apr 26, 2012
 *      Author: aguthals
 */

#ifndef ASSEMBLEDPEAK_H_
#define ASSEMBLEDPEAK_H_

#include "mzrange.h"
#include "spectrum.h"

using namespace specnets;

namespace abruijn
{
  class AssembledPeak : public specnets::MZRange
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    enum Symmetry
    {
      Symmetry_unknown = 0,
      Symmetry_B = 1,
      Symmetry_Y = 2
    };

    AssembledPeak();

    AssembledPeak(const AssembledPeak& other);

    AssembledPeak(const specnets::Spectrum& refSpec, unsigned int peakIdx);

    void clear();

    AssembledPeak& operator=(const AssembledPeak& other);

    void initialize(const specnets::Spectrum& refSpec, unsigned int peakIdx);

    void initializeEmpty(const string& refSpecID, float peakMass);

    virtual void addBinaryVersionInfo(map<string, unsigned short>& versions) const
    {
      versions[BIN_VERSION_ID] = BIN_VERSION;
      versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;
    }

    virtual bool saveToBinaryStream(FILE* fp) const;

    virtual bool loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions);

    inline const string& getSpecID() const
    {
      return m_specID;
    }

    inline const bool& isEndPt() const
    {
      return m_endPt;
    }

    inline void setEndPt(const bool& endPt)
    {
      m_endPt = endPt;
    }

    inline void setSymmetry(const Symmetry& newSymm)
    {
      m_BYsymmetry = newSymm;
    }

    inline const Symmetry& getSymmetry()
    {
      return m_BYsymmetry;
    }

  protected:

    string m_specID;
    bool m_endPt;
    Symmetry m_BYsymmetry;

  };

}

#endif /* ASSEMBLEDPEAK_H_ */
