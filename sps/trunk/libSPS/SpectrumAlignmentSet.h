/*
 * SpectrumAlignmentSet.h
 *
 *  Created on: May 1, 2012
 *      Author: aguthals
 */

#ifndef SPECTRUMALIGNMENTSET_H_
#define SPECTRUMALIGNMENTSET_H_

#include "SpectrumAlignment.h"

#include <vector>

using namespace std;
using namespace specnets;

namespace abruijn
{

  class SpectrumAlignment;

  class SpectrumAlignmentSet
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    SpectrumAlignmentSet(unsigned int sz = 0);

    SpectrumAlignmentSet(const SpectrumAlignmentSet& other);

    SpectrumAlignmentSet &operator=(const SpectrumAlignmentSet &other);

    abruijn::SpectrumAlignment &operator[](unsigned int i)
    {
      return m_alignments[i];
    }

    const abruijn::SpectrumAlignment &operator[](unsigned int i) const
    {
      return m_alignments[i];
    }

    unsigned int size() const
    {
      return m_alignments.size();
    }

    void resize(unsigned int newSize)
    {
      m_alignments.resize(newSize);
    }

    void appendAlignmentSet(const SpectrumAlignmentSet& other);

    void push_back(const abruijn::SpectrumAlignment& other);

    virtual void addBinaryVersionInfo(map<string, unsigned short>& versions) const
    {
      versions[BIN_VERSION_ID] = BIN_VERSION;
      versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;
      SpectrumAlignment align;
      align.addBinaryVersionInfo(versions);
    }

    virtual bool saveToBinaryStream(FILE* fp) const;

    virtual bool loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions);

  protected:
    vector<abruijn::SpectrumAlignment> m_alignments;
  };
}

#endif /* SPECTRUMALIGNMENTSET_H_ */
