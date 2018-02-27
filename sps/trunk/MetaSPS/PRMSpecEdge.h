/*
 * PRMSpecEdge.h
 *
 *  Created on: Jun 23, 2011
 *      Author: aguthals
 */

#ifndef PRMSPECEDGE_H_
#define PRMSPECEDGE_H_

#include "SpectrumPair.h"
#include "Logger.h"
#include "spectrum.h"
#include "SpecSet.h"

namespace specnets
{
  class PRMSpecEdge : public SpectrumPair
  {
  public:

    static void appendEdges(PRMSpecEdge* edge1,
                            PRMSpecEdge* edge2,
                            PRMSpecEdge* outEdge1Plus2);

    int index;

    PRMSpecEdge(void) :
      SpectrumPair(), index(0)
    {
    }

    PRMSpecEdge(const PRMSpecEdge& copyFrom) :
      SpectrumPair(), index(0)
    {
      (*this) = copyFrom;
    }

    PRMSpecEdge(const int idx,
                const SpectrumPair& specPair,
                const SpecSet& spectra) :
      SpectrumPair(), index(0)
    {
      initializeEdge(idx, specPair, spectra);
    }

    float getMinScore() const
    {
      return min(score1, score2);
    }

    PRMSpecEdge& operator=(const PRMSpecEdge &other);

    PRMSpecEdge& operator=(const SpectrumPair &other);

    bool operator<(const PRMSpecEdge &other) const
    {
      if (getMinScore() < other.getMinScore())
        return true;
      return false;
    }

    bool operator>(const PRMSpecEdge &other) const
    {
      if (getMinScore() > other.getMinScore())
        return true;
      return false;
    }

    bool operator==(const PRMSpecEdge &other) const
    {
      if (getMinScore() == other.getMinScore())
        return true;
      return false;
    }

    void initializeEdge(const int idx,
                        const SpectrumPair& specPair,
                        const SpecSet& spectra);

    float getShift(int idx1, int idx2) const;

    float getReversedShift(int idx1, int idx2) const;

    void reverse(bool reverse1, bool reverse2);

    void reverse(int idx);

  };
}

#endif /* PRMSPECEDGE_H_ */
