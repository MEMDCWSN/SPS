/*
 * SpectrumAligner.h
 *
 *  Created on: Apr 30, 2012
 *      Author: aguthals
 */

#ifndef SPECTRUMALIGNER_H_
#define SPECTRUMALIGNER_H_

#include "SpectrumAlignment.h"
#include "Logger.h"

using namespace std;
using namespace specnets;

namespace abruijn
{

  class SpectrumAlignment;

  class SpectrumAligner : public SpectrumAlignment
  {
  public:

    SpectrumAligner(bool matchSuffixMasses = true);

    SpectrumAligner(const SpectrumAligner& other);

    SpectrumAligner &operator=(const SpectrumAligner &other);

    void setSpec1(const Spectrum& newSpec1);

    void setSpec2(const Spectrum& newSpec2);

    void setShift(const MZRange& shiftB);

    void setShift(pair<MZRange, MZRange>& initialMatchedPeaks);

    void computeBestShift(int maxNumMods = -1);

    void scoreOverlap(int maxNumMods = -1);

  protected:

    Spectrum m_spec1;
    Spectrum m_spec2;
    Spectrum m_spec2Rev;

  };

}

#endif /* SPECTRUMALIGNER_H_ */
