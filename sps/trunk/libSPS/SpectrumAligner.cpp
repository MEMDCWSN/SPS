/*
 * SpectrumAligner.cpp
 *
 *  Created on: Apr 30, 2012
 *      Author: aguthals
 */

#include "SpectrumAligner.h"

using namespace std;
using namespace specnets;

namespace abruijn
{

  SpectrumAligner::SpectrumAligner(bool matchSuffixMasses) :
    SpectrumAlignment()
  {
    m_matchSuffixMasses = matchSuffixMasses;
    m_spec1.resize(0);
    m_spec2.resize(0);
    m_spec2Rev.resize(0);
  }

  SpectrumAligner::SpectrumAligner(const SpectrumAligner& other)
  {
    this->operator =(other);
  }

  SpectrumAligner &SpectrumAligner::operator=(const SpectrumAligner &other)
  {
    m_spec1 = other.m_spec1;
    m_spec2 = other.m_spec2;
    m_spec2Rev = other.m_spec2Rev;
    SpectrumAlignment::operator=((SpectrumAlignment&)other);
    return *this;
  }

  void SpectrumAligner::setSpec1(const Spectrum& newSpec1)
  {
    m_spec1 = newSpec1;
  }

  void SpectrumAligner::setSpec2(const Spectrum& newSpec2)
  {
    m_spec2 = newSpec2;
  }

  void SpectrumAligner::setShift(const MZRange& shiftB)
  {
    m_shiftB = shiftB;
  }

  void SpectrumAligner::setShift(pair<MZRange, MZRange>& initialMatchedPeaks) {

  }

  void SpectrumAligner::computeBestShift(int maxNumMods)
  {

  }

  void SpectrumAligner::scoreOverlap(int maxNumMods)
  {

  }
}
