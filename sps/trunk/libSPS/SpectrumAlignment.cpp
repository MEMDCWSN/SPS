/*
 * SpectrumAligner.cpp
 *
 *  Created on: Apr 30, 2012
 *      Author: aguthals
 */

#include "SpectrumAlignment.h"

using namespace std;
using namespace specnets;

namespace abruijn
{

  const unsigned short SpectrumAlignment::BIN_VERSION = 1;
  const unsigned short SpectrumAlignment::BIN_SUBVERSION = 1;

  const string SpectrumAlignment::BIN_VERSION_ID =
      "SpectrumAlignment_binVersion";
  const string SpectrumAlignment::BIN_SUBVERSION_ID =
      "SpectrumAlignment_binSubVersion";

  SpectrumAlignment::SpectrumAlignment()
  {
    m_matchedPeaksB.resize(0);
    m_matchedPeaksY.resize(0);
    m_matchSuffixMasses = true;
    m_spec1ID = "";
    m_spec2ID = "";
    m_score1 = 0;
    m_score2 = 0;
    m_prob1 = 1.0;
    m_prob2 = 1.0;
    m_shiftB.set(0, 0, 0);
    m_shiftY.set(0, 0, 0);
    m_revSpec2 = false;
  }

  SpectrumAlignment::SpectrumAlignment(const SpectrumAlignment& other)
  {
    this->operator =(other);
  }

  SpectrumAlignment & SpectrumAlignment::operator=(const SpectrumAlignment &other)
  {
    m_matchedPeaksB = other.m_matchedPeaksB;
    m_matchedPeaksY = other.m_matchedPeaksY;
    m_matchSuffixMasses = other.m_matchSuffixMasses;
    m_spec1ID = other.m_spec1ID;
    m_spec2ID = other.m_spec2ID;
    m_score1 = other.m_score1;
    m_score2 = other.m_score2;
    m_prob1 = other.m_prob1;
    m_prob2 = other.m_prob2;
    m_shiftB = other.m_shiftB;
    m_shiftY = other.m_shiftY;
    m_revSpec2 = other.m_revSpec2;
    return *this;
  }

  MZRange SpectrumAlignment::getShift(const string& spec1ID,
                                      const string& spec2ID) const
  {
    if (spec1ID == m_spec1ID && spec2ID == m_spec2ID)
    {
      return MZRange(m_shiftB);
    }
    else if (spec1ID == m_spec2ID && spec2ID == m_spec1ID)
    {
      return MZRange(0 - m_shiftB.getMass(),
                     m_shiftB.getIntensity(),
                     m_shiftB.getTolerance());
    }
    else
    {
      ERROR_MSG("Found unknown IDs \'" << spec1ID << "\' and \'" << spec2ID << "\'");
      abort();
    }
  }

  bool SpectrumAlignment::saveToBinaryStream(FILE* fp) const
  {
    unsigned int count;

    count = fwrite(&m_matchSuffixMasses, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_matchSuffixMasses");
      return false;
    }

    count = fwrite(&m_revSpec2, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_revSpec2");
      return false;
    }

    count = fwrite(&m_score1, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_score1");
      return false;
    }

    count = fwrite(&m_score2, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_score2");
      return false;
    }

    count = fwrite(&m_prob1, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_prob1");
      return false;
    }

    count = fwrite(&m_prob2, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_prob2");
      return false;
    }

    double mass = m_shiftB.getMass();
    double intensity = m_shiftB.getIntensity();
    double tolerance = m_shiftB.getTolerance();

    count = fwrite(&mass, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the m_shiftB mass");
      return false;
    }

    count = fwrite(&intensity, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the m_shiftB intensity");
      return false;
    }

    count = fwrite(&tolerance, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the m_shiftB tolerance");
      return false;
    }

    mass = m_shiftY.getMass();
    intensity = m_shiftY.getIntensity();
    tolerance = m_shiftY.getTolerance();

    count = fwrite(&mass, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the m_shiftY mass");
      return false;
    }

    count = fwrite(&intensity, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the m_shiftY intensity");
      return false;
    }

    count = fwrite(&tolerance, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the m_shiftY tolerance");
      return false;
    }

    vector<string> labels(2);
    labels[0] = m_spec1ID;
    labels[1] = m_spec2ID;

    if (!writeStringsToBinaryStream(fp, labels))
    {
      ERROR_MSG("Could not write the spectrum IDs");
      return false;
    }

    unsigned int numMatchedB = m_matchedPeaksB.size();
    count = fwrite(&numMatchedB, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write numMatchedB");
      return false;
    }

    if (numMatchedB > 0)
    {
      unsigned int bufSize = numMatchedB * 6;
      double* peakBuf = (double*)malloc(sizeof(double) * bufSize);
      unsigned int idxUse = 0;
      for (list<pair<MZRange, MZRange> >::const_iterator mIt =
          m_matchedPeaksB.begin(); mIt != m_matchedPeaksB.end(); mIt++)
      {
        peakBuf[idxUse++] = mIt->first.getMass();
        peakBuf[idxUse++] = mIt->first.getIntensity();
        peakBuf[idxUse++] = mIt->first.getTolerance();
        peakBuf[idxUse++] = mIt->second.getMass();
        peakBuf[idxUse++] = mIt->second.getIntensity();
        peakBuf[idxUse++] = mIt->second.getTolerance();
      }
      count = fwrite(peakBuf, sizeof(double), bufSize, fp);
      if (count == 0)
      {
        ERROR_MSG("Could not write m_matchedPeaksB");
        return false;
      }
      free(peakBuf);
    }

    unsigned int numMatchedY = m_matchedPeaksY.size();
    count = fwrite(&numMatchedY, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write numMatchedY");
      return false;
    }

    if (numMatchedY > 0)
    {
      unsigned int bufSize = numMatchedY * 6;
      double* peakBuf = (double*)malloc(sizeof(double) * bufSize);
      unsigned int idxUse = 0;
      for (list<pair<MZRange, MZRange> >::const_iterator mIt =
          m_matchedPeaksY.begin(); mIt != m_matchedPeaksY.end(); mIt++)
      {
        peakBuf[idxUse++] = mIt->first.getMass();
        peakBuf[idxUse++] = mIt->first.getIntensity();
        peakBuf[idxUse++] = mIt->first.getTolerance();
        peakBuf[idxUse++] = mIt->second.getMass();
        peakBuf[idxUse++] = mIt->second.getIntensity();
        peakBuf[idxUse++] = mIt->second.getTolerance();
      }
      count = fwrite(peakBuf, sizeof(double), bufSize, fp);
      if (count == 0)
      {
        ERROR_MSG("Could not write m_matchedPeaksY");
        return false;
      }
      free(peakBuf);
    }

    return true;
  }

  bool SpectrumAlignment::loadFromBinaryStream(FILE* fp, map<string,
      unsigned short>& versions)
  {
    unsigned short version = versions[BIN_VERSION_ID];
    unsigned short subVersion = versions[BIN_SUBVERSION_ID];

    unsigned int count;

    count = fread(&m_matchSuffixMasses, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_matchSuffixMasses");
      return false;
    }

    count = fread(&m_revSpec2, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_revSpec2");
      return false;
    }

    count = fread(&m_score1, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_score1");
      return false;
    }

    count = fread(&m_score2, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_score2");
      return false;
    }

    count = fread(&m_prob1, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_prob1");
      return false;
    }

    count = fread(&m_prob2, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_prob2");
      return false;
    }

    double mass = 0;
    double intensity = 0;
    double tolerance = 0;

    count = fread(&mass, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the m_shiftB mass");
      return false;
    }

    count = fread(&intensity, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the m_shiftB intensity");
      return false;
    }

    count = fread(&tolerance, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the m_shiftB tolerance");
      return false;
    }
    m_shiftB.set(mass, intensity, tolerance);

    count = fread(&mass, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the m_shiftY mass");
      return false;
    }

    count = fread(&intensity, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the m_shiftY intensity");
      return false;
    }

    count = fread(&tolerance, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the m_shiftY tolerance");
      return false;
    }
    m_shiftY.set(mass, intensity, tolerance);

    vector<string> labels;
    if (!readStringsFromBinaryStream(fp, labels) || labels.size() != 2)
    {
      ERROR_MSG("Could not read the spectrum IDs");
      return false;
    }
    m_spec1ID = labels[0];
    m_spec2ID = labels[1];

    unsigned int numMatchedB = 0;
    count = fread(&numMatchedB, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read numMatchedB");
      return false;
    }
    m_matchedPeaksB.clear();
    pair<MZRange, MZRange> nextPair;

    if (numMatchedB > 0)
    {
      unsigned int bufSize = numMatchedB * 6;
      double* peakBuf = (double*)malloc(sizeof(double) * bufSize);
      count = fread(peakBuf, sizeof(double), bufSize, fp);
      if (count == 0)
      {
        ERROR_MSG("Could not read m_matchedPeaksB");
        return false;
      }
      unsigned int idxUse = 0;
      for (unsigned int i = 0; i < numMatchedB; i++)
      {
        nextPair.first.setMass(peakBuf[idxUse++]);
        nextPair.first.setIntensity(peakBuf[idxUse++]);
        nextPair.first.setTolerance(peakBuf[idxUse++]);
        nextPair.second.setMass(peakBuf[idxUse++]);
        nextPair.second.setIntensity(peakBuf[idxUse++]);
        nextPair.second.setTolerance(peakBuf[idxUse++]);
        m_matchedPeaksB.push_back(nextPair);
      }
      free(peakBuf);
    }

    unsigned int numMatchedY = 0;
    count = fread(&numMatchedY, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read numMatchedY");
      return false;
    }
    m_matchedPeaksY.clear();

    if (numMatchedY > 0)
    {
      unsigned int bufSize = numMatchedY * 6;
      double* peakBuf = (double*)malloc(sizeof(double) * bufSize);
      count = fread(peakBuf, sizeof(double), bufSize, fp);
      if (count == 0)
      {
        ERROR_MSG("Could not read m_matchedPeaksY");
        return false;
      }
      unsigned int idxUse = 0;
      for (unsigned int i = 0; i < numMatchedY; i++)
      {
        nextPair.first.setMass(peakBuf[idxUse++]);
        nextPair.first.setIntensity(peakBuf[idxUse++]);
        nextPair.first.setTolerance(peakBuf[idxUse++]);
        nextPair.second.setMass(peakBuf[idxUse++]);
        nextPair.second.setIntensity(peakBuf[idxUse++]);
        nextPair.second.setTolerance(peakBuf[idxUse++]);
        m_matchedPeaksY.push_back(nextPair);
      }
      free(peakBuf);
    }

    return true;
  }

}
