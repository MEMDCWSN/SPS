/*
 * SpectrumAlignmentSet.cpp
 *
 *  Created on: May 1, 2012
 *      Author: aguthals
 */

#include "SpectrumAlignmentSet.h"

using namespace std;
using namespace specnets;

namespace abruijn
{

  const unsigned short SpectrumAlignmentSet::BIN_VERSION = 1;
  const unsigned short SpectrumAlignmentSet::BIN_SUBVERSION = 1;

  const string SpectrumAlignmentSet::BIN_VERSION_ID =
      "SpectrumAlignmentSet_binVersion";
  const string SpectrumAlignmentSet::BIN_SUBVERSION_ID =
      "SpectrumAlignmentSet_binSubVersion";

  SpectrumAlignmentSet::SpectrumAlignmentSet(unsigned int sz)
  {
    m_alignments.resize(sz);
  }

  SpectrumAlignmentSet::SpectrumAlignmentSet(const SpectrumAlignmentSet& other)
  {
    this->operator =(other);
  }

  SpectrumAlignmentSet &SpectrumAlignmentSet::operator=(const SpectrumAlignmentSet &other)
  {
    m_alignments = other.m_alignments;
    return *this;
  }

  void SpectrumAlignmentSet::appendAlignmentSet(const SpectrumAlignmentSet& other)
  {

    unsigned int idxUse = size();
    resize(size() + other.size());

    for (unsigned int i = 0; i < other.size(); i++)
    {
      m_alignments[idxUse++] = other[i];
    }
  }

  void SpectrumAlignmentSet::push_back(const SpectrumAlignment& other)
  {
    m_alignments.push_back(other);
  }

  bool SpectrumAlignmentSet::saveToBinaryStream(FILE* fp) const
  {
    if (fp == 0)
    {
      return false;
    }

    unsigned int numAligns = size();
    unsigned int count = fwrite(&numAligns, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write numAligns");
      return false;
    }

    for (unsigned int i = 0; i < numAligns; i++)
    {
      if (!m_alignments[i].saveToBinaryStream(fp))
      {
        ERROR_MSG("Failed to save SpectrumAlignment " << i);
        return false;
      }
    }

    return true;
  }

  bool SpectrumAlignmentSet::loadFromBinaryStream(FILE* fp, map<string,
      unsigned short>& versions)
  {
    if (fp == 0)
    {
      return false;
    }

    unsigned short version = versions[BIN_VERSION_ID];
    unsigned short subVersion = versions[BIN_SUBVERSION_ID];

    unsigned int numAligns;
    unsigned int count = fread(&numAligns, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not fread numAligns");
      return false;
    }
    resize(numAligns);

    for (unsigned int i = 0; i < numAligns; i++)
    {
      if (!m_alignments[i].loadFromBinaryStream(fp, versions))
      {
        ERROR_MSG("Failed to load SpectrumAlignment " << i);
        return false;
      }
    }

    return true;
  }
}
