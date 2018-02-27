/*
 * SvmScaleParameterList.cpp
 *
 *  Created on: Jan 1, 2011
 *      Author: jsnedecor
 */

/**
 * Load in svm scaling parameters
 */

#include "SvmScaleParameterList.h"

namespace specnets
{
  //---------------------------------------------------------------------------
  SvmScaleParameterList::SvmScaleParameterList(void)
  {
    //initialize required header
    m_requiredHeader.push_back("Key");
    m_requiredHeader.push_back("StartRange");
    m_requiredHeader.push_back("EndRange");
    m_requiredHeader.push_back("Feature");
    m_requiredHeaderIndex.resize(m_requiredHeader.size());
  }

  //---------------------------------------------------------------------------
  SvmScaleParameterList::~SvmScaleParameterList(void)
  {
    // EMPTY
  }

  //---------------------------------------------------------------------------

  bool SvmScaleParameterList::loadSvmKeyRange(const char * keyFile)
  {
    if(!DelimitedTextReader::loadDelimitedFile(keyFile,"\t","#",m_header,m_lines,m_requiredHeader,m_requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open keyFile! " << keyFile);
      return false;
    }
    //fill key map
    for (unsigned int i = 0; i<m_lines.size(); i++)
    {
      unsigned int nameHeaderIndex = m_requiredHeaderIndex[3];
      unsigned int keyHeaderIndex = m_requiredHeaderIndex[0];
      unsigned int keyIndex;

      sscanf(m_lines[i][keyHeaderIndex].c_str(), "%u", &keyIndex);
      m_nameToKey[m_lines[i][nameHeaderIndex]] = keyIndex;
    }
    return true;
  }
  //---------------------------------------------------------------------------
  bool SvmScaleParameterList::getKeyIndexByName(string &name, unsigned int &outputIndex)
  {
    map<string, unsigned int>::const_iterator it;

    //DEBUG_VAR(m_nameToKey.size());
    it = m_nameToKey.find(name);

    if (it==m_nameToKey.end())
    {
      return false;
    }

    outputIndex = it->second;
    //DEBUG_VAR(name);
    //DEBUG_VAR(outputIndex);
    return true;
  }
  //---------------------------------------------------------------------------
  string SvmScaleParameterList::getKeyName(unsigned int lineIndex) const
  {
    unsigned int keyIndex = m_requiredHeaderIndex[3];

    return m_lines[lineIndex][keyIndex];
  }

  //---------------------------------------------------------------------------
  double SvmScaleParameterList::getStartRange(unsigned int lineIndex) const
  {
    unsigned int keyIndex = m_requiredHeaderIndex[1];

    double returnValue;
    sscanf(m_lines[lineIndex][keyIndex].c_str(), "%lf", &returnValue);
    return returnValue;
  }

  //---------------------------------------------------------------------------
  double SvmScaleParameterList::getEndRange(unsigned int lineIndex) const
  {
    unsigned int keyIndex = m_requiredHeaderIndex[2];

    double returnValue;
    sscanf(m_lines[lineIndex][keyIndex].c_str(), "%lf", &returnValue);
    return returnValue;
  }

  //---------------------------------------------------------------------------
  unsigned int SvmScaleParameterList::getKeyIndex(unsigned int lineIndex) const
  {
    unsigned int keyIndex = m_requiredHeaderIndex[0];

    unsigned int returnValue;
    sscanf(m_lines[lineIndex][keyIndex].c_str(), "%u", &returnValue);
    return returnValue;
  }
  //---------------------------------------------------------------------------
  unsigned int SvmScaleParameterList::size() const
  {
    return (unsigned int) m_nameToKey.size();
  }
}

