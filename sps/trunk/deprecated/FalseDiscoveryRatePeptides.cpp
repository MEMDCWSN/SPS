/*
 * FalseDiscoveryRatePeptide.cpp
 *
 *  Created on: Feb 14, 2011
 *      Author: jsnedecor
 */

#include "FalseDiscoveryRatePeptides.h"

namespace specnets
{
  //---------------------------------------------------------------------------

  FalseDiscoveryRatePeptides::FalseDiscoveryRatePeptides()
  {
    m_calculatedPValue = false;
  }

  //---------------------------------------------------------------------------

  FalseDiscoveryRatePeptides::~FalseDiscoveryRatePeptides()
  {
    //empty
  }
  //---------------------------------------------------------------------------

  void FalseDiscoveryRatePeptides::addFdrPeptide(FalseDiscoveryRatePeptide &input)
  {
    m_inputPeptides.push_back(input);
    m_calculatedPValue = false;
  }
  //---------------------------------------------------------------------------

  bool compareFDR(FalseDiscoveryRatePeptide i, FalseDiscoveryRatePeptide j)
  {
    return i.m_score > j.m_score;
  }
  //---------------------------------------------------------------------------

  bool FalseDiscoveryRatePeptides::calculatePValues()
  {
    sort(m_inputPeptides.begin(), m_inputPeptides.end(), compareFDR);

    unsigned int correctHits = 0;
    unsigned int incorrectHits = 0;

    m_pValues.resize(m_inputPeptides.size());

    for (int i = 0; i < m_inputPeptides.size(); i++)
    {
      if (m_inputPeptides[i].m_isDecoy)
      {
        incorrectHits++;
      }
      else
      {
        correctHits++;
      }
      m_pValues[i] = 1 - (((double)correctHits - (double)incorrectHits) / (double)correctHits);
    }
    m_calculatedPValue = true;
  }
  //---------------------------------------------------------------------------

  bool FalseDiscoveryRatePeptides::getPValueByScan(unsigned int scanNum,
                                                   double &score)
  {
    if (m_calculatedPValue)
    {
      for (int i = 0; i < m_inputPeptides.size(); i++)
      {
        if (m_inputPeptides[i].m_scanNum == scanNum)
        {
          score = m_inputPeptides[i].m_score;
          return true;
        }
      }
      return false;
    }
    else
    {
      WARN_MSG("Haven't calculated pvalues for FalseDiscoveryRate!");
      return false;
    }
  }
  //---------------------------------------------------------------------------

  bool FalseDiscoveryRatePeptides::getPValueByIndex(unsigned int indexNum,
                                                    double &score)
  {
    if (m_calculatedPValue)
    {
      if (m_inputPeptides.size() > indexNum)
      {
        score = m_pValues[indexNum];
        return true;
      }
      else
      {
        WARN_MSG("Index " << indexNum << " out of range of m_inputPeptides");
        return false;
      }
    }
    else
    {
      WARN_MSG("Haven't calculated pvalues for FalseDiscoveryRate!");
      return false;
    }
    return false;
  }
  //---------------------------------------------------------------------------

  bool FalseDiscoveryRatePeptides::getMQScoreByIndex(unsigned int indexNum, double &score)
  {
    score = m_inputPeptides[indexNum].m_score;
    return true;
  }

  //---------------------------------------------------------------------------

  InspectResultsLine * FalseDiscoveryRatePeptides::getResultByIndex(unsigned int indexNum)
  {
    return m_inputPeptides[indexNum].m_peptideResult;
  }

  //---------------------------------------------------------------------------

  bool FalseDiscoveryRatePeptides::getScanNumByIndex(unsigned int indexNum,
                                                     unsigned int &scanNum)
  {
    scanNum = m_inputPeptides[indexNum].m_scanNum;
    return true;
  }
  //---------------------------------------------------------------------------
  unsigned int FalseDiscoveryRatePeptides::size()
  {
    return m_inputPeptides.size();
  }

}

