/*
 * SvmModel.cpp
 *
 *  Created on: Jan 1, 2011
 *      Author: jsnedecor
 */

/**
 * Load in svm model parameters
 */

#include "SvmModel.h"

namespace specnets
{
  //---------------------------------------------------------------------------
  SvmModel::SvmModel(void)
  {
    m_modelKeyRanges = (SvmScaleParameterList*) NULL;
  }

  //---------------------------------------------------------------------------
  SvmModel::~SvmModel(void)
  {
    //EMPTY
  }

  //---------------------------------------------------------------------------

  SvmModel::SvmModel(SvmScaleParameterList * modelKeyRanges)
  {
    m_modelKeyRanges = modelKeyRanges;
  }

  //---------------------------------------------------------------------------
  void SvmModel::setKeyRanges(SvmScaleParameterList * modelKeyRanges)
  {
    m_modelKeyRanges = modelKeyRanges;
  }

  //---------------------------------------------------------------------------

  bool SvmModel::loadSvmModel(const char * modelFile)
  {
    ifstream modelHandle(modelFile, ios::binary);

    if (!modelHandle.is_open() || !modelHandle.good())
    {
      ERROR_MSG("Unable to read file! " << modelFile);
      return false;
    }

    string lineBuff;

    bool passedSvmMarker = false;

    int numSupportVector = -1;

    while (!passedSvmMarker)
    {
      DelimitedTextReader::getlineChomp(modelHandle, lineBuff);

      //skip empty lines
      if (lineBuff.compare("") == 0)
      {
        continue;
      }

      vector<string> headerTemp;

      stringSplit(lineBuff, headerTemp, " ");

      if (headerTemp[0].compare("rho") == 0)
      {
        sscanf(headerTemp[1].c_str(), "%lf", &m_rho);
      }
      else if (headerTemp[0].compare("gamma") == 0)
      {
        sscanf(headerTemp[1].c_str(), "%lf", &m_gamma);
      }
      else if (headerTemp[0].compare("total_sv") == 0)
      {
        sscanf(headerTemp[1].c_str(), "%d", &numSupportVector);
      }
      else if (headerTemp[0].compare("SV") == 0)
      {
        passedSvmMarker = true;
      }
    }

    while (!modelHandle.eof())
    {
      DelimitedTextReader::getlineChomp(modelHandle, lineBuff);

      //skip empty lines
      if (lineBuff.compare("") == 0)
      {
        continue;
      }

      vector<string> fields;
      stringSplit(lineBuff, fields, " ");

      double scalingFactor;
      sscanf(fields[0].c_str(), "%lf", &scalingFactor);
      m_scalingFactors.push_back(scalingFactor);

      map<unsigned int, double> currLine;

      for (unsigned int i = 1; i < fields.size(); i++)
      {
        //split scaling factor
        vector<string> currFactor;
        stringSplit(fields[i], currFactor, ":");

        unsigned int keyIndex;
        double keyScalingFactor;
        sscanf(currFactor[0].c_str(), "%u", &keyIndex);
        sscanf(currFactor[1].c_str(), "%lf", &keyScalingFactor);
        currLine[keyIndex] = keyScalingFactor;
      }
      m_model.push_back(currLine);
    }
    return true;
  }
  //---------------------------------------------------------------------------
  double SvmModel::getRho() const
  {
    return m_rho;
  }

  //---------------------------------------------------------------------------
  double SvmModel::getGamma() const
  {
    return m_gamma;
  }

  //---------------------------------------------------------------------------
  double SvmModel::getScalingFactor(unsigned int i) const
  {
    return m_scalingFactors[i];
  }
  //---------------------------------------------------------------------------

  bool SvmModel::getKeyScalingFactor(unsigned int i,
                                     unsigned int keyIndex,
                                     double &returnScalingFactor)
  {
    map<unsigned int, double>::const_iterator currMap;

    currMap = m_model[i].find(keyIndex);

    if (currMap == m_model[i].end())
    {
      return false;
    }
    else
    {
      returnScalingFactor = currMap->second;
      return true;
    }
  }

  //---------------------------------------------------------------------------
  unsigned int SvmModel::size() const
  {
    return (unsigned int)m_model.size();
  }
  //---------------------------------------------------------------------------
  void SvmModel::scaleVector(vector<float> &inputVector,
                             vector<float> &outputVector)
  {
    outputVector.resize(inputVector.size());

    for (unsigned i = 0; i < m_modelKeyRanges->size(); i++)
    {
      unsigned int currKeyIndex = m_modelKeyRanges->getKeyIndex(i) - 1; //key indexes are 1 indexes.

      //DEBUG_VAR(currKeyIndex);
      float oldValue = inputVector[currKeyIndex];
      float newValue;

      if (oldValue <= m_modelKeyRanges->getStartRange(i))
      {
        newValue = -1;
      }
      else if (oldValue >= m_modelKeyRanges->getEndRange(i))
      {
        newValue = 1;
      }
      else
      {
        newValue = -1 + 2 * ((oldValue - m_modelKeyRanges->getStartRange(i))
            / (m_modelKeyRanges->getEndRange(i)
                - m_modelKeyRanges->getStartRange(i)));
      }
      //DEBUG_VAR(oldValue);
      //DEBUG_VAR(newValue);
      outputVector[currKeyIndex] = newValue;
    }
  }
  //---------------------------------------------------------------------------
  double SvmModel::getSvmScore(vector<float> &inputVector)
  {
    double score = 0;

    for (unsigned int i=0; i < m_model.size(); i++)
    {
      //iterate through vector
      double innerProduct = 0;

      for (unsigned int j=0; j < inputVector.size(); j++)
      {
        if (i == 0)
        {
          //DEBUG_VAR(inputVector[j]);
        }

        double modelValue;
        if (getKeyScalingFactor(i,j+1,modelValue)) //key indices are from 1
        {
          double difference = inputVector[j] - modelValue;
          innerProduct += pow(difference,2);
        }
        else
        {
          innerProduct += pow((double)inputVector[j],2);
        }
      }

      innerProduct = exp(-1 * getGamma() * innerProduct);
      score += getScalingFactor(i) * innerProduct;
    }

    score -= getRho();
    return score;
  }
  //---------------------------------------------------------------------------
  bool SvmModel::getKeyIndexByName(string &name, unsigned int &outputIndex)
  {
    return m_modelKeyRanges->getKeyIndexByName(name, outputIndex);
  }
  //---------------------------------------------------------------------------
  unsigned int SvmModel::getKeyIndexSize(void)
  {
    return m_modelKeyRanges->size();
  }

}

