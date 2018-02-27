/*
 * MassErrorModel.cpp
 *
 *  Created on: Mar 4, 2014
 *      Author: aguthals
 */

#include <math.h>

#include "MassErrorModel.h"

namespace specnets
{
  const unsigned short MassErrorModel::BIN_VERSION = 1;
  const unsigned short MassErrorModel::BIN_SUBVERSION = 1;

  const string MassErrorModel::BIN_VERSION_ID = "MassErrorModel_binVersion";
  const string MassErrorModel::BIN_SUBVERSION_ID =
      "MassErrorModel_binSubVersion";

  const float MassErrorModel::BIN_WIDTH_PPM = 4.0;

  const float MassErrorModel::MAX_PPM_ERROR = 1000.0;

  float MassErrorModel::getPPMMassError(const Spectrum &msSpec,
                                        const Spectrum &prmSpec,
                                        const float prmMass1,
                                        const float srmMass1,
                                        const float prmMass2,
                                        const float srmMass2,
                                        const float aaMass,
                                        const float roundingFactor)
  {
    const float massBOffset =
        (msSpec.msFragType == Spectrum::FragType_ETD) ?
            18.0344 * roundingFactor : 1.0078 * roundingFactor;

    const float massYOffset =
        (msSpec.msFragType == Spectrum::FragType_ETD) ?
            1.9918 * roundingFactor : 19.0184 * roundingFactor;

    list<int> matchesB1;
    if (prmMass1 > 50.0)
    {
      msSpec.findPeaks(prmMass1 + massBOffset, 0.5, &matchesB1);
    }

    list<int> matchesB2;
    msSpec.findPeaks(prmMass2 + massBOffset, 0.5, &matchesB2);

    list<int> matchesY1;
    msSpec.findPeaks(srmMass1 + massYOffset, 0.5, &matchesY1);

    list<int> matchesY2;
    msSpec.findPeaks(srmMass2 + massYOffset, 0.5, &matchesY2);

    float minError = POS_INF;

    if ((prmMass1 > 50.0 && prmSpec.findPeaks(prmMass1, 0.5) < 0)
        || (prmMass2
            < (prmSpec.parentMass - (AAJumps::massMH * roundingFactor) - 50.0)
            && prmSpec.findPeaks(prmMass2, 0.5) < 0))
    {
      return minError;
    }

    if (prmMass1 <= 50.0)
    {
      for (list<int>::const_iterator b2It = matchesB2.begin();
          b2It != matchesB2.end(); b2It++)
      {
        int p2 = *b2It;
        const float massError = ((aaMass - abs(msSpec[p2][0] - massBOffset))
            / msSpec[p2][0]) / PPM_FACTOR;
        if (abs(massError) < abs(minError))
        {
          minError = massError;
        }
      }
    }

    for (list<int>::const_iterator b1It = matchesB1.begin();
        b1It != matchesB1.end(); b1It++)
    {
      int p1 = *b1It;

      for (list<int>::const_iterator b2It = matchesB2.begin();
          b2It != matchesB2.end(); b2It++)
      {
        int p2 = *b2It;

        const float avgMass = (msSpec[p2][0] + msSpec[p1][0]) / 2.0;
        const float massError = ((aaMass - abs(msSpec[p2][0] - msSpec[p1][0]))
            / avgMass) / PPM_FACTOR;

        if (abs(massError) < abs(minError))
        {
          minError = massError;
        }
      }
    }

    for (list<int>::const_iterator y1It = matchesY1.begin();
        y1It != matchesY1.end(); y1It++)
    {
      int p1 = *y1It;

      for (list<int>::const_iterator y2It = matchesY2.begin();
          y2It != matchesY2.end(); y2It++)
      {
        int p2 = *y2It;
        const float avgMass = (msSpec[p2][0] + msSpec[p1][0]) / 2.0;
        const float massError = ((aaMass - abs(msSpec[p2][0] - msSpec[p1][0]))
            / avgMass) / PPM_FACTOR;

        if (abs(massError) < abs(minError))
        {
          minError = massError;
        }
      }
    }

    return minError;
  }

  void MassErrorModel::getTDAHistogram(const MassErrorModel &targetDist,
                                       const MassErrorModel &decoyDist,
                                       OutputTable &tableOut)
  {
    if (targetDist.size() != decoyDist.size() || targetDist.size() == 0)
    {
      WARN_MSG("Invalid target/decoy distribution sizes " << targetDist.size() << "/" << decoyDist.size());
      return;
    }
    if (targetDist.getZeroBin() != decoyDist.getZeroBin())
    {
      WARN_MSG("Invalid target/decoy zero bins " << targetDist.getZeroBin() << "/" << decoyDist.getZeroBin());
      return;
    }
    tableOut.values.resize(targetDist.size() + 1);

    const int numCols = 3;
    int row = 0, col = 0;

    tableOut.values[row].resize(numCols);
    tableOut.values[row][col].second = true;
    tableOut.values[row][col++].first = "PPM Bin";
    tableOut.values[row][col].second = true;
    tableOut.values[row][col++].first = "Target Prob Density";
    tableOut.values[row][col].second = true;
    tableOut.values[row][col++].first = "Decoy Prob Density";

    col = 0;
    row++;
    for (int i = 0; i < targetDist.size(); i++)
    {
      tableOut.values[row].resize(numCols);

      tableOut.values[row][col].second = false;
      tableOut.values[row][col++].first =
          parseFloat(targetDist.getPPMCenterBin(i), 1);

      tableOut.values[row][col].second = false;
      tableOut.values[row][col++].first = parseDoubleSci(targetDist[i], 5);

      tableOut.values[row][col].second = false;
      tableOut.values[row][col++].first = parseDoubleSci(decoyDist[i], 5);

      col = 0;
      row++;
    }
  }

  float MassErrorModel::getEdgeScore(const Spectrum &msSpec,
                                     const Spectrum &prmSpec,
                                     const float prmMass1,
                                     const float srmMass1,
                                     const float prmMass2,
                                     const float srmMass2,
                                     const float aaMass,
                                     const float maxPPMError,
                                     const MassErrorModel &targetDist,
                                     const MassErrorModel &decoyDist,
                                     const float roundingFactor)
  {
    if (prmMass1 < 0 || prmMass2 < 0)
    {
      return getAbsentEdgeScore(targetDist, decoyDist);
    }
    /*else if (prmMass1 < 0 || prmMass2 < 0)
     {
     return getDangleEdgeScore(targetDist, decoyDist);
     }*/

    float ppmError = getPPMMassError(msSpec,
                                     prmSpec,
                                     prmMass1,
                                     srmMass1,
                                     prmMass2,
                                     srmMass2,
                                     aaMass,
                                     roundingFactor);

    int absoluteBinTarget = targetDist.getAbsoluteBin(ppmError);

    int absoluteBinDecoy = decoyDist.getAbsoluteBin(ppmError);

    if (abs(ppmError) > maxPPMError || absoluteBinTarget < 0
        || absoluteBinDecoy < 0)
    {
      return getAbsentEdgeScore(targetDist, decoyDist);
    }

    return log((double)(targetDist[absoluteBinTarget]
        / decoyDist[absoluteBinDecoy]));
  }

  MassErrorModel::MassErrorModel() :
      m_numerators(0), m_denomintor(0), m_probDensity(0), m_zeroBin(-1), m_edgesPresent(0), m_dangleEdgesPresent(0), m_dangleEdgeProb(0), m_edgeProb(0), m_binWidthPPM(BIN_WIDTH_PPM), m_maxPPMError(MAX_PPM_ERROR)
  {

  }

  MassErrorModel::MassErrorModel(const SpecSet &identifiedSpecsMS,
                                 const SpecSet &specsPRM,
                                 const AAJumps &aaJumps,
                                 const bool useDecoy) :
      m_numerators(0), m_denomintor(0), m_probDensity(0), m_zeroBin(-1), m_edgesPresent(0), m_edgeProb(0), m_dangleEdgesPresent(0), m_dangleEdgeProb(0), m_binWidthPPM(BIN_WIDTH_PPM), m_maxPPMError(MAX_PPM_ERROR)
  {
    initialize(identifiedSpecsMS, specsPRM, aaJumps, useDecoy);
  }

  MassErrorModel::~MassErrorModel()
  {
  }

  void MassErrorModel::initialize(const SpecSet &identifiedSpecsMS,
                                  const SpecSet &specsPRM,
                                  const AAJumps &aaJumps,
                                  const bool useDecoy)
  {
    const int numBins = (floatToInt(m_maxPPMError / m_binWidthPPM) * 2) - 1;

    m_numerators.assign(numBins, 0);
    m_denomintor = 2;
    m_edgesPresent = 1;
    m_dangleEdgesPresent = 1;
    m_probDensity.assign(numBins, 0);

    m_zeroBin = (numBins - 1) / 2;

    for (int i = 0; i < identifiedSpecsMS.size(); i++)
    {
      const Spectrum &inSpec = identifiedSpecsMS[i];

      if (inSpec.psmList.size() == 0)
      {
        continue;
      }

      const bool isDecoy = inSpec.psmList.front()->m_isDecoy;

      if (useDecoy && !isDecoy)
      {
        continue;
      }
      if ((!useDecoy) && (isDecoy || inSpec.psmList.front()->m_fdr > 0.01))
      {
        continue;
      }

      const string peptide = inSpec.psmList.front()->m_annotation;

      if (aaJumps.modsUsed == 0 && AAJumps::isModified(peptide))
      {
        continue;
      }

      addPPMMassErrors(inSpec, specsPRM[i], peptide, aaJumps);
    }

    m_edgeProb = ((long double)m_edgesPresent) / ((long double)m_denomintor);

    m_dangleEdgeProb = ((long double)m_dangleEdgesPresent)
        / ((long double)m_denomintor);

    for (int i = 0; i < m_numerators.size(); i++)
    {
      m_probDensity[i] = ((long double)m_numerators[i])
          / ((long double)m_denomintor);

      if (m_numerators[i] == 0)
      {
        m_probDensity[i] = 1 / ((float)numBins);
      }
    }
    DEBUG_VAR(m_denomintor);
    DEBUG_VAR(m_edgesPresent);
    DEBUG_VAR(m_edgeProb);
    DEBUG_VAR(m_dangleEdgesPresent);
    DEBUG_VAR(m_dangleEdgeProb);
  }

  bool MassErrorModel::saveToBinaryStream(FILE* fp) const
  {
    if (fp == 0)
    {
      return false;
    }

    unsigned int count;

    unsigned int numBins = m_numerators.size();

    count = fwrite(&numBins, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to # of bins");
      return false;
    }

    count = fwrite(&m_zeroBin, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save zero bin");
      return false;
    }

    count = fwrite(&m_binWidthPPM, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save PPM bin width");
      return false;
    }

    count = fwrite(&m_maxPPMError, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save max PPM error");
      return false;
    }

    count = fwrite(&m_denomintor, sizeof(long long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save denominator");
      return false;
    }

    count = fwrite(&m_edgesPresent, sizeof(long long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save # of edges");
      return false;
    }

    count = fwrite(&m_edgeProb, sizeof(long double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save edge probability");
      return false;
    }

    count = fwrite(&m_dangleEdgesPresent, sizeof(long long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save # of edges");
      return false;
    }

    count = fwrite(&m_dangleEdgeProb, sizeof(long double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to save edge probability");
      return false;
    }

    if (numBins == 0)
    {
      return true;
    }

    long long* numeratorBuf = (long long*)malloc(sizeof(long long) * numBins);

    for (unsigned int i = 0; i < numBins; i++)
    {
      numeratorBuf[i] = m_numerators[i];
    }
    count = fwrite(numeratorBuf, sizeof(long long), numBins, fp);
    free(numeratorBuf);
    if (count == 0)
    {
      ERROR_MSG("Failed to save numerators");
      return false;
    }

    long double* probBuf = (long double*)malloc(sizeof(long double) * numBins);

    for (unsigned int i = 0; i < numBins; i++)
    {
      probBuf[i] = m_probDensity[i];
    }
    count = fwrite(probBuf, sizeof(long double), numBins, fp);
    free(probBuf);
    if (count == 0)
    {
      ERROR_MSG("Failed to save probabilities");
      return false;
    }

    return true;
  }

  bool MassErrorModel::loadFromBinaryStream(FILE* fp,
                                            map<string, unsigned short>& versions)
  {
    if (fp == 0)
    {
      return false;
    }
    unsigned short version = versions[BIN_VERSION_ID];
    unsigned short subVersion = versions[BIN_SUBVERSION_ID];

    if (version > BIN_VERSION
        || (version == BIN_VERSION && subVersion > BIN_SUBVERSION))
    {
      ERROR_MSG("Unsupported MassErrorModel version " << version << "." << subVersion << " (this release supports up to " << BIN_VERSION << "." << BIN_SUBVERSION << "), you must obtain a later release of the code base to load this file.");
      return false;
    }

    unsigned int count;

    unsigned int numBins;

    count = fread(&numBins, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read # of bins");
      return false;
    }

    count = fread(&m_zeroBin, sizeof(int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read zero bin");
      return false;
    }

    count = fread(&m_binWidthPPM, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read PPM bin width");
      return false;
    }

    count = fread(&m_maxPPMError, sizeof(float), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read max PPM error");
      return false;
    }

    count = fread(&m_denomintor, sizeof(long long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read denominator");
      return false;
    }

    count = fread(&m_edgesPresent, sizeof(long long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read # of edges");
      return false;
    }

    count = fread(&m_edgeProb, sizeof(long double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read edge probability");
      return false;
    }

    count = fread(&m_dangleEdgesPresent, sizeof(long long), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read # of edges");
      return false;
    }

    count = fread(&m_dangleEdgeProb, sizeof(long double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read edge probability");
      return false;
    }

    m_numerators.resize(numBins);
    m_probDensity.resize(numBins);
    if (numBins == 0)
    {
      return true;
    }

    long long* numeratorBuf = (long long*)malloc(sizeof(long long) * numBins);
    count = fread(numeratorBuf, sizeof(long long), numBins, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read numerators");
      free(numeratorBuf);
      return false;
    }

    for (unsigned int i = 0; i < numBins; i++)
    {
      m_numerators[i] = numeratorBuf[i];
    }

    free(numeratorBuf);

    long double* probBuf = (long double*)malloc(sizeof(long double) * numBins);
    count = fread(probBuf, sizeof(long double), numBins, fp);
    if (count == 0)
    {
      ERROR_MSG("Failed to read probabilities");
      free(probBuf);
      return false;
    }

    for (unsigned int i = 0; i < numBins; i++)
    {
      m_probDensity[i] = probBuf[i];
    }

    free(probBuf);

    return true;
  }

  bool MassErrorModel::saveBinaryFile(const string& filename)
  {
    FILE* fp = fopen(filename.c_str(), "wb");
    if (fp == 0)
    {
      ERROR_MSG("Error opening \'" << filename << "\' for writing !");
      return false;
    }

    map<string, unsigned short> versions;
    versions[BIN_VERSION_ID] = BIN_VERSION;
    versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;

    if (!writeStringMapToBinaryStream<unsigned short>(fp, versions))
    {
      ERROR_MSG("Error saving version info for MassErrorModel");
      fclose(fp);
      return false;
    }

    if (!saveToBinaryStream(fp))
    {
      ERROR_MSG("Error saving MassErrorModel");
      fclose(fp);
      return false;
    }
    fclose(fp);
    return true;
  }

  bool MassErrorModel::loadBinaryFile(const string& filename)
  {
    FILE* fp = fopen(filename.c_str(), "rb");
    if (fp == 0)
    {
      ERROR_MSG("Error opening \'" << filename << "\' for reading !");
      return false;
    }

    map<string, unsigned short> versions;
    if (!readStringMapFromBinaryStream<unsigned short>(fp, versions))
    {
      ERROR_MSG("Error reading version info for MassErrorModel");
      fclose(fp);
      return false;
    }

    if (!loadFromBinaryStream(fp, versions))
    {
      ERROR_MSG("Error loading MassErrorModel'");
      fclose(fp);
      return false;
    }
    fclose(fp);
    return true;
  }

  int MassErrorModel::getAbsoluteBin(const float &ppmError) const
  {
    const float absError = abs(ppmError);

    const float halfBin = m_binWidthPPM / 2.0;

    const int relativeBin =
        (absError <= halfBin) ? 0 :
            floatToInt(floor((absError - halfBin) / m_binWidthPPM)) + 1;

    const int absoluteBin =
        (ppmError < 0) ? m_zeroBin - relativeBin : m_zeroBin + relativeBin;

    if (absoluteBin < 0 || absoluteBin >= m_numerators.size())
    {
      return -1;
    }
    return absoluteBin;
  }

  void MassErrorModel::addPPMError(const float &ppmError)
  {
    m_denomintor++;
    if (abs(ppmError) <= m_maxPPMError)
    {
      const int absoluteBin = getAbsoluteBin(ppmError);

      if (absoluteBin >= 0)
      {
        //DEBUG_VAR(absoluteBin);
        m_numerators[absoluteBin]++;
        m_edgesPresent++;
      }
    }
  }

  void MassErrorModel::addPPMMassErrors(const Spectrum &msSpec,
                                        const Spectrum &prmSpec,
                                        const string &peptide,
                                        const AAJumps &aaJumps)
  {
    vector<float> prmMasses0;
    // Get PRM masses, add zero mass
    aaJumps.getPRMMasses(peptide, prmMasses0, 0, 0, true);

    Spectrum msSpecCopy = msSpec;
    msSpecCopy.setPeakTolerance(0);
    msSpecCopy.rankFilterPeaks(10);

    const float pepMass = prmMasses0[prmMasses0.size() - 1];

    for (int mIdx = 0; mIdx < prmMasses0.size() - 1; mIdx++)
    {
      float prm1 = prmMasses0[mIdx];
      float srm1 = pepMass - prm1;
      float prm2 = prmMasses0[mIdx + 1];
      float srm2 = pepMass - prm2;
      float aaMass = prm2 - prm1;

      float ppmError = getPPMMassError(msSpecCopy,
                                       prmSpec,
                                       prm1,
                                       srm1,
                                       prm2,
                                       srm2,
                                       aaMass);

      if (ppmError > POS_INF / 10.0 && prm1 > 50.0
          && prm2 < prmSpec.parentMass - AAJumps::massMH - 50.0)
      {
        bool found1 = prmSpec.findPeaks(prm1, 0.5) >= 0;

        bool found2 = prmSpec.findPeaks(prm2, 0.5) >= 0;

        if ((found1 && (!found2)) || ((!found1) && found2))
        {
          m_dangleEdgesPresent++;
        }
      }
      addPPMError(ppmError);
    }
  }
}

