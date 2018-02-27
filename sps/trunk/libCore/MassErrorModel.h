/*
 * MassErrorModel.h
 *
 *  Created on: Mar 4, 2014
 *      Author: aguthals
 */

#ifndef MASSERRORMODEL_H_
#define MASSERRORMODEL_H_

#include <math.h>

#include "SpecSet.h"
#include "OutputTable.h"

namespace specnets
{

  class MassErrorModel
  {
  protected:

    vector<long long> m_numerators;

    long long m_denomintor;

    long long m_edgesPresent;

    long long m_dangleEdgesPresent;

    long double m_edgeProb;

    long double m_dangleEdgeProb;

    vector<long double> m_probDensity;

    int m_zeroBin;

    float m_binWidthPPM;

    float m_maxPPMError;

  public:

    static const float BIN_WIDTH_PPM;

    static const float MAX_PPM_ERROR;

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    static float getPPMMassError(const Spectrum &msSpec,
                                 const Spectrum &prmSpec,
                                 const float prmMass1,
                                 const float srmMass1,
                                 const float prmMass2,
                                 const float srmMass2,
                                 const float aaMass,
                                 const float roundingFactor = 1.0);

    static void getTDAHistogram(const MassErrorModel &targetDist,
                                const MassErrorModel &decoyDist,
                                OutputTable &tableOut);

    static float getAbsentEdgeScore(const MassErrorModel &targetDist,
                                    const MassErrorModel &decoyDist)
    {
      return log((1.0 - (double)targetDist.m_edgeProb)
          / (1.0 - (double)decoyDist.m_edgeProb));
    }

    static float getDangleEdgeScore(const MassErrorModel &targetDist,
                                    const MassErrorModel &decoyDist)
    {
      return log(((double)targetDist.m_dangleEdgeProb)
          / ((double)decoyDist.m_dangleEdgeProb));
    }

    static float getBaseEdgeScore(const MassErrorModel &targetDist,
                                  const MassErrorModel &decoyDist)
    {
      return log(((double)targetDist.m_edgeProb)
          / ((double)decoyDist.m_edgeProb));
    }

    static float getEdgeScore(const Spectrum &msSpec,
                              const Spectrum &prmSpec,
                              const float prmMass1,
                              const float srmMass1,
                              const float prmMass2,
                              const float srmMass2,
                              const float aaMass,
                              const float maxPPMError,
                              const MassErrorModel &targetDist,
                              const MassErrorModel &decoyDist,
                              const float roundingFactor = 1.0);

    MassErrorModel();

    MassErrorModel(const SpecSet &identifiedSpecsMS,
                   const SpecSet &specsPRM,
                   const AAJumps &aaJumps,
                   const bool useDecoy);

    ~MassErrorModel();

    void initialize(const SpecSet &identifiedSpecsMS,
                    const SpecSet &specsPRM,
                    const AAJumps &aaJumps,
                    const bool useDecoy);

    bool saveToBinaryStream(FILE* fp) const;

    bool loadFromBinaryStream(FILE* fp, map<string, unsigned short>& versions);

    bool saveBinaryFile(const string& filename);

    bool loadBinaryFile(const string& filename);

    unsigned int size() const
    {
      return m_probDensity.size();
    }

    int getAbsoluteBin(const float &ppmError) const;

    long double getProbDensity(const float &ppmError) const
    {
      int absoluteBin = getAbsoluteBin(ppmError);
      if (absoluteBin < 0)
      {
        return -1.0;
      }
      return m_probDensity[absoluteBin];
    }

    inline const int &getZeroBin() const
    {
      return m_zeroBin;
    }

    float getPPMCenterBin(const int &absoluteBin) const
    {
      float relativeBin = ((float)absoluteBin) - ((float)m_zeroBin);

      return relativeBin * m_binWidthPPM;
    }

    void addPPMError(const float &ppmError);

    inline long double &operator[](const int &absoluteBin)
    {
      return m_probDensity[absoluteBin];
    }

    inline const long double &operator[](const int &absoluteBin) const
    {
      return m_probDensity[absoluteBin];
    }

    void addPPMMassErrors(const Spectrum &msSpec,
                          const Spectrum &prmSpec,
                          const string &peptide,
                          const AAJumps &aaJumps);

  };
}

#endif /* MASSERRORMODEL_H_ */
