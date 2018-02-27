/*
 * AssembledPeak.cpp
 *
 *  Created on: Apr 26, 2012
 *      Author: aguthals
 */

#include "AssembledPeak.h"

using namespace std;
using namespace specnets;

namespace abruijn
{

  const unsigned short AssembledPeak::BIN_VERSION = 1;
  const unsigned short AssembledPeak::BIN_SUBVERSION = 1;

  const string AssembledPeak::BIN_VERSION_ID = "AssembledPeak_binVersion";
  const string AssembledPeak::BIN_SUBVERSION_ID = "AssembledPeak_binSubVersion";

  AssembledPeak::AssembledPeak(void) :
    m_specID(""), m_endPt(false), m_BYsymmetry(Symmetry_unknown), MZRange()
  {
  }

  AssembledPeak::AssembledPeak(const AssembledPeak& other) :
    m_specID(other.m_specID), m_endPt(other.m_endPt),
        m_BYsymmetry(other.m_BYsymmetry), MZRange((MZRange&)other)
  {
  }

  AssembledPeak::AssembledPeak(const Spectrum& refSpec, unsigned int peakIdx) :
    m_specID(""), m_endPt(false), m_BYsymmetry(Symmetry_unknown), MZRange()
  {
    this->initialize(refSpec, peakIdx);
  }

  void AssembledPeak::clear()
  {
    m_specID = "";
    m_endPt = false;
    this->setMass(0);
    this->setIntensity(0);
    this->setTolerance(0);
    m_BYsymmetry = Symmetry_unknown;
  }

  AssembledPeak& AssembledPeak::operator=(const AssembledPeak& other)
  {
    m_specID = other.m_specID;
    m_endPt = other.m_endPt;
    m_BYsymmetry = other.m_BYsymmetry;
    MZRange::operator =((MZRange&)other);
    return (*this);
  }

  void AssembledPeak::initialize(const Spectrum& refSpec, unsigned int peakIdx)
  {
    m_specID = refSpec.getUniqueID();
    this->setMass(refSpec[peakIdx][0]);
    this->setIntensity(refSpec[peakIdx][1]);
    this->setTolerance(refSpec.getTolerance(peakIdx));
    m_endPt = false;
    m_BYsymmetry = Symmetry_unknown;

    double massB0 = 0, massBk = refSpec.parentMass - AAJumps::massMH, massY0 =
        AAJumps::massH2O, massYk = refSpec.parentMass - AAJumps::massHion;

    if (this->operator ==(massB0) || this->operator ==(massBk))
    {
      m_endPt = true;
      m_BYsymmetry = Symmetry_B;
    }
    else if (this->operator ==(massY0) || this->operator ==(massYk))
    {
      m_endPt = true;
      m_BYsymmetry = Symmetry_Y;
    }
  }

  void AssembledPeak::initializeEmpty(const string& refSpecID, float peakMass)
  {
    m_specID = refSpecID;
    this->setMass(peakMass);
    this->setIntensity(-1.0);
    this->setTolerance(0);
    m_endPt = false;
    m_BYsymmetry = Symmetry_B;
  }

  bool AssembledPeak::saveToBinaryStream(FILE* fp) const
  {
    unsigned int count;
    double mass = getMass();
    double intensity = getIntensity();
    double tolerance = getTolerance();
    unsigned short symm = (unsigned short)m_BYsymmetry;

    count = fwrite(&m_endPt, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_endPt");
      return false;
    }

    count = fwrite(&symm, sizeof(unsigned short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_BYsymmetry");
      return false;
    }

    count = fwrite(&mass, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the mass");
      return false;
    }

    count = fwrite(&intensity, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the intensity");
      return false;
    }

    count = fwrite(&tolerance, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write the tolerance");
      return false;
    }

    vector<string> labels(1);
    labels[0] = m_specID;

    if (!writeStringsToBinaryStream(fp, labels))
    {
      ERROR_MSG("Could not write the spectrum ID");
      return false;
    }

    return true;
  }

  bool AssembledPeak::loadFromBinaryStream(FILE* fp,
                                           map<string, unsigned short>& versions)
  {
    unsigned short version = versions[BIN_VERSION_ID];
    unsigned short subVersion = versions[BIN_SUBVERSION_ID];

    unsigned int count;
    double mass, intensity, tolerance;
    unsigned short symm;

    count = fread(&m_endPt, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_endPt");
      return false;
    }

    count = fread(&symm, sizeof(unsigned short), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_BYsymmetry");
      return false;
    }

    if (symm == (unsigned short)Symmetry_unknown)
    {
      m_BYsymmetry = Symmetry_unknown;
    }
    else if (symm == (unsigned short)Symmetry_B)
    {
      m_BYsymmetry = Symmetry_B;
    }
    else if (symm == (unsigned short)Symmetry_Y)
    {
      m_BYsymmetry = Symmetry_Y;
    }
    else
    {
      ERROR_MSG("Found unknown m_BYsymmetry \'" << symm << "\'");
      return false;
    }

    count = fread(&mass, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the mass");
      return false;
    }
    setMass(mass);

    count = fread(&intensity, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the intensity");
      return false;
    }
    setIntensity(intensity);

    count = fread(&tolerance, sizeof(double), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read the tolerance");
      return false;
    }
    setTolerance(tolerance);

    vector<string> labels;

    if (!readStringsFromBinaryStream(fp, labels) || labels.size() != 1)
    {
      ERROR_MSG("Could not read the spectrum ID");
      return false;
    }
    m_specID = labels[0];

    return true;
  }

}
