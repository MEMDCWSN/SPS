/*
 * AbruijnNode.cpp
 *
 *  Created on: Apr 26, 2012
 *      Author: aguthals
 */

#include "AbruijnNode.h"
#include "AbruijnGraph.h"

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_set>
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_set>
#  include <unordered_map>
#endif

using namespace std;
using namespace specnets;

namespace abruijn
{

  const unsigned short AbruijnNode::BIN_VERSION = 1;
  const unsigned short AbruijnNode::BIN_SUBVERSION = 1;

  const string AbruijnNode::BIN_VERSION_ID = "AbruijnNode_binVersion";
  const string AbruijnNode::BIN_SUBVERSION_ID = "AbruijnNode_binSubVersion";

  AbruijnNode* AbruijnNode::castNodePtr(Node* nodePtr)
  {
    AbruijnNode* temp = dynamic_cast<AbruijnNode*> (nodePtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast Node \'" << nodePtr->toString() << "\' to an AbruijnNode");
      abort();
    }
    return temp;
  }

  AbruijnNode::AbruijnNode(void) :
    Node(), m_isCorrect(false), m_isGreen(false), m_isChimeric(false),
        m_isAnnotated(false), m_fPath(""), m_rPath(""), m_peaks(0)
  {
  }

  AbruijnNode::AbruijnNode(const AbruijnNode& other) :
    Node(), m_isCorrect(false), m_isGreen(false), m_isChimeric(false),
        m_isAnnotated(false), m_fPath(""), m_rPath(""), m_peaks(0)
  {
    this->operator =(other);
  }

  string AbruijnNode::getGraphvizLabel(void) const
  {
    ostringstream out;
    out << Node::getGraphvizLabel();
    for (unsigned int i = 0; i < size(); i++)
    {
      out << " : (" << this->operator [](i).getSpecID() << ","
          << this->operator [](i).getMass() << ")";
    }
    return out.str();
  }

  string AbruijnNode::getGraphvizFillColor(void) const
  {
    return (this->isGreen()) ? "green" : "red";
  }

  void AbruijnNode::initialize(unsigned int sz)
  {
    m_peaks.resize(sz);
    m_isCorrect = false;
    m_isChimeric = false;
    m_isAnnotated = false;
    m_fPath = "";
    m_rPath = "";
    m_isGreen = false;
  }

  bool AbruijnNode::haveRealPeak()
  {
    int sz = m_peaks.size();
    bool foundPeak = false;
    for (unsigned int j = 0; j < sz; j++)
    {
      float inten = m_peaks[j].getIntensity();
      const AssembledPeak& aPeak = m_peaks[j];
      const string& specID = m_peaks[j].getSpecID();
      if (m_peaks[j].getIntensity() > -0.1
          && m_peaks[j].getSpecID().find(AbruijnGraph::ENDPT_SPEC_ID)
              == string::npos
          && m_peaks[j].getSpecID().find(AbruijnGraph::STARTPT_SPEC_ID)
              == string::npos)
      {
        foundPeak = true;
        break;
      }
    }
    return foundPeak;
  }

  AbruijnNode & AbruijnNode::operator=(const AbruijnNode &other)
  {
    if (this == &other)
    {
      return *this;
    }

    Node::operator=((const Node&)other);

    if (size() != other.size())
    {
      m_peaks.resize(other.size());
    }

    for (unsigned int i = 0; i < size(); i++)
    {
      m_peaks[i] = other[i];
    }

    m_isCorrect = other.m_isCorrect;
    m_isChimeric = other.m_isChimeric;
    m_isAnnotated = other.m_isAnnotated;
    m_fPath = other.m_fPath;
    m_rPath = other.m_rPath;
    m_isGreen = other.m_isGreen;

    return *this;
  }

  void AbruijnNode::copy(Node& otherNode)
  {
    this->operator =(*castNodePtr(&otherNode));
  }

  bool AbruijnNode::saveToBinaryStream(FILE* fp) const
  {
    if (!Node::saveToBinaryStream(fp))
    {
      return false;
    }
    unsigned int count;

    count = fwrite(&m_isCorrect, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_isCorrect");
      return false;
    }

    count = fwrite(&m_isChimeric, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_isChimeric");
      return false;
    }

    count = fwrite(&m_isGreen, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_isGreen");
      return false;
    }

    count = fwrite(&m_isAnnotated, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write m_isAnnotated");
      return false;
    }

    unsigned int numPeaks = size();
    count = fwrite(&numPeaks, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not write numPeaks");
      return false;
    }

    for (unsigned int i = 0; i < numPeaks; i++)
    {
      if (!m_peaks[i].saveToBinaryStream(fp))
      {
        ERROR_MSG("Failed to save AssembledPeak " << i);
        return false;
      }
    }

    return true;
  }

  bool AbruijnNode::loadFromBinaryStream(FILE* fp,
                                         map<string, unsigned short>& versions)
  {
    if (!Node::loadFromBinaryStream(fp, versions))
    {
      return false;
    }
    unsigned int count;

    count = fread(&m_isCorrect, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_isCorrect");
      return false;
    }

    count = fread(&m_isChimeric, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_isChimeric");
      return false;
    }

    count = fread(&m_isGreen, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_isGreen");
      return false;
    }

    count = fread(&m_isAnnotated, sizeof(bool), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read m_isAnnotated");
      return false;
    }

    unsigned int numPeaks = 0;
    count = fread(&numPeaks, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      ERROR_MSG("Could not read numPeaks");
      return false;
    }
    m_peaks.resize(numPeaks);

    for (unsigned int i = 0; i < numPeaks; i++)
    {
      if (!m_peaks[i].loadFromBinaryStream(fp, versions))
      {
        ERROR_MSG("Failed to load AssembledPeak " << i);
        return false;
      }
    }
    //clearLabelIndex();

    return true;
  }

  bool AbruijnNode::isComposite() const
  {
    tr1::unordered_set<string> uniqueSpecs(size());
    for (unsigned int i = 0; i < m_peaks.size(); i++)
    {
      uniqueSpecs.insert(m_peaks[i].getSpecID());
    }
    return uniqueSpecs.size() < m_peaks.size();
  }

  bool AbruijnNode::isCompositeWithOther(const AbruijnNode& other) const
  {
    tr1::unordered_map<string, MZRange> uniqueSpecs(size() + other.size());
    for (unsigned int i = 0; i < m_peaks.size(); i++)
    {
      const string& specID = m_peaks[i].getSpecID();
      if (uniqueSpecs.count(specID) > 0)
      {
        return true;
      }
      uniqueSpecs[specID] = m_peaks[i].getTightBound();
    }
    for (unsigned int i = 0; i < other.m_peaks.size(); i++)
    {
      const string& specID = other.m_peaks[i].getSpecID();
      MZRange tightBound = other.m_peaks[i].getTightBound();
      if (uniqueSpecs.count(specID) > 0)
      {
        return true;
      }
      uniqueSpecs[specID] = tightBound;
    }
    return false;
  }

  void AbruijnNode::mergeWith(const AbruijnNode& other)
  {
    unsigned int idxUse = size();
    this->resize(size() + other.size());
    for (unsigned int i = 0; i < other.size(); i++)
    {
      m_peaks[idxUse++] = other[i];
    }
  }

  bool AbruijnNode::isEndPoint() const
  {
    for (unsigned int i = 0; i < m_peaks.size(); i++)
    {
      if (!m_peaks[i].isEndPt())
      {
        return false;
      }
    }
    return true;
  }

  bool AbruijnNode::isAllYSymmetric()
  {
    bool foundY = false;
    bool foundB = false;

    for (unsigned int i = 0; i < size(); i++)
    {
      if (m_peaks[i].getSymmetry() == AssembledPeak::Symmetry_B)
      {
        foundB = true;
      }
      else if (m_peaks[i].getSymmetry() == AssembledPeak::Symmetry_Y)
      {
        foundY = true;
      }
    }
    return foundY && !foundB;
  }
}

