/*
 * AbruijnNode.h
 *
 *  Created on: Apr 26, 2012
 *      Author: aguthals
 */

#ifndef ABRUIJNNODE_H_
#define ABRUIJNNODE_H_

#include "Node.h"
#include "AbruijnEdge.h"
#include "AssembledPeak.h"
#include "spectrum.h"
#include "SpecSet.h"

#include <vector>
#include <map>

using namespace std;
using namespace specnets;

namespace specnets
{
  class Node;
  class Edge;
}

namespace abruijn
{
  class AssembledPeak;
  class AbruijnEdge;

  class AbruijnNode : public specnets::Node
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    static AbruijnNode* castNodePtr(specnets::Node* nodePtr);

    bool m_isCorrect;
    bool m_isChimeric;
    bool m_isAnnotated;

    AbruijnNode();

    AbruijnNode(const AbruijnNode& other);

    virtual ~AbruijnNode()
    {
    }

    virtual string getGraphvizLabel(void) const;

    virtual string getGraphvizFillColor(void) const;

    void initialize(unsigned int sz = 0);

    inline void resize(unsigned int sz)
    {
      m_peaks.resize(sz);
    }

    inline void mergeColors(const AbruijnNode& other)
    {
      m_isGreen = m_isGreen || other.isGreen();
    }

    inline void loadPeak(const specnets::Spectrum& refSpec,
                         unsigned int peakIdx)
    {
      m_peaks.resize(1);
      m_peaks[0].initialize(refSpec, peakIdx);
      m_isGreen = false;
    }

    inline void addPeak(const specnets::Spectrum& refSpec, unsigned int peakIdx)
    {
      unsigned int newIdx = size();
      m_peaks.resize(newIdx + 1);
      m_peaks[newIdx].initialize(refSpec, peakIdx);
    }

    inline void loadEmptyPeak(const string& refSpecID, float prmMass)
    {
      m_peaks.resize(1);
      m_peaks[0].initializeEmpty(refSpecID, prmMass);
      m_isGreen = true;
    }

    inline void addEmptyPeak(const string& refSpecID, float prmMass)
    {
      unsigned int newIdx = size();
      m_peaks.resize(newIdx + 1);
      m_peaks[newIdx].initializeEmpty(refSpecID, prmMass);
      m_isGreen = true;
    }

    bool haveRealPeak();

    AbruijnNode &operator=(const AbruijnNode &other);

    virtual void copy(specnets::Node& otherNode);

    inline AssembledPeak &operator[](unsigned int i)
    {
      return m_peaks[i];
    }

    inline const AssembledPeak &operator[](unsigned int i) const
    {
      return m_peaks[i];
    }

    inline AssembledPeak* findSpectrum(const string& specID)
    {
      for (unsigned int i = 0; i < size(); i++)
      {
        if (m_peaks[i].getSpecID() == specID)
        {
          return &m_peaks[i];
        }
      }
      return (AssembledPeak*)0;
    }

    virtual void addBinaryVersionInfo(map<string, unsigned short>& versions) const
    {
      Node::addBinaryVersionInfo(versions);
      versions[BIN_VERSION_ID] = BIN_VERSION;
      versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;
      AssembledPeak pk;
      pk.addBinaryVersionInfo(versions);
    }

    virtual bool saveToBinaryStream(FILE* fp) const;

    virtual bool loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions);

    bool isComposite() const;

    bool isCompositeWithOther(const AbruijnNode& other) const;

    void mergeWith(const AbruijnNode& other);

    inline virtual bool isGreen() const
    {
      return m_isGreen;
    }

    inline virtual bool setGreen(const bool& isGreen)
    {
      m_isGreen = isGreen;
    }

    virtual bool isEndPoint() const;

    inline virtual unsigned int size() const
    {
      return m_peaks.size();
    }

    bool isAllYSymmetric();

    inline const string& getFPath() const
    {
      return m_fPath;
    }

    inline const string& getRPath() const
    {
      return m_rPath;
    }

    inline void setFPath(const string& fPath)
    {
      m_fPath = fPath;
    }

    inline void setRPath(const string& rPath)
    {
      m_rPath = rPath;
    }
    /*
     inline void addedEdge(AbruijnEdge* newEdge)
     {
     if (m_haveIndexedLabels)
     {
     if (newEdge->getTo() == this)
     {
     m_incomingLabels.insertWord(newEdge->getRLabelRef(), &newEdge);
     }
     else
     {
     m_outgoingLabels.insertWord(newEdge->getFLabelRef(), &newEdge);
     }
     }
     }

     inline void removingEdge(AbruijnEdge* dyingEdge)
     {
     if (m_haveIndexedLabels)
     {
     if (dyingEdge->getTo() == this)
     {
     m_incomingLabels.deleteWord(dyingEdge->getRLabelRef(), &dyingEdge);
     }
     else
     {
     m_outgoingLabels.deleteWord(dyingEdge->getFLabelRef(), &dyingEdge);
     }
     }
     }
     */

  protected:
    vector<AssembledPeak> m_peaks;
    bool m_isGreen;

    string m_fPath;
    string m_rPath;
  };
}

#endif /* ABRUIJNNODE_H_ */
