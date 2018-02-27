/*
 * AbruijnEdge.h
 *
 *  Created on: Apr 26, 2012
 *      Author: aguthals
 */

#ifndef ABRUIJNEDGE_H_
#define ABRUIJNEDGE_H_

#include "Edge.h"
#include "AbruijnNode.h"
#include "Logger.h"
#include "aminoacid.h"
#include "utils.h"

#include <string>
#include <list>
#include <set>

using namespace std;
using namespace specnets;

namespace specnets
{
  class Node;
  class Edge;
}

namespace abruijn
{
  class AAJumps;
  class AbruijnEdge;
  class AbruijnNode;

  class AbruijnEdge : public specnets::Edge
  {
  public:

    class MergeSizeCompare
    {
      bool m_lookReverse;
    public:
      MergeSizeCompare(const bool& lookReverse) :
        m_lookReverse(lookReverse)
      {
      }

      bool operator()(const AbruijnEdge* lhs, const AbruijnEdge* rhs) const
      {
        const string& lPath = (m_lookReverse) ? lhs->getRLabel()
            : lhs->getFLabel();
        const string& rPath = (m_lookReverse) ? rhs->getRLabel()
            : rhs->getFLabel();
        if (lPath.length() == rPath.length())
        {
          return lhs->getSpecIDs().size() >= rhs->getSpecIDs().size();
        }
        return lPath.length() >= rPath.length();
      }
    };

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    static AbruijnEdge* castEdgePtr(specnets::Edge* edgePtr);

    bool m_isCorrect;
    bool m_isAnnotated;

    // Default constructor. Should always be called by derived class constructors
    AbruijnEdge(void);

    // Copy constructor. Should always be called by derived class constructors
    AbruijnEdge(const AbruijnEdge& other);

    AbruijnEdge(const list<AbruijnEdge*>& edgeList, const bool& lookReverse);

    virtual string getGraphvizLabel(void) const;

    virtual string toString(void) const;

    AbruijnEdge &operator=(const AbruijnEdge &other);

    // copies another edge
    virtual void copy(Edge& otherEdge);

    unsigned int mergeEdges(const list<AbruijnEdge*>& edgeList,
                            const bool& lookReverse,
                            list<AbruijnEdge*>* usedEdges = 0);

    virtual void addBinaryVersionInfo(map<string, unsigned short>& versions) const
    {
      Edge::addBinaryVersionInfo(versions);
      versions[BIN_VERSION_ID] = BIN_VERSION;
      versions[BIN_SUBVERSION_ID] = BIN_SUBVERSION;
    }

    virtual bool saveToBinaryStream(FILE* fp) const;

    virtual bool loadFromBinaryStream(FILE* fp,
                                      map<string, unsigned short>& versions);

    inline const double& getMass() const
    {
      return m_mass;
    }

    inline const string& getFLabel() const
    {
      return m_fLabel;
    }

    inline const string& getJumpLabel() const
    {
      return m_jumpLabel;
    }

    inline const string& getRLabel() const
    {
      return m_rLabel;
    }

    inline const list<string>& getSpecIDs() const
    {
      return m_specIDs;
    }

    inline const double& getRWeight() const
    {
      return m_rWeight;
    }

    inline const bool& expandOnly() const
    {
      return m_expandOnly;
    }

    inline void setExpandOnly(const bool& expandOnly)
    {
      m_expandOnly = expandOnly;
    }

    template<class InputIterator> void loadIDs(InputIterator start,
                                               InputIterator end)
    {
      m_specIDs.clear();
      m_specIDs.insert(m_specIDs.end(), start, end);
    }

    inline void addID(const string& specID)
    {
      m_specIDs.push_back(specID);
    }

    inline const bool& isModified() const
    {
      return m_modified;
    }

    inline void setMass(double mass)
    {
      m_mass = mass;
    }

    inline void setFLabel(const string& fLabel)
    {
      m_fLabel = fLabel;
    }

    inline void setJumpLabel(const string& jumpLabel)
    {
      m_jumpLabel = jumpLabel;
    }

    inline void setRLabel(const string& rLabel)
    {
      m_rLabel = rLabel;
    }

    inline void setModified(const bool& modified)
    {
      m_modified = modified;
    }

    inline void setRWeight(const double& rWeight)
    {
      m_rWeight = rWeight;
    }

    inline void set(const double& mass,
                    const double& _weight,
                    const double& rWeight,
                    const string& jumpLabel,
                    const string& fLabel,
                    const string& rLabel,
                    const bool& modified,
                    const list<string>& specIDs)
    {
      setWeight(_weight);
      m_rWeight = rWeight;
      m_mass = mass;
      m_jumpLabel = jumpLabel;
      m_fLabel = fLabel;
      m_rLabel = rLabel;
      m_modified = modified;
      m_specIDs = specIDs;
    }

  protected:
    double m_mass;
    double m_rWeight;
    string m_jumpLabel;
    string m_fLabel;
    string m_rLabel;
    bool m_modified;
    list<string> m_specIDs;
    bool m_expandOnly;
  };
}

#endif /* ABRUIJNEDGE_H_ */
