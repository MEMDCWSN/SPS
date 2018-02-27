/*
 * JumpEdge.h
 *
 *  Created on: Jan 30, 2013
 *      Author: aguthals
 */

#ifndef JUMPEDGE_H_
#define JUMPEDGE_H_

#include "Edge.h"
#include "AbruijnEdge.h"
#include "Logger.h"
#include "aminoacid.h"
#include "utils.h"
#include "mzrange.h"

#include <string>
#include <list>

using namespace std;
using namespace specnets;

namespace specnets
{
  class Node;
  class Edge;
  class AAJumps;
}

namespace abruijn
{

  class JumpEdge : public specnets::Edge
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    static JumpEdge* castEdgePtr(specnets::Edge* edgePtr);

    // Default constructor. Should always be called by derived class constructors
    JumpEdge();

    // Copy constructor. Should always be called by derived class constructors
    JumpEdge(const JumpEdge& other);

    // Initializes edges from possibly modified AA jump label
    JumpEdge(const string& jumpLabel, MZRange* jumpMass = 0);

    virtual string getGraphvizLabel(void) const;

    virtual string toString(void) const;

    JumpEdge &operator=(const JumpEdge &other);

    JumpEdge &operator=(const AbruijnEdge &other);

    // copies another edge
    virtual void copy(Edge& otherEdge);

    inline void initialize()
    {
      m_jumpCodes.clear();
      m_jumpMass.set(0, 0, 0);
      m_rWeight;
      this->setWeight(0);
    }

    /**
     * Returns:
     *   0 - both edges are equal
     *   1 - this edge supports other, but not vis-versa
     *   2 - other edge supports this, but not vis-versa
     *   3 - edges are not equal or supporting
     */
    unsigned short compareTo(const JumpEdge& other) const;

    void loadLabel(const string& jumpLabel, MZRange* jumpMass = 0);

    string getLabel() const;

    inline double getMass() const
    {
      return m_jumpMass.getMass();
    }

    inline MZRange getMassRange() const
    {
      return MZRange(m_jumpMass);
    }

    inline unsigned int getLength() const
    {
      return m_jumpCodes.size();
    }

    inline list<int>::iterator begin()
    {
      return m_jumpCodes.begin();
    }

    inline list<int>::const_iterator begin() const
    {
      return m_jumpCodes.begin();
    }

    inline list<int>::iterator end()
    {
      return m_jumpCodes.end();
    }

    inline list<int>::const_iterator end() const
    {
      return m_jumpCodes.end();
    }

    inline list<int>::reverse_iterator rbegin()
    {
      return m_jumpCodes.rbegin();
    }

    inline list<int>::const_reverse_iterator rbegin() const
    {
      return m_jumpCodes.rbegin();
    }

    inline list<int>::reverse_iterator rend()
    {
      return m_jumpCodes.rend();
    }

    inline list<int>::const_reverse_iterator rend() const
    {
      return m_jumpCodes.rend();
    }

    inline list<int>::reference front()
    {
      return m_jumpCodes.front();
    }

    inline list<int>::const_reference front() const
    {
      return m_jumpCodes.front();
    }

    inline list<int>::reference back()
    {
      return m_jumpCodes.back();
    }

    inline list<int>::const_reference back() const
    {
      return m_jumpCodes.back();
    }

    inline void reverseJumpCodes()
    {
      m_jumpCodes.reverse();
    }

    inline void addJump(const unsigned int& jumpCode)
    {
      m_jumpCodes.push_back(jumpCode);
      m_jumpMass.setMass(m_jumpMass.getMass()
          + specnets::AAJumps::getGlobalJumps()[jumpCode]);
    }

    inline void addJumpFront(const unsigned int& jumpCode)
    {
      m_jumpCodes.push_front(jumpCode);
      m_jumpMass.setMass(m_jumpMass.getMass()
          + specnets::AAJumps::getGlobalJumps()[jumpCode]);
    }

    inline unsigned int popJumpsFront(const unsigned int& numPop)
    {
      for (unsigned int i = 0; i < numPop; i++)
      {
        m_jumpMass.setMass(m_jumpMass.getMass()
            - specnets::AAJumps::getGlobalJumps()[m_jumpCodes.front()]);
        m_jumpCodes.pop_front();
      }
      return m_jumpCodes.size();
    }

    inline void setRWeight(const double& rWeight)
    {
      m_rWeight = rWeight;
    }

    inline void addRWeight(const double& addRWeight)
    {
      m_rWeight += addRWeight;
    }

    inline const double& getRWeight() const
    {
      return m_rWeight;
    }

  protected:
    list<int> m_jumpCodes;
    MZRange m_jumpMass;
    double m_rWeight;
  };
}

#endif /* JUMPEDGE_H_ */
