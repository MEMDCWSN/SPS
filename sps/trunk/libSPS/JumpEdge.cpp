/*
 * AbruijnEdge.cpp
 *
 *  Created on: Apr 26, 2012
 *      Author: aguthals
 */

#include "JumpEdge.h"

#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace specnets;

namespace abruijn
{

  const unsigned short JumpEdge::BIN_VERSION = 1;
  const unsigned short JumpEdge::BIN_SUBVERSION = 1;

  const string JumpEdge::BIN_VERSION_ID = "JumpEdge_binVersion";
  const string JumpEdge::BIN_SUBVERSION_ID = "JumpEdge_binSubVersion";

  JumpEdge* JumpEdge::castEdgePtr(Edge* edgePtr)
  {
    JumpEdge* temp = dynamic_cast<JumpEdge*> (edgePtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast Edge \'" << edgePtr->toString()
          << "\' to a JumpEdge");

      abort();
    }
    return temp;
  }

  JumpEdge::JumpEdge(void) :
    Edge(), m_jumpCodes(), m_rWeight(0), m_jumpMass()
  {
  }

  JumpEdge::JumpEdge(const JumpEdge& other) :
    Edge(), m_jumpCodes(), m_rWeight(0), m_jumpMass()
  {
    this->operator =(other);
  }

  JumpEdge::JumpEdge(const string& jumpLabel, MZRange* jumpMass) :
    Edge(), m_jumpCodes(), m_rWeight(0), m_jumpMass()
  {
    this->loadLabel(jumpLabel, jumpMass);
  }

  string JumpEdge::getGraphvizLabel(void) const
  {
    ostringstream out;
    out << this->getLabel() << " (" << weight << "/" << m_rWeight << ")";
    return out.str();
  }

  string JumpEdge::toString(void) const
  {
    ostringstream out;
    out << Edge::toString() << " : " << getGraphvizLabel();
    return out.str();
  }

  JumpEdge & JumpEdge::operator=(const JumpEdge &other)
  {
    if (this == &other)
    {
      return *this;
    }

    Edge::operator=((const Edge&)other);
    m_jumpCodes = other.m_jumpCodes;
    m_rWeight = other.m_rWeight;
    m_jumpMass = other.m_jumpMass;

    return *this;
  }

  JumpEdge & JumpEdge::operator=(const AbruijnEdge &other)
  {
    MZRange jumpMass(other.getMass(), 0, 0.1);
    Edge::operator=((const Edge &)other);
    this->loadLabel(other.getJumpLabel(), &jumpMass);
    m_rWeight = other.getRWeight();

    return *this;
  }

  void JumpEdge::copy(Edge& otherEdge)
  {
    this->operator =(*castEdgePtr(&otherEdge));
  }

  /**
   * Returns:
   *   0 - both edges are equal
   *   1 - this edge supports other, but not vis-versa
   *   2 - other edge supports this, but not vis-versa
   *   3 - edges are not equal or supporting
   */
  unsigned short JumpEdge::compareTo(const JumpEdge& other) const
  {
    if (other.m_jumpCodes.size() == 0 && m_jumpCodes.size() == 0)
    {
      return (m_jumpMass == other.m_jumpMass) ? 0 : 3;
    }

    list<int>::const_iterator mIt = m_jumpCodes.begin();
    list<int>::const_iterator oIt = other.m_jumpCodes.begin();

    while (mIt != m_jumpCodes.end() && oIt != other.m_jumpCodes.end() && *mIt
        == *oIt)
    {
      mIt++;
      oIt++;
    }

    if (mIt != m_jumpCodes.end() && oIt != other.m_jumpCodes.end())
    {
      return 3;
    }
    else if (mIt == m_jumpCodes.end() && oIt != other.m_jumpCodes.end())
    {
      return 1;
    }
    else if (mIt != m_jumpCodes.end() && oIt == other.m_jumpCodes.end())
    {
      return 2;
    }
    else
    {
      return 0;
    }
  }

  void JumpEdge::loadLabel(const string& jumpLabel, MZRange* jumpMass)
  {
    m_jumpCodes.clear();
    if (jumpLabel.length() == 0)
    {
      if (jumpMass == 0)
      {
        ERROR_MSG("Found value for jump!!");
        abort();
      }
      m_jumpMass = *jumpMass;
      return;
    }
    vector<string> outputJumps;
    specnets::AAJumps::getSingleJumps(jumpLabel, outputJumps);
    specnets::AAJumps& globalJumps = specnets::AAJumps::getGlobalJumps();
    double totalMass = 0;

    for (unsigned int i = 0; i < outputJumps.size(); i++)
    {
      int code = globalJumps.massLookup[outputJumps[i]];
      totalMass += globalJumps[code];
      m_jumpCodes.push_back(code);
    }
    m_jumpMass.set(totalMass, 0, 0.0001);
  }

  string JumpEdge::getLabel() const
  {
    specnets::AAJumps& globalJumps = specnets::AAJumps::getGlobalJumps();
    string label = "";
    for (list<int>::const_iterator cIt = m_jumpCodes.begin(); cIt
        != m_jumpCodes.end(); cIt++)
    {
      label += globalJumps.jumpLabels[*cIt].first;
    }

    if (m_jumpCodes.size() == 0)
    {
      label = "[";
      label += parseFloat(m_jumpMass.getMass(), 2);
      label += "]";
    }
    return label;
  }
}

