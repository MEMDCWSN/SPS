/*
 * GFNode.cpp
 *
 *  Created on: Jul 22, 2013
 *      Author: aguthals
 */

#include "GFNode.h"
#include <stdlib.h>

using namespace std;
using namespace specnets;

namespace abruijn
{

  const unsigned short GFNode::BIN_VERSION = 1;
  const unsigned short GFNode::BIN_SUBVERSION = 1;

  const string GFNode::BIN_VERSION_ID = "GFNode_binVersion";
  const string GFNode::BIN_SUBVERSION_ID = "GFNode_binSubVersion";

  GFNode* GFNode::castNodePtr(Node* nodePtr)
  {
    GFNode* temp = dynamic_cast<GFNode*> (nodePtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast Node \'" << nodePtr->toString() << "\' to a GFNode");
      abort();
    }
    return temp;
  }

  GFNode::GFNode(void) :
    Node(), m_mass(0), m_score(0), m_isValid(true)
  {
  }

  GFNode::GFNode(const unsigned int &mass, const unsigned int &score) :
    m_mass(mass), m_score(score), m_isValid(true)
  {
  }

  GFNode::GFNode(const GFNode &other) :
    Node((const Node &)other), m_mass(other.m_mass), m_score(other.m_score),
        m_isValid(other.m_isValid)
  {
  }

  string GFNode::getGraphvizLabel(void) const
  {
    ostringstream out;
    out << Node::getGraphvizLabel();
    out << " : (" << m_mass << "," << m_score << ")";
    return out.str();
  }

  GFNode & GFNode::operator=(const GFNode &other)
  {
    if (this == &other)
    {
      return *this;
    }

    Node::operator=((const Node&)other);
    m_mass = other.m_mass;
    m_score = other.m_score;
    m_isValid = other.m_isValid;

    return *this;
  }

}
