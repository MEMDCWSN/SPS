/*
 * GFNode.h
 *
 *  Created on: Jul 22, 2013
 *      Author: aguthals
 */

#ifndef GFNODE_H_
#define GFNODE_H_

#include "Node.h"

using namespace std;
using namespace specnets;

namespace specnets
{
  class Node;
  class Edge;
}

namespace abruijn
{

  class GFNode : public specnets::Node
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    static GFNode* castNodePtr(specnets::Node *nodePtr);

    GFNode();

    GFNode(const unsigned int &mass, const unsigned int &score);

    GFNode(const GFNode &other);

    virtual ~GFNode()
    {
    }

    virtual string getGraphvizLabel(void) const;

    GFNode & operator=(const GFNode &other);

    inline void invalidate()
    {
      m_isValid = false;
    }

    inline const bool& isValid() const
    {
      return m_isValid;
    }

    inline void addScore(const unsigned int &newScore)
    {
      m_score += newScore;
    }

    inline const unsigned int& getScore() const
    {
      return m_score;
    }

  protected:
    unsigned int m_mass;
    unsigned int m_score;
    bool m_isValid;

  };
}

#endif /* GFNODE_H_ */
