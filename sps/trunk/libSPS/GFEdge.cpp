/*
 * GFEdge.cpp
 *
 *  Created on: Jul 22, 2013
 *      Author: aguthals
 */

#include "GFEdge.h"
#include <stdlib.h>

using namespace std;
using namespace specnets;

namespace abruijn
{

  const unsigned short GFEdge::BIN_VERSION = 1;
  const unsigned short GFEdge::BIN_SUBVERSION = 1;

  const string GFEdge::BIN_VERSION_ID = "GFEdge_binVersion";
  const string GFEdge::BIN_SUBVERSION_ID = "GFEdge_binSubVersion";

  GFEdge* GFEdge::castEdgePtr(Edge* edgePtr)
  {
    GFEdge* temp = dynamic_cast<GFEdge*> (edgePtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast Edge \'" << edgePtr->toString() << "\' to a GFEdge");

      abort();
    }
    return temp;
  }

  GFEdge::GFEdge() :
    Edge(), m_aaIndex(-1)
  {
  }

  GFEdge::GFEdge(const GFEdge &other) :
    Edge((const Edge&)other), m_aaIndex(other.m_aaIndex)
  {
  }

  GFEdge::GFEdge(const int &aaIndex) :
    Edge(), m_aaIndex(aaIndex)
  {
  }

  GFEdge & GFEdge::operator=(const GFEdge &other)
  {

    if (this == &other)
    {
      return *this;
    }

    Edge::operator=((const Edge&)other);
    m_aaIndex = other.m_aaIndex;

    return *this;
  }

}
