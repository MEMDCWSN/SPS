/*
 * GFEdge.h
 *
 *  Created on: Jul 22, 2013
 *      Author: aguthals
 */

#ifndef GFEDGE_H_
#define GFEDGE_H_

#include "Edge.h"

using namespace std;
using namespace specnets;

namespace abruijn
{
  class GFEdge : public specnets::Edge
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    static GFEdge* castEdgePtr(specnets::Edge* edgePtr);

    GFEdge();

    GFEdge(const GFEdge& other);

    GFEdge(const int &aaIndex);

    inline const int& getAA() const
    {
      return m_aaIndex;
    }

    GFEdge & operator=(const GFEdge &other);

  protected:

    int m_aaIndex;

  };
}

#endif /* GFEDGE_H_ */
