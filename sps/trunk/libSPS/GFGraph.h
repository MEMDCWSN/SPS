/*
 * GFGraph.h
 *
 *  Created on: Jul 22, 2013
 *      Author: aguthals
 */

#ifndef GFGRAPH_H_
#define GFGRAPH_H_

#include <vector>
#include <list>

#include "BaseGraph.h"
#include "GFTable.h"
#include "GFNode.h"
#include "GFEdge.h"

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

namespace specnets
{
  class Node;
  class Edge;
  class BaseGraph;
  class GFTable;
}

namespace abruijn
{
  class GFNode;
  class GFEdge;

  class GFGraph : public specnets::BaseGraph
  {
  public:
    // binary format version
    static const unsigned short BIN_VERSION;

    // binary format sub-version
    static const unsigned short BIN_SUBVERSION;

    // string identifier for version
    static const string BIN_VERSION_ID;

    // string identifier for sub-version
    static const string BIN_SUBVERSION_ID;

    // Upcasts a BaseGraph to a GFGraph
    static GFGraph* castGraphPtr(specnets::BaseGraph *graphPtr);

    // Default constructor
    GFGraph();

    // Construct from a GFTable
    GFGraph(const specnets::GFTable &other);

    void initialize();

    void initialize(const specnets::GFTable &other);

    /**
     * Clones a Node as an GFNode and returns a pointer to the new
     *   GFNode (as a Node pointer)
     * @param copyNode
     * @return new Node pointer
     */
    virtual specnets::Node* cloneNode(specnets::Node& copyNode) const;

    /**
     * Clones an Edge as an GFEdge and returns a pointer to the new
     *   GFEdge (as an Edge pointer)
     * @param copyEdge
     * @return new Edge pointer
     */
    virtual specnets::Edge* cloneEdge(specnets::Edge& copyEdge) const;

    /**
     * Creates a new GFNode and returns its pointer (as a Node pointer)
     * @return new Node pointer
     */
    virtual specnets::Node* createNode(void) const;

    /**
     * Creates a new GFEdge and returns its pointer (as an Edge pointer)
     * @return new Edge pointer
     */
    virtual specnets::Edge* createEdge(void) const;

    /** Intersects this graph with another. Intersection begins at the specified start masses in each graph
     *    and ends at which ever ending mass comes first. Any tailing, or non-overlapping, region that is
     *    encountered is appended to the end of the graph.
     * @param other GFGraph to intersect
     * @param myStartMass starting mass of the alignment in this graph
     * @param otherStartMass starting mass of the alignment in the other graph
     * @return false if the intersection is empty, otherwise true
     **/
    bool intersect(GFGraph &other,
                   const string &mySpecAlign,
                   const unsigned int &myStartMass,
                   const string &otherSpecAlign,
                   const unsigned int &otherStartMass);

    bool intersect(GFGraph &other, GFNode *myStartNode, GFNode *otherStartNode);

  protected:

    tr1::unordered_map<string, vector<set<GFNode*> > > m_massToNodes;
    set<GFNode*> m_endPtNodes;

  };
}

#endif /* GFGRAPH_H_ */
