/*
 * GFGraph.cpp
 *
 *  Created on: Jul 22, 2013
 *      Author: aguthals
 */

#include "GFGraph.h"
#include "aminoacid.h"
#include "Logger.h"

using namespace std;
using namespace specnets;

namespace abruijn
{
  class GFGraph;

  const unsigned short GFGraph::BIN_VERSION = 1;
  const unsigned short GFGraph::BIN_SUBVERSION = 1;

  const string GFGraph::BIN_VERSION_ID = "GFGraph_binVersion";
  const string GFGraph::BIN_SUBVERSION_ID = "GFGraph_binSubVersion";

  GFGraph* GFGraph::castGraphPtr(BaseGraph* graphPtr)
  {
    GFGraph* temp = dynamic_cast<GFGraph*> (graphPtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast BaseGraph \'" << graphPtr->toString() << "\' to a GFGraph");
      abort();
    }
    return temp;
  }

  GFGraph::GFGraph() :
    BaseGraph(), m_massToNodes(), m_endPtNodes()
  {
  }

  GFGraph::GFGraph(const GFTable &other) :
    BaseGraph(), m_massToNodes(), m_endPtNodes()
  {
    initialize(other);
  }

  void GFGraph::initialize()
  {
    m_massToNodes.rehash(0);
    m_endPtNodes.clear();
    BaseGraph::initialize(0, 0);
  }

  void GFGraph::initialize(const specnets::GFTable &other)
  {
    if (other.size() == 0)
    {
      m_massToNodes.rehash(0);
      m_endPtNodes.clear();
      BaseGraph::initialize(0, 0);
      return;
    }
    m_massToNodes.rehash(other.size());
    BaseGraph::initialize(other.size(), other.size());

    vector<set<GFNode*> > massToNodes(other.size());
    m_massToNodes[other.getUniqueID()] = massToNodes;

    vector<set<GFNode*> > &massToNodes_ref = m_massToNodes[other.getUniqueID()];

    list<pair<unsigned int, unsigned int> > nextCells;
    list<GFNode*> nextNodes;

    // Initialize the queue at mass 0
    for (unsigned int s = 0; s < other[0].size(); s++)
    {
      if (other[0][s].isValid == (char)1)
      {
        GFNode startNode(0, s);
        GFNode *addedNode =
            GFNode::castNodePtr(BaseGraph::addNode((Node &)startNode));
        massToNodes_ref[0].insert(addedNode);

        nextCells.push_back(pair<unsigned int, unsigned int> (0, s));
        nextNodes.push_back(addedNode);
      }
    }

    const vector<int> &aaMasses = other.getAAMasses();

    unsigned int curMass, curScore, m, s;
    GFNode *curNode;
    while (nextCells.size() > 0)
    {
      curMass = nextCells.front().first;
      curScore = nextCells.front().second;
      curNode = nextNodes.front();
      nextCells.pop_front();
      nextNodes.pop_front();

      for (unsigned int aa = 0; aa < aaMasses.size(); aa++)
      {
        m = curMass + aaMasses[aa];
        s = curScore + other.getScore(m);
        if (m >= other.size() || s >= other[m].size() || other[m][s].isValid
            == (char)0)
        {// skip nodes that don't exist in the GF table
          continue;
        }

        GFNode nextNode(m, s);
        GFNode *addedNode =
            GFNode::castNodePtr(BaseGraph::addNode((Node &)nextNode));
        massToNodes_ref[m].insert(addedNode);
        if (m >= other.getPmLowerBound())
        {
          m_endPtNodes.insert(addedNode);
        }

        nextCells.push_back(pair<unsigned int, unsigned int> (m, s));
        nextNodes.push_back(addedNode);

        GFEdge nextEdge(aa);
        GFEdge *addedEdge = GFEdge::castEdgePtr(BaseGraph::addEdge(curNode,
                                                                   addedNode,
                                                                   nextEdge));
      }
    }
  }

  Node* GFGraph::cloneNode(Node& copyNode) const
  {
    Node* copy = new GFNode(*(GFNode::castNodePtr(&copyNode)));
    return copy;
  }

  Edge* GFGraph::cloneEdge(Edge& copyEdge) const
  {
    Edge* copy = new GFEdge(*(GFEdge::castEdgePtr(&copyEdge)));
    return copy;
  }

  Node* GFGraph::createNode(void) const
  {
    Node* node = new GFNode();
    return node;
  }

  Edge* GFGraph::createEdge(void) const
  {
    Edge* edge = new GFEdge();
    return edge;
  }

  bool GFGraph::intersect(GFGraph &other,
                          const string &mySpecAlign,
                          const unsigned int &myStartMass,
                          const string &otherSpecAlign,
                          const unsigned int &otherStartMass)
  {
    if (m_massToNodes.count(mySpecAlign) == 0)
    {
      ERROR_MSG("Failed to locate spectrum \'" << mySpecAlign << "\' in GFGraph " << this->toString());
      abort();
    }
    if (other.m_massToNodes.count(otherSpecAlign) == 0)
    {
      ERROR_MSG("Failed to locate spectrum \'" << otherSpecAlign << "\' in GFGraph " << other.toString());
      abort();
    }
    vector<set<GFNode*> > &myNodes = m_massToNodes[mySpecAlign];
    vector<set<GFNode*> > &otherNodes = other.m_massToNodes[otherSpecAlign];

    if (myStartMass >= myNodes.size())
    {
      ERROR_MSG("Invalid start mass " << myStartMass);
      abort();
    }
    if (otherStartMass >= otherNodes.size())
    {
      ERROR_MSG("Invalid start mass " << otherStartMass);
      abort();
    }

    set<GFNode*> &myStartNodes = myNodes[myStartMass];
    set<GFNode*> &otherStartNodes = otherNodes[otherStartMass];
    bool foundIntersection = false;

    for (set<GFNode*>::iterator myIt = myStartNodes.begin(); myIt
        != myStartNodes.end(); myIt++)
    {
      for (set<GFNode*>::iterator oIt = otherStartNodes.begin(); oIt
          != otherStartNodes.end(); oIt++)
      {
        if (intersect(other, *myIt, *oIt))
        {
          foundIntersection = true;
        }
      }
    }

    return foundIntersection;
  }

  bool GFGraph::intersect(GFGraph &other,
                          GFNode *myStartNode,
                          GFNode *otherStartNode)
  {
    list<pair<GFNode*, GFNode*> > dfsStack;

    if (myStartNode->isValid() && otherStartNode->isValid())
    {
      dfsStack.push_back(pair<GFNode*, GFNode*> (myStartNode, otherStartNode));
    }

    GFNode *myNode, *otherNode, *nextNode, *nextNodeOther;
    GFEdge *nextEdge;
    pair<GFNode*, GFNode*> emptyNodePair((GFNode*)NULL, (GFNode*)NULL);
    vector<pair<GFNode*, GFNode*> > intersectEdges(0);

    /**
     * TODO: update spectrum->mass->node references with other graph's references
     */

    bool foundEndPt = false;
    while (dfsStack.size() > 0)
    {
      myNode = dfsStack.back().first;
      otherNode = dfsStack.back().second;
      dfsStack.pop_back();

      if ((!foundEndPt) && (m_endPtNodes.count(myNode) > 0
          || other.m_endPtNodes.count(otherNode) > 0))
      {
        foundEndPt = true;
      }

      // Index this graph's outgoing edges by AA jump index
      const IndexVector& myOutEdges = this->getOutEdges(myNode);

      intersectEdges.assign(intersectEdges.size(), emptyNodePair);
      for (unsigned int i = 0; i < myOutEdges.size(); i++)
      {
        if (myOutEdges[i] < 0)
        {
          continue;
        }
        nextEdge = GFEdge::castEdgePtr(this->getEdge(myOutEdges[i]));
        nextNode = GFNode::castNodePtr(this->getNode(nextEdge->getToIndex()));
        if (!nextNode->isValid())
        {
          continue;
        }
        if (nextEdge->getAA() >= intersectEdges.size())
        {
          intersectEdges.resize(nextEdge->getAA() + 1, emptyNodePair);
        }
        intersectEdges[nextEdge->getAA()].first = nextNode;
      }

      // Index the other graph's outgoing edges by AA jump index
      const IndexVector& otherOutEdges = other.getOutEdges(otherNode);
      for (unsigned int i = 0; i < otherOutEdges.size(); i++)
      {
        if (otherOutEdges[i] < 0)
        {
          continue;
        }
        nextEdge = GFEdge::castEdgePtr(other.getEdge(otherOutEdges[i]));
        nextNode = GFNode::castNodePtr(other.getNode(nextEdge->getToIndex()));
        if (!nextNode->isValid())
        {
          continue;
        }
        if (nextEdge->getAA() >= intersectEdges.size())
        {
          intersectEdges.resize(nextEdge->getAA() + 1, emptyNodePair);
        }
        intersectEdges[nextEdge->getAA()].second = nextNode;
      }

      // Invalidate nodes not in the intersection and continue exploring nodes that are
      for (unsigned int i = 0; i < intersectEdges.size(); i++)
      {
        nextNode = intersectEdges[i].first;
        nextNodeOther = intersectEdges[i].second;

        if (nextNodeOther == (GFNode*)NULL && nextNode != (GFNode*)NULL)
        {
          nextNode->invalidate();
          continue;
        }
        else if (nextNodeOther != (GFNode*)NULL && nextNode != (GFNode*)NULL)
        {
          nextNode->addScore(nextNodeOther->getScore());
          dfsStack.push_back(pair<GFNode*, GFNode*> (nextNode, nextNodeOther));
        }
        else
        {
          continue;
        }
      }
    }

    return foundEndPt;
  }
}
