/*
 * ContigNetwork.cpp
 *
 *  Created on: Jan 5, 2012
 *      Author: aguthals
 */

#include "ContigNetwork.h"

using namespace std;

namespace specnets
{

  ContigNetwork::ContigNetwork() :
    m_size(0), m_rootSpectra(0x0), m_edges(0x0), m_nodes(0x0), m_graph(0x0),
        m_pmTol(1.0), m_pkTol(0.5)
  {
  }

  ContigNetwork::ContigNetwork(SpecSet& _rootSpectra,
                               SpectrumPairSet& alignments,
                               float pkTol,
                               float pmTol) :
    m_size(0), m_rootSpectra(0x0), m_edges(0x0), m_nodes(0x0), m_graph(0x0),
        m_pmTol(pmTol), m_pkTol(pkTol)
  {
    initialize(_rootSpectra, alignments, pkTol, pmTol);
  }

  ContigNetwork::~ContigNetwork()
  {
    if (m_rootSpectra != 0x0)
    {
      delete m_rootSpectra;
    }

    if (m_nodes != 0x0)
    {
      eraseNodes();
      delete m_nodes;
    }

    if (m_edges != 0x0)
    {
      eraseEdges();
      delete m_edges;
    }

    if (m_graph != 0x0)
    {
      delete m_graph;
    }

    m_maxNodeIdx = 0;
    m_maxEdgeIdx = 0;
  }

  void ContigNetwork::initialize(SpecSet& _rootSpectra,
                                 SpectrumPairSet& alignments,
                                 float pkTol,
                                 float pmTol)
  {
    m_size = _rootSpectra.size();
    m_maxNodeIdx = 0;
    m_maxEdgeIdx = 0;

    if (!m_rootSpectra)
      m_rootSpectra = new SpecSet(_rootSpectra.size());

    m_rootSpectra->operator =(_rootSpectra);
    m_rootSpectra->setPeakTolerance(pkTol, false);
    m_rootSpectra->setParentMassTolerance(pmTol, false);
    if (!m_nodes)
      m_nodes = new map<int, Contig*> ;
    if (!m_graph)
      m_graph = new map<int, map<int, int> > ;
    if (!m_edges)
      m_edges = new map<int, PRMSpecEdge*> ;
    for (int i = 0; i < m_size; i++)
    {
      Contig newContig(i, m_rootSpectra);
      m_maxNodeIdx++;
      addNode(newContig);
    }

    for (int i = 0; i < alignments.size(); i++)
    {
      /*
       if (alignments[i].spec1 == 409 && alignments[i].spec2 == 1091)
       {
       DEBUG_VAR(alignments[i].spec1);
       DEBUG_VAR(alignments[i].spec2);
       DEBUG_VAR(alignments[i].shift1);
       DEBUG_VAR(alignments[i].shift2);
       DEBUG_VAR(alignments[i].score1);
       DEBUG_VAR(alignments[i].score2);
       DEBUG_VAR(alignments[i].spec2rev);

       cout << "\n";
       _rootSpectra[alignments[i].spec1].output(cout);

       cout << "\n";
       _rootSpectra[alignments[i].spec2].output(cout);
       }*/

      PRMSpecEdge newEdge(i, alignments[i], *m_rootSpectra);
      m_maxEdgeIdx++;
      addEdge(newEdge);
    }
    m_pkTol = pkTol;
    m_pmTol = pmTol;
  }

  int ContigNetwork::assembleIteratively(float minScore, int minNumMatchedPeaks)
  {

    //DEBUG_VAR(minNumMatchedPeaks);

    int numMerged = 0;

    for (map<int, Contig*>::iterator nodeIt = m_nodes->begin(); nodeIt
        != m_nodes->end(); nodeIt++)
    {
      nodeIt->second->setPeakTolerance(m_pkTol);
      nodeIt->second->setParentMassTol(m_pmTol);
    }

    PRMSpecEdge* edgeIt = getHighestScoringEdge();
    while (edgeIt != (PRMSpecEdge*)0 && edgeIt->getMinScore() > minScore)
    {

      int edgeIdx = edgeIt->index;
      int node1 = edgeIt->spec1;

      DEBUG_MSG("Merging contigs " << edgeIt->spec1 << " and " << edgeIt->spec2 << " w/ edge " << edgeIt->index);

      bool debug = false;//edgeIt->spec1 == 285 && edgeIt->spec2 == 366;

      if (!mergeEdge(edgeIdx, minNumMatchedPeaks))
      {
        WARN_MSG("Failed to merge contigs " << edgeIt->spec1 << " and " << edgeIt->spec2);
        removeEdge(edgeIdx);
        edgeIt = getHighestScoringEdge();

        if (debug)
        {
          break;
        }
        continue;
      }

      if (debug)
      {
        break;
      }

      numMerged++;

      rescoreConnectingEdges(node1);

      edgeIt = getHighestScoringEdge();
    }

    return numMerged;
  }

  void ContigNetwork::outputAbinfo(abinfo_t& outputAB)
  {
    outputAB.clear();

    for (map<int, Contig*>::iterator nodeIt = m_nodes->begin(); nodeIt
        != m_nodes->end(); nodeIt++)
    {
      if (nodeIt->second->assembledStars->size() > 0)
      {
        if (nodeIt->second->assembledStars->count(0) == 0)
        {
          ERROR_MSG("Abinfo does not have it's data at index 0!!!");
          abort();
        }
        outputAB[nodeIt->first] = (*nodeIt->second->assembledStars)[0];
      }

    }
  }

  void ContigNetwork::outputConsensusSpectra(SpecSet& outputSpecs)
  {
    outputSpecs.resize(m_size);

    for (unsigned int i = 0; i < m_size; i++)
    {

      if (m_nodes->count(i) == 0)
      {
        outputSpecs[i].resize(0);
        continue;
      }

      outputSpecs[i] = (Spectrum&)(*(*m_nodes)[i]);
    }
  }

  void ContigNetwork::outputComponents(vector<set<unsigned int> >& outputComps)
  {
    outputComps.resize(m_size);
    for (unsigned int i = 0; i < m_size; i++)
    {
      outputComps[i].clear();
      if (m_nodes->count(i) == 0)
      {
        continue;
      }
      Contig* nextContig = (*m_nodes)[i];
      for (map<int, pair<int, bool> >::iterator childIt = nextContig->childSpectra->begin(); childIt != nextContig->childSpectra->end(); childIt++) {
        outputComps[i].insert(childIt->first);
      }
    }
  }

  void ContigNetwork::outputReversed(vector<bool>& outputRev)
  {
    outputRev.resize(m_size);
    outputRev.assign(m_size, false);
    for (unsigned int i = 0; i < m_size; i++)
    {

      if (m_nodes->count(i) == 0)
      {
        continue;
      }
      Contig* ctig = (*m_nodes)[i];
      for (map<int, pair<int, bool> >::iterator childIt =
          ctig->childSpectra->begin(); childIt != ctig->childSpectra->end(); childIt++)
      {
        outputRev[childIt->first] = childIt->second.second;
      }
    }
  }

  bool ContigNetwork::mergeEdge(int edgeIdx, int minNumMatchedPeaks)
  {
    PRMSpecEdge* edge = getEdge(edgeIdx);
    if (edge == (PRMSpecEdge*)0)
    {
      return false;
    }

    int node1 = edge->spec1;
    int node2 = edge->spec2;

    if (edge->spec2rev)
    {
      reverseNode(node2);
    }
    if (edge->spec2rev)
    {
      ERROR_MSG("Contig " << node2 << " is still reversed!!");
      abort();
    }

    Contig* mergedContig = getNode(node1);
    Contig* contig2 = getNode(node2);

    if (!mergedContig->merge(edge,
                             contig2,
                             minNumMatchedPeaks,
                             m_pkTol,
                             m_pmTol))
    {
      return false;
    }
    mergedContig->setPeakTolerance(m_pkTol);
    mergedContig->setParentMassTol(m_pmTol);

    list<PRMSpecEdge> edgesToAdd;
    list<int> edgesToRemove;

    map<int, int>* node2Neighbors = &(*m_graph)[node2];
    for (map<int, int>::iterator neighIt = node2Neighbors->begin(); neighIt
        != node2Neighbors->end(); neighIt++)
    {
      int node3 = neighIt->first;

      if (node3 == node1)
      {
        continue;
      }

      int otherEdge3Idx = neighIt->second;
      PRMSpecEdge* otherEdge3Ptr = getEdge(otherEdge3Idx);
      PRMSpecEdge otherEdge3;

      PRMSpecEdge::appendEdges(edge, otherEdge3Ptr, &otherEdge3);

      if (!containsEdge(node1, node3))
      {
        edgesToAdd.push_back(otherEdge3);
      }
      else
      {
        PRMSpecEdge* origEdge = getEdge(node1, node3);

        pair<float, float> scoreOrigRes = getEdgeScore(*origEdge);
        pair<float, float> scoreOther3Res = getEdgeScore(otherEdge3);
        float scoreOrig = min(scoreOrigRes.first, scoreOrigRes.second);
        float scoreOther3 = min(scoreOther3Res.first, scoreOther3Res.second);

        if (scoreOrig < scoreOther3)
        {
          edgesToAdd.push_back(otherEdge3);
          edgesToRemove.push_back(origEdge->index);
        }
      }
    }

    contig2 = (Contig*)0;
    edge = 0;
    removeNode(node2);

    for (list<int>::iterator removeIt = edgesToRemove.begin(); removeIt
        != edgesToRemove.end(); removeIt++)
    {
      removeEdge(*removeIt);
    }

    for (list<PRMSpecEdge>::iterator addIt = edgesToAdd.begin(); addIt
        != edgesToAdd.end(); addIt++)
    {
      addIt->index = m_maxEdgeIdx;
      m_maxEdgeIdx++;
      addEdge(*addIt);
    }

    return true;
  }

  bool ContigNetwork::rescoreConnectingEdges(int nodeIdx)
  {
    if (!containsNode(nodeIdx))
    {
      ERROR_MSG("Node " << nodeIdx << " does not exist");
      return false;
    }

    map<int, int>* neighborRef = &(*m_graph)[nodeIdx];
    for (map<int, int>::iterator neighIt = neighborRef->begin(); neighIt
        != neighborRef->end(); neighIt++)
    {
      PRMSpecEdge* edge = getEdge(neighIt->second);
      pair<float, float> scoreOrigRes = getEdgeScore(*edge);
      edge->score1 = scoreOrigRes.first;
      edge->score2 = scoreOrigRes.second;
    }

    return true;
  }

  list<pair<int, int> >* ContigNetwork::getNeighborList(int metaContigIdx) const
  {
    list<pair<int, int> >* neighbors = new list<pair<int, int> > ;

    if (!containsNode(metaContigIdx))
    {
      return neighbors;
    }
    map<int, int>* neighborRef = &(*m_graph)[metaContigIdx];
    for (map<int, int>::iterator neighIt = neighborRef->begin(); neighIt
        != neighborRef->end(); neighIt++)
    {
      neighbors->push_back(pair<int, int> (neighIt->first, neighIt->second));
    }
    return neighbors;
  }

  bool ContigNetwork::containsNode(int metaContigIdx) const
  {
    return m_graph->count(metaContigIdx) > 0;
  }

  bool ContigNetwork::containsEdge(int node1, int node2) const
  {
    if ((!containsNode(node1)) || (!containsNode(node2))
        || (*m_graph)[node1].count(node2) == 0)
    {
      return false;
    }
    return true;
    /*

     if (m_edges->count((*m_graph)[node1][node2]) == 0) {
     ERROR_MSG("Graph is inconsistent");
     return false;
     }
     return true;
     */
  }

  bool ContigNetwork::containsEdge(int edgeIdx) const
  {
    if (m_edges->count(edgeIdx) == 0)
    {
      return false;
    }
    return true;
    /*
     PRMSpecEdge* edge = (*m_edges)[edgeIdx];

     if ((!containsNode(edge->spec1)) || (!containsNode(edge->spec2))
     || (*m_graph)[edge->spec1].count(edge->spec2) == 0) {
     ERROR_MSG("Graph is inconsistent");
     return false;
     }
     return true;
     */
  }

  Contig* ContigNetwork::getNode(int index) const
  {
    if (!containsNode(index))
    {
      ERROR_MSG("Node " << index << " does not exist");
      return (Contig*)0;
    }
    else
    {
      return (*m_nodes)[index];
    }
  }

  PRMSpecEdge* ContigNetwork::getEdge(int index) const
  {
    if (!containsEdge(index))
    {
      ERROR_MSG("Edge " << index << " does not exist");
      return (PRMSpecEdge*)0;
    }
    else
    {
      return (*m_edges)[index];
    }
  }

  PRMSpecEdge* ContigNetwork::getEdge(int node1, int node2) const
  {
    if (!containsEdge(node1, node2))
    {
      ERROR_MSG("Edge between nodes " << node1 << " and " << node2 << " does not exist");
      return (PRMSpecEdge*)0;
    }
    else
    {
      return (*m_edges)[(*m_graph)[node1][node2]];
    }
  }

  bool ContigNetwork::addNode(const Contig& newNode)
  {
    if (containsNode(newNode.index))
    {
      ERROR_MSG("Node " << newNode.index << " already exists");
      return false;
    }

    (*m_nodes)[newNode.index] = new Contig(newNode);
    //DEBUG_VAR((*m_nodes)[newNode.index]->childSpectra->size());
    map<int, int> subGraph;
    (*m_graph)[newNode.index] = subGraph;
    return true;
  }

  bool ContigNetwork::addEdge(const PRMSpecEdge& newEdge)
  {
    if (containsEdge(newEdge.spec1, newEdge.spec2))
    {
      ERROR_MSG("Edge between m_nodes " << newEdge.spec1 << " and " << newEdge.spec2 << " already exists");
      return false;
    }
    else if (containsEdge(newEdge.index))
    {
      ERROR_MSG("Edge at index " << newEdge.index << " already exists");
      return false;
    }

    addEdge(newEdge.spec1, newEdge.spec2, newEdge.index);

    (*m_edges)[newEdge.index] = new PRMSpecEdge(newEdge);

    return true;
  }

  bool ContigNetwork::removeNode(int nodeIdx)
  {

    if (!containsNode(nodeIdx))
    {
      ERROR_MSG("Node " << nodeIdx << " does not exist");
      return false;
    }

    map<int, int>* neighborRef = &(*m_graph)[nodeIdx];
    for (map<int, int>::iterator neighIt = neighborRef->begin(); neighIt
        != neighborRef->end(); neighIt++)
    {
      //removeEdge(nodeIdx, neighIt->first);

      int edgeIdx = neighIt->second;
      int node2 = neighIt->first;
      delete (*m_edges)[edgeIdx];

      m_edges->erase(edgeIdx);
      (*m_graph)[node2].erase(nodeIdx);
    }

    m_graph->erase(nodeIdx);

    delete (*m_nodes)[nodeIdx];

    m_nodes->erase(nodeIdx);

    return true;
  }

  bool ContigNetwork::removeEdge(int edgeIdx)
  {
    PRMSpecEdge* edgeToRemove = getEdge(edgeIdx);

    if (edgeToRemove == 0)
    {
      return false;
    }

    int node1 = edgeToRemove->spec1;
    int node2 = edgeToRemove->spec2;

    edgeToRemove = (PRMSpecEdge*)0;

    (*m_graph)[node1].erase(node2);
    (*m_graph)[node2].erase(node1);

    delete (*m_edges)[edgeIdx];

    m_edges->erase(edgeIdx);

    return true;
  }

  bool ContigNetwork::removeEdge(int node1, int node2)
  {
    PRMSpecEdge* edgeToRemove = getEdge(node1, node2);

    if (edgeToRemove == 0)
    {
      return false;
    }
    else
    {
      int edgeIdx = edgeToRemove->index;
      edgeToRemove = (PRMSpecEdge*)0;
      return removeEdge(edgeIdx);
    }
  }

  bool ContigNetwork::reverseNode(int metaContigIdx)
  {
    DEBUG_TRACE;
    if (!containsNode(metaContigIdx))
    {
      return false;
    }
    DEBUG_TRACE;
    (*m_nodes)[metaContigIdx]->reverse();
    DEBUG_TRACE;
    map<int, int>* neighborRef = &(*m_graph)[metaContigIdx];
    for (map<int, int>::iterator neighIt = neighborRef->begin(); neighIt
        != neighborRef->end(); neighIt++)
    {
      int edgeIdx = neighIt->second;
      PRMSpecEdge* connectingEdge = (*m_edges)[edgeIdx];
      connectingEdge->reverse(metaContigIdx);
    }
    DEBUG_TRACE;
    return true;
  }

  float ContigNetwork::getConsensusShift(PRMSpecEdge& edge,
                                         int nodeFrom,
                                         int nodeTo) const
  {

    Contig* node1 = getNode(nodeFrom);
    Contig* node2 = getNode(nodeTo);
    PRMSpecEdge* edgeUse = &edge;

    Contig node2Copy;
    PRMSpecEdge edgeCopy;

    if (edge.spec2rev)
    {
      node2Copy = *node2;
      node2Copy.reverse();
      node2 = &node2Copy;

      edgeCopy = edge;
      edgeCopy.reverse(nodeTo);
      edgeUse = &edgeCopy;
    }

    float totalScoreF, totalScoreR;
    float leftEdgeF = 0, rightEdgeF = 0;
    for (map<int, pair<float, float> >::iterator refIt =
        node1->rootRef->begin(); refIt != node1->rootRef->end(); refIt++)
    {
      if (refIt->first == edgeUse->spec1)
      {
        continue;
      }
      leftEdgeF = min(leftEdgeF, refIt->second.first);
    }
    for (map<int, pair<float, float> >::iterator refIt =
        node2->rootRef->begin(); refIt != node2->rootRef->end(); refIt++)
    {
      if (refIt->first == edgeUse->spec2)
      {
        continue;
      }
      rightEdgeF = min(rightEdgeF, refIt->second.first);
    }
    /*
     if (debug) {
     cout << "Root shift between " << idx1 << " and " << idx2 << " = "
     << fShift << ", " << rShift << "\n";
     cout << "consensusFShift = " << fShift << " - " << leftEdgeF << " + "
     << rightEdgeF << " = " << fShift - leftEdgeF + rightEdgeF
     << "\n";
     cout << "consensusRShift = " << rShift << " - " << leftEdgeF << " + "
     << rightEdgeR << " = " << rShift - leftEdgeF + rightEdgeR
     << "\n";
     }*/

    if (edgeUse->spec2 == nodeTo)
    {
      return edgeUse->getShift(nodeFrom, nodeTo) - leftEdgeF + rightEdgeF
          - node1->endGaps.first + node2->endGaps.first;
    }
    else
    {
      return edgeUse->getShift(nodeFrom, nodeTo) + leftEdgeF - rightEdgeF
          + node1->endGaps.first - node2->endGaps.first;
    }
  }

  pair<float, float> ContigNetwork::getEdgeScore(PRMSpecEdge& edge)
  {

    float shift = getConsensusShift(edge, edge.spec1, edge.spec2);

    Contig* node1 = getNode(edge.spec1);
    Contig* node2 = getNode(edge.spec2);

    Contig node2Copy;

    if (edge.spec2rev)
    {
      node2Copy = *node2;
      node2Copy.reverse();
      node2 = &node2Copy;
    }

    m_alignmentObj.setSpec1((Spectrum*)node1);
    m_alignmentObj.setSpec2((Spectrum*)node2);

    pair<int, pair<float, float> > alignRes =
        m_alignmentObj.getShiftScore(shift, m_pkTol, 0);

    return pair<float, float> (alignRes.second.first, alignRes.second.second);
  }

  PRMSpecEdge* ContigNetwork::getHighestScoringEdge()
  {

    if (m_edges->size() == 0)
    {
      return (PRMSpecEdge*)0;
    }

    map<int, PRMSpecEdge*>::iterator edgeIt = m_edges->begin();

    PRMSpecEdge* bestEdge = edgeIt->second;

    edgeIt++;

    for (; edgeIt != m_edges->end(); edgeIt++)
    {
      bestEdge = (edgeIt->second->getMinScore() > bestEdge->getMinScore())
          ? edgeIt->second : bestEdge;
    }
    return bestEdge;
  }

  void ContigNetwork::addEdge(int spec1, int spec2, int edgeIdx)
  {

    map<int, int> emptySubGraph;
    if (m_graph->count(spec1) == 0)
    {
      (*m_graph)[spec1] = emptySubGraph;
      (*m_graph)[spec1][spec2] = edgeIdx;
    }
    else
    {
      if ((*m_graph)[spec1].count(spec2) > 0)
      {
        WARN_MSG("Over-writting edge between " << spec1 << " and " << spec2 << "(edge index from " << (*m_graph)[spec1][spec2] << " to " << edgeIdx << ")");
      }
      (*m_graph)[spec1][spec2] = edgeIdx;
    }

    if (m_graph->count(spec2) == 0)
    {
      (*m_graph)[spec2] = emptySubGraph;
      (*m_graph)[spec2][spec1] = edgeIdx;
    }
    else
    {
      if ((*m_graph)[spec2].count(spec1) > 0)
      {
        WARN_MSG("Over-writting edge between " << spec2 << " and " << spec1 << "(edge index from " << (*m_graph)[spec2][spec1] << " to " << edgeIdx << ")");
      }
      (*m_graph)[spec2][spec1] = edgeIdx;
    }
  }

  void ContigNetwork::eraseNodes()
  {
    for (map<int, Contig*>::iterator nodeIt = m_nodes->begin(); nodeIt
        != m_nodes->end(); nodeIt++)
    {
      delete nodeIt->second;
      nodeIt->second = 0;
    }
  }

  void ContigNetwork::eraseEdges()
  {
    for (map<int, PRMSpecEdge*>::iterator eIt = m_edges->begin(); eIt
        != m_edges->end(); eIt++)
    {
      delete eIt->second;
      eIt->second = 0;
    }
  }

  int ContigNetwork::getNewNodeIndex()
  {
    int idx = m_maxNodeIdx;
    m_maxNodeIdx++;
    return idx;
  }

  int ContigNetwork::getNewEdgeIndex()
  {
    int idx = m_maxEdgeIdx;
    m_maxEdgeIdx++;
    return idx;
  }

}
