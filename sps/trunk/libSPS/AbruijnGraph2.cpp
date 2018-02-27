/*
 * AbruijnGraph2.cpp
 *
 *  Created on: May 1, 2012
 *      Author: aguthals
 */

#include "AbruijnGraph2.h"
#include "TrieMap.h"
#include "UnionFind.h"

#include <limits>
#include <algorithm>

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

namespace abruijn
{
  class AbruijnGraph2;

  AbruijnGraph2::ParallelEdge::ParallelEdge() :
    m_jumpCode(-1), m_edgeIds(0)
  {
  }

  AbruijnGraph2::ParallelEdge::ParallelEdge(const JumpEdge& newEdge) :
    m_jumpCode(-1), m_edgeIds(1)
  {
    this->operator =(newEdge);
  }

  AbruijnGraph2::ParallelEdge::ParallelEdge(const ParallelEdge& other) :
    m_jumpCode(other.m_jumpCode), m_edgeIds(other.m_edgeIds)
  {
  }

  AbruijnGraph2::ParallelEdge& AbruijnGraph2::ParallelEdge::operator=(const JumpEdge &other)
  {
    if (other.getLength() != 1)
    {
      ERROR_MSG("Parallel edges can only have length 1");
      abort();
    }
    m_edgeIds.resize(1);
    m_edgeIds[0] = other.getIndex();
    m_jumpCode = other.front();
    return *this;
  }

  AbruijnGraph2::ParallelEdge& AbruijnGraph2::ParallelEdge::operator=(const AbruijnGraph2::ParallelEdge &other)
  {
    if (this == &other)
    {
      return *this;
    }

    m_edgeIds = other.m_edgeIds;
    m_jumpCode = other.m_jumpCode;
    return *this;
  }

  AbruijnGraph2::ParallelPath::ParallelPath() :
    m_path(0), m_nodes(0), m_outgoingEdges(0), m_pathNodes(), m_fWeight(0),
        m_rWeight(0), m_outgoingEdgeIDs()
  {
  }

  AbruijnGraph2::ParallelPath::ParallelPath(const ParallelPath& other) :
    m_path(other.m_path), m_nodes(other.m_nodes),
        m_outgoingEdges(other.m_outgoingEdges), m_pathNodes(other.m_pathNodes),
        m_fWeight(other.m_fWeight), m_rWeight(other.m_rWeight),
        m_outgoingEdgeIDs(other.m_outgoingEdgeIDs)
  {

  }

  AbruijnGraph2::ParallelPath::ParallelPath(const AbruijnGraph2& abG,
                                            const AbruijnNode& start) :
    m_path(0), m_nodes(0), m_outgoingEdges(0), m_pathNodes(), m_fWeight(0),
        m_rWeight(0), m_outgoingEdgeIDs()
  {
    initialize(abG, start);
  }

  void AbruijnGraph2::ParallelPath::initialize(const AbruijnGraph2& abG,
                                               const AbruijnNode& start)
  {
    clear();
    m_nodes.resize(1);
    m_nodes[0] = start;

    const IndexVector& outEdges = abG.getOutEdges(start.getIndex());
    for (unsigned long i = 0; i < outEdges.size(); i++)
    {
      if (outEdges[i] < 0)
      {
        continue;
      }
      JumpEdge* abEdge = JumpEdge::castEdgePtr(abG.getEdge(outEdges[i]));

      m_outgoingEdges.push_back(*abEdge);
    }
  }

  AbruijnGraph2::ParallelPath &AbruijnGraph2::ParallelPath::operator=(const ParallelPath &other)
  {
    if (this == &other)
    {
      return *this;
    }
    m_path = other.m_path;
    m_nodes = other.m_nodes;
    m_outgoingEdges = other.m_outgoingEdges;
    m_outgoingEdgeIDs = other.m_outgoingEdgeIDs;
    m_pathNodes = other.m_pathNodes;
    m_fWeight = other.m_fWeight;
    m_rWeight = other.m_rWeight;
    return *this;
  }

  void AbruijnGraph2::ParallelPath::advance(const AbruijnGraph2& abG,
                                            const JumpEdge& advanceEdge)
  {
    if (AbruijnGraph2::DEBUG_EXPANSION)
    {
      DEBUG_MSG("Advancing path " << toString() << " by " << advanceEdge.toString());
      for (unsigned int i = 0; i < m_outgoingEdges.size(); i++)
      {
        DEBUG_VAR(m_outgoingEdges[i].toString());
      }
    }

    if (advanceEdge.getLength() == 0)
    {
      m_path.resize(m_path.size() + 1);
      m_path[m_path.size() - 1] = advanceEdge;
      m_nodes.resize(m_nodes.size() + 1);
      // add mass edge
    }
    else
    {
      unsigned int idxUse = m_path.size();
      m_path.resize(m_path.size() + advanceEdge.getLength());
      for (list<int>::const_iterator cIt = advanceEdge.begin(); cIt
          != advanceEdge.end(); cIt++)
      {
        m_path[idxUse++].m_jumpCode = *cIt;
      }
      m_nodes.resize(m_nodes.size() + advanceEdge.getLength());
    }

    vector<JumpEdge> oldOutgoingEdges(m_outgoingEdges);
    m_outgoingEdges.resize(0);
    m_outgoingEdgeIDs.clear();

    if (AbruijnGraph2::DEBUG_EXPANSION)
    {
      DEBUG_MSG("New path: " << toString());
    }

    for (vector<JumpEdge>::const_iterator eIt = oldOutgoingEdges.begin(); eIt
        != oldOutgoingEdges.end(); eIt++)
    {
      if (m_pathNodes.count(eIt->getToNodeIndex()) > 0)
      {
        continue;
      }
      if (advanceEdge.getLength() == 0)
      {
        if (eIt->compareTo(advanceEdge) != 0)
        {
          continue;
        }

        if (!addNode(abG, m_nodes.size() - 1, eIt->getToNodeIndex()))
          continue;

        addEdge(m_path.size() - 1, *eIt);
      }

      if (eIt->getLength() == 0)
      {
        continue;
      }
      const unsigned short comparVal = eIt->compareTo(advanceEdge);

      if (comparVal >= 3)
      {
        continue;
      }
      else if (comparVal == 2)
      {
        if (m_outgoingEdgeIDs.count(eIt->getIndex()) == 0)
        {
          m_outgoingEdgeIDs.insert(eIt->getIndex());
          JumpEdge tempEdge2(*eIt);
          tempEdge2.popJumpsFront(advanceEdge.getLength());
          m_outgoingEdges.push_back(tempEdge2);
        }
        continue;
      }
      else
      {
        if (!addNode(abG, m_nodes.size() - advanceEdge.getLength()
            + eIt->getLength() - 1, eIt->getToNodeIndex()))
          continue;

        addEdge(m_path.size() - advanceEdge.getLength() + eIt->getLength() - 1,
                *eIt);
      }
    }

    if (AbruijnGraph2::DEBUG_EXPANSION)
    {
      for (unsigned int i = 0; i < m_outgoingEdges.size(); i++)
      {
        DEBUG_VAR(m_outgoingEdges[i].toString());
      }
    }
  }

  bool AbruijnGraph2::ParallelPath::operator==(const ParallelPath& other) const
  {
    if (m_path != other.m_path)
    {
      return false;
    }

    if (!std::includes(m_pathNodes.begin(),
                       m_pathNodes.end(),
                       other.m_pathNodes.begin(),
                       other.m_pathNodes.end()))
    {
      return false;
    }

    return true;
  }

  bool AbruijnGraph2::ParallelPath::addNode(const AbruijnGraph2& abG,
                                            const unsigned int& pathIdx,
                                            const long& node)
  {
    AbruijnNode& abNode = m_nodes.at(pathIdx);
    AbruijnNode* other = abG.getNode(node);

    if (abNode.isCompositeWithOther(*other))
    {
      return false;
    }

    if (AbruijnGraph2::DEBUG_EXPANSION)
    {
      DEBUG_MSG("Adding node " << node << " at index " << pathIdx);
    }

    abNode.mergeWith(*other);
    m_pathNodes.insert(node);

    list<pair<long, unsigned int> > backBFS;
    backBFS.push_back(pair<long, unsigned int> (node, pathIdx));

    while (backBFS.size() > 0)
    {
      long curNode = backBFS.front().first;
      unsigned int curIdx = backBFS.front().second;
      backBFS.pop_front();

      if (curIdx > 0)
      {
        JumpEdge backPath;
        if (m_path[curIdx - 1].m_jumpCode < 0)
        {
          MZRange massR =
              abG.getEdge(m_path[curIdx - 1].m_edgeIds.front())->getMassRange();
          backPath.loadLabel("", &massR);
        }
        else
        {
          for (unsigned int i = 0; i < curIdx; i++)
          {
            backPath.addJumpFront(m_path[i].m_jumpCode);
          }
        }

        const IndexVector& inEdges = abG.getInEdges(curNode);
        for (unsigned int i = 0; i < inEdges.size(); i++)
        {
          if (inEdges[i] < 0)
          {
            continue;
          }
          JumpEdge* abEdge = abG.getEdge(inEdges[i]);
          JumpEdge revEdge(*abEdge);
          revEdge.reverseJumpCodes();

          if (revEdge.compareTo(backPath) > 1)
          {
            continue;
          }
          addEdge(curIdx - 1, *abEdge);

          unsigned int jLength = (revEdge.getLength() == 0) ? 1
              : revEdge.getLength();
          AbruijnNode& myAbNode = m_nodes.at(curIdx - jLength);
          AbruijnNode* otherNode = abG.getNode(revEdge.getFromNodeIndex());

          if (m_pathNodes.count(revEdge.getFromNodeIndex()) > 0
              || myAbNode.isCompositeWithOther(*otherNode))
          {
            continue;
          }

          myAbNode.mergeWith(*otherNode);
          m_pathNodes.insert(otherNode->getIndex());

          backBFS.push_back(pair<long, unsigned int> (revEdge.getFromNodeIndex(),
                                                      curIdx - jLength));
        }
      }

      JumpEdge forPath;
      if (curIdx < m_path.size())
      {
        if (m_path[curIdx].m_jumpCode < 0)
        {
          MZRange massR =
              abG.getEdge(m_path[curIdx].m_edgeIds.front())->getMassRange();
          forPath.loadLabel("", &massR);
        }
        else
        {
          for (unsigned int i = curIdx; i < m_path.size(); i++)
          {
            forPath.addJump(m_path[i].m_jumpCode);
          }
        }
      }

      if (AbruijnGraph2::DEBUG_EXPANSION)
      {
        DEBUG_VAR(forPath.toString());
      }

      const IndexVector& outEdges = abG.getOutEdges(curNode);
      for (unsigned long i = 0; i < outEdges.size(); i++)
      {
        if (outEdges[i] < 0)
        {
          continue;
        }
        JumpEdge* abEdge = JumpEdge::castEdgePtr(abG.getEdge(outEdges[i]));
        const unsigned short nextCompareVal = abEdge->compareTo(forPath);

        if (AbruijnGraph2::DEBUG_EXPANSION)
        {
          DEBUG_MSG("Considering (" << nextCompareVal << ") : " << abEdge->toString());
        }

        if (nextCompareVal >= 3)
        {
          continue;
        }
        else if (nextCompareVal == 2)
        {
          if (m_outgoingEdgeIDs.count(abEdge->getIndex()) == 0)
          {
            if (AbruijnGraph2::DEBUG_EXPANSION)
            {
              DEBUG_MSG("Adding to outgoing set");
            }
            m_outgoingEdgeIDs.insert(abEdge->getIndex());
            JumpEdge tempEdge2(*abEdge);
            tempEdge2.popJumpsFront(forPath.getLength());
            m_outgoingEdges.push_back(tempEdge2);
          }
        }
        else
        {
          unsigned int jLength = (abEdge->getLength() == 0) ? 1
              : abEdge->getLength();

          addEdge(curIdx + jLength - 1, *abEdge);

          AbruijnNode& myAbNode = m_nodes.at(curIdx + jLength);
          AbruijnNode* otherNode = abG.getNode(abEdge->getToNodeIndex());

          if (m_pathNodes.count(abEdge->getToNodeIndex()) > 0
              || myAbNode.isCompositeWithOther(*otherNode))
          {
            continue;
          }

          myAbNode.mergeWith(*otherNode);
          m_pathNodes.insert(otherNode->getIndex());

          backBFS.push_back(pair<long, unsigned int> (otherNode->getIndex(),
                                                      curIdx + jLength));
        }
      }
    }

    return true;
  }

  void AbruijnGraph2::ParallelPath::addEdge(const unsigned int& pathIdx,
                                            const JumpEdge& newEdge)
  {

    ParallelEdge& pEdge = m_path.at(pathIdx);
    if (pEdge.addEdge(newEdge))
    {
      m_fWeight += newEdge.getWeight();
      m_rWeight += newEdge.getRWeight();
      if (AbruijnGraph2::DEBUG_EXPANSION)
      {
        DEBUG_MSG("Added edge " << newEdge.toString() << " at index " << pathIdx);
      }
    }
  }

  string AbruijnGraph2::ParallelPath::toString() const
  {
    string st;
    specnets::AAJumps& jumps = specnets::AAJumps::getGlobalJumps();
    for (unsigned int i = 0; i < m_path.size(); i++)
    {
      if (m_path[i].m_jumpCode < 0)
      {
        st += "[]";
      }
      else
      {
        st += jumps.getLabel(m_path[i].m_jumpCode);
      }
    }
    return st;
  }

  const unsigned short AbruijnGraph2::BIN_VERSION = 1;
  const unsigned short AbruijnGraph2::BIN_SUBVERSION = 1;

  const string AbruijnGraph2::BIN_VERSION_ID = "AbruijnGraph_binVersion";
  const string AbruijnGraph2::BIN_SUBVERSION_ID = "AbruijnGraph_binSubVersion";

  const unsigned int AbruijnGraph2::MAX_NUM_JUMPS = 2;
  const unsigned int AbruijnGraph2::MAX_NUM_MODS_PER_JUMP = 0;

  const double AbruijnGraph2::MAX_JUMP_MASS = 320.0;

  const string AbruijnGraph2::STARTPT_SPEC_ID = "startPtSpec";

  const string AbruijnGraph2::ENDPT_SPEC_ID = "endPtSpec";

  unsigned int AbruijnGraph2::NUM_PATHS_PER_NODE = 10;

  bool AbruijnGraph2::DEBUG_EXPANSION = false;

  bool AbruijnGraph2::REVERSE_STARS = true;

  bool AbruijnGraph2::ENFORCE_B_ENDPTS = true;

  double AbruijnGraph2::GAP_PENALTY = 0.10;

  double AbruijnGraph2::PEAK_TOLERANCE = 0.04;

  AbruijnGraph2* AbruijnGraph2::castGraphPtr(BaseGraph* graphPtr)
  {
    AbruijnGraph2* temp = dynamic_cast<AbruijnGraph2*> (graphPtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast BaseGraph \'" << graphPtr->toString() << "\' to an AbruijnGraph2");
      abort();
    }
    return temp;
  }

  AbruijnGraph2* AbruijnGraph2::castNodePtr(Node* nodePtr)
  {
    AbruijnGraph2* temp = dynamic_cast<AbruijnGraph2*> (nodePtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast Node \'" << nodePtr->toString() << "\' to an AbruijnGraph2");
      abort();
    }
    return temp;
  }

  AbruijnGraph2::AbruijnGraph2() :
    BaseGraph(), m_globalJumps(specnets::AAJumps::getGlobalJumps()),
        m_consensusPath(), m_reversed(false), m_nodeLookup(),
        m_nodesToExpand(), m_startPtAssembly()
  {
    initialize();
  }

  AbruijnGraph2::AbruijnGraph2(const AbruijnGraph2& other) :
    BaseGraph(), m_consensusPath(), m_reversed(other.m_reversed),
        m_nodeLookup(), m_globalJumps(specnets::AAJumps::getGlobalJumps()),
        m_nodesToExpand(), m_startPtAssembly()
  {
    this->operator =(other);
  }

  AbruijnGraph2::AbruijnGraph2(const pair<pair<vector<int> , vector<int> > ,
                                   vector<pair<vector<int> , vector<double> > > >& inAbinfo,
                               const SpecSet& assembledSpecs,
                               pair<pair<vector<int> , vector<int> > , vector<
                                   pair<vector<int> , vector<double> > > >* outAbinfo) :
    m_globalJumps(specnets::AAJumps::getGlobalJumps())
  {
    reSequencePPaths(inAbinfo, assembledSpecs, outAbinfo);
  }

  void AbruijnGraph2::initialize()
  {
    BaseGraph::initialize();
    m_consensusPath.clear();
    m_reversed = false;
    m_nodeLookup.clear();
    m_nodesToExpand.clear();
    m_startPtAssembly.clear();
    m_label = "AbruijnGraph2";
  }

  void AbruijnGraph2::compress(vector<long>* outputNewNodeIdxs,
                               vector<long>* outputNewEdgeIdxs)
  {
    vector<long>* newNodeIdxs = (outputNewNodeIdxs == 0) ? 0
        : outputNewNodeIdxs;
    vector<long>* newEdgeIdxs = (outputNewEdgeIdxs == 0) ? new vector<long> ()
        : outputNewEdgeIdxs;

    BaseGraph::compress(newNodeIdxs, newEdgeIdxs);

    if (outputNewEdgeIdxs == 0)
      delete newEdgeIdxs;
  }

  void AbruijnGraph2::reSequencePPaths(const pair<pair<vector<int> ,
                                           vector<int> > , vector<pair<vector<
                                           int> , vector<double> > > >& abinfo,
                                       const SpecSet& assembledSpecs,
                                       pair<pair<vector<int> , vector<int> > ,
                                           vector<pair<vector<int> , vector<
                                               double> > > >* outAbinfo)
  {
    if (ENFORCE_B_ENDPTS)
    {
      DEBUG_MSG("Re-sequencing with b end-points ...");
      m_reSequencePPaths(abinfo, assembledSpecs, REVERSE_STARS, true, outAbinfo);

      DEBUG_VAR(m_consensusPath.toString());
      DEBUG_VAR(m_consensusPath.getWeight());

      return;
    }

    pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> , vector<
        double> > > > yAbinfo;

    DEBUG_MSG("Re-sequencing with y end-points ...");
    m_reSequencePPaths(abinfo, assembledSpecs, REVERSE_STARS, false, &yAbinfo);

    if (m_consensusPath.getLength() == 0)
    {
      return;
    }

    double yPathScore = m_consensusPath.getWeight();

    AbruijnGraph2 yGraph(*this);

    DEBUG_MSG("Re-sequencing with b end-points ...");
    m_reSequencePPaths(abinfo, assembledSpecs, REVERSE_STARS, true, outAbinfo);

    double bPathScore = m_consensusPath.getWeight();

    DEBUG_VAR(yPathScore);
    DEBUG_VAR(bPathScore);

    if (yPathScore > bPathScore + 0.001)
    {
      DEBUG_MSG("Choosing y end-points ...");
      this->copy((BaseGraph&)yGraph);
      if (outAbinfo != 0)
      {
        *outAbinfo = yAbinfo;
      }
    }
    else
    {
      DEBUG_MSG("Choosing b end-points ...");
    }

  }

  Node* AbruijnGraph2::cloneNode(Node& copyNode) const
  {
    Node* copy = new AbruijnNode(*(AbruijnNode::castNodePtr(&copyNode)));
    return copy;
  }

  Edge* AbruijnGraph2::cloneEdge(Edge& copyEdge) const
  {
    Edge* copy = new JumpEdge(*(JumpEdge::castEdgePtr(&copyEdge)));
    return copy;
  }

  Node* AbruijnGraph2::createNode(void) const
  {
    Node* node = new AbruijnNode();
    return node;
  }

  Edge* AbruijnGraph2::createEdge(void) const
  {
    Edge* edge = new JumpEdge();
    return edge;
  }

  AbruijnGraph2 & AbruijnGraph2::operator=(const AbruijnGraph2 &other)
  {
    if (this == &other)
    {
      return *this;
    }
    BaseGraph::operator=((const BaseGraph&)other);

    m_internalCopy(other);

    return *this;
  }

  void AbruijnGraph2::copy(BaseGraph& otherGraph)
  {
    this->operator =(*castGraphPtr(&otherGraph));
  }

  void AbruijnGraph2::addBinaryVersionInfo(map<string, unsigned short>& versions) const
  {
    BaseGraph::addBinaryVersionInfo(versions);
    versions[BIN_VERSION_ID] = BIN_VERSION;
    versions[BIN_SUBVERSION_ID] = BIN_VERSION;
    AbruijnNode nd;
    nd.addBinaryVersionInfo(versions);
    JumpEdge eg;
    eg.addBinaryVersionInfo(versions);
  }

  void AbruijnGraph2::outputConsensusSpectrum(Spectrum& outputSpec) const
  {
    outputSpec.resize(m_consensusPath.m_nodes.size());
    if (outputSpec.size() == 0)
    {
      WARN_MSG("No consensus present, skipping call to outputConsensusSpectrum!!");
      return;
    }
    unsigned int idxUse = 0;
    vector<int> nodeToIdx(m_consensusPath.m_nodes.size(), -1);
    for (unsigned int i = 0; i < m_consensusPath.m_nodes.size(); i++)
    {
      const AbruijnNode* abNode = &m_consensusPath.m_nodes[i];
      if (DEBUG_EXPANSION)
      {
        DEBUG_VAR(abNode->getGraphvizLabel());
      }
      if (abNode->size() == 0)
      {
        continue;
      }
      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Found peak");
      }
      outputSpec[idxUse][1] = 0;
      nodeToIdx[i] = idxUse;
      for (unsigned int j = 0; j < abNode->size(); j++)
      {
        if ((*abNode)[j].getIntensity() > -0.1
            && (*abNode)[j].getSpecID().find(ENDPT_SPEC_ID) == string::npos
            && (*abNode)[j].getSpecID().find(STARTPT_SPEC_ID) == string::npos)
        {
          outputSpec[idxUse][1] += (*abNode)[j].getIntensity();
        }
      }
      idxUse++;
    }
    outputSpec.resize(idxUse);
    outputSpec[0][0] = 0;
    double totMass = 0;

    for (unsigned int i = 1; i < m_consensusPath.m_nodes.size(); i++)
    {
      const ParallelEdge& pEdge = m_consensusPath.m_path[i - 1];
      double jumpMass;
      if (pEdge.m_jumpCode < 0)
      {
        if (pEdge.m_edgeIds.size() == 0)
        {
          ERROR_MSG("Found empty label-free edge!!!");
          abort();
        }
        jumpMass = getEdge(pEdge.m_edgeIds.front())->getMass();
      }
      else
      {
        jumpMass = m_globalJumps[pEdge.m_jumpCode];
      }

      if (nodeToIdx[i] >= 0)
      {
        outputSpec[nodeToIdx[i]][0] = jumpMass + totMass;
      }
      totMass += jumpMass;
    }
    outputSpec.setCharge(1);
    outputSpec.setParentMass(totMass + specnets::AAJumps::massMH);
    outputSpec.sortPeaks();
  }

  void AbruijnGraph2::m_reSequencePPaths(const pair<pair<vector<int> , vector<
                                             int> > , vector<pair<vector<int> ,
                                             vector<double> > > >& inAbinfo,
                                         const SpecSet& assembledSpecs,
                                         const bool reverseStars,
                                         const bool useBEndpts,
                                         pair<
                                             pair<vector<int> , vector<int> > ,
                                             vector<pair<vector<int> , vector<
                                                 double> > > >* outAbinfo)
  {
    // Clear out old graph
    initialize();
    if (outAbinfo != 0)
    {
      outAbinfo->first.first.resize(0);
      outAbinfo->first.second.resize(0);
      outAbinfo->second.resize(0);
    }

    // Index assembled specs w/ their reversed flags
    map<unsigned int, bool> uAssembledSpecs;
    for (unsigned int i = 0; i < inAbinfo.first.first.size(); i++)
    {
      if (inAbinfo.first.first[i] >= 0)
      {
        uAssembledSpecs[inAbinfo.first.first[i]] = inAbinfo.first.second[i]
            == 1;
      }
    }

    DEBUG_VAR(uAssembledSpecs.size());

    // Copy spectra so we can add endpoints, set tolerance, etc.
    SpecSet newAssembledSpecs(uAssembledSpecs.size(), true);

    // Average score in each spectrum
    vector<double> avgSpecScore(uAssembledSpecs.size(), 0);

    vector<vector<bool> > usedPeaks(uAssembledSpecs.size());
    // which peaks were previously assembled
    vector<vector<bool> > alignedPeaks(uAssembledSpecs.size());
    // which peaks we have connected with jump edges
    vector<vector<set<int> > > connectedPeaks(uAssembledSpecs.size());

    set<string> assembldeSpecIDs;

    // map global spectrum index to local
    map<int, unsigned int> specIdxToLoc;
    unsigned int idxUse = 0;
    double totalIntensity = 0, specIntensity = 0;
    double numIntensity = 0;
    for (map<unsigned int, bool>::const_iterator sIt = uAssembledSpecs.begin(); sIt
        != uAssembledSpecs.end(); sIt++)
    {
      newAssembledSpecs[idxUse] = assembledSpecs[sIt->first];
      newAssembledSpecs[idxUse].setPeakTolerance(AbruijnGraph2::PEAK_TOLERANCE);

      if (reverseStars && sIt->second)
      {
        newAssembledSpecs[idxUse].reverse(0);
      }
      newAssembledSpecs[idxUse].scan = sIt->first + 1;
      newAssembledSpecs[idxUse].setReversed(sIt->second);

      if (useBEndpts)
      {
        // Add b end-points and remove y end-points
        newAssembledSpecs[idxUse].addZPMpeaks(-1.0, 0, false);
        newAssembledSpecs[idxUse].maximizeZPMpeaks(-1.0, 0, false);
        int y0Idx =
            newAssembledSpecs[idxUse].findPeaks(specnets::AAJumps::massH2O);
        int
            ykIdx =
                newAssembledSpecs[idxUse].findPeaks(newAssembledSpecs[idxUse].parentMass
                    - specnets::AAJumps::massHion);
        if (y0Idx >= 0)
        {
          newAssembledSpecs[idxUse].removePeak(y0Idx);
          ykIdx--;
        }
        if (ykIdx >= 0)
        {
          newAssembledSpecs[idxUse].removePeak(ykIdx);
        }
      }
      else
      {
        // Add y end-points and remove b end-points
        newAssembledSpecs[idxUse].addZPMpeaks(-1.0, 0, true);
        newAssembledSpecs[idxUse].maximizeZPMpeaks(-1.0, 0, true);
        int b0Idx = newAssembledSpecs[idxUse].findPeaks(0);
        int
            bkIdx =
                newAssembledSpecs[idxUse].findPeaks(newAssembledSpecs[idxUse].parentMass
                    - specnets::AAJumps::massMH);
        if (b0Idx >= 0)
        {
          newAssembledSpecs[idxUse].removePeak(b0Idx);
          bkIdx--;
        }
        if (bkIdx >= 0)
        {
          newAssembledSpecs[idxUse].removePeak(bkIdx);
        }
      }

      newAssembledSpecs[idxUse].consolidatePeaks();

      specIdxToLoc[sIt->first] = idxUse;
      alignedPeaks[idxUse].resize(newAssembledSpecs[idxUse].size(), false);
      connectedPeaks[idxUse].resize(newAssembledSpecs[idxUse].size());
      assembldeSpecIDs.insert(newAssembledSpecs[idxUse].getUniqueID());

      list<int> peaksToRemove;
      Spectrum& spec = newAssembledSpecs[idxUse];

      specIntensity = 0;
      for (unsigned int i = 1; i < spec.size() - 1; i++)
      {
        specIntensity += spec[i][1];
        numIntensity += 1.0;
        avgSpecScore[idxUse] += 1.0;
      }
      totalIntensity += specIntensity;
      avgSpecScore[idxUse] = specIntensity / avgSpecScore[idxUse];

      usedPeaks[idxUse].resize(spec.size(), false);
      idxUse++;
    }
    newAssembledSpecs.index();

    double endptScore = totalIntensity / numIntensity;
    for (unsigned int i = 0; i < newAssembledSpecs.size(); i++)
    {
      Spectrum& spec = newAssembledSpecs[i];
      spec[0][1] = endptScore;
      spec[spec.size() - 1][1] = endptScore;
    }

    bool foundGlue = false;
    set<unsigned int> foundSpecs;
    m_nodeLookup.rehash(newAssembledSpecs.size());

    // Make an AbruijnNode for each consensus vertex
    for (unsigned int i = 0; i < inAbinfo.second.size(); i++)
    {
      AbruijnNode newNode;
      foundGlue = (foundGlue || inAbinfo.second[i].first.size() > 1);
      for (unsigned int j = 0; j < inAbinfo.second[i].first.size(); j++)
      {
        int specIdx = specIdxToLoc[inAbinfo.second[i].first[j]];
        int peakIdx =
            newAssembledSpecs[specIdx].findPeaks(inAbinfo.second[i].second[j]);
        if (peakIdx < 0 || usedPeaks[specIdx][peakIdx])
          continue;

        usedPeaks[specIdx][peakIdx] = true;

        foundSpecs.insert(specIdx);
        newNode.addPeak(newAssembledSpecs[specIdx], peakIdx);
        alignedPeaks[specIdx][peakIdx] = true;
      }
      newNode.setGreen(false);
      AbruijnNode* addedNode = addAbruijnNode(newNode);
    }

    if (!foundGlue)
    {
      WARN_MSG("Found no glues, cannot build abruijn graph!!");
      initialize();
      return;
    }

    for (set<unsigned int>::const_iterator specIdxIt = foundSpecs.begin(); specIdxIt
        != foundSpecs.end(); specIdxIt++)
    {
      unsigned int specIdx = *specIdxIt;
      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Processing spectrum " << newAssembledSpecs[specIdx].scan << " ...");
        newAssembledSpecs[specIdx].output(cerr);
      }
      // Make an AbruijnNode for every other peak
      for (unsigned int i = 0; i < newAssembledSpecs[specIdx].size(); i++)
      {
        if (!alignedPeaks[specIdx][i])
        {
          AbruijnNode tempNode;
          tempNode.loadPeak(newAssembledSpecs[specIdx], i);
          AbruijnNode* newNode = addAbruijnNode(tempNode);
        }
      }
      // Connect nodes with AA jumps in the same spectrum
      for (unsigned int i = 0; i < newAssembledSpecs[specIdx].size(); i++)
      {
        for (unsigned int j = i + 1; j < newAssembledSpecs[specIdx].size(); j++)
        {
          list<JumpEdge*> addedEdges;
          addAbruijnEdges(newAssembledSpecs[specIdx],
                          i,
                          j,
                          avgSpecScore[specIdx],
                          &addedEdges);

          if (addedEdges.size() > 0)
          {
            connectedPeaks[specIdx][i].insert(j);
          }
        }
      }
    }

    unsigned int lastIdx = 0;
    // Connect consecutive aligned peaks with label-free edges (if no edge already exists)
    for (unsigned int i = 1; i < inAbinfo.second.size(); i++)
    {
      const pair<vector<int> , vector<double> >& lastAbVert =
          inAbinfo.second[lastIdx];

      //DEBUG_MSG(lastIdx << " -> " << i);

      map<unsigned int, unsigned int> locAssembledPeaks;
      for (unsigned int j = 0; j < inAbinfo.second[i].first.size(); j++)
      {
        unsigned int specIdx = specIdxToLoc[inAbinfo.second[i].first[j]];
        locAssembledPeaks[specIdx]
            = newAssembledSpecs[specIdx].findClosest(inAbinfo.second[i].second[j]);

        //DEBUG_MSG("Spectrum " << specIdx << ", mass " << inAbinfo.second[i].second[j] << ", peak " << locAssembledPeaks[specIdx]);
      }
      for (unsigned int j = 0; j < lastAbVert.first.size(); j++)
      {
        unsigned int specIdx = specIdxToLoc[lastAbVert.first[j]];
        if (locAssembledPeaks.count(specIdx) == 0)
        {
          continue;
        }
        unsigned int peakIdxFrom =
            newAssembledSpecs[specIdx].findClosest(lastAbVert.second[j]);
        unsigned int peakIdxTo = locAssembledPeaks[specIdx];
        if (peakIdxFrom >= peakIdxTo)
        {
          continue;
          /*
           DEBUG_MSG("Spectrum " << specIdx << ", mass " << lastAbVert.second[j] << ", peak " << peakIdxFrom);
           newAssembledSpecs[specIdx].output(cerr);
           abort();*/
        }
        MZRange peakFrom =
            newAssembledSpecs[specIdx].getPeak(peakIdxFrom).getTightBound();
        MZRange peakTo =
            newAssembledSpecs[specIdx].getPeak(peakIdxTo).getTightBound();
        AbruijnNode
            * nodeFrom =
                m_nodeLookup[newAssembledSpecs[specIdx].getUniqueID()][peakFrom].first;
        AbruijnNode
            * nodeTo =
                m_nodeLookup[newAssembledSpecs[specIdx].getUniqueID()][peakTo].first;
        if (connectedPeaks[specIdx][peakIdxFrom].count(peakIdxTo) == 0)
        {
          if (DEBUG_EXPANSION)
          {
            DEBUG_MSG("Adding label-free edge in " << newAssembledSpecs[specIdx].scan << ": " << peakFrom.getMass() << " -> " << peakTo.getMass());
          }
          addLabelFreeAbEdge(nodeFrom,
                             nodeTo,
                             peakTo.getMass() - peakFrom.getMass(),
                             peakTo.getIntensity(),
                             peakFrom.getIntensity());
        }
      }
      lastIdx = i;
    }

    DEBUG_MSG("Adding endpoint edges ...");
    addEndpointEdges(assembldeSpecIDs, newAssembledSpecs);

    DEBUG_MSG("Merging edges ...");
    mergeAllParallelLabelFreeEdgesByMass();
    mergeAllParallelEdges();

    DEBUG_VAR(this->numNodes());
    DEBUG_VAR(this->numEdges());

    /*
     DEBUG_MSG("Adding missing glues ...");
     for (unsigned int i = 0; i < newAssembledSpecs.size(); i++)
     {
     if (foundSpecs.count(i) == 0)
     {
     continue;
     }

     for (unsigned int j = i + 1; j < newAssembledSpecs.size(); j++)
     {
     if (foundSpecs.count(j) == 0)
     {
     continue;
     }
     addMissingGlues(newAssembledSpecs[i], newAssembledSpecs[j]);
     }
     }
     */

    if (!computeHeaviestParallelPathDAG())
    {
      ERROR_MSG("Graph contains a cycle!!");
      abort();
    }

    if (outAbinfo == 0)
    {
      return;
    }

    // format output abinfo
    outAbinfo->second.resize(0);
    map<int, int> filteredSpecs;

    for (unsigned int i = 0; i < m_consensusPath.m_nodes.size(); i++)
    {
      AbruijnNode* abNode = &m_consensusPath.m_nodes[i];
      if (abNode->size() > 0)
      {
        pair<vector<int> , vector<double> > outVerts;
        for (unsigned int j = 0; j < abNode->size(); j++)
        {
          if ((*abNode)[j].getIntensity() > -0.1
              && (*abNode)[j].getSpecID().find(ENDPT_SPEC_ID) == string::npos
              && (*abNode)[j].getSpecID().find(STARTPT_SPEC_ID) == string::npos)
          {
            Spectrum* spec =
                newAssembledSpecs.getIndex((*abNode)[j].getSpecID());
            filteredSpecs[spec->scan - 1] = (spec->isReversed()) ? 1 : 0;
            outVerts.first.push_back(spec->scan - 1);
            outVerts.second.push_back((*abNode)[j].getMass());
          }
        }
        outAbinfo->second.push_back(outVerts);
      }
    }

    outAbinfo->first.first.resize(filteredSpecs.size());
    outAbinfo->first.second.resize(filteredSpecs.size());
    idxUse = 0;
    for (map<int, int>::const_iterator mIt = filteredSpecs.begin(); mIt
        != filteredSpecs.end(); mIt++)
    {
      outAbinfo->first.first[idxUse] = mIt->first;
      outAbinfo->first.second[idxUse++] = mIt->second;
    }
  }

  void AbruijnGraph2::m_internalCopy(const AbruijnGraph2 &other)
  {
    m_startPtAssembly = other.m_startPtAssembly;
    m_reversed = other.m_reversed;
    m_consensusPath = other.m_consensusPath;

    m_nodeLookup = other.m_nodeLookup;

    for (tr1::unordered_map<string, map<MZRange, pair<AbruijnNode*, TwoValues<
        bool> > > >::iterator sIt = m_nodeLookup.begin(); sIt
        != m_nodeLookup.end(); sIt++)
    {

      for (map<MZRange, pair<AbruijnNode*, TwoValues<bool> > >::iterator mIt =
          sIt->second.begin(); mIt != sIt->second.end(); mIt++)
      {
        mIt->second.first = getNode(mIt->second.first->getIndex());
      }
    }

    m_nodesToExpand.clear();
    for (set<AbruijnNode*>::const_iterator sIt = other.m_nodesToExpand.begin(); sIt
        != other.m_nodesToExpand.end(); sIt++)
    {
      m_nodesToExpand.insert(getNode((*sIt)->getIndex()));
    }

  }

  void AbruijnGraph2::addSpectrum(const Spectrum& prmSpec, bool allEndPts)
  {
    if (m_nodeLookup.count(prmSpec.getUniqueID()) >= 0)
    {
      WARN_MSG("Found duplicate spectrum \'" << prmSpec.getUniqueID() << "\', skipping call to addSpectrum!");
      return;
    }

    //DEBUG_VAR(prmSpec.getUniqueID());

    double avgScore = 0.0;
    double numScore = 0.0;

    for (unsigned int i = 0; i < prmSpec.size(); i++)
    {
      // Create a node for every peak in the spectrum
      AbruijnNode tempNode;
      tempNode.loadPeak(prmSpec, i);
      if (allEndPts)
      {
        tempNode[0].setEndPt(true);
      }
      AbruijnNode* newNode = addAbruijnNode(tempNode);

      if (i > 0 && i < prmSpec.size() - 1)
      {
        avgScore += prmSpec[i][1];
        numScore += 1.0;
      }
    }

    avgScore = avgScore / numScore;

    for (unsigned int i = 0; i < prmSpec.size(); i++)
    {
      for (unsigned int j = i + 1; j < prmSpec.size(); j++)
      {
        addAbruijnEdges(prmSpec, i, j, avgScore);
      }
    }
  }

  void AbruijnGraph2::addGlues(const SpectrumAlignment& prmAlign)
  {
    string spec1ID(prmAlign.getSpec1ID());
    string spec2ID(prmAlign.getSpec2ID());

    if (m_nodeLookup.count(spec1ID) == 0)
    {
      ERROR_MSG("Cannot find spectrum \'" << spec1ID << "\'");
      abort();
    }
    if (m_nodeLookup.count(spec2ID) == 0)
    {
      ERROR_MSG("Cannot find spectrum \'" << spec2ID << "\'");
      abort();
    }

    list<pair<MZRange, MZRange> > matchedPeaksB;
    list<pair<MZRange, MZRange> > matchedPeaksY;
    prmAlign.outputMatchedPeaks(matchedPeaksB, matchedPeaksY);

    if (matchedPeaksB.size() == 0)
    {
      WARN_MSG("Skipping alignment between \'" << spec1ID << "\' and \'" << spec2ID << "\' with 0 matched peaks");
      return;
    }

    list<pair<MZRange, MZRange> >::const_iterator mIt;

    for (mIt = matchedPeaksB.begin(); mIt != matchedPeaksB.end(); mIt++)
    {
      if (m_nodeLookup[spec1ID].count(mIt->first) == 0)
      {
        ERROR_MSG("Cannot find peak " << mIt->first.getMass() << " +/- " << mIt->first.getTolerance() << " in spectrum \'" << spec1ID << "\'");
        abort();
      }
      if (m_nodeLookup[spec2ID].count(mIt->second) == 0)
      {
        ERROR_MSG("Cannot find peak " << mIt->second.getMass() << " +/- " << mIt->second.getTolerance() << " in spectrum \'" << spec2ID << "\'");
        abort();
      }
    }

    AbruijnNode* nodePeak1;
    AbruijnNode* nodePeak2;

    for (mIt = matchedPeaksB.begin(); mIt != matchedPeaksB.end(); mIt++)
    {
      nodePeak1 = m_nodeLookup[spec1ID][mIt->first.getTightBound()].first;
      nodePeak2 = m_nodeLookup[spec2ID][mIt->second.getTightBound()].first;

      AssembledPeak* aPeak1 = nodePeak1->findSpectrum(spec1ID);
      AssembledPeak* aPeak2 = nodePeak2->findSpectrum(spec2ID);

      if (aPeak1 == 0 || aPeak2 == 0)
      {
        ERROR_MSG("Could not locate spectra \'" << spec1ID << "\' and \'" << spec2ID << "\'");
        abort();
      }

      aPeak1->setSymmetry(AssembledPeak::Symmetry_B);
      aPeak2->setSymmetry(AssembledPeak::Symmetry_B);

      if (nodePeak1 == nodePeak2)
      {
        continue;
      }
      glueNodes(nodePeak1, &nodePeak2);
    }

    for (mIt = matchedPeaksY.begin(); mIt != matchedPeaksY.end(); mIt++)
    {
      nodePeak1 = m_nodeLookup[spec1ID][mIt->first.getTightBound()].first;
      nodePeak2 = m_nodeLookup[spec2ID][mIt->second.getTightBound()].first;

      AssembledPeak* aPeak1 = nodePeak1->findSpectrum(spec1ID);
      AssembledPeak* aPeak2 = nodePeak2->findSpectrum(spec2ID);

      if (aPeak1 == 0 || aPeak2 == 0)
      {
        ERROR_MSG("Could not locate spectra \'" << spec1ID << "\' and \'" << spec2ID << "\'");
        abort();
      }

      if (aPeak1->getSymmetry() == AssembledPeak::Symmetry_unknown)
        aPeak1->setSymmetry(AssembledPeak::Symmetry_Y);
      if (aPeak2->getSymmetry() == AssembledPeak::Symmetry_unknown)
        aPeak2->setSymmetry(AssembledPeak::Symmetry_Y);
    }
  }

  AbruijnNode* AbruijnGraph2::addAbruijnNode(AbruijnNode& copyNode)
  {
    AbruijnNode* newNode =
        AbruijnNode::castNodePtr(BaseGraph::addNode((Node&)copyNode));

    if (DEBUG_EXPANSION)
    {
      DEBUG_MSG("Created node " << newNode << " - " << newNode->getGraphvizLabel());
    }

    for (unsigned int i = 0; i < newNode->size(); i++)
    {
      AssembledPeak* aPeak = &newNode->operator [](i);
      const string& specID = aPeak->getSpecID();

      if (m_nodeLookup.count(specID) == 0)
      {
        map<MZRange, pair<AbruijnNode*, TwoValues<bool> > > subMap;
        subMap[aPeak->getTightBound()]
            = pair<AbruijnNode*, TwoValues<bool> > (newNode,
                                                    TwoValues<bool> (false,
                                                                     false));
        m_nodeLookup[specID] = subMap;
      }
      else
      {
        m_nodeLookup[specID][aPeak->getTightBound()].first = newNode;
      }
    }
    return newNode;
  }

  JumpEdge* AbruijnGraph2::addLabelFreeAbEdge(AbruijnNode* from,
                                              AbruijnNode* to,
                                              const double& mass,
                                              const double& weight,
                                              const double& rWeight)
  {
    MZRange massR(mass, 0, AbruijnGraph2::PEAK_TOLERANCE);
    JumpEdge newEdge("", &massR);

    JumpEdge* addedEdge =
        JumpEdge::castEdgePtr(this->addEdge(from, to, newEdge));
    return addedEdge;
  }

  void AbruijnGraph2::recoverSourceNodeScores()
  {
    for (unsigned long i = 0; i <= this->maxNodeIdx(); i++)
    {
      Node* node = this->lookupNode(i);
      if (node == 0)
        continue;

      AbruijnNode* abNode = AbruijnNode::castNodePtr(node);
      double totalMissedScore = 0;
      double totalMissedScoreR = 0;
      for (unsigned int j = 0; j < abNode->size(); j++)
      {
        if ((*abNode)[j].getSpecID().find(AbruijnGraph2::STARTPT_SPEC_ID)
            != string::npos
            || (*abNode)[j].getSpecID().find(AbruijnGraph2::ENDPT_SPEC_ID)
                != string::npos || (*abNode)[j].getIntensity() < -0.1)
        {
          continue;
        }

        if (!m_nodeLookup[(*abNode)[j].getSpecID()][(*abNode)[j].getTightBound()].second[0])
        {
          totalMissedScore += (*abNode)[j].getIntensity();
          if (DEBUG_EXPANSION)
            DEBUG_MSG("Adding missing score " << (*abNode)[j].getIntensity() << " from peak " << (*abNode)[j].getSpecID() << " / " << (*abNode)[j].getMass() << " to node " << abNode->getIndex());
        }
        if (!m_nodeLookup[(*abNode)[j].getSpecID()][(*abNode)[j].getTightBound()].second[1])
        {
          totalMissedScoreR += (*abNode)[j].getIntensity();
          if (DEBUG_EXPANSION)
            DEBUG_MSG("Adding missing reverse score " << (*abNode)[j].getIntensity() << " from peak " << (*abNode)[j].getSpecID() << " / " << (*abNode)[j].getMass() << " to node " << abNode->getIndex());
        }
      }

      if (totalMissedScore > 0)
      {
        const IndexVector& inEdges = this->getInEdges(node);
        for (unsigned long j = 0; j < inEdges.size(); j++)
        {
          if (inEdges[j] < 0)
            continue;
          JumpEdge* abEdge = getEdge(inEdges[j]);
          if (DEBUG_EXPANSION)
            DEBUG_MSG("Adding " << totalMissedScore << " to " << abEdge->toString());

          abEdge->setWeight(abEdge->getWeight() + totalMissedScore);
        }
      }
      if (totalMissedScoreR > 0)
      {
        const IndexVector& outEdges = this->getOutEdges(node);
        for (unsigned long j = 0; j < outEdges.size(); j++)
        {
          if (outEdges[j] < 0)
            continue;
          JumpEdge* abEdge = getEdge(outEdges[j]);
          if (DEBUG_EXPANSION)
            DEBUG_MSG("Adding " << totalMissedScoreR << " to " << abEdge->toString());

          abEdge->setRWeight(abEdge->getRWeight() + totalMissedScoreR);
        }
      }
    }
  }

  void AbruijnGraph2::reverseEdgeWeights()
  {
    for (unsigned long i = 0; i <= this->maxEdgeIdx(); i++)
    {
      Edge* edge = this->lookupEdge(i);
      if (edge == 0)
        continue;

      JumpEdge* abEdge = JumpEdge::castEdgePtr(edge);
      const double weight = abEdge->getWeight();
      abEdge->setWeight(abEdge->getRWeight());
      abEdge->setRWeight(weight);
    }
  }

  void AbruijnGraph2::addEndpointEdges(const set<string>& assembledSpecIDs,
                                       const SpecSet& assmebledSpecs)
  {
    BaseGraph startPtGraph;
    map<string, Node*> startNodeLookup;
    map<Node*, string> startSpecLookup;
    map<Edge*, double> startJumpMassLookup;

    BaseGraph endPtGraph;
    map<string, Node*> endNodeLookup;
    map<Node*, string> endSpecLookup;
    map<Edge*, double> endJumpMassLookup;

    Node* startBegin = 0;
    Node* endBegin = 0;

    for (set<string>::const_iterator sIt = assembledSpecIDs.begin(); sIt
        != assembledSpecIDs.end(); sIt++)
    {
      const Spectrum* spec = assmebledSpecs.getIndex(*sIt);

      if (m_nodeLookup.count(*sIt) == 0)
      {
        continue;
      }

      Node* newNode = startPtGraph.addNode();
      startNodeLookup[*sIt] = newNode;
      startSpecLookup[newNode] = *sIt;

      newNode = endPtGraph.addNode();
      endNodeLookup[*sIt] = newNode;
      endSpecLookup[newNode] = *sIt;
    }

    for (set<string>::const_iterator sIt = assembledSpecIDs.begin(); sIt
        != assembledSpecIDs.end(); sIt++)
    {
      const Spectrum* spec = assmebledSpecs.getIndex(*sIt);
      if (m_nodeLookup.count(*sIt) == 0)
      {
        continue;
      }

      Node* fromNode = startNodeLookup[*sIt];

      for (unsigned int i = 0; i < spec->size(); i++)
      {
        MZRange peak = spec->getPeak(i).getTightBound();
        if (m_nodeLookup[*sIt].count(peak) == 0)
        {
          continue;
        }
        AbruijnNode* abNode = m_nodeLookup[*sIt][peak].first;
        for (unsigned int j = 0; j < abNode->size(); j++)
        {
          if ((*abNode)[j].getIntensity() < -0.1)
          {
            continue;
          }
          const string& nextID = (*abNode)[j].getSpecID();
          if (nextID == *sIt)
          {
            continue;
          }

          double shift = peak.getMass() - (*abNode)[j].getMass();
          Node* toNode = startNodeLookup[nextID];

          if (toNode->getIndex() == fromNode->getIndex())
          {
            continue;
          }

          unsigned int
              toSpecIdx =
                  assmebledSpecs.getIndex(nextID)->findClosest((*abNode)[j].getMass());

          Edge* newEdge = startPtGraph.addEdge(fromNode, toNode);
          newEdge->setWeight((double)(toSpecIdx + 1));
          startJumpMassLookup[newEdge] = shift;
          newEdge = startPtGraph.addEdge(toNode, fromNode);
          newEdge->setWeight((double)(i + 1));
          startJumpMassLookup[newEdge] = 0.0 - shift;

          if (shift > 0)
          {
            if (startBegin == 0 || startBegin == toNode)
            {
              startBegin = fromNode;
            }
          }
          else
          {
            if (startBegin == 0 || startBegin == fromNode)
            {
              startBegin = toNode;
            }
          }
        }
      }
      for (int i = spec->size() - 1; i >= 0; i--)
      {
        MZRange peak = spec->getPeak(i).getTightBound();
        if (m_nodeLookup[*sIt].count(peak) == 0)
        {
          continue;
        }
        AbruijnNode* abNode = m_nodeLookup[*sIt][peak].first;
        float baseDiff = spec->parentMass - peak.getMass();
        Node* fromNode = endNodeLookup[*sIt];

        for (unsigned int j = 0; j < abNode->size(); j++)
        {
          if ((*abNode)[j].getIntensity() < -0.1)
          {
            continue;
          }
          const string& nextID = (*abNode)[j].getSpecID();
          const Spectrum* spec2 = assmebledSpecs.getIndex(nextID);
          if (nextID == *sIt)
          {
            continue;
          }

          double destDiff = spec2->parentMass - (*abNode)[j].getMass();
          double shift = destDiff - baseDiff;
          Node* toNode = endNodeLookup[nextID];

          if (toNode->getIndex() == fromNode->getIndex())
          {
            continue;
          }

          unsigned int toSpecIdx = spec2->findClosest((*abNode)[j].getMass());

          Edge* newEdge = endPtGraph.addEdge(fromNode, toNode);
          newEdge->setWeight((double)(spec2->size() - toSpecIdx));
          endJumpMassLookup[newEdge] = shift;
          newEdge = endPtGraph.addEdge(toNode, fromNode);
          newEdge->setWeight((double)(spec->size() - i));
          endJumpMassLookup[newEdge] = 0.0 - shift;
          if (shift > 0)
          {
            if (endBegin == 0 || endBegin == toNode)
            {
              endBegin = fromNode;
            }
          }
          else
          {
            if (endBegin == 0 || endBegin == fromNode)
            {
              endBegin = toNode;
            }
            //DEBUG_MSG("Adding edge " << newEdge->toString() << " w/ shift " << 0.0 - shift);
          }
        }
      }
    }

    if (startBegin == 0)
    {
      ERROR_MSG("Failed to locate begin node in start-pt graph!!");
      abort();
    }
    if (endBegin == 0)
    {
      ERROR_MSG("Failed to locate begin node in end-pt graph!!");
      abort();
    }

    Tree startPtTree;
    Tree endPtTree;

    startPtGraph.getLightestPaths(startBegin, startPtTree);
    endPtGraph.getLightestPaths(endBegin, endPtTree);

    Spectrum startPtSpec;
    startPtSpec.fileName = AbruijnGraph2::STARTPT_SPEC_ID;
    startPtSpec.scan = 1;
    map<float, set<AbruijnNode*> > nodesToMergeStart;

    m_startPtAssembly.clear();

    for (Tree::iterator tIt = startPtTree.begin(); tIt != startPtTree.end(); tIt++)
    {
      const string& nextSpecID = startSpecLookup[tIt->first];

      if (m_startPtAssembly.count(nextSpecID) > 0)
      {
        continue;
      }
      const Spectrum* spec2 = assmebledSpecs.getIndex(nextSpecID);

      double jumpMass = 0;
      Tree::iterator curNode = tIt;
      while (curNode->first != startBegin)
      {
        jumpMass += startJumpMassLookup[curNode->second.second];
        curNode
            = startPtTree.find(startPtGraph.getNode(curNode->second.second->getFromNodeIndex()));
      }

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Mass " << jumpMass << ": " << nextSpecID);
      }
      m_startPtAssembly[nextSpecID] = jumpMass;

      MZRange massB02 = spec2->getPeak(0);
      MZRange massY02 =
          spec2->getPeak(spec2->findClosest(specnets::AAJumps::massH2O));

      AbruijnNode* spec2B0Node =
          m_nodeLookup[nextSpecID][massB02.getTightBound()].first;

      massB02.setMass(jumpMass);

      if (massB02 == 0)
      {
        massB02.setMass(0);
      }
      else if (jumpMass < 0)
      {
        for (unsigned int i = 0; i < startPtSpec.size(); i++)
        {
          startPtSpec[i][0] += 0.0 - jumpMass;
        }
        massB02.setMass(0);
      }

      int peakIdx = startPtSpec.findPeaks(jumpMass, massB02.getTolerance());
      if (peakIdx >= 0)
      {
        startPtSpec[peakIdx][1] += massB02.getIntensity();
        nodesToMergeStart[startPtSpec[peakIdx][0]].insert(spec2B0Node);
      }
      else
      {
        startPtSpec.insertPeak(&massB02);
        set<AbruijnNode*> newSet;
        newSet.insert(spec2B0Node);
        nodesToMergeStart[massB02.getMass()] = newSet;
      }
    }

    startPtSpec.parentCharge = 2;
    startPtSpec.parentMass = startPtSpec[startPtSpec.size() - 1][0]
        + specnets::AAJumps::massMH;

    Spectrum endPtSpec;
    endPtSpec.fileName = AbruijnGraph2::ENDPT_SPEC_ID;
    endPtSpec.scan = 1;
    map<float, set<AbruijnNode*> > nodesToMergeEnd;

    for (Tree::iterator tIt = endPtTree.begin(); tIt != endPtTree.end(); tIt++)
    {
      const string& nextSpecID = endSpecLookup[tIt->first];
      const Spectrum* spec2 = assmebledSpecs.getIndex(nextSpecID);

      double jumpMass = 0;
      Tree::iterator curNode = tIt;
      while (curNode->first != endBegin)
      {
        jumpMass += endJumpMassLookup[curNode->second.second];
        curNode
            = endPtTree.find(endPtGraph.getNode(curNode->second.second->getFromNodeIndex()));
      }

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Mass " << jumpMass << ": " << nextSpecID);
      }

      MZRange massBk2 = spec2->getPeak(spec2->findClosest(spec2->parentMass
          - specnets::AAJumps::massMH));
      MZRange massYk2 = spec2->getPeak(spec2->findClosest(spec2->parentMass
          - specnets::AAJumps::massHion));

      AbruijnNode* spec2BkNode =
          m_nodeLookup[nextSpecID][massBk2.getTightBound()].first;

      massBk2.setMass(jumpMass);

      if (massBk2 == 0)
      {
        massBk2.setMass(0);
      }
      else if (jumpMass < 0)
      {
        for (unsigned int i = 0; i < endPtSpec.size(); i++)
        {
          endPtSpec[i][0] += 0.0 - jumpMass;
        }
        massBk2.setMass(0);
      }

      int peakIdx = endPtSpec.findPeaks(jumpMass, massBk2.getTolerance());
      if (peakIdx >= 0)
      {
        endPtSpec[peakIdx][1] += massBk2.getIntensity();
        nodesToMergeEnd[endPtSpec[peakIdx][0]].insert(spec2BkNode);
      }
      else
      {
        endPtSpec.insertPeak(&massBk2);
        set<AbruijnNode*> newSet;
        newSet.insert(spec2BkNode);
        nodesToMergeEnd[massBk2.getMass()] = newSet;
      }
    }

    endPtSpec.parentCharge = 2;
    endPtSpec.parentMass = endPtSpec[endPtSpec.size() - 1][0]
        + specnets::AAJumps::massMH;

    set<AbruijnNode*> removedNodes;
    list<int> removedPeaks;

    // Make an AbruijnNode for every peak
    for (unsigned int i = 0; i < startPtSpec.size(); i++)
    {
      list<AbruijnNode*>
          mergeNodes(nodesToMergeStart[startPtSpec[i][0]].begin(),
                     nodesToMergeStart[startPtSpec[i][0]].end());

      //DEBUG_VAR(mergeNodes.size());

      bool foundReal = false;
      for (list<AbruijnNode*>::iterator sIt = mergeNodes.begin(); sIt
          != mergeNodes.end(); sIt++)
      {
        if (removedNodes.count(*sIt) == 0)
        {
          foundReal = true;
          removedNodes.insert(*sIt);
        }
        else
        {
          *sIt = 0;
        }
      }
      if (!foundReal)
      {
        removedPeaks.push_back(i);
        continue;
      }

      AbruijnNode tempNode;
      tempNode.loadPeak(startPtSpec, i);
      AbruijnNode* newNode = addAbruijnNode(tempNode);

      for (list<AbruijnNode*>::iterator sIt = mergeNodes.begin(); sIt
          != mergeNodes.end(); sIt++)
      {
        if (*sIt != 0 && !newNode->isCompositeWithOther(*(*sIt)))
        {
          glueNodes(newNode, &(*sIt));
        }
      }

      if (DEBUG_EXPANSION)
      {
        DEBUG_VAR(newNode->getGraphvizLabel());
      }
    }
    startPtSpec.removePeaks(removedPeaks);

    removedPeaks.clear();
    // Connect nodes with AA jumps in the same spectrum
    /*
     for (unsigned int i = 0; i < startPtSpec.size(); i++)
     {
     for (unsigned int j = i + 1; j < startPtSpec.size(); j++)
     {
     addAbruijnEdges(startPtSpec, i, j);
     }
     }*/

    // Make an AbruijnNode for every peak
    for (unsigned int i = 0; i < endPtSpec.size(); i++)
    {
      list<AbruijnNode*> mergeNodes(nodesToMergeEnd[endPtSpec[i][0]].begin(),
                                    nodesToMergeEnd[endPtSpec[i][0]].end());

      //DEBUG_VAR(mergeNodes.size());

      bool foundReal = false;
      for (list<AbruijnNode*>::iterator sIt = mergeNodes.begin(); sIt
          != mergeNodes.end(); sIt++)
      {
        if (removedNodes.count(*sIt) == 0)
        {
          foundReal = true;
          removedNodes.insert(*sIt);
        }
        else
        {
          *sIt = 0;
        }
      }
      if (!foundReal)
      {
        removedPeaks.push_back(i);
        continue;
      }

      AbruijnNode tempNode;
      tempNode.loadPeak(endPtSpec, i);
      AbruijnNode* newNode = addAbruijnNode(tempNode);

      for (list<AbruijnNode*>::iterator sIt = mergeNodes.begin(); sIt
          != mergeNodes.end(); sIt++)
      {
        if (*sIt != 0 && !newNode->isCompositeWithOther(*(*sIt)))
        {
          glueNodes(newNode, &(*sIt));
        }
      }

      if (DEBUG_EXPANSION)
      {
        DEBUG_VAR(newNode->getGraphvizLabel());
      }
    }
    endPtSpec.removePeaks(removedPeaks);

    if (DEBUG_EXPANSION)
    {
      DEBUG_MSG("Start-pt spectrum:");
      startPtSpec.output(cerr);
      cerr << "\n";
      DEBUG_MSG("End-pt spectrum:");
      endPtSpec.output(cerr);
    }

    // Connect nodes with AA jumps in the same spectrum
    /*
     for (unsigned int i = 0; i < endPtSpec.size(); i++)
     {
     for (unsigned int j = i + 1; j < endPtSpec.size(); j++)
     {
     addAbruijnEdges(endPtSpec, i, j);
     }
     }*/
  }

  void AbruijnGraph2::addAbruijnEdges(const Spectrum& prmSpec,
                                      const unsigned int& peakIdxFrom,
                                      const unsigned int& peakIdxTo,
                                      const double& avgSpecScore,
                                      list<JumpEdge*>* addedEdges)
  {
    if (addedEdges != 0)
    {
      addedEdges->clear();
    }

    if (peakIdxTo <= peakIdxFrom)
    {
      ERROR_MSG("Invalid indices for jump \'" << peakIdxFrom << " -> " << peakIdxTo << "\'");
      abort();
    }
    const string specID = prmSpec.getUniqueID();
    MZRange peakFrom = prmSpec.getPeak(peakIdxFrom);
    MZRange peakTo = prmSpec.getPeak(peakIdxTo);

    if (m_nodeLookup.count(specID) == 0
        || m_nodeLookup[specID].count(peakFrom.getTightBound()) == 0
        || m_nodeLookup[specID].count(peakTo.getTightBound()) == 0)
    {
      ERROR_MSG("Cannot find AbruijnNodes for peaks " << peakIdxFrom << " and " << peakIdxTo << " in spectrum " << specID);
      abort();
    }

    pair<AbruijnNode*, TwoValues<bool> >& pairTo =
        m_nodeLookup[specID][peakTo.getTightBound()];
    pair<AbruijnNode*, TwoValues<bool> >& pairFrom =
        m_nodeLookup[specID][peakFrom.getTightBound()];
    AbruijnNode* nodeFrom =
        m_nodeLookup[specID][peakFrom.getTightBound()].first;
    AbruijnNode* nodeTo = pairTo.first;

    list<pair<string, double> > foundJumps;
    vector<string> jumpPrefixes;
    vector<string> jumpSuffixes;
    bool foundRedundancy;

    double foundJumpMass = prmSpec[peakIdxTo][0] - prmSpec[peakIdxFrom][0];
    double jumpTol = prmSpec.getTolerance(peakIdxTo)
        + prmSpec.getTolerance(peakIdxFrom);

    if (foundJumpMass > MAX_JUMP_MASS)
      return;

    // match peak mass difference to the mass of an amino acid combination
    m_globalJumps.findJumpsWLabels(foundJumpMass,
                                   jumpTol,
                                   foundJumps,
                                   -1,
                                   AbruijnGraph2::MAX_NUM_JUMPS,
                                   AbruijnGraph2::MAX_NUM_MODS_PER_JUMP);

    /*
     if (DEBUG_EXPANSION)
     {
     DEBUG_MSG("foundJumpMass = " << foundJumpMass << ", jumpTol = " << jumpTol << ", foundJumps.size() = " << foundJumps.size());
     }
     */
    set<unsigned int> sizeSet;
    unsigned int maxLength = 0;
    for (list<pair<string, double> >::iterator jIt = foundJumps.begin(); jIt
        != foundJumps.end(); jIt++)
    {
      sizeSet.insert(jIt->first.length());
      maxLength = max(maxLength, (unsigned int)jIt->first.length());
    }

    // iterate over all possible AA combinations
    for (list<pair<string, double> >::iterator jIt = foundJumps.begin(); jIt
        != foundJumps.end(); jIt++)
    {
      const unsigned int numJumps = specnets::AAJumps::getNumJumps(jIt->first);
      double jumpMass = m_globalJumps.getModPeptideMass(jIt->first);
      bool modified = jIt->first.find("(") != string::npos;
      string noModJump = specnets::AAJumps::stripMods(jIt->first);
      double edgeWeight = prmSpec[peakIdxTo][1] - abs(jumpMass - foundJumpMass);
      double rEdgeWeight = prmSpec[peakIdxFrom][1] - abs(jumpMass
          - foundJumpMass);

      if (numJumps > 1)
      {
        // don't add redundant edges (no "PETI" edge if already exists a "PE" or "P" edge)
        specnets::AAJumps::getPrefixJumps(jIt->first, jumpPrefixes);
        foundRedundancy = false;
        // look for redundant jumps originating from source vertex
        for (unsigned int pIt = 0; pIt < jumpPrefixes.size() - 1; pIt++)
        {
          double intermedPeak = prmSpec[peakIdxFrom][0]
              + m_globalJumps.getModPeptideMass(jumpPrefixes[pIt]);
          if (prmSpec.findPeaks(intermedPeak, prmSpec.getTolerance(peakIdxFrom))
              >= 0)
          {
            foundRedundancy = true;
            break;
          }
        }

        if (foundRedundancy)
        {
          continue;
        }

        edgeWeight -= avgSpecScore * GAP_PENALTY * (double)(numJumps - 1);
        rEdgeWeight -= avgSpecScore * GAP_PENALTY * (double)(numJumps - 1);
      }

      //DEBUG_MSG("Adding edge " << i << " -> " << j << " (" << jIt->first << ")");

      // add intermediate nodes for jumps of length > 1
      // addExpandedPath(from, to, prmSpec[j][1], jIt->first);

      if (modified)
      {
        ERROR_MSG("Adding modified edge " << jIt->first);
        abort();
      }

      list<string> tempList;
      tempList.push_back(prmSpec.getUniqueID());
      JumpEdge newEdge(jIt->first);
      newEdge.setWeight(edgeWeight);
      newEdge.setRWeight(rEdgeWeight);
      JumpEdge* addedEdge = JumpEdge::castEdgePtr(this->addEdge(nodeFrom,
                                                                nodeTo,
                                                                newEdge));
      pairTo.second[0] = true;
      pairFrom.second[1] = true;

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding edge " << addedEdge->toString());
      }

      if (addedEdges != 0)
      {
        addedEdges->push_back(addedEdge);
      }
    }
  }

  void AbruijnGraph2::mergeParallelEdges(AbruijnNode* from,
                                         AbruijnNode* to,
                                         bool mergeMods)
  {
    list<list<JumpEdge*> > parallelSets;

    const IndexVector& outEdges = this->getOutEdges(from->getIndex());
    for (unsigned long i = 0; i < outEdges.size(); i++)
    {
      if (outEdges[i] < 0)
        continue;

      JumpEdge* edge = getEdge(outEdges[i]);
      if (edge->getToNodeIndex() != to->getIndex())
      {
        continue;
      }

      bool foundSame = false;
      for (list<list<JumpEdge*> >::iterator lIt = parallelSets.begin(); lIt
          != parallelSets.end(); lIt++)
      {
        if (edge->compareTo(*(lIt->front())) == 0)
        {
          foundSame = true;
          lIt->push_back(edge);
          break;
        }
      }

      if (!foundSame)
      {
        list<JumpEdge*> newL;
        newL.push_back(edge);
        parallelSets.push_back(newL);
      }
    }

    list<JumpEdge*> edgesToRemove;
    list<JumpEdge> edgesToAdd;
    for (list<list<JumpEdge*> >::const_iterator lIt = parallelSets.begin(); lIt
        != parallelSets.end(); lIt++)
    {
      if (lIt->size() <= 1)
      {
        continue;
      }

      JumpEdge templateEdge(*(lIt->front()));
      templateEdge.setWeight(0);
      templateEdge.setRWeight(0);
      for (list<JumpEdge*>::const_iterator eIt = lIt->begin(); eIt
          != lIt->end(); eIt++)
      {
        templateEdge.addWeight((*eIt)->getWeight());
        templateEdge.addRWeight((*eIt)->getRWeight());
      }

      edgesToRemove.insert(edgesToRemove.end(), lIt->begin(), lIt->end());
      edgesToAdd.push_back(templateEdge);
    }

    for (list<JumpEdge*>::iterator eIt = edgesToRemove.begin(); eIt
        != edgesToRemove.end(); eIt++)
    {
      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Removing edge " << (*eIt)->getFromNodeIndex() << " -> " << (*eIt)->getToNodeIndex() << " index=" << (*eIt)->getIndex() << ", jump=" << (*eIt)->getMass() << ", weight=" << (*eIt)->getWeight());
      }
      this->removeEdge((Edge**)(&(*eIt)));
    }

    for (list<JumpEdge>::iterator eIt = edgesToAdd.begin(); eIt
        != edgesToAdd.end(); eIt++)
    {
      JumpEdge* addedEdge =
          JumpEdge::castEdgePtr(this->addEdge(from, to, *eIt));

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding edge " << from->getIndex() << " -> " << to->getIndex() << " index=" << addedEdge->getIndex() << ", jump=" << addedEdge->getMass() << ", weight=" << addedEdge->getWeight());
      }
    }
  }

  void AbruijnGraph2::mergeAllParallelEdges(bool mergeMods)
  {
    set<long> seenNeij;
    for (unsigned long i = 0; i <= this->maxNodeIdx(); i++)
    {
      Node* from = this->lookupNode(i);
      if (from == 0)
        continue;

      const IndexVector& outEdges = this->getOutEdges(from->getIndex());
      seenNeij.clear();
      for (unsigned long j = 0; j < outEdges.size(); j++)
      {
        if (outEdges[j] < 0 || seenNeij.count(outEdges[j]) > 0)
          continue;
        seenNeij.insert(outEdges[j]);

        Node* to = this->getNode(this->getEdge(outEdges[j])->getToNodeIndex());

        mergeParallelEdges(AbruijnNode::castNodePtr(from),
                           AbruijnNode::castNodePtr(to),
                           mergeMods);
      }
    }
  }

  void AbruijnGraph2::mergeParallelLabelFreeEdgesByMass(AbruijnNode* from,
                                                        AbruijnNode* to)
  {
    map<MZRange, list<JumpEdge*> > parallelSets;
    set<MZRange> jumpsToMerge;

    const IndexVector& outEdges = this->getOutEdges(from->getIndex());
    for (unsigned long i = 0; i < outEdges.size(); i++)
    {
      if (outEdges[i] < 0)
        continue;

      Edge* edge = this->getEdge(outEdges[i]);
      if (edge->getToNodeIndex() != to->getIndex())
      {
        continue;
      }
      JumpEdge* abEdge = JumpEdge::castEdgePtr(edge);
      float jumpMass = abEdge->getMass();
      double jumpWeight = abEdge->getWeight();

      if (abEdge->getLength() > 0)
      {
        continue;
      }
      //DEBUG_MSG(from->getIndex() << " -> " << to->getIndex() << " : jumpMass = " << jumpMass);
      MZRange jump(jumpMass, jumpWeight, AbruijnGraph2::PEAK_TOLERANCE);

      if (parallelSets.count(jump) > 0)
      {
        parallelSets[jump].push_back(abEdge);
        jumpsToMerge.insert(jump);
      }
      else
      {
        list<JumpEdge*> tempList;
        tempList.push_back(abEdge);
        parallelSets[jump] = tempList;
      }
      //DEBUG_TRACE;
    }

    if (jumpsToMerge.size() == 0)
    {
      return;
    }
    list<JumpEdge*> edgesToRemove;
    for (set<MZRange>::const_iterator lIt = jumpsToMerge.begin(); lIt
        != jumpsToMerge.end(); lIt++)
    {
      list<JumpEdge*>& edgeList = parallelSets[*lIt];
      JumpEdge templateEdge(*edgeList.front());
      double totalWeight = 0;
      double totalRWeight = 0;

      for (list<JumpEdge*>::iterator eIt = edgeList.begin(); eIt
          != edgeList.end(); eIt++)
      {
        edgesToRemove.push_back(*eIt);
        totalWeight += (*eIt)->getWeight();
        totalRWeight += (*eIt)->getRWeight();
      }
      templateEdge.setWeight(totalWeight);
      templateEdge.setRWeight(totalRWeight);

      Edge* addedEdge = this->addEdge(from, to, templateEdge);
    }
    for (list<JumpEdge*>::iterator eIt = edgesToRemove.begin(); eIt
        != edgesToRemove.end(); eIt++)
    {
      this->removeEdge((Edge**)(&(*eIt)));
    }
  }

  void AbruijnGraph2::mergeAllParallelLabelFreeEdgesByMass()
  {
    set<long> seenNeij;
    for (unsigned long i = 0; i <= this->maxNodeIdx(); i++)
    {
      Node* from = this->lookupNode(i);
      if (from == 0)
        continue;

      const IndexVector& outEdges = this->getOutEdges(from->getIndex());
      seenNeij.clear();
      for (unsigned long j = 0; j < outEdges.size(); j++)
      {
        if (outEdges[j] < 0 || seenNeij.count(outEdges[j]) > 0)
          continue;
        seenNeij.insert(outEdges[j]);

        Node* to = this->getNode(this->getEdge(outEdges[j])->getToNodeIndex());

        mergeParallelLabelFreeEdgesByMass(AbruijnNode::castNodePtr(from),
                                          AbruijnNode::castNodePtr(to));
      }
    }
  }

  bool AbruijnGraph2::computeHeaviestParallelPathDAG()
  {
    compress();

    vector<vector<ParallelPath*> > heaviestPaths(this->numNodes());
    vector<unsigned int> processedOutEdges(this->numNodes(), 0);

    list<Node*> topOrder;
    if (!getTopologicalOrderingDAG(topOrder))
    {
      ERROR_MSG("Graph contains a cycle!!!!");
      return false;
    }

    for (list<Node*>::iterator topIt = topOrder.begin(); topIt
        != topOrder.end(); topIt++)
    {
      vector<ParallelPath*> ppaths;
      AbruijnNode* abNode = AbruijnNode::castNodePtr(*topIt);
      ppaths.push_back(new ParallelPath(*this, *abNode));

      if (AbruijnGraph2::DEBUG_EXPANSION)
      {
        DEBUG_MSG("Processing paths ending at " << abNode->getGraphvizLabel());
      }

      const IndexVector& inEdges = getInEdges((*topIt)->getIndex());
      for (unsigned int i = 0; i < inEdges.size(); i++)
      {
        if (inEdges[i] < 0)
        {
          continue;
        }
        JumpEdge* abEdge = getEdge(inEdges[i]);

        vector<ParallelPath*>& prevPpaths =
            heaviestPaths[abEdge->getFromNodeIndex()];
        for (vector<ParallelPath*>::const_iterator ppIt = prevPpaths.begin(); ppIt
            != prevPpaths.end(); ppIt++)
        {
          ParallelPath* newPath = new ParallelPath(**ppIt);
          newPath->advance(*this, *abEdge);
          ppaths.push_back(newPath);
        }

        processedOutEdges[abEdge->getFromNodeIndex()]++;
        if (this->getNumOutEdges(getNode(abEdge->getFromNodeIndex()))
            == processedOutEdges[abEdge->getFromNodeIndex()])
        {
          for (vector<ParallelPath*>::iterator ppIt = prevPpaths.begin(); ppIt
              != prevPpaths.end(); ppIt++)
          {
            delete *ppIt;
            *ppIt = 0;
          }
          prevPpaths.resize(0);
        }
      }

      sort(ppaths.begin(), ppaths.end(), ParallelPath::WeightComparator);

      vector<ParallelPath*>& rankedPPaths = heaviestPaths[(*topIt)->getIndex()];
      rankedPPaths.resize(0);
      for (unsigned int i = 0; i < ppaths.size(); i++)
      {
        if (rankedPPaths.size() == 0)
        {
          rankedPPaths.push_back(ppaths[i]);
          continue;
        }
        if (rankedPPaths.size() == AbruijnGraph2::NUM_PATHS_PER_NODE
            || rankedPPaths.back()->operator ==(*ppaths[i]))
        {
          delete ppaths[i];
          ppaths[i] = 0;
          continue;
        }
        rankedPPaths.push_back(ppaths[i]);
      }
    }

    ParallelPath* heaviestPath = 0;
    double maxWeight = 0 - numeric_limits<double>::max();
    for (unsigned long i = 0; i < heaviestPaths.size(); i++)
    {
      if (heaviestPaths[i].size() == 0)
      {
        continue;
      }
      if (heaviestPaths[i].front()->getWeight() > maxWeight)
      {
        heaviestPath = heaviestPaths[i].front();
        maxWeight = heaviestPath->getWeight();
      }
    }

    m_consensusPath = *heaviestPath;
    heaviestPath = 0;

    for (unsigned long i = 0; i < heaviestPaths.size(); i++)
    {
      if (heaviestPaths[i].size() == 0)
      {
        continue;
      }
      for (unsigned int j = 0; j < heaviestPaths[i].size(); j++)
      {
        delete heaviestPaths[i][j];
        heaviestPaths[i][j] = 0;
      }
      heaviestPaths[i].resize(0);
    }
    return true;
  }

  /*

   bool AbruijnGraph2::computeHeaviestPathDAG()
   {
   compress();

   Node* start = this->addNode();
   Node* end = this->addNode();

   for (unsigned long i = 0; i <= this->maxNodeIdx(); i++)
   {
   Node* node = this->lookupNode(i);
   if (node == 0)
   continue;
   AbruijnNode* abNode = AbruijnNode::castNodePtr(node);
   if (!abNode->isGreen() && abNode->getIndex() != start->getIndex()
   && abNode->getIndex() != end->getIndex())
   {
   Edge* newEdge = this->addEdge(start, node);
   newEdge = this->addEdge(node, end);
   }
   }

   DEBUG_MSG("Computing heaviest path with sink edge weights ...");
   recoverSourceNodeScores();
   Path completePath;
   pair<bool, double> pathScore = this->getHeaviestPathDAG(start,
   end,
   completePath);

   if ((!pathScore.first) || completePath.size() <= 2)
   {
   this->removeNode(&start);
   this->removeNode(&end);
   return false;
   }

   DEBUG_MSG("Computing heaviest path with source edge weights ...");
   reverseEdgeWeights();
   Path completePathR;
   pair<bool, double> pathScoreR = this->getHeaviestPathDAG(start,
   end,
   completePathR);

   this->removeNode(&start);
   this->removeNode(&end);

   if (pathScoreR.second > pathScore.second + 0.0001)
   {
   DEBUG_MSG("Choosing source edge weights ...");
   pathScore = pathScoreR;
   completePath = completePathR;
   }
   else
   {
   DEBUG_MSG("Choosing sink edge weights ...");
   reverseEdgeWeights();
   }

   m_consensusPath.resize(completePath.size() - 2);
   m_consensusPathVerts.resize(m_consensusPath.size() + 1);
   unsigned int idxUse = 0;
   for (unsigned int i = 1; i < completePath.size() - 1; i++)
   {
   m_consensusPathVerts[idxUse]
   = AbruijnNode::castNodePtr(this->getNode(completePath[i]->getFromNodeIndex()));
   m_consensusPath[idxUse++] = completePath[i];
   }
   m_consensusPathVerts[idxUse]
   = AbruijnNode::castNodePtr(this->getNode(completePath[completePath.size()
   - 2]->getToNodeIndex()));

   //DEBUG_VAR(pathScore.second);
   m_consensusPathScore = pathScore.second;
   return true;
   }
   */

  void AbruijnGraph2::removeSymmetricYNodes()
  {
    list<Node*> nodesToRemove;
    for (unsigned long i = 0; i <= this->maxNodeIdx(); i++)
    {
      Node* node = this->lookupNode(i);
      if (node == 0)
        continue;
      AbruijnNode* abNode = AbruijnNode::castNodePtr(node);
      if (!abNode->isAllYSymmetric())
      {
        continue;
      }
      nodesToRemove.push_back(node);
      for (unsigned int i = 0; i < abNode->size(); i++)
      {
        AssembledPeak* aPeak = &abNode->operator [](i);
        MZRange tightBound = aPeak->getTightBound();
        string specID(aPeak->getSpecID());

        if (m_nodeLookup.count(specID) > 0
            && m_nodeLookup[specID].count(tightBound) > 0
            && m_nodeLookup[specID][tightBound].first == abNode)
        {
          m_nodeLookup[specID].erase(tightBound);
          if (m_nodeLookup[specID].size() == 0)
          {
            m_nodeLookup.erase(specID);
          }
        }
      }
    }

    for (list<Node*>::iterator nIt = nodesToRemove.begin(); nIt
        != nodesToRemove.end(); nIt++)
    {
      this->removeNode(&(*nIt));
    }
  }

  void AbruijnGraph2::removeSourceGreenNodes()
  {
    list<Node*> nodesToRemove;
    for (unsigned long i = 0; i <= this->maxNodeIdx(); i++)
    {
      Node* node = this->lookupNode(i);
      if (node == 0)
        continue;
      AbruijnNode* abNode = AbruijnNode::castNodePtr(node);
      if ((abNode->isGreen()) && ((this->getNumInEdges(node) == 0)
          || (this->getNumOutEdges(node) == 0)))
      {
        nodesToRemove.push_back(node);
        for (unsigned int i = 0; i < abNode->size(); i++)
        {
          AssembledPeak* aPeak = &abNode->operator [](i);
          MZRange tightBound = aPeak->getTightBound();
          string specID(aPeak->getSpecID());

          if (m_nodeLookup.count(specID) > 0
              && m_nodeLookup[specID].count(tightBound) > 0
              && m_nodeLookup[specID][tightBound].first == abNode)
          {
            m_nodeLookup[specID].erase(tightBound);
            if (m_nodeLookup[specID].size() == 0)
            {
              m_nodeLookup.erase(specID);
            }
          }
        }
      }
    }

    for (list<Node*>::iterator nIt = nodesToRemove.begin(); nIt
        != nodesToRemove.end(); nIt++)
    {
      this->removeNode(&(*nIt));
    }
  }

  void AbruijnGraph2::glueNodes(AbruijnNode* node1, AbruijnNode** node2)
  {
    if ((!this->containsNode(node1)) || (!this->containsNode(*node2)))
    {
      ERROR_MSG("Cannot find one of \'" << node1->toString() << "\' or \'" << (*node2)->toString() << "\'");
      abort();
    }
    if (node1->getIndex() == (*node2)->getIndex())
    {
      return;
    }
    if (node1->isCompositeWithOther(**node2))
    {
      ERROR_MSG("Caught composite glue of nodes " << node1->getGraphvizLabel() << " and " << (*node2)->getGraphvizLabel());
      abort();
    }

    if (DEBUG_EXPANSION)
    {
      DEBUG_MSG("Gluing \'" << node1->getGraphvizLabel() << "\' and \'" << (*node2)->getGraphvizLabel() << "\'");
    }

    node1->mergeColors(**node2);

    unsigned int size1 = node1->size();
    unsigned int size2 = (*node2)->size();
    node1->resize(size1 + size2);

    for (unsigned int i = 0; i < size2; i++)
    {
      (*node1)[i + size1] = (**node2)[i];
    }

    const IndexVector& inEdges = this->getInEdges((*node2)->getIndex());
    for (unsigned long i = 0; i < inEdges.size(); i++)
    {
      if (inEdges[i] < 0)
        continue;

      Edge* edge = this->getEdge(inEdges[i]);
      Node* fromNode = this->getNode(edge->getFromNodeIndex());

      JumpEdge* abEdge = JumpEdge::castEdgePtr(this->addEdge(fromNode,
                                                             node1,
                                                             *edge));

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding edge " << abEdge->getFromNodeIndex() << " -> " << abEdge->getToNodeIndex() << " index=" << abEdge->getIndex());
      }
    }

    const IndexVector& outEdges = this->getOutEdges((*node2)->getIndex());
    for (unsigned long i = 0; i < outEdges.size(); i++)
    {
      if (outEdges[i] < 0)
        continue;

      Edge* edge = this->getEdge(outEdges[i]);
      Node* toNode = this->getNode(edge->getToNodeIndex());
      JumpEdge* abEdge = JumpEdge::castEdgePtr(this->addEdge(node1,
                                                             toNode,
                                                             *edge));

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding edge " << abEdge->getFromNodeIndex() << " -> " << abEdge->getToNodeIndex() << " index=" << abEdge->getIndex());
      }
    }

    for (unsigned int i = 0; i < (**node2).size(); i++)
    {
      AssembledPeak* aPeak = &(**node2)[i];
      const string& specID = aPeak->getSpecID();
      m_nodeLookup[specID][aPeak->getTightBound()].first = node1;
    }

    if (m_nodesToExpand.count(*node2) > 0)
    {
      m_nodesToExpand.erase(*node2);
    }
    this->removeNode((Node**)node2);
  }
}
