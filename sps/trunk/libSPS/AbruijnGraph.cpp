/*
 * AbruijnGraph.cpp
 *
 *  Created on: May 1, 2012
 *      Author: aguthals
 */

#include "AbruijnGraph.h"
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
  class AbruijnGraph;

  const unsigned short AbruijnGraph::BIN_VERSION = 1;
  const unsigned short AbruijnGraph::BIN_SUBVERSION = 1;

  const string AbruijnGraph::BIN_VERSION_ID = "AbruijnGraph_binVersion";
  const string AbruijnGraph::BIN_SUBVERSION_ID = "AbruijnGraph_binSubVersion";

  const unsigned int AbruijnGraph::MAX_NUM_JUMPS = 2;
  const unsigned int AbruijnGraph::MAX_NUM_MODS_PER_JUMP = 0;

  const double AbruijnGraph::MAX_JUMP_MASS = 320.0;

  const string AbruijnGraph::STARTPT_SPEC_ID = "startPtSpec";

  const string AbruijnGraph::ENDPT_SPEC_ID = "endPtSpec";

  unsigned int AbruijnGraph::PATH_EXPAND_LIMIT = 5;

  bool AbruijnGraph::DEBUG_EXPANSION = false;

  bool AbruijnGraph::REVERSE_STARS = true;

  bool AbruijnGraph::ENFORCE_B_ENDPTS = true;

  double AbruijnGraph::GAP_PENALTY = 0.30;

  AbruijnGraph* AbruijnGraph::castGraphPtr(BaseGraph* graphPtr)
  {
    AbruijnGraph* temp = dynamic_cast<AbruijnGraph*> (graphPtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast BaseGraph \'" << graphPtr->toString() << "\' to an AbruijnGraph");
      abort();
    }
    return temp;
  }

  AbruijnGraph* AbruijnGraph::castNodePtr(Node* nodePtr)
  {
    AbruijnGraph* temp = dynamic_cast<AbruijnGraph*> (nodePtr);
    if (temp == 0)
    {
      ERROR_MSG("Failed to cast Node \'" << nodePtr->toString() << "\' to an AbruijnGraph");
      abort();
    }
    return temp;
  }

  string AbruijnGraph::nodeSetToString(const set<AbruijnNode*>& inputNodeSet)
  {
    list<long> nodeIDs;
    for (set<AbruijnNode*>::const_iterator nIt = inputNodeSet.begin(); nIt
        != inputNodeSet.end(); nIt++)
    {
      nodeIDs.push_back((*nIt)->getIndex());
    }
    nodeIDs.sort();
    stringstream buff;
    for (list<long>::const_iterator nIt = nodeIDs.begin(); nIt != nodeIDs.end(); nIt++)
    {
      buff << *nIt << ":";
    }
    return buff.str();
  }

  AbruijnGraph::AbruijnGraph() :
    BaseGraph(), m_globalJumps(specnets::AAJumps::getGlobalJumps()),
        m_consensusPath(0), m_consensusPathVerts(0), m_reversed(false),
        m_nodeLookup(), m_nodesToExpand(), m_startPtAssembly(),
        m_consensusPathScore(0)
  {
    initialize();
  }

  AbruijnGraph::AbruijnGraph(const AbruijnGraph& other) :
    BaseGraph(), m_consensusPath(other.m_consensusPath.size()),
        m_consensusPathVerts(other.m_consensusPathVerts.size()),
        m_reversed(other.m_reversed), m_nodeLookup(),
        m_globalJumps(specnets::AAJumps::getGlobalJumps()), m_nodesToExpand(),
        m_startPtAssembly(), m_consensusPathScore(other.m_consensusPathScore)
  {
    this->operator =(other);
  }

  AbruijnGraph::AbruijnGraph(const pair<pair<vector<int> , vector<int> > ,
                                 vector<pair<vector<int> , vector<double> > > >& inAbinfo,
                             const SpecSet& assembledSpecs,
                             const float& peakTol,
                             pair<pair<vector<int> , vector<int> > , vector<
                                 pair<vector<int> , vector<double> > > >* outAbinfo) :
    m_globalJumps(specnets::AAJumps::getGlobalJumps())
  {
    reSequencePPaths(inAbinfo, assembledSpecs, peakTol, outAbinfo);
  }

  void AbruijnGraph::initialize()
  {
    BaseGraph::initialize();
    m_consensusPath.resize(0);
    m_consensusPathVerts.resize(0);
    m_reversed = false;
    m_nodeLookup.clear();
    m_nodesToExpand.clear();
    m_startPtAssembly.clear();
    m_label = "AbruijnGraph";
    m_consensusPathScore = 0;
  }

  void AbruijnGraph::compress(vector<long>* outputNewNodeIdxs,
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

  void AbruijnGraph::reSequencePPaths(const pair<
                                          pair<vector<int> , vector<int> > ,
                                          vector<pair<vector<int> , vector<
                                              double> > > >& abinfo,
                                      const SpecSet& assembledSpecs,
                                      const float& peakTol,
                                      pair<pair<vector<int> , vector<int> > ,
                                          vector<pair<vector<int> , vector<
                                              double> > > >* outAbinfo)
  {
    if (ENFORCE_B_ENDPTS)
    {
      DEBUG_MSG("Re-sequencing with b end-points ...");
      m_reSequencePPaths(abinfo,
                         assembledSpecs,
                         peakTol,
                         REVERSE_STARS,
                         true,
                         outAbinfo);

      DEBUG_VAR(m_consensusPathScore);

      return;
    }

    pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> , vector<
        double> > > > yAbinfo;

    DEBUG_MSG("Re-sequencing with y end-points ...");
    m_reSequencePPaths(abinfo,
                       assembledSpecs,
                       peakTol,
                       REVERSE_STARS,
                       false,
                       &yAbinfo);

    if (m_consensusPath.size() == 0)
    {
      return;
    }

    double yPathScore = m_consensusPathScore;

    AbruijnGraph yGraph(*this);

    DEBUG_MSG("Re-sequencing with b end-points ...");
    m_reSequencePPaths(abinfo,
                       assembledSpecs,
                       peakTol,
                       REVERSE_STARS,
                       true,
                       outAbinfo);

    double bPathScore = m_consensusPathScore;

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

  Node* AbruijnGraph::cloneNode(Node& copyNode) const
  {
    Node* copy = new AbruijnNode(*(AbruijnNode::castNodePtr(&copyNode)));
    return copy;
  }

  Edge* AbruijnGraph::cloneEdge(Edge& copyEdge) const
  {
    Edge* copy = new AbruijnEdge(*(AbruijnEdge::castEdgePtr(&copyEdge)));
    return copy;
  }

  Node* AbruijnGraph::createNode(void) const
  {
    Node* node = new AbruijnNode();
    return node;
  }

  Edge* AbruijnGraph::createEdge(void) const
  {
    Edge* edge = new AbruijnEdge();
    return edge;
  }

  AbruijnGraph & AbruijnGraph::operator=(const AbruijnGraph &other)
  {
    if (this == &other)
    {
      return *this;
    }
    BaseGraph::operator=((const BaseGraph&)other);

    m_internalCopy(other);

    return *this;
  }

  void AbruijnGraph::copy(BaseGraph& otherGraph)
  {
    this->operator =(*castGraphPtr(&otherGraph));
  }

  void AbruijnGraph::addBinaryVersionInfo(map<string, unsigned short>& versions) const
  {
    BaseGraph::addBinaryVersionInfo(versions);
    versions[BIN_VERSION_ID] = BIN_VERSION;
    versions[BIN_SUBVERSION_ID] = BIN_VERSION;
    AbruijnNode nd;
    nd.addBinaryVersionInfo(versions);
    AbruijnEdge eg;
    eg.addBinaryVersionInfo(versions);
  }

  void AbruijnGraph::outputConsensusSpectrum(Spectrum& outputSpec) const
  {
    outputSpec.resize(m_consensusPathVerts.size());
    if (outputSpec.size() == 0)
    {
      WARN_MSG("No consensus present, skipping call to outputConsensusSpectrum!!");
      return;
    }
    unsigned int idxUse = 0;
    map<Node*, unsigned int> nodeToIdx;
    for (unsigned int i = 0; i < m_consensusPathVerts.size(); i++)
    {
      AbruijnNode* abNode = m_consensusPathVerts[i];
      if (DEBUG_EXPANSION)
      {
        DEBUG_VAR(abNode->getGraphvizLabel());
      }
      if (!abNode->haveRealPeak())
      {
        continue;
      }
      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Found peak");
      }
      outputSpec[idxUse][1] = 0;
      nodeToIdx[(Node*)abNode] = idxUse;
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

    for (unsigned int i = 0; i < m_consensusPath.size(); i++)
    {
      AbruijnEdge* abEdge = AbruijnEdge::castEdgePtr(m_consensusPath[i]);
      if (DEBUG_EXPANSION)
      {
        DEBUG_VAR(abEdge->toString());
      }

      Node* nodeTo = this->getNode(abEdge->getToNodeIndex());

      if (nodeToIdx.count(nodeTo) > 0)
      {
        outputSpec[nodeToIdx[nodeTo]][0] = abEdge->getMass() + totMass;
      }
      totMass += abEdge->getMass();
    }
    outputSpec.setCharge(1);
    outputSpec.setParentMass(totMass + specnets::AAJumps::massMH);
    outputSpec.sortPeaks();
  }

  void AbruijnGraph::m_reSequencePPaths(const pair<pair<vector<int> , vector<
                                            int> > , vector<pair<vector<int> ,
                                            vector<double> > > >& inAbinfo,
                                        const SpecSet& assembledSpecs,
                                        const float& peakTol,
                                        const bool reverseStars,
                                        const bool useBEndpts,
                                        pair<pair<vector<int> , vector<int> > ,
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
      newAssembledSpecs[idxUse].setPeakTolerance(peakTol);

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
          list<AbruijnEdge*> addedEdges;
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

    this->expandForwardReversePaths(peakTol);

    if (!computeHeaviestPathDAG())
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

    for (unsigned int i = 0; i < m_consensusPathVerts.size(); i++)
    {
      AbruijnNode* abNode = m_consensusPathVerts[i];
      if (abNode->haveRealPeak())
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

  void AbruijnGraph::m_internalCopy(const AbruijnGraph &other)
  {
    m_startPtAssembly = other.m_startPtAssembly;
    m_reversed = other.m_reversed;
    m_consensusPathScore = other.m_consensusPathScore;

    m_consensusPath.resize(other.m_consensusPath.size());
    for (unsigned int i = 0; i < other.m_consensusPath.size(); i++)
    {
      m_consensusPath[i]
          = BaseGraph::getEdge(other.m_consensusPath[i]->getIndex());
    }

    m_consensusPathVerts.resize(other.m_consensusPathVerts.size());
    for (unsigned int i = 0; i < other.m_consensusPathVerts.size(); i++)
    {
      m_consensusPathVerts[i]
          = getNode(other.m_consensusPathVerts[i]->getIndex());
    }

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

  void AbruijnGraph::addSpectrum(const Spectrum& prmSpec, bool allEndPts)
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

  void AbruijnGraph::addGlues(const SpectrumAlignment& prmAlign)
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

  AbruijnNode* AbruijnGraph::addAbruijnNode(AbruijnNode& copyNode)
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

  AbruijnEdge* AbruijnGraph::addLabelFreeAbEdge(AbruijnNode* from,
                                                AbruijnNode* to,
                                                const double& mass,
                                                const double& weight,
                                                const double& rWeight)
  {
    AbruijnEdge newEdge;
    list<string> tempList;
    newEdge.set(mass, weight, rWeight, "", "", "", false, tempList);

    AbruijnEdge* addedEdge = AbruijnEdge::castEdgePtr(this->addEdge(from,
                                                                    to,
                                                                    newEdge));
    return addedEdge;
  }

  void AbruijnGraph::recoverSourceNodeScores()
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
        if ((*abNode)[j].getSpecID().find(AbruijnGraph::STARTPT_SPEC_ID)
            != string::npos
            || (*abNode)[j].getSpecID().find(AbruijnGraph::ENDPT_SPEC_ID)
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
          AbruijnEdge* abEdge = getEdge(inEdges[j]);
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
          AbruijnEdge* abEdge = getEdge(outEdges[j]);
          if (DEBUG_EXPANSION)
            DEBUG_MSG("Adding " << totalMissedScoreR << " to " << abEdge->toString());

          abEdge->setRWeight(abEdge->getRWeight() + totalMissedScoreR);
        }
      }
    }
  }

  void AbruijnGraph::reverseEdgeWeights()
  {
    for (unsigned long i = 0; i <= this->maxEdgeIdx(); i++)
    {
      Edge* edge = this->lookupEdge(i);
      if (edge == 0)
        continue;

      AbruijnEdge* abEdge = AbruijnEdge::castEdgePtr(edge);
      const double weight = abEdge->getWeight();
      abEdge->setWeight(abEdge->getRWeight());
      abEdge->setRWeight(weight);
    }
  }

  void AbruijnGraph::addEndpointEdges(const set<string>& assembledSpecIDs,
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
    startPtSpec.fileName = AbruijnGraph::STARTPT_SPEC_ID;
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
    endPtSpec.fileName = AbruijnGraph::ENDPT_SPEC_ID;
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

  void AbruijnGraph::addAbruijnEdges(const Spectrum& prmSpec,
                                     const unsigned int& peakIdxFrom,
                                     const unsigned int& peakIdxTo,
                                     const double& avgSpecScore,
                                     list<AbruijnEdge*>* addedEdges)
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
    AbruijnNode* nodeFrom = pairFrom.first;
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
                                   MAX_NUM_JUMPS,
                                   MAX_NUM_MODS_PER_JUMP);

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
      double edgeWeight = prmSpec[peakIdxTo][1] - abs(jumpMass - foundJumpMass)
          - (((double)numJumps) / 10.0);
      double rEdgeWeight = prmSpec[peakIdxFrom][1] - abs(jumpMass
          - foundJumpMass) - (((double)numJumps) / 10.0);

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
      AbruijnEdge newEdge;
      newEdge.set(foundJumpMass,
                  edgeWeight,
                  rEdgeWeight,
                  jIt->first,
                  noModJump,
                  reverseString(noModJump),
                  modified,
                  tempList);

      AbruijnEdge* addedEdge = AbruijnEdge::castEdgePtr(this->addEdge(nodeFrom,
                                                                      nodeTo,
                                                                      newEdge));
      pairTo.second[0] = true;
      pairFrom.second[1] = true;

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding edge " << nodeFrom->getIndex() << " -> " << nodeTo->getIndex() << " \'" << addedEdge->getFLabel() << "\' index=" << addedEdge->getIndex() << ", jump=" << foundJumpMass << ", weight=" << addedEdge->getWeight());
      }

      if (addedEdges != 0)
      {
        addedEdges->push_back(addedEdge);
      }
    }
  }

  void AbruijnGraph::expandForwardReversePaths(const float& peakTol)
  {
    DEBUG_MSG("Merging edges ...");
    mergeAllParallelLabelFreeEdgesByMass(peakTol);
    mergeAllParallelEdges();

    DEBUG_VAR(this->numNodes());
    DEBUG_VAR(this->numEdges());

    DEBUG_MSG("Merging forward paths ...");
    processAllPaths(false, false);
    DEBUG_MSG("Merging reverse paths ...");
    processAllPaths(true, false);

    DEBUG_VAR(this->numNodes());
    DEBUG_VAR(this->numEdges());

    DEBUG_MSG("Expanding forward paths ...");
    processAllPaths(false, true);
    //DEBUG_MSG("Merging edges ...");
    //mergeAllParallelLabelFreeEdgesByMass(peakTol);
    //mergeAllParallelEdges();
    DEBUG_MSG("Expanding reverse paths ...");
    processAllPaths(true, true);

    DEBUG_VAR(this->numNodes());
    DEBUG_VAR(this->numEdges());

    for (unsigned long i = 0; i <= this->maxEdgeIdx(); i++)
    {
      Edge* edge = this->lookupEdge(i);
      if (edge == 0)
      {
        continue;
      }
      AbruijnEdge* abEdge = AbruijnEdge::castEdgePtr(edge);
      if (abEdge->expandOnly())
      {
        this->removeEdge(&edge);
      }
    }
    //DEBUG_MSG("Merging edges ...");
    //mergeAllParallelLabelFreeEdgesByMass(peakTol);
    //mergeAllParallelEdges();
  }

  typedef sps::tuple<AbruijnNode*, list<AbruijnEdge*> , unsigned int>
      PathSource;

  void AbruijnGraph::expandPaths(AbruijnNode* start, const bool& goReverse)
  {
    list<PathSource> nodeQueue;
    PathSource nextSource;
    nextSource.m0 = start;
    nextSource.m1.clear();
    nextSource.m2 = 0;

    nodeQueue.push_front(nextSource);

    if (DEBUG_EXPANSION)
      DEBUG_VAR(start->getIndex());

    tr1::unordered_set<string> outgoingLabels;
    list<AbruijnEdge*> prefixEdges;
    list<AbruijnEdge*> filteredPrefixEdges;
    set<AbruijnNode*> mergeNodeSet;
    set<AbruijnNode*> destNodes;
    list<AbruijnEdge*> destEdges;

    tr1::unordered_map<string, AbruijnNode*> createdNodes;
    list<AbruijnNode*> nodesToMerge;

    TrieMap<AbruijnEdge*> prefixMap;
    AbruijnEdge mergeEdge;

    while (nodeQueue.size() > 0)
    {
      //DEBUG_TRACE;
      nextSource = nodeQueue.front();
      //DEBUG_TRACE;
      nodeQueue.pop_front();
      //DEBUG_TRACE;

      AbruijnNode* source = nextSource.m0;
      const string& pathForward = source->getFPath();
      const string& pathReverse = source->getRPath();
      set<AbruijnEdge*> edgesToRemove(nextSource.m1.begin(),
                                      nextSource.m1.end());
      const unsigned int pathLength = nextSource.m2;

      const string& pathUse = (goReverse) ? pathReverse : pathForward;
      const string& pathDontUse = (goReverse) ? pathForward : pathReverse;

      if (DEBUG_EXPANSION)
      {
        DEBUG_VAR(source->getGraphvizLabel());
        DEBUG_MSG("pathUse = " << pathUse << ", pathDontUse = " << pathDontUse << ", pathLength = " << pathLength);
      }

      outgoingLabels.clear();
      outgoingLabels.rehash((goReverse) ? this->getNumInEdges(source)
          : this->getNumOutEdges(source));

      const IndexVector& outEdges = (goReverse) ? this->getInEdges(source)
          : this->getOutEdges(source);

      set<AbruijnNode*> skipNodes;

      //DEBUG_TRACE;

      for (unsigned long i = 0; i < outEdges.size(); i++)
      {
        if (outEdges[i] < 0)
        {
          continue;
        }
        AbruijnEdge* abEdge =
            AbruijnEdge::castEdgePtr(this->getEdge(outEdges[i]));

        // do not expand label free edges
        const string& jumpLabel = abEdge->getJumpLabel();
        const unsigned int jumpLen = jumpLabel.length();
        const string& label = (goReverse) ? abEdge->getRLabel()
            : abEdge->getFLabel();
        const unsigned int extPathLen = pathLength
            + specnets::AAJumps::getNumJumps(label);

        if (DEBUG_EXPANSION)
          DEBUG_MSG("Considering - extPathLen = " << extPathLen << ", " << abEdge->toString());

        if (jumpLen == 0)
        {
          continue;
        }
        if (extPathLen > PATH_EXPAND_LIMIT)
        {
          skipNodes.insert((goReverse) ? getNode(abEdge->getFromNodeIndex())
              : getNode(abEdge->getToNodeIndex()));
          continue;
        }

        const string revExt = (goReverse) ? abEdge->getFLabel().substr(jumpLen)
            : abEdge->getRLabel().substr(jumpLen);

        if ((!source->isGreen()) || (revExt == pathDontUse))
        {
          if (DEBUG_EXPANSION)
            DEBUG_MSG("Indexing - revExt = " << revExt << ", " << abEdge->toString());
          prefixMap.insertWord(label, &abEdge);
          if (isPrefix(pathUse, label))
          {
            outgoingLabels.insert(label);
          }
        }
      }

      createdNodes.rehash(outgoingLabels.size());

      for (tr1::unordered_set<string>::const_iterator labIt =
          outgoingLabels.begin(); labIt != outgoingLabels.end(); labIt++)
      {
        prefixMap.findPrefixes(*labIt, prefixEdges);

        /*if (source->getIndex() == 315 || source->getIndex() == 316)
         {
         DEBUG_MSG("Evaluating \'" << *labIt << "\', numPrefixes=" << prefixEdges.size());
         }*/
        if (prefixEdges.size() <= 1
            || mergeEdge.mergeEdges(prefixEdges,
                                    goReverse,
                                    &filteredPrefixEdges) <= 1)
        {
          continue;
        }
        filteredPrefixEdges = prefixEdges;

        unsigned int minPathLength = AbruijnGraph::MAX_NUM_JUMPS;//numeric_limits<unsigned int>::max();
        set<unsigned int> jumpLengths;
        mergeNodeSet.clear();
        for (list<AbruijnEdge*>::iterator eIt = filteredPrefixEdges.begin(); eIt
            != filteredPrefixEdges.end(); eIt++)
        {
          if (DEBUG_EXPANSION)
          {
            DEBUG_MSG("Prefix " << (*eIt)->toString());
          }
          jumpLengths.insert((unsigned int)(*eIt)->getJumpLabel().length());

          AbruijnNode* dest = (goReverse) ? getNode((*eIt)->getFromNodeIndex())
              : getNode((*eIt)->getToNodeIndex());
          mergeNodeSet.insert(dest);
          minPathLength = min(minPathLength,
                              (unsigned int)(*eIt)->getJumpLabel().length());
        }

        const string nodeSetID = nodeSetToString(mergeNodeSet);
        if (createdNodes.count(nodeSetID) > 0 && jumpLengths.size() == 1)
        {
          AbruijnNode* toNode = createdNodes[nodeSetID];
          Edge* addedEdge = (goReverse) ? this->addEdge(toNode,
                                                        source,
                                                        mergeEdge)
              : this->addEdge(source, toNode, mergeEdge);

          if (goReverse)
          {
            toNode->setFPath("");
          }
          else
          {
            toNode->setRPath("");
          }

          for (list<AbruijnEdge*>::iterator eIt = filteredPrefixEdges.begin(); eIt
              != filteredPrefixEdges.end(); eIt++)
          {
            const string prefixLabel = (goReverse) ? (*eIt)->getRLabel()
                : (*eIt)->getFLabel();
            if (prefixLabel == *labIt)
            {
              edgesToRemove.insert(*eIt);
            }
          }

          if (DEBUG_EXPANSION)
          {
            DEBUG_MSG("Adding " << addedEdge->toString());
          }
          continue;
        }

        for (list<AbruijnEdge*>::iterator eIt = filteredPrefixEdges.begin(); eIt
            != filteredPrefixEdges.end(); eIt++)
        {
          const string& prefixLabel = (goReverse) ? (*eIt)->getRLabel()
              : (*eIt)->getFLabel();

          AbruijnNode* toNode = (goReverse)
              ? getNode((*eIt)->getFromNodeIndex())
              : getNode((*eIt)->getToNodeIndex());

          if (specnets::AAJumps::getNumJumps((*eIt)->getJumpLabel())
              > minPathLength)
          {
            this->injectExpandedPath(*eIt, minPathLength, goReverse);

            AbruijnNode* newToNode = (goReverse)
                ? getNode((*eIt)->getFromNodeIndex())
                : getNode((*eIt)->getToNodeIndex());
          }
        }

        const string nextUsePath = labIt->substr(minPathLength);
        //nextUsePath.erase(0, minPathLength);
        const string nextDontUsePath =
            specnets::AAJumps::reversePeptide(labIt->substr(0, minPathLength))
                + pathDontUse;

        const string& nextFPath = (goReverse) ? nextDontUsePath : nextUsePath;
        const string& nextRPath = (goReverse) ? nextUsePath : nextDontUsePath;

        //DEBUG_VAR(*labIt);
        //DEBUG_VAR(mergeNodeSet.size());

        //DEBUG_VAR(nextFPath);
        //DEBUG_VAR(nextRPath);

        // Edges that need to be removed following expansion of the merged node
        list<AbruijnEdge*> edgesToRemoveNext;

        AbruijnNode* addedNode = mergeNodes(filteredPrefixEdges,
                                            source,
                                            nextFPath,
                                            nextRPath,
                                            *labIt,
                                            edgesToRemove,
                                            edgesToRemoveNext,
                                            goReverse);

        if (addedNode == 0)
        {
          continue;
        }

        if (!addedNode->isGreen())
        {
          createdNodes[nodeSetID] = addedNode;
        }

        nextSource.m0 = addedNode;
        nextSource.m1.clear();
        nextSource.m1.insert(nextSource.m1.begin(),
                             edgesToRemoveNext.begin(),
                             edgesToRemoveNext.end());

        nextSource.m2 = pathLength + minPathLength;

        if (DEBUG_EXPANSION)
          DEBUG_MSG("Continuing to " << addedNode->getGraphvizLabel() << " with \'" << nextUsePath << "\' / \'" << nextDontUsePath << "\'");

        nodeQueue.push_front(nextSource);
      }
      prefixMap.clear();
      createdNodes.clear();

      for (set<AbruijnEdge*>::iterator eIt = edgesToRemove.begin(); eIt
          != edgesToRemove.end(); eIt++)
      {
        this->removeExpandedPath((AbruijnEdge**)&(*eIt));
      }

      for (set<AbruijnNode*>::iterator nIt = skipNodes.begin(); nIt
          != skipNodes.end(); nIt++)
      {
        if (goReverse)
          mergeParallelEdges(*nIt, source);
        else
          mergeParallelEdges(source, *nIt);
      }
      skipNodes.clear();
    }
  }

  void AbruijnGraph::mergePaths(AbruijnNode* start, const bool& goReverse)
  {
    list<PathSource> nodeQueue;
    PathSource nextSource;
    nextSource.m0 = start;
    nextSource.m1.clear();
    nextSource.m2 = 0;

    nodeQueue.push_front(nextSource);

    if (DEBUG_EXPANSION)
      DEBUG_VAR(start->getIndex());

    TrieMap<AbruijnEdge*> prefixMap;
    tr1::unordered_set<string> outgoingLabels;
    list<AbruijnEdge*> prefixEdges;
    set<AbruijnNode*> mergedNodes;
    tr1::unordered_set<string> mergedSpecIDs;
    list<list<AbruijnNode*> > nodeSetsToMerge;

    while (nodeQueue.size() > 0)
    {
      //DEBUG_TRACE;
      AbruijnNode* source = nodeQueue.front().m0;
      //DEBUG_TRACE;
      nodeQueue.pop_front();
      //DEBUG_TRACE;

      if (m_nodesToExpand.count(source) == 0)
      {
        continue;
      }
      else
      {
        m_nodesToExpand.erase(source);
      }

      if (DEBUG_EXPANSION)
      {
        DEBUG_VAR(source->getGraphvizLabel());
      }

      mergedNodes.clear();
      nodeSetsToMerge.clear();

      const IndexVector& outEdges = (goReverse) ? this->getInEdges(source)
          : this->getOutEdges(source);
      outgoingLabels.rehash(outEdges.size());

      for (unsigned long i = 0; i < outEdges.size(); i++)
      {
        if (outEdges[i] < 0)
        {
          continue;
        }
        AbruijnEdge* abEdge =
            AbruijnEdge::castEdgePtr(this->getEdge(outEdges[i]));

        // do not expand label free edges
        const unsigned int jumpLen = abEdge->getJumpLabel().length();
        const string& label = (goReverse) ? abEdge->getRLabel()
            : abEdge->getFLabel();

        if (jumpLen == 0)
        {
          continue;
        }

        if (DEBUG_EXPANSION)
          DEBUG_MSG("Indexing - " << abEdge->toString());

        prefixMap.insertWord(label, &abEdge);
        outgoingLabels.insert(label);
      }

      for (tr1::unordered_set<string>::const_iterator labIt =
          outgoingLabels.begin(); labIt != outgoingLabels.end(); labIt++)
      {
        if (DEBUG_EXPANSION)
        {
          DEBUG_MSG("Checking " << *labIt << " ...");
        }

        if (prefixMap.findWord(*labIt, &prefixEdges) <= 1)
        {
          continue;
        }

        list<AbruijnNode*> nodesToMerge;
        mergedSpecIDs.clear();
        mergedSpecIDs.rehash(prefixEdges.size());
        prefixEdges.sort(AbruijnEdge::MergeSizeCompare(goReverse));

        for (list<AbruijnEdge*>::iterator eIt = prefixEdges.begin(); eIt
            != prefixEdges.end(); eIt++)
        {
          if (DEBUG_EXPANSION)
          {
            DEBUG_MSG("Prefix " << (*eIt)->toString());
          }
          AbruijnNode* dest = (goReverse) ? getNode((*eIt)->getFromNodeIndex())
              : getNode((*eIt)->getToNodeIndex());

          if (mergedNodes.count(dest) > 0)
          {
            continue;
          }

          bool foundSameSpec = false;
          for (unsigned int i = 0; i < dest->size(); i++)
          {
            const string& specID = (*dest)[i].getSpecID();
            if (mergedSpecIDs.count(specID) > 0)
            {
              foundSameSpec = true;
              break;
            }
            mergedSpecIDs.insert(specID);
          }
          if (foundSameSpec)
          {
            continue;
          }
          mergedNodes.insert(dest);
          nodesToMerge.push_back(dest);
        }

        if (nodesToMerge.size() > 0)
        {
          nodeSetsToMerge.push_back(nodesToMerge);
        }
      }
      prefixMap.clear();
      outgoingLabels.clear();

      for (list<list<AbruijnNode*> >::iterator setIt = nodeSetsToMerge.begin(); setIt
          != nodeSetsToMerge.end(); setIt++)
      {
        list<AbruijnNode*>::iterator nIt = setIt->begin();
        AbruijnNode* nextDest = *nIt;
        for (nIt++; nIt != setIt->end(); nIt++)
        {
          glueNodes(nextDest, &(*nIt));
        }

        if (goReverse)
        {
          mergeParallelEdges(nextDest, source);
        }
        else
        {
          mergeParallelEdges(source, nextDest);
        }

        if (setIt->size() > 1)
        {
          if (DEBUG_EXPANSION)
          {
            DEBUG_MSG("Continuing to " << nextDest->getGraphvizLabel());
          }

          if (m_nodesToExpand.count(nextDest) > 0)
          {
            m_nodesToExpand.erase(nextDest);
          }

          nextSource.m0 = nextDest;
          nodeQueue.push_front(nextSource);
        }
      }
    }
  }

  class AbDistCompare
  {
    bool m_reverse;
  public:
    AbDistCompare(bool sortReverse = false) :
      m_reverse(sortReverse)
    {
    }
    bool operator()(const pair<double, AbruijnNode*>& lhs, const pair<double,
        AbruijnNode*>& rhs) const
    {
      return (m_reverse) ? (lhs.first >= rhs.first) : (lhs.first <= rhs.first);
    }
  };

  void AbruijnGraph::processAllPaths(const bool& goReverse,
                                     const bool& expandYes)
  {
    if (this->numNodes() == 0)
    {
      WARN_MSG("Calling processAllPaths on an empty graph!!");
      return;
    }
    /**
     * To minimize running time, need to expand nodes furthest away from mass 0 first
     */
    //this->compress();

    m_nodesToExpand.clear();

    list<pair<double, AbruijnNode*> > rankedDistNodes;
    for (unsigned long nodeIdx = 0; nodeIdx <= this->maxNodeIdx(); nodeIdx++)
    {
      Node* node = lookupNode(nodeIdx);
      if (node == 0)
      {
        continue;
      }
      AbruijnNode* abNode = AbruijnNode::castNodePtr(node);
      m_nodesToExpand.insert(abNode);
      for (unsigned int i = 0; i < abNode->size(); i++)
      {
        const string& specID = (*abNode)[i].getSpecID();
        if (specID.find(AbruijnGraph::STARTPT_SPEC_ID) != string::npos
            || specID.find(AbruijnGraph::ENDPT_SPEC_ID) != string::npos)
        {
          continue;
        }
        if (m_startPtAssembly.count(specID) == 0)
        {
          ERROR_MSG("Could not find spectrum " << specID << " in start-pt assembly!");
          abort();
        }
        double startDist = (*abNode)[i].getMass() + m_startPtAssembly[specID];
        rankedDistNodes.push_back(pair<double, AbruijnNode*> (startDist, abNode));
        break;
      }
    }

    rankedDistNodes.sort(AbDistCompare(!goReverse));

    // expand paths from all nodes
    for (list<pair<double, AbruijnNode*> >::iterator nIt =
        rankedDistNodes.begin(); nIt != rankedDistNodes.end(); nIt++)
    {
      if (m_nodesToExpand.count(nIt->second) == 0)
        continue;

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Calling expand/merge on " << nIt->second->getIndex() << " with distance " << nIt->first);
      }

      if (expandYes)
        expandPaths(nIt->second, goReverse);
      else
        mergePaths(nIt->second, goReverse);
    }
    m_nodesToExpand.clear();
    /*for (unsigned int i = 0; i < MAX_NUM_JUMPS - 1; i++)
     {
     removeSourceGreenNodes();
     }*/
  }

  void AbruijnGraph::mergeParallelEdges(AbruijnNode* from,
                                        AbruijnNode* to,
                                        bool mergeMods)
  {
    map<string, list<AbruijnEdge*> > parallelSets;
    set<string> labelsToMerge;

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
      AbruijnEdge* abEdge = AbruijnEdge::castEdgePtr(edge);
      string jumpLabel = abEdge->getJumpLabel();

      if (jumpLabel.length() == 0)
      {
        continue;
      }

      if (mergeMods)
      {
        jumpLabel = specnets::AAJumps::stripMods(jumpLabel);
      }

      if (parallelSets.count(jumpLabel) > 0)
      {
        parallelSets[jumpLabel].push_back(abEdge);
        labelsToMerge.insert(jumpLabel);
      }
      else
      {
        list<AbruijnEdge*> tempList;
        tempList.push_back(abEdge);
        parallelSets[jumpLabel] = tempList;
      }
    }

    if (labelsToMerge.size() == 0)
    {
      return;
    }

    list<AbruijnEdge*> edgesToRemove;
    list<AbruijnEdge> edgesToAdd;
    for (set<string>::const_iterator lIt = labelsToMerge.begin(); lIt
        != labelsToMerge.end(); lIt++)
    {
      list<AbruijnEdge*>& edgeList = parallelSets[*lIt];
      if (edgeList.size() <= 1)
        continue;

      AbruijnEdge templateEdge;
      unsigned int numMerged = templateEdge.mergeEdges(edgeList, false);
      edgesToRemove.insert(edgesToRemove.end(),
                           edgeList.begin(),
                           edgeList.end());

      templateEdge.setJumpLabel(*lIt);
      edgesToAdd.push_back(templateEdge);
    }

    for (list<AbruijnEdge*>::iterator eIt = edgesToRemove.begin(); eIt
        != edgesToRemove.end(); eIt++)
    {
      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Removing edge " << (*eIt)->getFromNodeIndex() << " -> " << (*eIt)->getToNodeIndex() << " \'" << (*eIt)->getFLabel() << "\' index=" << (*eIt)->getIndex() << ", jump=" << (*eIt)->getMass() << ", weight=" << (*eIt)->getWeight());
      }
      this->removeEdge((Edge**)(&(*eIt)));
    }

    for (list<AbruijnEdge>::iterator eIt = edgesToAdd.begin(); eIt
        != edgesToAdd.end(); eIt++)
    {
      AbruijnEdge* addedEdge = AbruijnEdge::castEdgePtr(this->addEdge(from,
                                                                      to,
                                                                      *eIt));

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding edge " << from->getIndex() << " -> " << to->getIndex() << " \'" << addedEdge->getFLabel() << "\' index=" << addedEdge->getIndex() << ", jump=" << addedEdge->getMass() << ", weight=" << addedEdge->getWeight());
      }
    }
  }

  void AbruijnGraph::mergeAllParallelEdges(bool mergeMods)
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

  void AbruijnGraph::pruneParallelEdgesByMass(AbruijnNode* from,
                                              AbruijnNode* to,
                                              float peakTol)
  {
    map<MZRange, AbruijnEdge*> heavEdges;
    list<AbruijnEdge*> edgesToRemove;

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
      AbruijnEdge* abEdge = AbruijnEdge::castEdgePtr(edge);
      float jumpMass = abEdge->getMass();
      double jumpWeight = abEdge->getWeight();

      MZRange jump(jumpMass, jumpWeight, peakTol / 2.0);

      map<MZRange, AbruijnEdge*>::iterator findIt = heavEdges.find(jump);

      if (findIt == heavEdges.end())
      {
        heavEdges[jump] = abEdge;
      }
      else if (findIt->first.getIntensity() < jumpWeight)
      {
        edgesToRemove.push_back(findIt->second);
        heavEdges.erase(jump);
        heavEdges[jump] = abEdge;
      }
      else
      {
        edgesToRemove.push_back(abEdge);
      }
    }

    if (edgesToRemove.size() == 0)
    {
      return;
    }

    for (list<AbruijnEdge*>::iterator eIt = edgesToRemove.begin(); eIt
        != edgesToRemove.end(); eIt++)
    {
      this->removeEdge((Edge**)(&(*eIt)));
    }
  }

  void AbruijnGraph::pruneAllParallelEdgesByMass(float peakTol)
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

        pruneParallelEdgesByMass(AbruijnNode::castNodePtr(from),
                                 AbruijnNode::castNodePtr(to),
                                 peakTol);
      }
    }
  }

  void AbruijnGraph::mergeParallelLabelFreeEdgesByMass(AbruijnNode* from,
                                                       AbruijnNode* to,
                                                       float peakTol)
  {
    map<MZRange, list<AbruijnEdge*> > parallelSets;
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
      AbruijnEdge* abEdge = AbruijnEdge::castEdgePtr(edge);
      string jumpLabel = abEdge->getJumpLabel();
      float jumpMass = abEdge->getMass();
      double jumpWeight = abEdge->getWeight();

      if (jumpLabel.length() > 0)
      {
        continue;
      }
      //DEBUG_MSG(from->getIndex() << " -> " << to->getIndex() << " : jumpMass = " << jumpMass);
      MZRange jump(jumpMass, jumpWeight, peakTol / 2.0);

      if (parallelSets.count(jump) > 0)
      {
        parallelSets[jump].push_back(abEdge);
        jumpsToMerge.insert(jump);
      }
      else
      {
        list<AbruijnEdge*> tempList;
        tempList.push_back(abEdge);
        parallelSets[jump] = tempList;
      }
      //DEBUG_TRACE;
    }

    if (jumpsToMerge.size() == 0)
    {
      return;
    }
    list<AbruijnEdge*> edgesToRemove;
    for (set<MZRange>::const_iterator lIt = jumpsToMerge.begin(); lIt
        != jumpsToMerge.end(); lIt++)
    {
      list<AbruijnEdge*>& edgeList = parallelSets[*lIt];
      AbruijnEdge templateEdge(*edgeList.front());
      double totalWeight = 0;
      double totalRWeight = 0;

      for (list<AbruijnEdge*>::iterator eIt = edgeList.begin(); eIt
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
    for (list<AbruijnEdge*>::iterator eIt = edgesToRemove.begin(); eIt
        != edgesToRemove.end(); eIt++)
    {
      this->removeEdge((Edge**)(&(*eIt)));
    }
  }

  void AbruijnGraph::mergeAllParallelLabelFreeEdgesByMass(float peakTol)
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
                                          AbruijnNode::castNodePtr(to),
                                          peakTol);
      }
    }
  }

  bool AbruijnGraph::computeHeaviestPathDAG()
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

  void AbruijnGraph::removeExpandedPath(AbruijnEdge** pathPartEdge)
  {
    AbruijnNode* from = getNode((*pathPartEdge)->getFromNodeIndex());
    AbruijnNode* to = getNode((*pathPartEdge)->getToNodeIndex());
    AbruijnNode* nextNode;

    if (DEBUG_EXPANSION)
    {
      DEBUG_MSG("Removing " << (*pathPartEdge)->toString());
    }
    this->removeEdge((Edge**)pathPartEdge);

    list<AbruijnNode*> nodesToCheck;
    nodesToCheck.push_back(from);
    while (nodesToCheck.size() > 0)
    {
      nextNode = nodesToCheck.front();
      nodesToCheck.pop_front();

      if ((!nextNode->isGreen()) || this->getNumOutEdges(nextNode) > 0)
      {
        continue;
      }
      const IndexVector& inEdges = this->getInEdges(nextNode);
      for (unsigned int i = 0; i < inEdges.size(); i++)
      {
        if (inEdges[i] < 0)
          continue;

        nodesToCheck.push_back(getNode(getEdge(inEdges[i])->getFromNodeIndex()));
      }
      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Removing " << nextNode->getGraphvizLabel());
      }
      if (m_nodesToExpand.count(nextNode) > 0)
      {
        m_nodesToExpand.erase(nextNode);
      }
      this->removeNode((Node**)&nextNode);
    }

    nodesToCheck.clear();
    nodesToCheck.push_back(to);
    while (nodesToCheck.size() > 0)
    {
      nextNode = nodesToCheck.front();
      nodesToCheck.pop_front();

      if ((!nextNode->isGreen()) || this->getNumInEdges(nextNode) > 0)
      {
        continue;
      }
      const IndexVector& outEdges = this->getOutEdges(nextNode);
      for (unsigned int i = 0; i < outEdges.size(); i++)
      {
        if (outEdges[i] < 0)
          continue;

        nodesToCheck.push_back(getNode(getEdge(outEdges[i])->getToNodeIndex()));
      }
      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Removing " << nextNode->getGraphvizLabel());
      }
      if (m_nodesToExpand.count(nextNode) > 0)
      {
        m_nodesToExpand.erase(nextNode);
      }
      this->removeNode((Node**)&nextNode);
    }
  }

  void AbruijnGraph::injectExpandedPath(AbruijnEdge* templateEdge,
                                        const unsigned int& step,
                                        const bool& goReverse)
  {
    const string& jumpLabel = templateEdge->getJumpLabel();
    unsigned int numJumps = specnets::AAJumps::getNumJumps(jumpLabel);

    if (numJumps <= 1)
    {
      return;
    }

    //DEBUG_MSG("Expanding " << (*templateEdge)->toString());

    AbruijnNode* from = getNode(templateEdge->getFromNodeIndex());
    AbruijnNode* to = getNode(templateEdge->getToNodeIndex());
    bool modified = false;//jumpLabel.find('(') != string::npos;

    unsigned int jumpLen = jumpLabel.length();

    const string leftJump = (goReverse) ? jumpLabel.substr(0, jumpLen - step)
        : jumpLabel.substr(0, step);
    const string rightJump = (goReverse) ? jumpLabel.substr(jumpLen - step,
                                                            step)
        : jumpLabel.substr(step, jumpLen - step);

    const string& leftFJump = templateEdge->getFLabel();
    const string rightFJump = leftFJump.substr(leftJump.length());

    const string& rightRJump = templateEdge->getRLabel();
    const string leftRJump = rightRJump.substr(rightJump.length());

    const double partWeight = templateEdge->getWeight() / ((double)numJumps);
    const double remainingWeight = templateEdge->getWeight() - partWeight;
    const double partRWeight = templateEdge->getRWeight() / ((double)numJumps);
    const double remainingRWeight = templateEdge->getRWeight() - partRWeight;

    AbruijnNode* newNode = AbruijnNode::castNodePtr(this->addNode());
    newNode->setFPath(rightFJump);
    newNode->setRPath(leftRJump);

    AbruijnEdge* nextEdge;

    if (goReverse)
    {
      this->moveEdge(templateEdge->getIndex(),
                     newNode->getIndex(),
                     to->getIndex());
      nextEdge = AbruijnEdge::castEdgePtr(addEdge(from, newNode));

      nextEdge->setJumpLabel(leftJump);
      templateEdge->setJumpLabel(rightJump);
      nextEdge->setFLabel(leftFJump);
      templateEdge->setFLabel(rightFJump);
      nextEdge->setRLabel(leftRJump);
      templateEdge->setRLabel(rightRJump);
    }
    else
    {
      this->moveEdge(templateEdge->getIndex(),
                     from->getIndex(),
                     newNode->getIndex());
      nextEdge = AbruijnEdge::castEdgePtr(addEdge(newNode, to));

      nextEdge->setJumpLabel(rightJump);
      templateEdge->setJumpLabel(leftJump);
      nextEdge->setFLabel(rightFJump);
      templateEdge->setFLabel(leftFJump);
      nextEdge->setRLabel(rightRJump);
      templateEdge->setRLabel(leftRJump);
    }

    nextEdge->setModified(modified);
    nextEdge->setWeight(remainingWeight);
    templateEdge->setWeight(partWeight);
    nextEdge->setRWeight(remainingRWeight);
    templateEdge->setRWeight(partRWeight);
    nextEdge->setMass(m_globalJumps.getModPeptideMass(nextEdge->getJumpLabel()));
    templateEdge->setMass(m_globalJumps.getModPeptideMass(templateEdge->getJumpLabel()));

    const list<string>& specIDs = templateEdge->getSpecIDs();
    nextEdge->loadIDs(specIDs.begin(), specIDs.end());

    map<string, double> basePeakMasses;
    double leftMass = m_globalJumps.getModPeptideMass(leftJump);

    for (unsigned int i = 0; i < from->size(); i++)
    {
      basePeakMasses[(*from)[i].getSpecID()] = (*from)[i].getMass();
    }

    for (list<string>::const_iterator sIt = specIDs.begin(); sIt
        != specIDs.end(); sIt++)
    {
      if (basePeakMasses.count(*sIt) == 0)
      {
        ERROR_MSG("Could not find spectrum " << *sIt << " in origin node " << from->getGraphvizLabel());
        DEBUG_VAR(templateEdge->getIndex());
        abort();
      }
      double pkMass = basePeakMasses[*sIt] + leftMass;
      newNode->addEmptyPeak(*sIt, pkMass);
    }
    if (DEBUG_EXPANSION)
    {
      DEBUG_MSG("Moved " << templateEdge->toString());
      DEBUG_MSG("Added " << nextEdge->toString());
    }
  }
  /*
   void AbruijnGraph::getEndPointSpectra(map<string, MZRange>& leftAlignedShifts,
   Spectrum& outputStartPtSpec,
   Spectrum& outputEndPtSpec,
   SpectrumAlignmentSet& outputStartPtAligns,
   SpectrumAlignmentSet& outputEndPtAligns)
   {
   outputStartPtSpec.resize(0);
   outputEndPtSpec.resize(0);
   leftAlignedShifts.clear();

   outputStartPtAligns.resize(0);
   outputEndPtAligns.resize(0);

   map<string, map<string, pair<MZRange, MZRange> > > assembly;
   map<string, pair<MZRange, MZRange> > subAssembly;
   map<string, pair<MZRange, MZRange> > refEndPts;
   pair<MZRange, MZRange> shift;

   list<pair<MZRange, MZRange> > matchedPeaksB;
   list<pair<MZRange, MZRange> > matchedPeaksY;
   for (unsigned int i = 0; i < m_alignments.size(); i++)
   {
   const string& spec1ID = m_alignments[i].getSpec1ID();
   const string& spec2ID = m_alignments[i].getSpec2ID();

   Spectrum* spec1 = m_assembledSpecs.getIndex(spec1ID);
   Spectrum* spec2 = m_assembledSpecs.getIndex(spec2ID);

   m_alignments[i].outputMatchedPeaks(matchedPeaksB, matchedPeaksY);

   MZRange& spec1F = matchedPeaksB.front().first;
   MZRange& spec1B = matchedPeaksB.back().first;
   MZRange& spec2F = matchedPeaksB.front().second;
   MZRange& spec2B = matchedPeaksB.back().second;

   int endIdx1 = spec1->findPeaks(spec1->parentMass
   - specnets::AAJumps::massMH);
   if (endIdx1 < 0)
   {
   ERROR_MSG("Failed to locate endpoint for spectrum " << spec1ID);
   abort();
   }
   int endIdx2 = spec2->findPeaks(spec2->parentMass
   - specnets::AAJumps::massMH);
   if (endIdx2 < 0)
   {
   ERROR_MSG("Failed to locate endpoint for spectrum " << spec2ID);
   abort();
   }

   int startIdx1 = spec1->findPeaks(0);
   if (startIdx1 < 0)
   {
   ERROR_MSG("Failed to locate endpoint for spectrum " << spec1ID);
   abort();
   }
   int startIdx2 = spec2->findPeaks(0);
   if (startIdx2 < 0)
   {
   ERROR_MSG("Failed to locate endpoint for spectrum " << spec2ID);
   abort();
   }

   if (refEndPts.count(spec1ID) == 0)
   {
   shift.first = spec1->getPeak(startIdx1);
   shift.second = spec1->getPeak(endIdx1);
   refEndPts[spec1ID] = shift;
   }

   if (refEndPts.count(spec2ID) == 0)
   {
   shift.first = spec1->getPeak(startIdx2);
   shift.second = spec1->getPeak(endIdx2);
   refEndPts[spec2ID] = shift;
   }

   shift.first.set(spec1F.getMass() - spec2F.getMass(),
   (*spec2)[startIdx2][1],
   spec2->getTolerance(startIdx2)
   + spec1->getTolerance(startIdx1));

   float shiftBack = spec1B.getMass() - spec2B.getMass();
   shiftBack = shiftBack + (*spec2)[endIdx2][0] - (*spec1)[endIdx1][0];
   shift.second.set(shiftBack,
   (*spec2)[endIdx2][1] + (*spec1)[endIdx1][1],
   spec2->getTolerance(endIdx2)
   + spec1->getTolerance(endIdx1));

   //DEBUG_MSG(spec1ID << " - " << spec2ID << ": " << shift.first.toString() << " | " << shift.second.toString());

   if (assembly.count(spec1ID) == 0)
   {
   subAssembly.clear();
   subAssembly[spec2ID] = shift;
   assembly[spec1ID] = subAssembly;
   }
   else
   {
   assembly[spec1ID][spec2ID] = shift;
   }

   shift.first.setMass(0.0 - shift.first.getMass());
   shift.second.setMass(0.0 - shift.second.getMass());

   if (assembly.count(spec2ID) == 0)
   {
   subAssembly.clear();
   subAssembly[spec1ID] = shift;
   assembly[spec2ID] = subAssembly;
   }
   else
   {
   assembly[spec2ID][spec1ID] = shift;
   }
   }

   subAssembly.clear();

   const string rootSpecID = m_assembledSpecs[0].getUniqueID();
   set<string> visited;
   visited.insert(rootSpecID);
   list<pair<string, pair<MZRange, MZRange> > > queue;
   pair<string, pair<MZRange, MZRange> > next;
   next.first = rootSpecID;
   next.second.first.setMass(0);
   next.second.second.setMass(0);
   queue.push_back(next);

   while (queue.size() > 0)
   {
   next = queue.front();
   queue.pop_front();
   subAssembly[next.first] = next.second;

   float shiftF = next.second.first.getMass();
   float shiftB = next.second.second.getMass();

   map<string, pair<MZRange, MZRange> >& neijList = assembly[next.first];

   for (map<string, pair<MZRange, MZRange> >::iterator nIt =
   neijList.begin(); nIt != neijList.end(); nIt++)
   {
   if (visited.count(nIt->first) > 0)
   {
   continue;
   }

   visited.insert(nIt->first);
   next.first = nIt->first;
   next.second.first = nIt->second.first;
   next.second.second = nIt->second.second;
   next.second.first += shiftF;
   next.second.second += shiftB;
   queue.push_back(next);
   }
   }

   //DEBUG_VAR(rootSpecID);

   float minShiftFront = 0;
   float minShiftBack = 0;
   for (map<string, pair<MZRange, MZRange> >::const_iterator aIt =
   subAssembly.begin(); aIt != subAssembly.end(); aIt++)
   {
   if (aIt->second.first.getMass() < minShiftFront)
   {
   minShiftFront = aIt->second.first.getMass();
   }
   if (aIt->second.second.getMass() < minShiftBack)
   {
   minShiftBack = aIt->second.second.getMass();
   }

   //DEBUG_MSG(aIt->first << ": " << aIt->second.first.toString() << " | " << aIt->second.second.toString());
   }

   list<MZRange> newPeaksStart;
   list<MZRange> newPeaksEnd;

   outputStartPtSpec.fileName = "startPtSpec";
   outputStartPtSpec.scan = 1;
   outputStartPtSpec.resize(0);

   outputEndPtSpec.fileName = "endPtSpec";
   outputEndPtSpec.scan = 1;
   outputEndPtSpec.resize(0);

   outputStartPtAligns.resize(subAssembly.size());
   outputEndPtAligns.resize(subAssembly.size());
   unsigned int idxUse = 0;

   for (map<string, pair<MZRange, MZRange> >::const_iterator aIt =
   subAssembly.begin(); aIt != subAssembly.end(); aIt++)
   {
   MZRange newPeakFront(aIt->second.first);
   newPeakFront.setMass(newPeakFront.getMass() - minShiftFront);
   newPeakFront.setIntensity(refEndPts[aIt->first].first.getIntensity());
   newPeakFront.setTolerance(refEndPts[aIt->first].first.getTolerance());

   MZRange newPeakBack(aIt->second.second);
   newPeakBack.setMass(newPeakBack.getMass() - minShiftBack);
   newPeakBack.setIntensity(refEndPts[aIt->first].second.getIntensity());
   newPeakBack.setTolerance(refEndPts[aIt->first].second.getTolerance());

   newPeaksStart.push_back(newPeakFront);
   newPeaksEnd.push_back(newPeakBack);

   Spectrum* spec = m_assembledSpecs.getIndex(aIt->first);

   int endIdx =
   spec->findPeaks(spec->parentMass - specnets::AAJumps::massMH);

   int startIdx = spec->findPeaks(0);

   outputStartPtAligns[idxUse].setSpec1ID(aIt->first);
   outputStartPtAligns[idxUse].setSpec2ID(outputStartPtSpec.getUniqueID());
   outputStartPtAligns[idxUse].clearMP();
   MZRange startPeak((*spec)[startIdx][0],
   (*spec)[startIdx][1],
   spec->getTolerance(startIdx));
   outputStartPtAligns[idxUse].addMatchedPeaksB(startPeak, newPeakFront);

   outputEndPtAligns[idxUse].setSpec1ID(aIt->first);
   outputEndPtAligns[idxUse].setSpec2ID(outputEndPtSpec.getUniqueID());
   MZRange endPeak((*spec)[endIdx][0],
   (*spec)[endIdx][1],
   spec->getTolerance(endIdx));
   outputEndPtAligns[idxUse].addMatchedPeaksB(endPeak, newPeakBack);
   idxUse++;

   leftAlignedShifts[aIt->first] = newPeakFront;
   }

   //outputStartPtSpec.insertPeaks(newPeaksStart);
   outputStartPtSpec.parentCharge = 1;

   for (list<MZRange>::iterator pIt = newPeaksStart.begin(); pIt
   != newPeaksStart.end(); pIt++)
   {
   int pkIdx = outputStartPtSpec.findPeaks(pIt->getMass());
   if (pkIdx >= 0)
   {
   outputStartPtSpec[pkIdx][1] += pIt->getIntensity();
   }
   else
   {
   outputStartPtSpec.insertPeak(&(*pIt));
   }
   }

   outputStartPtSpec.parentMass = outputStartPtSpec[outputStartPtSpec.size()
   - 1][0];

   //outputEndPtSpec.insertPeaks(newPeaksEnd);
   outputEndPtSpec.parentCharge = 1;

   for (list<MZRange>::iterator pIt = newPeaksEnd.begin(); pIt
   != newPeaksEnd.end(); pIt++)
   {
   int pkIdx = outputEndPtSpec.findPeaks(pIt->getMass());
   if (pkIdx >= 0)
   {
   outputEndPtSpec[pkIdx][1] += pIt->getIntensity();
   }
   else
   {
   outputEndPtSpec.insertPeak(&(*pIt));
   }
   }

   outputEndPtSpec.parentMass = outputEndPtSpec[outputEndPtSpec.size() - 1][0];

   }
   */

  void AbruijnGraph::removeSymmetricYNodes()
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

  void AbruijnGraph::removeSourceGreenNodes()
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

  void AbruijnGraph::glueNodes(AbruijnNode* node1, AbruijnNode** node2)
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

      AbruijnEdge* abEdge = AbruijnEdge::castEdgePtr(this->addEdge(fromNode,
                                                                   node1,
                                                                   *edge));

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding edge " << abEdge->getFromNodeIndex() << " -> " << abEdge->getToNodeIndex() << " \'" << abEdge->getFLabel() << "\' index=" << abEdge->getIndex());
      }
    }

    const IndexVector& outEdges = this->getOutEdges((*node2)->getIndex());
    for (unsigned long i = 0; i < outEdges.size(); i++)
    {
      if (outEdges[i] < 0)
        continue;

      Edge* edge = this->getEdge(outEdges[i]);
      Node* toNode = this->getNode(edge->getToNodeIndex());
      AbruijnEdge* abEdge = AbruijnEdge::castEdgePtr(this->addEdge(node1,
                                                                   toNode,
                                                                   *edge));

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding edge " << abEdge->getFromNodeIndex() << " -> " << abEdge->getToNodeIndex() << " \'" << abEdge->getFLabel() << "\' index=" << abEdge->getIndex());
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

  AbruijnNode* AbruijnGraph::mergeNodes(const list<AbruijnEdge*>& edgeList,
                                        AbruijnNode* source,
                                        const string& fLabel,
                                        const string& rLabel,
                                        const string& pathLabel,
                                        set<AbruijnEdge*>& edgesToRemove,
                                        list<AbruijnEdge*>& edgesToRemoveNext,
                                        const bool& lookReverse)
  {
    if (edgeList.size() <= 1)
    {
      return 0;
    }

    AbruijnEdge mergedEdge(edgeList, lookReverse);
    if (lookReverse)
    {
      mergedEdge.setRLabel(pathLabel);
      mergedEdge.setFLabel(fLabel);
    }
    else
    {
      mergedEdge.setFLabel(pathLabel);
      mergedEdge.setRLabel(rLabel);
    }

    set<AbruijnNode*> destNodes;
    for (list<AbruijnEdge*>::const_iterator eIt = edgeList.begin(); eIt
        != edgeList.end(); eIt++)
    {
      AbruijnNode* dest = (lookReverse) ? getNode((*eIt)->getFromNodeIndex())
          : getNode((*eIt)->getToNodeIndex());
      destNodes.insert(dest);
      const string& jumpLabel = (lookReverse) ? (*eIt)->getRLabel()
          : (*eIt)->getFLabel();
      if (jumpLabel == pathLabel)
      {
        edgesToRemove.insert(*eIt);
      }
    }

    AbruijnNode* firstNode = *(destNodes.begin());
    if (destNodes.size() == 1 && ((lookReverse && firstNode->getRPath()
        == rLabel) || (!lookReverse && firstNode->getFPath() == fLabel)))
    {
      Edge* addedEdge = (lookReverse) ? this->addEdge(firstNode,
                                                      source,
                                                      mergedEdge)
          : this->addEdge(source, firstNode, mergedEdge);

      edgesToRemove.insert(edgeList.begin(), edgeList.end());

      if (DEBUG_EXPANSION)
      {
        DEBUG_MSG("Adding " << addedEdge->toString());
      }
      return 0;
    }

    if (DEBUG_EXPANSION)
    {
      DEBUG_MSG("Expanding " << pathLabel << " from " << source->getGraphvizLabel());
      for (list<AbruijnEdge*>::const_iterator eIt = edgeList.begin(); eIt
          != edgeList.end(); eIt++)
      {
        string eIDs = "";
        for (list<string>::const_iterator idIt = (*eIt)->getSpecIDs().begin(); idIt
            != (*eIt)->getSpecIDs().end(); idIt++)
        {
          eIDs += *idIt + " : ";
        }
        if (lookReverse)
        {
          DEBUG_MSG("\'" << (*eIt)->getRLabel() << "\' -> " << (*eIt)->getFromNodeIndex() << " (" << eIDs << ") -> " << getNode((*eIt)->getFromNodeIndex())->getGraphvizLabel());
        }
        else
        {
          DEBUG_MSG("\'" << (*eIt)->getFLabel() << "\' -> " << (*eIt)->getToNodeIndex() << " (" << eIDs << ") -> " << getNode((*eIt)->getToNodeIndex())->getGraphvizLabel());
        }
      }
    }

    set<AbruijnEdge*> sourceEdges(edgeList.begin(), edgeList.end());
    set<string> edgeSpecs(mergedEdge.getSpecIDs().begin(),
                          mergedEdge.getSpecIDs().end());

    unsigned int newSize = 0;
    bool isGreen = (lookReverse) ? rLabel.length() > 0 : fLabel.length() > 0;

    for (set<AbruijnNode*>::const_iterator nIt = destNodes.begin(); nIt
        != destNodes.end(); nIt++)
    {
      newSize += (*nIt)->size();
    }

    AbruijnNode* newNode = AbruijnNode::castNodePtr(this->addNode());
    newNode->setFPath(fLabel);
    newNode->setRPath(rLabel);

    AbruijnEdge* addedMergedEdge = (lookReverse)
        ? AbruijnEdge::castEdgePtr(addEdge(newNode, source, mergedEdge))
        : AbruijnEdge::castEdgePtr(addEdge(source, newNode, mergedEdge));

    if (DEBUG_EXPANSION)
    {
      DEBUG_MSG("Adding " << addedMergedEdge->toString());
    }

    newNode->resize(newSize);
    newNode->setGreen(isGreen);
    unsigned int idxUse = 0;
    //unsigned int numIDBins = edgeSpecIDs.size() * 2;

    //DEBUG_VAR(source->getIndex());

    //unsigned int debugIdx = 733;

    set<string> usedSpecIDs;
    //set<AbruijnNode*> newInNodes;
    //set<AbruijnNode*> newOutNodes;

    for (set<AbruijnNode*>::iterator nIt = destNodes.begin(); nIt
        != destNodes.end(); nIt++)
    {
      for (unsigned int i = 0; i < (*nIt)->size(); i++)
      {
        const string& specID = (*nIt)->operator [](i).getSpecID();
        if (usedSpecIDs.count(specID) == 0)
          (*newNode)[idxUse++].operator =((*nIt)->operator [](i));
        usedSpecIDs.insert(specID);
      }

      if (!lookReverse)
      {
        const IndexVector& outEdges = this->getOutEdges(*nIt);
        for (unsigned long i = 0; i < outEdges.size(); i++)
        {
          if (outEdges[i] < 0)
            continue;

          //DEBUG_VAR(outEdges[i]);
          AbruijnEdge* abEdge = this->getEdge(outEdges[i]);
          AbruijnNode* dest = this->getNode(abEdge->getToNodeIndex());

          if (dest->getIndex() == newNode->getIndex())
          {
            continue;
          }

          const unsigned int abJumpLen = abEdge->getJumpLabel().length();
          AbruijnEdge* addedEdge = 0;

          // add outgoing label free edges to red nodes only
          if (abJumpLen == 0 && (!isGreen))
          {
            addedEdge
                = AbruijnEdge::castEdgePtr(addEdge(newNode, dest, *abEdge));
            if (DEBUG_EXPANSION)
            {
              string eIDs = "";
              for (list<string>::const_iterator idIt =
                  abEdge->getSpecIDs().begin(); idIt
                  != abEdge->getSpecIDs().end(); idIt++)
              {
                eIDs += *idIt + " : ";
              }
              DEBUG_MSG("Copying edge " << abEdge->toString() << " to " << addedEdge->toString() << " (" << eIDs << ")");
            }
            continue;
          }
          else if (abJumpLen == 0)
          {
            continue;
          }

          const bool isPrefixFor = isPrefix(abEdge->getFLabel(), fLabel);
          const bool isExtenFor = isPrefix(fLabel, abEdge->getFLabel());

          if ((!isPrefixFor) && (!isExtenFor))
          {
            continue;
          }

          const string revExt = abEdge->getRLabel().substr(abJumpLen);
          const bool isPrefixRev = isPrefix(revExt, rLabel);
          //const bool isExtenRev = isPrefix(rLabel, revExt);

          if (!isPrefixRev)
          {
            continue;
          }

          addedEdge = AbruijnEdge::castEdgePtr(addEdge(newNode, dest, *abEdge));
          addedEdge->setExpandOnly(false);
          //newOutNodes.insert(dest);

          const string revJump =
              specnets::AAJumps::reversePeptide(abEdge->getJumpLabel());

          if (isGreen)
            addedEdge->setRLabel(revJump + rLabel);

          if (DEBUG_EXPANSION)
          {
            string eIDs = "";
            for (list<string>::const_iterator idIt =
                abEdge->getSpecIDs().begin(); idIt
                != abEdge->getSpecIDs().end(); idIt++)
            {
              eIDs += *idIt + " : ";
            }
            DEBUG_MSG("Copying edge " << abEdge->toString() << " to " << addedEdge->toString() << " (" << eIDs << ")");
          }

          if (isPrefixFor && (!isExtenFor))
          {
            addedEdge->setExpandOnly(true);
            edgesToRemoveNext.push_back(addedEdge);
            if (DEBUG_EXPANSION)
            {
              DEBUG_MSG("Registered for deletion: " << addedEdge->toString());
            }
          }
        }
      }

      const IndexVector& inEdges = this->getInEdges(*nIt);

      for (unsigned long i = 0; i < inEdges.size(); i++)
      {
        if (inEdges[i] < 0)
          continue;

        //DEBUG_VAR(inEdges[i]);
        AbruijnEdge* abEdge = this->getEdge(inEdges[i]);
        AbruijnNode* dest = this->getNode(abEdge->getFromNodeIndex());

        if (dest->getIndex() == newNode->getIndex())
        {
          continue;
        }

        if ((!lookReverse) && dest->getIndex() == source->getIndex())
        {
          continue;
        }

        const unsigned int abJumpLen = abEdge->getJumpLabel().length();
        AbruijnEdge* addedEdge = 0;

        // add outgoing label free edges to red nodes only
        if (abJumpLen == 0 && (!isGreen) && lookReverse)
        {
          addedEdge = AbruijnEdge::castEdgePtr(addEdge(dest, newNode, *abEdge));
          if (DEBUG_EXPANSION)
          {
            string eIDs = "";
            for (list<string>::const_iterator idIt =
                abEdge->getSpecIDs().begin(); idIt
                != abEdge->getSpecIDs().end(); idIt++)
            {
              eIDs += *idIt + " : ";
            }
            DEBUG_MSG("Copying edge " << abEdge->toString() << " to " << addedEdge->toString() << " (" << eIDs << ")");
          }
          continue;
        }
        else if (abJumpLen == 0)
        {
          continue;
        }

        const bool isPrefixRev = isPrefix(abEdge->getRLabel(), rLabel);
        const bool isExtenRev = isPrefix(rLabel, abEdge->getRLabel());

        if ((!isPrefixRev) && (!isExtenRev))
        {
          continue;
        }

        const string forExt = abEdge->getFLabel().substr(abJumpLen);
        const bool isPrefixFor = isPrefix(forExt, fLabel);
        //const bool isExtenFor = isPrefix(fLabel, forExt);

        if (!isPrefixFor)
        {
          continue;
        }

        if (!lookReverse)
        {
          bool foundSameSpec = false;
          for (list<string>::const_iterator idIt = abEdge->getSpecIDs().begin(); idIt
              != abEdge->getSpecIDs().end(); idIt++)
          {
            if (edgeSpecs.count(*idIt) > 0)
            {
              foundSameSpec = true;
              break;
            }
          }
          if (foundSameSpec)
            continue;
        }

        addedEdge = AbruijnEdge::castEdgePtr(addEdge(dest, newNode, *abEdge));
        //newInNodes.insert(dest);

        if ((!lookReverse) || isGreen)
          addedEdge->setFLabel(abEdge->getJumpLabel() + fLabel);

        if ((!lookReverse) && isPrefixRev)
        {
          addedEdge->setRLabel(rLabel);
        }

        if (DEBUG_EXPANSION)
        {
          string eIDs = "";
          for (list<string>::const_iterator idIt = abEdge->getSpecIDs().begin(); idIt
              != abEdge->getSpecIDs().end(); idIt++)
          {
            eIDs += *idIt + " : ";
          }
          DEBUG_MSG("Copying edge " << abEdge->toString() << " to " << addedEdge->toString() << " (" << eIDs << ")");
        }

        if (isPrefixRev && (!isExtenRev))
        {
          addedEdge->setExpandOnly(true);
          if (lookReverse)
          {
            edgesToRemoveNext.push_back(addedEdge);
            if (DEBUG_EXPANSION)
            {
              DEBUG_MSG("Registered for deletion: " << addedEdge->toString());
            }
          }
          else
          {
            if (DEBUG_EXPANSION)
            {
              DEBUG_MSG("Registered for deletion after reverse expand: " << addedEdge->toString());
            }
          }
        }
      }
    }
    (*newNode).resize(idxUse);

    return newNode;
  }
}
