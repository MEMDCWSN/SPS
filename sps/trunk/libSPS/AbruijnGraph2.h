/*
 * AbruijnGraph22.h
 *
 *  Created on: May 1, 2012
 *      Author: aguthals
 */

#ifndef ABRUIJNGRAPH2_H_
#define ABRUIJNGRAPH2_H_

#include "BaseGraph.h"
#include "aminoacid.h"
#include "AbruijnNode.h"
#include "SpectrumAlignment.h"
#include "SpectrumAlignmentSet.h"
#include "SpecSet.h"
#include "spectrum.h"
#include "JumpEdge.h"

#include <vector>
#include <string>
#include <set>
#include <stdlib.h>

using namespace std;
using namespace specnets;

namespace specnets
{
  class Node;
  class Edge;
  class BaseGraph;
}

namespace abruijn
{
  class AbruijnNode;
  class SpectrumAlignment;
  class SpectrumAlignmentSet;

  typedef list<pair<AbruijnNode*, AbruijnNode*> > AbruijnAlignment;

  class AbruijnGraph2 : public specnets::BaseGraph,
                        public specnets::Node
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    static const unsigned int MAX_NUM_JUMPS;

    static const unsigned int MAX_NUM_MODS_PER_JUMP;

    static const double MAX_JUMP_MASS;

    static const string STARTPT_SPEC_ID;

    static const string ENDPT_SPEC_ID;

    static unsigned int NUM_PATHS_PER_NODE;

    static bool DEBUG_EXPANSION;

    static bool REVERSE_STARS;

    static bool ENFORCE_B_ENDPTS;

    static double GAP_PENALTY;

    static double PEAK_TOLERANCE;

  private:

    class ParallelEdge
    {
    public:

      ParallelEdge();

      ParallelEdge(const JumpEdge& newEdge);

      ParallelEdge(const ParallelEdge& other);

      ParallelEdge &operator=(const JumpEdge &other);

      ParallelEdge &operator=(const ParallelEdge &other);

      inline bool operator==(const ParallelEdge &other) const
      {
        return m_jumpCode == other.m_jumpCode;
      }

      inline bool addEdge(const JumpEdge& newEdge)
      {
        if (!std::binary_search(m_edgeIds.begin(),
                                m_edgeIds.end(),
                                newEdge.getIndex()))
        {
          m_edgeIds.push_back(newEdge.getIndex());
          std::sort(m_edgeIds.begin(), m_edgeIds.end());
          return true;
        }
        return false;
      }

      int m_jumpCode;
      vector<long> m_edgeIds;
    };

    class ParallelPath
    {

    public:

      ParallelPath();

      ParallelPath(const ParallelPath& other);

      ParallelPath(const AbruijnGraph2& abG, const AbruijnNode& start);

      inline void clear()
      {
        m_path.resize(0);
        m_nodes.resize(0);
        m_pathNodes.clear();
        m_fWeight = 0;
        m_rWeight = 0;
        m_outgoingEdgeIDs.clear();

        m_outgoingEdges.resize(0);
      }

      void initialize(const AbruijnGraph2& abG, const AbruijnNode& start);

      bool operator==(const ParallelPath &other) const;

      ParallelPath &operator=(const ParallelPath &other);

      void advance(const AbruijnGraph2& abG, const JumpEdge& advanceEdge);

      bool merge(const AbruijnGraph2& abG, ParallelPath& other);

      bool addNode(const AbruijnGraph2& abG,
                   const unsigned int& pathIdx,
                   const long& node);

      void addEdge(const unsigned int& pathIdx, const JumpEdge& newEdge);

      inline unsigned int getLength() const
      {
        return m_path.size();
      }

      inline double getWeight() const
      {
        return max(m_fWeight, m_rWeight);
      }

      string toString() const;

      vector<ParallelEdge> m_path;
      vector<AbruijnNode> m_nodes;
      vector<JumpEdge> m_outgoingEdges;
      std::set<long> m_outgoingEdgeIDs;

      std::set<long> m_pathNodes;
      double m_fWeight;
      double m_rWeight;

      static inline bool WeightComparator(const ParallelPath* lhs,
                                          const ParallelPath* rhs)
      {
        return lhs->getWeight() > rhs->getWeight();
      }
    };

  public:

    static AbruijnGraph2* castGraphPtr(specnets::BaseGraph* graphPtr);

    static AbruijnGraph2* castNodePtr(specnets::Node* nodePtr);

    AbruijnGraph2();

    AbruijnGraph2(const AbruijnGraph2& other);

    AbruijnGraph2(const pair<pair<vector<int> , vector<int> > , vector<pair<
                      vector<int> , vector<double> > > >& inAbinfo,
                  const SpecSet& assembledSpecs,
                  pair<pair<vector<int> , vector<int> > , vector<pair<vector<
                      int> , vector<double> > > >* outAbinfo = 0);

    virtual ~AbruijnGraph2(void)
    {
    }

    void initialize();

    virtual void compress(vector<long>* outputNewNodeIdxs = 0,
                          vector<long>* outputNewEdgeIdxs = 0);

    virtual void clear(void)
    {
      initialize();
    }

    void
    reSequencePPaths(const pair<pair<vector<int> , vector<int> > , vector<pair<
                         vector<int> , vector<double> > > >& abinfo,
                     const SpecSet& assembledSpecs,
                     pair<pair<vector<int> , vector<int> > , vector<pair<
                         vector<int> , vector<double> > > >* outAbinfo = 0);

    virtual specnets::Node* cloneNode(specnets::Node& copyNode) const;

    virtual specnets::Edge* cloneEdge(specnets::Edge& copyEdge) const;

    virtual specnets::Node* createNode(void) const;

    virtual specnets::Edge* createEdge(void) const;

    inline AbruijnNode* getNode(unsigned long nodeIndex) const
    {
      return AbruijnNode::castNodePtr(BaseGraph::getNode(nodeIndex));
    }

    inline JumpEdge* getEdge(unsigned long edgeIndex) const
    {
      return JumpEdge::castEdgePtr(BaseGraph::getEdge(edgeIndex));
    }

    AbruijnGraph2 &operator=(const AbruijnGraph2 &other);

    virtual void copy(specnets::BaseGraph& otherGraph);

    virtual void
    addBinaryVersionInfo(map<string, unsigned short>& versions) const;

    void outputConsensusSpectrum(Spectrum& outputSpec) const;

  protected:

    //SpecSet m_assembledSpecs;
    //SpectrumAlignmentSet m_alignments;

    ParallelPath m_consensusPath;
    tr1::unordered_map<string, map<MZRange,
        pair<AbruijnNode*, TwoValues<bool> > > > m_nodeLookup;

    set<AbruijnNode*> m_nodesToExpand;
    map<string, double> m_startPtAssembly;

    specnets::AAJumps& m_globalJumps;

    bool m_reversed;

    void
    m_reSequencePPaths(const pair<pair<vector<int> , vector<int> > , vector<
                           pair<vector<int> , vector<double> > > >& abinfo,
                       const SpecSet& assembledSpecs,
                       const bool reverseStars,
                       const bool useBEndpts,
                       pair<pair<vector<int> , vector<int> > , vector<pair<
                           vector<int> , vector<double> > > >* outAbinfo = 0);

    void m_internalCopy(const AbruijnGraph2 &other);

    void addSpectrum(const specnets::Spectrum& prmSpec, bool allEndPts = false);

    inline void addSpectra(const specnets::SpecSet& prmSpectra)
    {
      for (unsigned int i = 0; i < prmSpectra.size(); i++)
      {
        addSpectrum(prmSpectra[i]);
      }
    }

    void addGlues(const SpectrumAlignment& prmAlign);

    inline void addGlues(const SpectrumAlignmentSet& prmAligns)
    {
      for (unsigned int i = 0; i < prmAligns.size(); i++)
      {
        addGlues(prmAligns[i]);
      }
    }

    AbruijnNode* addAbruijnNode(AbruijnNode& copyNode);

    JumpEdge* addLabelFreeAbEdge(AbruijnNode* from,
                                 AbruijnNode* to,
                                 const double& mass,
                                 const double& weight,
                                 const double& rWeight);

    void recoverSourceNodeScores();

    void reverseEdgeWeights();

    void addEndpointEdges(const set<string>& assembledSpecIDs,
                          const SpecSet& assmebledSpecs);

    void addAbruijnEdges(const Spectrum& prmSpec,
                         const unsigned int& peakIdxFrom,
                         const unsigned int& peakIdxTo,
                         const double& avgSpecScore,
                         list<JumpEdge*>* addedEdges = 0);

    /**
     * Combines edges with the same label between the same two nodes (adds together weight)
     */
    void mergeParallelEdges(AbruijnNode* from, AbruijnNode* to, bool mergeMods =
        false);

    /*
     * Calls mergeParallelEdges on all connected node pairs
     */
    void mergeAllParallelEdges(bool mergeMods = false);

    /**
     * Combines edges with no label but the same mass between the same two nodes (adds together weight)
     */
    void mergeParallelLabelFreeEdgesByMass(AbruijnNode* from, AbruijnNode* to);

    /**
     * Calls mergeParallelLabelFreeEdgesByMass on all connected node pairs
     */
    void mergeAllParallelLabelFreeEdgesByMass();

    //bool computeHeaviestPathDAG();

    bool computeHeaviestParallelPathDAG();

    /*
     void getEndPointSpectra(map<string, specnets::MZRange>& leftAlignedShifts,
     specnets::Spectrum& outputStartPtSpec,
     specnets::Spectrum& outputEndPtSpec,
     SpectrumAlignmentSet& outputStartPtAligns,
     SpectrumAlignmentSet& outputEndPtAligns);
     */

    void removeSymmetricYNodes();

    void removeSourceGreenNodes();

    void glueNodes(AbruijnNode* node1, AbruijnNode** node2);

  };
}

#endif /* ABRUIJNGRAPH2_H_ */
