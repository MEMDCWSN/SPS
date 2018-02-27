/*
 * AbruijnGraph.h
 *
 *  Created on: May 1, 2012
 *      Author: aguthals
 */

#ifndef ABRUIJNGRAPH_H_
#define ABRUIJNGRAPH_H_

#include "BaseGraph.h"
#include "aminoacid.h"
#include "AbruijnNode.h"
#include "AbruijnEdge.h"
#include "SpectrumAlignment.h"
#include "SpectrumAlignmentSet.h"
#include "SpecSet.h"
#include "spectrum.h"
#include "JumpEdge.h"

#include <vector>
#include <string>
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
  class AbruijnEdge;
  class SpectrumAlignment;
  class SpectrumAlignmentSet;

  typedef list<pair<AbruijnNode*, AbruijnNode*> > AbruijnAlignment;

  class AbruijnGraph : public specnets::BaseGraph,
                       public specnets::Node
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

    // Maximum number of AA jumps allowed per abruijn edge
    static const unsigned int MAX_NUM_JUMPS;

    // Maximum number of PTMs allowed per abruijn edge
    static const unsigned int MAX_NUM_MODS_PER_JUMP;

    // Maximum jump mass
    static const double MAX_JUMP_MASS;

    // Universal identifier given the "start-point" spectrum
    static const string STARTPT_SPEC_ID;

    // Universal identifier given the "end-point" spectrum
    static const string ENDPT_SPEC_ID;

    // Maximum number of consecutive AA that are expanded from every node
    static unsigned int PATH_EXPAND_LIMIT;

    // Whether or not to print verbose debug messages
    static bool DEBUG_EXPANSION;

    // If true, then star spectra should be reversed if specified by abinfo
    static bool REVERSE_STARS;

    // If true, enforce that b-endpoints be used as endpoints. Otherwise, chose the configuration
    // that yields the highest scoring heaviest path
    static bool ENFORCE_B_ENDPTS;

    // Universal penalty for multi-jump edge scores
    static double GAP_PENALTY;

    // Upcasts a BaseGraph to an AbruijnGraph
    static AbruijnGraph* castGraphPtr(specnets::BaseGraph* graphPtr);

    // Upcasts a Node to an AbruijnNode
    static AbruijnGraph* castNodePtr(specnets::Node* nodePtr);

    // Creates a string identifier from a set of AbruijnNodes
    static string nodeSetToString(const set<AbruijnNode*>& inputNodeSet);

    // Default constructor
    AbruijnGraph();

    // Copy constructor
    AbruijnGraph(const AbruijnGraph& other);

    /**
     * Backwards-compatible copy constructor for old abinfo format. This will construct
     *   the abruijn graph, expand parallel paths, and find the consensus.
     * @param inAbinfo
     * @param assembledSpecs set of star spectra referenced by abinfo
     * @param peakTol peak toelrance of star spectra
     * @param outAbinfo if non-NULL, this will be filled with the resulting contig structure
     */
    AbruijnGraph(const pair<pair<vector<int> , vector<int> > , vector<pair<
                     vector<int> , vector<double> > > >& inAbinfo,
                 const SpecSet& assembledSpecs,
                 const float& peakTol,
                 pair<pair<vector<int> , vector<int> > , vector<pair<
                     vector<int> , vector<double> > > >* outAbinfo = 0);

    virtual ~AbruijnGraph(void)
    {
    }

    /**
     * Clears out all class data structures
     */
    void initialize();

    /**
     * Compresses all node and edge indices
     * @param outputNewNodeIdxs if specified, this will contain a mapping
     *   of old node indices to new indices
     * @param outputNewEdgeIdxs if specified, this will contain a mapping
     *   of old edge indices to new indices
     */
    virtual void compress(vector<long>* outputNewNodeIdxs = 0,
                          vector<long>* outputNewEdgeIdxs = 0);

    virtual void clear(void)
    {
      initialize();
    }

    /**
     * Has the same function as the constructor
     */
    void
    reSequencePPaths(const pair<pair<vector<int> , vector<int> > , vector<pair<
                         vector<int> , vector<double> > > >& abinfo,
                     const SpecSet& assembledSpecs,
                     const float& peakTol,
                     pair<pair<vector<int> , vector<int> > , vector<pair<
                         vector<int> , vector<double> > > >* outAbinfo = 0);

    /**
     * Clones a Node as an AbruijnNode and returns a pointer to the new
     *   AbruijnNode (as a Node pointer)
     * @param copyNode
     * @return new Node pointer
     */
    virtual specnets::Node* cloneNode(specnets::Node& copyNode) const;

    /**
     * Clones an Edge as an AbruijnEdge and returns a pointer to the new
     *   AbruijnEdge (as an Edge pointer)
     * @param copyEdge
     * @return new Edge pointer
     */
    virtual specnets::Edge* cloneEdge(specnets::Edge& copyEdge) const;

    /**
     * Creates a new AbruijnNode and returns its pointer (as a Node pointer)
     * @return new Node pointer
     */
    virtual specnets::Node* createNode(void) const;

    /**
     * Creates a new AbruijnEdge and returns its pointer (as an Edeg pointer)
     * @return new Edge pointer
     */
    virtual specnets::Edge* createEdge(void) const;

    /**
     * Returns the AbruijnNode pointer for a given node index
     * @param nodeIndex
     * @return AbruijnNode pointer
     */
    inline AbruijnNode* getNode(unsigned long nodeIndex) const
    {
      return AbruijnNode::castNodePtr(BaseGraph::getNode(nodeIndex));
    }

    /**
     * Returns the AbruijnEdge pointer for a given edge index
     * @param nodeIndex
     * @return AbruijnEdge pointer
     */
    inline AbruijnEdge* getEdge(unsigned long edgeIndex) const
    {
      return AbruijnEdge::castEdgePtr(BaseGraph::getEdge(edgeIndex));
    }

    /**
     * Copies another AbruijnGraph and returns this graph's reference
     * @param other
     * @return reference to this graph
     */
    AbruijnGraph &operator=(const AbruijnGraph &other);

    /**
     * Copies another BaseGraph, which is upcasted to an AbruijnGraph
     * @param otherGraph
     */
    virtual void copy(specnets::BaseGraph& otherGraph);

    /**
     * Adds the current binary format version and subversion to a map
     */
    virtual void
    addBinaryVersionInfo(map<string, unsigned short>& versions) const;

    /**
     * Outputs the consensus path as a Soectru
     * @param outputSpec
     */
    void outputConsensusSpectrum(Spectrum& outputSpec) const;

    //void appendGraph(const AbruijnGraph& otherGraph);

    //void removeEmptyNodes(void);
    /*
     void loadAssembly(const specnets::SpecSet& prmSpectra,
     const SpectrumAlignmentSet& prmAligns,
     bool addEndpointGlues = true);

     //string getConsensusSeq(const bool reportMods = false) const;
     * */

  protected:

    //SpecSet m_assembledSpecs;
    //SpectrumAlignmentSet m_alignments;

    // Heaviest path through the A-Bruijn graph
    Path m_consensusPath;

    // Nodes connected by the heaviest path
    vector<AbruijnNode*> m_consensusPathVerts;

    // Score of the heaviest path
    double m_consensusPathScore;

    // Whether or not this contig has been reversed
    bool m_reversed;

    // Set of amino acid masses used to impose jump edges in star spectra
    specnets::AAJumps& m_globalJumps;

    // Lookups AbruijnNodes by their assembled star peaks
    /**
     * Spectrum ID -> Peak Mass -> (first)  AbruijnNode
     *                             (second) [0] true if this peak's score was added to the incoming AbruijnEdge(s)
     *                                      [1] true if this peak's score was added to the reverse score of outgoing AbruijnEdge(s)
     */
    tr1::unordered_map<string, map<MZRange,
        pair<AbruijnNode*, TwoValues<bool> > > > m_nodeLookup;

    // REMAINING USED FOR INTERIOR PROCESSING

    // The set of nodes that require expansion
    set<AbruijnNode*> m_nodesToExpand;

    // The shift of each star spectrum wrt the left-most star
    map<string, double> m_startPtAssembly;

    void
    m_reSequencePPaths(const pair<pair<vector<int> , vector<int> > , vector<
                           pair<vector<int> , vector<double> > > >& abinfo,
                       const SpecSet& assembledSpecs,
                       const float& peakTol,
                       const bool reverseStars,
                       const bool useBEndpts,
                       pair<pair<vector<int> , vector<int> > , vector<pair<
                           vector<int> , vector<double> > > >* outAbinfo = 0);

    void m_internalCopy(const AbruijnGraph &other);

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

    AbruijnEdge* addLabelFreeAbEdge(AbruijnNode* from,
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
                         list<AbruijnEdge*>* addedEdges = 0);

    void expandForwardReversePaths(const float& peakTol);

    void expandPaths(AbruijnNode* source, const bool& goReverse);

    void mergePaths(AbruijnNode* start, const bool& goReverse);

    void processAllPaths(const bool& goReverse, const bool& expandYes);

    /**
     * Combines edges with the same label between the same two nodes (adds together weight)
     */
    void mergeParallelEdges(AbruijnNode* from, AbruijnNode* to, bool mergeMods =
        false);

    /*
     * Calls mergeParallelEdges on all connected node pairs
     */
    void mergeAllParallelEdges(bool mergeMods = false);

    /*
     * Prunes any edge between two nodes that is not the maximal scoring edge over all others
     *   with the same jump mass
     */
    void pruneParallelEdgesByMass(AbruijnNode* from,
                                  AbruijnNode* to,
                                  float peakTol);

    /*
     * Calls pruneParallelEdgesByMass on all connected node pairs
     */
    void pruneAllParallelEdgesByMass(float peakTol);

    /**
     * Combines edges with no label but the same mass between the same two nodes (adds together weight)
     */
    void mergeParallelLabelFreeEdgesByMass(AbruijnNode* from,
                                           AbruijnNode* to,
                                           float peakTol);

    /**
     * Calls mergeParallelLabelFreeEdgesByMass on all connected node pairs
     */
    void mergeAllParallelLabelFreeEdgesByMass(float peakTol);

    bool computeHeaviestPathDAG();

    //void assembleConsensus(const bool& addEndpointGlues);

    void removeExpandedPath(AbruijnEdge** pathPartEdge);

    void injectExpandedPath(AbruijnEdge* templateEdge,
                            const unsigned int& step,
                            const bool& goReverse);
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

    AbruijnNode* mergeNodes(const list<AbruijnEdge*>& edgeList,
                            AbruijnNode* source,
                            const string& fLabel,
                            const string& rLabel,
                            const string& pathLabel,
                            set<AbruijnEdge*>& edgesToRemove,
                            list<AbruijnEdge*>& edgesToRemoveNext,
                            const bool& lookReverse);

  };
}

#endif /* ABRUIJNGRAPH_H_ */
