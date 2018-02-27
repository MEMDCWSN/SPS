/*
 * ContigNetwork.h
 *
 *  Created on: Jan 5, 2012
 *      Author: aguthals
 */

#ifndef CONTIGNETWORK_H_
#define CONTIGNETWORK_H_

using namespace std;

#include "Contig.h"
#include "prm_alignment.h"
#include "PRMSpecEdge.h"

namespace specnets
{

  /**
   * This class was designed to merge contig PRM spectra into meta-contigs,
   * but it can be used to merge any group of aligned PRM spectra into contigs
   */
  class ContigNetwork
  {
  public:

    // Should not be needed but ExecAssembly needs them
    float m_pkTol;
    float m_pmTol;

    /**
     * Default constructor. Must call initialize() after calling this
     */
    ContigNetwork(void);

    /**
     * Primary constructor, calls initialize with the same parameters.
     * @param _rootSpectra PRM spectra that will be merged into contigs. Each of these is
     *   initialized to a Contig containing one spectrum.
     * @param alignments alignments of _m_rootSpectra, each of which will impose an alignment
     *   edge between contigs.
     */
    ContigNetwork(SpecSet& _rootSpectra,
                  SpectrumPairSet& alignments,
                  float pkTol = 0.5,
                  float pmTol = 1.0);

    /**
     * Deconstructor. De-allocates all class variables stored on the heap
     */
    virtual ~ContigNetwork(void);

    /**
     * Initializes the m_graph to contain a node for evey contig and an edge for every
     *   alignment. Also allocates memory for class data structures.
     * @param _rootSpectra PRM spectra that will be merged into contigs. Each of these is
     *   initialized to a Contig containing one spectrum.
     * @param alignments alignments of _m_rootSpectra, each of which will impose an
     *   alignment edge between contigs.
     */
    virtual void initialize(SpecSet& _rootSpectra,
                            SpectrumPairSet& alignments,
                            float pkTol = 0.5,
                            float pmTol = 1.0);

    /**
     * Iteratively (1) merges the highest scoring edge and (2) updates the scores of
     *   neighboring m_edges until the highest scoring m_edges has score < minScore
     * @param minScore minimum allowable edge/alignment score
     * @param minNumMatchedPeaks minimum allowable number matching peaks between contigs
     * @return the number of m_edges that were merged
     */
    virtual int assembleIteratively(float minScore, int minNumMatchedPeaks);

    /**
     * Outputs abinfo of contigs assembling at least 2 spectra or more.
     * @param outputAB
     */
    virtual void outputAbinfo(abinfo_t& outputAB);

    /**
     * Outputs consensus spectra of contigs w/ matching indices
     * @param outputSpecs
     */
    virtual void outputConsensusSpectra(SpecSet& outputSpecs);

    virtual void outputComponents(vector<set<unsigned int> >& outputComps);

    /**
     * Outputs a reversed flag for eeach meta-contig w/ matching indices
     * @param outputRev
     */
    virtual void outputReversed(vector<bool>& outputRev);

    /**
     * Merges the edge at edgeIdx. This will reverse a node connected to the edge if required
     *   by the alignment (by calling reverseNode()). Then it will merge the two Contigs (by
     *   calling Contig::merge()), update neighboring m_edges (without rescoring them), and
     *   remove one of the two merged Contigs (the other one becoming the merged Contig)
     * @param edgeIdx Index of edge to be merged. This edge will no longer exist if the
     *   function is successful.
     * @return true if an edge at edgeIdx exists and it was merged, false if not
     */
    virtual bool mergeEdge(int edgeIdx, int minNumMatchedPeaks);

    /**
     * Re-assigns the alignment scores of all m_edges connected to the Contig at nodeIdx.
     *   Alignment scores are computed by computeEdgeScore()
     * @param nodeIdx index of Contig to re-score m_edges to.
     * @return true if the Contig at nodeIdx exists and m_edges were re-scored, false if not
     */
    virtual bool rescoreConnectingEdges(int nodeIdx);

    /**
     * Generates a list of all neighbors to a given Contig node along with the indices of all
     *   connecting m_edges. The returned data structure must be de-allocated by the caller.
     * @param metaContigIdx index of Contig node
     * @return pointer to a list of all neighbors to metaContigIdx. Each element contains
     *   first=index of neighboring node and second=index of connecting edge. This must be
     *   de-allocated by the caller via "delete"
     */
    virtual list<pair<int, int> >* getNeighborList(int metaContigIdx) const;

    /**
     * Returns true if a node with index metaContigIdx exists in the m_graph, false if not
     */
    virtual bool containsNode(int metaContigIdx) const;

    /**
     * Returns true if an edge between node1 and node2 exists in the m_graph, false if not
     */
    virtual bool containsEdge(int node1, int node2) const;

    /**
     * Returns true if an edge with index edgeIdx exists in the m_graph, false if not
     */
    virtual bool containsEdge(int edgeIdx) const;

    /**
     * Returns a pointer to the Contig at index in the m_graph.
     * @param index
     * @return Pointer to Contig at that index. If Contig does not exist, An error message
     *   is printed and 0 is returned.
     */
    virtual Contig* getNode(int index) const;

    /**
     * Returns a pointer to the edge at index in the m_graph.
     * @param index
     * @return Pointer to edge at that index. If edge does not exist, An error message
     *   is printed and 0 is returned.
     */
    virtual PRMSpecEdge* getEdge(int index) const;

    /**
     * Returns a pointer to the edge between node1 and node2.
     * @param node1
     * @param node2
     * @return Pointer to the edge. If the edge does not exist, An error message
     *   is printed and 0 is returned.
     */
    virtual PRMSpecEdge* getEdge(int node1, int node2) const;

    /**
     * Adds a Contig to the m_graph as a new node. This fails if a node at newNode.index
     *   already exists.
     * @param newNode Contig that will be copied into the m_graph
     * @return true if node was added, false if not
     */
    virtual bool addNode(const Contig& newNode);

    /**
     * Adds an edge to the m_graph. This fails if an edge at newEdge.index
     *   already exists ir an edge between m_nodes newEdge.spec1 and newEdge.spec2 already exists
     * @param newEdge Edge that will be copied into the m_graph
     * @return true if edge was added, false if not
     */
    virtual bool addEdge(const PRMSpecEdge& newEdge);

    /**
     * Removes the node at nodeIdx from the m_graph
     * @param nodeIdx
     * @return true if node existed and was removed, false if not
     */
    virtual bool removeNode(int nodeIdx);

    /**
     * Removes the edge at edgeIdx from the m_graph
     * @param edgeIdx
     * @return true if the edge existed and was removed, false if not
     */
    virtual bool removeEdge(int edgeIdx);

    /**
     * Removes the edge between node1 and node2 from the m_graph
     * @param node1
     * @param node2
     * @return true if the edge existed and was removed, false if not
     */
    virtual bool removeEdge(int node1, int node2);

    /**
     * Reverses the Contig at index metaContigIdx. Also updates all connecting m_edges
     * @param metaContigIdx
     * @return true if node existed and was reversed, false if not
     */
    virtual bool reverseNode(int metaContigIdx);

    /**
     * Computes the shift between the consensus PRM spectra of m_nodes connected by a given edge
     * @param edge edge between 2 Contig m_nodes
     * @param nodeFrom either edge->spec1 or edge->spec2
     * @param nodeTo either edge->spec2 or edge->spec1
     * @return shift of Contig at nodeTo wrt Contig at nodeFrom
     */
    virtual float
    getConsensusShift(PRMSpecEdge& edge, int nodeFrom, int nodeTo) const;

    /**
     * Computes the score of the alignment between the consensus PRM spectra of Contig m_nodes
     *   connected by an input edge.
     * @param edge pointer to edge between 2 m_nodes
     * @return pair.first=score of overlap with edge.spec1, pair-second=score of overlap with edge.spec2
     */
    virtual pair<float, float> getEdgeScore(PRMSpecEdge& edge);

    /**
     * Returns a pointer to the highest scoring edge in the m_graph
     */
    virtual PRMSpecEdge* getHighestScoringEdge();

    /**
     * Gets a new unique node index
     */
    virtual int getNewNodeIndex();

    /**
     * Gets a new unique edge index
     */
    virtual int getNewEdgeIndex();

  protected:

    virtual void addEdge(int spec1, int spec2, int edgeIdx);

    void eraseNodes();

    void eraseEdges();

    // Keep track of maximum node and edge indices so none overlap
    int m_maxNodeIdx;
    int m_maxEdgeIdx;

    /**
     * Number of meta-contigs
     */
    int m_size;

    /**
     * Alignment object for re-scoring m_edges
     */
    PRMAlignment m_alignmentObj;

    /**
     * SPS contigs
     */
    SpecSet* m_rootSpectra;

    /**
     * Holds edge references for easy lookup
     */
    map<int, PRMSpecEdge*>* m_edges;

    /**
     * Assembled meta-contigs
     */
    map<int, Contig*>* m_nodes;

    /**
     * key = c1 meta-contig index
     *   value->key = c2 meta-contig index
     *     value->value = index of edge connecting c1 and c2
     */
    map<int, map<int, int> >* m_graph;
  };

}

#endif /* CONTIGNETWORK_H_ */
