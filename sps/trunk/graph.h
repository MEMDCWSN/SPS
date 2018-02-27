#ifndef GRAPH_H
#define GRAPH_H

#include "spectrum.h"
#include "aminoacid.h"
#include "SpectrumPairSet.h"
#include "label.h"
#include <vector>
#include <string>

namespace specnets
{
	using namespace std;

	/**
	 * Class for abstract graphs
	 *
	 */
	class Graph {
	public:
		/**
		 * Sets of destination vertices as edge indices.
		 */
		vector<vector<int> > vNext;

		/**
		 * Vertex scores
		 */
		vector<float> vScores;

		/**
		 * Edges as a set of (from,to) vertex indices pairs.
		 */
		vector<TwoValues<int> > edges;

		/**
		 * Edge scores
		 */
		vector<float> eScores;

		/**
		 * Outputs the graph in Graphviz for graphical rendering
		 *
		 *@param Output file name
		 */
		void output_graphviz(char *filename);

		/**
		 * TODO: add description
		 *
		 *@param aligns
		 */
		void build(SpectrumPairSet &aligns);

		/**
		 * TODO: add description
		 *
		 *@param aligns
		 */
		void add(SpectrumPairSet &aligns);

		/**
		 * TODO: add description
		 *
		 *@param pathLength
		 *@param minNumPaths
		 *@param vSets
		 *@param numPaths
		 */
		void findConnectedPairs(unsigned int pathLength, unsigned int minNumPaths,
				vector<vector<int> > &vSets, vector<int> &numPaths);

		/**
		 * Returns the number of vertices in the graph
		 *
		 *@return Number of vertices in the graph
		 */
		unsigned int numVertices() {
			return vNext.size();
		}

		/**
		 * Returns the number of edges in the graph
		 *
		 *@return Number of edges in the graph
		 */
		unsigned int numEdges() {
			return edges.size();
		}
	};

	/**
	 * Class for abstract graphs with additional mass spec info
	 */
	class MSGraph: public Graph {
	public:

		/**
		 * Amino acid masses per edge.
		 */
		vector<float> eMasses;

		/**
		 * Vertex labels.
		 */
		vector<string> vLabels;

		/**
		 * Peak labels.
		 */
		vector<PeakLabel> vPeakLabels;

		/**
		 * 1-based index of the edge in the heaviest path; >0 -> draw edge as bold.
		 */
		vector<short> edgeInPath;

		MSGraph() {
			eMasses.resize(0);
			vLabels.resize(0);
		}

		/**
		 * Returns the number of vertices in the graph
		 *
		 *@return Number of vertices in the graph
		 */
		unsigned int numVerts() {
			return vNext.size();
		}

		/**
		 * Builds graph from a spectrum by connecting pair of consecutive peaks with an edge.
		 *
		 *@param spec Input spectrum
		 */
		void ConnectConsecutive(Spectrum &spec);

		/**
		 * Builds graph from a spectrum by connecting pair of peaks whose
		 *   masses differ by a mass in jumps (within mass tolerance peakTol).
		 *
		 *@param spec Input spectrum
		 *@param jumps Set of masses used to connect peaks in the spectrum
		 *@param peakTol Tolerance for mass errors
		 */
		void ConnectJumps(Spectrum &spec, AAJumps &jumps, float peakTol);

		/**
		 * TODO: add description
		 *
		 */
		void consolidateEdgeScores();

		/**
		 * TODO: add description
		 *
		 *@param path
		 *@param useVertexScores
		 *@param resultSpec
		 *@param pathVertIdx
		 *@return
		 */
		float heaviestPath(MSGraph &path, bool useVertexScores = true,
				Spectrum *resultSpec = 0, vector<int> *pathVertIdx = 0);

		//	void findModPairs(SpecSet &specSet, float pmTol, list<vector<int> > &vSets);

		/**
		 * TODO: add description
		 *
		 *@param components
		 */
		void findConnected(vector<MSGraph> &components);

		/**
		 * TODO: add description
		 *
		 *@param parentMasses
		 *@param pmTol
		 *@param vSets
		 */
		void findModPairs(vector<int> &parentMasses, float pmTol,
				list<vector<int> > &vSets);

		/**
		 * TODO: add description
		 *
		 *@param filename the name of the file to output the graph to.
		 */
		void output_graphviz(char *filename);

		/**
		 * TODO: add description
		 *
		 *@param filename
		 */
		void output_clusters(char *filename);

		/**
		 * TODO: add description
		 *
		 *@param counts
		 */
		void info_heaviestPath(vector<float> &counts);
	};
}

#endif
