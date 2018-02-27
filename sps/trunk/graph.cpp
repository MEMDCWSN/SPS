#include "graph.h"
#include "SetMerger.h"
#include "clusters.h"
#include "Logger.h"
#include <cmath>
#include <iostream>
#include <fstream>

namespace specnets
{
	void Graph::build(SpectrumPairSet &aligns) {

		// Initialize edges/eScores + find maximum vertex index
		edges.resize(aligns.size());   eScores.resize(aligns.size());
		int maxVertexIdx=0;
		for(unsigned int i=0; i<aligns.size(); i++) {
			edges[i].set(aligns[i].spec1,aligns[i].spec2);
			eScores[i] = aligns[i].score1 + aligns[i].score2;
			maxVertexIdx=max(maxVertexIdx,aligns[i].spec1);
			maxVertexIdx=max(maxVertexIdx,aligns[i].spec2);
		}

		// initialize sizes in vNext
		vNext.resize(maxVertexIdx+1);   vScores.resize(maxVertexIdx+1);
		vector<int> countNext(vNext.size());
		for(unsigned int i=0; i<vNext.size(); i++) { countNext[i]=0; vScores[i]=0; }
		for(unsigned int i=0; i<edges.size(); i++) { countNext[edges[i][0]]++; }
		for(unsigned int i=0; i<vNext.size(); i++) { vNext[i].resize(countNext[i]); }

		// Fill vNext
		for(unsigned int i=0; i<vNext.size(); i++) { countNext[i]=0; }
		for(unsigned int i=0; i<edges.size(); i++) {
	//		vNext[edges[i][0]][countNext[edges[i][0]]] = edges[i][1];
			vNext[edges[i][0]][countNext[edges[i][0]]] = i;
			countNext[edges[i][0]]++;
		}
	}

	void Graph::add(SpectrumPairSet &aligns) {
		unsigned int newEdgeBase = edges.size();

		// Initialize edges/eScores + find maximum vertex index
		edges.resize(edges.size()+aligns.size());   eScores.resize(edges.size()+aligns.size());
		int maxVertexIdx=vNext.size();
		for(unsigned int i=0; i<aligns.size(); i++) {
			edges[newEdgeBase+i].set(aligns[i].spec1,aligns[i].spec2);
			eScores[i] = aligns[i].score1 + aligns[i].score2;
			maxVertexIdx=max(maxVertexIdx,aligns[i].spec1);
			maxVertexIdx=max(maxVertexIdx,aligns[i].spec2);
		}

		// Adapt sizes in vNext
		unsigned int newVertexBase = vNext.size();
		vNext.resize(maxVertexIdx+1);   vScores.resize(maxVertexIdx+1);
		vector<int> countNext(vNext.size());
		for(unsigned int i=0; i<newVertexBase; i++) countNext[i]=vNext[i].size();
		for(unsigned int i=newVertexBase; i<vNext.size(); i++) { countNext[i]=0; vScores[i]=0; }
		for(unsigned int i=newEdgeBase; i<edges.size(); i++) { countNext[edges[i][0]]++; }
		for(unsigned int i=0; i<vNext.size(); i++) { vNext[i].resize(countNext[i]); }

		// Fill vNext
		for(unsigned int i=0; i<vNext.size(); i++) { countNext[i]=0; }
		for(unsigned int i=0; i<edges.size(); i++) {
	//		vNext[edges[i][0]][countNext[edges[i][0]]] = edges[i][1];
			vNext[edges[i][0]][countNext[edges[i][0]]] = i;
			countNext[edges[i][0]]++;
		}
	}

	void Graph::findConnectedPairs(unsigned int pathLength, unsigned int minNumPaths,
									vector<vector<int> > &vSets, vector<int> &numPaths) {
		unsigned int numVerts = vNext.size();

		// Get sets of predecessors per vertex
		vector<vector<int> > vPreds(numVerts);
		vector<int> predCounts(numVerts);	for(unsigned int i=0; i<numVerts; i++) predCounts[i]=0;
		for(unsigned int i=0; i<edges.size(); i++) predCounts[edges[i][1]]++;
		for(unsigned int i=0; i<numVerts; i++) { vPreds[i].resize(predCounts[i]); predCounts[i]=0; }
		for(unsigned int i=0; i<edges.size(); i++) { vPreds[edges[i][1]][predCounts[edges[i][1]]++]=i; } //cerr<<"Pred of "<<edges[i][1]<<": ["<<edges[i][0]<<","<<edges[i][1]<<"]\n"; }
	/*	vector<vector<int> > vPreds(vNext.size());
		vector<int> predCounts(numVerts);	for(unsigned int i=0; i<numVerts; i++) predCounts[i]=0;
		for(unsigned int i=0; i<numVerts; i++)
			for(unsigned int j=0; j<vNext[i].size(); j++) predCounts[vNext[i][j]]++;
		for(unsigned int i=0; i<numVerts; i++) { vPreds[i].resize(predCounts[i]); predCounts[i]=0; }
		for(unsigned int i=0; i<numVerts; i++)
			for(unsigned int j=0; j<vNext[i].size(); j++)
				{ int next=vNext[i][j]; vPreds[next][predCounts[next]] = i; predCounts[next]++; }
	*/
		// Initialize component structures
		vector<vector<vector<TwoValues<int> > > > pathCounts;  // vertex idx [0], all preds [1], pathLengths [2], (pred idx,path count) [3]
		pathCounts.resize(numVerts);

		// Find connected pairs by looking backwards from each vertex
		unsigned int plIdx,       // path length index
								 predIdx,     // predecessor index (in the set of predecessors)
								 curVPred,    // predecessor index (in the set of vertices)
								 pathPredIdx; // index of path precessor vertex
		vector<vector<int> > curCounts; // Temporary values of predecessor path counts for current vertex. Zero-th column contains sum of all other lines (1st dimension in the array)
		curCounts.resize(numVerts); for(unsigned int i=0; i<numVerts; i++) curCounts[i].resize(pathLength+1);

		for(unsigned int v=0; v<numVerts; v++) {
			// initialize structures
			for(unsigned int i=0; i<numVerts; i++)
				for(plIdx=0; plIdx<=pathLength; plIdx++) curCounts[i][plIdx]=0;

	//cerr<<"--- Current vertex "<<v<<", preds = ";
	//for(unsigned int i=0; i<vPreds[v].size(); i++) cerr<<vPreds[v][i]<<" "; cerr<<endl;

			// Count number of paths from all current predecessors
			for(predIdx=0; predIdx<vPreds[v].size(); predIdx++) {
				curVPred = edges[vPreds[v][predIdx]][0];
				curCounts[curVPred][1] = 1;   curCounts[curVPred][0] = 1;
				for(unsigned int i=0; i<pathCounts[curVPred].size(); i++) {
					for(plIdx=2; plIdx<=pathLength; plIdx++) {  // If a vertex reaches a predecessor with a path of length k then it reached the current vertex with a path of length k+1
						curCounts[ pathCounts[curVPred][i][plIdx-1][0] ][plIdx] += pathCounts[curVPred][i][plIdx-1][1];
						curCounts[ pathCounts[curVPred][i][plIdx-1][0] ][0] += pathCounts[curVPred][i][plIdx-1][1];
					}
				}
			}

	//cerr<<"curCounts transpose, size "<<pathLength+1<<" x "<<numVerts<<endl;
	//for(unsigned int i=0; i<=pathLength; i++)
	//	{ for(unsigned int j=0; j<numVerts; j++) cerr<<curCounts[j][i]<<" "; cerr<<endl; }

			pathPredIdx=0;
			for(unsigned int i=0; i<numVerts; i++) if(curCounts[i][0]>0) pathPredIdx++;

			// Build persistent structures
			pathCounts[v].resize(pathPredIdx);   pathPredIdx=0;
			for(predIdx=0,pathPredIdx=0; predIdx<numVerts; predIdx++) {
				if(curCounts[predIdx][0]>0) {
					pathCounts[v][pathPredIdx].resize(pathLength+1);  // entry 0 contains sums across all path lengths
					for(plIdx=0; plIdx<=pathLength; plIdx++)
						{ pathCounts[v][pathPredIdx][plIdx].set(predIdx, curCounts[predIdx][plIdx] ); }
					pathPredIdx++;
				}
			}

	//cerr<<"pathCounts["<<v<<"] transpose, size "<<pathLength+1<<" x "<<pathCounts[v].size()<<endl;
	//for(unsigned int i=0; i<=pathLength; i++)
	//	{ for(unsigned int j=0; j<pathCounts[v].size(); j++) cerr<<"("<<pathCounts[v][j][i][0]<<","<<pathCounts[v][j][i][1]<<"] "; cerr<<endl; }

		}
		// Build vSets using two vectors
		//   successors from vertex1 with path length <= pathLength
		//   predecessors from vertex2 with path length <= pathLength
		//   vSets has all vertices that occur in both sets in a path of length pathLength
		pathPredIdx=0;   vector<bool> hasPairs(numVerts);
		for(unsigned int v=0; v<numVerts; v++)
			for(predIdx=0; predIdx<pathCounts[v].size(); predIdx++)
	//			if(pathCounts[v][predIdx][pathLength][1]>0 and pathCounts[v][predIdx][0][1]>=(int)minNumPaths) {  // minNumPaths includes all paths
				if(pathCounts[v][predIdx][pathLength][1]>=(int)minNumPaths)  // minNumPaths for pathLength only
					{ pathPredIdx++; hasPairs[v]=true; }
		vSets.resize(pathPredIdx);   numPaths.resize(pathPredIdx);   pathPredIdx=0;

	//cout << "Number of sets: " << vSets.size() << endl;

		TwoValues<int> pair;           // Pair of spectra that meet the pathLength/minNumPaths conditions
		vector<vector<bool> > succs(pathLength);  // Vertices found to be successors of pair[0]
		for(unsigned int i=0; i<pathLength; i++) succs[i].resize(numVerts);
		vector<list<int> > toProcess(pathLength);  // Stack of vertices left to process (per path length from base vertex)
		unsigned int commonCount;           // Number of vertices that are successors of pair[0] and predecessors of pair[1]
		vector<bool> common(vNext.size());  // Vertices that are successors of pair[0] and predecessors of pair[1]
		for(unsigned int i=1; i<numVerts; i++) common[i]==false;

		for(unsigned int v=0; v<numVerts; v++) {
			if(!hasPairs[v]) continue;
			for(predIdx=0; predIdx<pathCounts[v].size(); predIdx++)
	//			if(pathCounts[v][predIdx][pathLength][1]>0 and pathCounts[v][predIdx][0][1]>=(int)minNumPaths) {  // minNumPaths includes all paths
				if(pathCounts[v][predIdx][pathLength][1]>=(int)minNumPaths) {  // minNumPaths for pathLength only
					pair.set(pathCounts[v][predIdx][0][0],v);
					numPaths[pathPredIdx]=pathCounts[v][predIdx][0][1];

					// Find successors of pair[0] at most pathLength away
					toProcess[0].clear();
					for(unsigned int i=0; i<vNext[pair[0]].size(); i++) { toProcess[0].insert(toProcess[0].begin(),edges[vNext[pair[0]][i]][1]); }
					for(unsigned int i=0; i<pathLength; i++) for(unsigned int j=0; j<numVerts; j++) succs[i][j]=false;   //succs[pair[0]]=true;
					for(plIdx=0; plIdx<pathLength; plIdx++) {
						if(plIdx<pathLength-1) toProcess[plIdx+1].clear();
						for(list<int>::iterator iter=toProcess[plIdx].begin(); iter!=toProcess[plIdx].end(); iter++) {
							if(succs[plIdx][*iter]) continue; else succs[plIdx][*iter]=true;
							if(plIdx<pathLength-1)
								for(unsigned int i=0; i<vNext[*iter].size(); i++)
									{ toProcess[plIdx+1].insert(toProcess[plIdx+1].begin(),edges[vNext[*iter][i]][1]); }
						}
					}

					// Find common successors/predecessors and build vSets accordingly
					commonCount=0;    common[pair[0]]=true;   common[pair[1]]=true;
					for(unsigned int i=1; i<pathCounts[v].size(); i++) {
						for(unsigned int j=0; j<pathLength; j++)
							if(succs[pathLength-1-j][pathCounts[v][i][j][0]] and !common[pathCounts[v][i][j][0]])
								{ commonCount++; common[pathCounts[v][i][j][0]]=true; break; }
					}

					vSets[pathPredIdx].resize(commonCount+2);   commonCount=2;
					vSets[pathPredIdx][0]=pair[0];              common[pair[0]]=false;
					vSets[pathPredIdx][1]=pair[1];              common[pair[1]]=false;
	//cerr<<"vSets["<<pathPredIdx<<"], size() = "<<vSets[pathPredIdx].size()<<": "<<pair[0]<<" "<<pair[1]<<" ";
					for(unsigned int i=1; i<numVerts; i++)
						if(common[i]) {
							vSets[pathPredIdx][commonCount++]=i;
							common[i]=false;
	//cerr<<i<<" ";
						}
	//cerr<<endl;
	/*				for(unsigned int i=0; i<pathCounts[v].size(); i++) {
						if(pathCounts[v][i][0][0]==pair[0] or pathCounts[v][i][0][0]==pair[1]) continue;
						for(unsigned int j=0; j<pathLength; j++)
							if(succs[pathLength-1-j][pathCounts[v][i][j][0]]) {
								vSets[pathPredIdx][commonCount++]=pathCounts[v][i][j][0];
								break;
							}
					}
	*/
					if(pathPredIdx>0 and pathPredIdx%1000==0)
					  DEBUG_MSG("Generated vSets["<<pathPredIdx<<"]...");
					pathPredIdx++;
				}
		}
	}

	void Graph::output_graphviz(char *filename){
		ofstream out(filename, ios::binary);
		if(!out) { ERROR_MSG("Error opening "<<filename<<"!"); return; }
		out << "digraph G {\n";

		for(int eIdx=0; eIdx<(int)edges.size(); eIdx++)
			out << edges[eIdx][0]<<" -> "<<edges[eIdx][1]<<" [label = \""<<eScores[eIdx]<<"\"];\n";

		out << "}\n";
		out.close();
	}

	void MSGraph::ConnectConsecutive(Spectrum &spec){
		unsigned int numPeaks = spec.size();
		vNext.resize(numPeaks);     vScores.resize(numPeaks);
		if(numPeaks>0) { edges.resize(numPeaks-1); eMasses.resize(numPeaks-1); eScores.resize(numPeaks-1); }
		else { edges.resize(0); eMasses.resize(0); eScores.resize(0); return; }

		// Connect spectrum peaks
		for(unsigned int p=0; p<numPeaks-1; p++) {
			vNext[p].resize(1);    vNext[p][0]=p;     vScores[p]=spec[p][1];
			edges[p].set(p,p+1);   eMasses[p]=spec[p+1][0]-spec[p][0];
			eScores[p] = spec[p+1][1];   if(p==0) eScores[p]+=vScores[p];
		}
		vNext[numPeaks-1].resize(0);  vScores[numPeaks-1]=spec[numPeaks-1][1];
	}

	void MSGraph::ConnectJumps(Spectrum &spec, AAJumps &jumps, float peakTol){
		unsigned int numPeaks = spec.size(), numEdges=0;
		if(numPeaks==0) return;

		vector<int> inDegree(numPeaks); for(unsigned int i=0;i<numPeaks;i++) inDegree[i]=0;
		vNext.resize(numPeaks);     vScores.resize(numPeaks);
		edges.resize(numPeaks*(numPeaks-1));   eMasses.resize(numPeaks*(numPeaks-1));   eScores.resize(numPeaks*(numPeaks-1));

		// Connect spectrum peaks
		for(unsigned int p=0; p<numPeaks-1; p++) {
			unsigned int jumpIdx=0, nextPeak=p+1, numSuccs=0;
			vNext[p].resize(numPeaks-1-p);
			while (nextPeak<numPeaks) {
			// Find smallest jump that connects the two peaks
				while (jumpIdx>0 and spec[nextPeak][0]-spec[p][0]-jumps[jumpIdx]<-peakTol) jumpIdx--;
				if(spec[nextPeak][0]-spec[p][0]-jumps[jumpIdx]<-peakTol) // Peaks are too close, no jump applies
					{ nextPeak++; continue; }
				while (jumpIdx<jumps.size() and spec[nextPeak][0]-spec[p][0]-jumps[jumpIdx]>peakTol) jumpIdx++;
				if(jumpIdx>jumps.size()) break; // Peaks are too far apart, no jump can connect them

				if(abs(spec[nextPeak][0]-spec[p][0]-jumps[jumpIdx])<=peakTol) {
					vNext[p][numSuccs++] = numEdges;
					edges[numEdges].set(p,nextPeak);   inDegree[nextPeak]++;
					eMasses[numEdges] = spec[nextPeak][0]-spec[p][0];  // eMasses[numEdges] = jumps[jumpIdx];
					eScores[numEdges] = spec[nextPeak][1];  if(inDegree[p]==0) eScores[numEdges]+=spec[p][1];
					numEdges++;
				}
				nextPeak++;
			}
			vNext[p].resize(numSuccs);   vScores[p]=spec[p][1];
		}
		edges.resize(numEdges);   eMasses.resize(numEdges);   eScores.resize(numEdges);
		vNext[numPeaks-1].resize(0);  vScores[numPeaks-1]=spec[numPeaks-1][1];
	}

	//
	//  consolidateEdgeScores - replaces each edge's score by eScore*vScore[destination vertex].
	//    A bonus of vScore[source vertex] is given when the source vertex is a source (no incoming edges)
	//
	void MSGraph::consolidateEdgeScores() {
		vector<int> numPreds(vNext.size());

		for(unsigned int i=0; i<numPreds.size(); i++) numPreds[i]=0;
		for(unsigned int i=0; i<edges.size(); i++) numPreds[edges[i][1]]++;
		for(unsigned int i=0; i<edges.size(); i++) {
			eScores[i] = eScores[i]*max((double)vScores[edges[i][1]],0.01);
			if(numPreds[edges[i][0]]==0) eScores[i]+=vScores[edges[i][0]];
		}
	}

	//
	//  heaviestPath - finds the heaviest path in the graph and returns it in path.
	//    return value is the score of the heaviest path.
	//
	//  NOTE: Only considers weights on edges. Use consolidateEdgeScores to include
	//         vertex weights into edge weights.
	//
	float MSGraph::heaviestPath(MSGraph &path, bool useVertexScores, Spectrum *resultSpec, vector<int> *pathVertIdx) {
		vector<TwoValues<int> > predsToProcess(vNext.size()); // Number of predecessor edges left to process, Total number of predecessor edges
		vector<TwoValues<int> > bestPred(vNext.size());  // Best predecessor edge index (pos.0), Path lenght (col.1)
		vector<float> bestPathScore(vNext.size());       // Score of the highest scoring path ending at the vertex
		list<int> sources;                               // List of source vertices awaiting processing
		vector<float> eScoresHP(edges.size());           // Edge scores - may or may not combine edge/vertex scores
		unsigned int processedVertices = 0,              // Number of processed vertices
			bestPathIdx = 0;                             // Index of the vertex where the highest scoring path ends

		if(edges.size()==0) {
			path.edges.resize(0);   path.eScores.resize(0);   path.eMasses.resize(0);
			path.vNext.resize(0);   path.vScores.resize(0);
			if(resultSpec) resultSpec->resize(0);   if(pathVertIdx) pathVertIdx->resize(0);
			return 0;
		}

		// Initialize structures
		for(unsigned int i=0; i<vNext.size(); i++) {
			predsToProcess[i].set(0,0);   bestPred[i].set(-1,0);   bestPathScore[i] = 0;
		}
		for(unsigned int i=0; i<edges.size(); i++)
			if(edges[i][0]!=edges[i][1])   // Ignore one vertex cycles
				{ predsToProcess[edges[i][1]][0]++; predsToProcess[edges[i][1]][1]++; }
		for(unsigned int i=0; i<vNext.size(); i++)
			if(predsToProcess[i][0]==0) { sources.push_front(i); processedVertices++; }
		if(edgeInPath.size()!=edges.size()) edgeInPath.resize(edges.size());
		if(!useVertexScores) for(unsigned int i=0; i<edges.size(); i++) { eScoresHP[i] = eScores[i]; edgeInPath[i]=0; }
		else for(unsigned int i=0; i<edges.size(); i++) {
			eScoresHP[i] = max((double)vScores[edges[i][1]],0.01);
	//		eScoresHP[i] = eScores[i]*max((double)vScores[edges[i][1]],0.01);
			if(predsToProcess[edges[i][0]][1]==0) eScoresHP[i]+=vScores[edges[i][0]];
			edgeInPath[i]=0;
		}

		while (sources.size()>0) {
			int curVert=sources.front();   sources.pop_front();
			if(bestPathScore[curVert]>bestPathScore[bestPathIdx]) bestPathIdx=curVert;

			// Propagate paths to all successors of curVert
			for(unsigned int succIdx=0; succIdx<vNext[curVert].size(); succIdx++) {
				int edgeIdx = vNext[curVert][succIdx], nextVert = edges[edgeIdx][1];
				if(edges[edgeIdx][0]==edges[edgeIdx][1]) continue;  // avoid single vertex cycles
				if(predsToProcess[nextVert][0]==0) continue;        // nextPred was processed prematurely to break a cycle
				if(bestPathScore[curVert]+eScoresHP[edgeIdx] > bestPathScore[nextVert]) {
					bestPred[nextVert].set(edgeIdx,bestPred[curVert][1]+1);
					bestPathScore[nextVert] = bestPathScore[curVert]+eScoresHP[edgeIdx];
				}
				predsToProcess[nextVert][0]--;
				if(predsToProcess[nextVert][0]==0) { sources.push_front(nextVert); processedVertices++; }
			}

			if(sources.size()==0 and processedVertices!=vNext.size()) {
				float minPredRatio=1;  int minPredIdx=-1; // minPredRatio is percentage of unresolved preds and minPredIdx is vertex with smallest minPredRatio
				for(unsigned int i=0; i<predsToProcess.size(); i++)
					if(predsToProcess[i][0]>0 and (((float)predsToProcess[i][0])/((float)predsToProcess[i][1]))<=minPredRatio)
						{ minPredRatio = (((float)predsToProcess[i][0])/((float)predsToProcess[i][1])); minPredIdx=i; }
				DEBUG_MSG("WARNING: heaviestPath may perform incorrectly if the graph contains cycles! Unblocked vertex "<<minPredIdx<<" with "<<predsToProcess[minPredIdx][0]<<" predecessors unresolved (minPredRatio = "<<minPredRatio<<").");
				sources.push_front(minPredIdx); processedVertices++;   predsToProcess[minPredIdx][0]=0;
			}
		}

		// Build heaviest path graph
		int pathPivot = bestPred[bestPathIdx][0], pathVert=0, pathEdge=0;
		if(pathPivot<0) return -1;
		int curVert = edges[pathPivot][0], succVert = edges[pathPivot][1];

		path.edges.resize(bestPred[bestPathIdx][1]);    path.eScores.resize(path.edges.size());   path.eMasses.resize(path.edges.size());
		path.vNext.resize(bestPred[bestPathIdx][1]+1);  path.vScores.resize(path.vNext.size());
		if(vLabels.size()>0) path.vLabels.resize(path.vNext.size());
		if(vPeakLabels.size()>0) path.vPeakLabels.resize(path.vNext.size());
		if(resultSpec) resultSpec->resize(path.vNext.size());
		if(pathVertIdx) pathVertIdx->resize(path.vNext.size());

		path.vNext[pathVert].resize(0);                 path.vScores[pathVert] = vScores[succVert];
		if(vLabels.size()>0) path.vLabels[pathVert] = vLabels[succVert];
		if(vPeakLabels.size()>0) path.vPeakLabels[pathVert] = vPeakLabels[succVert];
		if(resultSpec) (*resultSpec)[path.vNext.size()-1].set(eMasses[pathPivot],path.vScores[pathVert]);
		if(pathVertIdx) (*pathVertIdx)[path.vNext.size()-1] = succVert;
		pathVert++;
		while(pathPivot>=0) {
			edgeInPath[pathPivot] = path.vNext.size()-pathVert;
	//cerr<<"pathPivot = "<<pathPivot<<", current edge = ("<<curVert<<","<<edges[pathPivot][1]<<"), bestPred = ["<<bestPred[curVert][0]<<","<<bestPred[curVert][1]<<"]\n";
			path.vNext[pathVert].resize(1);               path.vNext[pathVert][0]=pathEdge;
			path.vScores[pathVert] = vScores[curVert];
			if(vLabels.size()>0) path.vLabels[pathVert] = vLabels[curVert];
			if(vPeakLabels.size()>0) path.vPeakLabels[pathVert] = vPeakLabels[curVert];

			path.edges[pathEdge].set(pathVert,pathVert-1);
			path.eScores[pathEdge] = eScores[pathPivot];  path.eMasses[pathEdge] = eMasses[pathPivot];

			if(pathVertIdx) (*pathVertIdx)[path.vNext.size()-1-pathVert] = curVert;

			pathVert++;   pathEdge++;   pathPivot = bestPred[curVert][0];   if (pathPivot>=0) curVert = edges[pathPivot][0];
			if(resultSpec) {
				if (pathPivot>=0) (*resultSpec)[path.vNext.size()-pathVert].set(eMasses[pathPivot],path.vScores[pathVert-1]);
				else (*resultSpec)[path.vNext.size()-pathVert].set(0,path.vScores[pathVert-1]); }
		}

		if(resultSpec) {
			for(unsigned int i=1;i<resultSpec->size();i++) (*resultSpec)[i][0]+=(*resultSpec)[i-1][0];  // Cumulative sum of the path edge masses
			resultSpec->parentMass = (*resultSpec)[resultSpec->size()-1][0]+AAJumps::massMH;
		}
		return bestPathScore[bestPathIdx];
	}

	void MSGraph::findConnected(vector<MSGraph> &components) {
		SetMerger helper(vNext.size());
		for(unsigned int v=0; v<vNext.size(); v++) helper.createset(v);  // Initialize with one-vertex sets
		for(unsigned int e=0; e<edges.size(); e++)
			helper.merge(helper.membership[edges[e][0]],helper.membership[edges[e][1]]);
		helper.compressSetIndices();

		components.resize(helper.numSets);   vector<int> edgeCount(helper.numSets), vCount(helper.numSets);
		vector<int> vRename(vNext.size());  // Index conversion table from current vertex indices to indices in the new components
		// Count #vertices/#edges per component
		for(unsigned int cIdx=0; cIdx<components.size(); cIdx++) { edgeCount[cIdx]=0; vCount[cIdx]=0; }
		for(unsigned int e=0; e<edges.size(); e++) edgeCount[helper.membership[edges[e][0]]]++;
		for(unsigned int v=0; v<vNext.size(); v++) vRename[v]=vCount[helper.membership[v]]++;

		// Allocate space for the new graphs' data
		unsigned int setIdx;
		for(setIdx=0; setIdx<components.size(); setIdx++) {
			components[setIdx].edges.resize(edgeCount[setIdx]);
			if(eScores.size()>0) components[setIdx].eScores.resize(edgeCount[setIdx]);
			if(eMasses.size()>0) components[setIdx].eMasses.resize(edgeCount[setIdx]);
			if(edgeInPath.size()>0) components[setIdx].edgeInPath.resize(edgeCount[setIdx]);
			edgeCount[setIdx]=0;

			components[setIdx].vNext.resize(vCount[setIdx]);
			if(vScores.size()>0) components[setIdx].vScores.resize(vCount[setIdx]);
			if(vLabels.size()>0) components[setIdx].vLabels.resize(vCount[setIdx]);
			if(vPeakLabels.size()>0) components[setIdx].vPeakLabels.resize(vCount[setIdx]);
		}

		// Copy vertex/edges data to new graphs
		for(unsigned int v=0; v<vNext.size(); v++) {
			setIdx = helper.membership[v]; int vIdx = vRename[v];
			if(vScores.size()>0) components[setIdx].vScores[vIdx] = vScores[v];
			if(vLabels.size()>0) components[setIdx].vLabels[vIdx] = vLabels[v];
			if(vPeakLabels.size()>0) components[setIdx].vPeakLabels[vIdx] = vPeakLabels[v];

			// copy edges
			components[setIdx].vNext[vIdx].resize(vNext[v].size());  // Vertex degree does not change.
			for(unsigned int e=0; e<vNext[v].size(); e++) {
				int eIdx=edgeCount[setIdx]++;
				components[setIdx].edges[eIdx].set(vRename[edges[vNext[v][e]][0]],vRename[edges[vNext[v][e]][1]]);
				if(eScores.size()>0) components[setIdx].eScores[eIdx] = eScores[vNext[v][e]];
				if(eMasses.size()>0) components[setIdx].eMasses[eIdx] = eMasses[vNext[v][e]];
				if(edgeInPath.size()>0) components[setIdx].edgeInPath[eIdx] = edgeInPath[vNext[v][e]];
			}
		}
	}


	//
	//  findModPairs - find pairs of spectra A,B such that:
	//    C=succ(A)=pred(B), D=succ(A)=pred(B) and PM(D)+PM(A)=PM(C)+PM(D)
	//    Every pair A,B generates with entry in vSets where the first two elements
	//    are A,B and the remaining elements are all possible C,D vertices.
	//
	//  NOTE: Assumes that specSet[i].parentMass<=specSet[j].parentMass, for all i<=j
	//
	void MSGraph::findModPairs(vector<int> &parentMasses, float pmTol, list<vector<int> > &vSets) {
		unsigned int numVerts = vNext.size();

		int maxMassInt=parentMasses[0]; for(unsigned int i=1; i<numVerts; i++) if(maxMassInt<parentMasses[i]) maxMassInt=parentMasses[i];
		unsigned int pmTolInt = (int)ceil(10*pmTol);   maxMassInt += pmTolInt+1;

		vector<bool> existPM(maxMassInt);    // [i]==1 if spectra with this PM exist between A and B
		for(unsigned int i=0; i<maxMassInt; i++) existPM[i]=false;

		vector<bool> succs(numVerts);  // [i]==1 if i is a sucessor of the current vertex
		vector<bool> processed(numVerts);  // [i]==1 if i has already been considered as a pair (of baseV)


		// Get sets of predecessors per vertex
		vector<vector<int> > vPreds(numVerts);
		vector<int> predCounts(numVerts);	for(unsigned int i=0; i<numVerts; i++) predCounts[i]=0;
		for(unsigned int i=0; i<edges.size(); i++) predCounts[edges[i][1]]++;
		for(unsigned int i=0; i<numVerts; i++) { vPreds[i].resize(predCounts[i]); predCounts[i]=0; }
		for(unsigned int i=0; i<edges.size(); i++) { vPreds[edges[i][1]][predCounts[edges[i][1]]++]=i; } //cerr<<"Pred of "<<edges[i][1]<<": ["<<edges[i][0]<<","<<edges[i][1]<<"]\n"; }


		unsigned int baseV, succVidx, succV, pairIdx, pairV;
		int otherPM;
		list<int> siblings;
		vector<int> curSet;
		for(baseV=0; baseV<numVerts; baseV++) {
			// Initialize successor structures
			for(unsigned int i=0; i<numVerts; i++) processed[i]=false;
			for(succVidx=0; succVidx<vNext[baseV].size(); succVidx++) {
				succV = edges[vNext[baseV][succVidx]][1];    succs[succV]=true;
	//cerr << "succ "<<succV<<", "<<parentMasses[succV]<<", "<<existPM[parentMasses[succV]]<<endl;
			}

			// Find pairs
			for(succVidx=0; succVidx<vNext[baseV].size(); succVidx++) {
				succV = edges[vNext[baseV][succVidx]][1];

				for(pairIdx=0; pairIdx<vNext[succV].size(); pairIdx++) {
					pairV = edges[vNext[succV][pairIdx]][1];
					if(processed[pairV]) continue;
	//cerr << "Testing a pair: baseV = "<<baseV<<", succV = "<<succV<<", pairV = "<<pairV<<"\n";
					// Initialize data structures
					unsigned int siblingIdx, siblingV;
					siblings.clear();   // List of matched siblings
					// Find the parent masses of every eligible sibling
					for(siblingIdx=0; siblingIdx<vPreds[pairV].size(); siblingIdx++) {
						siblingV = edges[vPreds[pairV][siblingIdx]][0];   if(!succs[siblingV]) continue;
						for(unsigned int i=max(0,parentMasses[siblingV]-(int)pmTolInt); i<=parentMasses[siblingV]+pmTolInt; i++) existPM[i]=true;
					}

					for(siblingIdx=0; siblingIdx<vPreds[pairV].size(); siblingIdx++) {
						siblingV = edges[vPreds[pairV][siblingIdx]][0];
						otherPM = parentMasses[baseV]+parentMasses[pairV]-parentMasses[siblingV];
	//cerr<<"siblingV = "<<siblingV<<", "<<parentMasses[baseV]<<", "<<parentMasses[pairV]<<", "<<parentMasses[succV]<<", "<<otherPM<<"... ";
	//					if(!succs[siblingV] or !existPM[otherPM]) { cerr<<"not added: "<<succs[siblingV]<<"/"<<otherPM<<"/"<<existPM[(int)round(10*otherPM)]<<"\n"; continue; } // There is no adequate sibling
						if(!succs[siblingV] or !existPM[otherPM]) continue; // There is no adequate sibling
						siblings.insert(siblings.begin(),siblingV);
	//cerr<<"added.\n";
					}

					if(siblings.size()>1) {  // Report new pair (baseV,pairV)
						list<int>::iterator iter=siblings.begin();    curSet.resize(siblings.size()+2);
						curSet[0]=baseV;  curSet[1]=pairV;
						for(unsigned int i=2; i<curSet.size(); i++)	{ curSet[i]=*iter; iter++; }
						vSets.insert(vSets.begin(),curSet);
					}

					processed[pairV]=true;
					for(siblingIdx=0; siblingIdx<vPreds[pairV].size(); siblingIdx++) {
						siblingV = edges[vPreds[pairV][siblingIdx]][0];   if(!succs[siblingV]) continue;
						for(unsigned int i=max(0,parentMasses[siblingV]-(int)pmTolInt); i<=parentMasses[siblingV]+pmTolInt; i++) existPM[i]=false;
					}
				}
			}

			// reset successor structures
			for(succVidx=0; succVidx<vNext[baseV].size(); succVidx++) {
				succV = edges[vNext[baseV][succVidx]][1];    succs[succV]=false;
			}
		}
	}

	void MSGraph::output_graphviz(char *filename){
		char *sBuf = (char *)malloc(1024);   if(sBuf==(char *)0) { ERROR_MSG("ERROR: Out of memory (MSGraph::output_graphviz)"); return; }
		ofstream out(filename, ios::binary);
		if(!out) { ERROR_MSG("Error opening "<<filename<<"!"); return; }
		out << "digraph G {\n";

		vector<string> *labels;
		if(vLabels.size()>0) labels=&vLabels; else {
			labels = new vector<string>;   labels->resize(vNext.size());
			for(unsigned int i=0; i<vNext.size(); i++) { sprintf(sBuf,"%d",i); (*labels)[i] = string(sBuf); }
		}

		for(unsigned int i=0; i<vPeakLabels.size(); i++) {
			out << (*labels)[i];
			if(!vPeakLabels[i].isEmpty()) {
				out << "\t[";
				if(vPeakLabels[i].nb()>0 and vPeakLabels[i].ny()==0) out<<" shape=diamond,color=green,style=filled, ";
				if(vPeakLabels[i].nb()==0 and vPeakLabels[i].ny()>0) out<<" shape=triangle,color=yellow,style=filled, ";
				if(vPeakLabels[i].nb()>0 and vPeakLabels[i].ny()>0) {
					out<<" shape=polygon,sides=4,";
					if(vPeakLabels[i].nb()>vPeakLabels[i].ny()) out<<"color=green,style=filled,";
					if(vPeakLabels[i].nb()<vPeakLabels[i].ny()) out<<"color=yellow,style=filled,";
				}
				out << " label=\"" << (*labels)[i] <<" : "<< vPeakLabels[i].asString() << "\" ]\n";
			} else out << "\t[label=\""<<(*labels)[i]<<" : "<<vScores[i]<<"\"];\n";
		}

	//fprintf(fid,'\t %d -> %d [ color = green, style = bold, label = %.1f ];\n', aligns(idxOk(i),[1 2 3])); end;

	//	for(unsigned int i=0; i<vLabels.size(); i++) out << i << " [ label = \""<<vLabels[i]<<"\"]\n";
		string styleStr;
		for(unsigned int eIdx=0; eIdx<edges.size(); eIdx++) {
	//		out << edges[eIdx][0]<<" -> "<<edges[eIdx][1]<<" [label = "<<eMasses[eIdx]<<"];\n";
			if(eIdx<edgeInPath.size() and edgeInPath[eIdx]) { styleStr.assign(",style=bold,color=red"); sprintf(sBuf,"%d: ",edgeInPath[eIdx]); } else { styleStr.assign(""); sBuf[0]=0; }
			if(eIdx<eMasses.size())
				out << (*labels)[edges[eIdx][0]]<<" -> "<<(*labels)[edges[eIdx][1]]<<" [label = \"("<<sBuf<<eMasses[eIdx]<<"/"<<eScores[eIdx]<<")\""<<styleStr<<"];\n";
			else out << (*labels)[edges[eIdx][0]]<<" -> "<<(*labels)[edges[eIdx][1]]<<" [label = \"("<<sBuf<<eScores[eIdx]<<" )\""<<styleStr<<"];\n";
		}

		out << "}\n";
		out.close();

		if(labels!=&vLabels) delete labels;
		free(sBuf);
	}

	//
	//  NOTE: Assumes the graph to be of the format constructed and returned by heaviestPath
	//
	void MSGraph::output_clusters(char *filename){
		unsigned int numPeaks = vNext.size();
		Clusters clst;

	//	clst.specIdx.resize(1);    clst.specIdx[0].resize(1);   clst.specIdx[0][0]=0;   // This special case only contains consensus PRMs
	//	clst.shifts.resize(1);     clst.shifts[0].resize(1);    clst.shifts[0][0].set(0,0);
		clst.specIdx.resize(1);    clst.specIdx[0].resize(0);
		clst.shifts.resize(1);     clst.shifts[0].resize(0);
		clst.endpoints.resize(1);  clst.endpoints[0].resize(0);
		clst.consensus.resize(1);  clst.consensus[0].resize(numPeaks);

		float curPeakMass=0;
		for(unsigned int i=0; i<numPeaks-1; i++) {
			clst.consensus[0][i].set(curPeakMass,vScores[numPeaks-1-i]);
			curPeakMass += eMasses[vNext[numPeaks-1-i][0]];
		}
		clst.consensus[0][numPeaks-1].set(curPeakMass,vScores[0]);
		clst.consensus[0].parentMass = curPeakMass+19;

		clst.Save(filename);
	}

	//
	// info_heaviestPath - counts =
	//   (%path_vertices with #b>#y, %path_vertices with #b<#y, %path_vertices endpoints>0 and with #b == #y == 0,
	//   #b in path, total #b, #y in path, total #y, #endpoint peaks in path, total #peaks in path) -> counts positions 0-8
	//
	void MSGraph::info_heaviestPath(vector<float> &counts) {
		counts.resize(9); for(unsigned int i=0; i<9; i++) counts[i]=0;
		if(vPeakLabels.size()==0) return;

		vector<bool> vertexInPath(vNext.size());
		for(unsigned int i=0; i<vNext.size(); i++)
			{ vertexInPath[i]=false;  counts[4]+=vPeakLabels[i].nb();  counts[6]+=vPeakLabels[i].ny(); }

		for(unsigned int e=0; e<edges.size(); e++)
			if(edgeInPath[e]) { vertexInPath[edges[e][0]]=true; vertexInPath[edges[e][1]]=true; }

		float vCount=0; // # vertices in heaviest path
		for(unsigned int i=0; i<vNext.size(); i++)
			if(vertexInPath[i]) {
				vCount++;
				if(vPeakLabels[i].nb()>vPeakLabels[i].ny()) counts[0]++;
				if(vPeakLabels[i].nb()<vPeakLabels[i].ny()) counts[1]++;
				if(vPeakLabels[i].ne()>0 and vPeakLabels[i].nb()==0 and vPeakLabels[i].ny()==0) counts[2]++;
				counts[3]+=vPeakLabels[i].nb();   counts[5]+=vPeakLabels[i].ny();   counts[7]+=vPeakLabels[i].ne();
				counts[8]+=vPeakLabels[i].nb()+vPeakLabels[i].ny()+vPeakLabels[i].no()+vPeakLabels[i].nm()+vPeakLabels[i].ne();
			}
		counts[0] = counts[0]/vCount;   counts[1] = counts[1]/vCount;   counts[2] = counts[2]/vCount;
	}
}
