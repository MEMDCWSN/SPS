#include "inputParams.h"
#include "alignment_scoring.h"
#include "batch.h"
#include "filters.h"
#include "abruijn.h"
#include "graph.h"
#include "SetMerger.h"
//#include "spectral_pairs.h"
#include "SpectralPairs.h"
#include "SpectrumPairSet.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;
using namespace specnets;

/*
  masab - Indentifies sets of spectra connected by spectral alignment, builds an
    A-Bruijn graph for every connected component and finds a consensus de novo
    sequence for each component.

  INPUT:
    INPUT_SPECS     - Input deconvolved spectra (mostly b or mostly y but not both)
	INPUT_ALIGNS    - Set of detected spectral alignment pairs
	INPUT_ALIGNSPA  - Set of detected partial-overlap spectral alignment pairs
  ?  INPUT_LABELS    -
    SPEC_TYPE_MSMS  - Spectrum type: 0 (PRM), 1 (MS/MS)
	TOLERANCE_PEAK  - Peak mass tolerance (in Daltons, defaul 0.5 Da)
    TOLERANCE_PM    - Parent mass tolerance (in Daltons, default 1 Da)
	MIN_RATIO       - Minimum matched score in aligned spectra to accept a spectra alignment
	PENALTY_PTM     - Score penalty for matching a PTM in a spectral alignment
	PENALTY_SAME_VERTEX - Penalty for using complementary peaks (i.e. b/y) in the same alignment.
	                        Usually '-infinity' , i.e. -1000000000
	MAX_AA_JUMP     - Largest mass jump in the ABruijn graphs (usually 2).
    EDGE_SCORE_TYPE -
	            EST_EDGE_MULT = 0,       // Edge multiplicity
	            EST_EDGE_SCORES = 1,     // Edge scores: each edge collects its score from the
	                                                     destination peak in the original spectrum (default)
	            EST_ABVERTEX_SCORES = 2; // Vertex scores: each edge collects its score from the destination A-Bruijn vertex
	                                     //   Edge multiplicity is added to vertex scores to help distinguish between edges
    GRAPH_TYPE - ABruijn graph edges are derived by gluing
                0 - path graphs from sets of consecutive matched peaks
                1 - path graphs (each spectrum is represented as a path through _all_ of its peaks)
                2 - spectrum graphs (one per spectrum) - default option
	PATH_MIN_PEAKS - Minimum number of peaks for a spectrum to be part of a protein contig sequence
	PATH_MIN_SPECS - Minimum number of spectra (observing PATH_MIN_PEAKS) for a valid protein contig sequence

  OUTPUT:
	OUTPUT_SPECS  -
	OUTPUT_MODPOS -
	OUTPUT_RATIOS   -
	OUTPUT_RATIOSPA -
*/

// Items that need fixing to avoid having to split aligns into ASP/PA
// - MSGraph g;   g.build(alignsASP);   g.add(cAligns[cIdx]); - should be able to initialize using alignsPA
// - VertexSet::addEndpointEdges - only input is alignsASP?
// - SetMerger::splitSet(cIdx,usedSpecs) - ASP/PA cAligns/cAligns_idx are class members

void SplitASPPA(SpectrumPairSet &pairs, float pmTol, float maxModMass, SpectrumPairSet &pairsASP, SpectrumPairSet &pairsPA) {
	pairsASP.resize(pairs.size());   pairsPA.resize(pairs.size());
	unsigned int idxPair, idxASP=0, idxPA=0;
	for(idxPair=0; idxPair<pairs.size(); idxPair++)
		if((fabs(pairs[idxPair].shift1)<=pmTol and fabs(pairs[idxPair].shift2)<=maxModMass) or (fabs(pairs[idxPair].shift1)<=maxModMass and fabs(pairs[idxPair].shift2)<=pmTol))
			pairsASP[idxASP++]=pairs[idxPair]; else pairsPA[idxPA++]=pairs[idxPair];
	pairsASP.resize(idxASP);   pairsPA.resize(idxPA);
}

int main(int argc, char **argv){
    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("masab.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "masab.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs;   paramStrs.resize(6);
	paramStrs[0] = "OUTPUT_SPECS";
	paramStrs[1] = "OUTPUT_MODPOS";
	paramStrs[2] = "PENALTY_PTM";
	paramStrs[3] = "PENALTY_SAME_VERTEX";
	paramStrs[4] = "GRAPH_TYPE";
	paramStrs[5] = "INPUT_ALIGNS";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"masab.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_ALIGNS, OUTPUT_SPECS, OUTPUT_MODPOS, PENALTY_PTM, PENALTY_SAME_VERTEX, GRAPH_TYPE\n";
		return -1;
	}

	const char *alignsFN = params.getValue("INPUT_ALIGNS");
	const char *resultsFN = params.getValue("OUTPUT_SPECS");
	const char *modPosFN = params.getValue("OUTPUT_MODPOS");
	float penalty_ptm = (float) params.getValueDouble("PENALTY_PTM");
	float penalty_sameVert = (float) params.getValueDouble("PENALTY_SAME_VERTEX");
	int graphType = params.getValueInt("GRAPH_TYPE");

	const char *labelsFN = params.getValue("INPUT_LABELS");
	if(params.paramPresent("AMINO_ACID_MASSES")) { AAJumps jumps(-1); jumps.loadJumps(params.getValue("AMINO_ACID_MASSES"),true); }
	int maxAAjump = params.paramPresent("MAX_AA_JUMP")?params.getValueInt("MAX_AA_JUMP"):0;
	float maxModMass = params.paramPresent("MAX_MOD_MASS")?(float) params.getValueDouble("MAX_MOD_MASS"):100.0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1.0;
    short edgeScoreType = (short)(params.paramPresent("EDGE_SCORE_TYPE")?params.getValueInt("EDGE_SCORE_TYPE"):1);
    unsigned int minMatchedPeaks = params.paramPresent("MIN_MATCHED_PEAKS")?(int) params.getValueInt("MIN_MATCHED_PEAKS"):1;
    unsigned int minEdgesToComponent = params.paramPresent("MIN_EDGES_TO_COMPONENT")?(int) params.getValueInt("MIN_EDGES_TO_COMPONENT"):1;
    unsigned int pathMinSpecs = (unsigned int)(params.paramPresent("PATH_MIN_SPECS")?params.getValueInt("PATH_MIN_SPECS"):1);
    short pathMinPeaks = (short)(params.paramPresent("PATH_MIN_PEAKS")?params.getValueInt("PATH_MIN_PEAKS"):1);
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;

	bool noSequencing = params.paramPresent("NO_SEQUENCING")?(bool)params.getValueInt("NO_SEQUENCING"):false;
	bool addEndpoints = params.paramPresent("ADD_ENDPOINTS")?(bool)params.getValueInt("ADD_ENDPOINTS"):true;
	const char *wholeABFN = params.getValue("OUTPUT_COMPLETE_ABRUIJN");

	SpecSet specSet;   short loadOk=0;
	if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
	else if(params.paramPresent("INPUT_SPECS_PKLBIN")) loadOk=specSet.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS_PKLBIN"));
	if (loadOk<=0 or specSet.size()==0) return -1;

	SpectrumPairSet aligns;
	SpectrumPairSet tmpAlignsPA;   SpectrumPairSet tmpAlignsASP;  // Quick-fix variables used below
	if (params.paramPresent("INPUT_ALIGNS")) {
		if (aligns.loadFromBinaryFile(alignsFN)==0) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
        else  cout << "Loading aligns complete. Num pairs: " << aligns.size() << "\n";
	}

	vector<SpectrumPeakLabels> labels(0), newLabels;
	if(strlen(labelsFN)!=0) {
		int res = LoadLabels(specSet, labelsFN, labels);
		if(res<=0) { cerr<<"Error reading labels file "<<labelsFN<<".\n"; return -1; }
	}

	if (addEndpoints) {
		// Make sure every spectrum has zero/parentMass-19 nodes
		if(!labels.empty()) for(unsigned int i=0; i<specSet.size(); i++) specSet[i].addZPMpeaks(peakTol,true,&labels[i]);
		else for(unsigned int i=0; i<specSet.size(); i++) specSet[i].addZPMpeaks(peakTol,ionOffset,true);
	}

	// Separate aligns into components
//	vector<SpectrumPairSet> cAligns;    vector<vector<int> > cAligns_idx;  // To allow going back from components to original order
	SetMerger components(specSet.size());
	components.createSets(specSet.size(),2,tmpAlignsASP,aligns);
	components.splitAligns(tmpAlignsASP,aligns);
//	components.splitAligns(aligns,cAligns,cAligns_idx);
	cout<<"Got "<<components.size()<<" component(s).\n";

	vector<MSGraph> spectrumGraphs(specSet.size());   // Used to build an MSGraph for each spectrum
	AAJumps jumps2(2);

	// Split aligns per component and process each component
	char sBuf[1024]; AAJumps jumps(1);
	unsigned int numElemsInSets=components.numElemsInSets();  // Maximum possible number of components (all singletons)
	vector<vector<float> > cStats(numElemsInSets);  for(unsigned int cIdx=0;cIdx<cStats.size();cIdx++) { cStats[cIdx].resize(9); for(unsigned int i=0;i<9;i++) cStats[cIdx][i]=0; }  // Absolute value 9 depends on MSGraph::info_heaviestPath
	Clusters           pathSpectra;   pathSpectra.resize(numElemsInSets);  // Keep the de-novo reconstructed heaviestPath sequences as spectra in a Cluster variable
	vector<list<int> > cSpectra(numElemsInSets);  // Lists of spectrum indices per component
	vector<bool>       specFlipped(specSet.size()); for(unsigned int i=0;i<specFlipped.size();i++) specFlipped[i]=false;

	// Keeps track of which spectrum peaks were matched (dim.3) in each ABruijn vertex (dim.2)
	//   from the de novo sequence (heaviest path) in each component (dim.1)
	vector<vector<list<TwoValues<int> > > > abVertices(numElemsInSets);
	for(unsigned int i=0;i<abVertices.size();i++) abVertices[i].resize(0);

	// Keeps track of which spectrum peaks were matched (dim.3) for all
	//   ABruijn vertices (dim.2) in each component (dim.1)
	vector<vector<list<TwoValues<int> > > > abVerticesAll(numElemsInSets);
	for(unsigned int i=0;i<abVertices.size();i++) abVerticesAll[i].resize(0);

    vector<vector<short> > abCounts(numElemsInSets);   // Records info on number of vertices and edges per ABruijn graph
    for(unsigned int i=0; i<abCounts.size(); i++) { abCounts[i].resize(2); abCounts[i][0]=0; abCounts[i][1]=0; }

//	for(unsigned int cIdx=20; cIdx<21; cIdx++) {
	unsigned int prevNumSpecs;
    for(unsigned int cIdx=0; cIdx<components.size(); cIdx++) {
		cout << "Processing component "<<cIdx<<endl;

		prevNumSpecs = 0;
		while(prevNumSpecs!=components.sets[cIdx].size()) {  // Iterates ABruijn/sequencing/split until all
			                                                 //   remaining spectra match the best ABruijn path
			prevNumSpecs = components.sets[cIdx].size();

			cout << "  - Spectrum indices: "; for(list<int>::iterator iter=components.sets[cIdx].begin();iter!=components.sets[cIdx].end();iter++) cout<<*iter<<" "; cout<<endl; cout.flush();
			cout << "  - Component defined by "<< components.cAlignsPA[cIdx].size() << " pairs...\n"; cout.flush();
			vector<vector<TwoValues<int> > > matches;

			//
			// Choose consensus orientations for pairwise alignments - a spectrum can only be
			//   used as-is or reversed, not both. Also determines set of matched peaks for
			//   the consensus orientations.
			//
			vector<float> modPos;
	//		SplitPairs3(specSet, cAlignsASP[cIdx], components.cAlignsPA[cIdx], peakTol, maxAAjump, penalty_sameVert, penalty_ptm, matches, matchesPA, specFlipped, modPos, false, labelsP);
			SplitPairs(specSet, components.cAlignsPA[cIdx], peakTol, pmTol, maxAAjump, maxModMass, penalty_sameVert, penalty_ptm, matches, specFlipped, modPos, minMatchedPeaks, minEdgesToComponent, false, NULL, &labels);
if(cIdx==62)
	for(unsigned int i=0; i<components.cAlignsPA[cIdx].size(); i++) {
		cerr<<"*** Pair ("<<components.cAlignsPA[cIdx][i].spec1<<","<<components.cAlignsPA[cIdx][i].spec2<<"), shifts "<<components.cAlignsPA[cIdx][i].shift1<<"/"<<components.cAlignsPA[cIdx][i].shift2<<": ";
		for(unsigned int j=0; j<matches[i].size(); j++) cerr<<"("<<matches[i][j][0]<<","<<matches[i][j][1]<<")";
		cerr<<endl;
	}
			SplitASPPA(components.cAlignsPA[cIdx],pmTol,maxModMass,tmpAlignsASP,tmpAlignsPA);

/*			vector<bool> present(specSet.size()); for(unsigned int i=0; i<present.size();i++) present[i]=false;
			for(unsigned int i=0;i<components.cAlignsPA[cIdx].size();i++)
				if(matches[i].size()>0) {
					if(!present[components.cAlignsPA[cIdx][i].spec1]) { cSpectra[cIdx].push_back(components.cAlignsPA[cIdx][i].spec1); present[components.cAlignsPA[cIdx][i].spec1]=true; }
					if(!present[components.cAlignsPA[cIdx][i].spec2]) { cSpectra[cIdx].push_back(components.cAlignsPA[cIdx][i].spec2); present[components.cAlignsPA[cIdx][i].spec2]=true; }
				}
*/
			//
			// Maximize endpoint scores for the matched spectra
			//
			list<int>::iterator sicIter;
			if(addEndpoints)
				for(sicIter=cSpectra[cIdx].begin();sicIter!=cSpectra[cIdx].end();sicIter++)
					specSet[*sicIter].maximizeZPMpeaks(peakTol,true);

	#ifdef DBG_MASAB
			// Add labels to the split pairs graph for graphviz output
			MSGraph g;   g.build(tmpAlignsASP);   g.add(tmpAlignsPA);
			g.vLabels.resize(specSet.size());
			if(specSet.size()==specSet.size())
				for(unsigned int i=0; i<specSet.size(); i++) {
					if(specFlipped[i]) sprintf(sBuf,"v%d_R",i); else sprintf(sBuf,"v%d",i);
					g.vLabels[i] = string((const char *)sBuf);
				}
			else for(unsigned int i=0; i<specSet.size(); i++) {
					sprintf(sBuf,"v%d",i);    g.vLabels[2*i] = string((const char *)sBuf);
					sprintf(sBuf,"v%d_R",i);  g.vLabels[2*i+1] = string((const char *)sBuf);
				 }
			sprintf(sBuf,"split_pairs_graph_%d.txt",cIdx);  g.output_graphviz(sBuf);
	#endif

			//
			// Build spectrum graphs for the spectra in this component
			//
			if(graphType>0)
				for(unsigned int i=0; i<components.cAlignsPA[cIdx].size(); i++) {
					if(spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].numVerts()==0) {
						if(graphType==1) spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].ConnectConsecutive(specSet[components.cAlignsPA[cIdx][i].spec1]);
						else spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].ConnectJumps(specSet[components.cAlignsPA[cIdx][i].spec1],jumps2,peakTol);
//sprintf(sBuf,"spectrum_graph_%d.txt",components.cAlignsPA[cIdx][i].spec1);
//spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].output_graphviz(sBuf);
					}
					if(spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].numVerts()==0) {
						if(graphType==1) spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].ConnectConsecutive(specSet[components.cAlignsPA[cIdx][i].spec2]);
						else spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].ConnectJumps(specSet[components.cAlignsPA[cIdx][i].spec2],jumps2,peakTol);
//sprintf(sBuf,"spectrum_graph_%d.txt",components.cAlignsPA[cIdx][i].spec2);
//spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].output_graphviz(sBuf);
					}
				}

			//
			// Build A-Bruijn graph
			//
			VertexSet vSet(specSet,2048);
			vector<bool> usedSpectra(specSet.size());  for(unsigned int i=0; i<specSet.size(); i++) usedSpectra[i]=false;
			for(unsigned int i=0; i<matches.size(); i++)
				if(matches[i].size()>0) {
					if(graphType==0)  // adds edges for all consecutive matched spectrum peaks
						vSet.addGlues(components.cAlignsPA[cIdx][i].spec1,components.cAlignsPA[cIdx][i].spec2,matches[i]);
					else  // add edges later
						vSet.addGlues(components.cAlignsPA[cIdx][i].spec1,components.cAlignsPA[cIdx][i].spec2,matches[i],&spectrumGraphs);
					usedSpectra[components.cAlignsPA[cIdx][i].spec1] = true;
					usedSpectra[components.cAlignsPA[cIdx][i].spec2] = true;
				}
			list<int> includedSpecs;     includedSpecs.clear();
			for(unsigned int i=0; i<specSet.size(); i++)
				if(usedSpectra[i]) includedSpecs.push_back(i);

			int specCount=0; for(unsigned int i=0;i<usedSpectra.size();i++) if(usedSpectra[i]) specCount++;
			if(specCount==0) { cout<<"  - component is empty, ABruijn graph not built\n"; cout.flush(); continue; }
			cout<<"  - ABruijn graph built on "<<specCount<<" spectra.\n";

			if(specCount>10000) {
				cout << "Spectrum count per component exceeds maximum of 7500 spectra! Skipping component...\n";
				continue;
			}

			//
			// Find/count/list/split composite vertices
			//
			int compositeVertexCount=0;  list<int> compositeSet;
			for(unsigned int i=0; i<vSet.vertices.size(); i++)
				if(vSet.vertices[i].size()>0 and vSet.vertices[i].compositeVertex)
					{ compositeVertexCount++; compositeSet.push_front(i);
cerr<<">>>>> composite vertex "<<i<<": ";
for(list<TwoValues<int> >::iterator iter=vSet.vertices[i].specPeaks.begin(); iter!=vSet.vertices[i].specPeaks.end(); iter++)
	cerr<<"("<<(*iter)[0]<<","<<(*iter)[1]<<")";
cerr<<endl;
					}
			cout<<"  - Abruijn graph contains "<<compositeVertexCount<<" composite vertices: ";
			list<int>::iterator iter=compositeSet.begin();
			for(; iter!=compositeSet.end(); iter++) cout<<(*iter)<<" ";
			if(compositeVertexCount>0) {
				cout<<"-> splitting...";cout.flush();
				vSet.splitComposite(spectrumGraphs,peakTol,&usedSpectra);
				cout<<"done.\n";
			} else cout<<endl;
			cout.flush();

			//
			// Add spectrum graph edges to ABruijn graph (may connect non-aligned peaks, e.g. if
			//   missing in all but one of the aligned spectra)
			//
			if(graphType==1) { vSet.addEdges(spectrumGraphs,&usedSpectra);  vSet.consolidatePaths(); }
			if(graphType==2) { vSet.addEdges(spectrumGraphs,&usedSpectra); }
// Changed 2010/09/14			if(not addEndpoints) vSet.addEndpointEdges(tmpAlignsASP,matches,modPos,jumps2,peakTol);
			if(addEndpoints) vSet.addEndpointEdges(tmpAlignsASP,matches,modPos,jumps2,peakTol);

			//
			//  ABruijn path-finishing procedures should go here.
			//

	#ifdef DBG_MASAB
			vSet.outputGraph(labels);
	#endif

			//
			// Create ABruijn MSGraph from VertexSet for b-ion and y-ion endpoints; choose option with highest-scoring path
			//
			VertexSet copy(vSet), *vSetP;
			MSGraph *abg, abgB, abgY, path;
			vector<int> vSet_indexB, vSet_indexY, *vSet_indexP; // Correspondences between vertex indices in the ABruijn graph and simplified graph
			vector<int> pathVertsIdxB, pathVertsIdxY, *pathVertsIdx;  // Indices of the vertices in the heaviest path
			Spectrum specWithBep, specWithYep;
			float scoreWithBep, scoreWithYep;

			if(addEndpoints) vSet.removeEndpoints(false,peakTol);  // Remove y-mass endpoints added above
			//
			//  ABruijn endpoint-finishing procedures should go here.
			//
			vSet.buildGraph(abgB,jumps,peakTol,vSet_indexB,labels,edgeScoreType);
//char filename[2048];   sprintf(filename,"component_%d.txt",cIdx+1);
//abgB.output_graphviz(filename);
			if(noSequencing) {
				pathVertsIdxB.resize(0);
				if(wholeABFN) vSet.getMatchedPeaks(pathVertsIdxB,abVerticesAll[cIdx]);
				continue;
			}
			scoreWithBep = abgB.heaviestPath(path,false,&specWithBep,&pathVertsIdxB);

			if(addEndpoints) copy.removeEndpoints(true,peakTol);  // Remove b-mass endpoints added above
			//
			//  ABruijn endpoint-finishing procedures should go here.
			//
			copy.buildGraph(abgY,jumps,peakTol,vSet_indexY,labels,edgeScoreType);
//sprintf(filename,"component_%d_y.txt",cIdx+1);
//abgY.output_graphviz(filename);
			scoreWithYep = abgY.heaviestPath(path,false,&specWithYep,&pathVertsIdxY);

//			scoreWithBep = max(scoreWithBep,scoreWithYep)+1;  // Force B endpoints

			if(scoreWithBep>=scoreWithYep) { abg=&abgB; pathSpectra.consensus[cIdx]=specWithBep; vSet_indexP = &vSet_indexB; vSetP = &vSet; pathVertsIdx=&pathVertsIdxB; cerr<<"  - Selected B endpoints.\n";}
			  else { abg=&abgY; pathSpectra.consensus[cIdx]=specWithYep; vSet_indexP = &vSet_indexY; vSetP = &copy; pathVertsIdx=&pathVertsIdxY; cerr<<"  - Selected Y endpoints.\n";}
			for(unsigned int vIdx=0; vIdx<pathVertsIdx->size(); vIdx++) (*pathVertsIdx)[vIdx] = (*vSet_indexP)[(*pathVertsIdx)[vIdx]];  // Convert simplified graph vertex indices to ABruijn vertex indices.

	#ifdef DBG_MASAB
	//		sprintf(sBuf,"graph_ma_%d.txt",cIdx);
	//		vSetP->output_graphviz_ma(sBuf, *pathVertsIdx);
	#endif
			vSetP->getMatchedPeaks(*pathVertsIdx,abVertices[cIdx]);
			pathSpectra.set_endpoints(cIdx,specSet,abVertices[cIdx],peakTol,pmTol);

			if(wholeABFN) {
				pathVertsIdx->resize(0);
				vSetP->getMatchedPeaks(*pathVertsIdx,abVerticesAll[cIdx]);
			}

			abg->info_heaviestPath(cStats[cIdx]);   abCounts[cIdx][0]=abg->numVertices();   abCounts[cIdx][1]=abg->numEdges();
	//		cout<<"  - Heaviest path stats: ["; for(unsigned int i=0;i<cStats[cIdx].size();i++) {cout<<cStats[cIdx][i]; if(i<cStats[cIdx].size()-1) cout<<", "; } cout<<"]\n";

	#ifdef DBG_MASAB
			sprintf(sBuf,"graph_%d.txt",cIdx); abg->output_graphviz(sBuf);
	#endif

			//
			// Find spectra with not enough ABruijn vertices on the heaviest path (if any)
			//   and remove them from the current component. Unused spectra may define new (leftover)
			//   components if connected by at least one edge (ie, at least two spectra in a
			//   connected component).
			//
			list<int> usedSpecs;     usedSpecs.clear();
			list<int> unusedSpecs;   unusedSpecs.clear();
			vector<short> numPathPeaks(specSet.size());   for(unsigned int i=0; i<specSet.size(); i++) numPathPeaks[i]=0;
			for(unsigned int i=0; i<abVertices[cIdx].size(); i++) {
//cerr<<"AB vertex "<<i<<": ";
				list<TwoValues<int> >::iterator vNext; // Used to remove remaining composite vertices
				for(list<TwoValues<int> >::iterator vIter=abVertices[cIdx][i].begin(); vIter!=abVertices[cIdx][i].end(); ) {
//cerr<<"("<<(*vIter)[0]<<","<<(*vIter)[1]<<")";
					vNext=vIter;  if(vNext!=abVertices[cIdx][i].end()) vNext++;
					if(vNext!=abVertices[cIdx][i].end() and (*vIter)[0]==(*vNext)[0]) {
						cerr<<"  - ERROR: inconsistent contig vertex containing ("<<(*vIter)[0]<<","<<(*vIter)[1]<<") and ("<<(*vNext)[0]<<","<<(*vNext)[1]<<")"; cerr.flush();
						vIter=abVertices[cIdx][i].erase(vIter);  // Remove inconsistent spectrum/peaks
						vIter=abVertices[cIdx][i].erase(vIter);
					} else {
						if(++numPathPeaks[(*vIter)[0]]>=pathMinPeaks) usedSpecs.push_back((*vIter)[0]);
						vIter++;
					}
				}
//cerr<<endl;
			}
			usedSpecs.sort();     usedSpecs.unique();

			// Prune ABruijn vertices: remove peaks from spectra without enough matches to the consensus path
			unsigned int maxPeaksPerVertex=0;
			for(unsigned int i=0; i<abVertices[cIdx].size(); i++) {
				for(list<TwoValues<int> >::iterator vIter=abVertices[cIdx][i].begin(); vIter!=abVertices[cIdx][i].end(); )
					if(numPathPeaks[(*vIter)[0]]>=pathMinPeaks) vIter++;
					else { unusedSpecs.push_back((*vIter)[0]); vIter=abVertices[cIdx][i].erase(vIter); }
				if(abVertices[cIdx][i].size()>maxPeaksPerVertex) maxPeaksPerVertex=abVertices[cIdx][i].size();
				if(abVertices[cIdx][i].size()>usedSpecs.size()) {
					cerr<<"  - ERROR: inconsistent contig vertex: "; cerr.flush();
					for(list<TwoValues<int> >::iterator vIter=abVertices[cIdx][i].begin(); vIter!=abVertices[cIdx][i].end(); vIter++)
						{ cerr<<"("<<(*vIter)[0]<<","<<(*vIter)[1]<<")"; cerr.flush(); }
					cerr<<"\n";
				}
			}
			unusedSpecs.sort();   unusedSpecs.unique();

			// Output complete set of vertices
			if(wholeABFN) {
				for(unsigned int i=0; i<abVerticesAll[cIdx].size(); i++)
					for(list<TwoValues<int> >::iterator vIter=abVerticesAll[cIdx][i].begin(); vIter!=abVerticesAll[cIdx][i].end(); )
						if(numPathPeaks[(*vIter)[0]]>=pathMinPeaks) vIter++; else vIter=abVerticesAll[cIdx][i].erase(vIter);
			}

			if(usedSpecs.size()>=pathMinSpecs)  {
				if (usedSpecs.size()<components.sets[cIdx].size()-1) {  // Whenever there are at least 2 unused spectra
					cout<<"  - Keeping "<<usedSpecs.size()<<" spectra; number of components: "<<components.size()<<" -> ";
					components.splitSet(cIdx,usedSpecs);
					cout<<components.size()<<"\n"; cout.flush();
				} else if (usedSpecs.size()==components.sets[cIdx].size()-1) components.removeElements(cIdx,unusedSpecs);
				if(maxPeaksPerVertex>usedSpecs.size()) {
					cerr<<"  - ERROR: inconsistent contig containing a vertex with "<<maxPeaksPerVertex<<" spectrum peaks but only "<<usedSpecs.size()<<" used spectra (contig deleted)\n";
					pathSpectra.consensus[cIdx].resize(0);
					abVertices[cIdx].resize(0);
				}
			} else {
				if (includedSpecs.size()>0) {  // Whenever there is at least 1 included spectrum
					cout<<"  - Keeping "<<components.sets[cIdx].size()-includedSpecs.size()<<" spectra; number of components: "<<components.size()<<" -> ";
					components.splitSet(cIdx,includedSpecs);
					cout<<components.size()<<"\n"; cout.flush();
				}
				cout<<"  - Only "<<usedSpecs.size()<<" spectra with at least "<<pathMinPeaks<<" peaks in the consensus path - component deleted.\n";
				pathSpectra.consensus[cIdx].resize(0);  // Remove de novo sequences for poor contigs (not enough spectra with enough matched peaks)
			}

			cout.flush();   cerr.flush();

			if(components.cAlignsPA[cIdx].size()==0 and components.cAlignsASP[cIdx].size()==0) { // No pairs left
				prevNumSpecs = components.sets[cIdx].size();
				cout<<"  --> No pairs left for second iteration - keeping ABruijn path\n";
			}
		} // while (prevNumSpecs!=components.sets[cIdx].size())
	}
	// Resize down to the final number of resulting connected components
	abCounts.resize(components.size());
	cStats.resize(components.size());
	pathSpectra.resize(components.size());
	cSpectra.resize(components.size());
	abVertices.resize(components.size());
	abVerticesAll.resize(components.size());

	components.saveas_binListArray("components.bla");
	Save_binArray("component_stats.bna",cStats);
	Save_binListArray<int,list<int>,list<int>::iterator>("component_spectra.bla",cSpectra);
	Save_abinfo("component_info.bin",specSet, components.sets, specFlipped, abVertices);
	if(params.paramPresent("OUTPUT_CLUSTERS")) pathSpectra.Save(params.getValue("OUTPUT_CLUSTERS"));
	else pathSpectra.Save("path_spectra_as_cluster.txt");
	pathSpectra.consensus.SaveSpecSet_pklbin(params.getValue("OUTPUT_SPECS"));
	Save_binArray("abCounts.bin", abCounts);
	if(params.paramPresent("OUTPUT_COMPLETE_ABRUIJN")) Save_abinfo(wholeABFN,specSet, components.sets, specFlipped, abVerticesAll);

	return(0);
}
