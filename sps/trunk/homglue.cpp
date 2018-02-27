#include "inputParams.h"
#include "alignment_modmut.h"
#include "spectral_pairs.h"
#include "abruijn.h"
#include "graph.h"
#include "SetMerger.h"
#include "SpectrumPairSet.h"
#include "SpectralPairs.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;
using namespace specnets;

/*
  homglue - Glues spectra into ABruijn graphs if matched to the overlapping
  protein locations. If ClustalW protein sequence alignments are available, also
  transfers spectrum/protein-K matches to a reference protein-R if there is a
  ClustalW alignment between protein-K and protein-R.

  INPUT:
    INPUT_SPECS     - Input deconvolved spectra (mostly b or mostly y but not both)
    INPUT_FASTA     - Protein sequence database
    INPUT_MATCHED_PROTS     - Binary file with protein matches per spectrum: protein index and number of mods
    INPUT_MATCHED_PEAKS_IDX - Binary pklbin file with sets of matched masses between each spectrum/protein
    SPEC_TYPE_MSMS  - Spectrum type: 0 (PRM), 1 (MS/MS), default is PRM
	TOLERANCE_PEAK  - Peak mass tolerance (in Daltons, default 0.5 Da)
    TOLERANCE_PM    - Parent mass tolerance (in Daltons, default 1 Da)
	MAX_AA_JUMP     - Largest mass jump in the ABruijn graphs (default is 2)
	MIN_CONTIG_SET  - Minimum number of spectra per ABruijn graph (default is 1)
    EDGE_SCORE_TYPE -
	            EST_EDGE_MULT = 0,       // Edge multiplicity
	            EST_EDGE_SCORES = 1,     // Edge scores: each edge collects its score from the
	                                                     destination peak in the original spectrum (DEFAULT)
	            EST_ABVERTEX_SCORES = 2; // Vertex scores: each edge collects its score from the destination A-Bruijn vertex
	                                     //   Edge multiplicity is added to vertex scores to help distinguish between edges
    GRAPH_TYPE - ABruijn graph edges are derived by gluing
                0 - path graphs from sets of consecutive matched peaks
                1 - path graphs (each spectrum is represented as a path through
                      _all_ of its peaks) - DEFAULT option
                2 - spectrum graphs (one per spectrum)
	PATH_MIN_PEAKS - Minimum number of peaks for a spectrum to be part of a protein contig sequence
	PATH_MIN_SPECS - Minimum number of spectra (observing PATH_MIN_PEAKS) for a valid protein contig sequence

  OUTPUT:
	OUTPUT_SPECS - Highest-scoring paths in the resulting ABruijn graphs
	OUTPUT_CSV   - Matches between the ABruijn highest-scoring paths and the protein sequences
*/
int main(int argc, char **argv){

	// Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("homglue.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "homglue.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs;   paramStrs.resize(5);
	paramStrs[0] = "OUTPUT_SPECS";
	paramStrs[1] = "OUTPUT_CSV";
	paramStrs[2] = "INPUT_MATCHED_PROTS";
	paramStrs[3] = "INPUT_MATCHED_PEAKS_IDX";
	paramStrs[4] = "INPUT_FASTA";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"homglue.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: OUTPUT_SPECS, OUTPUT_CSV, INPUT_FASTA, INPUT_MATCHED_PROTS, INPUT_MATCHED_PEAKS_IDX\n";
		return -1;
	}

	int graphType = params.paramPresent("GRAPH_TYPE")?params.getValueInt("GRAPH_TYPE"):1;
    short edgeScoreType = (short)(params.paramPresent("EDGE_SCORE_TYPE")?params.getValueInt("EDGE_SCORE_TYPE"):1);
	int minContigSet = params.paramPresent("MIN_CONTIG_SET")?params.getValueInt("MIN_CONTIG_SET"):1;
	int maxAAjump = params.paramPresent("MAX_AA_JUMP")?params.getValueInt("MAX_AA_JUMP"):2;
	float maxModMass = params.paramPresent("MAX_MOD_MASS")?(float) params.getValueDouble("MAX_MOD_MASS"):100.0;
	float peakTol = params.getValue("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.getValue("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;
	bool specNets = params.paramPresent("INPUT_SPECNETS_ALIGNS") and params.paramPresent("INPUT_SPECNETS_MATCHED");
	float penalty_ptm = params.getValue("PENALTY_PTM")?(float)params.getValueDouble("PENALTY_PTM"):0;
	float penalty_sameVert = params.getValue("PENALTY_SAME_VERTEX")?(float)params.getValueDouble("PENALTY_SAME_VERTEX"):-1000000;
    unsigned int cIdx, pairIdx, peakIdx, specIdx;

	SpecSet specSet;   short loadOk=0;
	if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
	else if(params.paramPresent("INPUT_SPECS_PKLBIN")) loadOk=specSet.loadPklBin(params.getValue("INPUT_SPECS_PKLBIN"));
	if (loadOk<=0 or specSet.size()==0) return -1;

	//
	// Process Contig <-> Protein Alignments
	//
    DB_fasta db;
    int count = db.Load(params.getValue("INPUT_FASTA"));  if(count<=0) { cerr<<"Error reading "<<params.getValue("INPUT_FASTA")<<"!\n"; return(-1); }
    for(unsigned int specIdx=0;specIdx<count;specIdx++) db.getMassesSpec(specIdx);

    AMM_match_set matchSet(db,specSet.size());  // Matches between SPS contigs and either homolog(loaded) or reference(transfered) proteins
    if(not matchSet.LoadMatches(params.getValue("INPUT_MATCHED_PROTS"),params.getValue("INPUT_MATCHED_PEAKS_IDX")))
      return -1;

	// Load clustalw protein <-> protein alignments (if any)
	if(params.paramPresent("INPUT_HOMOLOGIES")) {
		vector<TwoValues<int> > matchedProts;
		vector<vector<TwoValues<int> > > matchedAAs;
		if(not Load_clustalw_multiple(params.getValue("INPUT_HOMOLOGIES"), db, matchedProts, matchedAAs, 1)) return -1;

		// Use clustalw alignments to transitively import contig matches to the reference protein
		for(unsigned int matchIdx=0; matchIdx<matchedProts.size(); matchIdx++)
			matchSet.MergeIntoReference(matchedProts[matchIdx][0],matchedProts[matchIdx][1],matchedAAs[matchIdx]);
	}

	// Output modified spectrum/protein matches (after transferring some matches to the reference protein)
	if(params.paramPresent("OUTPUT_MATCHES_REF")) {
		matchSet.matchedPeaks.SaveSpecSet_pklbin((string(params.getValue("OUTPUT_MATCHES_REF"))+string("_midx.pklbin")).c_str());
	    Save_binArray((string(params.getValue("OUTPUT_MATCHES_REF"))+string("_mp.bin")).c_str(), matchSet.specProtMatches);
	}

	vector<bool> specFlipped(specSet.size());
	for(unsigned int specIdx=0; specIdx<specFlipped.size(); specIdx++)
		specFlipped[specIdx]=matchSet.specProtMatches[specIdx][0]==-1?false:(bool)matchSet.specProtMatches[specIdx][2];

	//
	// Sequence resulting contigs
	//
	SpectrumPairSet pairsASP(0);
	SpectrumPairSet pairsPA(0);
    vector<vector<TwoValues<int> > > pairMatchedPeaks(0),pairMatchedPeaksASP(0);
    if(not params.paramPresent("SKIP_GLUES") or params.getValueInt("SKIP_GLUES")!=1)
    	matchSet.GenerateGlues(pairsPA,pairMatchedPeaks);

	// Load pairwise alignments from spectral networks (if available)
	if(specNets) {
//		if (!Load_results_bin(params.getValue("INPUT_SPECNETS_ALIGNS"), pairsASP)) { cerr << "Error reading "<<params.getValue("INPUT_SPECNETS_ALIGNS")<<"!\n"; return -1; }
//		cout << "Loading spectral networks pairs complete (" << pairsASP.size() << " pairs)\n";
		if (!pairsASP.loadFromBinaryFile(params.getValue("INPUT_SPECNETS_ALIGNS"))) { cerr << "Error reading "<<params.getValue("INPUT_SPECNETS_ALIGNS")<<"!\n"; return -1; }
		cout << "Loading spectral networks pairs complete (" << pairsASP.size() << " pairs)\n";

		SpecSet snMatched;
	    if(snMatched.loadPklBin(params.getValue("INPUT_SPECNETS_MATCHED"))<=0)
			{ cerr<<"Error loading "<<params.getValue("INPUT_SPECNETS_MATCHED")<<"\n"; return -1; }
	    cout << "Loading spectral network matched peaks complete (" << snMatched.size() << " pairs)\n";
		pairMatchedPeaksASP.resize(snMatched.size());
		for(pairIdx=0; pairIdx<snMatched.size(); pairIdx++) {
			pairMatchedPeaksASP[pairIdx].resize(snMatched[pairIdx].size());
			for(peakIdx=0; peakIdx<snMatched[pairIdx].size(); peakIdx++) {
				pairMatchedPeaksASP[pairIdx][peakIdx].set((int)snMatched[pairIdx][peakIdx][0],(int)snMatched[pairIdx][peakIdx][1]);
			}
		}
	}

	if(params.paramPresent("INPUT_SPECNETS_ALIGNS_ALL")) {
		SpectrumPairSet pairsASP_all(0);
//		if (!Load_results_bin(params.getValue("INPUT_SPECNETS_ALIGNS_ALL"), pairsASP_all)) { cerr << "Error reading "<<params.getValue("INPUT_SPECNETS_ALIGNS_ALL")<<"!\n"; return -1; }
//		cout << "Loading full set of spectral networks pairs complete (" << pairsASP_all.size() << " pairs)\n";
		if (!pairsASP_all.loadFromBinaryFile(params.getValue("INPUT_SPECNETS_ALIGNS_ALL"))) { cerr << "Error reading "<<params.getValue("INPUT_SPECNETS_ALIGNS_ALL")<<"!\n"; return -1; }
		cout << "Loading full set of spectral networks pairs complete (" << pairsASP_all.size() << " pairs)\n";

		// Find edges between unidentified spectra
		SetMerger components_annot(specSet.size());
		components_annot.createSets(specSet.size(),minContigSet,pairsASP,pairsPA);
		components_annot.compressSetIndices();
		unsigned int aspCount, szKept=0;
		for(pairIdx=0; pairIdx<pairsASP_all.size(); pairIdx++) {
			if(components_annot.membership[pairsASP_all[pairIdx].spec1]==components_annot.membership[pairsASP_all[pairIdx].spec2]) continue;  // Already in the same component
			if(matchSet.matchedPeaks[pairsASP_all[pairIdx].spec1].size()>0 or matchSet.matchedPeaks[pairsASP_all[pairIdx].spec2].size()>0) continue;  // At least one spectrum is annotated
			pairsASP_all[szKept++]=pairsASP_all[pairIdx];
		}
		pairsASP_all.resize(szKept);

		// Add edges between unidentified spectra to global set of spectral networks edges
		aspCount = pairsASP.size();
		pairsASP.resize(aspCount+szKept);   pairMatchedPeaksASP.resize(aspCount+szKept);
		for(pairIdx=0; pairIdx<szKept; pairIdx++) {
			pairsASP[aspCount+pairIdx] = pairsASP_all[pairIdx];
			pairMatchedPeaksASP[aspCount+pairIdx].resize(0);
		}
	}

#ifdef DBG_HOMGLUE
    cerr<<"Got "<<pairsPA.size()<<" pairs of glued spectra\n";
    ofstream dbg_output;   dbg_output.open("dbg_glues.txt", ios::binary);
    for(unsigned int i=0; i<pairsPA.size(); i++) {
    	dbg_output<<"Pair "<<pairsPA[i].spec1<<" / "<<pairsPA[i].spec2<<":\n";
    	for(unsigned int j=0; j<pairMatchedPeaks[i].size(); j++)
    		dbg_output<<" peaks ("<<pairMatchedPeaks[i][j][0]<<","<<pairMatchedPeaks[i][j][1]<<"), masses ("<<specSet[pairsPA[i].spec1][pairMatchedPeaks[i][j][0]][0]<<","<<specSet[pairsPA[i].spec2][pairMatchedPeaks[i][j][1]][0]<<")\n";
    	dbg_output<<endl;

    }
    dbg_output.close();
#endif

    // Find connected components (aligned spectra)
	SetMerger components(specSet.size());
	components.createSets(specSet.size(),minContigSet,pairsASP,pairsPA);
	if(specNets) { // Keep only annotated components (including singletons) and unidentified networks
		for(unsigned int setIdx=0; setIdx<components.sets.size(); setIdx++)
			if(components.sets[setIdx].size()==1 and matchSet.matchedPeaks[components.sets[setIdx].front()].size()==0)
				components.removeSet(setIdx,false);
		components.compressSetIndices();
	}
	components.splitAligns(pairsASP,pairsPA);
	cout<<"Got "<<components.size()<<" component(s), minimum size was "<<minContigSet<<".\n";
	cout.flush();

	// Sequence ABruijn graphs
	AAJumps jumps(maxAAjump);   char sBuf[1024];
	unsigned int numElemsInSets=components.numElemsInSets();  // Maximum possible number of components (all singletons)
	vector<MSGraph> spectrumGraphs(specSet.size());
	Clusters           pathSpectra;   pathSpectra.resize(numElemsInSets);  // Keep the de-novo reconstructed heaviestPath sequences as spectra in a Cluster variable
	vector<AMM_match_spec> matches(numElemsInSets);   // Matches between sequenced CSPS contigs and protein sequences
	AMM_match_set csps_matches(db,numElemsInSets);    // Container for all the matches between csps contigs and the protein sequence database

	vector<vector<float> > cStats(numElemsInSets);  for(unsigned int cIdx=0;cIdx<cStats.size();cIdx++) { cStats[cIdx].resize(9); for(unsigned int i=0;i<9;i++) cStats[cIdx][i]=0; }  // Absolute value 9 depends on MSGraph::info_heaviestPath
	vector<list<int> > cSpectra(numElemsInSets);  // Lists of spectrum indices per component
	vector<vector<list<TwoValues<int> > > > abVertices(numElemsInSets); for(unsigned int i=0;i<abVertices.size();i++) abVertices[i].resize(0);   // Keeps track of which spectrum peaks were matched (dim.3) in each ABruijn vertex (dim.2) in each component (dim.1)
    vector<vector<short> > abCounts(numElemsInSets);   // Records info on number of vertices and edges per ABruijn graph
    for(unsigned int i=0; i<abCounts.size(); i++) { abCounts[i].resize(2); abCounts[i][0]=0; abCounts[i][1]=0; }

	for(cIdx=0; cIdx<components.size(); cIdx++) {
		AMM_match curMatch;        // Resulting path-protein match
		TwoValues<float> curPairMasses;    TwoValues<int> curPairIdx;
		if(components.sets[cIdx].size()==1) {  // No gluing/resequencing necessary, just reproduce input contig/protein match
			matches[cIdx].init(1,0);           specIdx = components.sets[cIdx].front();
			curMatch.orientationPRMs=0;        curMatch.modSize=0;    curMatch.matchScore=0;
			curMatch.matchedIndices.clear();   curMatch.matchedMasses.clear();
			curMatch.proteinIdx = matchSet.specProtMatches[specIdx][0];

			for(peakIdx=0; peakIdx<matchSet.matchedPeaks[specIdx].size(); peakIdx++) {
				curPairIdx.set((int)round(matchSet.matchedPeaks[specIdx][peakIdx][0]),(int)round(matchSet.matchedPeaks[specIdx][peakIdx][1]));
				curMatch.matchedIndices.push_front(curPairIdx);
				curPairMasses.set(specSet[specIdx][curPairIdx[0]][0],db.masses[curMatch.proteinIdx][curPairIdx[1]][0]);
				curMatch.matchedMasses.push_front(curPairMasses);
				curMatch.matchScore += specSet[specIdx][curPairIdx[0]][1];
			}

			if(not curMatch.matchedIndices.empty()) {
				curMatch.aaStart = curMatch.matchedIndices.back()[1];
				curMatch.aaEnd   = curMatch.matchedIndices.front()[1];
			}
			matches[cIdx].addMatch(0,curMatch);
			csps_matches.set(cIdx,curMatch);
			pathSpectra.consensus[cIdx] = specSet[cIdx];
			continue;
		}

		cout << "Processing component "<<cIdx<<" ("<<components.sets[cIdx].size()<<" elements)\n"; cout.flush();
		cout << "  - spectrum indices: "; for(list<int>::iterator iter=components.sets[cIdx].begin(); iter!=components.sets[cIdx].end();iter++) cout<<(*iter)<<" "; cout<<"\n"; cout.flush();

		SpectrumPairSet &curAligns    = components.cAlignsPA[cIdx];
		vector<int> &curAlignsIdx     = components.cAlignsPA_idx[cIdx];
		SpectrumPairSet &curAlignsASP = components.cAlignsASP[cIdx];
		vector<int> &curAlignsASPIdx  = components.cAlignsASP_idx[cIdx];

		//
		// Build spectrum graphs for the spectra in this component
		//
		if(graphType>0) {
			// alignsPA
			for(unsigned int i=0; i<curAligns.size(); i++) {
//if(curAligns[i].spec1>spectrumGraphs.size() or curAligns[i].spec2>spectrumGraphs.size()) {
//	cerr<<"ERROR: pair ("<<curAligns[i].spec1<<","<<curAligns[i].spec2<<") with only "<<spectrumGraphs.size()<<" spectra!\n"; exit(-1);
//}
				if(spectrumGraphs[curAligns[i].spec1].numVerts()==0) {
					if(graphType==1) spectrumGraphs[curAligns[i].spec1].ConnectConsecutive(specSet[curAligns[i].spec1]);
					else spectrumGraphs[curAligns[i].spec1].ConnectJumps(specSet[curAligns[i].spec1],jumps,peakTol);
//					sprintf(sBuf,"spectrum_graph_%d.txt",curAligns[i].spec1); spectrumGraphs[curAligns[i].spec1].output_graphviz(sBuf);
				}
				if(spectrumGraphs[curAligns[i].spec2].numVerts()==0) {
					if(graphType==1) spectrumGraphs[curAligns[i].spec2].ConnectConsecutive(specSet[curAligns[i].spec2]);
					else spectrumGraphs[curAligns[i].spec2].ConnectJumps(specSet[curAligns[i].spec2],jumps,peakTol);
//					sprintf(sBuf,"spectrum_graph_%d.txt",curAligns[i].spec2); spectrumGraphs[curAligns[i].spec2].output_graphviz(sBuf);
				}
			}
			// alignsASP
			for(unsigned int i=0; i<curAlignsASP.size(); i++) {
//if(curAlignsASP[i].spec1>spectrumGraphs.size() or curAlignsASP[i].spec2>spectrumGraphs.size()) {
//	cerr<<"ERROR: pair ("<<curAlignsASP[i].spec1<<","<<curAlignsASP[i].spec2<<") with only "<<spectrumGraphs.size()<<" spectra!\n"; exit(-1);
//}
				if(spectrumGraphs[curAlignsASP[i].spec1].numVerts()==0) {
					if(graphType==1) spectrumGraphs[curAlignsASP[i].spec1].ConnectConsecutive(specSet[curAlignsASP[i].spec1]);
					else spectrumGraphs[curAlignsASP[i].spec1].ConnectJumps(specSet[curAlignsASP[i].spec1],jumps,peakTol);
//					sprintf(sBuf,"spectrum_graph_%d.txt",curAlignsASP[i].spec1); spectrumGraphs[curAlignsASP[i].spec1].output_graphviz(sBuf);
				}
				if(spectrumGraphs[curAlignsASP[i].spec2].numVerts()==0) {
					if(graphType==1) spectrumGraphs[curAlignsASP[i].spec2].ConnectConsecutive(specSet[curAlignsASP[i].spec2]);
					else spectrumGraphs[curAlignsASP[i].spec2].ConnectJumps(specSet[curAlignsASP[i].spec2],jumps,peakTol);
//					sprintf(sBuf,"spectrum_graph_%d.txt",curAlignsASP[i].spec2); spectrumGraphs[curAlignsASP[i].spec2].output_graphviz(sBuf);
				}
			}
		}

		//
		// Build A-Bruijn graph
		//
		VertexSet vSet(specSet,8192);
		vector<bool> usedSpectra(specSet.size());  for(unsigned int i=0; i<specSet.size(); i++) usedSpectra[i]=false;
		for(pairIdx=0; pairIdx<curAlignsIdx.size(); pairIdx++) {
			if(graphType==0) vSet.addGlues(curAligns[pairIdx].spec1,curAligns[pairIdx].spec2,pairMatchedPeaks[curAlignsIdx[pairIdx]]);
			else vSet.addGlues(curAligns[pairIdx].spec1,curAligns[pairIdx].spec2,pairMatchedPeaks[curAlignsIdx[pairIdx]],&spectrumGraphs);
			usedSpectra[curAligns[pairIdx].spec1] = true;
			usedSpectra[curAligns[pairIdx].spec2] = true;
		}
		if(curAlignsASP.size()>0 and pairMatchedPeaksASP[curAlignsASPIdx[0]].empty()) {  // Need to compute matched peaks for unidentified spectral network
			vector<vector<TwoValues<int> > > matches;    vector<float> modPos;
			SplitPairs(specSet, curAlignsASP, peakTol, pmTol, maxAAjump, maxModMass, penalty_sameVert, penalty_ptm, matches, specFlipped, modPos);

		  for(pairIdx=0; pairIdx<curAlignsASP.size(); pairIdx++) {
				if(graphType==0) vSet.addGlues(curAlignsASP[pairIdx].spec1,curAlignsASP[pairIdx].spec2,matches[pairIdx]);
				else vSet.addGlues(curAlignsASP[pairIdx].spec1,curAlignsASP[pairIdx].spec2,matches[pairIdx],&spectrumGraphs);
				usedSpectra[curAlignsASP[pairIdx].spec1] = true;
				usedSpectra[curAlignsASP[pairIdx].spec2] = true;
			}
		} else {
			for(pairIdx=0; pairIdx<curAlignsASP.size(); pairIdx++) {
				if(graphType==0) vSet.addGlues(curAlignsASP[pairIdx].spec1,curAlignsASP[pairIdx].spec2,pairMatchedPeaksASP[curAlignsASPIdx[pairIdx]]);
				else vSet.addGlues(curAlignsASP[pairIdx].spec1,curAlignsASP[pairIdx].spec2,pairMatchedPeaksASP[curAlignsASPIdx[pairIdx]],&spectrumGraphs);
				usedSpectra[curAlignsASP[pairIdx].spec1] = true;
				usedSpectra[curAlignsASP[pairIdx].spec2] = true;
			}
		}
		int specCount=0; for(unsigned int i=0;i<usedSpectra.size();i++) if(usedSpectra[i]) specCount++;
		cout<<"  - ABruijn graph built on "<<specCount<<" spectra\n";

		// Variables used when renaming vertices in the ABruijn PRM graph and outputting it to graphviz
		MSGraph abg;
		MSGraph path;
		vector<int> vSet_index;    // Correspondences between vertex indices in the ABruijn graph and simplified graph
		vector<int> pathVertsIdx;  // Indices of the vertices in the heaviest path

		//
		// Find/count/list/split composite vertices
		//
		int compositeVertexCount=0;  list<int> compositeSet;
		for(unsigned int i=0; i<vSet.vertices.size(); i++)
			if(vSet.vertices[i].size()>0 and vSet.vertices[i].compositeVertex)
				{ compositeVertexCount++; compositeSet.push_front(i); }
		cout<<"  - Abruijn graph contains "<<compositeVertexCount<<" composite vertices: ";
		list<int>::iterator iter=compositeSet.begin();
		for(; iter!=compositeSet.end(); iter++) cout<<(*iter)<<" ";
		if(compositeVertexCount>0)
			{ cout<<"-> splitting...";cout.flush(); vSet.splitComposite(spectrumGraphs,peakTol,&usedSpectra); cout<<"done.\n";}
		else
			cout<<endl;
		cout.flush();

		//
		// Add spectrum graph edges to ABruijn graph
		//
		if(graphType==1) { vSet.addEdges(spectrumGraphs,&usedSpectra);  vSet.consolidatePaths(); }
		if(graphType==2) { vSet.addEdges(spectrumGraphs,&usedSpectra); }
//		vSet.addEndpointEdges(components.cAlignsASP[cIdx],matches,modPos,jumps,peakTol);

		vSet.removeEndpoints(false,peakTol);
		vector<SpectrumPeakLabels> foo_labels;
		vSet.buildGraph(abg,jumps,peakTol,vSet_index,foo_labels, edgeScoreType);
//char filename[2048];   sprintf(filename,"graphs/component_%d.txt",cIdx+1);
//abg.output_graphviz(filename);
		curMatch.matchScore = abg.heaviestPath(path,false,&pathSpectra.consensus[cIdx],&pathVertsIdx);
		for(unsigned int vIdx=0; vIdx<pathVertsIdx.size(); vIdx++) pathVertsIdx[vIdx] = vSet_index[pathVertsIdx[vIdx]];  // Convert simplified graph vertex indices to ABruijn vertex indices.
		vSet.getMatchedPeaks(pathVertsIdx,abVertices[cIdx]);

		//
		// Get correspondence between protein(s) and heaviest path
		//
		matches[cIdx].init(1,0);
		curMatch.orientationPRMs=0;        curMatch.modSize=0;
		curMatch.matchedIndices.clear();   curMatch.matchedMasses.clear();
cerr<<" -->> abVertices["<<cIdx<<"].size() = "<<abVertices[cIdx].size()<<endl;
cerr<<" -->> matches["<<cIdx<<"] = "<<matches[cIdx].matches[0].size()<<" hits\n";
		if(abVertices[cIdx].size()>0) {
/* Commented on Feb 24, 2010 and replaced with code compatible with Spectral Networks components
			curMatch.proteinIdx = matchSet.specProtMatches[abVertices[cIdx][0].front()[0]][0];  // Just take the first matched protein - most likely all spectra in the ABruijn component were matched to the same protein (except for spectral networks inputs)
			for(unsigned int vIdx=0; vIdx<abVertices[cIdx].size(); vIdx++) {
				for(list<TwoValues<int> >::iterator peakIter = abVertices[cIdx][vIdx].begin(); peakIter!=abVertices[cIdx][vIdx].end(); peakIter++) {
					for(peakIdx=0; peakIdx<matchSet.matchedPeaks[(*peakIter)[0]].size(); peakIdx++)
						if(matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][0]==(*peakIter)[1]) break;   // Found a peak in an ABruijn vertex that was initially matched to the protein sequence
                                                                                                       // Possible problem: different peaks in the same ABruijn vertex may match different protein locations - just taking the first match may not be the best option

					if(peakIdx<matchSet.matchedPeaks[(*peakIter)[0]].size() and matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][0]==(*peakIter)[1]) {
						curPairIdx.set(vIdx,(int)round(matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][1]));
						curMatch.matchedIndices.push_front(curPairIdx);
						curPairMasses.set(pathSpectra.consensus[cIdx][vIdx][0],db.masses[curMatch.proteinIdx][curPairIdx[1]][0]);
						curMatch.matchedMasses.push_front(curPairMasses);
						break;
					}
				}
			}
*/
			// Select the protein with the highest number of matched peaks
			vector<unsigned int> protMatches(db.size());
			unsigned int vIdx, bestCount=0, bestIdx=-1;
			int protIdx;
			list<TwoValues<int> >::iterator peakIter;
			for(protIdx=0; protIdx<(int)protMatches.size(); protIdx++) protMatches[protIdx]=0;
			for(vIdx=0; vIdx<abVertices[cIdx].size(); vIdx++)
				for(peakIter = abVertices[cIdx][vIdx].begin(); peakIter!=abVertices[cIdx][vIdx].end(); peakIter++) {
					protIdx = matchSet.specProtMatches[(*peakIter)[0]][0];
					if(protIdx>=0) {
						protMatches[protIdx]++;
						if(bestCount<protMatches[protIdx]) { bestCount=protMatches[protIdx]; bestIdx=protIdx; }
					}
				}
			curMatch.proteinIdx = bestIdx;

			// Determine ABruijn-vertex/protein matches from spectrum-peak/protein matches
			if(pathSpectra.consensus[cIdx].size()>0) {
				for(vIdx=0; vIdx<abVertices[cIdx].size(); vIdx++) {
					for(peakIter = abVertices[cIdx][vIdx].begin(); peakIter!=abVertices[cIdx][vIdx].end(); peakIter++) {
						if(matchSet.specProtMatches[(*peakIter)[0]][0]==curMatch.proteinIdx) { // must match same protein
							for(peakIdx=0; peakIdx<matchSet.matchedPeaks[(*peakIter)[0]].size(); peakIdx++)
								if(matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][0]==(*peakIter)[1]) break;   // Found a peak in an ABruijn vertex that was initially matched to the protein sequence
																											   // Possible problem: different peaks in the same ABruijn vertex may match different protein locations - just taking the first match may not be the best option
							if(peakIdx<matchSet.matchedPeaks[(*peakIter)[0]].size() and matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][0]==(*peakIter)[1]) {
								curPairIdx.set(vIdx,(int)round(matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][1]));
								curMatch.matchedIndices.push_front(curPairIdx);
								curPairMasses.set(pathSpectra.consensus[cIdx][vIdx][0],db.masses[curMatch.proteinIdx][curPairIdx[1]][0]);
								curMatch.matchedMasses.push_front(curPairMasses);
								break;
							}
						}
					}
				}
			}

			if(not curMatch.matchedIndices.empty()) {
				curMatch.aaStart = curMatch.matchedIndices.back()[1];
				curMatch.aaEnd   = curMatch.matchedIndices.front()[1];
			}
cerr<<" -->A> matches["<<cIdx<<"] = "<<matches[cIdx].matches[0].size()<<" hits\n";
			matches[cIdx].addMatch(0,curMatch);
cerr<<" -->B> matches["<<cIdx<<"] = "<<matches[cIdx].matches[0].size()<<" hits\n";
			csps_matches.set(cIdx,curMatch);
		}

		abg.info_heaviestPath(cStats[cIdx]);   abCounts[cIdx][0]=abg.numVertices();   abCounts[cIdx][1]=abg.numEdges();
#ifdef DBG_HOMGLUE
		cout<<"  - Heaviest path stats: ["; for(unsigned int i=0;i<cStats[cIdx].size();i++) {cout<<cStats[cIdx][i]; if(i<cStats[cIdx].size()-1) cout<<", "; } cout<<"]\n";
		sprintf(sBuf,"graph_%d.txt",cIdx); abg.output_graphviz(sBuf);
		sprintf(sBuf,"graph_ma_%d.txt",cIdx);   vSet.output_graphviz_ma(sBuf, *pathVertsIdx);
#endif

		//
		// Find spectra with no ABruijn vertices on the heaviest path (if any)
		//   and remove them from the current component
		//
		list<int> usedSpecs;   usedSpecs.clear();
		for(unsigned int i=0; i<abVertices[cIdx].size(); i++)
			for(list<TwoValues<int> >::iterator vIter=abVertices[cIdx][i].begin(); vIter!=abVertices[cIdx][i].end(); vIter++)
				usedSpecs.push_back((*vIter)[0]);
		usedSpecs.sort();   usedSpecs.unique();
		if(usedSpecs.size()>0 and usedSpecs.size()<components.sets[cIdx].size()-1) {  // Whenever there are at least 2 unused spectra
			cout<<"  - Keeping "<<usedSpecs.size()<<" spectra; number of components: "<<components.size()<<" -> ";
			components.splitSet(cIdx,usedSpecs);
			cout<<components.size()<<"\n";
		}
cerr<<" -->C> matches["<<cIdx<<"] = "<<matches[cIdx].matches[0].size()<<" hits\n";
		cout.flush();   cerr.flush();
	}

	// Resize down to the final number of resulting connected components
	abCounts.resize(components.size());
	cStats.resize(components.size());
	pathSpectra.resize(components.size());
	cSpectra.resize(components.size());
	abVertices.resize(components.size());
    matches.resize(components.size());
    csps_matches.resize(components.size());

    //
    // Output results
    //
//    cout << "Done searching. Beginning output... "; cout.flush();
    BufferedLineReader blr;
    if(params.paramPresent("INPUT_SPECS_NAMES")) blr.Load(params.getValue("INPUT_SPECS_NAMES"));
    char sep=';';
    ofstream output;   output.open(params.getValue("OUTPUT_CSV"), ios::binary);
//    output << "FASTA file: " << params.getValue("INPUT_FASTA") << endl;
    output << "SPS contig"<<sep<<"CSPS index"<<sep<<"#mods"<<sep<<"Direction"<<sep<<"Match score"<<sep<<"Protein index"<<sep<<"Match start"<<sep<<"Match end"<<sep<<"Match length"<<sep<<"Matched sequence"<<sep<<"Matched peaks"<<sep<<"Protein reference\n";
    for(int i=0; i<components.size(); i++) {
//    	cout<<" --> Outputting component "<<i<<" with "<<matches[i].matches.size()<<" matches"; cout.flush();
//    	if(not matches[i].matches.empty() and not matches[i].matches[0].empty()) cout<<" ("<<matches[i].matches[0].front().matchedIndices.size()<<" peak/protein matches)\n"; else cout<<"\n"; cout.flush();
//    	output<< "Spectrum " << i << endl;
    	if(i<blr.size()) output << blr.getline(i); output << sep;
    	matches[i].output_csv(output,db,peakTol,sep,0,i);
//    	output << endl;
    }
    output.close();

	// Output modified spectrum/protein matches (after transferring some matches to the reference protein)
	if(params.paramPresent("OUTPUT_MATCHES_CSPS")) {
		csps_matches.matchedPeaks.SaveSpecSet_pklbin((string(params.getValue("OUTPUT_MATCHES_CSPS"))+string("_midx.pklbin")).c_str());
	    Save_binArray((string(params.getValue("OUTPUT_MATCHES_CSPS"))+string("_mp.bin")).c_str(), csps_matches.specProtMatches);
	}

    components.saveas_binListArray("components.bla");
	Save_binArray("component_stats.bna",cStats);
	Save_binListArray<int,list<int>,list<int>::iterator>("component_spectra.bla",cSpectra);
	Save_abinfo("component_info.bin",specSet, components.sets, specFlipped, abVertices);
	pathSpectra.Save("path_spectra_as_cluster.txt");
	pathSpectra.consensus.SaveSpecSet_pklbin(params.getValue("OUTPUT_SPECS"));

//    cout << "done.\n";
}

