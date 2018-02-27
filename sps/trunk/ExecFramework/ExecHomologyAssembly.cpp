// Header Includes
#include "ExecHomologyAssembly.h"

// Module Includes
#include "Logger.h"
//#include "FileUtils.h"

// SpecNets Includes
#include "alignment_modmut.h"
#include "abruijn.h"
#include "graph.h"
#include "SetMerger.h"
#include "SpectrumPairSet.h"
#include "SpectralPairs.h"

#define DEBUG_HOMOLOGY 0

using namespace std;
using namespace specnets;

namespace specnets
{
	ExecHomologyAssembly::ExecHomologyAssembly(void) :
		m_spectra(0x0), m_db(0x0), ownInput(true),
		m_contigs(0x0), m_cspsMatchedIndices(0x0), m_cspsMatchedProts(0x0), ownOutput(true)
	{
		m_name = "ExecHomologyAssembly";
		m_type = "ExecHomologyAssembly";
	}

	// -------------------------------------------------------------------------

	ExecHomologyAssembly::ExecHomologyAssembly(const ParameterList & inputParams) :
		ExecBase(inputParams),
		m_spectra(0x0), m_db(0x0), ownInput(true),
		m_contigs(0x0), m_cspsMatchedIndices(0x0), m_cspsMatchedProts(0x0), ownOutput(true)

	{
		m_name = "ExecHomologyAssembly";
		m_type = "ExecHomologyAssembly";
	}

	// -------------------------------------------------------------------------

	ExecHomologyAssembly::ExecHomologyAssembly(const ParameterList & inputParams,
			SpecSet * spectra,
			DB_fasta * db,
			SpecSet * contigSpectra,
			SpecSet * cspsMatchedIndices,
			vector<vector<int> > * cspsMatchedProts) :
	  ExecBase(inputParams),
			m_spectra(spectra),
			m_db(db),
			ownInput(false),
			m_contigs(contigSpectra),
			m_cspsMatchedIndices(cspsMatchedIndices),
			m_cspsMatchedProts(cspsMatchedProts),
			ownOutput(false)
	{
		m_name = "ExecHomologyAssembly";
		m_type = "ExecHomologyAssembly";
	}

	// -------------------------------------------------------------------------

	ExecHomologyAssembly::~ExecHomologyAssembly(void)
	{
		if (ownInput)
		{
			if(m_spectra) delete m_spectra;
			if(m_db) delete m_db;
		}
		if (ownOutput)
		{
			if(m_contigs) delete m_contigs;
			if(m_cspsMatchedIndices) delete m_cspsMatchedIndices;
			if(m_cspsMatchedProts) delete m_cspsMatchedProts;
		}
	}

	// -------------------------------------------------------------------------

	ExecBase * ExecHomologyAssembly::clone(const ParameterList & inputParams) const
	{
		return new ExecHomologyAssembly(inputParams);
	}

	// -------------------------------------------------------------------------

	bool ExecHomologyAssembly::invoke(void)
	{
		if(!m_spectra or m_spectra->size()==0) {
			ERROR_MSG("ERROR: empty set of input spectra");
			return false;
		}

		if(!m_db or m_db->size()==0) {
			ERROR_MSG("ERROR: empty database");
			return false;
		}

//    Save_abinfo(m_params.getValue("OUTPUT_CSPS_ABRUIJN").c_str(), *m_spectra, m_components.sets, m_specFlipped, m_abVertices);

    float pmTol = (float) m_params.getValueDouble("TOLERANCE_PM");
    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");

    unsigned int startIdx = m_params.exists("IDX_START") ? max(0,m_params.getValueInt("IDX_START")) : 0;
    unsigned int endIdx = m_params.exists("IDX_END") ? min(m_spectra->size()-1, (unsigned int) m_params.getValueInt("IDX_END")) : m_spectra->size()-1;
    unsigned int startIdx_db = m_params.exists("IDX_START_DB") ? max(0,m_params.getValueInt("IDX_START_DB")) :0;
    unsigned int endIdx_db = m_params.exists("IDX_END_DB") ? min(m_db->size()-1, (unsigned int) m_params.getValueInt("IDX_END_DB")) : m_db->size()-1;
    if(startIdx > m_spectra->size()) return false;
    if(startIdx_db > m_db->size()) return false;

  	int graphType = m_params.exists("GRAPH_TYPE")?m_params.getValueInt("GRAPH_TYPE"):2;
  	short edgeScoreType = (short)(m_params.exists("EDGE_SCORE_TYPE")?m_params.getValueInt("EDGE_SCORE_TYPE"):1);
  	int minContigSet = m_params.exists("MIN_CONTIG_SET")?m_params.getValueInt("MIN_CONTIG_SET"):1;
  	int maxAAjump = m_params.exists("MAX_AA_JUMP")?m_params.getValueInt("MAX_AA_JUMP"):2;
  	float penalty_ptm = m_params.exists("PENALTY_PTM")?(float)m_params.getValueDouble("PENALTY_PTM"):0;
  	float penalty_sameVert = m_params.exists("PENALTY_SAME_VERTEX")?(float)m_params.getValueDouble("PENALTY_SAME_VERTEX"):-1000000;
    float maxModMass = m_params.exists("MAX_MOD_MASS") ? (float)max(0.0,m_params.getValueDouble("MAX_MOD_MASS")) : 372.2;
    int specType = m_params.exists("SPEC_TYPE_MSMS") ? m_params.getValueInt("SPEC_TYPE_MSMS")>0 : 0;
    float ionOffset = specType ? AAJumps::massHion : 0;
    unsigned int cIdx, pairIdx, peakIdx, specIdx;

  	bool specNets = m_params.exists("INPUT_SPECNETS_ALIGNS") and m_params.exists("INPUT_SPECNETS_MATCHED");

    DEBUG_TRACE;

  	//
  	// Process Contig <-> Protein Alignments
  	//
    AMM_match_set matchSet(*m_db,m_spectra->size());
    matchSet.SetMatches(m_spectra);

    DEBUG_TRACE;

    for(unsigned int i=0;i<m_db->size();i++) m_db->getMassesSpec(i);

//for(unsigned int matchIdx=0; matchIdx<m_psmSet->size(); matchIdx++) {
//	cerr<<" -- ExecHomologyAssembly: m_psmSet["<<matchIdx<<"][0] = "<< (*m_psmSet)[matchIdx][0]<<endl;
//}

  	// Load clustalw protein <-> protein alignments (if any)
  	if(m_params.exists("INPUT_HOMOLOGIES")) {
  		vector<TwoValues<int> > matchedProts;
  		vector<vector<TwoValues<int> > > matchedAAs;
  		if(not Load_clustalw_multiple(m_params.getValue("INPUT_HOMOLOGIES").c_str(), *m_db, matchedProts, matchedAAs)) {
  			DEBUG_MSG("ERROR loading "<<m_params.getValue("INPUT_HOMOLOGIES"));
  			return false;
  		}

  		// Use clustalw alignments to transitively import contig matches to the reference protein
  		for(unsigned int matchIdx=0; matchIdx<matchedProts.size(); matchIdx++) {
cerr<<" -- ExecHomologyAssembly: matchedProts["<<matchIdx<<"][0/1] = "<<matchedProts[matchIdx][0]<<"/"<<matchedProts[matchIdx][1]<<", matchedAAs["<<matchIdx<<"].size() = "<<matchedAAs[matchIdx].size()<<endl;
  			matchSet.MergeIntoReference(matchedProts[matchIdx][0],matchedProts[matchIdx][1],matchedAAs[matchIdx]);
  		}
  	}

  	// Output modified spectrum/protein matches (after transferring some matches to the reference protein)
  	if(m_params.exists("OUTPUT_MATCHES_REF_MIDX"))
  		matchSet.matchedPeaks.savePklBin(m_params.getValue("OUTPUT_MATCHES_REF_MIDX").c_str());
  	if(m_params.exists("OUTPUT_MATCHES_REF_MP"))
  		Save_binArray(m_params.getValue("OUTPUT_MATCHES_REF_MP").c_str(), matchSet.specProtMatches);

  	m_specFlipped.resize(m_spectra->size());
  	for(unsigned int specIdx=0; specIdx<m_specFlipped.size(); specIdx++)
  		m_specFlipped[specIdx]=matchSet.specProtMatches[specIdx][0]==-1?false:(bool)matchSet.specProtMatches[specIdx][2];

  	//
  	// Sequence resulting contigs
  	//
  	SpectrumPairSet pairsASP(0);
  	SpectrumPairSet pairsPA(0);
      vector<vector<TwoValues<int> > > pairMatchedPeaks(0),pairMatchedPeaksASP(0);
      if(not m_params.exists("SKIP_GLUES") or m_params.getValueInt("SKIP_GLUES")!=1)
      	matchSet.GenerateGlues(pairsPA,pairMatchedPeaks);

  	// Load pairwise alignments from spectral networks (if available)
  	if(specNets) {
  //		if (!Load_results_bin(m_params.getValue("INPUT_SPECNETS_ALIGNS"), pairsASP)) { cerr << "Error reading "<<m_params.getValue("INPUT_SPECNETS_ALIGNS")<<"!\n"; return -1; }
  //		cout << "Loading spectral networks pairs complete (" << pairsASP.size() << " pairs)\n";
  		if (!pairsASP.loadFromBinaryFile(m_params.getValue("INPUT_SPECNETS_ALIGNS"))) { ERROR_MSG("Error reading "<<m_params.getValue("INPUT_SPECNETS_ALIGNS")<<"!"); return -1; }
  		DEBUG_MSG("Loading spectral networks pairs complete (" << pairsASP.size() << " pairs)");

  		SpecSet snMatched;
  	    if(snMatched.loadPklBin(m_params.getValue("INPUT_SPECNETS_MATCHED").c_str())<=0)
  			{ ERROR_MSG("Error loading "<<m_params.getValue("INPUT_SPECNETS_MATCHED")); return -1; }
  	    DEBUG_MSG("Loading spectral network matched peaks complete (" << snMatched.size() << " pairs)");
  		pairMatchedPeaksASP.resize(snMatched.size());
  		for(pairIdx=0; pairIdx<snMatched.size(); pairIdx++) {
  			pairMatchedPeaksASP[pairIdx].resize(snMatched[pairIdx].size());
  			for(peakIdx=0; peakIdx<snMatched[pairIdx].size(); peakIdx++) {
  				pairMatchedPeaksASP[pairIdx][peakIdx].set((int)snMatched[pairIdx][peakIdx][0],(int)snMatched[pairIdx][peakIdx][1]);
  			}
  		}
  	}

  	if(m_params.exists("INPUT_SPECNETS_ALIGNS_ALL")) {
  		SpectrumPairSet pairsASP_all(0);
  //		if (!Load_results_bin(m_params.getValue("INPUT_SPECNETS_ALIGNS_ALL"), pairsASP_all)) { cerr << "Error reading "<<m_params.getValue("INPUT_SPECNETS_ALIGNS_ALL")<<"!\n"; return -1; }
  //		cout << "Loading full set of spectral networks pairs complete (" << pairsASP_all.size() << " pairs)\n";
  		if (!pairsASP_all.loadFromBinaryFile(m_params.getValue("INPUT_SPECNETS_ALIGNS_ALL"))) { ERROR_MSG("Error reading "<<m_params.getValue("INPUT_SPECNETS_ALIGNS_ALL")<<"!"); return -1; }
  		DEBUG_MSG("Loading full set of spectral networks pairs complete (" << pairsASP_all.size() << " pairs)");

  		// Find edges between unidentified spectra
  		SetMerger components_annot(m_spectra->size());
  		components_annot.createSets(m_spectra->size(),minContigSet,pairsASP,pairsPA);
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
      ERROR_MSG("Got "<<pairsPA.size()<<" pairs of glued spectra");
      ofstream dbg_output;   dbg_output.open("dbg_glues.txt", ios::binary);
      for(unsigned int i=0; i<pairsPA.size(); i++) {
      	dbg_output<<"Pair "<<pairsPA[i].spec1<<" / "<<pairsPA[i].spec2<<":\n";
      	for(unsigned int j=0; j<pairMatchedPeaks[i].size(); j++)
      		dbg_output<<" peaks ("<<pairMatchedPeaks[i][j][0]<<","<<pairMatchedPeaks[i][j][1]<<"), masses ("<<(*m_spectra)[pairsPA[i].spec1][pairMatchedPeaks[i][j][0]][0]<<","<<(*m_spectra)[pairsPA[i].spec2][pairMatchedPeaks[i][j][1]][0]<<")\n";
      	dbg_output<<endl;

      }
      dbg_output.close();
  #endif

      // Find connected components (aligned spectra)
  	m_components.resize(m_spectra->size());
  	m_components.createSets(m_spectra->size(),minContigSet,pairsASP,pairsPA);
  	if(specNets) { // Keep only annotated components (including singletons) and unidentified networks
  		for(unsigned int setIdx=0; setIdx<m_components.sets.size(); setIdx++)
  			if(m_components.sets[setIdx].size()==1 and matchSet.matchedPeaks[m_components.sets[setIdx].front()].size()==0)
  				m_components.removeSet(setIdx,false);
  		m_components.compressSetIndices();
  	}
  	m_components.splitAligns(pairsASP,pairsPA);
  	DEBUG_MSG("Got "<<m_components.size()<<" component(s), minimum size was "<<minContigSet<<".");

  	// Sequence ABruijn graphs
  	AAJumps jumps(maxAAjump);   char sBuf[1024];
  	unsigned int numElemsInSets=m_components.numElemsInSets();  // Maximum possible number of components (all singletons)
  	vector<MSGraph> spectrumGraphs(m_spectra->size());
  	Clusters           pathSpectra;   pathSpectra.resize(numElemsInSets);  // Keep the de-novo reconstructed heaviestPath sequences as spectra in a Cluster variable
  	vector<AMM_match_spec> matches(numElemsInSets);   // Matches between sequenced CSPS contigs and protein sequences
  	AMM_match_set csps_matches(*m_db,numElemsInSets);    // Container for all the matches between csps contigs and the protein sequence database

    DEBUG_TRACE;

  	vector<vector<float> > cStats(numElemsInSets);  for(unsigned int cIdx=0;cIdx<cStats.size();cIdx++) { cStats[cIdx].resize(9); for(unsigned int i=0;i<9;i++) cStats[cIdx][i]=0; }  // Absolute value 9 depends on MSGraph::info_heaviestPath
  	vector<list<int> > cSpectra(numElemsInSets);  // Lists of spectrum indices per component
  	m_abVertices.resize(numElemsInSets); for(unsigned int i=0;i<m_abVertices.size();i++) m_abVertices[i].resize(0);   // Keeps track of which spectrum peaks were matched (dim.3) in each ABruijn vertex (dim.2) in each component (dim.1)
      vector<vector<short> > abCounts(numElemsInSets);   // Records info on number of vertices and edges per ABruijn graph
      for(unsigned int i=0; i<abCounts.size(); i++) { abCounts[i].resize(2); abCounts[i][0]=0; abCounts[i][1]=0; }

  	for(cIdx=0; cIdx<m_components.size(); cIdx++) {

      //DEBUG_VAR(cIdx);

  		AMM_match curMatch;        // Resulting path-protein match
  		TwoValues<float> curPairMasses;    TwoValues<int> curPairIdx;
  		if(m_components.sets[cIdx].size()==1) {  // No gluing/resequencing necessary, just reproduce input contig/protein match
        if (DEBUG_HOMOLOGY) DEBUG_TRACE;
  			matches[cIdx].init(1,0);           specIdx = m_components.sets[cIdx].front();
  			curMatch.orientationPRMs=0;        curMatch.modSize=0;    curMatch.matchScore=0;
  			curMatch.matchedIndices.clear();   curMatch.matchedMasses.clear();
  			curMatch.proteinIdx = matchSet.specProtMatches[specIdx][0];
        if (DEBUG_HOMOLOGY) DEBUG_VAR(matchSet.matchedPeaks.size());
        if (DEBUG_HOMOLOGY) DEBUG_VAR(matchSet.matchedPeaks[specIdx].size());
  			for(peakIdx = 0; peakIdx < matchSet.matchedPeaks[specIdx].size(); peakIdx++) {
          if (DEBUG_HOMOLOGY) DEBUG_VAR(peakIdx);
          float peak1 = matchSet.matchedPeaks[specIdx][peakIdx][0];
          float peak2 = matchSet.matchedPeaks[specIdx][peakIdx][1];
          if (DEBUG_HOMOLOGY) DEBUG_VAR(peak1);
          if (DEBUG_HOMOLOGY) DEBUG_VAR(peak2);
  				curPairIdx.set( (int)round(peak1), (int)round(peak2) );
          if (DEBUG_HOMOLOGY) DEBUG_VAR(curPairIdx[0]);
          if (DEBUG_HOMOLOGY) DEBUG_VAR(curPairIdx[1]);
  				curMatch.matchedIndices.push_front(curPairIdx);
          if (DEBUG_HOMOLOGY) DEBUG_VAR(specIdx);
          if (DEBUG_HOMOLOGY) DEBUG_VAR(curMatch.proteinIdx);
          if (DEBUG_HOMOLOGY) DEBUG_VAR(m_spectra->size());
          if (DEBUG_HOMOLOGY) DEBUG_VAR((*m_spectra)[specIdx].size());
          int index1 = (*m_spectra)[specIdx][curPairIdx[0]][0];
          if (DEBUG_HOMOLOGY) DEBUG_VAR(index1);
          if (DEBUG_HOMOLOGY) DEBUG_VAR(curMatch.proteinIdx)
          if (DEBUG_HOMOLOGY) DEBUG_VAR(curPairIdx[1]);
          if (DEBUG_HOMOLOGY) DEBUG_VAR(m_db->masses.size());
          if (DEBUG_HOMOLOGY) DEBUG_VAR(m_db->masses[curMatch.proteinIdx].size());
          int index2 = m_db->masses[curMatch.proteinIdx][curPairIdx[1]][0];
          if (DEBUG_HOMOLOGY) DEBUG_VAR(index2);
  				curPairMasses.set(index1,index2);
          if (DEBUG_HOMOLOGY) DEBUG_TRACE;
  				curMatch.matchedMasses.push_front(curPairMasses);
  				curMatch.matchScore += (*m_spectra)[specIdx][curPairIdx[0]][1];
          if (DEBUG_HOMOLOGY) DEBUG_TRACE;
  			}

        if (DEBUG_HOMOLOGY) DEBUG_TRACE;
  			if(not curMatch.matchedIndices.empty()) {
  				curMatch.aaStart = curMatch.matchedIndices.back()[1];
  				curMatch.aaEnd   = curMatch.matchedIndices.front()[1];
  			}
        if (DEBUG_HOMOLOGY) DEBUG_TRACE;
  			matches[cIdx].addMatch(0,curMatch);
  			csps_matches.set(cIdx,curMatch);
  			pathSpectra.consensus[cIdx] = (*m_spectra)[specIdx];
  			continue;
  		}

  		DEBUG_MSG("Processing component "<<cIdx<<" ("<<m_components.sets[cIdx].size()<<" elements)\n");
  		stringstream aux; aux <<"  - spectrum indices: "; for(list<int>::iterator iter=m_components.sets[cIdx].begin(); iter!=m_components.sets[cIdx].end();iter++) aux<<(*iter)<<" "; DEBUG_MSG(aux.str());

  		SpectrumPairSet &curAligns    = m_components.cAlignsPA[cIdx];
  		vector<int> &curAlignsIdx     = m_components.cAlignsPA_idx[cIdx];
  		SpectrumPairSet &curAlignsASP = m_components.cAlignsASP[cIdx];
  		vector<int> &curAlignsASPIdx  = m_components.cAlignsASP_idx[cIdx];

  		//
  		// Build spectrum graphs for the spectra in this component
  		//
  		if(graphType>0) {
  			// alignsPA
  			for(unsigned int i=0; i<curAligns.size(); i++) {
  				if(spectrumGraphs[curAligns[i].spec1].numVerts()==0) {
  					if(graphType==1) spectrumGraphs[curAligns[i].spec1].ConnectConsecutive((*m_spectra)[curAligns[i].spec1]);
  					else spectrumGraphs[curAligns[i].spec1].ConnectJumps((*m_spectra)[curAligns[i].spec1],jumps,peakTol);
  				}
  				if(spectrumGraphs[curAligns[i].spec2].numVerts()==0) {
  					if(graphType==1) spectrumGraphs[curAligns[i].spec2].ConnectConsecutive((*m_spectra)[curAligns[i].spec2]);
  					else spectrumGraphs[curAligns[i].spec2].ConnectJumps((*m_spectra)[curAligns[i].spec2],jumps,peakTol);
  				}
  			}
  			// alignsASP
  			for(unsigned int i=0; i<curAlignsASP.size(); i++) {
  				if(spectrumGraphs[curAlignsASP[i].spec1].numVerts()==0) {
  					if(graphType==1) spectrumGraphs[curAlignsASP[i].spec1].ConnectConsecutive((*m_spectra)[curAlignsASP[i].spec1]);
  					else spectrumGraphs[curAlignsASP[i].spec1].ConnectJumps((*m_spectra)[curAlignsASP[i].spec1],jumps,peakTol);
  				}
  				if(spectrumGraphs[curAlignsASP[i].spec2].numVerts()==0) {
  					if(graphType==1) spectrumGraphs[curAlignsASP[i].spec2].ConnectConsecutive((*m_spectra)[curAlignsASP[i].spec2]);
  					else spectrumGraphs[curAlignsASP[i].spec2].ConnectJumps((*m_spectra)[curAlignsASP[i].spec2],jumps,peakTol);
  				}
  			}
  		}

  		//
  		// Build A-Bruijn graph
  		//
  		VertexSet vSet(*m_spectra,8192);
  		vector<bool> usedSpectra(m_spectra->size());  for(unsigned int i=0; i<m_spectra->size(); i++) usedSpectra[i]=false;
  		for(pairIdx=0; pairIdx<curAlignsIdx.size(); pairIdx++) {
  			if(graphType==0) vSet.addGlues(curAligns[pairIdx].spec1,curAligns[pairIdx].spec2,pairMatchedPeaks[curAlignsIdx[pairIdx]]);
  			else vSet.addGlues(curAligns[pairIdx].spec1,curAligns[pairIdx].spec2,pairMatchedPeaks[curAlignsIdx[pairIdx]],&spectrumGraphs);
  			usedSpectra[curAligns[pairIdx].spec1] = true;
  			usedSpectra[curAligns[pairIdx].spec2] = true;
  		}
  		if(curAlignsASP.size()>0 and pairMatchedPeaksASP[curAlignsASPIdx[0]].empty()) {  // Need to compute matched peaks for unidentified spectral network
  			vector<vector<TwoValues<int> > > matches;    vector<float> modPos;
  			SplitPairs(*m_spectra, curAlignsASP, peakTol, pmTol, maxAAjump, maxModMass, penalty_sameVert, penalty_ptm, matches, m_specFlipped, modPos);

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
  		DEBUG_MSG("  - ABruijn graph built on "<<specCount<<" spectra");

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
  	  stringstream aux2;
  		aux2<<"  - Abruijn graph contains "<<compositeVertexCount<<" composite vertices: ";
  		list<int>::iterator iter=compositeSet.begin();
  		for(; iter!=compositeSet.end(); iter++) aux2<<(*iter)<<" ";
  		if(compositeVertexCount>0)
  			{ aux2<<"-> splitting..."; vSet.splitComposite(spectrumGraphs,peakTol,&usedSpectra); aux2<<"done.";}
  		DEBUG_MSG(aux2.str());

  		//
  		// Add spectrum graph edges to ABruijn graph
  		//
  		if(graphType==1) { vSet.addEdges(spectrumGraphs,&usedSpectra);  vSet.consolidatePaths(); }
  		if(graphType==2) { vSet.addEdges(spectrumGraphs,&usedSpectra); }
  //		vSet.addEndpointEdges(m_components.cAlignsASP[cIdx],matches,modPos,jumps,peakTol);

  		vSet.removeEndpoints(false,peakTol);
  		vector<SpectrumPeakLabels> foo_labels;
  		vSet.buildGraph(abg,jumps,peakTol,vSet_index,foo_labels, edgeScoreType);
  //char filename[2048];   sprintf(filename,"graphs/component_%d.txt",cIdx+1);
  //abg.output_graphviz(filename);
  		curMatch.matchScore = abg.heaviestPath(path,false,&pathSpectra.consensus[cIdx],&pathVertsIdx);
  		for(unsigned int vIdx=0; vIdx<pathVertsIdx.size(); vIdx++) pathVertsIdx[vIdx] = vSet_index[pathVertsIdx[vIdx]];  // Convert simplified graph vertex indices to ABruijn vertex indices.
  		vSet.getMatchedPeaks(pathVertsIdx,m_abVertices[cIdx]);

  		//
  		// Get correspondence between protein(s) and heaviest path
  		//
  		matches[cIdx].init(1,0);
  		curMatch.orientationPRMs=0;        curMatch.modSize=0;
  		curMatch.matchedIndices.clear();   curMatch.matchedMasses.clear();
  		if(m_abVertices[cIdx].size()>0) {
  			// Select the protein with the highest number of matched peaks
  			vector<unsigned int> protMatches(m_db->size());
  			unsigned int vIdx, bestCount=0, bestIdx=0;
  			int protIdx;
  			list<TwoValues<int> >::iterator peakIter;
  			for(protIdx=0; protIdx<(int)protMatches.size(); protIdx++) protMatches[protIdx]=0;
  			for(vIdx=0; vIdx<m_abVertices[cIdx].size(); vIdx++)
  				for(peakIter = m_abVertices[cIdx][vIdx].begin(); peakIter!=m_abVertices[cIdx][vIdx].end(); peakIter++) {
  					protIdx = matchSet.specProtMatches[(*peakIter)[0]][0];
  					if(protIdx>=0) {
  						protMatches[protIdx]++;
  						if(bestCount<protMatches[protIdx]) { bestCount=protMatches[protIdx]; bestIdx=protIdx; }
  					}
  				}
  			curMatch.proteinIdx = bestIdx;

  			// Determine ABruijn-vertex/protein matches from spectrum-peak/protein matches
  			if(pathSpectra.consensus[cIdx].size()>0) {
  				for(vIdx=0; vIdx<m_abVertices[cIdx].size(); vIdx++) {
  					for(peakIter = m_abVertices[cIdx][vIdx].begin(); peakIter!=m_abVertices[cIdx][vIdx].end(); peakIter++) {
  						if(matchSet.specProtMatches[(*peakIter)[0]][0]==curMatch.proteinIdx) { // must match same protein
  							for(peakIdx=0; peakIdx<matchSet.matchedPeaks[(*peakIter)[0]].size(); peakIdx++)
  								if(matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][0]==(*peakIter)[1]) break;   // Found a peak in an ABruijn vertex that was initially matched to the protein sequence
  																											   // Possible problem: different peaks in the same ABruijn vertex may match different protein locations - just taking the first match may not be the best option
  							if(peakIdx<matchSet.matchedPeaks[(*peakIter)[0]].size() and matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][0]==(*peakIter)[1]) {
  								curPairIdx.set(vIdx,(int)round(matchSet.matchedPeaks[(*peakIter)[0]][peakIdx][1]));
  								curMatch.matchedIndices.push_front(curPairIdx);
  								curPairMasses.set(pathSpectra.consensus[cIdx][vIdx][0],m_db->masses[curMatch.proteinIdx][curPairIdx[1]][0]);
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
  			matches[cIdx].addMatch(0,curMatch);
  			csps_matches.set(cIdx,curMatch);
  		}

  		abg.info_heaviestPath(cStats[cIdx]);   abCounts[cIdx][0]=abg.numVertices();   abCounts[cIdx][1]=abg.numEdges();
  #ifdef DBG_HOMGLUE
      stringstream aux3;
  		aux3<<"  - Heaviest path stats: ["; for(unsigned int i=0;i<cStats[cIdx].size();i++) {aux3<<cStats[cIdx][i]; if(i<cStats[cIdx].size()-1) aux3<<", "; } aux<<"]"; DEBUG_MSG(aux3.str());
  		sprintf(sBuf,"graph_%d.txt",cIdx); abg.output_graphviz(sBuf);
  		sprintf(sBuf,"graph_ma_%d.txt",cIdx);   vSet.output_graphviz_ma(sBuf, *pathVertsIdx);
  #endif

  		//
  		// Find spectra with no ABruijn vertices on the heaviest path (if any)
  		//   and remove them from the current component
  		//
  		list<int> usedSpecs;   usedSpecs.clear();
  		for(unsigned int i=0; i<m_abVertices[cIdx].size(); i++)
  			for(list<TwoValues<int> >::iterator vIter=m_abVertices[cIdx][i].begin(); vIter!=m_abVertices[cIdx][i].end(); vIter++)
  				usedSpecs.push_back((*vIter)[0]);
  		usedSpecs.sort();   usedSpecs.unique();
  		if(usedSpecs.size()>0 and usedSpecs.size()<m_components.sets[cIdx].size()-1) {  // Whenever there are at least 2 unused spectra
  			DEBUG_MSG("  - Keeping "<<usedSpecs.size()<<" spectra; number of components: "<<m_components.size()<<" -> ");
  			m_components.splitSet(cIdx,usedSpecs);
  			DEBUG_MSG(m_components.size());
  		}
  	}

  	// Resize down to the final number of resulting connected components
  	abCounts.resize(m_components.size());
  	cStats.resize(m_components.size());
  	pathSpectra.resize(m_components.size());
  	cSpectra.resize(m_components.size());
  	m_abVertices.resize(m_components.size());
  	matches.resize(m_components.size());
  	csps_matches.resize(m_components.size());

  	//
  	// Output results
  	//
  	if(m_params.exists("OUTPUT_CSV")) {
  		BufferedLineReader blr;
  		if(m_params.exists("INPUT_SPECS_NAMES")) blr.Load(m_params.getValue("INPUT_SPECS_NAMES").c_str());
  		char sep=';';
  		ofstream output;   output.open(m_params.getValue("OUTPUT_CSV").c_str());
  		output << "SPS contig"<<sep<<"CSPS index"<<sep<<"#mods"<<sep<<"Direction"<<sep<<"Match score"<<sep<<"Protein index"<<sep<<"Match start"<<sep<<"Match end"<<sep<<"Match length"<<sep<<"Matched sequence"<<sep<<"Matched peaks"<<sep<<"Protein reference\n";
  		for(int i=0; i<m_components.size(); i++) {
  			if(i<blr.size()) output << blr.getline(i); output << sep;
  			matches[i].output_csv(output,*m_db,peakTol,sep,0,i);
  		}
  		output.close();
  	}

  	// Output modified spectrum/protein matches (after transferring some matches to the reference protein)
  	if(m_params.exists("OUTPUT_CLUSTERS")) pathSpectra.Save(m_params.getValue("OUTPUT_CLUSTERS").c_str());

  	if(m_contigs) {
  		m_contigs->resize(pathSpectra.consensus.size());
  		for (unsigned int i=0; i < m_contigs->size(); i++)
  			(*m_contigs)[i] = pathSpectra.consensus[i];
  	}

  	if(m_cspsMatchedIndices) {
  		m_cspsMatchedIndices->resize(csps_matches.matchedPeaks.size());
  		for (unsigned int i=0; i < m_cspsMatchedIndices->size(); i++)
  			(*m_cspsMatchedIndices)[i] = csps_matches.matchedPeaks[i];
  	}

  	if(m_cspsMatchedProts) {
  		m_cspsMatchedProts->resize(csps_matches.specProtMatches.size());
  		for (unsigned int i=0; i < m_cspsMatchedProts->size(); i++) {
				(*m_cspsMatchedProts)[i] = csps_matches.specProtMatches[i];
/*  			m_cspsMatchedIndices[i].resize(csps_matches.specProtMatches[i].size());
  			for (unsigned int j=0; j < m_cspsMatchedProts->size(); j++)
  				(*m_cspsMatchedProts)[i][j] = csps_matches.specProtMatches[i][j];
*/
  		}
  	}

  	// Deprecated output code
  	// if(m_params.exists("OUTPUT_MATCHES_CSPS")) {
  	//   csps_matches.matchedPeaks.savePklBin((string(m_params.getValue("OUTPUT_MATCHES_CSPS"))+string("_midx.pklbin")).c_str());
  	//   Save_binArray((string(m_params.getValue("OUTPUT_MATCHES_CSPS"))+string("_mp.bin")).c_str(), csps_matches.specProtMatches);
  	// }
  	// pathSpectra.consensus.savePklBin(m_params.getValue("OUTPUT_SPECS"));
  	// m_components.saveas_binListArray("components.bla");
  	// Save_binArray("component_stats.bna",cStats);
  	// Save_binListArray<int,list<int>,list<int>::iterator>("component_spectra.bla",cSpectra);

    return true;
	}

	// -------------------------------------------------------------------------

	bool ExecHomologyAssembly::loadInputData(void)
	{
		if(ownInput) {
    	if(!m_spectra) m_spectra = new SpecSet;
    	if(!m_db) m_db = new DB_fasta;
    }
	  m_spectra->resize(0);

  	if(ownOutput) {
    	if(!m_contigs) m_contigs = new SpecSet;
    	if(!m_cspsMatchedIndices) m_cspsMatchedIndices = new SpecSet;
    	if(!m_cspsMatchedProts) m_cspsMatchedProts = new vector<vector<int> >;
    }
  	m_contigs->resize(0);
  	m_cspsMatchedIndices->resize(0);
  	m_cspsMatchedProts->resize(0);

		if(m_params.exists("AMINO_ACID_MASSES")) {
			AAJumps tmpJumps(-1);
			tmpJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(),true); // Set global defaults for amino acid masses
		}

    if (!m_params.exists("INPUT_SPECS_PKLBIN") or (m_spectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str())<=0 or m_spectra->size()==0)) {
    	ERROR_MSG("Error reading input spectra from "<< m_params.getValue("INPUT_SPECS_PKLBIN"));
    	return false;
    }

    if (!m_params.exists("INPUT_FASTA")) {
			ERROR_MSG("Parameters are incomplete. INPUT_FASTA is missing.");
      return false;
    } else if (m_db->Load(m_params.getValue("INPUT_FASTA").c_str())<=0) {
			ERROR_MSG("Error reading database sequences from "<<m_params.getValue("INPUT_FASTA"));
			return false;
		}

/*
    if (!m_params.exists("INPUT_MATCHED_PEAKS_IDX") or (m_matchedIndices->LoadSpecSet_pklbin(m_params.getValue("INPUT_MATCHED_PEAKS_IDX").c_str())<=0 or m_matchedIndices->size()==0)) {
    	ERROR_MSG("Error reading input matched peaks from "<< m_params.getValue("INPUT_MATCHED_PEAKS_IDX"));
    	return false;
    }
*/
/*
    if (!m_params.exists("INPUT_MATCHED_PROTS") or Load_binArray<int>(m_params.getValue("INPUT_MATCHED_PROTS").c_str(),*m_matchedProts)<=0) {
    	ERROR_MSG("Error reading input matched peaks from "<< m_params.getValue("INPUT_MATCHED_PROTS"));
    	return false;
    }
*/

    return true;
	}

	// -------------------------------------------------------------------------

	bool ExecHomologyAssembly::saveOutputData(void)
	{
		if(m_contigs and m_params.exists("OUTPUT_SPECS"))
			m_contigs->savePklBin(m_params.getValue("OUTPUT_SPECS").c_str());

		if(m_cspsMatchedIndices and m_params.exists("OUTPUT_MATCHES_CSPS_MIDX"))
			m_cspsMatchedIndices->savePklBin(m_params.getValue("OUTPUT_MATCHES_CSPS_MIDX").c_str());

		if(m_cspsMatchedProts and m_params.exists("OUTPUT_MATCHES_CSPS_MP"))
			Save_binArray(m_params.getValue("OUTPUT_MATCHES_CSPS_MP").c_str(), *m_cspsMatchedProts);

    if(m_spectra and m_params.exists("OUTPUT_CSPS_ABRUIJN"))
      Save_abinfo(m_params.getValue("OUTPUT_CSPS_ABRUIJN").c_str(), *m_spectra, m_components.sets, m_specFlipped, m_abVertices);

		return true;
	}

	// -------------------------------------------------------------------------

	bool ExecHomologyAssembly::saveInputData(std::vector<std::string> & filenames)
	{
		return false;
	}

	// -------------------------------------------------------------------------

	bool ExecHomologyAssembly::loadOutputData(void)
	{
		return false;
	}

	// -------------------------------------------------------------------------

  vector<ExecBase*> const & ExecHomologyAssembly::split(int numSplit)
	{
		m_subModules.resize(0);
		return m_subModules;
	}

	// -------------------------------------------------------------------------

	bool ExecHomologyAssembly::merge(void)
	{
		return false;
	}

	// -------------------------------------------------------------------------

	bool ExecHomologyAssembly::validateParams(std::string & error)
	{
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
//    VALIDATE_PARAM_EXIST("IDX_START");
//    VALIDATE_PARAM_EXIST("IDX_END");
//    VALIDATE_PARAM_EXIST("IDX_START_DB");
//    VALIDATE_PARAM_EXIST("IDX_END_DB");
//    VALIDATE_PARAM_EXIST("GRAPH_TYPE");
//    VALIDATE_PARAM_EXIST("EDGE_SCORE_TYPE");
//    VALIDATE_PARAM_EXIST("MIN_CONTIG_SET");
//    VALIDATE_PARAM_EXIST("MAX_AA_JUMP");
//    VALIDATE_PARAM_EXIST("PENALTY_PTM");
//    VALIDATE_PARAM_EXIST("PENALTY_SAME_VERTEX");
//    VALIDATE_PARAM_EXIST("MAX_MOD_MASS");
//    VALIDATE_PARAM_EXIST("SPEC_TYPE_MSMS");

    m_isValid = true;
    return true;
	}

	// -------------------------------------------------------------------------


}
