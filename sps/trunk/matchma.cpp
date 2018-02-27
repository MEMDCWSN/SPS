// Main file for matchma - Match Multiple Alignments against a
//   database of fasta protein sequences and find the best match

#include "db_fasta.h"
#include "clusters.h"
#include "alignment_modmut.h"
#include "spectrum.h"
#include "inputParams.h"

#include <iostream>
#include <fstream>
#include <deque>
using namespace std;
using namespace specnets;

int main(int argc, char **argv) {
    DB_fasta db;
    vector<DB_fasta> exonAlleles(2);
    Clusters clusters;
    SpecSet testSpecs;
    int modIdx, specIdx, protIdx, exonIdx;

	InputParams params; 	vector<const char *> paramStrs;   paramStrs.resize(3);
	paramStrs[0] = "INPUT_FASTA";
	paramStrs[1] = "OUTPUT_CSV";
	paramStrs[2] = "ENFORCE_ENDPEAKS";

	if(argc==1)	params.readParams("matchma.params"); else params.readParams(argv[1]);
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"matchma.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_FASTA, OUTPUT_CSV, ENFORCE_ENDPEAKS\n";
		return -1;
	}

	const char *fastaFN = params.getValue("INPUT_FASTA");
	const char *outputFN = params.getValue("OUTPUT_CSV");
	bool enforceEndpeaks = params.getValueInt("ENFORCE_ENDPEAKS")>0;
	unsigned int startIdx = (params.paramPresent("IDX_START")?max(0,params.getValueInt("IDX_START")):0);
	unsigned int endIdx = (params.paramPresent("IDX_END")?max(0,params.getValueInt("IDX_END")):0);
	unsigned int startIdx_db = (params.paramPresent("IDX_START_DB")?max(0,params.getValueInt("IDX_START_DB")):0);
	unsigned int endIdx_db = (params.paramPresent("IDX_END_DB")?max(0,params.getValueInt("IDX_END_DB")):0);
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float maxModMass=372.2;   if(params.paramPresent("MAX_MOD_MASS")) maxModMass = (float)max(0.0,params.getValueDouble("MAX_MOD_MASS"));
	float minModMass=-372.2;  if(params.paramPresent("MIN_MOD_MASS")) minModMass = max(minModMass,(float)params.getValueDouble("MIN_MOD_MASS"));
	int showMatchedPRMs=0;	  if (params.paramPresent("SHOW_MATCHED_PRMS")) showMatchedPRMs=params.getValueInt("SHOW_MATCHED_PRMS");
	int maxNumMods = (params.paramPresent("MAX_NUM_MODS")?max(0,params.getValueInt("MAX_NUM_MODS")):1);
	int minNumMatchPeaks = (params.paramPresent("MIN_NUM_MATCH_PEAKS")?max(0,params.getValueInt("MIN_NUM_MATCH_PEAKS")):4);
	int matchesPerMod = (params.paramPresent("MATCHES_PER_MOD")?max(0,params.getValueInt("MATCHES_PER_MOD")):2);
	if(params.paramPresent("AMINO_ACID_MASSES")) { AAJumps jumps(-1); jumps.loadJumps(params.getValue("AMINO_ACID_MASSES"),true); }
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;

  // Read FASTA sequences
    int count = db.Load(fastaFN);  if(count<=0) { cerr<<"Error reading "<<fastaFN<<"!\n"; return(-1); }
    if(not params.paramPresent("IDX_END_DB")) endIdx_db = db.size()-1; else endIdx_db = min(endIdx_db,db.size()-1);

  // Read single spectrum files + cluster information
	if(params.paramPresent("INPUT_CLUSTERS")) count=clusters.Load(params.getValue("INPUT_CLUSTERS"));
	else if(params.paramPresent("INPUT_CLUSTERS_PKLBIN")) count=clusters.Load_pklbin(params.getValue("INPUT_CLUSTERS_PKLBIN"));
	if (count<=0 or clusters.size()==0) return -1;
    if(not params.paramPresent("IDX_END")) endIdx = clusters.size()-1; else endIdx = min(endIdx,clusters.size()-1);

  // Read consensus MA spectra
    Spectrum tmpSpec;                 tmpSpec.peakList.reserve(1024);
    Spectrum cSpec;                   cSpec.peakList.reserve(1024);
    Spectrum cSpecRev;                cSpecRev.peakList.reserve(1024);
    vector<AMM_match_spec> matches(clusters.size());
    AMM_match curMatch;
    SpecSet matchedIdx(clusters.size()), selectedSpec(clusters.size());   AMM_match *bestMatch;
    vector<vector<int> > proteinMatches(clusters.size());  // Best match[i]: protein / #mods / b/y-match (col 0/1/2)
    for(specIdx=0; specIdx<proteinMatches.size(); specIdx++)
    	{ proteinMatches[specIdx].resize(3); proteinMatches[specIdx][0]=-1; proteinMatches[specIdx][1]=-1; proteinMatches[specIdx][2]=-1; }

    for(specIdx=startIdx; specIdx<=endIdx; specIdx++) {
    	matches[specIdx].init(matchesPerMod,maxNumMods);
    	cout << "Spectrum " << specIdx; cout.flush();
        clusters.getSpecIfB(specIdx,cSpec,peakTol);	    clusters.getSpecIfY(specIdx,cSpecRev,peakTol,ionOffset);

        cout << ": matching as b..."; cout.flush();
        curMatch.orientationPRMs=0;
        for(protIdx=startIdx_db; protIdx<=endIdx_db; protIdx++) {
            curMatch.proteinIdx=protIdx;
            scoreOverlapAMM(cSpec, db.getMassesSpec(protIdx), maxNumMods, minNumMatchPeaks, matches[specIdx], curMatch, pmTol, peakTol, 57, maxModMass, minModMass, enforceEndpeaks);
        }

        cout << " done. Matching as y..."; cout.flush();
        curMatch.orientationPRMs=1;
		for(protIdx=startIdx_db; protIdx<=endIdx_db; protIdx++) {
			curMatch.proteinIdx=protIdx;
			scoreOverlapAMM(cSpecRev, db.getMassesSpec(protIdx), maxNumMods, minNumMatchPeaks, matches[specIdx], curMatch, pmTol, peakTol, 57, maxModMass, minModMass, enforceEndpeaks);
		}
        cout << " done.\n"; cout.flush();

        // Get indices of matched peaks for best match
        unsigned int modIdx=0;
        while(modIdx<=maxNumMods and matches[specIdx].matches[modIdx].empty()) modIdx++;
        if(modIdx<=maxNumMods) {
	        bestMatch=&(matches[specIdx].matches[modIdx].back());   proteinMatches[specIdx][1]=modIdx++;
cerr<<"  - Best match with "<<modIdx-1<<" mods and score "<<bestMatch->matchScore<<"\n";
	        for(; modIdx<=maxNumMods; modIdx++)
	        	if(not matches[specIdx].matches[modIdx].empty() and matches[specIdx].matches[modIdx].back().matchScore > bestMatch->matchScore+0.0001)
	        		{ bestMatch = &(matches[specIdx].matches[modIdx].back()); proteinMatches[specIdx][1]=modIdx;
cerr<<"  - Best match with "<<modIdx<<" mods and score "<<bestMatch->matchScore<<"\n";
	        		}
	        proteinMatches[specIdx][0]=bestMatch->proteinIdx;   proteinMatches[specIdx][2]=bestMatch->orientationPRMs;
	        if(bestMatch->orientationPRMs==0) selectedSpec[specIdx]=cSpec; else selectedSpec[specIdx]=cSpecRev;
	        matchedIdx[specIdx].resize(bestMatch->matchedIndices.size());
cerr<<"  - Selected match to protein "<<bestMatch->proteinIdx<<", #mods = "<<modIdx<<", score "<<bestMatch->matchScore<<", #peaks = "<<bestMatch->matchedIndices.size()<<"\n";
	        deque<TwoValues<int> >::reverse_iterator peaksIter = bestMatch->matchedIndices.rbegin();
	        for(unsigned int peakIdx=0; peakIdx<matchedIdx[specIdx].size(); peakIdx++, peaksIter++)
	        	matchedIdx[specIdx][peakIdx].set((float)(*peaksIter)[0],(float)(*peaksIter)[1]);
        }
    }

    // Output results
    AMM_match *prevExon;
    cout << "Done searching. Beginning output... "; cout.flush();
	char sep=';';
    ofstream output;   output.open(outputFN, ios::binary);  // output.open("matchma_debug_output.txt", ios::binary);
    output << "FASTA file: " << fastaFN << endl;
    output << "Direction"<<sep<<"Match score"<<sep<<"Protein index"<<sep<<"Match start"<<sep<<"Match end"<<sep<<"Matched sequence"<<sep<<"Matched peaks"<<sep<<"Protein reference\n";
    for(int i=startIdx; i<=endIdx; i++) {
    	output<< "Spectrum " << i << endl; matches[i].output_csv(output,db,peakTol,sep,showMatchedPRMs); output << endl;
    }
    output.close();

    if(params.paramPresent("OUTPUT_MATCHED_PEAKS_IDX"))
    	matchedIdx.SaveSpecSet_pklbin(params.getValue("OUTPUT_MATCHED_PEAKS_IDX"));

    if(params.paramPresent("OUTPUT_MATCHED_PROTS"))
    	Save_binArray(params.getValue("OUTPUT_MATCHED_PROTS"), proteinMatches);

    if(params.paramPresent("OUTPUT_MATCHED_SPECS"))
    	selectedSpec.SaveSpecSet_pklbin(params.getValue("OUTPUT_MATCHED_SPECS"));

    cout << "done.\n";
}
