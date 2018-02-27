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

int main(int argc, char **argv) {
    DB_fasta db;
    vector<DB_fasta> exonAlleles(2);
    Clusters clusters;
    SpecSet testSpecs;
    int modIdx, specIdx, protIdx, exonIdx;

	InputParams params; 	vector<char *> paramStrs;   paramStrs.resize(4);
	paramStrs[0] = "INPUT_CLUSTERS";
	paramStrs[1] = "INPUT_FASTA";
	paramStrs[2] = "OUTPUT_CSV";
	paramStrs[3] = "ENFORCE_ENDPEAKS";

	if(argc==1)	params.readParams("matchma.params"); else params.readParams(argv[1]);
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"matchma.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_CLUSTERS, INPUT_FASTA, OUTPUT_CSV, ENFORCE_ENDPEAKS\n";
		return -1;
	}

	char *clustersFN = params.getValue(paramStrs[0]);
	char *fastaFN = params.getValue(paramStrs[1]);
	char *outputFN = params.getValue(paramStrs[2]);
	bool enforceEndpeaks = params.getValueInt("ENFORCE_ENDPEAKS")>0;
	unsigned int startIdx = (params.paramPresent("IDX_START")?max(0,params.getValueInt("IDX_START")):0);
	unsigned int endIdx = (params.paramPresent("IDX_END")?max(0,params.getValueInt("IDX_END")):0);
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float maxModMass=372.2;   if(params.paramPresent("MAX_MOD_MASS")) maxModMass = (float)max(0.0,params.getValueDouble("MAX_MOD_MASS"));
	float minModMass=-372.2;  if(params.paramPresent("MIN_MOD_MASS")) minModMass = max(minModMass,(float)params.getValueDouble("MIN_MOD_MASS"));
	int showMatchedPRMs=0;	  if (params.paramPresent("SHOW_MATCHED_PRMS")) showMatchedPRMs=params.getValueInt("SHOW_MATCHED_PRMS");
	int maxNumMods = (params.paramPresent("MAX_NUM_MODS")?max(0,params.getValueInt("MAX_NUM_MODS")):1);
	int minNumMatchPeaks = (params.paramPresent("MIN_NUM_MATCH_PEAKS")?max(0,params.getValueInt("MIN_NUM_MATCH_PEAKS")):4);
	int matchesPerMod = (params.paramPresent("MATCHES_PER_MOD")?max(0,params.getValueInt("MATCHES_PER_MOD")):2);
	if(params.paramPresent("AMINO_ACID_MASSES")) { AAJumps jumps(-1); jumps.loadJumps(params.getValue("AMINO_ACID_MASSES"),true); }

  // Read FASTA sequences
    int count = db.Load(fastaFN);  if(count<=0) { cerr<<"Error reading "<<fastaFN<<"!\n"; return(-1); }

  // Read single spectrum files + cluster information
    count = clusters.Load(clustersFN);   if(count<=0) { cerr<<"Error reading "<<clustersFN<<"!\n"; return(-1); }
	if(not params.paramPresent("IDX_END")) endIdx = clusters.size()-1; else endIdx = min(endIdx,clusters.size()-1);

  // Read consensus MA spectra
    Spectrum tmpSpec;                 tmpSpec.peakList.reserve(1024);
    Spectrum cSpec;                   cSpec.peakList.reserve(1024);
    Spectrum cSpecRev;                cSpecRev.peakList.reserve(1024);
    vector<AMM_match_spec> matches(clusters.size());
    AMM_match curMatch;
    vector<vector<vector<AMM_match> > > exonBestMatches; // Used to keep best matches per exon type (e.g., VMax)
						//  first dimension is exon type
						//  second dimension is number of mods (zero to maxNumMods)
						//  third dimension is number of spectrum peaks
    SpecSet matchedIdx(clusters.size()), selectedSpec(clusters.size());   AMM_match *bestMatch;
    vector<vector<int> > proteinMatches(clusters.size());  // Best match[i]: protein / #mods / b/y-match (col 0/1/2)
    for(specIdx=0; specIdx<proteinMatches.size(); specIdx++)
    	{ proteinMatches[specIdx].resize(3); proteinMatches[specIdx][0]=-1; proteinMatches[specIdx][1]=-1; proteinMatches[specIdx][2]=-1; }

    count = exonAlleles[0].Load("test_alleles1.fasta");   if(count<=0) { cerr<<"Error reading test_alleles1.fasta!\n"; return(-1); }
    count = exonAlleles[1].Load("test_alleles2.fasta");   if(count<=0) { cerr<<"Error reading test_alleles2.fasta!\n"; return(-1); }
    exonBestMatches.resize(2);
    for(specIdx=startIdx; specIdx<=endIdx; specIdx++) {
    	matches[specIdx].init(matchesPerMod,maxNumMods);
    	cout << "Spectrum " << specIdx; cout.flush();
        clusters.getSpecIfB(specIdx,cSpec);	    clusters.getSpecIfY(specIdx,cSpecRev);

        cout << ": matching as b..."; cout.flush();
        curMatch.orientationPRMs=0;
        for(protIdx=0; protIdx<db.size(); protIdx++) {
            curMatch.proteinIdx=protIdx;
//            scoreOverlapAMM(cSpec, db.getMassesSpec(protIdx), maxNumMods, minNumMatchPeaks, matches[specIdx], curMatch, pmTol, peakTol, 57, maxModMass, minModMass, enforceEndpeaks);
        }

/*
		// Standard spectrum/protein matches
        cout << " done. Matching as y..."; cout.flush();
        curMatch.orientationPRMs=1;
		for(protIdx=0; protIdx<db.size(); protIdx++) {
			curMatch.proteinIdx=protIdx;
	            scoreOverlapAMM(cSpecRev, db.getMassesSpec(protIdx), maxNumMods, minNumMatchPeaks, matches[specIdx], curMatch, pmTol, peakTol, 57, maxModMass, minModMass, enforceEndpeaks);
		}
*/
        // Multi-exon matches
        cout << " done. Matching as y..."; cout.flush();
        curMatch.orientationPRMs=1;
    	for(exonIdx=0; exonIdx<exonAlleles.size(); exonIdx++) {
    	    exonBestMatches[exonIdx].resize(maxNumMods+1);
			for(modIdx=0; modIdx<=maxNumMods; modIdx++) {
				exonBestMatches[exonIdx][modIdx].resize(cSpecRev.size());
				for(unsigned int peakIdx=0; peakIdx<cSpecRev.size(); peakIdx++)
					exonBestMatches[exonIdx][modIdx][peakIdx].reset();
			}

			for(protIdx=0; protIdx<exonAlleles[exonIdx].size(); protIdx++) {
				curMatch.proteinIdx=protIdx;
		        curMatch.exonAllele.set(exonIdx,protIdx);
cerr<<"matching ("<<curMatch.exonAllele[0]<<","<<curMatch.exonAllele[1]<<")\n";
				if(exonIdx==0)
					scoreOverlapAMMme(cSpecRev, exonAlleles[exonIdx].getMassesSpec(protIdx), maxNumMods, minNumMatchPeaks, matches[specIdx], curMatch, 0, &exonBestMatches[exonIdx], pmTol, peakTol, 57, maxModMass, minModMass, enforceEndpeaks);
				else
					scoreOverlapAMMme(cSpecRev, exonAlleles[exonIdx].getMassesSpec(protIdx), maxNumMods, minNumMatchPeaks, matches[specIdx], curMatch, &exonBestMatches[exonIdx-1], &exonBestMatches[exonIdx], pmTol, peakTol, 57, maxModMass, minModMass, enforceEndpeaks);
			}
/*
for(modIdx=0; modIdx<=maxNumMods; modIdx++) {
	cerr<<modIdx<<" mods: ";
	for(unsigned int peakIdx=0; peakIdx<cSpecRev.size(); peakIdx++)
		cerr<<exonBestMatches[exonIdx][modIdx][peakIdx].matchScore<<"\t";
	cerr<<endl;
}
*/
    	}

        cout << " done.\n"; cout.flush();

/*
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
*/
    }

    // Output results
    AMM_match *prevExon;
    cout << "Done searching. Beginning output... "; cout.flush();
	char sep=';';
    ofstream output;   output.open(outputFN);  // output.open("matchma_debug_output.txt");
    output << "FASTA file: " << fastaFN << endl;
    output << "Direction"<<sep<<"Match score"<<sep<<"Protein index"<<sep<<"Match start"<<sep<<"Match end"<<sep<<"Matched sequence"<<sep<<"Matched peaks"<<sep<<"Protein reference\n";
//    for(int i=startIdx; i<=endIdx; i++) {
//    	output<< "Spectrum " << i << endl; matches[i].output_csv(output,db,peakTol,sep,showMatchedPRMs); output << endl;
//    }
    for(specIdx=startIdx; specIdx<=endIdx; specIdx++) {
    	output<< "Spectrum " << specIdx << endl;
        for(unsigned int modIdx=0; modIdx<=maxNumMods; modIdx++) {
        	output << "Num mod/muts = " << modIdx << endl;
            deque<AMM_match>::reverse_iterator p = matches[specIdx].matches[modIdx].rbegin();
cerr<<"output ("<<p->exonAllele[0]<<","<<p->exonAllele[1]<<")\n";
            while(p!=matches[specIdx].matches[modIdx].rend()) {
            	prevExon = p->output_csv(output,exonAlleles[p->exonAllele[0]],peakTol,sep,showMatchedPRMs);
            	while(prevExon) {
cerr<<" -- dereferencing to <"<<prevExon->exonAllele[0]<<","<<prevExon->exonAllele[1]<<")\n";
            		output<<"+++"; output.flush();
            		prevExon = prevExon->output_csv(output,exonAlleles[prevExon->exonAllele[0]],peakTol,sep,showMatchedPRMs);
            	}
            	p++;
            }
        }
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
