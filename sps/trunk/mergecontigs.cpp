#include "spectrum.h"
#include "clusters.h"
#include "inputParams.h"
#include "batch.h"
#include "alignment_scoring.h"
#include "denovo.h"

// TODO:
//  - Finish FindConnectorSpectra (in batch.cpp)
//  - Add reversal of contig2 to FindConnectorSpectra
//  - Add IDX_START, IDX_END functionality for possible grid parallelization
//  - Add final step to select which pairs of contigs are connected by which connectors
//     1. Find best pair of contigs c1,c2, shift(c1,c2) + set of connectors C (with corresponding shifts)
//     2. Remove (c_i,c_j) pairs where c_i==c1 or c_j==c2. Remove c \in C from set of available connectors
//     3. Update connected contigs' scores to reflect removal of all connectors in C
//     4. Repeat from 1 until no contigs are significantly connected
//


int main(int argc, char **argv) {
	InputParams params; 	vector<char *> paramStrs;   paramStrs.resize(3);
	paramStrs[0] = "INPUT_CLUSTERS";
	paramStrs[1] = "INPUT_SPECS";
	paramStrs[2] = "OUTPUT_ALIGNS";
	
	if(argc==1)	params.readParams("mergecontigs.params"); else params.readParams(argv[1]);
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"mergecontigs.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_CLUSTERS, INPUT_SPECS, OUTPUT_ALIGNS\n";
		return -1;
	}

	char *clustersFN = params.getValue(paramStrs[0]);
	char *specSetFN = params.getValue(paramStrs[1]);
	char *alignsFN = params.getValue(paramStrs[2]);

	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0.4;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	int mergedPeakCount = params.paramPresent("MERGED_PEAK_COUNT")?params.getValueInt("MERGED_PEAK_COUNT"):2;

    Clusters contigs;
	cout << "Reading contigs... "; cerr.flush();
    contigs.Load(clustersFN);   if(contigs.size()<=0) { cerr<<"Error reading "<<clustersFN<<"!\n"; return(-1); }
	cout << "done: "<<contigs.size()<<" contigs read.\n";

	SpecSet specSet;
    specSet.LoadSpecSet_pkl(specSetFN);   if (specSet.size()==0) return -1;
    cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";
	
	list<Results_CS> results;
	list<vector<float> > ratios;
	findConnectorSpectra(contigs, specSet, peakTol, minRatio, 57, mergedPeakCount, results, ratios);

	// Save results
	Save_resultsCS(alignsFN, results, ';');

	// Save ratios
	if(params.paramPresent("OUTPUT_RATIOS"))
		Save_binArray(params.getValue("OUTPUT_RATIOS"), ratios);
		
	return 0;
}
