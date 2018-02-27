#include "inputParams.h"
#include "alignment_scoring.h"
#include "batch.h"
#include "hash.h"
#include "setmerger.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

int main(int argc, char **argv){
    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("batch_clst.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "batch_clst.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<char *> paramStrs;   paramStrs.resize(1);
//	paramStrs[0] = "INPUT_SPECS";
	paramStrs[0] = "OUTPUT_PAIRS";
	if(!params.confirmParams(paramStrs) or !(params.paramPresent("INPUT_SPECS") or params.paramPresent("INPUT_SPECS_PKLBIN"))) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"batch_asp.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_SPECS, OUTPUT_PAIRS\n";
		return -1;
	}
	
//	char *specSetFN = params.getValue(paramStrs[0]);
	char *pairsFN = params.getValue(paramStrs[0]);
	
	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float resolution = params.paramPresent("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
	unsigned int minClusterSize = params.paramPresent("MIN_CLUSTER_SIZE")?(unsigned int)params.getValueInt("MIN_CLUSTER_SIZE"):3;
	unsigned int specIdx, peakIdx;
    
    SpecSet specSet, specSetTopK;  short loadOk;
    if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
    else loadOk=specSet.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS_PKLBIN"));
    if (loadOk<=0 or specSet.size()==0) return -1;
    cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";
    cout << "Finding all clusters with at least " << minClusterSize << " spectra.\n";

	// Filter only top k peaks
	specSetTopK.resize(specSet.size());
	for(specIdx=0; specIdx<specSet.size(); specIdx++) specSet[specIdx].selectTopK(3,&specSetTopK[specIdx]);

    list<vector<unsigned int> > pairs;
	HashTable index;
	index.Build(specSetTopK, resolution, peakTol, 0);
	index.QuerySet(specSetTopK, pairs, 0, false);
	cout<<"Hashing returns "<<pairs.size()<<" pairs.\n";

	// Compute normalized dot products
	vector<TwoValues<float> > totIntensity(specSet.size());  // Total intensity (pos.0), Sqrt of the sum of squared intensities (pos.1)
	for(specIdx=0; specIdx<specSet.size(); specIdx++) {  
		totIntensity[specIdx].set(0,0);
		for(peakIdx=0; peakIdx<specSet[specIdx].size(); peakIdx++) {  // Compute norm to convert spectrum to a unit-length vector
			totIntensity[specIdx][0]+=specSet[specIdx][peakIdx][1];
			totIntensity[specIdx][1]+=(specSet[specIdx][peakIdx][1]*specSet[specIdx][peakIdx][1]);
		}
		totIntensity[specIdx][1] = sqrt(totIntensity[specIdx][1]);
		
		// smooth spectra here if it becomes a requirement to match 2+ peaks within peakTol Da windows
		
		for(peakIdx=0; peakIdx<specSet[specIdx].size(); peakIdx++)  // Convert peak intensities to log-values to enable ScoreOverlap6 (intensity sums instead of products)
			specSet[specIdx][peakIdx][1] = log10(specSet[specIdx][peakIdx][1]);
	}

	list<vector<unsigned int> >::iterator pairsIter=pairs.begin();
	vector<short> idx1(2048), idx2(2048), idxMatched1(2048), idxMatched2(2048);
	list<vector<float> > ratios;   vector<float> curRatio(2);
	float angleCos, intensity1, intensity2, totalMatched1, totalMatched2;
	while(pairsIter!=pairs.end()) {
		if(fabs(specSet[(*pairsIter)[0]].parentMass-specSet[(*pairsIter)[1]].parentMass)>pmTol)
			{ pairsIter = pairs.erase(pairsIter); continue; }
			
		FindMatchPeaksAll(specSet[(*pairsIter)[0]], specSet[(*pairsIter)[1]], 0, peakTol, idx1, idx2);
		ScoreOverlap6(specSet[(*pairsIter)[0]], idx1, specSet[(*pairsIter)[1]], idx2, 0, 
					 peakTol, idxMatched1, idxMatched2, 3*peakTol, NULL);

		angleCos=0; intensity1=0; intensity2=0; totalMatched1=0; totalMatched2=0;
		for(unsigned int i=0; i<idxMatched1.size(); i++) {
			intensity1 = pow(10,specSet[(*pairsIter)[0]][idxMatched1[i]][1]);  totalMatched1 += intensity1;
			intensity2 = pow(10,specSet[(*pairsIter)[1]][idxMatched2[i]][1]);  totalMatched2 += intensity2;
			angleCos+=intensity1*intensity2;
		}
		angleCos = angleCos/(totIntensity[(*pairsIter)[0]][1]*totIntensity[(*pairsIter)[1]][1]);
		curRatio[0] = totalMatched1/totIntensity[(*pairsIter)[0]][0];
		curRatio[1] = totalMatched2/totIntensity[(*pairsIter)[1]][0];
		
		if(angleCos>.8 and min(curRatio[0],curRatio[1])>=minRatio) { ratios.push_back(curRatio); pairsIter++; } else {
//cerr<<"("<<(*pairsIter)[0]<<","<<(*pairsIter)[1]<<"): angleCos = "<<angleCos<<", matched intensities = "<<pairsIter->score1<<"/"<<pairsIter->score2<<endl;
			pairsIter = pairs.erase(pairsIter);
		}
	}
	cout<<"Filtering retains "<<pairs.size()<<" pairs.\n";

	if(params.paramPresent("OUTPUT_CLUSTERS")) {  // Separate aligns into components
		SetMerger clusters(specSet.size());
		for(unsigned int i=0; i<specSet.size(); i++) clusters.createset(i);
		for(pairsIter=pairs.begin(); pairsIter!=pairs.end(); pairsIter++) 
			clusters.merge(clusters.membership[(*pairsIter)[0]],clusters.membership[(*pairsIter)[1]]);
		if(minClusterSize>1) clusters.removeSmallSets(minClusterSize);
		clusters.compressSetIndices();
		clusters.saveas_binListArray(params.getValue("OUTPUT_CLUSTERS"));

		// Remove pairs/ratios in removed clusters
		list<vector<float> >::iterator ratiosIter = ratios.begin();   pairsIter=pairs.begin();
		while(pairsIter!=pairs.end()) {
			if(clusters.membership[(*pairsIter)[0]]<0 or clusters.membership[(*pairsIter)[1]]<0)
				{ pairsIter = pairs.erase(pairsIter);   ratiosIter = ratios.erase(ratiosIter); }
			else { pairsIter++;   ratiosIter++; }
		}

		cout<<"Found "<<clusters.numSets<<" clusters with >="<<minClusterSize<<" spectra ("<<pairs.size()<<" pairs)\n";
	}
	    	
    // Output results
	Save_binArray<unsigned int>(pairsFN, pairs);
	if(params.paramPresent("OUTPUT_RATIOS")) Save_binArray<float>(params.getValue("OUTPUT_RATIOS"), ratios);
    
    return(0);
}
