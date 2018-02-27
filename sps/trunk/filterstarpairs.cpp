#include "inputParams.h"
#include "alignment_scoring.h"
#include "spectral_pairs.h"
#include "batch.h"
#include "filters.h"
//#include "abruijn.h"
//#include "graph.h"
#include "SetMerger.h"
#include "SpectralPairs.h"
#include "SpectrumPair.h"
#include "SpectrumPairSet.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;
using namespace specnets;

int FilterPairs(SpecSet &specSet, SpectrumPairSet &aligns,
		SetMerger &components, InputParams &params, bool isASP);

int main(int argc, char **argv){
    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("filterstarpairs.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "filterstarpairs.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs;   paramStrs.resize(3);
	paramStrs[0] = "OUTPUT_ALIGNS";
	paramStrs[1] = "PENALTY_PTM";
	paramStrs[2] = "PENALTY_SAME_VERTEX";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"filterstarpairs.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: OUTPUT_ALIGNS, PENALTY_PTM, PENALTY_SAME_VERTEX\n";
		return -1;
	}

	SpecSet specSet;   short loadOk=0;
	if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
	else if(params.paramPresent("INPUT_SPECS_PKLBIN")) loadOk=specSet.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS_PKLBIN"));
	if (loadOk<=0 or specSet.size()==0) return -1;

	const char *alignsFN = params.getValue("INPUT_ALIGNS");
	bool alignPA = params.paramPresent("PARTIAL_OVERLAPS")?(bool) params.getValueInt("PARTIAL_OVERLAPS"):0;
	SpectrumPairSet alignsASP;    SpectrumPairSet alignsPA;

	if (alignPA) {
		if (alignsPA.loadFromBinaryFile(alignsFN)==0) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
	    else  cout << "Loading aligns complete. Number of input pairs: " << alignsPA.size() << "\n";
	} else {
		if (alignsASP.loadFromBinaryFile(alignsFN)==0) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
	    else  cout << "Loading aligns complete. Number of input pairs: " << alignsASP.size() << "\n";
	}

	// Separate aligns into isolated networks (connected components)
	SetMerger components(specSet.size());
	components.createSets(specSet.size(),2,alignsASP,alignsPA);

	if(alignPA) return FilterPairs(specSet,alignsPA,components,params,false);
	else return FilterPairs(specSet,alignsASP,components,params,true);
}

int FilterPairs(SpecSet &specSet, SpectrumPairSet &aligns,
		SetMerger &components, InputParams &params, bool isASP) {

	float penalty_ptm = (float) params.getValueDouble("PENALTY_PTM");
	float penalty_sameVert = (float) params.getValueDouble("PENALTY_SAME_VERTEX");
	int maxAAjump = params.paramPresent("MAX_AA_JUMP")?params.getValueInt("MAX_AA_JUMP"):0;
	float maxModMass = params.paramPresent("MAX_MOD_MASS")?(float) params.getValueDouble("MAX_MOD_MASS"):100.0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1.0;
	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):-1;
	unsigned int minMatchedPeaks = params.paramPresent("MIN_MATCHED_PEAKS")?(int) params.getValueInt("MIN_MATCHED_PEAKS"):0;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;

	vector<SpectrumPairSet> cAligns;    vector<vector<int> > cAligns_idx;  // To allow going back from components to original order
	components.splitAligns(aligns,cAligns,cAligns_idx);

	vector<vector<float> > ratios(aligns.size());
	for(unsigned int i=0;i<ratios.size();i++) {ratios[i].resize(3); ratios[i][0]=0; ratios[i][1]=0; ratios[i][2]=0;}

	vector<float> maxSpecScores(specSet.size());
	vector<bool> specFlipped(specSet.size());
	for(unsigned int i=0;i<specSet.size();i++) {
		maxSpecScores[i]=0; for(unsigned int j=0;j<specSet[i].size();j++) maxSpecScores[i]+=specSet[i][j][1];
		specSet[i].addZPMpeaks(peakTol,ionOffset,true);
		specFlipped[i]=false;
	}

	vector<vector<TwoValues<int> > > matches;
	vector<vector<float> > cRatios;
	vector<float> modPos;
	SpecSet matchedPeaks(2*aligns.size());
	for(unsigned int cIdx=0; cIdx<cAligns.size(); cIdx++) {
		SplitPairs(specSet, cAligns[cIdx], peakTol, pmTol, maxAAjump, maxModMass, penalty_sameVert, penalty_ptm, matches, specFlipped, modPos, 0, 0, false, false, &cRatios);

		// Copy alignment statistics to the corresponding global pair indices (as opposed to per-connected-component pair indices)
		for(unsigned int pairIdx=0; pairIdx<cAligns[cIdx].size(); pairIdx++) {
			for(unsigned int i=0; i<3; i++) ratios[cAligns_idx[cIdx][pairIdx]][i] = cRatios[pairIdx][i];

			unsigned int mpIdx = 2*cAligns_idx[cIdx][pairIdx];
			matchedPeaks[mpIdx].resize(matches[pairIdx].size());
			matchedPeaks[mpIdx+1].resize(matches[pairIdx].size());
			unsigned int s1=cAligns[cIdx][pairIdx].spec1, s2=cAligns[cIdx][pairIdx].spec2;
//if(s1==0 and s2==1) {
//	cerr<<" * pairIdx = "<<cAligns_idx[cIdx][pairIdx]<<", mpIdx = "<<mpIdx<<"\n";
//}
			for(unsigned int matchIdx=0; matchIdx<matches[pairIdx].size(); matchIdx++) {
				matchedPeaks[mpIdx][matchIdx] = specSet[s1][matches[pairIdx][matchIdx][0]];
				matchedPeaks[mpIdx+1][matchIdx] = specSet[s2][matches[pairIdx][matchIdx][1]];
//if(s1==0 and s2==1) {
//	cerr<<" * "<<matches[pairIdx][matchIdx][0]<<"/("<<matchedPeaks[mpIdx][matchIdx][0]<<","<<matchedPeaks[mpIdx][matchIdx][1]<<") - "
//		<<matches[pairIdx][matchIdx][1]<<"/("<<matchedPeaks[mpIdx+1][matchIdx][0]<<","<<matchedPeaks[mpIdx+1][matchIdx][1]<<")\n";
//}
			}
		}
/*		for(unsigned int pairIdx=0; pairIdx<cAligns[cIdx].size(); pairIdx++) {
			float score1=0, score2=0;   int s1=cAligns[cIdx][pairIdx].spec1, s2=cAligns[cIdx][pairIdx].spec2;
			for(unsigned int i=0; i<matches[pairIdx].size(); i++) {
				if(matches[pairIdx][i][0]>=0) score1+=specSet[s1][matches[pairIdx][i][0]][1];
				if(matches[pairIdx][i][1]>=0) score2+=specSet[s2][matches[pairIdx][i][1]][1];
			}
			ratios[cAligns_idx[cIdx][pairIdx]][0] = score1/maxSpecScores[s1];
			ratios[cAligns_idx[cIdx][pairIdx]][1] = score2/maxSpecScores[s2];
			ratios[cAligns_idx[cIdx][pairIdx]][2] = matches[pairIdx].size();
		}
*/
	}

	// Keep only the spectral pairs matching >= MIN_RATIO percent intensity in both aligned spectra
	unsigned int keptPairs=0;
	SpectrumPair curPair;
	for(unsigned int idxRatio=0; idxRatio<ratios.size(); idxRatio++) {
		curPair = aligns[idxRatio];
		unsigned int extraMatchPeaks=0;     // ASP pairs always match either/both PRMs at 0/18 or/and PM-19/PM-1 so these have higher numbers of peaks that need to match
		if(abs(curPair.shift1)<=peakTol) extraMatchPeaks++;  // matching 0/18 does not count towards achieving minMatchedPeaks
		if(abs(curPair.shift2)<=peakTol) extraMatchPeaks++;  // matching PM-19/PM-1 does not count towards achieving minMatchedPeaks

//if(aligns[idxRatio].spec1==9630 and aligns[idxRatio].spec2==9642) {
//	cerr<<" ----- aligns[idxRatio].shift1 = "<<aligns[idxRatio].shift1<<", peakTol = "<<peakTol<<", (abs(aligns[idxRatio].shift1)<=peakTol) = "<<(abs(aligns[idxRatio].shift1)<=peakTol)<<"\n";
//	cerr<<" ----- idxRatio = "<<idxRatio<<", keptPairs = "<<keptPairs<<", ratios[idxRatio][2] = "<<ratios[idxRatio][2]<<", minMatchedPeaks = "<<minMatchedPeaks<<", extraMatchPeaks = "<<extraMatchPeaks<<", pass = "<<(ratios[idxRatio][2]>=minMatchedPeaks+extraMatchPeaks)<<"\n";
//}
		if(ratios[idxRatio][0]>=minRatio and ratios[idxRatio][1]>=minRatio and ratios[idxRatio][2]>=minMatchedPeaks+extraMatchPeaks) {
			aligns[keptPairs] = aligns[idxRatio];
			for(unsigned int pivot=0;pivot<3;pivot++)
				ratios[keptPairs][pivot]=ratios[idxRatio][pivot];
			keptPairs++;
		}
	}
	aligns.resize(keptPairs);   ratios.resize(keptPairs);
	cout<<"Retained "<<keptPairs<<" pairs.\n";

//	Save_results_bin(params.getValue("OUTPUT_ALIGNS"), aligns.size(), aligns.begin());
	aligns.saveToBinaryFile(params.getValue("OUTPUT_ALIGNS"));
	if(params.paramPresent("OUTPUT_RATIOS")) Save_binArray(params.getValue("OUTPUT_RATIOS"), ratios);
	if(params.paramPresent("OUTPUT_SPECS")) matchedPeaks.SaveSpecSet_pklbin(params.getValue("OUTPUT_SPECS"));

	return 0;
}
