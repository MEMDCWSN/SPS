#include "inputParams.h"
#include "alignment_scoring.h"
#include "batch.h"
#include "filters.h"
#include "PepnovoTags.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;
using namespace specnets;

int main(int argc, char **argv){
    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("filterpairs.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "filterpairs.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs;   paramStrs.resize(1);
	paramStrs[0] = "OUTPUT_ALIGNS";
	if(!params.confirmParams(paramStrs) or !(params.paramPresent("INPUT_SPECS") or params.paramPresent("INPUT_SPECS_PKLBIN"))) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"filterpairs.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: OUTPUT_ALIGNS\n";
		return -1;
	}

	const char *alignsFN = params.getValue("OUTPUT_ALIGNS");
	int startBaseIdx; if(params.paramPresent("IDX_START")) startBaseIdx = max(0,params.getValueInt("IDX_START"));  else startBaseIdx=0;
	int endBaseIdx;   if(params.paramPresent("IDX_END")) endBaseIdx = max(0,params.getValueInt("IDX_END")); else endBaseIdx=-1;

	int aaDiff;   double minShift, maxShift;   float minOverlap=0;   short minNumMatchedPeaks=0;
	if (not params.paramPresent("AA_DIFF_COUNT") or not params.paramPresent("MIN_SHIFT") or not params.paramPresent("MAX_SHIFT"))
	    { cerr << "ERROR: Parameters file is incomplete. One of the following is missing: AA_DIFF_COUNT, MIN_SHIFT, MAX_SHIFT\n"; return -1; }
	aaDiff = params.getValueInt("AA_DIFF_COUNT");
	minShift = params.getValueDouble("MIN_SHIFT");
	maxShift = params.getValueDouble("MAX_SHIFT");
	bool alignPA = params.paramPresent("PARTIAL_OVERLAPS")?(bool) params.getValueInt("PARTIAL_OVERLAPS"):0;

	if(alignPA) {
		minOverlap = params.paramPresent("MIN_OVERLAP_AREA")?(float) params.getValueDouble("MIN_OVERLAP_AREA"):0;
		minNumMatchedPeaks = params.paramPresent("MIN_NUM_MATCHED_PEAKS")?(short) params.getValueInt("MIN_NUM_MATCHED_PEAKS"):0;
	}

	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float resolution = params.paramPresent("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
	int useMinDist = params.paramPresent("USE_MIN_DIST_57")?(int) params.getValueInt("USE_MIN_DIST_57"):1;
    bool specTypeMSMS = params.paramPresent("SPEC_TYPE_MSMS")?(bool)params.getValueInt("SPEC_TYPE_MSMS"):false;
	float symmetryOffset = specTypeMSMS?2*AAJumps::massHion:0;
	float ionOffset = specTypeMSMS?AAJumps::massHion:0;

	vector<vector<TTag> > tags(0);
	int tagsMatchFlank = params.paramPresent("TAGS_MATCH_FLANK")?(int) params.getValueInt("TAGS_MATCH_FLANK"):1;
	int tagsMatchCount = params.paramPresent("TAGS_MATCH_COUNT")?(int) params.getValueInt("TAGS_MATCH_COUNT"):1;
	if(params.paramPresent("INPUT_TAGS"))
		for(int i=0;i<=3;i++) {
			stringstream s; s<<i;
			LoadTags(string(params.getValue("INPUT_TAGS"))+s.str()+".txt", tags, tagsMatchFlank);
		}

    SpecSet specSet;   short loadOk;
    if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
    else loadOk=specSet.loadPklBin(params.getValue("INPUT_SPECS_PKLBIN"));
    if (loadOk<=0 or specSet.size()==0) return -1;

    cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";
	if(endBaseIdx<0) endBaseIdx=specSet.size()-1;
	if(startBaseIdx<0 or startBaseIdx>=specSet.size()) { cerr<<"Invalid start index "<<startBaseIdx<<" ("<<specSet.size()<<" spectra)\n"; return -1; }

	if (endBaseIdx>=specSet.size()) cout << "Warning: IDX_END ("<<endBaseIdx<<") was lowered to the index of the last spectrum in the dataset ("<<specSet.size()-1<<")\n";
	endBaseIdx = min(specSet.size()-1,(unsigned int)endBaseIdx);
	cout<<"Computing alignments between spectra "<<startBaseIdx<<":"<<endBaseIdx<<" and all others in the dataset.\n";
    vector<int> baseSpecIdx(endBaseIdx-startBaseIdx+1);  for(int i=0; i<endBaseIdx-startBaseIdx+1; i++) baseSpecIdx[i]=startBaseIdx+i;

    list<Results_ASP> resultsASP;    list<Results_PA> resultsPA;
    list<TwoValues<float> > ratios;
	vector<TwoValues<float> > means;
	vector<float> varTerms;
	if(not alignPA) {
		getPairAlignsASP(specSet, baseSpecIdx, aaDiff, minShift, maxShift, pmTol, peakTol, minRatio, resultsASP, ratios, means, varTerms, tags, tagsMatchFlank, tagsMatchCount, resolution, symmetryOffset);
	    Save_results_bin(alignsFN,resultsASP.size(),resultsASP.begin());
	} else {
		list<vector<float> > alignStats;    alignStats.clear();
		vector<vector<float> > specStats;   specStats.resize(0);
		list<TwoValues<int> > numMatchedPeaks;
		AAJumps jumps(0);   //jumps.alljumps(1000.0,resolution);

		for(unsigned int i=0; i<specSet.size(); i++) specSet[i].addZPMpeaks(peakTol,ionOffset,true);
		getPairAlignsPA2(specSet, startBaseIdx, endBaseIdx, peakTol, pmTol, minRatio, minOverlap, minNumMatchedPeaks, jumps, minShift, resultsPA, ratios, numMatchedPeaks, means, varTerms, alignStats, specStats);
	    Save_results_bin(alignsFN,resultsPA.size(),resultsPA.begin());

	    const char *filename = params.getValue("OUTPUT_STATS");
		if(params.paramPresent("OUTPUT_STATS")) {
			char *fnBuffer = (char *)malloc(strlen(filename)+20);
			sprintf(fnBuffer,"%s_specStats.bin",filename);
			Save_binArray(fnBuffer, specStats);
			sprintf(fnBuffer,"%s_alignStats.bin",filename);
			Save_binArray(fnBuffer, alignStats);
			free(fnBuffer);
		}
	}

	if(params.paramPresent("OUTPUT_RATIOS")) {
	    cout << "Outputting ratios..."; cout.flush();
	    Save_binArray(params.getValue("OUTPUT_RATIOS"), ratios);
	    cout << "done.\n"; cout.flush();
	}

	if(params.paramPresent("OUTPUT_MEANS")) {
	    cout << "Outputting means..."; cout.flush();
	    Save_binArray(params.getValue("OUTPUT_MEANS"), means);
	    cout << "done.\n"; cout.flush();
	}

	if(params.paramPresent("OUTPUT_VARIANCE")) {
	    cout << "Outputting variance terms..."; cout.flush();
	    Save_binArray(params.getValue("OUTPUT_VARIANCE"), varTerms);
	    cout << "done.\n"; cout.flush();
	}

    return(0);
}
