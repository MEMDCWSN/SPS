#include "inputParams.h"
#include "alignment_scoring.h"
#include "batch.h"
#include "filters.h"
#include "dekel_align.h"
#include "PepnovoTags.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char **argv){
    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("batch_asp.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "batch_asp.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<char *> paramStrs;   paramStrs.resize(4);
	paramStrs[0] = "OUTPUT_ALIGNS";
	paramStrs[1] = "AA_DIFF_COUNT";
	paramStrs[2] = "MIN_SHIFT";
	paramStrs[3] = "MAX_SHIFT";
	if(!params.confirmParams(paramStrs) or !(params.paramPresent("INPUT_SPECS") or params.paramPresent("INPUT_SPECS_PKLBIN"))) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"batch_asp.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_SPECS, OUTPUT_ALIGNS, AA_DIFF_COUNT, MIN_SHIFT, MAX_SHIFT\n";
		return -1;
	}
	
	char *alignsFN = params.getValue("OUTPUT_ALIGNS");
	int startBaseIdx; if(params.paramPresent("IDX_START")) startBaseIdx = max(0,params.getValueInt("IDX_START"));  else startBaseIdx=0;
	int endBaseIdx;   if(params.paramPresent("IDX_END")) endBaseIdx = max(0,params.getValueInt("IDX_END")); else endBaseIdx=-1;
	int aaDiff = params.getValueInt("AA_DIFF_COUNT");
	double minShift = params.getValueDouble("MIN_SHIFT");
	double maxShift = params.getValueDouble("MAX_SHIFT");
	
	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float resolution = params.paramPresent("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
	int useMinDist = params.paramPresent("USE_MIN_DIST_57")?(int) params.getValueInt("USE_MIN_DIST_57"):1;
    float symmetryOffset = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?2:0):0;

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
    else loadOk=specSet.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS_PKLBIN"));
    if (loadOk<=0 or specSet.size()==0) return -1;
    
    cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";
	if(endBaseIdx<0) endBaseIdx=specSet.size()-1;
	if(startBaseIdx<0 or startBaseIdx>=specSet.size()) { cerr<<"Invalid start index "<<startBaseIdx<<" ("<<specSet.size()<<" spectra)\n"; return -1; }

	if (endBaseIdx>=specSet.size()) cout << "Warning: IDX_END ("<<endBaseIdx<<") was lowered to the index of the last spectrum in the dataset ("<<specSet.size()-1<<")\n";
	endBaseIdx = min(specSet.size()-1,(unsigned int)endBaseIdx);
	cout<<"Computing alignments between spectra "<<startBaseIdx<<":"<<endBaseIdx<<" and all others in the dataset.\n";
    vector<int> baseSpecIdx(endBaseIdx-startBaseIdx+1);  for(int i=0; i<endBaseIdx-startBaseIdx+1; i++) baseSpecIdx[i]=startBaseIdx+i;

    list<Results_ASP> results;
    list<TwoValues<float> > ratios;   list<TwoValues<float> > pvalues;   TwoValues<float> curPvalue;
	vector<TwoValues<float> > means;
	vector<float> varTerms;
    getPairAlignsASP(specSet, baseSpecIdx, aaDiff, minShift, maxShift, pmTol, peakTol, minRatio, results, ratios, means, varTerms, tags, tagsMatchFlank, tagsMatchCount, resolution, symmetryOffset);

	// Filter by p-value
	if(params.paramPresent("MAX_PVALUE")) {
	    cout << "Filtering by p-value <= "<<params.getValueDouble("MAX_PVALUE")<<"..."; cout.flush();
		float revPvalue = 1-(float)params.getValueDouble("MAX_PVALUE");
		vector<float> stddevs(varTerms.size());
		for(unsigned int specIdx=0; specIdx<stddevs.size(); specIdx++)
			stddevs[specIdx]= sqrt( (means[specIdx][0]/(means[specIdx][0]-1)) * (varTerms[specIdx]-means[specIdx][0]*means[specIdx][0]) );  // Unbiased Gaussian standard deviation
	
	    list<Results_ASP>::iterator iterResults = results.begin();
		list<TwoValues<float> >::iterator iterRatios = ratios.begin();
		while(iterResults != results.end()) {
		        curPvalue[0] = Utils::gaussiancdf(iterResults->score1,means[iterResults->spec1][0],stddevs[iterResults->spec1]);
		        curPvalue[1] = Utils::gaussiancdf(iterResults->score2,means[iterResults->spec2][0],stddevs[iterResults->spec2]);
			if(curPvalue[0]>=revPvalue and curPvalue[1]>=revPvalue)
				{ iterResults++; iterRatios++; curPvalue[0]=1-curPvalue[0]; curPvalue[1]=1-curPvalue[1]; pvalues.push_back(curPvalue); }
			else { iterResults=results.erase(iterResults); iterRatios=ratios.erase(iterRatios); }
		}
		cout << "done. Retained "<<results.size()<<" spectral pairs.\n"; cout.flush();
	}
    
    // Output results
    cout << "Outputing results..."; cout.flush();
    Save_resultsASPbin(alignsFN, results);
    cout << "done.\n"; cout.flush();
    
	if(params.paramPresent("OUTPUT_RATIOS")) {
	    cout << "Outputing ratios..."; cout.flush();
for(list<TwoValues<float> >::iterator it=ratios.begin(); it!=ratios.end(); it++) cerr<<"("<<(*it)[0]<<","<<(*it)[1]<<")"; cerr<<endl;
	    Save_binArray(params.getValue("OUTPUT_RATIOS"), ratios);
	    cout << "done.\n"; cout.flush();
	}

	if(params.paramPresent("OUTPUT_PVALUES")) {
	    cout << "Outputing pvalues..."; cout.flush();
	    Save_binArray(params.getValue("OUTPUT_PVALUES"), pvalues);
	    cout << "done.\n"; cout.flush();
	}

	if(params.paramPresent("OUTPUT_MEANS")) {
	    cout << "Outputing means..."; cout.flush();
	    Save_binArray(params.getValue("OUTPUT_MEANS"), means);
	    cout << "done.\n"; cout.flush();
	}

	if(params.paramPresent("OUTPUT_VARIANCE")) {
	    cout << "Outputing variance terms..."; cout.flush();
	    Save_binArray(params.getValue("OUTPUT_VARIANCE"), varTerms);
	    cout << "done.\n"; cout.flush();
	}

    return(0);
}
