#include "inputParams.h"
#include "batch.h"
#include "filters.h"

#include <fstream>
using namespace std;
using namespace specnets;

int main(int argc, char **argv){
	InputParams params;
	
	if(argc==1)	params.readParams("mergefilter.params"); else params.readParams(argv[1]);
	unsigned int startIdx = params.paramPresent("IDX_START")?(unsigned int)params.getValueInt("IDX_START"):0;	
	unsigned int endIdx = params.paramPresent("IDX_END")?(unsigned int)params.getValueInt("IDX_END"):0;
	bool alignPA = params.paramPresent("PARTIAL_OVERLAPS")?(bool) params.getValueInt("PARTIAL_OVERLAPS"):0;
	bool filterTrigs = params.paramPresent("FILTER_TRIGS")?(bool) params.getValueInt("FILTER_TRIGS"):0;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):-1;
	const char *indicesFN = params.paramPresent("OUTPUT_INDICES")?params.getValue("OUTPUT_INDICES"):(char *)0;
	float minPValue = params.paramPresent("MAX_PVALUE")?params.getValueDouble("MAX_PVALUE"):1.1;
	bool doFilter=(filterTrigs==true or minRatio>0 or minPValue<=1.01);
	if(not doFilter) cout<<"Skipping pairs filters, just merging grid files...\n";
	
	vector<Results_ASP> alignsASP;   vector<Results_PA> alignsPA;
	if (params.paramPresent("INPUT_ALIGNS")) {
		if(alignPA)	{ if(not Load_results_bin(params.getValue("INPUT_ALIGNS"), alignsPA)) { cerr << "Error reading "<<params.getValue("INPUT_ALIGNS")<<"!\n"; return -1; } }
		else { if(not Load_results_bin(params.getValue("INPUT_ALIGNS"), alignsASP)) { cerr << "Error reading "<<params.getValue("INPUT_ALIGNS")<<"!\n"; return -1; } }
	} else {
		if (params.paramPresent("INPUT_ALIGNS_MULTIPLE")) {
			if(alignPA)	{ if(not Load_results_bin_multiple(params.getValue("INPUT_ALIGNS_MULTIPLE"), alignsPA)) { cerr << "Error reading "<<params.getValue("INPUT_ALIGNS_MULTIPLE")<<"!\n"; return -1; } }
			else { if(not Load_results_bin_multiple(params.getValue("INPUT_ALIGNS_MULTIPLE"), alignsASP)) { cerr << "Error reading "<<params.getValue("INPUT_ALIGNS_MULTIPLE")<<"!\n"; return -1; } }
		} 
	}
	unsigned int numPairs = alignPA?alignsPA.size():alignsASP.size();
    cout << "Loading aligns complete. Num pairs: " << numPairs << "\n"; cout.flush();
    if(numPairs<=0) { cerr<<"ERROR: Cannot merge/filter with "<<numPairs<<" input pairs!\n"; return -1; }
    if(not params.paramPresent("IDX_END")) endIdx = numPairs;
	if(not doFilter) {
		if(alignPA) { if(params.paramPresent("OUTPUT_ALIGNS")) Save_results_bin(params.getValue("OUTPUT_ALIGNS"),alignsPA.size(),alignsPA.begin()); alignsPA.resize(0); }
		else { if(params.paramPresent("OUTPUT_ALIGNS")) Save_results_bin(params.getValue("OUTPUT_ALIGNS"),alignsASP.size(),alignsASP.begin()); alignsASP.resize(0); }
	}

	vector<TwoValues<float> > gaussianParams;
	unsigned int idxPair, idxLast;
	if(params.paramPresent("INPUT_MEANS_MULTIPLE") and params.paramPresent("INPUT_VARS_MULTIPLE")) {
		// Get means/stddevs from multiple mean and variance terms
		if(not LoadGaussianParams(params.getValue("INPUT_MEANS_MULTIPLE"), params.getValue("INPUT_VARS_MULTIPLE"), gaussianParams))
			{ cerr<<"Error loading score distribution's parameters from "<<params.getValue("INPUT_MEANS_MULTIPLE")<<" / "<<params.getValue("INPUT_VARS_MULTIPLE")<<"!\n"; return -1; }
		if(params.paramPresent("OUTPUT_MEANS_STDDEVS")) Save_binArray(params.getValue("OUTPUT_MEANS_STDDEVS"),gaussianParams);
		cout<<"Done estimating means/stddevs\n";
	}
	if(not doFilter) gaussianParams.resize(0);
		
	// Load and filter ratios
	vector<vector<float> > ratios;
	if(params.paramPresent("INPUT_RATIOS_MULTIPLE")) {
		if(not Load_binArray_multiple(params.getValue("INPUT_RATIOS_MULTIPLE"), ratios))
			{ cerr << "Error reading "<<params.getValue("INPUT_RATIOS_MULTIPLE")<<"!\n"; return -1; }
		if(ratios.size()>0 and ratios[0].size()<2) { cerr << "Error reading "<<params.getValue("INPUT_RATIOS_MULTIPLE")<<" (invalid format)!\n"; return -1; }
		if(params.paramPresent("OUTPUT_RATIOS")) Save_binArray(params.getValue("OUTPUT_RATIOS"),ratios);
		if(minRatio<=0) ratios.resize(0); 
	}
	cout<<"Done merging ratios\n";
	
	vector<unsigned int> idxKept;
	if(not doFilter) {
		if(params.paramPresent("OUTPUT_INDICES")) {
			ratios.resize(0);   idxKept.resize(numPairs);
			for(idxPair=0; idxPair<numPairs; idxPair++) idxKept[idxPair]=idxPair;
			Save_binArray(params.getValue("OUTPUT_INDICES"),idxKept);
		}
		return 0;
	}
	
	vector<TwoValues<float> > pvalues;
	if(alignPA) {
		FilterAligns(alignsPA,idxKept,pvalues,gaussianParams,ratios,minPValue,minRatio,pmTol,filterTrigs);
		if(params.paramPresent("OUTPUT_ALIGNS")) Save_results_bin(params.getValue("OUTPUT_ALIGNS"),alignsPA.size(),alignsPA.begin());
	} else {
		FilterAligns(alignsASP,idxKept,pvalues,gaussianParams,ratios,minPValue,minRatio,pmTol,filterTrigs);
		if(params.paramPresent("OUTPUT_ALIGNS")) Save_results_bin(params.getValue("OUTPUT_ALIGNS"),alignsASP.size(),alignsASP.begin());
	}
	if(params.paramPresent("OUTPUT_PVALUES")) Save_binArray(params.getValue("OUTPUT_PVALUES"),pvalues);
	if(params.paramPresent("OUTPUT_INDICES")) Save_binArray(params.getValue("OUTPUT_INDICES"),idxKept);
	return 0;
}

