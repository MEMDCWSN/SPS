#include "inputParams.h"
#include "batch.h"
#include "filters.h"

#include <fstream>
using namespace std;

int main(int argc, char **argv){
	InputParams params;
	
	if(argc==1)	params.readParams("grid_followup.params"); else params.readParams(argv[1]);
	unsigned int startIdx = params.paramPresent("IDX_START")?(unsigned int)params.getValueInt("IDX_START"):0;	
	unsigned int endIdx = params.paramPresent("IDX_END")?(unsigned int)params.getValueInt("IDX_END"):0;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0;
	char *indicesFN = params.paramPresent("OUTPUT_INDICES")?params.getValue("OUTPUT_INDICES"):(char *)0;

	vector<Results_ASP> aligns;   vector<unsigned int> selected;
	if (params.paramPresent("INPUT_ALIGNS")) {
		if(not Load_resultsASPbin(params.getValue("INPUT_ALIGNS"), aligns)) { cerr << "Error reading "<<params.getValue("INPUT_ALIGNS")<<"!\n"; return -1; }
	} else {
		if (params.paramPresent("INPUT_ALIGNS_MULTIPLE")) {
			if(not Load_resultsASPbin_multiple(params.getValue("INPUT_ALIGNS_MULTIPLE"), aligns)) { cerr << "Error reading "<<params.getValue("INPUT_ALIGNS_MULTIPLE")<<"!\n"; return -1; }
		} 
	}
    cout << "Loading aligns complete. Num pairs: " << aligns.size() << "\n"; cout.flush();
    if(not params.paramPresent("IDX_END")) endIdx = aligns.size();

	vector<TwoValues<float> > gaussianParams;
	vector<unsigned int> idxKept(0);
	unsigned int idxPair, idxLast;
	if(params.paramPresent("INPUT_MEANS_MULTIPLE") and params.paramPresent("INPUT_VARS_MULTIPLE")) {
		// Get means/stddevs from multiple mean and variance terms
		if(not LoadGaussianParams(params.getValue("INPUT_MEANS_MULTIPLE"), params.getValue("INPUT_VARS_MULTIPLE"), gaussianParams))
			{ cerr<<"Error loading p-value distributions from "<<params.getValue("INPUT_MEANS_MULTIPLE")<<" / "<<params.getValue("INPUT_VARS_MULTIPLE")<<"!\n"; return -1; }
		if(params.paramPresent("OUTPUT_MEANS_STDDEVS")) Save_binArray(params.getValue("OUTPUT_MEANS_STDDEVS"),gaussianParams);
		cout<<"Done estimating means/stddevs\n";
		
		// Compute p-values and filter aligns (if any alignments were input)
		if(aligns.size()>0) {
			vector<TwoValues<float> > pvalues(aligns.size());
			for(idxPair=0; idxPair<aligns.size(); idxPair++) {
				pvalues[idxPair][0] = 1 - Utils::gaussiancdf(aligns[idxPair].score1,gaussianParams[aligns[idxPair].spec1][0],gaussianParams[aligns[idxPair].spec1][1]);
				pvalues[idxPair][1] = 1 - Utils::gaussiancdf(aligns[idxPair].score2,gaussianParams[aligns[idxPair].spec2][0],gaussianParams[aligns[idxPair].spec2][1]);
			}
			cout<<"Done estimating p-values\n";

			if(params.paramPresent("OUTPUT_PVALUES")) Save_binArray(params.getValue("OUTPUT_PVALUES"),pvalues);
			if(params.paramPresent("MAX_PVALUE")) { // and params.paramPresent("OUTPUT_ALIGNS")) {
				idxLast=0;   idxKept.resize(aligns.size());   float minPValue = params.getValueDouble("MAX_PVALUE");
				for(idxPair=0; idxPair<aligns.size(); idxPair++)
					if(pvalues[idxPair][0]<=minPValue and pvalues[idxPair][1]<=minPValue)
						{ if(idxLast<idxPair) aligns[idxLast]=aligns[idxPair]; idxKept[idxLast]=idxPair; idxLast++; }
				aligns.resize(idxLast);   idxKept.resize(idxLast);
				cout<<"Got "<<idxLast<<" pairs with both pvalues <= "<<minPValue<<endl;
//			    Save_resultsASPbin(params.getValue("OUTPUT_ALIGNS"), aligns);
			}
		}
	}

	// Load and filter ratios
	vector<vector<float> > ratios;
	if(params.paramPresent("INPUT_RATIOS_MULTIPLE")) {
		if(not Load_binArray_multiple(params.getValue("INPUT_RATIOS_MULTIPLE"), ratios))
			{ cerr << "Error reading "<<params.getValue("INPUT_RATIOS_MULTIPLE")<<"!\n"; return -1; }
		if(ratios.size()>0 and ratios[0].size()<2) { cerr << "Error reading "<<params.getValue("INPUT_RATIOS_MULTIPLE")<<" (invalid format)!\n"; return -1; }
		if(params.paramPresent("OUTPUT_RATIOS")) Save_binArray(params.getValue("OUTPUT_RATIOS"),ratios);

		if(minRatio>0) {
			if(idxKept.size()==0)
				{ idxKept.resize(aligns.size()); for(idxPair=0;idxPair<idxKept.size();idxPair++) idxKept[idxPair]=idxPair; }
			idxLast=0;
			for(idxPair=0;idxPair<idxKept.size();idxPair++) 
				if(ratios[idxKept[idxPair]][0]>=minRatio and ratios[idxKept[idxPair]][1]>=minRatio) 
					{ aligns[idxLast]=aligns[idxKept[idxPair]]; idxKept[idxLast]=idxKept[idxPair]; idxLast++; }
		}
		aligns.resize(idxLast);   idxKept.resize(idxLast);
		cout<<"Got "<<idxLast<<" pairs with both ratios >= "<<minRatio<<endl;
	}
	
	if(params.paramPresent("OUTPUT_INDICES")) Save_binArray(params.getValue("OUTPUT_INDICES"),idxKept);
	if(params.paramPresent("OUTPUT_ALIGNS")) Save_resultsASPbin(params.getValue("OUTPUT_ALIGNS"), aligns);
/*
	FilterTriangles(aligns, startIdx, endIdx, pmTol, selected);

	cout << "Retained "<<aligns.size()<<" alignments. Saving output files..."; cout.flush();
	Save_resultsASPbin(outputFN, aligns);
	if(indicesFN) Save_binArray(indicesFN, selected);
	cout << "done.\n";
*/
	return 0;
}

/*
int main(int argc, char **argv){
	InputParams params;
	vector<char *> paramStrs;   paramStrs.resize(2);
	
	paramStrs[0] = "INPUT_ALIGNS";
	paramStrs[1] = "OUTPUT_ALIGNS";
	if(argc==1)	params.readParams("grid_followup.params"); else params.readParams(argv[1]);
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"grid_followup.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_ALIGNS, OUTPUT_ALIGNS\n";
		return -1;
	}
	char *inputFN = params.getValue(paramStrs[0]);
	char *outputFN = params.getValue(paramStrs[1]);
	unsigned int startIdx = params.paramPresent("IDX_START")?(unsigned int)params.getValueInt("IDX_START"):0;	
	unsigned int endIdx = params.paramPresent("IDX_END")?(unsigned int)params.getValueInt("IDX_END"):0;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	char *indicesFN = params.paramPresent("OUTPUT_INDICES")?params.getValue("OUTPUT_INDICES"):(char *)0;

	vector<float> means, vars;
	
	
	vector<Results_ASP> aligns;   vector<unsigned int> selected;
	if (params.paramPresent("INPUT_ALIGNS")) {
		if(not Load_resultsASP(params.getValue("INPUT_ALIGNS"), aligns)) { cerr << "Error reading "<<params.getValue("INPUT_ALIGNS")<<"!\n"; return -1; }
	} else {
		if (params.paramPresent("INPUT_ALIGNS_MULTIPLE")) {
			if(not Load_resultsASP_multiple(params.getValue("INPUT_ALIGNS_MULTIPLE"), aligns)) { cerr << "Error reading "<<params.getValue("INPUT_ALIGNS_MULTIPLE")<<"!\n"; return -1; }
		} else { cerr << "ERROR: Parameters file is incomplete - must specify either INPUT_ALIGNS or INPUT_ALIGNS_MULTIPLE\n"; return -1; }
	}
    cout << "Loading aligns complete. Num pairs: " << aligns.size() << "\n"; cout.flush();
    if(not params.paramPresent("IDX_END")) endIdx = aligns.size();

	FilterTriangles(aligns, startIdx, endIdx, pmTol, selected);

	cout << "Retained "<<aligns.size()<<" alignments. Saving output files..."; cout.flush();
	Save_resultsASPbin(outputFN, aligns);
	if(indicesFN) Save_binArray(indicesFN, selected);
	cout << "done.\n";
	
	return 0;
}
*/
