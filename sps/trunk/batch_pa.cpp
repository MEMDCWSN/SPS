#include "alignment_scoring.h"
#include "batch.h"
#include "filters.h"
#include "clusters.h"
#include "inputParams.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

int main(int argc, char **argv){
	const short MIN_AA_DIST = 0;
	InputParams params;
	vector<char *> paramStrs;   paramStrs.resize(1);

	paramStrs[0] = "OUTPUT_ALIGNS";
	if(argc==1)	params.readParams("batch_pa.params"); else params.readParams(argv[1]);
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"batch_pa.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: OUTPUT_ALIGNS\n";
		return -1;
	}

	char *clustersFN = params.getValue("INPUT_CLUSTERS");  bool procClusters=params.paramPresent("INPUT_CLUSTERS");
	char *alignsFN = params.getValue(paramStrs[0]);

	int startIdx = params.paramPresent("IDX_START")?params.getValueInt("IDX_START"):0;
	int endIdx = params.paramPresent("IDX_END")?params.getValueInt("IDX_END"):-1;
	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0.3;
	float minOverlap = params.paramPresent("MIN_OVERLAP_AREA")?(float) params.getValueDouble("MIN_OVERLAP_AREA"):0;
	short minNumMatchedPeaks = params.paramPresent("MIN_NUM_MATCHED_PEAKS")?(short) params.getValueInt("MIN_NUM_MATCHED_PEAKS"):0;
	float minAbsShift = params.paramPresent("MIN_SHIFT")?(float) params.getValueDouble("MIN_SHIFT"):0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float resolution = params.paramPresent("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;
	cout << "Minimum match ratio is " << minRatio << "\n";

    SpecSet specSet;   short loadOk=0;
    if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
    else loadOk=specSet.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS_PKLBIN"));
    if (loadOk<=0 or specSet.size()==0) { cerr<<"ERROR reading input spectra!\n"; return -1; }
	for(unsigned int i=0; i<specSet.size(); i++) specSet[i].addZPMpeaks(peakTol,ionOffset,true);

    if(endIdx<0) endIdx = specSet.size()-1;
    cout << "Loading specs complete. Num specs: " << specSet.size() <<", computing pairwise aligns for indices "<<startIdx<<" to "<<endIdx<<"\n";

//	AAJumps jumps3(6,resolution);   if(peakTol>0.1) jumps3.removeHigherJumps(250.0); else jumps3.removeHigherJumps(500.0);
	AAJumps jumps(0);   //jumps.alljumps(1000.0,resolution);
	cout<<"Enforcing amino acid mass jumps up to 1000.0 Da ("<<jumps.size()<<" different valid jumps)\n";

	list<Results_PA> results;
	list<TwoValues<float> > ratios;
	list<TwoValues<int> > numMatchedPeaks;
	vector<TwoValues<float> > meansClst, meansSpecs;
	vector<float> varTerms;
	list<vector<float> > alignStats;    alignStats.clear();
	vector<vector<float> > specStats;   specStats.resize(0);
	if(procClusters) {
cerr<<"Warning: MIN_AA_DIST hard coded to zero!\n";
		specSet.addZPMpeaks(0.5,ionOffset,false);   specSet.setResolution(resolution,true);

	    Clusters clusters;
	    clusters.Load(clustersFN);
		clusters.setResolution(resolution,true);

		getPairAlignsPAext(clusters,specSet,startIdx,endIdx,MIN_AA_DIST,pmTol/resolution,peakTol/resolution,minRatio,results,ratios, meansClst, meansSpecs);
//list<TwoValues<float> >::iterator r = ratios.begin();
//for(list<Results_PA>::const_iterator i=results.begin(); i!=results.end(); i++, r++)
//	cerr<<"results for specs ("<<(*i).spec1<<","<<(*i).spec2<<") has shifts ["<<(*i).shift1*resolution<<","<<(*i).shift2*resolution<<"], scores ("<<(*i).score1<<","<<(*i).score2<<"), ratios ("<<(*r)[0]<<","<<(*r)[1]<<")\n";
	} else {
		cout<<"Computing alignments between spectra "<<startIdx<<":"<<endIdx<<" and all others in the dataset.\n";
		getPairAlignsPA2(specSet, startIdx, endIdx, peakTol, pmTol, minRatio, minOverlap, minNumMatchedPeaks, jumps, minAbsShift, results, ratios, numMatchedPeaks, meansSpecs, varTerms, alignStats, specStats);
	}

    ofstream output;
    output.open(alignsFN, ios::binary);  output << results.size() << endl;
	for(list<Results_PA>::iterator i=results.begin(); i!=results.end(); i++)
		if(procClusters) { (*i).shift1*=resolution; (*i).shift2*=resolution; (*i).output(output,';'); }
		else (*i).output(output,';');
    output.close();

	char *filename = params.getValue("OUTPUT_RATIOS");
	if(params.paramPresent("OUTPUT_RATIOS")) {
	    output.open(filename, ios::binary);
	    for(list<TwoValues<float> >::iterator r = ratios.begin(); r!=ratios.end(); r++) output<<(*r)[0]<<" "<<(*r)[1]<<endl;
	    output.close();
	}

	filename = params.getValue("OUTPUT_NUM_MATCHED_PEAKS");
	if(params.paramPresent("OUTPUT_NUM_MATCHED_PEAKS")) {
	    output.open(filename, ios::binary);
	    for(list<TwoValues<int> >::iterator r = numMatchedPeaks.begin(); r!=numMatchedPeaks.end(); r++) output<<(*r)[0]<<" "<<(*r)[1]<<endl;
	    output.close();
	}

	if(params.paramPresent("OUTPUT_MATCHED_PEAKS")) {
		SpecSet matchedPeaks(4*results.size());  // Each pair generates 4 sets of peaks
		                                         //   Peaks matched in spectra 1/2 with shift1
		                                         //   Peaks matched in spectra 1/2 with shift2
		unsigned int indexMP=0;
		unsigned int maxNumPeaks=0; for(unsigned int i=0; i<specSet.size(); i++) maxNumPeaks=max(maxNumPeaks,specSet[i].size());
		vector<int> idx1(maxNumPeaks), idx2(maxNumPeaks);

		for(list<Results_PA>::iterator curPair=results.begin(); curPair!=results.end(); curPair++) {
			FindMatchPeaks(specSet[curPair->spec1], specSet[curPair->spec2], curPair->shift1, peakTol, idx1, idx2);
			matchedPeaks[indexMP] = specSet[curPair->spec1];    matchedPeaks[indexMP++].selectIndices(idx1);
			matchedPeaks[indexMP] = specSet[curPair->spec2];    matchedPeaks[indexMP++].selectIndices(idx2);
			FindMatchPeaks(specSet[curPair->spec1], specSet[curPair->spec2], curPair->shift2, peakTol, idx1, idx2);
			matchedPeaks[indexMP] = specSet[curPair->spec1];    matchedPeaks[indexMP++].selectIndices(idx1);
			matchedPeaks[indexMP] = specSet[curPair->spec2];    matchedPeaks[indexMP++].selectIndices(idx2);
		}
//		matchedPeaks.SaveSpecSet_pklbin(params.getValue("OUTPUT_MATCHED_PEAKS"));
		matchedPeaks.SaveSpecSet_pkl(params.getValue("OUTPUT_MATCHED_PEAKS"));
	}

	filename = params.getValue("OUTPUT_STATS");
	if(params.paramPresent("OUTPUT_STATS")) {
		char *fnBuffer = (char *)malloc(strlen(filename)+20);
		sprintf(fnBuffer,"%s_specStats.bin",filename);
		Save_binArray(fnBuffer, specStats);
		sprintf(fnBuffer,"%s_alignStats.bin",filename);
		Save_binArray(fnBuffer, alignStats);
		free(fnBuffer);
	}

	filename = params.getValue("OUTPUT_MEANS_CLUSTERS");
	if(params.paramPresent("OUTPUT_MEANS_CLUSTERS")) {
	    output.open(filename, ios::binary);
        for(unsigned int i=0; i<meansClst.size(); i++) output << meansClst[i][0] << " " << meansClst[i][1] << endl;
	    output.close();
	}

	filename = params.getValue("OUTPUT_MEANS");
	if(params.paramPresent("OUTPUT_MEANS")) {
	    output.open(filename, ios::binary);
        for(unsigned int i=0; i<meansSpecs.size(); i++) output << meansSpecs[i][0] << " " << meansSpecs[i][1] << endl;
	    output.close();
	}

	filename = params.getValue("OUTPUT_VARIANCE");
	if(params.paramPresent("OUTPUT_VARIANCE")) {
	    output.open(filename, ios::binary);
        for(unsigned int i=0; i<varTerms.size(); i++) output << varTerms[i] << endl;
	    output.close();
	}

    return(0);
}
