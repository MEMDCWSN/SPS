#include "inputParams.h"
#include "denovo.h"
#include "clusters.h"
#include "batch.h"
#include "msn.h"
#include "spectrum_scoring.h"
#include "utils.h"
//#include "mzxml.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>

using namespace std;

int main(int argc, char **argv){
	InputParams params;
	vector<char *> paramStrs(1);
	
	paramStrs[0] = "OUTPUT_SPECS";
	if(argc==1)	params.readParams("ms3den.params"); else params.readParams(argv[1]);
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"ms3den.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: OUTPUT_SPECS\n";
		return -1;
	}
	if(!params.paramPresent("INPUT_CLUSTERS") and !params.paramPresent("INPUT_ALIGNS")) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"ms3den.params"; else cerr<<argv[1];
		cerr << " is incomplete. At least one of the following must be specified: INPUT_ALIGNS, INPUT_CLUSTERS\n";
		return -1;
	}
	
	char *outputSpecsFN = params.getValue(paramStrs[0]);
	char *outputSpecsIdxFN = params.paramPresent("OUTPUT_SPECS_IDX")?params.getValue("OUTPUT_SPECS_IDX"):(char *)0;
	char *aaMassesFN = params.paramPresent("AMINO_ACID_MASSES")?params.getValue("AMINO_ACID_MASSES"):(char *)0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float resolution = params.paramPresent("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
	int maxAAjump = params.paramPresent("MAX_AA_JUMP")?params.getValueInt("MAX_AA_JUMP"):2;
	float mpPenalty = params.paramPresent("MISSING_PEAK_PENALTY")?(float) params.getValueDouble("MISSING_PEAK_PENALTY"):-1.0;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
    bool allSpecsMS2 = params.paramPresent("ALL_SPECS_MS2")?((bool) params.getValueInt("ALL_SPECS_MS2")):false;
	bool enforcePM = params.paramPresent("STRICT_PM_TOLERANCE")?((int) params.getValueInt("STRICT_PM_TOLERANCE")?1:0):1;
    unsigned int idxStart = params.paramPresent("IDX_START")?(unsigned int) params.getValueInt("IDX_START"):0,
                 idxEnd   = params.paramPresent("IDX_END")?(unsigned int) params.getValueInt("IDX_END"):0;
	float ionOffset = specType?AAJumps::massHion:0;
//	AAJumps jumps(maxAAjump,resolution,peakTol);	cout<<maxAAjump<<" amino acid jumps enforced in denovo interpretations. Parent mass tolerance"<<(enforcePM?" ":" not ")<<"strictly enforced.\n";
//	AAJumps jumps(1);	jumps.alljumps(1500,resolution,peakTol);  cout<<"All valid jumps of mass <=1500 Da accepted in denovo interpretations. Parent mass tolerance"<<(enforcePM?" ":" not ")<<"strictly enforced.\n";
	AAJumps jumps(-1);    if(aaMassesFN) jumps.loadJumps(aaMassesFN);    jumps.alljumps(750,resolution,peakTol);  
	cout<<"All valid jumps of mass <=750 Da ("<<jumps.size()<<" jumps, "; if(aaMassesFN) cout<<"AAmasses from "<<aaMassesFN; else cout<<"standard AAmasses"; cout<<") accepted in denovo interpretations. Parent mass tolerance"<<(enforcePM?" ":" not ")<<"strictly enforced.\n";
//for(unsigned int i=0; i<jumps.size(); i++) cerr<<jumps[i]<<endl;
	Spectrum consensus;

	cout << setprecision(6);   cerr << setprecision(6);

    SpecSet specSet, curMS3set;  int loadOk=0;
    if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
    else if(params.paramPresent("INPUT_SPECS_PKLBIN")) loadOk=specSet.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS_PKLBIN"));
    if (loadOk<=0 or specSet.size()==0) { cerr<<"Error loading spectra!\n"; return -1; }
    cout << "Loading specs complete. Num specs: " << specSet.size() <<"\n";

	MS2ScoringModel ms2model;
	MS3ScoringModel ms3modelB, ms3modelY;
	if(params.paramPresent("MS2_MODEL")) {
		cout << "Using MS/MS model "<<params.getValue("MS2_MODEL")<<endl;
		if(!ms2model.LoadModel(params.getValue("MS2_MODEL"))) return -1;
	}
	if(params.paramPresent("MS3B_MODEL")) {
//		cout << "Using MS/MS/MS/b-ion model "<<params.getValue("MS3B_MODEL")<<endl;
//		if(!ms3modelB.LoadModel(params.getValue("MS3B_MODEL"))) return -1;
	}
	if(params.paramPresent("MS3Y_MODEL")) {
		cout << "Using MS/MS/MS/y-ion model "<<params.getValue("MS3Y_MODEL")<<endl;
		if(!ms3modelY.LoadModel(params.getValue("MS3Y_MODEL"))) return -1;
		cout << "Using MS/MS/MS/b-ion model "<<params.getValue("MS3Y_MODEL")<<endl;
		if(!ms3modelB.LoadModel(params.getValue("MS3Y_MODEL"))) return -1;
	}
    
	Clusters ms3sets;
	if(params.paramPresent("INPUT_CLUSTERS"))
		if(ms3sets.Load(params.getValue("INPUT_CLUSTERS"))<=0) { cerr<<"ERROR opening "<<params.getValue("INPUT_CLUSTERS")<<"\n"; return -1; }
	
	vector<Results_ASP> aligns;    vector<list<int> > ms2sets;   unsigned int pivot;
	if(params.paramPresent("INPUT_ALIGNS")) {
		Load_resultsASPbin(params.getValue("INPUT_ALIGNS"), aligns);
//		get_ms3_sets_max(specSet, aligns, ms2sets, 5);    Save_binListArray<int,list<int>,list<int>::iterator>("ms2sets.bla", ms2sets);
		get_ms3_sets_sparse(specSet, aligns, ms2sets, pmTol, 5);    Save_binListArray<int,list<int>,list<int>::iterator>("ms2sets.bla", ms2sets);
		ms3sets.resize(ms2sets.size());
		for(unsigned int setIdx=0; setIdx<ms2sets.size(); setIdx++) {
			ms3sets.specIdx[setIdx].resize(ms2sets[setIdx].size());   pivot=0;
			for(list<int>::iterator iter=ms2sets[setIdx].begin(); iter!=ms2sets[setIdx].end(); iter++)
				ms3sets.specIdx[setIdx][pivot++] = *iter;
			ms2sets[setIdx].clear();
		}
	}

    cout << "Loading sets complete. Num sets: " << ms3sets.size() <<"\n";
	cout.flush();
    if(params.paramPresent("IDX_START")) idxStart = min(idxStart,ms3sets.size()-1);
    if(params.paramPresent("IDX_END")) idxEnd = min(idxEnd,ms3sets.size()-1); else idxEnd=ms3sets.size()-1;

	// Binomial scores
	vector<vector<float> > binomialProbsSignal, binomialProbsNoise, binomialScores;
	Utils::binomial(5,.8,binomialProbsSignal);  // Warning: constant binomial probability of signal = 0.8
	Utils::binomial(5,.04,binomialProbsNoise);  // Warning: constant binomial probability of noise  = 0.04
	Utils::logscores(binomialProbsSignal, binomialProbsNoise, binomialScores);
	
	vector<vector<unsigned int> > scoreDists(idxEnd-idxStart+1);   // Score distributions
	list<Spectrum> curSeqs, allTopSeqs;   vector<vector<unsigned int> > seqsIndex(ms3sets.size());
	for(unsigned int i=0; i<seqsIndex.size(); i++) seqsIndex[i].resize(0);
	for(unsigned int i=0; i<specSet.size(); i++) specSet[i].addZPMpeaks(peakTol,ionOffset,false);
	for(unsigned int setIdx=idxStart; setIdx<=idxEnd; setIdx++) {
		curMS3set.resize(ms3sets.specIdx[setIdx].size());
cerr<<"Set "<<setIdx<<" ("<<ms3sets.specIdx[setIdx].size()<<" spectra) : ";
		if(ms3sets.specIdx[setIdx].size()==0) continue;
		for(unsigned int specIdx=0; specIdx<ms3sets.specIdx[setIdx].size(); specIdx++) {
			if(ms3sets.specIdx[setIdx][specIdx]>=(int)specSet.size()) { cerr<<"ERROR: ms3 set contains index larger than number of spectra: "<<ms3sets.specIdx[setIdx][specIdx]<<" > "<<specSet.size()<<". Exiting...\n"; exit(-1); }
cerr<<"--- specIdx = "<<ms3sets.specIdx[setIdx][specIdx]<<", "<<specSet[ms3sets.specIdx[setIdx][specIdx]].size()<<" peaks\n";
			curMS3set[specIdx] = specSet[ms3sets.specIdx[setIdx][specIdx]];
		}
		
		denovo_ms3(curMS3set, allSpecsMS2, ionOffset, peakTol, pmTol, resolution, mpPenalty, ms2model, ms3modelB, ms3modelY, jumps, ms3sets.consensus[setIdx], &curSeqs, enforcePM);
//		denovo_ms3(curMS3set, allSpecsMS2, ionOffset, peakTol, pmTol, resolution, mpPenalty, ms2model, ms3modelB, ms3modelY, jumps, ms3sets.consensus[setIdx], &curSeqs, enforcePM, &binomialScores);
cerr<<"Size of consensus: "<<ms3sets.consensus[setIdx].size()<<"\n";
		if(outputSpecsIdxFN) {
			if(curSeqs.size()==0) { seqsIndex[setIdx].resize(1); seqsIndex[setIdx][0]=allTopSeqs.size()+1; allTopSeqs.push_back(ms3sets.consensus[setIdx]); }
			else {
				seqsIndex[setIdx].resize(curSeqs.size());   for(unsigned int i=0; i<seqsIndex[setIdx].size(); i++) seqsIndex[setIdx][i]=allTopSeqs.size()+i+1;
				allTopSeqs.splice(allTopSeqs.end(),curSeqs);
			}
		}
		
		// Save scored spectra for debugging purposes
		for(unsigned int specIdx=0; specIdx<ms3sets.specIdx[setIdx].size(); specIdx++)
			specSet[ms3sets.specIdx[setIdx][specIdx]] = curMS3set[specIdx];

//curMS3set[0].output(cout);
		time_t time_before, time_after;
		if(params.paramPresent("OUTPUT_SCORE_DISTS")) {
cerr<<"Computing distribution of peptide scores (parent mass "<<curMS3set[0].parentMass<<")... "; cerr.flush();
			time(&time_before);
			EstimateScoresDistribution(curMS3set[0], peakTol, pmTol, resolution, scoreDists[setIdx-idxStart], false);
			time(&time_after);
//		EstimateScoresDistribution(specSet[ms3sets.specIdx[setIdx][0]], .5, .05, .1, scoreDists[setIdx-idxStart], false);
cerr<<"done in "<<difftime(time_after, time_before)<<" seconds\n";
		}
	}

	if(outputSpecsIdxFN) {
		Save_binListArray<unsigned int,vector<unsigned int>,vector<unsigned int>::iterator>(outputSpecsIdxFN,seqsIndex);
		SpecSet allSeqs(allTopSeqs.size());   int specsIdx=0;
		list<Spectrum>::iterator seqsIter = allTopSeqs.begin();
		while(seqsIter!=allTopSeqs.end()) { allSeqs[specsIdx++]=*seqsIter; seqsIter=allTopSeqs.erase(seqsIter); }
		allSeqs.SaveSpecSet_pklbin(outputSpecsFN);
		ms3sets.consensus.SaveSpecSet_pklbin("debug_top_denovo_seq_only.pklbin");
	} else ms3sets.consensus.SaveSpecSet_pklbin(outputSpecsFN);  // Outputs only one denovo reconstruction

	if(params.paramPresent("OUTPUT_SCORED_SPECS")) specSet.SaveSpecSet_pklbin(params.getValue("OUTPUT_SCORED_SPECS"));
	if(params.paramPresent("OUTPUT_SCORE_DISTS")) Save_binArray(params.getValue("OUTPUT_SCORE_DISTS"), scoreDists);
	
    return(0);
}
