#include "inputParams.h"
#include "alignment_scoring.h"
#include "batch.h"
#include "spectral_pairs.h"
#include "SpectralPairs.h"
#include "SpectrumPairSet.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;
using namespace specnets;

int main(int argc, char **argv){

    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("starden.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "starden.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs;   paramStrs.resize(4);
	//	paramStrs[0] = "INPUT_SPECS";
	paramStrs[0] = "INPUT_ALIGNS";
	paramStrs[1] = "OUTPUT_SPECS";
	paramStrs[2] = "OUTPUT_MODPOS";
	paramStrs[3] = "PENALTY_SAME_VERTEX";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"starden.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_ALIGNS, OUTPUT_SPECS, OUTPUT_MODPOS, PENALTY_SAME_VERTEX\n";
		return -1;
	}

	const char *alignsFN = params.getValue("INPUT_ALIGNS");
	const char *resultsFN = params.getValue("OUTPUT_SPECS");
	const char *modPosFN = params.getValue("OUTPUT_MODPOS");

/*
 * WARNING: OUTPUT_MODPOS assumes the modification mass is positive and contains
 *   the location of the mass offset on the spectrum with LOWER parent mass.
 */

	float penalty_sameVert = (float) params.getValueDouble("PENALTY_SAME_VERTEX");

//	int aaDiff = params.getValueInt(paramStrs[4]);
//	double minShift = params.getValueDouble(paramStrs[5]);
//	double maxShift = params.getValueDouble(paramStrs[6]);

	float penalty_ptm=params.paramPresent("PENALTY_PTM")?(float)params.getValueDouble("PENALTY_PTM"):0.0;

	// If specified, indicates the minimum number peaks on each side of a modification to _not_ make it a terminal modification
	//   This is encoded by computing the average prm score per spectrum (at least 57 Da apart) and multiplying it by PENALTY_PTM_PEAKS to obtain the penalty (per spectrum)
	float penalty_ptm_peaks=params.paramPresent("PENALTY_PTM_PEAKS")?(float)params.getValueInt("PENALTY_PTM_PEAKS"):-1.0;

	int startBaseIdx = params.paramPresent("IDX_START")?params.getValueInt("IDX_START"):0;
	int endBaseIdx = params.paramPresent("IDX_END")?params.getValueInt("IDX_END"):-1;
	int startStarIdx = params.paramPresent("STARS_START")?params.getValueInt("STARS_START"):0;
	int endStarIdx = params.paramPresent("STARS_END")?params.getValueInt("STARS_END"):0;
	int maxAAjump = params.paramPresent("MAX_AA_JUMP")?params.getValueInt("MAX_AA_JUMP"):0;
//    vector<int> baseSpecIdx(endBaseIdx-startBaseIdx+1);  for(int i=0; i<endBaseIdx-startBaseIdx+1; i++) baseSpecIdx[i]=startBaseIdx+i;
//	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0;
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float resolution = params.paramPresent("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;
	bool alignPA = params.paramPresent("PARTIAL_OVERLAPS")?(bool) params.getValueInt("PARTIAL_OVERLAPS"):0;
	if(params.paramPresent("AMINO_ACID_MASSES")) { AAJumps jumps(-1); jumps.loadJumps(params.getValue("AMINO_ACID_MASSES"),true); }

	//    SpecSet specSet;
    SpecSet specSet;   short loadOk=0;   unsigned int specIdx;
    if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
    else if(params.paramPresent("INPUT_SPECS_PKLBIN")) loadOk=specSet.loadPklBin(params.getValue("INPUT_SPECS_PKLBIN"));
    if (loadOk<=0 or specSet.size()==0) return -1;
    cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";

	vector<float> penalties(specSet.size());
	vector<int> idx1dp(2048), idx2dp(2048); // Indices of matched peaks for partial overlaps (default sizes should avoid most resize() operations)
	float score;
    for(specIdx=0; specIdx<specSet.size(); specIdx++) {
		specSet[specIdx].parentMass+=2*ionOffset;   // Quick fix to support alignment of MS/MS spectra
		if(penalty_ptm_peaks>0) {
			score=ScoreOverlap6(specSet[specIdx],specSet[specIdx],0,peakTol,idx1dp,idx2dp,AAJumps::minAAmass);
			if(idx1dp.size()==0) penalties[specIdx] = 0;
			else penalties[specIdx] = -penalty_ptm_peaks*score/((float)idx1dp.size());
		}
		specSet[specIdx].addZPMpeaks(peakTol,ionOffset,true);
    }

	SpectrumPairSet aligns;   unsigned int idxPair;
	if(alignPA) {
		SpectrumPairSet alignsPA;
		if (!alignsPA.loadFromBinaryFile(alignsFN)) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
		aligns.resize(alignsPA.size());
		for(idxPair=0;idxPair<aligns.size();idxPair++) aligns[idxPair]=alignsPA[idxPair];  // Loses shift2 (not needed below)
	} else if (!aligns.loadFromBinaryFile(alignsFN)) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
    cout << "Loading aligns complete. Num pairs: " << aligns.size() << "\n";
    if(aligns.size()==0) { cout<<"No work to do... exiting.\n"; return 0; }
    if(endBaseIdx<=0) endBaseIdx = aligns.size()-1;

#ifdef DEBUG
	ofstream debug("dekel_align_debug.txt");
#endif

	SpecSet resultSpecsSet;
	int idxResultSpecs=0, spec1, spec2, idxMatch;   float shift, curPenalty;
    if(params.paramPresent("INPUT_SPECS_PAIRS")) {
    	if(not resultSpecsSet.loadPklBin(params.getValue("INPUT_SPECS_PAIRS"))
    			or resultSpecsSet.size()!=2*aligns.size())
    		{ cerr<<"ERROR loading "<<params.getValue("INPUT_SPECS_PAIRS")<<"!\n"; return -1; }
    	cout << "Loading matched peaks complete. Num pairs: " << resultSpecsSet.size()/2 << "\n";
    } else {
    	resultSpecsSet.specs.resize(2*aligns.size());  // 2 results per pair of spectra
    	vector<float> modPositions;   modPositions.resize(aligns.size());        // Keeps track of the mass where the mod was placed
    	cout <<"Computing pairs from "<<startBaseIdx<<" to "<<endBaseIdx<<"..."; cout.flush();
    	for (idxPair=startBaseIdx; idxPair<=endBaseIdx; idxPair++, idxResultSpecs+=2) {
    		if(idxPair%50000==0){ cout<<"Processing pairs "<<idxPair<<" -> "<<endBaseIdx<<", "<<100.0*((float)idxPair-(float)startBaseIdx)/((float)endBaseIdx-(float)startBaseIdx)<<"% completed\n"; cout.flush(); }
    		spec1 = aligns[idxPair].spec1;    spec2 = aligns[idxPair].spec2;   shift = aligns[idxPair].shift1;
    		resultSpecsSet[idxResultSpecs].copyNP(specSet[spec1]);
    		resultSpecsSet[idxResultSpecs+1].copyNP(specSet[spec2]);
    		if(penalty_ptm_peaks>0) curPenalty=penalties[spec1]+penalties[spec2];
    		else curPenalty=penalty_ptm;

//    		if(abs(shift)<=2*pmTol and abs(specSet[spec2].parentMass+shift-specSet[spec1].parentMass)<=2*pmTol) continue; // Same-peptide pairs can't be used to compute spectral stars
    		if(not alignPA or abs(shift)<=2*pmTol or abs(specSet[spec2].parentMass+shift-specSet[spec1].parentMass)<=2*pmTol) {  // Common start or common end
#ifdef DEBUG
    			debug << "\nCurrent pair is ("<<spec1+1<<", "<<spec2+1<<"): \n";
    			modPositions[idxPair]=SpectrumAlignment(&specSet[spec1],&specSet[spec2],peakTol,&resultSpecsSet[idxResultSpecs],&resultSpecsSet[idxResultSpecs+1],maxAAjump,penalty_sameVert,penalty_ptm,debug);
#else
//cerr<<"Spectrum 1:\n"; specSet[spec1].output(cerr);
//cerr<<"Spectrum 2:\n"; specSet[spec2].output(cerr);
//cerr<<"peakTol = "<<peakTol<<", maxAAjump = "<<maxAAjump<<", penalty_sameVert = "<<penalty_sameVert<<", curPenalty = "<<curPenalty<<"\n";
    			modPositions[idxPair]=SpectrumAlignment(&specSet[spec1],&specSet[spec2],peakTol,&resultSpecsSet[idxResultSpecs],&resultSpecsSet[idxResultSpecs+1],maxAAjump,penalty_sameVert,curPenalty);
//for(unsigned int matchIdx=0; matchIdx<resultSpecsSet[idxResultSpecs].size(); matchIdx++)
//	cerr<<" * ("<<resultSpecsSet[idxResultSpecs][matchIdx][0]<<","<<resultSpecsSet[idxResultSpecs][matchIdx][1]<<") - ("
//		<<resultSpecsSet[idxResultSpecs+1][matchIdx][0]<<","<<resultSpecsSet[idxResultSpecs+1][matchIdx][1]<<")\n";
#endif
    		} else {  // Partial overlap
    			ScoreOverlap6(specSet[spec1],specSet[spec2],shift,peakTol,idx1dp,idx2dp,AAJumps::minAAmass);
    			resultSpecsSet[idxResultSpecs].resize(idx1dp.size());
    			for(idxMatch=0;idxMatch<idx1dp.size();idxMatch++) resultSpecsSet[idxResultSpecs][idxMatch]=specSet[spec1][idx1dp[idxMatch]];
    			resultSpecsSet[idxResultSpecs+1].resize(idx2dp.size());
    			for(idxMatch=0;idxMatch<idx2dp.size();idxMatch++) resultSpecsSet[idxResultSpecs+1][idxMatch]=specSet[spec2][idx2dp[idxMatch]];
    			modPositions[idxPair]=0;
    		}
    	}
    	cout << "done\n"; cout.flush();
#ifdef DEBUG
    	debug.close();
#endif
    	// Quick fix to support alignment of MS/MS spectra
    	for(specIdx=0; specIdx<resultSpecsSet.size(); specIdx++)
    		resultSpecsSet[specIdx].parentMass-=2*ionOffset;

    	resultSpecsSet.SaveSpecSet_pklbin(resultsFN);

    	Save_binArray(modPosFN, modPositions);
    	ofstream outs("old_modPos.txt");
    	for(idxPair=0; idxPair<(int)modPositions.size(); idxPair++) outs<<modPositions[idxPair]<<endl;
    	outs.close();
    }

	if(params.paramPresent("OUTPUT_STARS") or params.paramPresent("OUTPUT_STARS_ALL")) {
		unsigned int numStars=0, starIdx=0;
		SpecSet stars, curSpecs;
		vector<unsigned short> numPairs(specSet.size());  for(specIdx=0; specIdx<numPairs.size(); specIdx++) numPairs[specIdx]=0;
		vector<unsigned int> index;
		if(endStarIdx<=0 or endStarIdx>=specSet.size()) endStarIdx=specSet.size()-1;
		cout<<"Computing stars with indices from "<<startStarIdx<<" to "<<endStarIdx<<"..."; cout.flush();
		for (idxPair=startBaseIdx; idxPair<=endBaseIdx; idxPair++) {
			spec1 = aligns[idxPair].spec1;    spec2 = aligns[idxPair].spec2;   shift = aligns[idxPair].shift1;
			if(abs(shift)<=2*pmTol and abs(specSet[spec2].parentMass+shift-specSet[spec1].parentMass)<=2*pmTol) continue; // Same-peptide pairs can't be used to compute spectral stars
			if((spec1<startStarIdx or spec1>endStarIdx) and (spec2<startStarIdx or spec2>endStarIdx)) continue; // Star-spectrum will be computed by another instance
			if(spec1>=startStarIdx and spec1<=endStarIdx){ if(numPairs[spec1]==0) numStars++; numPairs[spec1]++; }
			if(spec2>=startStarIdx and spec2<=endStarIdx){ if(numPairs[spec2]==0) numStars++; numPairs[spec2]++; }
		}
		stars.resize(numStars);   index.resize(numStars);
		cout<<numStars<<" stars in range..."; cout.flush();

		for(specIdx=startStarIdx; specIdx<=endStarIdx; specIdx++) {
			if(numPairs[specIdx]==0) continue;
			index[starIdx]=specIdx;
			if(starIdx % 100 ==0) cout<<"Computing star spectrum "<<starIdx+1<<"/"<<numStars<<"\n";

			// Find all spectral pairs for this spectrum
			curSpecs.resize(numPairs[specIdx]);   numPairs[specIdx]=0;   idxResultSpecs=0;
			for(idxPair=startBaseIdx; idxPair<=endBaseIdx; idxPair++) {
				spec1 = aligns[idxPair].spec1;    spec2 = aligns[idxPair].spec2;   shift = aligns[idxPair].shift1;
//				if(abs(shift)<=2*pmTol and abs(specSet[spec2].parentMass+shift-specSet[spec1].parentMass)<=2*pmTol) continue; // Same-peptide pairs can't be used to compute spectral stars
				if(spec1==(int)specIdx)
					curSpecs[numPairs[specIdx]++] = resultSpecsSet[2*idxPair];
				else if (spec2==(int)specIdx)
					curSpecs[numPairs[specIdx]++] = resultSpecsSet[2*idxPair+1];
				if(numPairs[specIdx]==curSpecs.size()) break;
			}
			if(numPairs[specIdx]==1) stars[starIdx]=curSpecs[0];
			else ComputeSpectralStars(curSpecs, stars[starIdx], peakTol, resolution);
			starIdx++;
		}

		cout << "done\n"; cout.flush();

		if(params.paramPresent("OUTPUT_STARS")) stars.SaveSpecSet_pklbin(params.getValue("OUTPUT_STARS"));
		if(params.paramPresent("OUTPUT_STARS_INDEX")) Save_binArray(params.getValue("OUTPUT_STARS_INDEX"), index);

		if(params.paramPresent("OUTPUT_STARS_ALL")) {
			// Fill-in spectra that were not converted to spectral stars
			for(starIdx=0; starIdx<index.size(); starIdx++) specSet[index[starIdx]]=stars[starIdx];
			specSet.SaveSpecSet_pklbin(params.getValue("OUTPUT_STARS_ALL"));
		}
//		if(params.paramPresent("OUTPUT_SPECS_STARS_COUNTS"))
//			Save_binArray(params.getValue("OUTPUT_SPECS_STARS_COUNTS"), numPairs);
	}

return 0;

	resultSpecsSet.loadPklBin("pairs.pklbin");
	specSet.resize(0);   // Free some memory before loading MS/MS spectra
	if(params.paramPresent("MATCH_MSMS_INTENSITY") and (params.paramPresent("MATCH_MSMS_RATIOS") or params.paramPresent("MIN_MATCH_MSMS"))) {
		float minRatio=params.paramPresent("MIN_MATCH_MSMS")?params.getValueDouble("MIN_MATCH_MSMS"):-1.0;
		vector<vector<float> > ms2Ratios(aligns.size());
		Spectrum masses1, masses2;   float norm1, norm2;
		vector<int> idx1, idx2;
		unsigned int idxPeak, pivot;

		if(not specSet.loadPklBin(params.getValue("MATCH_MSMS_INTENSITY")))
			{ cerr<<"ERROR loading "<<params.getValue("MATCH_MSMS_INTENSITY")<<" (expected .pklbin format)\n"; exit(-1); }
		vector<float> specIntensity(specSet.size());
		for(pivot=0;pivot<specSet.size();pivot++){
			specIntensity[pivot]=0;
			for(idxPeak=0;idxPeak<specSet[pivot].size();idxPeak++) specIntensity[pivot]+=specSet[pivot][idxPeak][1];
		}

		idxResultSpecs=0;
		for(idxPair=0;idxPair<aligns.size();idxPair++) {
			ms2Ratios[idxPair].resize(4);   // % ms2 intensity in spec1/2 (cols 0/1), dot-product beween the not-normalized/normalized vectors (cols 2/3)
//if(idxPair>1000) continue;
if(idxPair%1000==0) cerr << "[ms2ratios] Current pair is "<<idxPair<<"/"<<aligns.size()<<"\n";
			if(resultSpecsSet[idxResultSpecs].size()>0 and resultSpecsSet[idxResultSpecs].size()==resultSpecsSet[idxResultSpecs+1].size()) {
				for(pivot=0;pivot<ms2Ratios[idxPair].size();pivot++) ms2Ratios[idxPair][pivot]=0;

				// Get sets of masses to match
				masses1 = resultSpecsSet[idxResultSpecs];   masses2 = resultSpecsSet[idxResultSpecs+1];
				for(idxPeak=0;idxPeak<masses1.size();idxPeak++) {
					masses1[idxPeak][0]+=AAJumps::massHion; masses1[idxPeak][1]=0;
					masses2[idxPeak][0]+=AAJumps::massHion; masses2[idxPeak][1]=0;
				}

				// Find corresponding intensities in MS/MS spectra
				spec1 = aligns[idxPair].spec1;   spec2 = aligns[idxPair].spec2;
				FindMatchPeaksAll2(masses1,specSet[spec1],0,peakTol,idx1,idx2);
				for(pivot=0;pivot<idx1.size();pivot++) masses1[idx1[pivot]][1]+=specSet[spec1][idx2[pivot]][1]/specIntensity[spec1];
				norm1=0; for(pivot=0;pivot<masses1.size();pivot++)	{
					norm1+=masses1[pivot][1]*masses1[pivot][1];
					ms2Ratios[idxPair][0]+=masses1[pivot][1];
				}
				norm1=sqrt(norm1);
				FindMatchPeaksAll2(masses2,specSet[spec2],0,peakTol,idx1,idx2);
				for(pivot=0;pivot<idx1.size();pivot++) masses2[idx1[pivot]][1]+=specSet[spec2][idx2[pivot]][1]/specIntensity[spec2];
				norm2=0; for(pivot=0;pivot<masses2.size();pivot++)	{
					norm2+=masses2[pivot][1]*masses2[pivot][1];
					ms2Ratios[idxPair][1]+=masses2[pivot][1];
				}
				norm2=sqrt(norm2);

				// Compute dot-products between the vectors
//				float v;
//				for(idxPeak=0;idxPeak<masses1.size();idxPeak++) { v=masses1[pivot][1]-masses2[pivot][1]; ms2Ratios[idxPair][2]+=v*v; }
				for(idxPeak=0;idxPeak<masses1.size();idxPeak++) ms2Ratios[idxPair][2]+=masses1[idxPeak][1]*masses2[idxPeak][1];
				ms2Ratios[idxPair][3]= ms2Ratios[idxPair][2]/(norm1*norm2);
			} else for(pivot=0;pivot<ms2Ratios[idxPair].size();pivot++) ms2Ratios[idxPair][pivot]=-1.0;
			idxResultSpecs+=2;
		}
		if(params.paramPresent("MATCH_MSMS_RATIOS")) Save_binArray(params.getValue("MATCH_MSMS_RATIOS"),ms2Ratios);
	}

	return(0);
}
