#include "alignment_scoring.h"
#include "spectrum_scoring.h"
#include "batch.h"
#include "filters.h"
#include "inputParams.h"
#include "denovo.h"
#include "db_fasta.h"
#include "msn.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

void Cyclize(Spectrum &in, Spectrum &out, float ionOffset, float mergeTol=0.5, float resolution=0.1);

int main(int argc, char **argv){

//	unit_MatchSpecsToPeps();
//	return 0;

	InputParams params;
	vector<char *> paramStrs;   paramStrs.resize(2);

	bool ok;
	if(argc==1)	ok=params.readParams("cycseq.params"); else ok=params.readParams(argv[1]); if(!ok) return -1;
/*	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"cycseq.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_SPECS, OUTPUT_ALIGNS\n";
		return -1;
	}
*/
	char *specSetFN = params.getValue("INPUT_SPECS");
	char *alignsFN = params.getValue("OUTPUT_ALIGNS");

	int startIdx = params.paramPresent("IDX_START")?params.getValueInt("IDX_START"):0;
	int endIdx = params.paramPresent("IDX_END")?params.getValueInt("IDX_END"):-1;
	// Alignment parameters
	float minRatio = params.paramPresent("MIN_RATIO")?(float) params.getValueDouble("MIN_RATIO"):0;
	float minOverlap = params.paramPresent("MIN_OVERLAP_AREA")?(float) params.getValueDouble("MIN_OVERLAP_AREA"):0;
	short minNumMatchedPeaks = params.paramPresent("MIN_NUM_MATCHED_PEAKS")?(short) params.getValueInt("MIN_NUM_MATCHED_PEAKS"):0;
	float minAbsShift = params.paramPresent("MIN_SHIFT")?(float) params.getValueDouble("MIN_SHIFT"):0;
	// Mass accuracy parameters
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float pmTol = params.paramPresent("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	float resolution = params.paramPresent("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
	float msn_peakTol = params.paramPresent("MSN_TOLERANCE_PEAK")?(float) params.getValueDouble("MSN_TOLERANCE_PEAK"):0.5;
	float msn_pmTol = params.paramPresent("MSN_TOLERANCE_PM")?(float) params.getValueDouble("MSN_TOLERANCE_PM"):1;
	float msn_resolution = params.paramPresent("MSN_RESOLUTION")?(float) params.getValueDouble("MSN_RESOLUTION"):0.1;
	// Sequencing parameters
	float minAAmass = params.paramPresent("MIN_AA_MASS")?(float) params.getValueDouble("MIN_AA_MASS"):50;
	float maxAAmass = params.paramPresent("MAX_AA_MASS")?(float) params.getValueDouble("MAX_AA_MASS"):-1;
    int maxAAjump = params.paramPresent("MAX_AA_JUMP")?(int) params.getValueInt("MAX_AA_JUMP"):1;
	float minPercTopDenovo = params.paramPresent("MIN_PERC_TOP_DENOVO")?(float) params.getValueDouble("MIN_PERC_TOP_DENOVO"):0.75;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
    int mergeType = params.paramPresent("PEAK_MERGE_TYPE")?(int) params.getValueInt("PEAK_MERGE_TYPE"):2;
	float ionOffset = specType?AAJumps::massHion:0;
	cout << "Minimum match ratio is " << minRatio << "\n";

    SpecSet specSet;
	list<Results_PA> results;
	list<TwoValues<float> > ratios;
	list<TwoValues<int> > numMatchedPeaks;
	vector<TwoValues<float> > meansSpecs;
	vector<float> varTerms;
	list<vector<float> > alignStats;    alignStats.clear();
	vector<vector<float> > specStats;   specStats.resize(0);
	SpecSet mergedSpecs;       // Consensus spectra input to denovo sequencing
	float shift;
	unsigned int specIdx, peakIdx;

	// Compound intensity of neutral losses
	vector<float> newLosses(2), noLosses(0);
	newLosses[0]=17.0;  // -NH3
	newLosses[1]=18.0;  // -H20
//	newLosses[2]=28.0;  // a-ion
//  newLosses[3]=46.0;  // a-H20
	AAJumps jumpsLosses(0);    jumpsLosses.addJumps(newLosses);

	SpecSet correctPeps;   correctPeps.loadPklBin("correct_peps.pklbin");
	vector<int> idxPep, idxSeq;

	if(params.paramPresent("INPUT_SPECS")) {
	    if(specSet.LoadSpecSet_pkl(specSetFN)<=0) { cerr<<"ERROR reading "<<specSetFN<<"!\n"; return -1; }
		for(specIdx=0; specIdx<specSet.size(); specIdx++) specSet[specIdx].addZPMpeaks(peakTol,ionOffset,false,false);

	    if(endIdx<0) endIdx = specSet.size()-1;
	    cout << "Loading specs complete. Num specs: " << specSet.size() <<", computing pairwise aligns for indices "<<startIdx<<" to "<<endIdx<<"\n";
		AAJumps jumps(0);

		vector<vector<float> > shifts;   Load_binArray("shifts.bin",shifts);

/*
		SpecSet specSetUI(2);   specSetUI[0]=specSet[0];   // Unit intensities for counting number of matched peaks
		for(peakIdx=0; peakIdx<specSetUI[0].size(); peakIdx++) specSetUI[0][peakIdx][1]=1;
		specSetUI[1]=specSetUI[0];

		vector<TwoValues<float> > scoredShifts;
		CyclicAlign(specSetUI[0], specSetUI[0], minAbsShift, minAAmass, ionOffset, peakTol, resolution, 0, scoredShifts);
		Save_binArray("autoconvolution_num_peaks.bin", scoredShifts);
		CyclicAlign(specSet[0], specSet[0], minAbsShift, minAAmass, ionOffset, peakTol, resolution, 0, scoredShifts);
		Save_binArray("autoconvolution_expint.bin", scoredShifts);


		getPairAlignsPA2(specSetUI, 0, 0, peakTol, pmTol, minRatio, minOverlap, minNumMatchedPeaks, jumps, minAbsShift, results, ratios, numMatchedPeaks, meansSpecs, varTerms, alignStats, specStats);
		results.front().spec2=0;
		shift = max(results.front().shift1,results.front().shift2);  // Choose the positive shift
		if(params.paramPresent("FORCE_MS3_SHIFT")) shift = (float) params.getValueDouble("FORCE_MS3_SHIFT");


		SpecSet specSet2;   specSet2.resize(specSet.size());  for(specIdx=0;specIdx<specSet.size();specIdx++) specSet2[specIdx]=specSet[specIdx];

		// Convert spectrum to unit intensities - i.e. work with matched #peaks rather than intensities
//		for(specIdx=0; specIdx<2; specIdx++)
//			for(peakIdx=0; peakIdx<specSet2[specIdx].size(); peakIdx++)
//				specSet2[specIdx][peakIdx][1]=1;

		// Process input spectra to add neutral-loss intensities to main peaks + make spectra symmetric (peak masses)
		vector<vector<int> > neighLosses;   float prevIntensity;
		for(specIdx=0; specIdx<min(2,(int)specSet.size()); specIdx++) {
			getNeighborsL(specSet2[specIdx], jumpsLosses, peakTol, neighLosses);
			for(peakIdx=specSet2[specIdx].size()-1; peakIdx>0; peakIdx--) { prevIntensity=specSet2[specIdx][peakIdx][1];
				for(unsigned int lossIdx=0; lossIdx<neighLosses[peakIdx].size(); lossIdx++)
//					if(specSet2[specIdx][neighLosses[peakIdx][lossIdx]][1] < prevIntensity)
						specSet2[specIdx][peakIdx][1]+=specSet2[specIdx][neighLosses[peakIdx][lossIdx]][1];
			}
			specSet2[specIdx].makeSymmetric(2*ionOffset,peakTol);    // add missing symmetric masses
		}
		specSet2.SaveSpecSet_pkl("specSet2.pkl");

		// Make spectra symmetric (peak intensities)
		SpecSet specSet2r(2);
		specSet2[0].reverse(2*ionOffset, &specSet2r[0]);
		if(specSet.size()>1) specSet2[1].reverse(2*ionOffset, &specSet2r[1]);
		specSet2r.SaveSpecSet_pkl("specSet2r.pkl");

		// Add spectrum with reversed version
		specSet2[0].mergeCommon(specSet2r[0],0,0,peakTol,0);
		if(specSet.size()>1) specSet2[1].mergeCommon(specSet2r[1],0,0,peakTol,0);
		specSet2.SaveSpecSet_pkl("specSet2s.pkl");

		SpecSet specSetF(2);  // Retains only masses present in the symmetric-mass spectra
		specSet2[0].mergeCommon(specSet[0],&specSetF[0],0,peakTol,3);   specSetF[0].addZPMpeaks(peakTol, ionOffset, true);
		if(specSet.size()>1) specSet2[1].mergeCommon(specSet[1],&specSetF[1],0,peakTol,3);   specSetF[1].addZPMpeaks(peakTol, ionOffset, true);
		specSet2.SaveSpecSet_pkl("specSetF.pkl");

		// Use best shift to generate merged version of the spectra
		mergedSpecs.resize(3);
	//	specSet[0].mergeCommon(specSet[0],mergedSpecs[0],shift,peakTol,0);
	//	specSetF[0].mergeCommon(specSetF[0],&mergedSpecs[0],shift,peakTol,mergeType);
	//	specSet[0].mergeCommon(specSet[0],&mergedSpecs[0],shift,peakTol,mergeType);

		// MS3+MS4 using neutral losses
//		specSetF[0].mergeCommon(specSetF[0],&mergedSpecs[1],shift,peakTol/2,mergeType);
//		mergedSpecs[1].mergeCommon(specSetF[1],&mergedSpecs[0],shift,peakTol/2,mergeType);

		// MS3+MS4 _without_ using neutral losses
//		specSet[0].mergeCommon(specSet[0],&mergedSpecs[1],shift,peakTol/2,mergeType);

		// Get the MS4-filtered version of the MS3 spectrum
//		mergedSpecs[1].mergeCommon(specSet[1],&mergedSpecs[0],shift,peakTol/2,mergeType);
//		mergedSpecs[0].addZPMpeaks(peakTol, ionOffset, true);

		// Get the MS3-filtered version of the MS4 spectrum
//		specSet[1].mergeCommon(mergedSpecs[1],&mergedSpecs[0],-shift,peakTol/2,mergeType);

		// MS3 only & using neutral losses
//		specSetF[0].mergeCommon(specSetF[0],&mergedSpecs[0],shift,peakTol,mergeType);

		mergedSpecs.resize(3);   shift=shifts[0][0];
		specSet[0].mergeCommon(specSet[0],&mergedSpecs[0],shift,peakTol,mergeType);

//		CyclicAlign(mergedSpecs[0], mergedSpecs[0], minAbsShift, minAAmass, ionOffset, peakTol, resolution, 0, scoredShifts);
//		Save_binArray("shifts_merged_expint.bin", scoredShifts);

		// Select the peaks in the overlapped portion of the auto-convolution
		mergedSpecs[0].parentMass -= shift;
		for(peakIdx=0; peakIdx<mergedSpecs[0].size(); peakIdx++) mergedSpecs[0][peakIdx][0]-= shift;
		mergedSpecs[0].addZPMpeaks(peakTol, ionOffset, false, false);
*/

//		MS2ScoringModel model;   model.LoadModel("dancik_model.txt");

		mergedSpecs.resize(shifts.size());   TwoValues<float> window(-50,50), windowDa(-0.5,0.5);
		for(unsigned int shiftIdx=0; shiftIdx<shifts.size(); shiftIdx++) {
//		mergedSpecs.resize(1);
//		for(unsigned int shiftIdx=0; shiftIdx<1; shiftIdx++) {
			shift = shifts[shiftIdx][0];
			specSet[0].mergeCommon(specSet[0],&mergedSpecs[shiftIdx],shift,peakTol,mergeType);
			mergedSpecs[shiftIdx].parentMass -= shift;
			for(peakIdx=0; peakIdx<mergedSpecs[shiftIdx].size(); peakIdx++) {
				mergedSpecs[shiftIdx][peakIdx][0]-= shift;
//				mergedSpecs[shiftIdx][peakIdx][1] = sqrt(mergedSpecs[shiftIdx][peakIdx][1]);
			}

//			mergedSpecs[shiftIdx].selectTopK(2,windowDa);
//			mergedSpecs[shiftIdx].selectTopK(8,window);
			mergedSpecs[shiftIdx].selectTopK(6,window);

/*			ScoreSpectrum(mergedSpecs[shiftIdx],model,false);
			mergedSpecs[shiftIdx].roundMasses();
			window.set(-0.3,0.3);
			mergedSpecs[shiftIdx].selectTopK(1,window);
//			window.set(-50,50);
//			mergedSpecs[shiftIdx].selectTopK(10,window);
			for(peakIdx=0; peakIdx<mergedSpecs[shiftIdx].size(); peakIdx++)
				mergedSpecs[shiftIdx][peakIdx][0]+=1;
*/
			mergedSpecs[shiftIdx].addZPMpeaks(peakTol, ionOffset, false, false);
			mergedSpecs[shiftIdx].maximizeZPMpeaks(peakTol,ionOffset);
		}
		mergedSpecs.SaveSpecSet_pklbin("merged_specs.pklbin");
		mergedSpecs.SaveSpecSet_pkl("merged_specs.pkl");

/*		mergedSpecs.SaveSpecSet_mgf("merged_specs.mgf");
		system("/home/nbandeira/misc/Pepnovo_newest/PepNovo_bin -model_dir /home/nbandeira/misc/Pepnovo_newest/Models -model CID_IT_TRYP -prm -file merged_specs.mgf -use_spectrum_mz -use_spectrum_charge -digest NON_SPECIFIC > pepnovo_specs.prms");
		SpecSet pepnovoSpecs;   pepnovoSpecs.LoadSpecSet_prmsv3("pepnovo_specs.prms");
		if(pepnovoSpecs.size()!=mergedSpecs.size()) { cerr<<"ERROR running pepnovo!\n"; return -1; }
		Spectrum tmp;
		for(specIdx=0; specIdx<mergedSpecs.size(); specIdx++) {
			for(peakIdx=0; peakIdx<pepnovoSpecs[specIdx].size(); peakIdx++) pepnovoSpecs[specIdx][peakIdx][0]+=1.0;
			for(peakIdx=0; peakIdx<mergedSpecs[specIdx].size(); peakIdx++) mergedSpecs[specIdx][peakIdx][1]=0;
			MergeIntoReference(mergedSpecs[specIdx],pepnovoSpecs[specIdx],peakTol,0.1,tmp);
			mergedSpecs[specIdx]=tmp;
		}
*/
//		return 0;


//		mergedSpecs.LoadSpecSet_pklbin("merged_specs_matlab.pklbin");
//		for(specIdx=0; specIdx<mergedSpecs.size(); specIdx++)
//			mergedSpecs[specIdx].addZPMpeaks(peakTol, ionOffset, false, false);
//		cout << "Loaded "<<mergedSpecs.size()<<" spectra from merged_specs_matlab.pklbin\n";
	}
	cout<<"Done with consensus, starting sequencing...\n"; cerr.flush();

	// Use denovoLR to output all suboptimal de novo sequences
	AAJumps jumpsDenovo(-1);
	if(params.paramPresent("AMINO_ACID_MASSES")) jumpsDenovo.loadJumps(params.getValue("AMINO_ACID_MASSES"));
	else jumpsDenovo.getjumps(1);
	if(params.paramPresent("MAX_AA_JUMP")) jumpsDenovo.getjumps(maxAAjump,resolution,peakTol);

	vector<vector<int> > neighsL, neighsR;
	list<list<int> > seqs, curSeqs;
	list<Spectrum> seqsList, cur_seqsList;   list<float> seqScores, cur_seqScores;
	SpecSet denovoSpecs;
  	vector<TwoValues<float> > seqsPairs;

	if(params.paramPresent("INPUT_DENOVO")) {
		if(denovoSpecs.LoadSpecSet_pkl(params.getValue("INPUT_DENOVO"))<0) { cerr<<"ERROR reading "<<params.getValue("INPUT_DENOVO")<<"!\n"; return -1; }
		denovoSpecs.SaveSpecSet_pklbin("specs_denovo.pklbin");  // Just to facilitate running Matlab benchmarking scripts
	} else {
		if(params.paramPresent("INPUT_CONSENSUS")) {
			if(mergedSpecs.LoadSpecSet_pkl(params.getValue("INPUT_CONSENSUS"))<0) { cerr<<"ERROR reading "<<params.getValue("INPUT_CONSENSUS")<<"!\n"; return -1; }
			mergedSpecs[0].addZPMpeaks(peakTol,ionOffset,false,false);
			if(params.paramPresent("MIN_AA_MASS") or params.paramPresent("MAX_AA_MASS"))
				getNeighborsL(mergedSpecs[0], minAAmass, maxAAmass, peakTol, neighsL);
			else getNeighborsL(mergedSpecs[0], jumpsDenovo, peakTol, neighsL);
			denovo_LtoR(mergedSpecs[0], neighsL, seqsList, seqScores, minPercTopDenovo);
		} else {
			for(specIdx=0; specIdx<mergedSpecs.size(); specIdx++) {
				cout<<" ===> sequencing "<<mergedSpecs[specIdx].size()<<" peaks from spectrum "<<specIdx<<"/"<<mergedSpecs.size()<<", parent mass = "<<mergedSpecs[specIdx].parentMass<<", shift = "<<specSet[0].parentMass-mergedSpecs[specIdx].parentMass<<"..."; cout.flush();

				// Sequence the filtered MS3 spectrum
				//		getNeighborsL(mergedSpecs[0], jumpsDenovo, peakTol, neighsL);
				//		getNeighborsL(mergedSpecs[0], minAAmass, maxAAmass, peakTol, neighsL);
				//		denovo_LtoR(mergedSpecs[0], neighsL, seqs, minPercTopDenovo);

				// Sequence the filtered MS4 spectrum
				if(params.paramPresent("MIN_AA_MASS") or params.paramPresent("MAX_AA_MASS")) {
					getNeighborsL(mergedSpecs[specIdx], minAAmass, maxAAmass, peakTol, neighsL);
					getNeighborsR(mergedSpecs[specIdx], minAAmass, maxAAmass, peakTol, neighsR);
				} else {
					getNeighborsL(mergedSpecs[specIdx], jumpsDenovo, peakTol, neighsL);
					getNeighborsR(mergedSpecs[specIdx], jumpsDenovo, peakTol, neighsR);
/*for(unsigned int i=0; i<mergedSpecs[specIdx].size(); i++) {
	cerr<<"Peak "<<i<<", mass "<<mergedSpecs[specIdx][i][0]<<": ";
	for(unsigned int j=0; j<neighsL[i].size(); j++)
		cerr<<"L("<<neighsL[i][j]<<","<<mergedSpecs[specIdx][neighsL[i][j]][0]<<")";
	for(unsigned int j=0; j<neighsR[i].size(); j++)
		cerr<<"R("<<neighsR[i][j]<<","<<mergedSpecs[specIdx][neighsR[i][j]][0]<<")";
	cerr<<endl;
}*/
				}
				denovo(mergedSpecs[specIdx],neighsL,neighsR,peakTol,pmTol,ionOffset,minPercTopDenovo,cur_seqsList,cur_seqScores,false);
				cout<<" got "<<cur_seqsList.size()<<" seqs\n"; cout.flush();

				float matchScore, denovoScore;
				denovoScore=-1; for(list<float>::iterator iter=cur_seqScores.begin(); iter!=cur_seqScores.end();iter++) denovoScore=max(denovoScore,*iter);
				for(unsigned int pepIdx=0; pepIdx<correctPeps.size(); pepIdx++)
					if(fabs(mergedSpecs[specIdx].parentMass-correctPeps[pepIdx].parentMass)<=2*peakTol) {
//						FindMatchPeaks(correctPeps[pepIdx],mergedSpecs[specIdx],0,peakTol,idxPep,idxSeq);
						ScoreOverlap6(correctPeps[pepIdx],mergedSpecs[specIdx],0,peakTol,idxPep,idxSeq,56);
						matchScore=0; for(unsigned int i=0; i<idxSeq.size();i++) matchScore+=mergedSpecs[specIdx][idxSeq[i]][1];
						cout<<" -.-> Correct peptide "<<pepIdx+1<<" explains "<<idxPep.size()<<"/"<<correctPeps[pepIdx].size()<<" peaks and "<<matchScore<<" / "<<100*matchScore/denovoScore<<"% of top score in merged spec "<<specIdx+1<<" (shift "<<specSet[0].parentMass-mergedSpecs[specIdx].parentMass<<"), pep masses: ";

						unsigned int j=0;
						for(unsigned int i=0;i<correctPeps[pepIdx].size();i++) {
							while(j<idxPep.size() and idxPep[j]<i) j++;
							if(j<idxPep.size() and idxPep[j]==i) cout<<"."; else cout<<"X";
						}
						cout<<endl;
					}

				seqsList.splice(seqsList.end(),cur_seqsList);      cur_seqsList.clear();
				seqScores.splice(seqScores.end(),cur_seqScores);   cur_seqScores.clear();
			}
		}
	}
	cout<<"Finished de novo sequencing with "<<seqsList.size()<<" sequences\n";

	// ********************************************************************************
	//   Output de novo seqs
	// ********************************************************************************
	list<float>::iterator iterScores=seqScores.begin();    specIdx=0;
	denovoSpecs.resize(seqsList.size());    seqsPairs.resize(seqsList.size());
	for(list<Spectrum>::iterator iterSeqs=seqsList.begin(); iterSeqs!=seqsList.end(); iterSeqs++, iterScores++, specIdx++) {
		denovoSpecs[specIdx] = *iterSeqs;   seqsPairs[specIdx][0]=*iterScores;   seqsPairs[specIdx][1]=(float)specIdx;
	}
//	denovoSpecs.SaveSpecSet_pkl("specs_denovo.pkl");
	denovoSpecs.SaveSpecSet_pklbin("specs_denovo.pklbin");

	// ********************************************************************************
	//   Estimate accuracy of de novo sequencing
	// ********************************************************************************
	vector<vector<float> > denovo_stats(denovoSpecs.size());  // # de novo peaks, % matched correct peaks, de novo score, index of correct peptide
	vector<short> seqType(denovoSpecs.size());  // types: 0 (incorrect), 1 (correct), 2 (correct but shorter), 3 (correct but longer)
	int numMatched, bestPep;
	for(specIdx=0; specIdx<denovoSpecs.size(); specIdx++) {
		denovo_stats[specIdx].resize(3);  // Num peaks, % matched correct peaks, de novo score, index of best-matched correct peptide
		denovo_stats[specIdx][0]=denovoSpecs[specIdx].size();   denovo_stats[specIdx][2]=seqsPairs[specIdx][0];
		denovo_stats[specIdx][1]=-1;
		for(unsigned int pepIdx=0; pepIdx<correctPeps.size(); pepIdx++) {
			FindMatchPeaks(correctPeps[pepIdx],denovoSpecs[specIdx],0,peakTol,idxPep,idxSeq);
			float percMatched = float(idxPep.size())/float(correctPeps[pepIdx].size());
			if(percMatched > denovo_stats[specIdx][1])
				{ numMatched = (int)idxPep.size();   bestPep=(int)pepIdx;   denovo_stats[specIdx][1]=percMatched; }
		}

		seqType[specIdx]=0;   if(correctPeps.size()==0) { cerr<<"Skipping...\n"; continue; }
		if((int)denovo_stats[specIdx][0]==numMatched) {
			if(numMatched==correctPeps[bestPep].size())
				seqType[specIdx]=1; else seqType[specIdx]=2;
		} else if(numMatched==correctPeps[bestPep].size())
			        seqType[specIdx]=3;
	}
	Save_binArray("denovo_stats.bin",denovo_stats);
	Save_binArray("seqType.bin",seqType);

	sort(seqsPairs.begin(),seqsPairs.end());
	Save_binArray("indices_expint.bin", seqsPairs);  // Output the indices (pos 1) of the denovo sequences sorted by increasing score (pos 0)

	unsigned int bestCorrect=0, bestShort=0, bestLong=0;
	unsigned int countCorrect=0, countShort=0, countLong=0;
	for(unsigned int rank=1; rank<=seqsPairs.size(); rank++) {
		specIdx = (unsigned int)seqsPairs[seqsPairs.size()-rank][1];
		if(seqType[specIdx]>0) {
			if(seqType[specIdx]==1) { countCorrect++; if(bestCorrect==0) bestCorrect=rank; }
			if(seqType[specIdx]==2) { countShort++;   if(bestShort==0) bestShort=rank; }
			if(seqType[specIdx]==3) { countLong++;    if(bestLong==0) bestLong=rank; }
//			if(bestCorrect>0 and bestShort>0 and bestLong>0) break;
		}
	}
	cout<<" --- Ranks of best correct/short/longer de novo sequences: "<<bestCorrect<<" / "<<bestShort<<" / "<<bestLong<<"\n";
	cout<<" --- Number of correct/short/longer de novo sequences: "<<countCorrect<<" / "<<countShort<<" / "<<countLong<<"\n";

//	return 0;

	// Add the shift to every mass in every de novo reconstruction
	// Add the shift as a peak to every de novo reconstruction
	Spectrum newPeak; newPeak.resize(1); newPeak[0].set(ionOffset,1.0);
	SpecSet cyclicDenovo(denovoSpecs.size());
//	SpecSet ms3spec(1);   ms3spec[0]=specSet[0];   ms3spec.SaveSpecSet_pkl("ms3spec.pkl");
	unsigned int numPeaks;
	for(specIdx=0; specIdx<denovoSpecs.size(); specIdx++) {
		numPeaks = denovoSpecs[specIdx].size();
		denovoSpecs[specIdx].addZPMpeaks(peakTol,ionOffset,false,false);
		if(not params.paramPresent("INPUT_CONSENSUS") and not params.paramPresent("INPUT_DENOVO")) {
			shift = specSet[0].parentMass - denovoSpecs[specIdx].parentMass;
			for(peakIdx=0; peakIdx<numPeaks; peakIdx++) denovoSpecs[specIdx][peakIdx][0]+=shift;
			denovoSpecs[specIdx].parentMass += shift;
		}

//		denovoSpecs[specIdx].mergePeakList(newPeak.peakList, 0);
		Cyclize(denovoSpecs[specIdx],cyclicDenovo[specIdx],ionOffset,peakTol,resolution);
		cyclicDenovo[specIdx].addZPMpeaks(peakTol,ionOffset,false,false);
//		continue;

		// Replicate each the de novo peptide after itself to allow for matching MS4 spectra that cross start/end boundaries
		denovoSpecs[specIdx].resize(2*numPeaks);
		for(peakIdx=numPeaks; peakIdx<denovoSpecs[specIdx].size(); peakIdx++) {
			denovoSpecs[specIdx][peakIdx]=denovoSpecs[specIdx][peakIdx-numPeaks];
			denovoSpecs[specIdx][peakIdx][0] += (denovoSpecs[specIdx].parentMass-ionOffset);
		}
		denovoSpecs[specIdx].parentMass += (denovoSpecs[specIdx].parentMass-ionOffset);

		denovoSpecs[specIdx].mergePeakList(newPeak.peakList, 0);
	}
	denovoSpecs.SaveSpecSet_pklbin("denovo_db.pklbin");

	cout<<"Searching de novo sequences against the MS4 spectra...";   cout.flush();

	// Choose best de novo reconstruction by matching against MS3 spectra
//	vector<TwoValues<float> > scores(0);
	vector<TwoValues<float> > scores(4), scoresEmpty(0);
	// Dancik scores with probabilities 0.8,0.6,0.4,0.3, p(noise)=0.1
	scores[0][0]=2.1;   scores[0][1]=-2.2;   // b-ion     (p=0.8)
	scores[1][0]=1.9;   scores[1][1]=-1.1;   // b-ion C13 isotope (p=0.7)
	scores[2][0]=1.4;   scores[2][1]=-0.6;   // b-ion-H2O (p=0.6))
//	scores[3][0]=1.1;   scores[3][1]=-0.4;   // b-ion-NH3 (p=0.3)
	scores[3][0]=1.9;   scores[3][1]=-1.1;   // b-ion-NH3 (p=0.7)

	// Dancik scores with p(noise)=0.05
/*	scores[0][0]=2.6;   scores[0][1]=-1.2;   // b-ion     (p=0.7)
	scores[1][0]=2.0;   scores[1][1]=-0.5;   // b-ion C13 isotope (p=0.4)
	scores[2][0]=1.6;   scores[2][1]=-0.23;   // b-ion-H2O (p=0.25))
	scores[3][0]=1.4;   scores[3][1]=-0.17;   // b-ion-NH3 (p=0.2)
*/
	// Num matched peaks
//	scores[0][0]=1;   scores[0][1]=0;   // b-ion
//	scores[1][0]=0;   scores[1][1]=0;   // b-ion C13 isotope
//	scores[2][0]=0;   scores[2][1]=0;   // b-ion-H2O
//	scores[3][0]=0;   scores[3][1]=0;   // b-ion-NH3


//	vector<float> ionOffsets(1);
	vector<float> ionOffsets(4);
	ionOffsets[0]=0;
	ionOffsets[1]=AAJumps::massHion;
	ionOffsets[2]=-AAJumps::massH2O;
	ionOffsets[3]=-AAJumps::massNH3;
	vector<TwoValues<float> > pepMatchScores;
	vector<list<unsigned int> > pepMatchedSpecs;
	SpecSet querySpecs;   querySpecs.resize(specSet.size()-1);   //querySpecs[0] = specSet[1];
	for(specIdx=0; specIdx<querySpecs.size(); specIdx++) querySpecs[specIdx]=specSet[specIdx+1];

	// Collect match-scores statistics for discriminant analysis; one spectrum per line.
	//  V col 0 - (MS3) Percentage of explained intensity
	//  . col 1 - (MS3) Percentage of expected true fragments found in the spectrum
	//  V col 2 - (MS3) Average Dancik score
	//  V col 3 - (MSn) Average percentage of explained intensity (over all matched spectra)
	//  V col 4 - (MSn) Average Dancik scores (approximate, over all matched spectra)
	//  V col 5 - (MSn) Number of matched MSn spectra
	//  . col 6 - (MSn) Number of peaks in all matched MSn spectra
	//  V col 7 - Length of de novo sequence
	vector<vector<float> > matchStats(denovoSpecs.size());
	for(specIdx=0; specIdx<matchStats.size(); specIdx++) {
		matchStats[specIdx].resize(8);
		for(unsigned int i=0; i<8; i++) matchStats[specIdx][i]=0;
		matchStats[specIdx][7]=denovoSpecs[specIdx].size();
	}

	// Match to MS3 spectrum
	// Cyclize each denovo sequence
	SpecSet ms3spec;   ms3spec.loadPklBin("ms3ori.pklbin");
//	for(peakIdx=0;peakIdx<ms3spec[0].size();peakIdx++) ms3spec[0][peakIdx][0]*=0.9995;  // Convert to near-integer masses
	cyclicDenovo.SaveSpecSet_pklbin("cyclic_denovo.pklbin");
	cout<<" - done cyclizing de novo sequences\n";
/*
	// Score match of MS3 spectrum to each cyclized denovo sequence
	cout<<" - matching cyclic de novo sequences to ms3 spectrum (Dancik)..."; cout.flush();
	MatchSpecsToPeps(ms3spec, cyclicDenovo, peakTol, ionOffset, ionOffset, 0, scores, ionOffsets, pepMatchScores, pepMatchedSpecs,false);
	for(unsigned int matchIdx=0; matchIdx<pepMatchScores.size(); matchIdx++) {
		specIdx = (unsigned int)round(pepMatchScores[matchIdx][1]);
		matchStats[specIdx][2] = pepMatchScores[matchIdx][0] / matchStats[specIdx][7];
	}
	cout<<" - done\n";

	cout<<" - matching cyclic de novo sequences to ms3 spectrum (Intensity)..."; cout.flush();
	MatchSpecsToPeps(ms3spec, cyclicDenovo, peakTol, ionOffset, ionOffset, 0, scoresEmpty, ionOffsets, pepMatchScores, pepMatchedSpecs,false);
	for(unsigned int matchIdx=0; matchIdx<pepMatchScores.size(); matchIdx++) {
		specIdx = (unsigned int)round(pepMatchScores[matchIdx][1]);
		matchStats[specIdx][0] = pepMatchScores[matchIdx][0];
	}
	cout<<" - done\n";
*/
	// Choose best de novo reconstruction by matching against MSn spectra
	vector<float> scores_msn(3);
	scores_msn[0]=4;   // Score if matched b-ion
	scores_msn[1]=-3;  // Penalty if missing b-ion
	scores_msn[2]=-2;  // Penalty if noise peak

	// Match to MSn spectra
	cout<<" - matching de novo sequences to MSn spectra (Dancik)..."; cout.flush();
	MatchSpecsToPeps(querySpecs, denovoSpecs, msn_peakTol, ionOffset, ionOffset, 0, scores, ionOffsets, pepMatchScores, pepMatchedSpecs,true);
	for(unsigned int matchIdx=0; matchIdx<pepMatchScores.size(); matchIdx++) {
		specIdx = (unsigned int)round(pepMatchScores[matchIdx][1]);
		matchStats[specIdx][5] = pepMatchedSpecs[specIdx].size();
		if(matchStats[specIdx][5]==0) matchStats[specIdx][4] = 0;
		else matchStats[specIdx][4] = pepMatchScores[matchIdx][0] / (matchStats[specIdx][7]*matchStats[specIdx][5]);
	}
	cout<<" - done\n";
/*
	cout<<" - matching de novo sequences to MSn spectra (Intensity)..."; cout.flush();
	MatchSpecsToPeps(querySpecs, denovoSpecs, msn_peakTol, ionOffset, ionOffset, 0, scoresEmpty, ionOffsets, pepMatchScores, pepMatchedSpecs,true);
	for(unsigned int matchIdx=0; matchIdx<pepMatchScores.size(); matchIdx++) {
		specIdx = (unsigned int)round(pepMatchScores[matchIdx][1]);
		if(matchStats[specIdx][5]==0) matchStats[specIdx][3] = 0;
		else matchStats[specIdx][3] = pepMatchScores[matchIdx][0] / matchStats[specIdx][5];
	}
	cout<<" - done\n";
*/
//	MatchSpecsToPeps_old(querySpecs, denovoSpecs, msn_peakTol, ionOffset, ionOffset, scores_msn, ionOffsets, pepMatchScores, pepMatchedSpecs,false);

/*for(unsigned int i=0; i<pepMatchScores.size();i++) {
	specIdx = (unsigned int)round(pepMatchScores[i][1]);
	if(pepMatchedSpecs[specIdx].size()>0) {
		cerr<<"Peptide "<<specIdx+1<<" got score "<<pepMatchScores[i][0]<<" with matched spectra";
		for(list<unsigned int>::iterator pmsIter=pepMatchedSpecs[specIdx].begin();pmsIter!=pepMatchedSpecs[specIdx].end();pmsIter++)
			cerr<<" "<<*pmsIter;
		cerr<<"\n";
	}
}*/

	Save_binArray("matches_stats.bin", matchStats);  // Output the matches' statistics

	Save_binArray("matches_expint.bin", pepMatchScores);  // Output the indices of the denovo sequences sorted by increasing match score
	Save_binListArray<unsigned int,list<unsigned int>,list<unsigned int>::iterator >("matches_specidx.bla",pepMatchedSpecs);

	cout<<"done\n";   cerr.flush();

	/***************************************************************************
	 *  Output of miscellaneous information
	**************************************************************************/

    ofstream output;
	if(params.paramPresent("OUTPUT_ALIGNS")) {
	    output.open(alignsFN, ios::binary);  output << results.size() << endl;
		for(list<Results_PA>::iterator i=results.begin(); i!=results.end(); i++)
			(*i).output(output,';');
	    output.close();
	}

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

/*	if(params.paramPresent("OUTPUT_MATCHED_PEAKS")) {
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
	}*/

	filename = params.getValue("OUTPUT_STATS");
	if(params.paramPresent("OUTPUT_STATS")) {
		char *fnBuffer = (char *)malloc(strlen(filename)+20);
		sprintf(fnBuffer,"%s_specStats.bin",filename);
		Save_binArray(fnBuffer, specStats);
		sprintf(fnBuffer,"%s_alignStats.bin",filename);
		Save_binArray(fnBuffer, alignStats);
		free(fnBuffer);
	}

    return(0);
}

void Cyclize(Spectrum &in, Spectrum &out, float ionOffset, float mergeTol, float resolution) {
	out.copyNP(in);
	SpecSet variants(in.size());
	unsigned int pivot;
	for(pivot=0; pivot<in.size(); pivot++) {
		variants[pivot]=in;
		variants[pivot].rotate(ionOffset-in[pivot][0],ionOffset);  // Assumes no cterm H2O
	}
	MergeSpecs(variants,mergeTol,resolution,0,out);
}
