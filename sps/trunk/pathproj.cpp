#include "inputParams.h"
#include "alignment_scoring.h"
#include "spectral_pairs.h"
#include "batch.h"
#include "filters.h"
#include "abruijn.h"
#include "graph.h"
#include "setmerger.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;

int main(int argc, char **argv){
    // Get input parameters and check minimum parameter set
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("pathproj.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "pathproj.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs;   paramStrs.resize(6);
	paramStrs[0] = "INPUT_SPECS_PKLBIN";
	paramStrs[1] = "INPUT_ALIGNS";
	paramStrs[2] = "INPUT_ANNOTATED";
	paramStrs[3] = "OUTPUT_ANNOTINFO";
	paramStrs[4] = "MIN_PERC_EXPINT";
	paramStrs[5] = "MIN_PERC_TP";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"pathproj.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_SPECS_PKLBIN, INPUT_ALIGNS, INPUT_ANNOTATED, MIN_PERC_EXPINT, MIN_PERC_TP\n";
		return -1;
	}

	const char *specSetFN = params.getValue("INPUT_SPECS_PKLBIN");
	const char *alignsFN = params.getValue("INPUT_ALIGNS");

	const char *annotatedFN = params.getValue("INPUT_ANNOTATED");
	const char *annotinfoFN = params.getValue("OUTPUT_ANNOTINFO");
	float minExpInt = (float) params.getValueDouble("MIN_PERC_EXPINT");
	float minTP = (float) params.getValueDouble("MIN_PERC_TP");
	float peakTol = params.paramPresent("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	unsigned int minMatchedPeaks = params.paramPresent("MIN_MATCHED_PEAKS")?(int) params.getValueInt("MIN_MATCHED_PEAKS"):0;
/*	int maxAAjump = params.getValue("MAX_AA_JUMP")?params.getValueInt("MAX_AA_JUMP"):0;
	float pmTol = params.getValue("TOLERANCE_PM")?(float) params.getValueDouble("TOLERANCE_PM"):1;
	char *modPosFN = params.getValue("OUTPUT_MODPOS");
	float penalty_ptm = (float) params.getValueDouble("PENALTY_PTM");
	float penalty_sameVert = (float) params.getValueDouble("PENALTY_SAME_VERTEX");
*/
    SpecSet specSet;
    if (specSet.LoadSpecSet_pklbin(specSetFN)<=0) { cerr << "Error reading "<<specSetFN<<"!\n"; return -1; }
    else  cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";
    unsigned int numSpecs = specSet.size();
    SpecSet specSetRev(numSpecs),   // Reversed versions of every proximately-annotated spectrum
			specSet_proj(numSpecs); // Final projected version for identified/projected spectra

	vector<Results_ASP> aligns;
	if (params.paramPresent("INPUT_ALIGNS")) {
		if (Load_results_bin(alignsFN, aligns)==0) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
	    else  cout << "Loading aligns complete. Num pairs: " << aligns.size() << "\n";
	}

	SpecSet matchedPeaks(min((unsigned int)specSet.size(),(unsigned int)aligns.size()));  // Indices of matched peaks for propagated annotations
	unsigned int numMatchedPairs = 0;       // Number of pairs used for propagation

	vector<vector<int> > matchInfo(numSpecs);  // Proximate level of annotation (col.1), index of the spectrum where the annotation came from (col.2),
	                                           //  index of the first peak to the right of the annotation (col.3),
	                                           //  100*percentage of explained score for the selected annotation (col.4)
	                                           //  100*percentage of matched peptide breakpoints (either b or y but not both) (col.5)
	for (unsigned int i=0; i<numSpecs; i++) { matchInfo[i].resize(5); for (unsigned int j=0; j<5; j++) matchInfo[i][j]=0; }

	vector<bool> annotated(numSpecs);   vector<vector<short> > annotatedFromFile;
	int annotatedCount=0;
	if (params.paramPresent("INPUT_ANNOTATED")) {
		Load_binArray(annotatedFN, annotatedFromFile);
		if(annotatedFromFile.size()!=numSpecs or annotatedFromFile[0].size()!=2) { cerr<<"ERROR reading "<<annotatedFN<<"!\n"; return -1; }
		for(unsigned int i=0; i<numSpecs; i++) {
			specSet_proj[i]=specSet[i];
			if(annotatedFromFile[i].size()!=2) { cerr<<"ERROR reading "<<annotatedFN<<" at entry number "<<i<<"!\n"; return -1; }
			annotated[i]=(annotatedFromFile[i][0]>=minExpInt-0.0001 and annotatedFromFile[i][1]>=minTP-0.0001);
			if(annotated[i]) {
				matchInfo[i][3]=annotatedFromFile[i][0];
				matchInfo[i][4]=annotatedFromFile[i][1];
//				specSet_proj[i]=specSet[i];
				annotatedCount++;
			}
		}
		cout<<"Number of annotated spectra: "<<annotatedCount<<endl;
	}
#ifdef DEBUG
	ofstream debug("dekel_align_debug.txt");
#endif

	vector<list<int> > neighs(numSpecs);
	for(unsigned int i=0; i<aligns.size(); i++) { neighs[aligns[i].spec1].push_back(aligns[i].spec2); neighs[aligns[i].spec2].push_back(aligns[i].spec1); }

	vector<short> scheduled(numSpecs); for(unsigned int i=0; i<numSpecs; i++) scheduled[i]=0;  // Indicates the iteration on which a spectrum is to be processed
	vector<list<int> > annNeighs(numSpecs);  // Annotated neighbors to consider when annotating the current spectrum
	list<int> toProcess; toProcess.clear();  // Spectra eligible for annotation on the next iteration
	int annIter=1,              // Number of hops between current iteration and initially annotated spectra
	    lastInIteration=0,      // Last spectrum index in toProcess with the current annIter
	    lastInNextIteration=-1; // Last spectrum index in toProcess with the next annIter
	for(unsigned int specIdx=0; specIdx<numSpecs; specIdx++)
		if(annotated[specIdx])
			for(list<int>::iterator curNeigh=neighs[specIdx].begin(); curNeigh!=neighs[specIdx].end(); curNeigh++)
				if(not annotated[*curNeigh]) {
					if(scheduled[*curNeigh]==0) { scheduled[*curNeigh]=1; annNeighs[*curNeigh].clear(); toProcess.push_back(*curNeigh); lastInIteration = *curNeigh; }
					annNeighs[*curNeigh].push_back(specIdx);
				}

cerr<<"(before loop): 2016 sizes = ("<<specSet[2016].size()<<","<<specSet_proj[2016].size()<<"), 720 sizes = ("<<specSet[720].size()<<","<<specSet_proj[720].size()<<"), 2751 sizes = ("<<specSet[2751].size()<<","<<specSet_proj[2751].size()<<")\n";
cerr.flush();

	for(list<int>::iterator curSpec=toProcess.begin(); curSpec!=toProcess.end(); curSpec++) {
		Spectrum curProj, bestProj;
		float bestScore=0, curScore=0;
		list<int> neighToProcess;     // Index of the node being annotated (propagated onto); data structure is a side effect of using ProjectSpectrum: this always contains the single element *curSpec
		list<int> curDeltas, bestDeltas, curBestDeltas;
		vector<TwoValues<int> > matchesSpec;
		list<int>::iterator curNeigh;
		int bestNeigh, bestMatchedPeakCount=0;
		float maxSpecScore=0; for(unsigned int i=0; i<specSet[*curSpec].size(); i++) maxSpecScore+=specSet[*curSpec][i][1];

		bestNeigh = annNeighs[*curSpec].empty()?-1:annNeighs[*curSpec].front();  // initialize with first neighbor

		// Find the annotated neighbor with the best projection onto the current spectrum
		specSet[*curSpec].reverse(0, &specSetRev[*curSpec]);
		curDeltas.clear();   bestDeltas.clear();   curBestDeltas.clear();
		neighToProcess.clear();   neighToProcess.push_front(*curSpec);
		matchedPeaks[numMatchedPairs].resize(0);

		/* ProjectSpectrum is used to project the non-annotated spectrum (*curSpec) onto each
		 * annotated neighbor (curNeigh). At the end of the loop bestProj contains the modified
		 * spectrum of the annotated neighbor best-matched to *curSpec (masses modified by delta
		 * and peak scores increased by matched peaks in *curSpec)
		 */
		for(curNeigh=annNeighs[*curSpec].begin(); curNeigh!=annNeighs[*curSpec].end(); curNeigh++) {
			curDeltas.clear();   curScore=0;   matchesSpec.resize(0);
			ProjectSpectrum(specSet, specSet[*curNeigh], neighToProcess, curDeltas, curScore, curBestDeltas, peakTol, &curProj, &matchesSpec, minMatchedPeaks);
			curScore=0; for(unsigned int p=0; p<matchesSpec.size(); p++) curScore+=specSet[*curSpec][matchesSpec[p][1]][1];
			if(matchesSpec.size()>0 and curScore>bestScore) {
				bestScore=curScore; bestDeltas.clear(); bestDeltas.splice(bestDeltas.end(),curBestDeltas); bestNeigh = *curNeigh; bestProj = curProj; bestMatchedPeakCount = matchesSpec.size();
				specSet_proj[*curSpec] = specSet[*curSpec];
				specSet_proj[*curSpec].resize(matchesSpec.size());
				matchedPeaks[numMatchedPairs].resize(matchesSpec.size()); // Update sets of peaks matched by projections
				for(unsigned int p=0; p<matchesSpec.size(); p++) {
					specSet_proj[*curSpec][p] = specSet[*curSpec][matchesSpec[p][1]];
					matchedPeaks[numMatchedPairs][p].set(p,matchesSpec[p][0]);
//					matchedPeaks[numMatchedPairs][p].set(matchesSpec[p][1],matchesSpec[p][0]);
				}
			}

			curDeltas.clear();   curScore=0;   matchesSpec.resize(0);
			ProjectSpectrum(specSetRev, specSet[*curNeigh], neighToProcess, curDeltas, curScore, curBestDeltas, peakTol, &curProj, &matchesSpec, minMatchedPeaks);
			curScore=0; for(unsigned int p=0; p<matchesSpec.size(); p++) curScore+=specSetRev[*curSpec][matchesSpec[p][1]][1];
			if(matchesSpec.size()>0 and curScore>bestScore) {
				bestScore=curScore; bestDeltas.clear(); bestDeltas.splice(bestDeltas.end(),curBestDeltas); bestNeigh = *curNeigh; bestProj = curProj; bestMatchedPeakCount = matchesSpec.size();
				specSet_proj[*curSpec] = specSetRev[*curSpec];
				specSet_proj[*curSpec].resize(matchesSpec.size());
				matchedPeaks[numMatchedPairs].resize(matchesSpec.size()); // Update sets of peaks matched by projections
				for(unsigned int p=0; p<matchesSpec.size(); p++) {
					specSet_proj[*curSpec][p] = specSetRev[*curSpec][matchesSpec[p][1]];
					matchedPeaks[numMatchedPairs][p].set(p,matchesSpec[p][0]);
//					matchedPeaks[numMatchedPairs][p].set(matchesSpec[p][1],matchesSpec[p][0]);
				}
			}
		}

		scheduled[*curSpec] = false;
		if(bestProj.size()>0) {
			annotated[*curSpec] = true;
//			specSet[*curSpec] = bestProj;   // Propagate the annotation, keep all annotation peaks even if not from this spectrum (accumulated evidence)
			specSet[*curSpec] = specSet_proj[*curSpec];   // Propagate the annotation, keep only annotated peaks
if(*curSpec==720 or *curSpec==2751) {
	cerr<<"(*curSpec=="<<*curSpec<<"): 2016 sizes = ("<<specSet[2016].size()<<","<<specSet_proj[2016].size()<<"), 720 sizes = ("<<specSet[720].size()<<","<<specSet_proj[720].size()<<"), 2751 sizes = ("<<specSet[2751].size()<<","<<specSet_proj[2751].size()<<")\n";
	cerr<<"-- bestProj:\n"; bestProj.output(cerr);
	cerr<<"-- specSet:\n"; specSet[*curSpec].output(cerr);
	cerr<<"-- specSet_proj:\n"; specSet_proj[*curSpec].output(cerr);
	cerr.flush();
}

			// Statistics for best projection
			matchInfo[*curSpec][0] = annIter;              matchInfo[*curSpec][1] = bestNeigh;
			matchInfo[*curSpec][2] = bestDeltas.empty()?0:bestDeltas.front();
			matchInfo[*curSpec][3] = maxSpecScore>0?(int)round(10000*(bestScore/maxSpecScore)):0;
			matchInfo[*curSpec][4] = (int)round(10000*(((float)bestMatchedPeakCount)/((float)specSet[bestNeigh].size())));

			// Update set of pairs used for projections
			aligns[numMatchedPairs].spec1 = *curSpec;
			aligns[numMatchedPairs].spec2 = bestNeigh;
			aligns[numMatchedPairs].shift1 = matchInfo[*curSpec][2];
			aligns[numMatchedPairs].score1 = bestScore;
			aligns[numMatchedPairs].score2 = -1;   // missing value
			numMatchedPairs++;

			// Check if any of the current spectrum's neighbors needs processing
			for(curNeigh=neighs[*curSpec].begin(); curNeigh!=neighs[*curSpec].end(); curNeigh++) {
				if(not annotated[*curNeigh]) {
					if(scheduled[*curNeigh]==0) { scheduled[*curNeigh]=annIter+1; annNeighs[*curNeigh].clear(); toProcess.push_back(*curNeigh); lastInNextIteration = *curNeigh; }
					if(scheduled[*curNeigh]>annIter) annNeighs[*curNeigh].push_back(*curSpec);
				}
			}
		}

		if((*curSpec)==lastInIteration) { annIter++; lastInIteration = lastInNextIteration; lastInNextIteration=-1; }
	}

cerr<<"(after loop): 2016 sizes = ("<<specSet[2016].size()<<","<<specSet_proj[2016].size()<<"), 720 sizes = ("<<specSet[720].size()<<","<<specSet_proj[720].size()<<"), 2751 sizes = ("<<specSet[2751].size()<<","<<specSet_proj[2751].size()<<")\n";
cerr.flush();

	for(unsigned int i=0; i<numSpecs; i++) {
		if(matchInfo[i][0]>0) matchInfo[i][1]++; // Convert indices to 1-based instead of 0-based
//		if(not annotated[i]) specSet_proj[i].resize(0);
		if(specSet_proj[i].size()==0) specSet_proj[i]=specSet[i];  // Recover original spectra for unsuccessful projections
	}
	Save_binArray(annotinfoFN, matchInfo);

	// Output projected peptide annotations
	if (params.paramPresent("OUTPUT_SPECS_PROJ"))
		specSet_proj.SaveSpecSet_pklbin(params.getValue("OUTPUT_SPECS_PROJ"));

	// Output sets of peaks matched by projection
	if (params.paramPresent("OUTPUT_SPECS_MATCHIDX")) {
		matchedPeaks.resize(numMatchedPairs);
		matchedPeaks.SaveSpecSet_pklbin(params.getValue("OUTPUT_SPECS_MATCHIDX"));
	}

	// Output set of pairs used for projections
	if (params.paramPresent("OUTPUT_ALIGNS")) {
		aligns.resize(numMatchedPairs);
		Save_results_bin(params.getValue("OUTPUT_ALIGNS"),numMatchedPairs,aligns.begin());
	}

#ifdef DEBUG
	debug.close();
#endif

	return(0);

}
