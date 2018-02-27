#include "msn.h"
#include "alignment_scoring.h"
#include "denovo.h"
#include "aminoacid.h"
#include <vector>
#include <cmath>

namespace specnets
{

	void get_ms3_sets_max(SpecSet &specs, vector<Results_ASP> &aligns, vector<list<int> > &sets, int maxNumMS3) {
		unsigned int pivot, counter=0;
		int parent;
		vector<list<pair<float, int> > > setsAll;   pair<float, int> curPair;

		setsAll.resize(specs.size()); for(pivot=0; pivot<setsAll.size(); pivot++) setsAll[pivot].clear();
		for(pivot=0; pivot<aligns.size(); pivot++) {
			if(specs[aligns[pivot].spec1].parentMass > specs[aligns[pivot].spec2].parentMass) { parent=aligns[pivot].spec1; curPair.first=aligns[pivot].score1; curPair.second=aligns[pivot].spec2; }
			else { parent=aligns[pivot].spec2; curPair.first=aligns[pivot].score2; curPair.second=aligns[pivot].spec1; }
			if(specs[parent].parentMass-specs[curPair.second].parentMass<50) continue;
			setsAll[parent].push_back(curPair);   if(setsAll[parent].size()>counter) counter=setsAll[parent].size();
		}
		if(maxNumMS3<=0) maxNumMS3 = counter+1;  // Use as many child spectra as possible

		list<pair<float, int> >::reverse_iterator iter;   sets.resize(specs.size());
		for(pivot=0; pivot<setsAll.size(); pivot++) {
			sets[pivot].clear();             if(setsAll[pivot].empty()) continue;
			sets[pivot].push_back(pivot);    setsAll[pivot].sort();
			for(counter=0, iter=setsAll[pivot].rbegin(); (int)counter<=maxNumMS3 and iter!=setsAll[pivot].rend(); counter++, iter++)
				sets[pivot].push_back(iter->second);
		}
	}

	struct msn_SpecInfo { float parentMass, score; int index; };
	bool msn_cmp_SpecInfo(msn_SpecInfo &i1, msn_SpecInfo &i2)
	{ if(i1.parentMass==i2.parentMass) return i1.score<i2.score; else return i1.parentMass<i2.parentMass; }

	void get_ms3_sets_sparse(SpecSet &specs, vector<Results_ASP> &aligns,
														 vector<list<int> > &sets, float pmTol, unsigned int maxNumMS3) {
		unsigned int pivot, setIdx;
		int parent;
		vector<list<msn_SpecInfo> > setsAll;   msn_SpecInfo curInfo;

		setsAll.resize(specs.size()); for(pivot=0; pivot<setsAll.size(); pivot++) setsAll[pivot].clear();
		for(pivot=0; pivot<aligns.size(); pivot++) {
			if(specs[aligns[pivot].spec1].parentMass > specs[aligns[pivot].spec2].parentMass) { parent=aligns[pivot].spec1; curInfo.parentMass=specs[aligns[pivot].spec2].parentMass; curInfo.score=aligns[pivot].score1; curInfo.index=aligns[pivot].spec2; }
			else { parent=aligns[pivot].spec2; curInfo.parentMass=specs[aligns[pivot].spec1].parentMass; curInfo.score=aligns[pivot].score2; curInfo.index=aligns[pivot].spec1; }
			if(specs[parent].parentMass-specs[curInfo.index].parentMass<50) continue;
			setsAll[parent].push_back(curInfo);
		}

		sets.resize(specs.size());
		Spectrum tmpSpec;   vector<int> idxAll, idxMatched1, idxMatched2, specIndices;
		setIdx = 0;
		for(pivot=0; pivot<setsAll.size(); pivot++) {
			if(setsAll[pivot].empty()) continue;
			int numElems = setsAll[pivot].size(), curElem=0;

			setsAll[pivot].sort(msn_cmp_SpecInfo);
			tmpSpec.resize(numElems);   idxAll.resize(numElems);   specIndices.resize(numElems);
			idxMatched1.resize(0);      idxMatched2.resize(0);
			for(list<msn_SpecInfo>::iterator iter=setsAll[pivot].begin(); iter!=setsAll[pivot].end(); iter++) {
				tmpSpec[curElem][0] = iter->parentMass;    tmpSpec[curElem][1] = iter->score;
				specIndices[curElem] = iter->index;
				idxAll[curElem]=curElem;                   curElem++;
			}

			ScoreOverlap6mp(tmpSpec, idxAll, tmpSpec, idxAll, 0, pmTol, idxMatched1, idxMatched2, 50);

			sets[setIdx].clear();
			for(unsigned int i=0; i<idxMatched1.size(); i++) sets[setIdx].push_back(specIndices[idxMatched1[i]]);

			if(maxNumMS3>0 and maxNumMS3<sets[setIdx].size()) { sets[setIdx].clear(); continue; }

			sets[setIdx].push_front(pivot); setIdx++;
		}
		sets.resize(setIdx);
	}

	/*
	 * TODO:
	 *   - Why are matchesSame/matchesDiff necessary? Can't these be replaced by a
	 *       minNumMatches parameter to MergeSpecs?
	*/

	void MergeIntoReference(Spectrum &refSpec, Spectrum &newPeaks, float mergeTol,
									 float resolution, Spectrum &result);

		/* Inner workings, K MS3 spectra
		 *  - Computes all the possible alignments between MS2/MS3, MS3/MS3 spectra: LL, LR, RL, RR
		 *      and saves the resulting lists of matched peaks in matchesSame, matchesDiff
		 *  - For each of the 2^K possible configurations:
		 *    1) Determine the subsets of peaks that were matched at least once to one other spectrum in
		 *         this configuration (fills peakUsed by iterating over matchesSame/matchesDiff)
		 *    2) Call MergeSpecs to generate a consensus spectrum for all internal peaks
		 *    3) Instantiate the endpoints set and spectrum
		 *    4) Merge internal peaks into endpoints if pmTol<peakTol and vice-versa if pmTol>peakTol
		 *    5) Run denovo on the final spectrum assuming that all peaks are b/y (2 separate denovo runs),
		 *         denovo should ignore the score of all symmetric peaks
		 *    6) Keep the final spectrum with the best denovo score, in percentage of the total spectrum score
		*/
	void denovo_ms3(SpecSet &specsL, bool allSpecsMS2, float ionOffset, float peakTol, float pmTol, float resolution,
										float mpPenalty, MS2ScoringModel &ms2model, MS3ScoringModel &ms3modelB, MS3ScoringModel &ms3modelY,
										AAJumps &jumps, Spectrum &consensus, list<Spectrum> *allTopSeqs, bool enforcePM,
										vector<vector<float> > *binomialScores, float minPerc) {
		SpecSet specsR(specsL.size()),      // Contains the right-aligned variants of every spectrum
				specsLnorm(specsL.size()),  // Contains the normalized left-aligned variants of every spectrum
				specsRnorm(specsL.size());  // Contains the normalized right-aligned variants of every spectrum
		vector<vector<TwoValues<vector<int> > > > matchesSame;
											 // pos[i][j], i>j: matched peaks when spectra i,j are both aligned on the left  (lower triangular)
											 // pos[i][j], i<j: matched peaks when spectra i,j are both aligned on the right (upper triangular)
												 // on both variants, [i][j][0] - list of matched peaks in spectrum i
												 //                   [i][j][1] - list of matched peaks in spectrum j
												 // i==j always empty
		vector<vector<TwoValues<vector<int> > > > matchesDiff;
											 // pos[i][j], i>j: matched peaks when spectra i,j are aligned on the right/left  (lower triangular)
											 // pos[i][j], i<j: matched peaks when spectra i,j are aligned on the left/right (upper triangular)
												 // on both variants, [i][j][0] - list of matched peaks in spectrum i
												 //                   [i][j][1] - list of matched peaks in spectrum j
												 // i==j, i==0, j==0 always empty
		vector<vector<bool> > peakUsed;       // pos[i][j]: indicates whether peak j from spectrum i was matched at least once in the current configuration
		SpecSet endpoints(specsL.size());     // Contains the endpoint peaks for each configuration
		vector<TwoValues<float> > bounds(specsL.size());  // Contains each spectrum's endpoint for the current configuration (used in MergeSpecs)
		SpecSet specsToMerge(specsL.size()),        // Contains the spectra to merge for each configuration
						specsToMergeNorm(specsL.size());    // Contains the normalized spectra to merge
		Spectrum specIP,     // Used to hold the merged internal peaks
						 specEP,     // Used to hold the merged endpoints
						 specConf,   // Used to hold the final merged spectrum for the current configuration
						 specDenovo; // Used to hold the peaks selected by denovo
		vector<bool> leftAligned(specsL.size());    // Indicates whether the i-th spectrum is left-aligned in the current configuration
		unsigned int specsIdx, specsIdx2, peakIdx,  // Iterator vars for spectra and peaks
								iterIdx,                        // Iterator var over all possible configurations
								iterCount=1<<(specsL.size()-1); // Total number of configurations to try
		float mass,        // Temporary mass values
					curScore,    // Score of the denovo interpretation of the current configuration
					totCurScore, // Maximum score of the current configuration
					bestScore=0, // Score of the denovo interpretation of the best configuration
					symOffset=2*ionOffset;  // Offset such that b_i + y_{n-i} = parent mass + symOffset

		if(ms2model.probs.size()>0) ScoreSpectrum(specsL[0],ms2model,false);
	//cerr<<"specsL[0]:\n";
	//specsL[0].output(cerr);

		// Reduce relative weight of MS3 PRM scores
		for(specsIdx=1; specsIdx<specsL.size(); specsIdx++)
			for(peakIdx=0; peakIdx<specsL[specsIdx].size(); peakIdx++) specsL[specsIdx][peakIdx][1]/=2;

	//	for(specsIdx=0; specsIdx<specsL.size(); specsIdx++)
	//		{ specsLnorm[specsIdx]=specsL[specsIdx]; specsLnorm[specsIdx].normalize(); }

		// Initialize variables
		matchesSame.resize(specsL.size());   matchesDiff.resize(specsL.size());   peakUsed.resize(specsL.size());
		for(specsIdx=0; specsIdx<specsL.size(); specsIdx++) {
			matchesSame[specsIdx].resize(specsL.size());
			matchesDiff[specsIdx].resize(specsL.size());
	//		peakUsed[specsIdx].resize(specsL[specsIdx].size());
			if(specsIdx==0) {
				peakUsed[specsIdx].resize(specsL[specsIdx].size());
				endpoints[specsIdx].resize(2);  endpoints[0][0].set(ionOffset,0);
				for(peakIdx=0; peakIdx<specsL[0].size(); peakIdx++) endpoints[0][0][1]=max(endpoints[0][0][1],specsL[0][peakIdx][1]);
				endpoints[0][1].set(specsL[0].parentMass-AAJumps::massMH+ionOffset,endpoints[0][0][1]);
				specsR[specsIdx]=specsL[specsIdx];   // specsL[0] is a special case: base MS/MS spectrum doesn't need to be right-shifted
				bounds[0].set(endpoints[0][0][0],endpoints[0][1][0]);   // The first spectrum is always aligned on the left
			} else {
				endpoints[specsIdx].resize(1);   endpoints[specsIdx][0][1]=0;
				specsR[specsIdx]=specsL[specsIdx];
				if(ms3modelB.probs.size()>0) ScoreSpectrumMS3(specsL[specsIdx],specsL[0],ms3modelB,specsL[0].parentMass-specsL[specsIdx].parentMass,false);
				if(ms3modelY.probs.size()>0) ScoreSpectrumMS3(specsR[specsIdx],specsL[0],ms3modelY,specsL[0].parentMass-specsR[specsIdx].parentMass,false);
				mass = max((float)0,specsL[0].parentMass-specsL[specsIdx].parentMass);
				for(peakIdx=0; peakIdx<specsR[specsIdx].size(); peakIdx++) {
					specsR[specsIdx][peakIdx][0]+=mass;  // Shift all peak masses in specsR
	//				endpoints[specsIdx][0][1]=max(endpoints[specsIdx][0][1],specsR[specsIdx][peakIdx][1]);  // Set endpoint score to maximum score
					endpoints[specsIdx][0][1]+=specsR[specsIdx][peakIdx][1];  // Set endpoint score to total summed score
				}
			}
		}
		for(specsIdx=0; specsIdx<specsR.size(); specsIdx++)
			{ specsRnorm[specsIdx]=specsR[specsIdx]; specsRnorm[specsIdx].normalize(); }

		// Find matched peaks between all spectral pairs
		for(specsIdx=1; specsIdx<specsL.size(); specsIdx++)
			for(specsIdx2=0; specsIdx2<specsIdx; specsIdx2++) {
				// Both matches on the left (matchesSame[i][j], i>j, lower triangular)
				FindMatchPeaksAll(specsL[specsIdx], specsL[specsIdx2], 0, peakTol, matchesSame[specsIdx][specsIdx2][0], matchesSame[specsIdx][specsIdx2][1]);
				// Both matches on the right (matchesSame[i][j], i<j, upper triangular)
				FindMatchPeaksAll(specsR[specsIdx2], specsR[specsIdx], 0, peakTol, matchesSame[specsIdx2][specsIdx][0], matchesSame[specsIdx2][specsIdx][1]);
				if(specsIdx2>0) {
					// Spectra aligned left/right (matchesDiff[i][j], i>j, lower triangular)
					FindMatchPeaksAll(specsL[specsIdx], specsR[specsIdx2], 0, peakTol, matchesDiff[specsIdx][specsIdx2][0], matchesDiff[specsIdx][specsIdx2][1]);
					// Spectra aligned right/left (matchesDiff[i][j], i<j, upper triangular)
					FindMatchPeaksAll(specsL[specsIdx2], specsR[specsIdx], 0, peakTol, matchesDiff[specsIdx2][specsIdx][0], matchesDiff[specsIdx2][specsIdx][1]);
				}
			}

		// Iterate over all possible combinations
		specsToMerge[0] = specsL[0];  // The first spectrum never changes position
		specsToMergeNorm[0] = specsLnorm[0];
		specConf.copyNP(specsL[0]);   // Initialize all non-peak values in specConf
		specIP.copyNP(specsL[0]);     // Initialize all non-peak values in specIP
		for(peakIdx=0; peakIdx<specsL[0].size(); peakIdx++) peakUsed[0][peakIdx]=true;  // Always use all peaks from the main spectrum
		for(iterIdx=0; iterIdx<iterCount; iterIdx++) {
			// Initialize leftAligned to current configuration
			for(specsIdx=1; specsIdx<specsL.size(); specsIdx++)
				leftAligned[specsIdx]=iterIdx&(1<<(specsIdx-1));

			// Reset peakUsed
	//		for(specsIdx=0; specsIdx<specsL.size(); specsIdx++)   // Only use peaks from the main spectrum if they're matched to other spectra
			for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) {   // Always use all peaks from the main spectrum
				if(leftAligned[specsIdx]) peakUsed[specsIdx].resize(specsL[specsIdx].size()); else peakUsed[specsIdx].resize(specsR[specsIdx].size());
				for(peakIdx=0; peakIdx<peakUsed[specsIdx].size(); peakIdx++) peakUsed[specsIdx][peakIdx]=false;
			}
			// Find peaks that are matched at least once
			for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) {
				if(leftAligned[specsIdx])  // Set specsToMerge and endpoints assuming that all matching peaks are b-ions
	//				{ specsToMerge[specsIdx] = specsL[specsIdx];  endpoints[specsIdx][0][0]=specsL[specsIdx].parentMass-AAJumps::massMH+ionOffset; bounds[specsIdx].set(bounds[0][0],endpoints[specsIdx][0][0]); }
					{ specsToMerge[specsIdx] = specsL[specsIdx];  endpoints[specsIdx][0][0]=specsL[specsIdx].parentMass-(allSpecsMS2?(AAJumps::massMH):(AAJumps::massHion))+ionOffset; bounds[specsIdx].set(bounds[0][0],endpoints[specsIdx][0][0]); }
				else
					{ specsToMerge[specsIdx] = specsR[specsIdx];  endpoints[specsIdx][0][0]=specsL[0].parentMass-specsL[specsIdx].parentMass+ionOffset; bounds[specsIdx].set(endpoints[specsIdx][0][0],bounds[0][1]); }

				for(specsIdx2=0; specsIdx2<specsIdx; specsIdx2++) {
					vector<int> *matched, *matched2;  // Pointers to the vectors of matched peaks between specsIdx/specsIdx2 in the current configuration
					if(specsIdx2==0) {
						if(leftAligned[specsIdx]) { matched=&matchesSame[specsIdx][0][0]; matched2=&matchesSame[specsIdx][0][1]; }
						else { matched=&matchesSame[0][specsIdx][1]; matched2=&matchesSame[0][specsIdx][0]; }
					} else {
						if(leftAligned[specsIdx] and leftAligned[specsIdx2])
							{ matched=&matchesSame[specsIdx][specsIdx2][0]; matched2=&matchesSame[specsIdx][specsIdx2][1]; }
						if(not leftAligned[specsIdx] and not leftAligned[specsIdx2])
							{ matched=&matchesSame[specsIdx2][specsIdx][1]; matched2=&matchesSame[specsIdx2][specsIdx][0]; }
						if(leftAligned[specsIdx] and not leftAligned[specsIdx2])
							{ matched=&matchesDiff[specsIdx][specsIdx2][0]; matched2=&matchesDiff[specsIdx][specsIdx2][1]; }
						if(not leftAligned[specsIdx] and leftAligned[specsIdx2])
							{ matched=&matchesDiff[specsIdx2][specsIdx][1]; matched2=&matchesDiff[specsIdx2][specsIdx][0]; }
					}
					for(peakIdx=0; peakIdx<matched->size(); peakIdx++) peakUsed[specsIdx][(*matched)[peakIdx]]=true;
					for(peakIdx=0; peakIdx<matched2->size(); peakIdx++) peakUsed[specsIdx2][(*matched2)[peakIdx]]=true;
				}
			}

	//for(unsigned int i=0; i<specsToMerge.size(); i++)
	//	{ cerr<<"specsToMerge["<<i<<"]:\n"; specsToMerge[i].output(cerr); }

			// Merge all matched peaks from all spectra
	//		MergeSpecs(specsToMerge, peakTol, resolution, mpPenalty, specIP);
			if(binomialScores) MergeSpecs(specsToMerge, peakTol, resolution, mpPenalty, specIP, &bounds, binomialScores, &peakUsed);
			else MergeSpecs(specsToMerge, peakTol, resolution, mpPenalty, specIP, 0, 0, &peakUsed);
	//		MergeSpecs(specsToMerge, peakTol, resolution, mpPenalty, specIP, &bounds, binomialScores, &peakUsed);
	//		MergeSpecs(specsToMerge, peakTol, resolution, mpPenalty, specIP);
	//cerr<<"specIP:\n";
	//specIP.output(cerr);

	//cerr<<"Configuration: "; for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) cerr<<leftAligned[specsIdx]<<","; cerr<<" spec.size() = "<<specConf.size()<<":\n";
			// Merge endpoints into a single spectrum
	//		MergeSpecs(endpoints, pmTol, resolution, mpPenalty, specEP);
			MergeSpecs(endpoints, 2*peakTol, resolution, mpPenalty, specEP);
	//cerr<<"specEP:\n";
	//specEP.output(cerr);

			// Merge endpoints and internal peaks
			MergeIntoReference(specEP, specIP, peakTol, resolution, specConf);

			// Filter peaks with a negative score (it's always better to jump over anyway)
			specConf.filterLowIntensity(0);

			totCurScore=0; for(peakIdx=0; peakIdx<specConf.size(); peakIdx++) totCurScore+=specConf[peakIdx][1];

			list<Spectrum> seqs;
			vector<int> tmpIdx1, tmpIdx2;
	//cerr<<"specConf:\n";
	//specConf.output(cerr);

	//cerr<<"--- sequencing configuration ("; for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) cerr<<leftAligned[specsIdx]<<","; cerr<<"), num peaks = "<<specConf.size()<<", bestScore = "<<bestScore<<"\n";

	//char tmpFN[255]; sprintf(tmpFN,"specConf_%d.pkl",iterIdx);
	//SpecSet tmpSet(1); tmpSet[0]=specConf; tmpSet.SaveSpecSet_pkl(tmpFN);

			if(enforcePM) {
				if(specsL.size()>1) curScore=denovo_exactPM(specConf, jumps, peakTol, pmTol, symOffset, resolution, 1, seqs, false, bestScore);
				else curScore=denovo_exactPM(specsL[0], jumps, peakTol, pmTol, symOffset, resolution, 1, seqs, false, bestScore);
			} else {
				if(specsL.size()>1) denovo(specConf, jumps, peakTol, pmTol, symOffset, 1, seqs);
				else denovo(specsL[0], jumps, peakTol, pmTol, symOffset, 1, seqs);
				curScore=0;
				if(!seqs.empty()) {
					for(peakIdx=0; peakIdx<seqs.front().size(); peakIdx++) curScore+=seqs.front()[peakIdx][1];
	//				FindMatchPeaksAll(specConf, seqs.front(), 0, peakTol, tmpIdx1, tmpIdx2);
	//				for(peakIdx=0; peakIdx<tmpIdx1.size(); peakIdx++) curScore+=specConf[tmpIdx1[peakIdx]][1];
				}
			}

			if(seqs.size()>0) {
				// Best %-explained score
	//			if(curScore/totCurScore>bestScore) {
	//				bestScore = curScore/totCurScore;   consensus = specDenovo;
				// Highest summed score
				if(curScore>bestScore) {
					bestScore = curScore;   consensus = seqs.front();
	cerr<<"Best score set to "<<bestScore<<" for configuration ("; for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) cerr<<leftAligned[specsIdx]<<","; cerr<<"), "<<seqs.size()<<" possible de novo seqs\n";
					if(allTopSeqs!=0) { allTopSeqs->clear(); allTopSeqs->swap(seqs); }

					// Generate consensus of normalized spectra used to estimate distribution of peptide scores for current configuration
					// Merge normalized spectra
	/*				for(specsIdx=1; specsIdx<specsL.size(); specsIdx++)
						if(leftAligned[specsIdx])  // Set specsToMerge and endpoints assuming that all matching peaks are b-ions
							{ specsToMergeNorm[specsIdx] = specsLnorm[specsIdx]; }
						else
							{ specsToMergeNorm[specsIdx] = specsRnorm[specsIdx]; }
					MergeSpecs(specsToMergeNorm, resolution, resolution, specIP, &bounds, binomialScores, &peakUsed);
					MergeSpecs(specsToMerge, resolution, resolution, specsL[0], &bounds, binomialScores, &peakUsed);
	*/
	//				if(binomialScores) MergeSpecs(specsToMerge, peakTol, resolution, mpPenalty, specsL[0], &bounds, binomialScores, &peakUsed);
	//				else MergeSpecs(specsToMerge, peakTol, resolution, mpPenalty, specsL[0], 0, 0, &peakUsed);
	//				specsL[0].parentMass = specsLnorm[0].parentMass;   // Output consensus spectrum for DEBUG purposes
					if(specsL.size()>1) specsL[0] = specConf;
				}
			}
		}

	//	specsL[0] = specsLnorm[0];  // For DEBUG purposes only
	}

	/* MergeSpecs - Merges all the peaks from all spectra in specs. A peak with mass m is
	 *   created in consensus if it exists in at least one spectrum in specs; its intensity
	 *   is set to the weighted sum of all peak masses within tolerance. Weighted sums are
	 *   computed by multiplying the intensity to merge by a linearly decreasing factor of
	 *   the peak mass offset (i.e. a peak at m+mergeTol only contributes 10% of its intensity
	 *   to the final score of peak m).
	 *
	 *   mpPenalty     - Score increment if a peak is not found in a spectrum (e.g. -5)
	 *   peaksAvailable - If not NULL then peak i from spectrum j is only considered to
	 *                      exist if peakAvailable[j][i]==true. If this parameter is NULL (default)
	 *                      then all peaks in specs are used.
	 *   bounds         - pos[i] contains the min,max peak masses for spectrum i
	 *   binomialScores - pos[i][j] contains the binomial log-score(signal/noise) for observing j peaks out of n+1 spectra
	*/
	void MergeSpecs(SpecSet &specs, float mergeTol, float resolution, float mpPenalty,
					 Spectrum &consensus, vector<TwoValues<float> > *bounds, vector<vector<float> > *binomialScores,
									 vector<vector<bool> > *peakAvailable) {
		vector<unsigned int> indices(specs.size());     // Used to iterate over specs, one mass value at a time
		vector<unsigned int> peakMasses(specs.size());  // Used to keep track of peaks' masses according to the current resolution (to ignore under-resolution peak mass differences)
		Spectrum ms3scores;        // Used to keep track of the score contributions from the MS3 spectra
		float curMass;             // Current mass being added to consensus
		unsigned int curMassInt;   // Integer version of curMass, used to estimate whether peaks have the same mass (within resolution)
		unsigned int specsIdx,     // Iterates specs
								 pivotIdx,     // Iterates peaks in current spectrum while adding their contribution to the consensus
									 consensusIdx; // Iterates consensus
		bool iterate=false;        // Marks whether there are still peaks left to process
		TwoValues<unsigned int> peakCounts;  // For each consensus spectrum mass, counts the number of times that a peak (was found, could be found). Used for the binomial factor.
		float penaltyFactor;

		const float symmetryPM = consensus.parentMass-AAJumps::massHion*consensus.idDist;

		/* Inner workings
		 *  - Creates a peak in consensus for every peak mass in specs, equal peak masses only create one consensus peak
		 *  - Peak mass offsets (to the consensus peak) induce a linear multiplicative
		 *      penalty on each peak's score contribution, decreasing to zero as the offset
		 *      approaches mergeTol
		 *  - If peakAvailable is not NULL then it is used to check whether peak P from spectrum S
		 *      can be used by checking whether (*peakAvailable)[S][P]==true. Otherwise all
		 *      peaks are considered available.
		 *  - All fields in consensus are supposed to have been pre-initialized by a previous call to
		 *      Spectrum::copyNP(). This function only changes peak masses/intensities.
		*/

		unsigned int maxNumPeaks=0;   curMass=1000000; // Supremum of any MS/MS peak mass
		for(specsIdx=0; specsIdx<specs.size(); specsIdx++) {
			indices[specsIdx]=0;   maxNumPeaks+=specs[specsIdx].size();
			if(peakAvailable)
				for(;indices[specsIdx]<specs[specsIdx].size() and not (*peakAvailable)[specsIdx][indices[specsIdx]]; indices[specsIdx]++);
			if(indices[specsIdx]<specs[specsIdx].size()) {
				curMass=min(curMass,specs[specsIdx][indices[specsIdx]][0]);
				peakMasses[specsIdx]=(unsigned int)round(specs[specsIdx][indices[specsIdx]][0]/resolution);
				iterate=true;
			}
		}
		consensus.resize(maxNumPeaks);   ms3scores.resize(maxNumPeaks);

		consensusIdx=0;
		while(iterate) {
			curMassInt = (unsigned int)round(curMass/resolution);
			consensus[consensusIdx].set(((float)curMassInt)*resolution,0);  // Guarantees that used mass is rounded towards current resolution
			ms3scores[consensusIdx].set(((float)curMassInt)*resolution,0);  // ms3scores is set to have the exact same masses as consensus
			peakCounts.set(0,0);

			// Add intensities of all peaks within tolerance and estimate next curMass
			iterate=false;   curMass=1000000; // Supremum of any MS/MS peak mass
			for(specsIdx=0; specsIdx<specs.size(); specsIdx++) {
				if(specs[specsIdx].size()==0) continue; // Ignore empty spectra

				// Find peaks within tolerance of curMass and add their intensity to the current consensus peak
				for(pivotIdx = indices[specsIdx]>0?indices[specsIdx]-1:0; pivotIdx>0 and specs[specsIdx][pivotIdx][0]>=consensus[consensusIdx][0]-mergeTol; pivotIdx--);
				if(specs[specsIdx][pivotIdx][0]<consensus[consensusIdx][0]-mergeTol) pivotIdx++;
				bool peakFound=false;  // Set to true below if a peak exists in the current spectrum for the current mass
				float peakScore=0;
				for(; pivotIdx<specs[specsIdx].size() and specs[specsIdx][pivotIdx][0]<=consensus[consensusIdx][0]+mergeTol; pivotIdx++) {
					if(peakAvailable and not (*peakAvailable)[specsIdx][pivotIdx]) continue;
					penaltyFactor = 1 - 0.9*(fabs(consensus[consensusIdx][0]-specs[specsIdx][pivotIdx][0])/mergeTol);
					// Supporting variable to add only maximum scoring peak
					peakScore = max(peakScore,penaltyFactor * specs[specsIdx][pivotIdx][1]);
					// Adds the intensity of all peaks in range
	//				peakScore += penaltyFactor * specs[specsIdx][pivotIdx][1];
					if(specsIdx>0) ms3scores[consensusIdx][1] += penaltyFactor * specs[specsIdx][pivotIdx][1];  // NOTE: 2* factor assumes that all ms3 spectra are symmetric!
					peakFound=true;
				}
				// Adds only maximum scoring peak
				consensus[consensusIdx][1] += peakFound ? peakScore : mpPenalty;
	//			if(specsIdx>0) ms3scores[consensusIdx][1] += peakScore;

				if(bounds and consensus[consensusIdx][0]>(*bounds)[specsIdx][0]+mergeTol and consensus[consensusIdx][0]<(*bounds)[specsIdx][1]-mergeTol)
					{ peakCounts[1]++;   peakCounts[0] += peakFound?1:0; }

				// Advance indices for the spectra that set curMassInt
				if(indices[specsIdx]<specs[specsIdx].size() and peakMasses[specsIdx]>curMassInt)
					{ curMass = min(curMass,specs[specsIdx][indices[specsIdx]][0]); iterate=true; }
				else while(indices[specsIdx]<specs[specsIdx].size() and peakMasses[specsIdx]==curMassInt) {
					indices[specsIdx]++;
					if(peakAvailable)
						for(;indices[specsIdx]<specs[specsIdx].size() and not (*peakAvailable)[specsIdx][indices[specsIdx]]; indices[specsIdx]++);
					if(indices[specsIdx]<specs[specsIdx].size()) {
						curMass=min(curMass,specs[specsIdx][indices[specsIdx]][0]);
						peakMasses[specsIdx]=(unsigned int)round(specs[specsIdx][indices[specsIdx]][0]/resolution);
						iterate = true;
					}
				}
			}
			// Binomial premium/penalty for the number of observed/expected peaks
			float premium=0;
			if(bounds and binomialScores) {
	//			premium = log( (pow(pTrue,peakCounts[0])*pow(1-pTrue,peakCounts[1]-peakCounts[0])) / (pow(pNoise,peakCounts[0])*pow(1-pNoise,peakCounts[1]-peakCounts[0])) );
				if(peakCounts[1]>=1 and peakCounts[1]<=binomialScores->size() and peakCounts[0]<(*binomialScores)[peakCounts[1]-1].size())
					premium = (*binomialScores)[peakCounts[1]-1][peakCounts[0]];
				consensus[consensusIdx][1] += premium;   ms3scores[consensusIdx][1] += premium;
	//cerr << " -- counts ("<<peakCounts[0]<<","<<peakCounts[1]<<"), premium = "<<premium<<", consensus["<<consensusIdx<<"] = ("<<consensus[consensusIdx][0]<<","<<consensus[consensusIdx][1]<<")\n";
			}
			consensusIdx++;
		}
		consensus.resize(consensusIdx);    ms3scores.resize(consensusIdx);

		float maxScore;  // Symmetric y-ion for the current b-ion
		vector<int> matches;
		unsigned int matchesIdx;
		if(specs.size()>1)
			for(consensusIdx=0; consensusIdx<consensus.size(); consensusIdx++) {
				ms3scores.findMatches(symmetryPM-consensus[consensusIdx][0], mergeTol, matches);
	// TODO: mpPenalty here should be proportional to the number of spectra without peaks at this mass
				maxScore=mpPenalty; for(matchesIdx=0; matchesIdx<matches.size(); matchesIdx++) if(ms3scores[matches[matchesIdx]][1]>maxScore) maxScore=ms3scores[matches[matchesIdx]][1];
	//			consensus[consensusIdx][1] -= maxScore;  // Subtract evidence that complementary mass is b-ion
	//			consensus[consensusIdx][1] += ms3scores[consensusIdx][1]-maxScore;  // Subtract evidence that complementary mass is b-ion
			}

	}

	// MergeIntoReference - merges a set of of new peaks (newPeaks) into a reference spectrum (refSpec)
	//  For any given unique mass m (within resolution) existing in either/both refSpec/newPeaks:
	//   if it exists in refSpec then its intensity is incremented using all newPeaks within tolerance (after linear mass offset penalties are imposed)
	//   otherwise its intensity is exactly the same as in newPeaks.
	void MergeIntoReference(Spectrum &refSpec, Spectrum &newPeaks, float mergeTol, float resolution, Spectrum &result) {
		result.copyNP(refSpec);   result.resize(refSpec.size()+newPeaks.size());
		unsigned int refIdx=0, newIdx=0, resultIdx=0,
								 pivotIdx; // used to find all newPeaks within tolerance of refSpec[refIdx]
		float penaltyFactor;

		// On each iteration, the index to the smallest mass is incremented (either refIdx or newIdx)
		while(refIdx<refSpec.size() or newIdx<newPeaks.size()) {
			if(refIdx==refSpec.size() or (newIdx<newPeaks.size() and (int)round(newPeaks[newIdx][0]/resolution)<(int)round(refSpec[refIdx][0]/resolution)))
				{ result[resultIdx++] = newPeaks[newIdx++]; continue; }

			// Find first new peak within tolerance of reference peak
			for(pivotIdx = (newIdx>0?newIdx-1:0); pivotIdx>0 and newPeaks[pivotIdx][0]>=refSpec[refIdx][0]-mergeTol; pivotIdx--);
			if(newPeaks[pivotIdx][0]<refSpec[refIdx][0]-mergeTol) pivotIdx++;
			for(; pivotIdx<newPeaks.size() and newPeaks[pivotIdx][0]<=refSpec[refIdx][0]+mergeTol; pivotIdx++) {
				penaltyFactor = 1 - 0.9*(fabs(refSpec[refIdx][0]-newPeaks[pivotIdx][0])/mergeTol);
				refSpec[refIdx][1]+=newPeaks[pivotIdx][1]*penaltyFactor;
			}

			int roundedRefMass = (int)round(refSpec[refIdx][0]/resolution);
			while(newIdx<newPeaks.size() and (int)round(newPeaks[newIdx][0]/resolution)==roundedRefMass)
				newIdx++;
			result[resultIdx++] = refSpec[refIdx++];
		}
		result.resize(resultIdx);
	}
}
