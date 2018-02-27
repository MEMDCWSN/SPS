#include "msn.h"
#include "alignment_scoring.h"
#include "denovo.h"
#include <vector>

void MergeSpecs(SpecSet &specs, float mergeTol, float resolution, Spectrum &consensus, 
                 vector<vector<bool> > *peakAvailable=0);
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
void denovo_ms3(SpecSet &specsL, float ionOffset, float peakTol, float pmTol, float resolution, 
                  AAJumps &jumps, Spectrum &consensus) {
	SpecSet specsR(specsL.size());  // Contains the right-aligned variants of every spectrum
	vector<vector<TwoValues<vector<short> > > > matchesSame;  
	                   // pos[i][j], i>j: matched peaks when spectra i,j are both aligned on the left  (lower triangular)
	                   // pos[i][j], i<j: matched peaks when spectra i,j are both aligned on the right (upper triangular)
                       // on both variants, [i][j][0] - list of matched peaks in spectrum i
                       //                   [i][j][1] - list of matched peaks in spectrum j
                       // i==j always empty
	vector<vector<TwoValues<vector<short> > > > matchesDiff;
	                   // pos[i][j], i>j: matched peaks when spectra i,j are aligned on the right/left  (lower triangular)
	                   // pos[i][j], i<j: matched peaks when spectra i,j are aligned on the left/right (upper triangular)
                       // on both variants, [i][j][0] - list of matched peaks in spectrum i
                       //                   [i][j][1] - list of matched peaks in spectrum j
                       // i==j, i==0, j==0 always empty
	vector<vector<bool> > peakUsed;       // pos[i][j]: indicates whether peak j from spectrum i was matched at least once in the current configuration
	SpecSet endpoints(specsL.size());     // Contains the endpoint peaks for each configuration
	SpecSet specsToMerge(specsL.size());  // Contains the spectra to merge for each configuration
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
	      bestScore=0; // Score of the denovo interpretation of the best configuration
	
//	MergeIntoReference(specsL[0], specsL[1], peakTol, resolution, consensus);
//	MergeSpecs(specsL, peakTol, resolution, consensus);

	// Initialize variables
	matchesSame.resize(specsL.size());   matchesDiff.resize(specsL.size());   peakUsed.resize(specsL.size());
	for(specsIdx=0; specsIdx<specsL.size(); specsIdx++) { 
		matchesSame[specsIdx].resize(specsL[specsIdx].size());
		matchesDiff[specsIdx].resize(specsL[specsIdx].size());
		peakUsed[specsIdx].resize(specsL[specsIdx].size());
		if(specsIdx==0) {
			endpoints[specsIdx].resize(2);  endpoints[0][0].set(0,0);
			for(peakIdx=0; peakIdx<specsL[0].size(); peakIdx++) endpoints[0][0][1]=max(endpoints[0][0][1],specsL[0][peakIdx][1]);
			endpoints[0][1].set(specsL[0].parentMass-AAJumps::massMH+ionOffset,endpoints[0][0][1]);
			specsR[specsIdx]=specsL[specsIdx];   // specsL[0] is a special case: base MS/MS spectrum doesn't need to be right-shifted
		} else {
			endpoints[specsIdx].resize(1);   endpoints[specsIdx][0][1]=0;
			specsR[specsIdx]=specsL[specsIdx];   
			mass = max((float)0,specsL[0].parentMass-specsL[specsIdx].parentMass);
			for(peakIdx=0; peakIdx<specsR[specsIdx].size(); peakIdx++) {
				specsR[specsIdx][peakIdx][0]+=mass;  // Shift all peak masses in specsR
				endpoints[specsIdx][0][1]=max(endpoints[specsIdx][0][1],specsR[specsIdx][peakIdx][1]);  // Set endpoint score to maximum score
//				endpoints[specsIdx][0][1]+=specsR[specsIdx][peakIdx][1];  // Set endpoint score to total summed score
			}
		}
	}
	
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
	specConf.copyNP(specsL[0]);   // Initialize all non-peak values in specConf
	for(peakIdx=0; peakIdx<specsL[0].size(); peakIdx++) peakUsed[0][peakIdx]=true; // Always use all peaks from the main spectrum
	for(iterIdx=0; iterIdx<iterCount; iterIdx++) {
		// Initialize leftAligned to current configuration
		for(specsIdx=1; specsIdx<specsL.size(); specsIdx++)
			leftAligned[specsIdx]=iterIdx&(1<<(specsIdx-1));

		// Reset peakUsed
//		for(specsIdx=0; specsIdx<specsL.size(); specsIdx++)   // Only use peaks from the main spectrum if they're matched to other spectra
		for(specsIdx=1; specsIdx<specsL.size(); specsIdx++)   // Always use all peaks from the main spectrum
			for(peakIdx=0; peakIdx<specsL[specsIdx].size(); peakIdx++) peakUsed[specsIdx][peakIdx]=false;

		// Find peaks that are matched at least once
		for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) {
			if(leftAligned[specsIdx])  // Set specsToMerge and endpoints assuming that all matching peaks are b-ions
				{ specsToMerge[specsIdx] = specsL[specsIdx];  endpoints[specsIdx][0][0]=specsL[specsIdx].parentMass-AAJumps::massMH+ionOffset; }
			else 
				{ specsToMerge[specsIdx] = specsR[specsIdx];  endpoints[specsIdx][0][0]=specsL[0].parentMass-specsL[specsIdx].parentMass+ionOffset; }

			for(specsIdx2=0; specsIdx2<specsIdx; specsIdx2++) {
				vector<short> *matched, *matched2;  // Pointers to the vectors of matched peaks between specsIdx/specsIdx2 in the current configuration
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
		
		// Merge all matched peaks from all spectra
		MergeSpecs(specsToMerge, peakTol, resolution, specIP, &peakUsed);
//		MergeSpecs(specsToMerge, peakTol, resolution, specIP);
		
//cerr<<"Configuration: "; for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) cerr<<leftAligned[specsIdx]<<","; cerr<<" spec.size() = "<<specConf.size()<<":\n";
		// Merge endpoints into a single spectrum
		MergeSpecs(endpoints, pmTol, resolution, specEP);
//cerr<<"specEP:\n";
//specEP.output(cerr);
		
		// Merge endpoints and internal peaks
		MergeIntoReference(specEP, specIP, peakTol, resolution, specConf);

		totCurScore=0; for(peakIdx=0; peakIdx<specConf.size(); peakIdx++) totCurScore+=specConf[peakIdx][1];
		
		list<Spectrum> seqs;
		vector<short> tmpIdx1, tmpIdx2;
//cerr<<"specConf:\n";
//specConf.output(cerr); // <-- Use this to determine why endpoint scores look incorrect (only one was maximized?)
		denovo(specConf, jumps, peakTol, pmTol, 2*ionOffset, 1, seqs);
		if(seqs.size()>0) {
			specDenovo = seqs.front();
//cerr<<"specDenovo:\n";
//specDenovo.output(cerr);			
			FindMatchPeaksAll(specConf, specDenovo, 0, peakTol, tmpIdx1, tmpIdx2);
			curScore=0; for(peakIdx=0; peakIdx<tmpIdx1.size(); peakIdx++) curScore+=specConf[tmpIdx1[peakIdx]][1];
			if(curScore/totCurScore>bestScore) {
				bestScore = curScore/totCurScore;   consensus = specDenovo;
//cerr<<"Best score set to "<<bestScore<<" for configuration ("; for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) cerr<<leftAligned[specsIdx]<<","; cerr<<"):\n";
			}
		}
		
		// Change specConf assuming that all peaks are y-ions - what's the point? The best result is when all are b-ions!
		//   - reversing code missing
/*		denovo(specConf, jumps, peakTol, 2*ionOffset, 1, tmp, false);
		specDenovo = tmp.front();
		FindMatchPeaksAll(specConf, specDenovo, 0, peakTol, tmpIdx1, tmpIdx2);
		curScore=0; for(peakIdx=0; peakIdx<tmpIdx1.size(); peakIdx++) curScore+=specConf[tmpIdx1[peakIdx]][1];
		if(curScore/totCurScore>bestScore) {
			bestScore = curScore/totCurScore;   consensus = specDenovo;
cerr<<"Best score set to "<<bestScore<<" for configuration ("; for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) cerr<<leftAligned[specsIdx]<<","; cerr<<"):\n";
		}
*/

//cerr<<"Spectrum for configuration ("; for(specsIdx=1; specsIdx<specsL.size(); specsIdx++) cerr<<leftAligned[specsIdx]<<","; cerr<<"):\n";
//specConf.output(cerr);
	}
}

/* MergeSpecs - Merges all the peaks from all spectra in specs. A peak with mass m is
 *   created in consensus if it exists in at least one spectrum in specs; its intensity
 *   is set to the weighted sum of all peak masses within tolerance. Weighted sums are
 *   computed by multiplying the intensity to merge by a linearly decreasing factor of
 *   the peak mass offset (i.e. a peak at m+mergeTol only contributes 10% of its intensity
 *   to the final score of peak m).
 * 
 *   peaksAvailable - If not NULL then peak i from spectrum j is only considered to
 *                      exist if peakAvailable[j][i]==true. If this parameter is NULL (default)
 *                      then all peaks in specs are used.
*/
void MergeSpecs(SpecSet &specs, float mergeTol, float resolution, Spectrum &consensus, 
                 vector<vector<bool> > *peakAvailable) {
	vector<unsigned int> indices(specs.size());     // Used to iterate over specs, one mass value at a time
	vector<unsigned int> peakMasses(specs.size());  // Used to keep track of peaks' masses according to the current resolution (to ignore under-resolution peak mass differences)
	float curMass;             // Current mass being added to consensus
	unsigned int curMassInt;   // Integer version of curMass, used to estimate whether peaks have the same mass (within resolution)
	unsigned int specsIdx,     // Iterates specs
	             pivotIdx,     // Iterates peaks in current spectrum while adding their contribution to the consensus
                 consensusIdx; // Iterates consensus
	bool iterate=false;        // Marks whether there are still peaks left to process
	float penaltyFactor;
	
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
	consensus.resize(maxNumPeaks);
	
	consensusIdx=0;
	while(iterate) {
		curMassInt = (unsigned int)round(curMass/resolution);
		consensus[consensusIdx].set(((float)curMassInt)*resolution,0);  // Guarantees that used mass is rounded towards current resolution

		// Add intensities of all peaks within tolerance and estimate next curMass
		iterate=false;   curMass=1000000; // Supremum of any MS/MS peak mass
		for(specsIdx=0; specsIdx<specs.size(); specsIdx++) {
			if(specs[specsIdx].size()==0) continue; // Ignore empty spectra
			
			// Find peaks within tolerance of curMass and add their intensity to the current consensus peak
			for(pivotIdx = indices[specsIdx]>0?indices[specsIdx]-1:0; pivotIdx>0 and specs[specsIdx][pivotIdx][0]>=consensus[consensusIdx][0]-mergeTol; pivotIdx--);
			if(specs[specsIdx][pivotIdx][0]<consensus[consensusIdx][0]-mergeTol) pivotIdx++;
			for(; pivotIdx<specs[specsIdx].size() and specs[specsIdx][pivotIdx][0]<=consensus[consensusIdx][0]+mergeTol; pivotIdx++) {
				if(peakAvailable and not (*peakAvailable)[specsIdx][pivotIdx]) continue;
				penaltyFactor = 1 - 0.9*(fabs(consensus[consensusIdx][0]-specs[specsIdx][pivotIdx][0])/mergeTol);
				consensus[consensusIdx][1] += penaltyFactor * specs[specsIdx][pivotIdx][1];
			}

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
		consensusIdx++;
	}

	consensus.resize(consensusIdx);
}

// MergeIntoReference - merges a set of of new peaks (newPeaks) into a reference spectrum (refSpec)
//  For any given unique mass m (within resolution) existing in either/both refSpec/newPeaks:
//   if it exists in refSpec then its intensity is incremented using all newPeaks within tolerance (after linear mass offset penalties are imposed)
//   otherwise its intensity is exactly the same as in newPeaks.
void MergeIntoReference(Spectrum &refSpec, Spectrum &newPeaks, float mergeTol, float resolution, Spectrum &result) {
	result.resize(refSpec.size()+newPeaks.size());
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
