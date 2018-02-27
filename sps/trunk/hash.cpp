#include "hash.h"
#include <cmath>

/**
 * Maximum hashtable peak spectrum mass.
 */
static const float maxPeakMass = 4000;

inline void HashTable::SpectrumToIndices(Spectrum &spec, vector<unsigned int> &indices) {
	unsigned int specIdx, indicesIdx=0, intPeakMass, minBin;
	indices.resize(spec.size()*binsPerPeak);
	for(specIdx=0; specIdx<spec.size(); specIdx++) {
		intPeakMass = (unsigned int)round(spec[specIdx][0]/resolution);
		if(intPeakMass>intPeakTol) minBin=intPeakMass-intPeakTol; else minBin=0;
		for(unsigned int pivot=minBin; pivot<=intPeakMass+intPeakTol and pivot<table.size(); pivot++)
			indices[indicesIdx++]=pivot;
	}
	indices.resize(indicesIdx);
}

//
//  Build - Builds the hashtable index from a set of spectra (every peak from
//   every spectrum generates 1+2*peakTol/resolution entries in the table).
//  endIdx - 1+index of the last spectrum to enter in the table. Set to 0 to enter all spectra.
//
void HashTable::Build(SpecSet &specs, float in_resolution, float peakTol, unsigned int endIdx) {
	vector<unsigned int> counts;   // Counts of how many spectra will collide in each bin. Used for predimensioning the hash table.
	vector<unsigned int> indices;  // Bin indices for a given spectrum
	unsigned int specIdx, pivot, pairIdx;

	// Initialize data fields
	resolution = in_resolution;
	intPeakTol = (unsigned int)round(peakTol/resolution);
	binsPerPeak= 2*intPeakTol+1;
	table.resize((unsigned int)ceil(maxPeakMass/resolution)+intPeakTol);
	indices.reserve(20*binsPerPeak);   // Predimension to 20 peaks per spectrum

	// Count number of entries per bin and predimension table
	if(endIdx==0 or endIdx>specs.size()) endIdx=specs.size();
	counts.resize(table.size());   for(pivot=0; pivot<counts.size(); pivot++) counts[pivot]=0;
	for(specIdx=0; specIdx<endIdx; specIdx++) {
		SpectrumToIndices(specs[specIdx],indices);
		for(pivot=0; pivot<indices.size(); pivot++) counts[indices[pivot]]++;
	}
	for(pivot=0; pivot<counts.size(); pivot++) { table[pivot].resize(counts[pivot]); counts[pivot]=0; }

	// Enter spectrum indices on table
	for(specIdx=0; specIdx<endIdx; specIdx++) {
		SpectrumToIndices(specs[specIdx],indices);
		for(pivot=0; pivot<indices.size(); pivot++) {
			table[indices[pivot]][counts[indices[pivot]]] = specIdx;
			counts[indices[pivot]]++;
		}
	}
}

//
//  QuerySet - Queries a set of spectra on the hashtable and returns the indices of
//    all pairs of spectra with at least one peak hashing to a common bin.
//
//  startIdx - index of the first spectrum to query
//  returnAllPairs - if true, the pairs only outputs pairs (i,j) of spectrum indices
//                     for i<j, i spectrum index in the table, j queried spectrum index
//
unsigned int HashTable::QuerySet(SpecSet &specs, list<vector<unsigned int> > &pairs,
                     unsigned int startIdx, bool returnAllPairs) {
	unsigned int specIdx, peakIdx, pairIdx, bin, maxSpecIdx;
	vector<unsigned int> pair;   pair.resize(2);

	vector<bool> pairUsed(specs.size());  // UGLY and O(n^2)! Replace with a queue of pairs that allows for direct querying of queue elements

	for(specIdx=startIdx; specIdx<specs.size(); specIdx++) {
		for(pairIdx=0; pairIdx<pairUsed.size(); pairIdx++) pairUsed[pairIdx]=false;
		if(returnAllPairs) maxSpecIdx=specs.size(); else maxSpecIdx=specIdx;
		for(peakIdx=0; peakIdx<specs[specIdx].size(); peakIdx++) {
			bin = (unsigned int)round(specs[specIdx][peakIdx][0]/resolution); if(bin>table.size()) break;
			for(pairIdx=0; pairIdx<table[bin].size() and table[bin][pairIdx]<maxSpecIdx; pairIdx++)
				if(!pairUsed[table[bin][pairIdx]]) {
					pair[0] = table[bin][pairIdx];   pair[1] = specIdx;
					pairs.push_back(pair);
					pairUsed[pair[0]]=true;
				}
		}
	}
	return 0;
}
