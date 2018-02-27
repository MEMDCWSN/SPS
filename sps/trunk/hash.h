#ifndef HASH_H
#define HASH_H

#include <vector>
#include <list>

#include "spectrum.h"
#include "batch.h"


/**
 * TODO: add description
 */
class HashTable {

	/**
	 * Hash table.
	 */
	vector<vector<unsigned int> > table;

	/**
	 * Peak tolerance in terms of hash table bins.
	 */
	unsigned int intPeakTol;

	/**
	 * Number of bins per peak. 2*intPeakTol+1
	 */
	unsigned int binsPerPeak;

	/**
	 * Mass resolution used to index the hash table: bin = round(mass/resolution).
	 */
	float resolution;

	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param indices
	 */
	void SpectrumToIndices(Spectrum &spec, vector<unsigned int> &indices);

	//static const float maxPeakMass = 4000;  // 4000 Da maximum hashable spectrum peak mass
public:

	HashTable() {
		table.resize(0);
		binsPerPeak = 0;
		intPeakTol = 0;
	}

	/**
	 * TODO: add description
	 *
	 *@param specs
	 *@param in_resolution
	 *@param peakTol
	 *@param endIdx
	 */
	void Build(SpecSet &specs, float in_resolution, float peakTol,
			unsigned int endIdx);

	/**
	 * TODO: add description
	 *
	 *@param specs
	 *@param pairs
	 *@param startIdx
	 *@param returnAllPairs
	 *@return
	 */
	unsigned int QuerySet(SpecSet &specs, list<vector<unsigned int> > &pairs,
			unsigned int startIdx, bool returnAllPairs);
};

#endif
