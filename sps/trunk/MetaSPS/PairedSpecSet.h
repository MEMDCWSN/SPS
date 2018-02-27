/*
 * PairedSpecSet.h
 *
 *  Created on: Nov 11, 2011
 *      Author: aguthals
 *
 */

#ifndef PAIREDSPECSET_H_
#define PAIREDSPECSET_H_

#include "spectrum.h"
#include "SpecSet.h"

using namespace std;

namespace specnets {

class PairedSpecSet {
public:

	PairedSpecSet();

	PairedSpecSet(SpecSet* allSpectra, bool enforceDaTol = false,
					vector<vector<bool> >* reversedPeaks = 0, float peakTol = 0, float srmOffset = 0.0 - AAJumps::massH2O - AAJumps::massNH);

	~PairedSpecSet(void);

	void initialize(SpecSet* allSpectra, bool enforceDaTol = false, vector<vector<bool> >* reversedPeaks =
			0, float peakTol = 0, float srmOffset = 0.0 - AAJumps::massH2O - AAJumps::massNH);

	void mergePRMs(unsigned short finalStage = 5);

	/**
	 * For every peak in a CID spectrum this will: boost its score if it has a matching ETD PRM or
	 *   move it to a true PRM mass if it has a matching ETD SRM or keep it unchanged if it has no
	 *   matching ETD PRM/SRM. Does the same procedure for ETD peaks compared against CID spectra.
	 * @param boostedSpecs output SpecSet parallel to input SpecSet, except it has boosted PRMs
	 * @return
	 */
	void boostPRMs(SpecSet& boostedSpecs);

	/**
	 * Finds peaks with same mass in CID/ETD pairs or HCD/ETD pairs and translates
	 *   them and their complementary SRMs to the PRM mass in their merged spectrum with summed PRM+SRM
	 *   intensities (at a singular PRM peak mass)
	 */
	void mergePRMsStage1();

	/**
	 * Finds peaks with mass difference 15+18 in CID/ETD pairs or HCD/ETD pairs AND where
	 *   at least one spectrum (CID, HCD or ETD) contains a PRM for these SRMs, translates them and
	 *   their complementary PRMs to the PRM mass in their merged spectrum with summed PRM+SRM
	 *   intensities
	 */
	void mergePRMsStage2();

	/**
	 * Finds PRM / SRM pairs in CID/ETD pairs or HCD/ETD pairs and translates them to the PRM mass
	 *   in their merged spectrum with summed PRM+SRM intensities
	 */
	void mergePRMsStage3();

	/**
	 * Finds peaks with mass difference 15+18 in CID/ETD pairs or HCD/ETD pairs and translates them and
	 *   their complementary PRMs to the PRM mass in their merged spectrum with summed PRM+SRM
	 *   intensities
	 */
	void mergePRMsStage4();

	/**
	 * Copies all un-labeled left-over peaks from CID, HCD, and/or ETD spectra to their merged spectra
	 */
	void mergePRMsStage5();

	void getMergedSpectrum(Spectrum* toSpec);

protected:

	// Merged spectrum.
	Spectrum* mergedSpectrum;

	SpecSet* m_boostedSpectra;

	// Labels for peaks in mergedSet that were derived from multiple complementary peaks in CID, HCD, and ETD.
	// 0 - Unlabeled, 1 - PRM, 2 - SRM, 3 - End point
	vector<short>* mergedLabels;

	// CID, HCD, and ETD PRM spectra
	SpecSet* inputSpectra;
	vector<vector<bool> >* reversedPeaks;
	vector<vector<short> >* usedPeaks;

	bool enforceDaTolerance;

	float m_peakTol;

	float m_srmOffset;

	void
	initializeUsedPeaks(vector<vector<short> >* usedPeaks, SpecSet* spectra);

	void insertMergedLabel(vector<short>* peakLabels, int index, short label);

	void boostPRMsCID(int CIDidx, list<int>& ETDIdxs);

	void boostPRMsETD(int ETDidx, list<int>& CIDDIdxs);

	void
	mergePRMsPair(int CIDIdx, int ETDIdx);

	void mergeSRMsPair(int CIDIdx, int ETDIdx, list<int>* checkPRM = 0);

	void mergePRMsSRMsPair(int CIDIdx, int ETDIdx);

	void mergeLeftovers(int childIdx);

	void addNewPeaks(list<MZRange>* newPeaks, short label);

	void getUniqueCIDETDPairs(list<pair<int, int> >& outputPairs);

	void getUniqueHCDETDPairs(list<pair<int, int> >& outputPairs);

	void getComplementIdxs(int idx1, int idx2, list<int>& outputIdxs);

};
}

#endif /* PAIREDSPECSET_H_ */
