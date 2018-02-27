#ifndef MSN_H
#define MSN_H

#include "batch.h"
#include "spectrum_scoring.h"

namespace specnets
{
	using namespace std;

	/**
	 * TODO: add description
	 *
	 *@param specs
	 *@param aligns
	 *@param sets
	 *@param maxNumMS3
	 */
	void get_ms3_sets_max(SpecSet &specs, vector<Results_ASP> &aligns, vector<list<
			int> > &sets, int maxNumMS3 = -1);

	/**
	 * TODO: add description
	 *
	 *@param specs
	 *@param aligns
	 *@param sets
	 *@param pmTol
	 *@param maxNumMS3
	 */
	void get_ms3_sets_sparse(SpecSet &specs, vector<Results_ASP> &aligns, vector<
			list<int> > &sets, float pmTol, unsigned int maxNumMS3 = 0);

	/**
	 * TODO: add description
	 *
	 *@param specs
	 *@param mergeTol
	 *@param resolution
	 *@param mpPenalty
	 *@param consensus
	 *@param bounds
	 *@param binomialScores
	 *@param peakAvailable
	 */
	void MergeSpecs(SpecSet &specs, float mergeTol, float resolution,
			float mpPenalty, Spectrum &consensus,
			vector<TwoValues<float> > *bounds = 0,
			vector<vector<float> > *binomialScores = 0,
			vector<vector<bool> > *peakAvailable = 0);

	/**
	 * TODO: add description
	 *
	 *@param allSpecsMS2
	 *@param ionOffset
	 *@param peakTol
	 *@param pmTol
	 *@param resolution
	 *@param mpPenalty
	 *@param ms2model
	 *@param ms3model
	 *@param ms3modelY
	 *@param jumps
	 *@param consensus
	 *@param allTopSeqs
	 *@param enforcePM
	 *@param binomialScores
	 */
	void denovo_ms3(SpecSet &specs, bool allSpecsMS2, float ionOffset,
			float peakTol, float pmTol, float resolution, float mpPenalty,
			MS2ScoringModel &ms2model, MS3ScoringModel &ms3model,
			MS3ScoringModel &ms3modelY, AAJumps &jumps, Spectrum &consensus, list<
					Spectrum> *allTopSeqs = (list<Spectrum> *) 0, bool enforcePM =
					true, vector<vector<float> > *binomialScores = 0,
			float minPerc = 0);

	/**
	 * TODO: add description
	 *
	 *@param refSpec
	 *@param newPeaks
	 *@param mergeTol
	 *@param resolution
	 *@param result
	 */
	void MergeIntoReference(Spectrum &refSpec, Spectrum &newPeaks, float mergeTol,
			float resolution, Spectrum &result);
}

#endif
