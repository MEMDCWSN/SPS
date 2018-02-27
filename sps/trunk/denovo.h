#ifndef DENOVO_H
#define DENOVO_H

#include "spectrum.h"
#include "aminoacid.h"

namespace specnets
{
	using namespace std;

		//void denovo(Spectrum &spec, AAJumps &jumps, float tolerance, float pmOffset, float minScorePerc,
	//                list<Spectrum> &seqs, bool symmetric=false);


	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param jumps
	 *@param peakTol
	 *@param neighs
	 */
	void getNeighborsL(Spectrum &spec, AAJumps &jumps, float peakTol, vector<
			vector<int> > &neighs);

	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param minPeakDist
	 *@param peakTol
	 *@param neighs
	 */
	void getNeighborsL(Spectrum &spec, float minPeakDist, float peakTol, vector<
			vector<int> > &neighs);

	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param minPeakDist
	 *@param maxPeakDist
	 *@param peakTol
	 *@param neighs
	 */
	void getNeighborsL(Spectrum &spec, float minPeakDist, float maxPeakDist,
			float peakTol, vector<vector<int> > &neighs);

	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param jumps
	 *@param peakTol
	 *@param neighs
	 */
	void getNeighborsR(Spectrum &spec, AAJumps &jumps, float peakTol, vector<
			vector<int> > &neighs);

	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param minPeakDist
	 *@param maxPeakDist
	 *@param peakTol
	 *@param neighs
	 */
	void getNeighborsR(Spectrum &spec, float minPeakDist, float maxPeakDist,
			float peakTol, vector<vector<int> > &neighs);

	/**
	 * TODO: add description
	 *
	 *@param inSpec
	 *@param jumps
	 *@param peakTol
	 *@param pmTol
	 *@param pmOffset
	 *@param minScorePerc
	 *@param seqs
	 *@param enforcePM
	 */
	void denovo(Spectrum &inSpec, AAJumps &jumps, float peakTol, float pmTol,
			float pmOffset, float minScorePerc, list<Spectrum> &seqs,
			bool enforcePM = false);

	/**
	 * TODO: add description
	 *
	 *@param inSpec
	 *@param neighsL
	 *@param neighsR
	 *@param peakTol
	 *@param pmTol
	 *@param ionOffset
	 *@param minScorePerc
	 *@param seqs
	 *@param seqScores
	 *@param ctermH2O
	 *@param enforcePM
	 */
	void denovo(Spectrum &inSpec, vector<vector<int> > &neighsL,
			vector<vector<int> > &neighsR, float peakTol, float pmTol,
			float ionOffset, float minScorePerc, list<Spectrum> &seqs,
			list<float> &seqScores, bool ctermH2O = true, bool enforcePM = false);

	/**
	 * TODO: add description
	 *
	 *@param inSpec
	 *@param jumps
	 *@param peakTol
	 *@param pmTol
	 *@param pmOffset
	 *@param resolution
	 *@param minScorePerc
	 *@param seqs
	 *@param seqsExactMasses
	 *@param minAbsScore
	 *@return
	 */
	float denovo_exactPM(Spectrum &inSpec, AAJumps &jumps, float peakTol,
			float pmTol, float pmOffset, float resolution, float minScorePerc,
			list<Spectrum> &seqs, bool seqsExactMasses = false, float minAbsScore =
					0);

	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param neighs
	 *@param idxMatched
	 *@return
	 */
	float denovo_LtoR(Spectrum &spec, vector<vector<int> > &neighs,
			vector<int> &idxMatched);
	//float denovo_LtoR(Spectrum &spec, vector<vector<int> > &neighs, list<list<int> > &seqs, float minPercMaxScore);

	/**
	 * TODO: add description
	 *
	 *@param spec
	 *@param neighs
	 *@param seqs
	 *@param seqScores
	 *@param minPercMaxScore
	 *@return
	 */
	float denovo_LtoR(Spectrum &spec, vector<vector<int> > &neighs,
			list<Spectrum> &seqs, list<float> &seqScores, float minPercMaxScore);

	/**
	 * Code to estimate the distribution of peptide scores over a given spectrum.
	 *
	 *@param spec
	 *@param peakTol
	 *@param pmTol
	 *@param resolution
	 *@param hist
	 *@param addIntensities
	 */
	void EstimateScoresDistribution(Spectrum &spec, float peakTol, float pmTol,
			float resolution, vector<unsigned int> &hist, bool addIntensities =
					true);

	/**
	 *
	 * TODO: add description
	 *
	 *@param spec
	 *@param bins
	 *@param peakTol
	 *@param resolution
	 */

	void peaksToMassBins(Spectrum &spec, Spectrum &bins, float peakTol,
			float resolution);
}
#endif
