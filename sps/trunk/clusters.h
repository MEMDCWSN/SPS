#ifndef CLUSTERS_H
#define CLUSTERS_H

#include "SpecSet.h"
#include "spectrum.h"
#include "abruijn.h"
#include <vector>

namespace specnets
{
	using namespace std;
	/**
	 * TODO: add description
	 *
	 */
	class Clusters {
	public:

		/**
		 * Indices of the spectra in each cluster (there should be some companion
		 * .pkl spectra file).
		 */
		vector<vector<int> > specIdx;

		/**
		 * Shifts (left and right) of each spectrum in each cluster.
		 */
		vector<vector<TwoValues<float> > > shifts;

		/**
		 * Consensus spectrum information per cluster - Consensus PRMs.
		 */
		SpecSet consensus;

		/**
		 * Endpoint types: 0=internal, 1=b0/bN, 2=y0/yN
		 */
		vector<vector<short> > endpoints;

		/**
		 * Default constructor
		 */
		Clusters();

		/**
		 * Initialize clusters from a set of spectra
		 *
		 *@param spectra Spectra used to initialize the clusters
		 */
		Clusters(SpecSet &spectra);

		/**
		 * Copies another Clusters object
		 */
		Clusters& operator =(const Clusters &other);

		/**
		 * TODO: add description
		 *
		 *@param newSize
		 *@return
		 */
		unsigned int resize(unsigned int newSize) {
			specIdx.resize(newSize);
			shifts.resize(newSize);
			consensus.resize(newSize);
			endpoints.resize(newSize);
			//for(unsigned int i=0;i<newSize;i++) { specIdx[i].resize(0); shifts[i].resize(0); consensus[i].resize(0); endpoints[i].resize(0); }
			return consensus.size();
		}

		/**
		 * TODO: add description
		 *
		 *@return
		 */
		unsigned int size() {
			return consensus.size();
		}

		/**
		 * TODO: add description
		 *
		 *@param newResolution
		 *@param enforceRounding
		 */
		void setResolution(float newResolution, bool enforceRounding);

		/**
     * Initializes this set of clusters from the abinfo_t data structure (see abruijn.h)
     *
     *@param contigSpectra Set of contig spectra
     *@param componentSets abinfo correspnding to contigSpectra that details which spectra and peaks
     *  were assembled into which abruijn vertices
     *@param starSpectra Set of star spectra that were assembled into contig spectra
     *  WARNING: addZPMpeaks will be called for starSpectra
     *@param pkTol peak tolerance in Da
     *@param pmTol parent mass tolerance in Da
     */
    void initializeFromAbinfo(SpecSet& contigSpectra,
                              abinfo_t& componentSets,
                              SpecSet& starSpectra,
                              float pkTol,
                              float pmTol);

		/**
		 * Determines the set of endpoints from a list of ABruijn vertices.
		 *
		 *@param idx Index of the cluster/consensus
		 *@param specs Set of assembled spectra
		 *@param abVertices Set of spectrum peaks (list, dim.2; each pair is (spectrum index, peak index)) in each ABruijn vertex (vector, dim.1)
		 *@param peakTol Peak mass tolerance (in Daltons)
		 *@param pmTol Parent mass tolerance (in Daltons)
		 */
		void set_endpoints(int idx, SpecSet &specs, vector<list<TwoValues<int> > > &abVertices, float peakTol, float pmTol);

		/**
		 * Generates the cluster/consensus spectrum assuming endpoints were b0/bN in assembled spectra
		 *
		 *@param idx Index of the cluster/consensus
		 *@param putHere Output spectrum
		 *@param peakTol Peak mass tolerance (in Daltons)
		 */
		void getSpecIfB(int idx, Spectrum &putHere, float peakTol=0, float ionOffset=0);

		/**
		 * Generates the cluster/consensus spectrum assuming endpoints were y0/yN in assembled spectra
		 *
		 *@param idx Index of the cluster/consensus
		 *@param putHere Output spectrum
		 *@param peakTol Peak mass tolerance (in Daltons)
		 */
		void getSpecIfY(int idx, Spectrum &putHere, float peakTol=0, float ionOffset=0);

		/**
		 * Loads clusters from a text file
		 *
		 *@param filename Name of the clusters file
		 *@return number of clusters (-1 if error)
		 */
		int Load(const char *contigFilename,
             const char *spectrumPklbinFilename = 0x0,
             const char *psmFilename = 0x0);

		/**
		 * Loads clusters from a .pklbin file
		 *
		 *@param filename Name of the clusters file
		 *@return number of clusters (-1 if error)
		 */
		int Load_pklbin(const char *filename,
                              const char *psmFilename = 0x0);


		/**
		 * Saves clusters to a text file
		 *
		 *@param filename Name of the clusters file
		 *@return number of clusters (-1 if error)
		 */
		int Save(const char *contigFilename,
			const char *spectrumPklbinFilename = 0x0,
			const char *psmFilename = 0x0);
	};
}

#endif
