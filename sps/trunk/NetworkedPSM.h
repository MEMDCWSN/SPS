/*
 * NetworkedPSM.h
 *
 *  Created on: Sep 23, 2015
 *      Author: isna
 */

#ifndef NETWORKEDPSM_H_
#define NETWORKEDPSM_H_

#include <string>
#include <vector>

#include "spectrum.h"

using namespace std;

namespace specnets {

	class NetworkedPSM {
		public:
			NetworkedPSM();
			NetworkedPSM(string annotation, int seedIndex);
			NetworkedPSM(string annotation, int seedIndex, int layer, float matchScore, int parent, float diffMass, float alignScore);
			virtual ~NetworkedPSM();

			NetworkedPSM* propagateToNeighborSpectrum(int layer, int parent, float diffMass, float alignScore, Spectrum & neighSpec, float fragTol);

			string toString();
			int compareTo(NetworkedPSM p);

			int getLayer(){ return m_layer; }
			int getSeedIndex(){ return m_seedSpecIndex; }
			string getAnnotation(){ return m_annotation; }
			float getAlignmentScore() { return m_matchScore; }//m_alignmentScore;

		private:
			string 			m_annotation;
			int 			m_modCount;
			float 			m_absModSum;
			int 			m_layer;
			float			m_matchScore; //peptide vs. PRM spectrum

			int				m_seedSpecIndex;	//seed index, zero-based
			int				m_parentSpecIndex;  //precedence spec index
			float			m_deltaFromParent;  //this - precedence
			float			m_alignmentScore;	//aligngf pvalue

			char 			m_wrappingNtermAA;
			char 			m_wrappingCtermAA;
			string 			m_stripSeq;
			vector<float> 	m_prmSpec;
			vector<float> 	m_modMass;
			float			m_neutralMass;

			void construct();
		//	vector<float> getModifiedPrmSpec(int modSite, float modMass);
		//	string getModifiedAnnotation(int modSite, float modMass);

			void getSubPrmSpec(vector<float> & subPrms, vector<float> & subMods, int start, int end);
			string getSubModPept(int start, int end, int modSite, float modMass);

			NetworkedPSM getBestAlignment(Spectrum & neighSpec, bool ntermRemovable, bool ctermRemovable, float fragTol);

			string getUnmodifiedPept();
			string getPeptIncludingNtermAA();
			string getPeptIncludingCtermAA();
			string getSeqIncludingNtermAA();
			string getSeqIncludingCtermAA();
	};

} /* namespace specnets */
#endif /* NETWORKEDPSM_H_ */
