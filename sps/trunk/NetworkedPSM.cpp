/*
 * NetworkedPSM.cpp
 *
 *  Created on: Sep 23, 2015
 *      Author: isna
 */

#include "NetworkedPSM.h"
#include "aminoacid.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sstream>

namespace specnets {

	static const double  isotope_space = 1.00235;
/*
	static const double	Hydrogen = 1.007825035;
	static const double	Oxygen = 15.99491463;
	static const double	H2O = Hydrogen*2 + Oxygen;
	static const double	Proton = 1.007276035;
	static const double	H2Oplus = H2O+Proton;//*/

	static const float min_distance = 54;
	static const float null_score = -10000;
	static const float ptm_penalty = 5;

	static float getAAMass(char aa){
		static float org_aa_masses[]={
			71.03711, 0, 103.00919, 115.02694, 129.04259,
			147.06841, 57.02146, 137.05891, 113.08406, 0,
			128.09496, 113.08406, 131.04049, 114.04293, 0,
			97.05276, 128.05858, 156.10111, 87.03203, 101.04768,
			0, 99.06841, 186.07931, 0, 163.06333, 0};

		return org_aa_masses[aa-'A'];
	}

	static float getMatchScore(Spectrum spec, vector<float> prms, float fragTol){

		if( prms.size() == 0 ) return null_score;

		float score = 0;

		int start = 0;
		for( int p=0; p<prms.size()-1; p++ ){
			for (int i=start; i<spec.size(); i++){
				if( spec[i][0] < prms[p]-fragTol ) continue;
				else if( spec[i][0] > prms[p]+fragTol ) {
					start = i;
					break;
				}
				score += spec[i][1];
			}
		}
		score += prms[prms.size()-1];//consider mod penalty

		return score;
	}

	static vector<float> getModPrms(vector<float> prms, vector<float> mods, int site, float delta){

		float prevM = (site==0)? 0 : prms[site-1];
		if( (prms[site]-prevM) + delta < min_distance ) return vector<float>(0);

		vector<float> mp(prms.size());
		for(int i=0; i<site; i++) mp[i] = prms[i];
		for(int i=site; i<prms.size()-1; i++) {
			mp[i] = prms[i]+delta;
		}

		if( round(delta) != 0 ){
			if( mods[site] != 0 ){ //introduce new mod into already modified site
				if( round(mods[site]+delta) == 0 ) {
					mp[prms.size()-1] = prms[prms.size()-1] + ptm_penalty;
				}
			}
			else mp[prms.size()-1] = prms[prms.size()-1] - ptm_penalty; //introduce new mod into new site
		}

		return mp;
	}

	static pair<int, float> localizeMOD(Spectrum spec, vector<float> prms, vector<float> mods, float delta, float fragTol){

		pair<int, float> site_score;
		int site = -1;
		double max = null_score;

		for(int i=0; i<prms.size(); i++){
			vector<float> mp = getModPrms(prms, mods, i, delta);
			if( mp.size() == 0 ) continue;
			float score = getMatchScore(spec, mp, fragTol);
			if( max < score ) {
				max = score;
				site = i;
			}
		}
		site_score.first = site;
		site_score.second = max;

		return site_score;
	}

	NetworkedPSM::NetworkedPSM() : m_layer(-1), m_matchScore(null_score), m_absModSum(0)
	{
	}

	NetworkedPSM::~NetworkedPSM()
	{
	}

	NetworkedPSM::NetworkedPSM(string annotation, int seedIndex) :
			m_annotation(annotation), m_seedSpecIndex(seedIndex), m_layer(0),
			m_matchScore(0), m_parentSpecIndex(-1), m_deltaFromParent(0), m_alignmentScore(0)
	{ // this is for seed id
		construct();
	}

	NetworkedPSM::NetworkedPSM(string annotation, int seedIndex, int layer, float matchScore, int parent, float diffMass, float alignScore) :
				m_annotation(annotation), m_seedSpecIndex(seedIndex), m_layer(layer),
				m_matchScore(matchScore), m_parentSpecIndex(parent), m_deltaFromParent(diffMass), m_alignmentScore(alignScore)
	{
		construct();
	}

	void NetworkedPSM::construct(){

		m_wrappingNtermAA = m_annotation[0];
		m_stripSeq = "";
		m_wrappingCtermAA = m_annotation[m_annotation.length()-1];

		m_modCount = 0;
		string 			delta;
		vector<int> 	msite;
		vector<string> 	msdelta;

		bool prevAA = true;
		for(int i=2; i<m_annotation.length()-2; i++){
			if( isalpha(m_annotation[i]) ){
				m_stripSeq += m_annotation[i];
				if(!prevAA){
					msdelta.push_back(delta);
					delta.clear();
				}
				prevAA = true;
			}
			else if(prevAA) {
				m_modCount++;
				if( m_stripSeq.length() == 0 ) msite.push_back(0);
				else msite.push_back(m_stripSeq.length()-1);
				prevAA = false;
				delta += m_annotation[i];
			}
			else if(!prevAA) {
				if( m_annotation[i] == '+' || m_annotation[i] == '-' ){
					msdelta.push_back(delta);
					delta.clear();
					msite.push_back(msite[msite.size()-1]);
				}
				delta += m_annotation[i];
			}
		}
		if(!prevAA){
			msdelta.push_back(delta);
		}

		m_prmSpec.resize(m_stripSeq.length());
		m_modMass.resize(m_stripSeq.length());
		m_absModSum = 0;

		if( m_modCount > 0 ){
			for(int k=0; k<msite.size(); k++){
				m_modMass[msite[k]] += atof(msdelta[k].c_str());
				m_absModSum += m_modMass[msite[k]];
			}
		}
		m_absModSum = fabs(m_absModSum);

		float prm_mass = 0;
		for(int i=0; i<m_stripSeq.length(); i++){
		  prm_mass += getAAMass(m_stripSeq[i]) + m_modMass[i];
		//  prm_mass += org_aa_masses[m_stripSeq[i]-'A'] + m_modMass[i];
		  m_prmSpec[i] = prm_mass;
		}
		m_neutralMass = prm_mass;
		m_prmSpec[m_prmSpec.size()-1] = -m_modCount * ptm_penalty;
	}

	string NetworkedPSM::toString(){
		stringstream ss;
		ss	<<	m_layer	<<"\t"<<
				m_matchScore <<"\t"<<
				m_annotation <<"\t"<<
				m_modCount <<"\t"<<
				(m_seedSpecIndex+1) <<"\t"<<
				(m_parentSpecIndex+1) <<"\t"<<
				m_deltaFromParent <<"\t"<<
				m_alignmentScore;
		return ss.str();
	}

	string NetworkedPSM::getUnmodifiedPept(){
		stringstream ss;
		ss	<<	m_wrappingNtermAA << "." << m_stripSeq << "." << m_wrappingCtermAA;
		return ss.str();
	}

	string NetworkedPSM::getSeqIncludingNtermAA(){
		stringstream ss;
		ss	<<	"*." << m_wrappingNtermAA << m_stripSeq << "." << m_wrappingCtermAA;
		return ss.str();
	}
	string NetworkedPSM::getSeqIncludingCtermAA(){
		stringstream ss;
		ss	<<	m_wrappingNtermAA << "." << m_stripSeq << m_wrappingCtermAA << ".*";
		return ss.str();
	}

	string NetworkedPSM::getPeptIncludingNtermAA(){
		stringstream ss;
		ss	<<	"*." << m_wrappingNtermAA << m_annotation.substr(2);
		return ss.str();
	}
	string NetworkedPSM::getPeptIncludingCtermAA(){
		stringstream ss;
		ss	<<	m_annotation.substr(0, m_annotation.length()-2) << m_wrappingCtermAA << ".*";
		return ss.str();
	}

	int NetworkedPSM::compareTo(NetworkedPSM p){

		if( m_matchScore > p.m_matchScore ) return -1;
		else if( m_matchScore < p.m_matchScore ) return 1;

		if( m_absModSum < p.m_absModSum ) return -1;
		else if( m_absModSum > p.m_absModSum ) return 1;

		return 0;
	}

	void NetworkedPSM::getSubPrmSpec(vector<float> & subPrms, vector<float> & subMods, int start, int end) {

		subPrms.resize(end-start);
		subMods.resize(end-start);

		float shiftMass = (start==0)? 0 : m_prmSpec[start-1];
		int modC = 0;

		for(int i=start; i<end-1; i++) subPrms[i-start] = m_prmSpec[i]-shiftMass;
		for(int i=start; i<end; i++) {
			subMods[i-start] = m_modMass[i];
			if( m_modMass[i] != 0 ) modC++;
		}
		subPrms[end-start-1] = -modC*ptm_penalty;
	}

	string NetworkedPSM::getSubModPept(int start, int end, int modSite, float modMass){

		string subSeq = m_stripSeq.substr(start, end-start);
		char prev = (start==0)? m_wrappingNtermAA : m_stripSeq[start-1];
		char next = (end==m_stripSeq.length())? m_wrappingCtermAA : m_stripSeq[end];

		stringstream ss;

		ss << prev << ".";
		for(int i=0; i<subSeq.length(); i++){
			ss << subSeq[i];

			float m = m_modMass[i+start];
			if( i == modSite ) m += modMass;

			if( round(m) != 0 ) {
				if( m > 0 ) ss << "+";
				char temp[10];
				sprintf(temp, "%.3f", m);
				ss << temp;
			}
		}
		ss << "." << next;

		return ss.str();
	}


	NetworkedPSM* NetworkedPSM::propagateToNeighborSpectrum(int layer, int parent, float observedMassDiff, float alignScore, Spectrum & neighSpec, float fragTol){
		if( m_layer == -1 ) cout << "WHAT!!! IMPOSSIBLE" << endl;

		NetworkedPSM bestPSM = getBestAlignment(neighSpec, true, true, fragTol);

		vector<NetworkedPSM> extraPool;
		if( isalpha(m_wrappingNtermAA) && ( -5 < observedMassDiff-getAAMass(m_wrappingNtermAA) ) )
			extraPool.push_back(NetworkedPSM(getPeptIncludingNtermAA(), m_seedSpecIndex).getBestAlignment(neighSpec, false, false, fragTol));

		if( isalpha(m_wrappingCtermAA) && ( -5 < observedMassDiff-getAAMass(m_wrappingCtermAA) ) )
			extraPool.push_back(NetworkedPSM(getPeptIncludingCtermAA(), m_seedSpecIndex).getBestAlignment(neighSpec, false, false, fragTol));

		if( m_modCount != 0 ) {

			if( m_absModSum < 5 ) {
				extraPool.push_back(NetworkedPSM(getUnmodifiedPept(), m_seedSpecIndex).getBestAlignment(neighSpec, false, false, fragTol));
			}
			else if( 50 < m_absModSum ){
				if( isalpha(m_wrappingNtermAA) )
					extraPool.push_back(NetworkedPSM(getSeqIncludingNtermAA(), m_seedSpecIndex).getBestAlignment(neighSpec, false, true, fragTol));

				if( isalpha(m_wrappingCtermAA) )
					extraPool.push_back(NetworkedPSM(getSeqIncludingCtermAA(), m_seedSpecIndex).getBestAlignment(neighSpec, true, false, fragTol));
			}
		}

		for(int i=0; i<extraPool.size(); i++){
			if( bestPSM.compareTo(extraPool[i]) == 1 ) bestPSM = extraPool[i];
		}

		return new NetworkedPSM( bestPSM.m_annotation, bestPSM.m_seedSpecIndex, layer, bestPSM.m_matchScore, parent, observedMassDiff, alignScore );
	}

	NetworkedPSM NetworkedPSM::getBestAlignment(Spectrum & neighSpec, bool ntermRemovable, bool ctermRemovable, float fragTol){

		float massDiff = neighSpec.parentMass - (m_neutralMass + AAJumps::massMH);

		vector<vector<float> > subPrmSpectra;
		vector<vector<float> > subModArrays;
		vector<pair<int, int> > subPeptSites;
		vector<float > subDeltas;

		subPrmSpectra.push_back(m_prmSpec); // add normal form
		subModArrays.push_back(m_modMass);
		subPeptSites.push_back(pair<int, int>(0, m_stripSeq.length()));
		subDeltas.push_back(massDiff);

		if( ntermRemovable ) {

			int subStart = 0;
			float subDelta = massDiff + m_prmSpec[subStart], prevDelta = massDiff;

			while( fabs(subDelta) < fabs(prevDelta) ){
				subStart++;

				vector<float> subPrms, subMods;
				getSubPrmSpec(subPrms, subMods, subStart, m_stripSeq.length());

				subPrmSpectra.push_back(subPrms);
				subModArrays.push_back(subMods);
				subPeptSites.push_back(pair<int, int>(subStart, m_stripSeq.length()));
				subDeltas.push_back(subDelta);

				prevDelta = subDelta;
				subDelta = massDiff + m_prmSpec[subStart];
			}
		}

		if( ctermRemovable ) {

			int subEnd = m_stripSeq.length()-2;
			float subDelta = massDiff + (m_neutralMass-m_prmSpec[subEnd]), prevDelta = massDiff;

			while( fabs(subDelta) < fabs(prevDelta) ){
				subEnd--;

				vector<float> subPrms, subMods;
				getSubPrmSpec(subPrms, subMods, 0, subEnd+2);

				subPrmSpectra.push_back(subPrms);
				subModArrays.push_back(subMods);
				subPeptSites.push_back(pair<int, int>(0, subEnd+2));
				subDeltas.push_back(subDelta);

				prevDelta = subDelta;
				subDelta = massDiff + (m_neutralMass-m_prmSpec[subEnd]);
			}
		}

//		int corr[] = {0, -1, -2, -3, 1, 2, 3};
		int corr[] = {0};

		float maxScore = null_score;
		int modSite = -1;
		float fixedDelta = 0;
		int	subIndex = -1;

		for(int i=0; i<1; i++) {
			for(int k=0; k<subPrmSpectra.size(); k++){

				float corrDelta = subDeltas[k] + corr[i]*isotope_space;
				pair<int, float> match= localizeMOD(neighSpec, subPrmSpectra[k], subModArrays[k], corrDelta, fragTol);
				if( match.first != -1 ){
					if( maxScore < match.second || ( maxScore == match.second && fabs(corrDelta) < fabs(fixedDelta) ) ){
						maxScore = match.second;
						modSite = match.first;
						fixedDelta = corrDelta;
						subIndex = k;
					}
				}
			}
		}

		if( modSite == -1 ) return NetworkedPSM();

		NetworkedPSM result = NetworkedPSM( getSubModPept(subPeptSites[subIndex].first, subPeptSites[subIndex].second, modSite, fixedDelta), m_seedSpecIndex);
		result.m_matchScore = maxScore;
		return result;
	}



	/*
	NetworkedPSM NetworkedPSM::propagateToNeighborSpectrum(int layer, int parent, float diffMass, float alignScore, Spectrum & neighSpec, float fragTol) {

		if( round(diffMass) == 0 )
			return NetworkedPSM( m_annotation, layer, getMatchScore(neighSpec, m_prmSpec, fragTol), parent, diffMass, alignScore );


		static int corr[] = {0, -1, 1, -2, -3};
		double maxScore = nullScore;
		int site = -1;
		float newDelta = -1;

	//	for(int corr=0; corr>-3; corr--) {
	//		float corrDelta = diffMass + corr*IsotopeSpace;
		for(int i=0; i<5; i++){
			float corrDelta = diffMass + corr[i]*IsotopeSpace;
			for(int i=0; i<m_modMass.size(); i++){
				double score = getMatchScore(neighSpec, getModifiedPrmSpec(i, corrDelta), fragTol);

				if( maxScore < score ) {
					maxScore = score;
					site = i;
					newDelta = corrDelta;
				}
			}
		}

		if( site == -1 ) return NetworkedPSM();
		return NetworkedPSM( getModifiedAnnotation(site, newDelta), layer, maxScore, parent, diffMass, alignScore );
	}

	vector<float> NetworkedPSM::getModifiedPrmSpec(int modSite, float modMass){

		if( round(m_prmSpec[modSite] + modMass) < 0 ) {
			vector<float> modifiedPS(0);
			return modifiedPS;
		}

		vector<float> modifiedPS(m_prmSpec.size());

		for(int i=0; i<modSite; i++) modifiedPS[i] = m_prmSpec[i];
		for(int i=modSite; i<m_prmSpec.size()-1; i++) {
			modifiedPS[i] = m_prmSpec[i]+modMass;
		}

		if( round(modMass) != 0 ){
			if( m_modMass[modSite] != 0 ){
				if( round(m_modMass[modSite]+modMass) == 0 ) {
					modifiedPS[m_prmSpec.size()-1] = m_prmSpec[m_prmSpec.size()-1] + ptmPenalty;
				}
			}
			else  modifiedPS[m_prmSpec.size()-1] = m_prmSpec[m_prmSpec.size()-1] - ptmPenalty;
		}

		return modifiedPS;
	}

	string NetworkedPSM::getModifiedAnnotation(int modSite, float modMass){

		stringstream ss;

		ss << m_annotation[0] << ".";
		for(int i=0; i<m_stripSeq.length(); i++){
			ss << m_stripSeq[i];

			float m = m_modMass[i];
			if( i == modSite ) m += modMass;

			if( round(m) != 0 ) {
				if( m > 0 ) ss << "+";
				char temp[10];
				sprintf(temp, "%.3f", m);
				ss << temp;
			}
		}
		ss << "." << m_annotation[m_annotation.length()-1];

		return ss.str();
	}//*/



} /* namespace specnets */





