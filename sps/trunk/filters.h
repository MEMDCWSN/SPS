#ifndef FILTER_H
/** TODO: add description */
#define FILTER_H

#include "spectrum.h"
#include "batch.h"

#include <cmath>
#include <iostream>

namespace specnets
{
	using namespace std;
	/**
	 * TODO: add description
	 *
	 *@param edge1
	 *@param edge2
	 *@param edge3
	 *@param pmTol
	 *@return
	 */
	inline bool ValidateTriangle(Results_SP &edge1, Results_SP &edge2, Results_SP &edge3, float pmTol) {
		return true;
	}

	/**
	 * TODO: add description
	 *
	 *@param edge1
	 *@param edge2
	 *@param edge3
	 *@param pmTol
	 *@return
	 */
	inline bool ValidateTriangle(Results_ASP &edge1, Results_ASP &edge2, Results_ASP &edge3, float pmTol) {
		Results_PA e1,e2,e3;

		e1.spec1=edge1.spec1;  e1.spec2=edge1.spec2;  e1.shift1=edge1.shift1; e1.shift2=0;
		e2.spec1=edge1.spec1;  e2.spec2=edge1.spec2;  e2.shift1=edge1.shift1; e2.shift2=0;
		e3.spec1=edge1.spec1;  e3.spec2=edge1.spec2;  e3.shift1=edge1.shift1; e3.shift2=0;
		return ValidateTriangle(e1,e2,e3,pmTol);
	}

	/**
	 * TODO: add description
	 *
	 *@param edge1
	 *@param edge2
	 *@param edge3
	 *@param pmTol
	 *@return
	 */
	inline bool ValidateTriangle(Results_PA &edge1, Results_PA &edge2, Results_PA &edge3, float pmTol) {

	/* Considers 3 possibilities for c1: edge1, edge2, edge3:
	 *          c1
	 *   O---------------->O
	 *    \------>O------>/
	 *       c2      c3
	 *  Note that one of the 3 possibilities _has_ to correspond to c1.
	*/
		Results_PA *short1, *short2;   // Pointers to the 2 short edges (c2,c3 above)
		float longShift,               // Value of the long shift (c1 above)
					curShift;

		// Find longShift, short1, short2
		longShift=abs(edge1.shift1)>abs(edge1.shift2)?edge1.shift1:edge1.shift2; short1=&edge2; short2=&edge3;
		curShift=abs(edge2.shift1)>abs(edge2.shift2)?edge2.shift1:edge2.shift2;
		if(abs(curShift)>abs(longShift)) { longShift=curShift; short1=&edge1; short2=&edge3; }
		curShift=abs(edge3.shift1)>abs(edge3.shift2)?edge3.shift1:edge3.shift2;
		if(abs(curShift)>abs(longShift)) { longShift=curShift; short1=&edge1; short2=&edge2; }

		// Test c1 = c2+c3
		if( abs(longShift - short1->shift1 - short2->shift1)<=pmTol or
				abs(longShift - short1->shift1 - short2->shift2)<=pmTol or
				abs(longShift - short1->shift2 - short2->shift1)<=pmTol or
				abs(longShift - short1->shift2 - short2->shift2)<=pmTol) return true;

		return false;
	}

	/**
	 * TODO: add description
	 *
	 *@param aligns
	 *@param idxStart
	 *@param idxEnd
	 *@param pmTol
	 *@param selectedIdx
	 *@return
	 */
	 /*
	 * FilterTriangles - Retains only edges forming a triangle (see ValidateTriangle
	 *   for an illustration of a triangle).
	 *
	 * Note: Assumes that:
	 *         - aligns is sorted by (.spec1,.spec2)
	 *         - contains no repeated edges
	 *         - for every aligns[i], aligns[i].spec1 < aligns[i].spec2
	 */
	template<class T> unsigned int FilterTriangles(vector<T> &aligns,
			unsigned int idxStart, unsigned int idxEnd, float pmTol, vector<
					unsigned int> &selectedIdx) {
		vector<bool> valid(aligns.size());  // Indicates whether an edge is part of some triangle
		vector<TwoValues<unsigned int> > adj;  // Start/end indices (in aligns) for each vertex's adjacency list
		unsigned int maxVertex=0,   // Largest vertex index in aligns
								 idxMain,       // Index of the edge being processed
								 idxPair;       // Index of the tentative pair for idxMain

		// Initialize adjacency lists
		for(idxMain=0; idxMain<aligns.size(); idxMain++) {
			valid[idxMain] = false;
			if((unsigned int)aligns[idxMain].spec1>maxVertex) maxVertex=aligns[idxMain].spec1;
			if((unsigned int)aligns[idxMain].spec2>maxVertex) maxVertex=aligns[idxMain].spec2;
		}
		adj.resize(maxVertex+1);
		for(unsigned int vertexIdx=0; vertexIdx<=maxVertex; vertexIdx++) adj[vertexIdx].set(1,0);  // Adjacency list is empty if start>end
		for(idxMain=0; idxMain<aligns.size(); idxMain++) {
			if(idxMain==0 or aligns[idxMain-1].spec1!=aligns[idxMain].spec1) adj[aligns[idxMain].spec1][0]=idxMain;
			adj[aligns[idxMain].spec1][1]=idxMain;
		}

	//cerr << "First edge ("<<aligns[0].spec1<<","<<aligns[0].spec2<<")\n";

		// Find triangles - for each edge (v1,v2), find a pairing edge (v1,v3) from the same vertex
		//  and look for a connector edge between the two destination vertices (v2,v3).
	//	pmTol *= 2;  // *2 = Once for c2 and once for c3
		TwoValues<unsigned int> searchBounds;  // lower/upper limits for the binary search for a connector edge
		unsigned int idxMiddle;                // index used in the binary search - middle index between searchBounds[0]/searchBounds[1]
		unsigned int targetVertex;             // vertex to find in the adjacency list
		idxEnd=min(idxEnd,(unsigned int)aligns.size());   idxEnd = max((unsigned int)0,idxEnd);   idxStart=max(idxStart,(unsigned int)0);

	//cerr<<"Looking for triangles with an edge between indices "<<idxStart<<" and "<<idxEnd<<"\n"; cerr.flush();
		for(idxMain=idxStart; idxMain<=idxEnd; idxMain++) {
			for(idxPair=adj[aligns[idxMain].spec1][0]; idxPair<adj[aligns[idxMain].spec1][1]; idxPair++) {
	//cerr << "Edge indices ("<<idxMain<<","<<idxPair<<"), edges ("<<aligns[idxMain].spec1<<"->"<<aligns[idxMain].spec2<<"), ("<<aligns[idxPair].spec1<<"->"<<aligns[idxPair].spec2<<")\n";
				if(idxPair==idxMain) continue;
				if(idxPair<idxMain) { searchBounds=adj[aligns[idxPair].spec2]; targetVertex=aligns[idxMain].spec2; }
				else { searchBounds=adj[aligns[idxMain].spec2]; targetVertex=aligns[idxPair].spec2; }

	//cerr << "Looking for targetVertex "<<aligns[idxPair].spec2<<" between indices "<<searchBounds[0]<<":"<<searchBounds[1]<<endl;
				// Binary search for targetVertex in the adjacency defined by searchBounds
				while(searchBounds[0]<=searchBounds[1]) {
	//cerr << " --- searchBounds = ("<<searchBounds[0]<<","<<searchBounds[1]<<")\n";
					if((unsigned int)aligns[searchBounds[1]].spec2==targetVertex) searchBounds[0]=searchBounds[1];
					if((unsigned int)aligns[searchBounds[0]].spec2==targetVertex or searchBounds[1]-searchBounds[0]<2) break;
					idxMiddle = max(searchBounds[0]+1,(unsigned int)round((searchBounds[0]+searchBounds[1])/2.0));
					if((unsigned int)aligns[idxMiddle].spec2<=targetVertex) searchBounds[0]=idxMiddle; else searchBounds[1]=idxMiddle;
				}
				if((unsigned int)aligns[searchBounds[0]].spec2!=targetVertex) continue;  // Connector edge not found

	//cerr<<"Found targetVertex at "<<searchBounds[0]<<", edge ("<<aligns[searchBounds[0]].spec1<<","<<aligns[searchBounds[0]].spec2<<")\n";

				if(ValidateTriangle(aligns[idxMain], aligns[idxPair], aligns[searchBounds[0]], pmTol)) {
					valid[idxMain] = true;   valid[idxPair] = true;   valid[searchBounds[0]] = true;
	//cerr<<"Validated triangle: idxMain="<<idxMain<<", idxPair="<<idxPair<<", searchBounds[0]="<<searchBounds[0]<<"\n"; cerr.flush();
				}
			}

			if((idxMain%1000)==0)
				{ cout<<"Done processing edge "<<idxMain<<" ("<<((double)idxMain-idxStart)/((double)idxEnd-idxStart)<<"% completed)\n"; cout.flush(); }
		}

		// Remove edges that do not participate in a triangle
		unsigned int idxKept=0; // Index to the last kept alignment
		for(idxMain=0; idxMain<aligns.size(); idxMain++)
			if(valid[idxMain]) aligns[idxKept++] = aligns[idxMain];
		aligns.resize(idxKept);   selectedIdx.resize(idxKept);   idxKept=0;
		for(idxMain=0; idxMain<valid.size(); idxMain++)
			if(valid[idxMain]) selectedIdx[idxKept++]=idxMain;

		return idxKept;  // Number of retained alignments
	}

	/**
	 * TODO: add description
	 *
	 *@param aligns
	 *@param idxKept
	 *@param pvalues
	 *@param gaussianParams
	 *@param ratios
	 *@param minPValue
	 *@param minRatio
	 *@param pmTol
	 *@param filterTrigs
	 *@return
	 */
	template<class T> unsigned int FilterAligns(vector<T> &aligns, vector<
			unsigned int> &idxKept, vector<TwoValues<float> > &pvalues, vector<
			TwoValues<float> > &gaussianParams, vector<vector<float> > &ratios,
			float minPValue, float minRatio, float pmTol, bool filterTrigs) {
		// Compute pvalues and filter by ratios and p-values
		unsigned int idxPair, idxLast=0;   idxKept.resize(aligns.size());   pvalues.resize(aligns.size());
		for(idxPair=0; idxPair<aligns.size(); idxPair++) {
			pvalues[idxPair][0] = 1 - Utils::gaussiancdf(aligns[idxPair].score1,gaussianParams[aligns[idxPair].spec1][0],gaussianParams[aligns[idxPair].spec1][1]);
			pvalues[idxPair][1] = 1 - Utils::gaussiancdf(aligns[idxPair].score2,gaussianParams[aligns[idxPair].spec2][0],gaussianParams[aligns[idxPair].spec2][1]);

			if(pvalues[idxPair][0]<=minPValue and pvalues[idxPair][1]<=minPValue and
					(ratios.size()==0 or (ratios[idxPair][0]>=minRatio and ratios[idxPair][1]>=minRatio)) )
				{ if(idxLast<idxPair) aligns[idxLast]=aligns[idxPair]; idxKept[idxLast]=idxPair; idxLast++; }
		}
		aligns.resize(idxLast);   idxKept.resize(idxLast);

		vector<unsigned int> idxKeptT;
		if(filterTrigs) {
			FilterTriangles(aligns, 0, aligns.size(),pmTol, idxKeptT);
			idxLast=0;
			for(unsigned int idxPair=0; idxPair<idxKeptT.size(); idxPair++) {
				if(idxLast<idxPair) { aligns[idxLast]=aligns[idxPair]; idxKept[idxLast]=idxKeptT[idxPair]; }
				idxLast++;
			}
			aligns.resize(idxLast);   idxKept.resize(idxLast);
		}
		return idxKept.size();
	}

	/**
	 * TODO: add description
	 *
	 *@param specSet
	 *@param results
	 *@param ratioType
	 *@param ratios
	 */
	void FilterMatchRatioASP(SpecSet &specSet, vector<Results_ASP> &results,
			short ratioType, vector<TwoValues<float> > &ratios);

	/**
	 * TODO: add description
	 *
	 *@param specSet
	 *@param results
	 *@param scoreMeans is mean (col.0) and sample size (col.1)
	 */void ComputeScoreMeans(SpecSet &specSet, vector<Results_ASP> &results, vector<
			TwoValues<float> > &scoreMeans);
}
#endif
