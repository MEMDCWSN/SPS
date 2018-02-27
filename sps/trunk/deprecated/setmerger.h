#ifndef SETMERGER_H
#define SETMERGER_H

#include <vector>
#include <list>
#include <iostream>
#include "batch.h"
using namespace std;

/**
 * Helper class used when merging sets.
 */
class SetMerger {
	vector<int> freePosList;
	int firstFreePos;
	void freePos(int pos) {
		freePosList[pos] = firstFreePos;
		firstFreePos = pos;
	}
	void reassignSet(int setIdx1, int setIdx2);
public:

	/**
	 * TODO: add description
	 */
	vector<int> membership;

	/**
	 * TODO: add description
	 */
	vector<list<int> > sets;

	/**
	 * TODO: add description
	 */
	unsigned int numSets;

	/**
	 * TODO: add description
	 */
	vector<vector<Results_ASP> > cAlignsASP;

	/**
	 * TODO: add description
	 */
	vector<vector<Results_PA> > cAlignsPA;

	/**
	 * To allow going back from components to original order of alignments.
	 */
	vector<vector<int> > cAlignsASP_idx;

	/**
	 * TODO: add description
	 */
	vector<vector<int> > cAlignsPA_idx;

	/**
	 * TODO: add description
	 *
	 *@param count
	 */
	SetMerger(int count = 0) {
		resize(count);
	}

	/**
	 * TODO: add description
	 *
	 *@param eltIdx
	 *@return
	 */
	int createset(int eltIdx) {
		if (firstFreePos < 0) {
			cerr
					<< "SetMerger::create_set() - Not enough memory to create set for element "
					<< eltIdx << "!\n";
			exit(-1);
		}
		int setIdx = firstFreePos;
		firstFreePos = freePosList[firstFreePos];
		addElement(setIdx, eltIdx);
		numSets++;
		return setIdx;
	}

	/**
	 * TODO: add description
	 *
	 *@param setIdx
	 *@param eltIdx
	 */
	void addElement(int setIdx, int eltIdx) {
		sets[setIdx].push_front(eltIdx);
		membership[eltIdx] = setIdx;
	}

	/**
	 * TODO: add description
	 *
	 *@param setIdx
	 *@param eltIndices
	 */
	void removeElements(int setIdx, list<int> &eltIndices);

	/**
	 * TODO: add description
	 *
	 *@param maxNumElems
	 *@param minSetSize
	 *@param pairsASP
	 *@param pairsPA
	 */
	void createSets(unsigned int maxNumElems, unsigned int minSetSize, vector<
			Results_ASP> &pairsASP, vector<Results_PA> &pairsPA);
	//	void splitAligns(vector<Results_ASP> &pairsASP, vector<Results_PA> &pairsPA, vector<int> *indicesASP=0, vector<int> *indicesPA=0);

	/**
	 * TODO: add description
	 *
	 *@param pairsASP
	 *@param pairsPA
	 */
	void splitAligns(vector<Results_ASP> &pairsASP, vector<Results_PA> &pairsPA);

	/**
	 * TODO: add description
	 *
	 *@param pairs
	 *@param cAligns
	 *@param cAligns_idx
	 *@return
	 */
	template<class T> void splitAligns(vector<T> &pairs,
			vector<vector<T> > &cAligns, vector<vector<int> > &cAligns_idx);

	/**
	 * TODO: add description
	 *
	 *@param other
	 *@return
	 */
	bool spliceSet(SetMerger &other);


	/**
	 * TODO: add description
	 *
	 *@param numMembers
	 */
	void resize(unsigned int numMembers);

	/**
	 * TODO: add description
	 *
	 *@return
	 */
	unsigned int size() {
		return numSets;
	}

	/**
	 * TODO: add description
	 *
	 *@return
	 */
	unsigned int numElemsInSets() {
		unsigned int c = 0;
		for (unsigned int i = 0; i < membership.size(); i++)
			if (membership[i] >= 0)
				c++;
		return c;
	}

	/**
	 * TODO: add description
	 *
	 *@param setIdx1
	 *@param setIdx2
	 */
	void merge(int setIdx1, int setIdx2);

	/**
	 * TODO: add description
	 *
	 *@param setIdx
	 *@param elemsToKeep
	 */
	void splitSet(int setIdx, list<int> &elemsToKeep);

	/**
	 * Removes set with index setIdx.If removeElmtsOnly is true then
	 * the set is replaced with an empty set and its former elements become unassigned;
	 * otherwise the set is also deleted.
	 *
	 *@param setIdx Set index
	 *@param removeElmtsOnly If true then set setIdx becomes and empty set, otherwise the set is deleted.
	 */
	void removeSet(unsigned int setIdx, bool removeElmtsOnly=true);

	/**
	 * Removes all sets with less than minSetSize elements.
	 *
	 *@param minSetSize Minimum number of elements to retain a set
	 */
	void removeSmallSets(unsigned int minSetSize);

	/**
	 * TODO: add description
	 */
	void compressSetIndices();

	/**
	 * TODO: add description
	 *
	 *@param filename
	 *@return
	 */
	int saveas_binListArray(const char *filename);
};

/**
 * TODO: add description
 */
template<class T> class SetMerger2 {
public:

	/**
	 * TODO: add description
	 */
	int foo;
};

#endif
