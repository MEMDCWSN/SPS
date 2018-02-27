#ifndef SETMERGER_H
#define SETMERGER_H

// Module Includes
#include "SpectrumPairSet.h"

// System Includes
#include <iostream>
#include <list>
#include <stdlib.h>
#include <vector>

/**
 * Helper class used when merging sets.
 */
namespace specnets
{
  class SetMerger
  {
  public:

    /**
     * TODO: add description
     *
     *@param count
     */
    SetMerger(int count = 0);

    /**
     * TODO: add description
     *
     *@param eltIdx
     *@return
     */
    int createset(int eltIdx);

    /**
     * TODO: add description
     *
     *@param setIdx
     *@param eltIdx
     */
    void addElement(int setIdx, int eltIdx);

    /**
     * TODO: add description
     *
     *@param setIdx
     *@param eltIndices
     */
    void removeElements(int setIdx, std::list<int> &eltIndices);

    /**
     * TODO: add description
     *
     *@param maxNumElems
     *@param minSetSize
     *@param pairsASP
     *@param pairsPA
     */
    void createSets(unsigned int maxNumElems,
                    unsigned int minSetSize,
                    SpectrumPairSet &pairsASP,
                    SpectrumPairSet &pairsPA);
    //	void splitAligns(vector<Results_ASP> &pairsASP, vector<Results_PA> &pairsPA, vector<int> *indicesASP=0, vector<int> *indicesPA=0);

    /**
     * TODO: add description
     *
     *@param pairsASP
     *@param pairsPA
     */
    void splitAligns(SpectrumPairSet &pairsASP, SpectrumPairSet &pairsPA);

    /**
     * TODO: add description
     *
     *@param pairs
     *@param cAligns
     *@param cAligns_idx
     *@return
     */
    void splitAligns(SpectrumPairSet &pairs,
                     std::vector<SpectrumPairSet> &cAligns,
                     std::vector<std::vector<int> > &cAligns_idx);

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
    unsigned int size();

    /**
     * TODO: add description
     *
     *@return
     */
    unsigned int numElemsInSets();

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
    void splitSet(int setIdx, std::list<int> &elemsToKeep);

    /**
     * Removes set with index setIdx.If removeElmtsOnly is true then
     * the set is replaced with an empty set and its former elements become unassigned;
     * otherwise the set is also deleted.
     *
     *@param setIdx Set index
     *@param removeElmtsOnly If true then set setIdx becomes and empty set, otherwise the set is deleted.
     */
    void removeSet(unsigned int setIdx, bool removeElmtsOnly = true);

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
    int saveas_binListArray(char *filename);

    /**
     * TODO: add description
     */
    std::vector<int> membership;

    /**
     * TODO: add description
     */
    std::vector<std::list<int> > sets;

    /**
     * TODO: add description
     */
    unsigned int numSets;

    /**
     * TODO: add description
     */
    std::vector<SpectrumPairSet> cAlignsASP;

    /**
     * TODO: add description
     */
    std::vector<SpectrumPairSet> cAlignsPA;

    /**
     * To allow going back from components to original order of alignments.
     */
    std::vector<std::vector<int> > cAlignsASP_idx;

    /**
     * TODO: add description
     */
    std::vector<std::vector<int> > cAlignsPA_idx;

  private:
    /**
     * TODO: add description
     */
    void freePos(int pos);

    /**
     * TODO: add description
     */
    void reassignSet(int setIdx1, int setIdx2);

    std::vector<int> freePosList;
    int firstFreePos;

  };

} // namespace specnets

#endif
