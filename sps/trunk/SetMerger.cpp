// Header Include
#include "SetMerger.h"
#include "utils.h"

// System Includes
#include <iostream>

using namespace std;

namespace specnets
{
  //---------------------------------------------------------------------------
  SetMerger::SetMerger(int count)
  {
    resize(count);
  }

  //---------------------------------------------------------------------------
  int SetMerger::createset(int eltIdx)
  {
    if (firstFreePos < 0)
    {
      std::cerr
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

  //---------------------------------------------------------------------------
  void SetMerger::addElement(int setIdx, int eltIdx)
  {
    sets[setIdx].push_front(eltIdx);
    membership[eltIdx] = setIdx;
  }

  //---------------------------------------------------------------------------
  void SetMerger::removeElements(int setIdx, list<int> &eltIndices)
  {
    if (setIdx < 0 or setIdx >= (int)sets.size())
      return;
    eltIndices.sort();

    list<int>::iterator setIter = sets[setIdx].begin(), rmIter =
        eltIndices.begin();
    while (setIter != sets[setIdx].end() and rmIter != eltIndices.end())
      if (*setIter == *rmIter)
      {
        setIter = sets[setIdx].erase(setIter);
        rmIter++;
      }
      else
      {
        if (*setIter < *rmIter)
          setIter++;
        else
          rmIter++;
      }
    for (rmIter = eltIndices.begin(); rmIter != eltIndices.end(); rmIter++)
      if (*rmIter >= 0 and *rmIter < (int)membership.size())
        membership[*rmIter] = -1;
  }

  //---------------------------------------------------------------------------
  void SetMerger::createSets(unsigned int maxNumElems,
                             unsigned int minSetSize,
                             SpectrumPairSet &pairsASP,
                             SpectrumPairSet &pairsPA)
  {
    unsigned int pivot;
    resize(maxNumElems);
    for (pivot = 0; pivot < maxNumElems; pivot++)
      createset(pivot);
    for (pivot = 0; pivot < pairsASP.size(); pivot++)
    {
      if ((unsigned int)pairsASP[pivot].spec1 < maxNumElems
          and (unsigned int)pairsASP[pivot].spec2 < maxNumElems)
        merge(membership[pairsASP[pivot].spec1],
              membership[pairsASP[pivot].spec2]);
    }

    for (pivot = 0; pivot < pairsPA.size(); pivot++)
      if ((unsigned int)pairsPA[pivot].spec1 < maxNumElems
          and (unsigned int)pairsPA[pivot].spec2 < maxNumElems)
        merge(membership[pairsPA[pivot].spec1],
              membership[pairsPA[pivot].spec2]);

    removeSmallSets(minSetSize);
    compressSetIndices();
  }

  //---------------------------------------------------------------------------
  void SetMerger::splitAligns(SpectrumPairSet &pairsASP,
                              SpectrumPairSet &pairsPA)
  {
    // Count number of alignments per component
    vector<int> cAlignsASP_counts(sets.size());
    for (unsigned int i = 0; i < cAlignsASP_counts.size(); i++)
      cAlignsASP_counts[i] = 0;
    vector<int> cAlignsPA_counts(sets.size());
    for (unsigned int i = 0; i < cAlignsPA_counts.size(); i++)
      cAlignsPA_counts[i] = 0;
    for (unsigned int i = 0; i < pairsASP.size(); i++)
      if (pairsASP[i].spec1 >= 0 and pairsASP[i].spec1 < membership.size()
          and membership[pairsASP[i].spec1] >= 0
          and membership[pairsASP[i].spec1] < sets.size())
        cAlignsASP_counts[membership[pairsASP[i].spec1]]++;
    for (unsigned int i = 0; i < pairsPA.size(); i++)
      if (pairsPA[i].spec1 >= 0 and pairsPA[i].spec1 < membership.size()
          and membership[pairsPA[i].spec1] >= 0
          and membership[pairsPA[i].spec1] < sets.size())
        cAlignsPA_counts[membership[pairsPA[i].spec1]]++;

    // Resize data structures
    cAlignsASP.resize(sets.size());
    cAlignsASP_idx.resize(sets.size());
    cAlignsPA.resize(sets.size());
    cAlignsPA_idx.resize(sets.size());
    for (unsigned int i = 0; i < cAlignsASP.size(); i++)
    {
      cAlignsASP[i].resize(cAlignsASP_counts[i]);
      cAlignsASP_idx[i].resize(cAlignsASP_counts[i]);
      cAlignsASP_counts[i] = 0;
    }
    for (unsigned int i = 0; i < cAlignsPA.size(); i++)
    {
      cAlignsPA[i].resize(cAlignsPA_counts[i]);
      cAlignsPA_idx[i].resize(cAlignsPA_counts[i]);
      cAlignsPA_counts[i] = 0;
    }

    // Split the alignments per component
    for (unsigned int i = 0; i < pairsASP.size(); i++)
      if (pairsASP[i].spec1 >= 0 and pairsASP[i].spec1 < membership.size()
          and membership[pairsASP[i].spec1] >= 0
          and membership[pairsASP[i].spec1] < sets.size())
      {
        int cIdx = membership[pairsASP[i].spec1], pairIdx =
            cAlignsASP_counts[cIdx]++;
        cAlignsASP[cIdx][pairIdx] = pairsASP[i];
        cAlignsASP_idx[cIdx][pairIdx] = i;
      }

    for (unsigned int i = 0; i < pairsPA.size(); i++)
      if (pairsPA[i].spec1 >= 0 and pairsPA[i].spec1 < membership.size()
          and membership[pairsPA[i].spec1] >= 0
          and membership[pairsPA[i].spec1] < sets.size())
      {
        int cIdx = membership[pairsPA[i].spec1], pairIdx =
            cAlignsPA_counts[cIdx]++;
        cAlignsPA[cIdx][pairIdx] = pairsPA[i];
        cAlignsPA_idx[cIdx][pairIdx] = i;
      }
  }

  //---------------------------------------------------------------------------
  void SetMerger::splitAligns(SpectrumPairSet &pairs,
                              vector<SpectrumPairSet> &cAligns,
                              vector<vector<int> > &cAligns_idx)
  {
    // Count number of alignments per component
    vector<int> cAligns_counts(sets.size());
    for (unsigned int i = 0; i < cAligns_counts.size(); i++)
      cAligns_counts[i] = 0;
    for (unsigned int i = 0; i < pairs.size(); i++)
      if (pairs[i].spec1 >= 0 && pairs[i].spec1 < membership.size()
          && membership[pairs[i].spec1] >= 0 && membership[pairs[i].spec1]
          < sets.size())
        cAligns_counts[membership[pairs[i].spec1]]++;

    // Resize data structures
    cAligns.resize(sets.size());
    cAligns_idx.resize(sets.size());
    for (unsigned int i = 0; i < cAligns.size(); i++)
    {
      cAligns[i].resize(cAligns_counts[i]);
      cAligns_idx[i].resize(cAligns_counts[i]);
      cAligns_counts[i] = 0;
    }
    // Split the alignments per component
    for (unsigned int i = 0; i < pairs.size(); i++)
      if (pairs[i].spec1 >= 0 and pairs[i].spec1 < membership.size()
          and membership[pairs[i].spec1] >= 0 and membership[pairs[i].spec1]
          < sets.size())
      {
        int cIdx = membership[pairs[i].spec1], pairIdx = cAligns_counts[cIdx]++;
        cAligns[cIdx][pairIdx] = pairs[i];
        cAligns_idx[cIdx][pairIdx] = i;
      }
  }

  //---------------------------------------------------------------------------
  bool SetMerger::spliceSet(SetMerger &other)
  {
    unsigned int pivot;
    if (membership.size() != other.membership.size())
    {
      cerr
          << "SetMerger::spliceSet(): Cannot merge sets with different numbers of elements!\n";
      return false;
    }
    for (pivot = 0; pivot < membership.size(); pivot++)
      if (min(membership[pivot], other.membership[pivot]) >= 0)
      {
        cerr
            << "SetMerger::spliceSet(): The same element is assigned to different sets!\n";
        return false;
      }

    if (other.sets.size() != other.numSets)
      other.compressSetIndices();

    unsigned int numOldElems = sets.size(), numNewElems = other.sets.size();
    numSets = numOldElems + numNewElems;
    sets.resize(numSets);
    cAlignsASP.resize(numSets);
    cAlignsPA.resize(numSets);
    cAlignsASP_idx.resize(numSets);
    cAlignsPA_idx.resize(numSets);

    for (pivot = 0; pivot < membership.size(); pivot++)
      membership[pivot] = max(membership[pivot], other.membership[pivot]);
    for (unsigned int setIdx = numOldElems, otherIdx = 0; setIdx < numSets; setIdx++, otherIdx++)
    {
      sets[setIdx].swap(other.sets[otherIdx]);

      cAlignsASP[setIdx].resize(other.cAlignsASP[otherIdx].size());
      cAlignsASP_idx[setIdx].resize(other.cAlignsASP[otherIdx].size());
      for (pivot = 0; pivot < other.cAlignsASP[otherIdx].size(); pivot++)
      {
        cAlignsASP[setIdx][pivot] = other.cAlignsASP[otherIdx][pivot];
        cAlignsASP_idx[setIdx][pivot] = other.cAlignsASP_idx[otherIdx][pivot];
      }
      cAlignsPA[setIdx].resize(other.cAlignsPA[otherIdx].size());
      cAlignsPA_idx[setIdx].resize(other.cAlignsPA[otherIdx].size());
      for (pivot = 0; pivot < other.cAlignsPA[otherIdx].size(); pivot++)
      {
        cAlignsPA[setIdx][pivot] = other.cAlignsPA[otherIdx][pivot];
        cAlignsPA_idx[setIdx][pivot] = other.cAlignsPA_idx[otherIdx][pivot];
      }
    }
    other.resize(0);
    return true;
  }

  //---------------------------------------------------------------------------
  void SetMerger::resize(unsigned int numMembers)
  {
    membership.resize(numMembers);
    sets.resize(numMembers);
    freePosList.resize(numMembers);
    for (unsigned int i = 0; i < numMembers; i++)
    {
      membership[i] = -1;
      sets[i].clear();
      if (i < numMembers - 1)
        freePosList[i] = i + 1;
      else
        freePosList[i] = -1;
    }
    if (numMembers > 0)
      firstFreePos = 0;
    else
      firstFreePos = -1;
    numSets = 0;
    cAlignsASP.resize(numSets);
    cAlignsPA.resize(numSets);
    cAlignsASP_idx.resize(numSets);
    cAlignsPA_idx.resize(numSets);
  }

  //---------------------------------------------------------------------------
  unsigned int SetMerger::size()
  {
    return numSets;
  }

  //---------------------------------------------------------------------------
  unsigned int SetMerger::numElemsInSets()
  {
    unsigned int c = 0;
    for (unsigned int i = 0; i < membership.size(); i++)
      if (membership[i] >= 0)
        c++;
    return c;
  }

  //---------------------------------------------------------------------------
  void SetMerger::merge(int setIdx1, int setIdx2)
  {
    if (setIdx1 == setIdx2 or min(setIdx1, setIdx2) < 0)
      return;
    reassignSet(setIdx1, setIdx2);
    numSets--;
  }

  //---------------------------------------------------------------------------
  void SetMerger::splitSet(int setIdx, list<int> &elemsToKeep)
  {
    if (setIdx < 0 or setIdx >= (int)sets.size())
      return;

    // Find the list of elements to remove
    list<int> otherElems;
    list<int>::iterator setIter = sets[setIdx].begin(), keepIter =
        elemsToKeep.begin();
    while (setIter != sets[setIdx].end() and keepIter != elemsToKeep.end())
      if (*setIter == *keepIter)
      {
        setIter++;
        keepIter++;
      }
      else
      {
        if (*setIter < *keepIter)
        {
          otherElems.push_back(*setIter);
          setIter++;
        }
        else
          keepIter++;
      }
    while (setIter != sets[setIdx].end())
    {
      otherElems.push_back(*setIter);
      setIter++;
    }
    removeElements(setIdx, otherElems);

    // Find the subset of alignments between the unused elements
    SpectrumPairSet ncPairsASP(cAlignsASP[setIdx].size());
    vector<unsigned int> ncPairsASP_idx(cAlignsASP[setIdx].size());
    SpectrumPairSet ncPairsPA(cAlignsPA[setIdx].size());
    vector<unsigned int> ncPairsPA_idx(cAlignsPA[setIdx].size());
    unsigned int ncIdx = 0, keepIdx = 0;
    for (unsigned int i = 0; i < cAlignsASP[setIdx].size(); i++)
    {
      if (membership[cAlignsASP[setIdx][i].spec1] >= 0
          and membership[cAlignsASP[setIdx][i].spec2] >= 0)
      {
        cAlignsASP[setIdx][keepIdx] = cAlignsASP[setIdx][i];
        cAlignsASP_idx[setIdx][keepIdx] = cAlignsASP_idx[setIdx][i];
        keepIdx++;
      }

      // Check if the alignment is between unused spectra
      for (setIter = otherElems.begin(); setIter != otherElems.end()
          and *setIter < cAlignsASP[setIdx][i].spec1; setIter++)
        ;
      if (setIter == otherElems.end() or *setIter > cAlignsASP[setIdx][i].spec1)
        continue;
      for (setIter = otherElems.begin(); setIter != otherElems.end()
          and *setIter < cAlignsASP[setIdx][i].spec2; setIter++)
        ;
      if (setIter == otherElems.end() or *setIter > cAlignsASP[setIdx][i].spec2)
        continue;
      ncPairsASP[ncIdx] = cAlignsASP[setIdx][i];
      ncPairsASP_idx[ncIdx++] = cAlignsASP_idx[setIdx][i];
    }
    cAlignsASP[setIdx].resize(keepIdx);
    cAlignsASP_idx[setIdx].resize(keepIdx);
    keepIdx = 0;
    ncPairsASP.resize(ncIdx);
    ncIdx = 0;
    for (unsigned int i = 0; i < cAlignsPA[setIdx].size(); i++)
    {
      if (membership[cAlignsPA[setIdx][i].spec1] >= 0
          and membership[cAlignsPA[setIdx][i].spec2] >= 0)
      {
        cAlignsPA[setIdx][keepIdx] = cAlignsPA[setIdx][i];
        cAlignsPA_idx[setIdx][keepIdx] = cAlignsPA_idx[setIdx][i];
        keepIdx++;
      }

      // Check if the alignment is between unused spectra
      for (setIter = otherElems. begin(); setIter != otherElems.end()
          and *setIter < cAlignsPA[setIdx][i].spec1; setIter++)
        ;
      if (setIter == otherElems.end() or *setIter > cAlignsPA[setIdx][i].spec1)
        continue;
      for (setIter = otherElems.begin(); setIter != otherElems.end()
          and *setIter < cAlignsPA[setIdx][i].spec2; setIter++)
        ;
      if (setIter == otherElems.end() or *setIter > cAlignsPA[setIdx][i].spec2)
        continue;
      ncPairsPA[ncIdx] = cAlignsPA[setIdx][i];
      ncPairsPA_idx[ncIdx++] = cAlignsPA_idx[setIdx][i];
    }
    cAlignsPA[setIdx].resize(keepIdx);
    cAlignsPA_idx[setIdx].resize(keepIdx);
    ncPairsPA.resize(ncIdx);
    if (ncPairsASP.size() + ncPairsPA.size() == 0)
      return;

    // Create new sets from the other eleme nts SetMerger newSets (membership.size());
    SetMerger newSets(membership.size());
    newSets.createSets(membership.size(), 2, ncPairsASP, ncPairsPA);
    newSets.splitAligns(ncPairsASP, ncPairsPA);
    for (unsigned int i = 0; i < newSets.size(); i++)
    {
      for (unsigned int j = 0; j < newSets.cAlignsASP_idx[i].size(); j++)
        newSets.cAlignsASP_idx[i][j]
            = ncPairsASP_idx[newSets.cAlignsASP_idx[i][j]];
      for (unsigned int j = 0; j < newSets.cAlignsPA_idx[i].size(); j++)
        newSets.cAlignsPA_idx[i][j]
            = ncPairsPA_idx[newSets.cAlignsPA_idx[i][j]];
    }

    // Merge new sets with the current sets
    spliceSet( newSets);
  }

  //---------------------------------------------------------------------------
  //
  // removeSet - Removes set with index setIdx. If removeElmtsOnly is true then
  //  the set is replaced with an empty set and its former elements become unassigned;
  //  otherwise the set is also deleted.
  //
  //---------------------------------------------------------------------------
  void SetMerger::removeSet(unsigned int setIdx, bool removeElmtsOnly)
  {
    if (setIdx >= sets.size())
      return;
    for (list<int>::iterator iter = sets[setIdx].begin(); iter
        != sets[setIdx].end(); iter++)
      membership[*iter] = -1;
    sets[setIdx].clear();
    if (not removeElmtsOnly)
    {
      freePos(setIdx);
      numSets--;
    }
  }

  //---------------------------------------------------------------------------
  //
  // removeSmall Sets - removes all sets with less than minSe tSize elements.
  //
  //---------------------------------------------------------------------------
  void SetMerger::removeSmallSets(unsigned int minSetSize)
  {
    for (unsigned int i = 0; i < sets.size(); i++)
      if (sets[i].size() > 0 and sets[i].size() < minSetSize)
      {
        for (list<int>::iterator iter = sets[i].begin(); iter != sets[i].end(); iter++)
          membership[*iter] = -1;
        sets[i].clear();
        freePos(i);
        numSets--;
      }
  }

  //---------------------------------------------------------------------------
  //
  //  compressSetIndices - makes sure all sets are on indices from zero to numSets-1
  //
  //---------------------------------------------------------------------------
  void SetMerger::compressSetIndices()
  {
    if (numSets == 0)
      return;
    unsigned int i;
    numSets = 0;
    list<unsigned int> usedIndices;
    usedIndices.clear();
    for (i = 0; i < sets.size(); i++)
      if (sets[i].size() > 0)
      {
        usedIndices.push_front(i);
        numSets++;
        sets[i].sort();
      }

    // Use higher index sets to fill in gaps between lower index sets
    for (i = 0; usedIndices.size() > 0 and i == usedIndices.front(); i++, usedIndices.pop_front())
      ;
    for (; usedIndices.size() > 0 and i < usedIndices.front(); i++)
    {
      if (sets[i].size() == 0)
      {
        reassignSet(i, usedIndices.front());
        sets[usedIndices.front()].clear();
        usedIndices.pop_front();
      }
    }
    sets.resize(numSets);

    // Reset freePosList
    if (freePosList.size() > (unsigned int)numSets)
    {
      for (i = 0; i < freePosList.size() - numSets - 1; i++)
        freePosList[i] = i + 1;
      freePosList[freePosList.size() - numSets - 1] = -1;
      firstFreePos = 0;
    }
    else
      firstFreePos = -1;
  }

  //---------------------------------------------------------------------------

  int SetMerger::saveas_binListArray(char *filename)
  {
    return Save_binListArray<int, list<int> , list<int>::iterator> (filename,
                                                                    sets);
  }

  //---------------------------------------------------------------------------
  void SetMerger::freePos(int pos)
  {
    freePosList[pos] = firstFreePos;
    firstFreePos = pos;
  }

  //---------------------------------------------------------------------------
  void SetMerger::reassignSet(int setIdx1, int setIdx2)
  {
    for (list<int>::iterator iter = sets[setIdx2].begin(); iter
        != sets[setIdx2].end(); iter++)
      membership[*iter] = setIdx1;
    sets[setIdx1].splice(sets[setIdx1].begin(), sets[setIdx2]);
    freePos(setIdx2);
  }

} // namespace specnets

