/*
 * Contig.cpp
 *
 *  Created on: Jun 23, 2011
 *      Author: aguthals
 */

#include "Contig.h"
#include "Logger.h"
#include "ExecFramework/ExecAssembly.h"
#include "prm_alignment.h"

namespace specnets
{
  Contig::Contig(void) :
    index(0), reversed(false), assembledStars(0x0), innerEdges(0x0),
        rootRef(0x0), childSpectra(0x0), childSpecSet(0x0), Spectrum()
  {
    create();
  }

  Contig::Contig(int idx, SpecSet* contigs, abinfo_t* _assembledStars) :
    index(idx), reversed(false), assembledStars(0x0), innerEdges(0x0),
        rootRef(0x0), childSpectra(0x0), childSpecSet(0x0), Spectrum()
  {
    create();
    initialize(idx, contigs, _assembledStars);
  }

  Contig::Contig(const Contig& other) :
    index(other.index), reversed(false), assembledStars(0x0), innerEdges(0x0),
        rootRef(0x0), childSpectra(0x0), childSpecSet(0x0), Spectrum()
  {
    create();
    copy(other);
  }

  Contig::~Contig(void)
  {
    //DEBUG_TRACE;
    // DEBUG_VAR(assembledStars);
    if (assembledStars)
      delete assembledStars;
    assembledStars = 0;
    //DEBUG_VAR(assembledStars);
    if (innerEdges)
      delete innerEdges;
    innerEdges = 0;
    if (rootRef)
      delete rootRef;
    rootRef = 0;
    if (childSpectra)
      delete childSpectra;
    childSpectra = 0;
    if (childSpecSet)
      delete childSpecSet;
    childSpecSet = 0;
  }

  Contig& Contig::operator=(const Contig &other)
  {
    //DEBUG_TRACE;
    if (this == &other)
    {
      return (*this);
    }
    copy(other);
    return (*this);
  }

  void Contig::initialize(int idx, SpecSet* contigs, abinfo_t* _assembledStars)
  {
    //DEBUG_TRACE;
    if (idx < 0 || idx >= contigs->size())
    {
      ERROR_MSG("Invalid index " << idx
          << " for meta-contig initialization from specset of size "
          << contigs->size());
      return;
    }
    if (_assembledStars != 0 && _assembledStars->count(idx) == 0)
    {
      ERROR_MSG("Invalid index " << idx
          << " not found in abinfo for meta-contig initialization");
      return;
    }
    Spectrum::operator =((*contigs)[idx]);
    scan = idx;

    index = idx;
    reversed = false;
    assembledStars->clear();
    if (_assembledStars != 0)
    {
      (*assembledStars)[0] = (*_assembledStars)[index];
    }
    innerEdges->clear();
    rootRef->clear();
    (*rootRef)[index] = pair<float, float> (0, 0);
    childSpectra->clear();
    (*childSpectra)[index] = pair<int, bool> (0, false);
    childSpecSet->resize(1);
    (*childSpecSet)[0] = (*contigs)[idx];
    endGaps = pair<float, float> (0, 0);
  }

  void Contig::copy(const Contig& other)
  {
    //DEBUG_TRACE;
    Spectrum::operator =((Spectrum&)other);
    index = other.index;
    reversed = other.reversed;
    endGaps = other.endGaps;
    assembledStars->operator =(*other.assembledStars);
    //assembledStars->clear();
    //assembledStars->insert(other.assembledStars->begin(),
    //                       other.assembledStars->end());

    //innerEdges->clear();

    //innerEdges->insert(innerEdges->begin(),
    //                   other.innerEdges->begin(),
    //                   other.innerEdges->end());

    innerEdges->operator =(*other.innerEdges);

    //rootRef->clear();

    //rootRef->insert(other.rootRef->begin(), other.rootRef->end());

    rootRef->operator =(*other.rootRef);

    //childSpectra->clear();

    //childSpectra->insert(other.childSpectra->begin(), other.childSpectra->end());

    childSpectra->operator =(*other.childSpectra);

    childSpecSet->operator =(*other.childSpecSet);
  }

  void Contig::reverse(void)
  {
    Spectrum::reverse(0.0 - AAJumps::massH2O);
    reversed = !reversed;

    float temp = endGaps.first;
    endGaps.first = endGaps.second;
    endGaps.second = temp;

    bool rev1;
    for (map<int, pair<int, bool> >::iterator childIt = childSpectra->begin(); childIt
        != childSpectra->end(); childIt++)
    {
      (*childSpecSet)[childIt->second.first].reverse(0.0 - AAJumps::massH2O);
      rev1 = childIt->second.second;
      childIt->second.second = !rev1;
    }

    for (list<PRMSpecEdge>::iterator innerIt = innerEdges->begin(); innerIt
        != innerEdges->end(); innerIt++)
    {
      innerIt->reverse(true, true);
    }

    for (map<int, pair<float, float> >::iterator refIt = rootRef->begin(); refIt
        != rootRef->end(); refIt++)
    {
      temp = refIt->second.second;
      refIt->second.second = refIt->second.first;
      refIt->second.first = temp;
    }
  }

  unsigned int Contig::assembleConsensus(int minNumMatchedPeaks,
                                         float pkTol,
                                         float pmTol)
  {
    SpecSet inputSpecs(childSpecSet->size());
    inputSpecs = *childSpecSet;

    vector<int> idxRef(inputSpecs.size());
    for (map<int, pair<int, bool> >::iterator childIt = childSpectra->begin(); childIt
        != childSpectra->end(); childIt++)
    {
      idxRef[childIt->second.first] = childIt->first;
      /*
       cout << "\n";
       DEBUG_VAR(childIt->first);
       DEBUG_VAR(childIt->second.first);
       DEBUG_VAR((*rootRef)[childIt->first].first);
       inputSpecs[childIt->second.first].output(cout);
       cout << "\n";
       */
    }

    inputSpecs.setPeakTolerance(pkTol, false);

    SpectrumPairSet inputPairs(innerEdges->size());
    DEBUG_VAR(inputSpecs.size());
    DEBUG_VAR(inputPairs.size());
    int idxUse = 0;
    for (list<PRMSpecEdge>::iterator innerIt = innerEdges->begin(); innerIt
        != innerEdges->end(); innerIt++)
    {
      inputPairs[idxUse] = (SpectrumPair &)(*innerIt);
      inputPairs[idxUse].spec1
          = (*childSpectra)[inputPairs[idxUse].spec1].first;
      inputPairs[idxUse].spec2
          = (*childSpectra)[inputPairs[idxUse].spec2].first;
      inputPairs[idxUse].shift2
          = inputSpecs[inputPairs[idxUse].spec1].parentMass * 1000.0;
      /*
       DEBUG_MSG("inputPairs[idxUse].spec1 = " << inputPairs[idxUse].spec1 << " (" << idxRef[inputPairs[idxUse].spec1] << ")");
       DEBUG_MSG("inputPairs[idxUse].spec2 = " << inputPairs[idxUse].spec2 << " (" << idxRef[inputPairs[idxUse].spec2] << ")");
       DEBUG_VAR(inputPairs[idxUse].shift1);
       DEBUG_VAR(inputPairs[idxUse].shift2);
       */
      ++idxUse;
    }

    if (inputSpecs.size() - 1 != inputPairs.size())
    {
      abort();
    }

    ParameterList inputParams;
    inputParams.addIfDoesntExist("PENALTY_PTM", "-200");
    inputParams.addIfDoesntExist("PENALTY_SAME_VERTEX", "-1000000");
    inputParams.addIfDoesntExist("GRAPH_TYPE", "2");
    inputParams.addIfDoesntExist("MAX_AA_JUMP", "2");
    inputParams.addIfDoesntExist("MAX_MOD_MASS", "100.0");
    inputParams.addIfDoesntExist("TOLERANCE_PEAK", parseFloat(pkTol, 5));
    inputParams.addIfDoesntExist("TOLERANCE_PM", parseFloat(pmTol, 5));
    inputParams.addIfDoesntExist("MIN_MATCHED_PEAKS",
                                 parseInt(minNumMatchedPeaks).c_str());
    inputParams.addIfDoesntExist("MIN_EDGES_TO_COMPONENT", "0");
    inputParams.addIfDoesntExist("PATH_MIN_SPECS", "2");
    inputParams.addIfDoesntExist("PATH_MIN_PEAKS",
                                 parseInt(minNumMatchedPeaks).c_str());
    inputParams.addIfDoesntExist("SPEC_TYPE_MSMS", "0");
    inputParams.addIfDoesntExist("NO_SEQUENCING", "0");
    inputParams.addIfDoesntExist("ADD_ENDPOINTS", "0");
    inputParams.addIfDoesntExist("OUTPUT_COMPLETE_ABRUIJN", "");
    inputParams.addIfDoesntExist("EDGE_SCORE_TYPE", "1");
    inputParams.addIfDoesntExist("IGNORE_REVERSALS", "1");

    //DEBUG_VAR(minNumMatchedPeaks);

    Clusters outputClusters;
    abinfo_t outputAbinfo;

    ExecAssembly assemblyObj(inputParams,
                             &inputSpecs,
                             &inputPairs,
                             &outputClusters,
                             &outputAbinfo);

    //Logger& currentLogger = Logger::getDefaultLogger();
    //Logger::setDefaultLogger(Logger::getLogger(1));

    //Logger::setDefaultLogger(Logger::getLogger(2));

    assemblyObj.invoke();

    //Logger::setDefaultLogger(currentLogger);

    unsigned int numComponents = outputClusters.consensus.size();
    if (outputClusters.consensus.size() == 0
        || outputClusters.consensus[0].size() == 0)
    {
      WARN_MSG("ExecAssembly RETURNED 0 COMPONENTS, MERGING FOR COMPONENT " << index
          << " FAILED!\n");
      //abort();
      return 0;
    }
    unsigned int myScan = scan;
    Spectrum::operator =(outputClusters.consensus[0]);
    scan = myScan;
    assembledStars->clear();
    (*assembledStars)[0] = outputAbinfo[0];

    for (unsigned int i = 0; i < (*assembledStars)[0].first.first.size(); i++)
    {
      int locIdx = (*assembledStars)[0].first.first[i];
      (*assembledStars)[0].first.first[i] = idxRef[locIdx];

      if ((*assembledStars)[0].first.second[i] == 1)
      {
        abort();
      }
    }
    for (unsigned int i = 0; i < (*assembledStars)[0].second.size(); i++)
    {
      for (unsigned int j = 0; j < (*assembledStars)[0].second[i].first.size(); j++)
      {
        int locIdx = (*assembledStars)[0].second[i].first[j];
        (*assembledStars)[0].second[i].first[j] = idxRef[locIdx];
      }
    }

    parentMass = peakList.back()[0] + AAJumps::massMH;
    parentCharge = 1;
    parentMZ = parentMass;

    setPeakTolerance(pkTol);
    setParentMassTol(pmTol);

    float minLeftEdgeF = 0;
    float maxRightEdge = peakList.back()[0];
    for (map<int, pair<int, bool> >::iterator childIt = childSpectra->begin(); childIt
        != childSpectra->end(); childIt++)
    {

      float FFShift = (childIt->first == index) ? 0
          : (*rootRef)[childIt->first].first;
      minLeftEdgeF = min(minLeftEdgeF, FFShift);
      float rightEdge = FFShift
          + (*childSpecSet)[childIt->second.first].back()->operator [](0);
      maxRightEdge = max(maxRightEdge, rightEdge);
      //if (debug) {
      //	cout << "\nContig " << *nodeIt << "(shift=" << FFShift << "):\n";
      //	oriented[*nodeIt].output(cout);
      //}
    }

    maxRightEdge -= minLeftEdgeF;

    if ((*assembledStars)[0].second.size() < 3)
    {
      ERROR_MSG("Resulting ab components are too small!");
      abort();
    }

    if ((*assembledStars)[0].second[1].first.size() == 0)
    {
      ERROR_MSG("No assembled peaks for second abruijn vertex!");
      abort();
    }

    pair<vector<int> , vector<double> >* first_mass =
        &(*assembledStars)[0].second[1];

    int locIdxF = first_mass->first.front();

    endGaps.first = first_mass->second.front() + (*rootRef)[locIdxF].first
        - minLeftEdgeF;
    endGaps.second = maxRightEdge - back()->operator [](0) - endGaps.first;

    /*
     if (debug) {
     cout << "\ngapF = " << gapF << " = " << first_mass.second.front()
     << " + " << contigIndent[locIdxF][0] << "\n";
     cout << "gapR = " << gapR << " = " << maxRightEdge << " - "
     << putSpec.peakList.back()[0] << " - " << gapF << "\n";
     cout << "First mass from " << locIdxF << "\n";
     cout << "minLeftEdgeF " << minLeftEdgeF << "\n";
     }
     */

    return numComponents;
  }

  bool Contig::merge(PRMSpecEdge* inEdge,
                     Contig* inOther,
                     int minNumMatchedPeaks,
                     float pktol,
                     float pmTol)
  {

    //DEBUG_VAR(minNumMatchedPeaks);
    Contig* other = inOther;
    Contig tempOther;

    PRMSpecEdge* edge = inEdge;

    Contig myCopy(*this);

    if (edge->spec2rev)
    {
      ERROR_MSG("Cannot merge with a reversed edge!");
      abort();
      return false;
    }

    for (map<int, pair<float, float> >::iterator refIt =
        other->rootRef->begin(); refIt != other->rootRef->end(); refIt++)
    {
      (*myCopy.rootRef)[refIt->first]
          = pair<float, float> (edge->getShift(myCopy.index, other->index)
              + refIt->second.first, edge->getReversedShift(myCopy.index,
                                                            other->index)
              + refIt->second.second);
    }

    PRMSpecEdge bestConn;
    int numMP = -1;
    pair<int, pair<float, float> > shiftScore;
    PRMAlignment alignmentObj;
    alignmentObj.spec2rev = false;

    for (map<int, pair<int, bool> >::iterator childIt =
        myCopy.childSpectra->begin(); childIt != myCopy.childSpectra->end(); childIt++)
    {

      alignmentObj.spec1 = childIt->first;
      alignmentObj.setSpec1(&(*myCopy.childSpecSet)[childIt->second.first]);

      for (map<int, pair<int, bool> >::iterator otherIt =
          other->childSpectra->begin(); otherIt != other->childSpectra->end(); otherIt++)
      {

        alignmentObj.spec2 = otherIt->first;
        alignmentObj.setSpec2(&(*other->childSpecSet)[otherIt->second.first]);
        alignmentObj.shift1 = (*myCopy.rootRef)[otherIt->first].first
            - (*myCopy.rootRef)[childIt->first].first;
        alignmentObj.shift2 = (*myCopy.rootRef)[otherIt->first].second
            - (*myCopy.rootRef)[childIt->first].second;

        shiftScore = alignmentObj.getShiftScore(alignmentObj.shift1, pktol, 0);
        /*
         DEBUG_VAR(pktol);
         DEBUG_VAR(alignmentObj.spec1);
         DEBUG_VAR(alignmentObj.spec2);
         DEBUG_VAR((*myCopy.rootRef)[otherIt->first].first);
         DEBUG_VAR((*myCopy.rootRef)[childIt->first].first);
         DEBUG_VAR(alignmentObj.shift1);
         DEBUG_VAR(shiftScore.first);
         */
        if (shiftScore.first > numMP)
        {
          numMP = shiftScore.first;
          bestConn = (SpectrumPair&)alignmentObj;
        }
      }
    }

    if (numMP < 0)
    {
      ERROR_MSG("Could not find overlapping contigs when merging meta-contigs " << index << " and " << other->index);
    }

    //DEBUG_MSG("Best node pair: " << bestConn.spec1 << " and " << bestConn.spec2 << "(" << numMP << " matching peaks)");

    myCopy.innerEdges->push_back(bestConn);
    myCopy.innerEdges->insert(myCopy.innerEdges->end(),
                              other->innerEdges->begin(),
                              other->innerEdges->end());

    int curSize = myCopy.childSpecSet->size();

    myCopy.childSpecSet->resize(myCopy.childSpecSet->size()
        + other->childSpectra->size());
    for (map<int, pair<int, bool> >::iterator otherIt =
        other->childSpectra->begin(); otherIt != other->childSpectra->end(); otherIt++)
    {
      int newIdx = otherIt->second.first + curSize;
      (*myCopy.childSpectra)[otherIt->first]
          = pair<int, bool> (newIdx, otherIt->second.second);
      (*myCopy.childSpecSet)[newIdx]
          = (*other->childSpecSet)[otherIt->second.first];
    }

    if (myCopy.assembleConsensus(minNumMatchedPeaks, pktol, pmTol) == 1)
    {
      copy(myCopy);
      return true;
    }
    else
    {
      return false;
    }
  }

  void Contig::create(void)
  {
    if (!assembledStars)
      assembledStars = new abinfo_t;
    if (!innerEdges)
      innerEdges = new list<PRMSpecEdge> ;
    if (!rootRef)
      rootRef = new map<int, pair<float, float> > ;
    if (!childSpectra)
      childSpectra = new map<int, pair<int, bool> > ;
    if (!childSpecSet)
      childSpecSet = new SpecSet;

    assembledStars->clear();
    innerEdges->clear();
    rootRef->clear();
    childSpectra->clear();
    childSpecSet->resize(0);
    endGaps = pair<float, float> (0, 0);
    index = 0;
    reversed = false;
  }
}
