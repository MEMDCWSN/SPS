#include "abruijn.h"

using namespace std;
using namespace specnets;

namespace specnets
{

  // Helper structure used in VertexSet::consolidatePaths
  // This structure is used to ensure that paths are converted to vertices in
  // reverse order - remove last vertex in the path and iterate back to first vertex
  struct vscp_index
  {
    int specIdx, peakIdx, vertexIdx;
  };
  bool operator<(const vscp_index &a, const vscp_index &b)
  {
    return a.specIdx < b.specIdx or (a.specIdx == b.specIdx and a.peakIdx
        < b.peakIdx);
  }

  Vertex::Vertex(const Vertex &other)
  { // copy constructor
    reset();
    specPeaks.assign(other.specPeaks.begin(), other.specPeaks.end());
    outMatchEdges.assign(other.outMatchEdges.begin(), other.outMatchEdges.end());
    outEdges.assign(other.outEdges.begin(), other.outEdges.end());
    compositeVertex = other.compositeVertex;
  }

  //
  // addPeak(specIdx, peakIdx) - adds a new peak to the vertex
  //
  void Vertex::addPeak(TwoValues<int> specPeak)
  {
    int lastSpecIdx = -1;
    list<TwoValues<int> >::iterator pivot = specPeaks.begin();
    while (pivot != specPeaks.end() and (*pivot) < specPeak)
    {
      lastSpecIdx = (*pivot)[0];
      pivot++;
    }
    if (pivot == specPeaks.end() or (*pivot)[0] != specPeak[0] or (*pivot)[1]
        != specPeak[1])
    { // Avoid repeated insertions of the same peak
      specPeaks.insert(pivot, specPeak);
      if (lastSpecIdx == specPeak[0] or (pivot != specPeaks.end()
          and (*pivot)[0] == specPeak[0]))
        compositeVertex = true;
    }
  }

  //
  //  addMatchEdge - adds a new 'outgoing' match edge (from the current vertex)
  //
  void Vertex::addMatchEdge(float edgeMass, Match from, Match to)
  {
    outMatchEdges.push_front(MatchEdge(edgeMass, from, to));
  }

  //
  //  addEdge - adds a new 'outgoing' simple edge (from the current vertex)
  //
  void Vertex::addEdge(float edgeMass,
                       float edgeScore,
                       int specIdx,
                       int peakIdx1,
                       int peakIdx2)
  {
    outEdges.push_front(SimpleEdge(edgeMass,
                                   edgeScore,
                                   specIdx,
                                   peakIdx1,
                                   peakIdx2));
  }

  //
  //  merge - merge this vertex with withVertex
  //
  void Vertex::merge(Vertex &withVertex)
  {
    specPeaks.merge(withVertex.specPeaks);
    outMatchEdges.splice(outMatchEdges.begin(), withVertex.outMatchEdges);
    outEdges.splice(outEdges.begin(), withVertex.outEdges);
    compositeVertex = compositeVertex or withVertex.compositeVertex;
    withVertex.reset();
  }

  //
  //  reset - Resets a vertex to empty (unsused) status
  //
  void Vertex::reset()
  {
    specPeaks.clear();
    outMatchEdges.clear();
    outEdges.clear();
    compositeVertex = false;
  }

  //
  //  replaceEdge - Helper function for consolidatePaths: replaces every edge to toSpecIdx/toPeakIdx
  //    with an edge to newPeakIdx and increments the edge mass by addEdgeMass
  //
  bool Vertex::replaceEdge(int toSpecIdx,
                           int toPeakIdx,
                           int newPeakIdx,
                           float addEdgeMass,
                           float addEdgeScore)
  {
    bool found = false;
    list<SimpleEdge>::iterator iter = outEdges.begin();
    for (; iter != outEdges.end(); iter++)
      if ((*iter).specIdx == toSpecIdx and (*iter).peaks[1] == toPeakIdx)
      {
        iter->peaks[1] = newPeakIdx;
        iter->edgeMass += addEdgeMass;
        iter->edgeScore += addEdgeScore;
        found = true;
      }
    return found;
  }

  VertexSet::VertexSet(const VertexSet &other)
  { // copy constructor
    freeVertices.resize(other.freeVertices.size());
    freeVertices.assign(other.freeVertices.begin(), other.freeVertices.end());
    firstFreeVertex = other.firstFreeVertex;
    peakToVertex.resize(other.peakToVertex.size());
    for (unsigned int i = 0; i < peakToVertex.size(); i++)
    {
      peakToVertex[i].resize(other.peakToVertex[i].size());
      peakToVertex[i].assign(other.peakToVertex[i].begin(),
                             other.peakToVertex[i].end());
    }
    vertices.resize(other.vertices.size());
    vertices.assign(other.vertices.begin(), other.vertices.end());
    numUsedVertices = other.numUsedVertices;
    specSet = other.specSet;
  }

  int VertexSet::mergeVertices(Match m)
  {
    int v1 = peakToVertex[m.specPeak1[0]][m.specPeak1[1]], v2 =
        peakToVertex[m.specPeak2[0]][m.specPeak2[1]], maxV = max(v1, v2), minV =
        min(v1, v2), v;

    //cerr<<" **** m=("<<m.specPeak1[0]<<","<<m.specPeak1[1]<<")/("<<m.specPeak2[0]<<","<<m.specPeak2[1]<<"), (v1,v2) = ("<<v1<<","<<v2<<"), maxV = "<<maxV; cerr.flush();
    //if(maxV>=0) cerr <<" (vertices[maxV].size() = "<<vertices[maxV].size()<<")";
    //cerr<<", minV = "<<minV<<"\n"; cerr.flush();

    if (maxV == minV and maxV > -1)
    {
      //cerr << "Vertex is already the same\n";
      return maxV;
    }
    if (maxV == -1)
    {
      //cerr << "New vertex created\n";
      v = getFreeVertex();
      peakToVertex[m.specPeak1[0]][m.specPeak1[1]] = v;
      vertices[v].addPeak(m.specPeak1);
      peakToVertex[m.specPeak2[0]][m.specPeak2[1]] = v;
      vertices[v].addPeak(m.specPeak2);
      numUsedVertices++;
      return v;
    }
    else
    {
      if (minV > -1)
      { // merge the two vertices in vertices[minV]
        //cerr << "Merging two separate vertices\n";
        list<TwoValues<int> >::iterator pivot1 =
            vertices[minV].specPeaks.begin(), pivot2 =
            vertices[maxV].specPeaks.begin();
        while (pivot1 != vertices[minV].specPeaks.end() and pivot2
            != vertices[maxV].specPeaks.end())
        {
          if ((*pivot1)[0] == (*pivot2)[0])
          {
            //					cout << "Composite vertex created by merging ["<<(*pivot1)[0]<<","<<(*pivot1)[1]<<"],["<<(*pivot2)[0]<<","<<(*pivot2)[1]<<"]. Vertices "<<minV<<" and "<<maxV<<" merged into "<<minV<<endl;
            //					cout << " --- mergeVertices called for ["<<m.specPeak1[0]<<","<<m.specPeak1[1]<<"] and ["<<m.specPeak2[0]<<","<<m.specPeak2[1]<<"]\n"; cout.flush();
            vertices[minV].compositeVertex = true;
            break;
          }
          if ((*pivot1) < (*pivot2))
            pivot1++;
          else
            pivot2++;
        }
        // Go through all peaks in vertices[maxV] and reset their vertex idx to minV
        list<TwoValues<int> >::iterator pivot =
            vertices[maxV].specPeaks.begin();

        while (pivot != vertices[maxV].specPeaks.end())
        {
          peakToVertex[(*pivot)[0]][(*pivot)[1]] = minV;
          pivot++;
        }

        vertices[minV].merge(vertices[maxV]);
        vertices[maxV].reset();
        addFreeVertex(maxV);
        numUsedVertices--;
        return minV;
      }
      else
      { // Merging minV into maxV
        //cerr << "Importing new peak into existing vertex\n";
        if (v1 == -1)
        {
          peakToVertex[m.specPeak1[0]][m.specPeak1[1]] = maxV;
          vertices[maxV].addPeak(m.specPeak1);
        }
        else
        {
          peakToVertex[m.specPeak2[0]][m.specPeak2[1]] = maxV;
          vertices[maxV].addPeak(m.specPeak2);
        }
        return maxV;
      }
    }
  }

  inline int VertexSet::addVertex(int specIdx, int peakIdx)
  {
    int v = peakToVertex[specIdx][peakIdx];
    if (v == -1)
    {
      v = getFreeVertex();
      peakToVertex[specIdx][peakIdx] = v;
      vertices[v].addPeak(TwoValues<int> (specIdx, peakIdx));
      numUsedVertices++;
    }
    return v;
  }

  void VertexSet::addGlues(int specIdx1,
                           int specIdx2,
                           vector<TwoValues<int> > &matches,
                           vector<MSGraph> *spectrumGraphs,
                           bool addMissingEdges)
  {
    float edgeMass1, edgeMass2;
    Match from, to;
    unsigned int i = 0;
    int fromV, toV;

    while (i < matches.size() and (matches[i][0] == -1 or matches[i][1] == -1))
    {
      cerr << "WARNING: missing a match in VertexSet::addGlues: (" << specIdx1
          << "," << matches[i][0] << ") to (" << specIdx2 << ","
          << matches[i][1] << ")\n";
      i++;
    }
    if (i < matches.size())
    {
      if (matches[i][0] >= ((*specSet)[specIdx1]).size() or matches[i][1]
          >= ((*specSet)[specIdx2]).size())
      {
        cerr << "ERROR: (VertexSet::addGlues) Gluing peaks (" << specIdx1
            << "," << matches[i][0] << ")/(" << specIdx2 << ","
            << matches[i][1] << ") but spectra (" << specIdx1 << ","
            << specIdx2 << ") only have (" << ((*specSet)[specIdx1]).size()
            << "," << ((*specSet)[specIdx2]).size() << ") peaks!\n";
        return;
      }
      from.set(specIdx1, matches[i][0], specIdx2, matches[i][1]);
      fromV = mergeVertices(from);
      //cerr << "First vertex index: "<<fromV<<" ["<<from.specPeak1[0]<<","<<from.specPeak1[1]<<","<<from.specPeak2[0]<<","<<from.specPeak2[1]<<"], size = "<<vertices[fromV].size()<<endl;
      i++;
      while (i < matches.size())
      {
        while (i < matches.size() and (matches[i][0] == -1 or matches[i][1]
            == -1))
        {
          cerr << "WARNING: missing a match in VertexSet::addGlues: ("
              << specIdx1 << "," << matches[i][0] << ") to (" << specIdx2
              << "," << matches[i][1] << ")\n";
          i++;
        }
        if (i == matches.size())
          break;
        if (matches[i][0] >= ((*specSet)[specIdx1]).size() or matches[i][1]
            >= ((*specSet)[specIdx2]).size())
        {
          cerr << "ERROR: (VertexSet::addGlues) Gluing peaks (" << specIdx1
              << "," << matches[i][0] << ")/(" << specIdx2 << ","
              << matches[i][1] << ") but spectra (" << specIdx1 << ","
              << specIdx2 << ") only have (" << ((*specSet)[specIdx1]).size()
              << "," << ((*specSet)[specIdx2]).size() << ") peaks!\n";
          return;
        }
        to.set(specIdx1, matches[i][0], specIdx2, matches[i][1]);
        toV = mergeVertices(to);
        //cerr << "Next vertex index: "<<toV<<" ["<<to.specPeak1[0]<<","<<to.specPeak1[1]<<","<<to.specPeak2[0]<<","<<to.specPeak2[1]<<"], size = "<<vertices[toV].size()<<endl;

        // add edges
        /*		edgeMass1 = (*specSet)[to.specPeak1[0]][to.specPeak1[1]][0]-(*specSet)[from.specPeak1[0]][from.specPeak1[1]][0];
         edgeMass2 = (*specSet)[to.specPeak2[0]][to.specPeak2[1]][0]-(*specSet)[from.specPeak2[0]][from.specPeak2[1]][0];
         vertices[fromV].addMatchEdge(edgeMass1, from, to);
         vertices[fromV].addMatchEdge(edgeMass2, from, to);
         */

        if (spectrumGraphs != 0)
        {
          if (addMissingEdges)
            cerr << "addMissingEdges not implemented yet!\n";

        }
        else
        {
          edgeMass1 = (*specSet)[to.specPeak1[0]][to.specPeak1[1]][0]
              - (*specSet)[from.specPeak1[0]][from.specPeak1[1]][0];
          edgeMass2 = (*specSet)[to.specPeak2[0]][to.specPeak2[1]][0]
              - (*specSet)[from.specPeak2[0]][from.specPeak2[1]][0];
          vertices[fromV].addEdge(edgeMass1,
                                  0,
                                  from.specPeak1[0],
                                  from.specPeak1[1],
                                  to.specPeak1[1]);
          vertices[fromV].addEdge(edgeMass2,
                                  0,
                                  from.specPeak2[0],
                                  from.specPeak2[1],
                                  to.specPeak2[1]);
        }

        from = to;
        fromV = toV;
        i++;
      }
    }
  }

  void VertexSet::addEdges(vector<MSGraph> &spectrumGraphs,
                           vector<bool> *usedSpectra)
  {
    int fromV, toV;
    for (unsigned int s = 0; s < spectrumGraphs.size(); s++)
    { // over all spectrum graphs
      if (usedSpectra and !(*usedSpectra)[s])
        continue;
      for (unsigned int e = 0; e < spectrumGraphs[s].edges.size(); e++)
      { // over all edges
        fromV = addVertex(s, spectrumGraphs[s].edges[e][0]);
        toV = addVertex(s, spectrumGraphs[s].edges[e][1]);
        vertices[fromV].addEdge(spectrumGraphs[s].eMasses[e],
                                spectrumGraphs[s].eScores[e],
                                s,
                                spectrumGraphs[s].edges[e][0],
                                spectrumGraphs[s].edges[e][1]);
      }
    }
  }

  //
  //  addEndpointEdges - Adds edges between endpoints of paired spectra to account for
  //    the placement of the modification.
  //
  void VertexSet::addEndpointEdges(SpectrumPairSet &aligns,
                                   vector<vector<TwoValues<int> > > &matches,
                                   vector<float> &modPos,
                                   AAJumps &jumps,
                                   float maxModMass,
                                   float peakTol,
                                   float pmTol)
  {
    unsigned int lowPM, highPM, fromV, toV;
    bool isModPair;
    unsigned int specIdxL, specIdxR, peakIdxL, peakIdxR;
    float jumpMassLeft, jumpMassRight;
    for (unsigned int pIdx = 0; pIdx < aligns.size(); pIdx++)
    {
      isModPair = false;
      jumpMassLeft = 0;
      jumpMassRight = 0;
      specIdxL = aligns[pIdx].spec1; // First assume that spec1 is not aligned to the right of spec2
      specIdxR = aligns[pIdx].spec2;
      //cerr<<"pair "<<pIdx<<" ["<<specIdxL<<","<<specIdxR<<"]"; cerr.flush();
      if ((*specSet)[specIdxL].size() < 4 or (*specSet)[specIdxR].size() < 4)
      {
        //cerr<<" -> spec sizes "<<(*specSet)[specIdxL].size()<<" / "<<(*specSet)[specIdxR].size()<<"\n"; cerr.flush();
        continue; // Spectra must have at least b0/y0/b(n-1)/y(n-1)
      }

      // modification pair?
      if ((fabs(aligns[pIdx].shift1) <= pmTol and fabs(aligns[pIdx].shift2)
          <= maxModMass and fabs(aligns[pIdx].shift2) > pmTol)
          or (fabs(aligns[pIdx].shift1) <= maxModMass
              and fabs(aligns[pIdx].shift2) <= pmTol
              and fabs(aligns[pIdx].shift1) > pmTol))
      {
        jumpMassLeft
            = max(fabs(aligns[pIdx].shift1), fabs(aligns[pIdx].shift2));
        jumpMassRight = jumpMassLeft; // Left/Right decision is made below based on modPos[pIdx]
        //cerr<<", mod="<<jumpMassLeft<<", modPos[pIdx]="<<modPos[pIdx]<<", pm(1)="<<(*specSet)[specIdxL].parentMass<<", pm(2)="<<(*specSet)[specIdxR].parentMass; cerr.flush();
        isModPair = true;
      }

      if (isModPair)
      {
        if ((*specSet)[specIdxL].parentMass > (*specSet)[specIdxR].parentMass)
        {
          // modPos assumes alignment from low-mass to high-mass spectrum
          specIdxL = aligns[pIdx].spec2;
          specIdxR = aligns[pIdx].spec1;
        }
      }
      else
      { // Partial overlaps - set both jumpMassLeft and jumpMassRight
        if (matches[pIdx].empty())
        {
          //cerr<<" -> MATCHES: ";
          //for(unsigned int i=0; i<matches[pIdx].size(); i++)
          //	cerr<<"("<<matches[pIdx][i][0]<<","<<matches[pIdx][i][1]<<")";
          //cerr<<endl;
          continue;
        }
        peakIdxL = matches[pIdx][0][0];
        peakIdxR = matches[pIdx][0][1];
        jumpMassLeft = (*specSet)[specIdxR][peakIdxR][0]
            - (*specSet)[specIdxL][peakIdxL][0];
        // Set jumpMassLeft to closest shift
        if (fabs(jumpMassLeft - aligns[pIdx].shift1) < fabs(jumpMassLeft
            - aligns[pIdx].shift2))
        {
          jumpMassLeft = aligns[pIdx].shift1;
          //cerr<<", left(1)="<<jumpMassLeft; cerr.flush();
        }
        else
        {
          jumpMassLeft = aligns[pIdx].shift2;
          //cerr<<", left(2)="<<jumpMassLeft; cerr.flush();
        }

        jumpMassRight = (*specSet)[specIdxR].parentMass + jumpMassLeft
            - (*specSet)[specIdxL].parentMass;
        //cerr<<", right="<<jumpMassRight; cerr.flush();
      }

      if (jumps.isValid(fabs(jumpMassLeft), peakTol))
      {
        // if isModPair then the mod/sequence-extension must be at the start to generate an edge
        if (!isModPair or (isModPair and modPos[pIdx] <= peakTol))
        {
          //cerr<<" +Start-edge"; cerr.flush();
          for (unsigned int p = 0; p < 2; p++)
            if (jumpMassLeft > -0.00001)
            {
              fromV = addVertex(specIdxL, p);
              toV = addVertex(specIdxR, p);
              if (fromV == toV)
                continue; // Don't add self cycles
              vertices[fromV].addEdge(jumpMassLeft, (*specSet)[specIdxL][p][1]
                  + (*specSet)[specIdxR][p][1], specIdxR, p, p);
            }
            else
            {
              fromV = addVertex(specIdxR, p);
              toV = addVertex(specIdxL, p);
              if (fromV == toV)
                continue; // Don't add self cycles
              vertices[fromV].addEdge(-jumpMassLeft, (*specSet)[specIdxL][p][1]
                  + (*specSet)[specIdxR][p][1], specIdxL, p, p);
            }
        }
      }
      //else cerr<<"(jumpMassLeft invalid)"; cerr.flush();

      if (jumps.isValid(fabs(jumpMassRight), pmTol))
      {
        // if isModPair then the mod/sequence-extension must be at the end to generate an edge
        if (!isModPair or (isModPair and modPos[pIdx]
            >= (*specSet)[specIdxL].parentMass - pmTol))
        {
          unsigned int szL = (*specSet)[specIdxL].size(), szR =
              (*specSet)[specIdxR].size(); // In this case the edge is from smaller PM to larger PM
          //cerr<<"pair ["<<specIdxL<<","<<specIdxR<<"], pMasses = ("<<(*specSet)[specIdxL].parentMass<<","<<(*specSet)[specIdxR].parentMass<<"), endpoint edge added at the end: jumpMassRight = "<<jumpMassRight<<endl;
          //cerr<<" +End-edge"; cerr.flush();
          for (unsigned int p = 1; p < 3; p++)
            if (jumpMassRight > -0.00001)
            {
              fromV = addVertex(specIdxL, szL - p);
              toV = addVertex(specIdxR, szR - p);
              if (fromV == toV)
                continue; // Don't add self cycles
              vertices[fromV].addEdge(jumpMassRight,
                                      (*specSet)[specIdxL][szL - p][1]
                                          + (*specSet)[specIdxR][szR - p][1],
                                      specIdxR,
                                      szL - p,
                                      szR - p);
            }
            else
            {
              fromV = addVertex(specIdxR, szR - p);
              toV = addVertex(specIdxL, szL - p);
              if (fromV == toV)
                continue; // Don't add self cycles
              vertices[fromV].addEdge(-jumpMassRight,
                                      (*specSet)[specIdxL][szL - p][1]
                                          + (*specSet)[specIdxR][szR - p][1],
                                      specIdxL,
                                      szR - p,
                                      szL - p);
            }
        }
      }
      //else cerr<<"(jumpMassRight invalid)"; cerr.flush();
      //cerr<<endl;
    } // for
  } // VertexSet::addEndpointEdges

  //
  //  removeEndpoints - Removes endpoints of b/y type (typeB=true/false) from the
  //     graph, including associated edges.
  //
  void VertexSet::removeEndpoints(bool typeB, float peakTol)
  {
    list<TwoValues<int> >::iterator peakIter;
    list<SimpleEdge>::iterator edgeIter;

    vector<int> inDegree(vertices.size());
    for (unsigned int vIdx = 0; vIdx < inDegree.size(); vIdx++)
      inDegree[vIdx] = 0;
    for (unsigned int vIdx = 0; vIdx < vertices.size(); vIdx++)
      for (edgeIter = vertices[vIdx].outEdges.begin(); edgeIter
          != vertices[vIdx].outEdges.end(); edgeIter++)
        inDegree[peakToVertex[edgeIter->specIdx][edgeIter->peaks[1]]]++;

    unsigned int otherV, specIdx, peakIdx;
    float peakMass, parentMass;
    bool erase;
    // Remove all edges to endpoints before removing endpoint vertices
    for (unsigned int vIdx = 0; vIdx < vertices.size(); vIdx++)
    {
      edgeIter = vertices[vIdx].outEdges.begin();
      while (edgeIter != vertices[vIdx].outEdges.end())
      {
        specIdx = edgeIter->specIdx;
        peakIdx = edgeIter->peaks[1];
        erase = false;
        peakMass = (*specSet)[specIdx][peakIdx][0], parentMass
            = (*specSet)[specIdx].parentMass;
        if (typeB)
        {
          if (peakMass <= peakTol or abs(peakMass + AAJumps::massMH
              - parentMass) <= peakTol)
            erase = true;
        }
        else
        {
          if (abs(peakMass - AAJumps::massH2O) <= peakTol or abs(peakMass
              + AAJumps::massHion - parentMass) <= peakTol)
            erase = true;
        }
        if (not erase)
          edgeIter++;
        else
        {
          otherV = peakToVertex[specIdx][peakIdx];
          if (--inDegree[otherV] == 0)
          {
            peakToVertex[specIdx][peakIdx] = -1;
            releaseVertex(otherV);
            if (vIdx == otherV)
              break;
          }
          edgeIter = vertices[vIdx].outEdges.erase(edgeIter);
        }
      }
    }

    // Remove endpoint peaks
    for (unsigned int vIdx = 0; vIdx < vertices.size(); vIdx++)
    {
      peakIter = vertices[vIdx].specPeaks.begin();
      while (peakIter != vertices[vIdx].specPeaks.end())
      {
        specIdx = (*peakIter)[0];
        peakIdx = (*peakIter)[1];
        erase = false;
        peakMass = (*specSet)[specIdx][peakIdx][0], parentMass
            = (*specSet)[specIdx].parentMass;
        if (typeB)
        {
          if (peakMass <= peakTol or abs(peakMass + AAJumps::massMH
              - parentMass) <= peakTol)
            erase = true;
        }
        else
        {
          if (abs(peakMass - AAJumps::massH2O) <= peakTol or abs(peakMass
              + AAJumps::massHion - parentMass) <= peakTol)
            erase = true;
        }
        if (not erase)
          peakIter++;
        else
        {
          peakToVertex[specIdx][peakIdx] = -1;
          peakIter = vertices[vIdx].specPeaks.erase(peakIter);
        }
      }
      if (vertices[vIdx].specPeaks.size() == 0)
        releaseVertex(vIdx);
    }

  }

  //
  //  scv_addEdges  - splitComposite auxiliary function; adds edges to the pred/succs/intn/comp
  //    lists of adjacent edges
  //
  inline void VertexSet::scv_addEdges(list<SimpleEdge> &edges, vector<
      TwoValues<bool> > &vertexAdjChanged)
  {
    unsigned int fromV, toV;
    list<SimpleEdge>::iterator iter = edges.begin();
    for (; iter != edges.end(); iter++)
    {
      fromV = addVertex(iter->specIdx, iter->peaks[0]);
      toV = addVertex(iter->specIdx, iter->peaks[1]);
      if (toV >= scv_predEdges.size() or fromV >= scv_predEdges.size())
      {
        int prevNumVerts = scv_predEdges.size();
        scv_predEdges.resize(prevNumVerts + 2048);
        scv_succEdges.resize(prevNumVerts + 2048);
        scv_intnEdges.resize(prevNumVerts + 2048);
        scv_compEdges.resize(prevNumVerts + 2048);
        vertexAdjChanged.resize(prevNumVerts + 2048);
        for (unsigned int i = prevNumVerts - 1; i < vertexAdjChanged.size(); i++)
          vertexAdjChanged[i].set(false, false);
      }

      if (!(vertices[fromV].compositeVertex or vertices[toV].compositeVertex))
        continue; // only care about edges involving composite vertices
      //cerr<<"  --- edge ("<<fromV<<","<<toV<<")/("<<iter->specIdx<<","<<iter->edgeMass<<","<<iter->peaks[0]<<","<<iter->peaks[1]<<") has compositeVertex ("<<vertices[fromV].compositeVertex<<","<<vertices[toV].compositeVertex<<")\n";

      if (fromV == toV)
      {
        scv_intnEdges[fromV].push_back(*iter);
        continue;
      }
      else if (vertices[fromV].compositeVertex
          and vertices[toV].compositeVertex)
      {
        scv_compEdges[fromV].push_back(*iter);
        scv_compEdges[toV].push_back(*iter);
        continue; // edges between different composite vertices may become composite/non-composite edges in the future
      }

      if (vertices[fromV].compositeVertex)
      {
        scv_succEdges[fromV].push_back(*iter);
        vertexAdjChanged[fromV][1] = true;
      }
      else
      {
        scv_predEdges[toV].push_back(*iter);
        vertexAdjChanged[toV][0] = true;
      }
    }
    //for(unsigned int i=0;i<scv_compVerts.size();i++)
    //	cerr<<"  --- vertex "<<scv_compVerts[i]<<"/"<<i<<" has incident edge list sizes ("<<scv_predEdges[scv_compVerts[i]].size()<<","<<scv_succEdges[scv_compVerts[i]].size()<<","<<scv_intnEdges[scv_compVerts[i]].size()<<","<<scv_compEdges[scv_compVerts[i]].size()<<")\n";

  }

  //
  //  scv_findSplitEdge  - splitComposite auxiliary function; finds the best split edge
  //    for a given composite vertex cvIdx and a pred/succ direction (dir=0/1)
  //
  inline void VertexSet::scv_findSplitEdge(int cvIdx,
                                           int direction,
                                           float peakTol)
  { // direction 0/1=in/out edges
    int vIdx = scv_compVerts[cvIdx];
    list<SimpleEdge> &edgesList = (direction == 0 ? scv_predEdges[vIdx]
        : scv_succEdges[vIdx]);
    if (edgesList.size() == 0)
    {
      scv_splitEdges[cvIdx][direction].clear();
      return;
    }
    vector<TwoValues<float> > edges(edgesList.size());
    list<SimpleEdge>::iterator iter = edgesList.begin();
    for (unsigned int i = 0; i < edges.size(); i++)
    {
      edges[i].set(peakToVertex[iter->specIdx][iter->peaks[direction]],
                   iter->edgeMass);
      iter++;
    }
    sort(edges.begin(), edges.end());

    //cerr<<"  - scv_findSplitEdge: picking best split edge out of "<<edges.size()<<" spectrum peak edges\n";

    int hmVert = -1, hmMult = -1, // Highest multiplicity vertex index and split edge multiplicity
        curVert, curMult; //   same for current vertex
    unsigned int eIdx = 0;
    float hmMass = 0.0, curMass, // Highest multiplicity and current edge masses
        hmMinMass, hmMaxMass, curMinMass, curMaxMass; // Limit bounds on the masses of the selected edges
    while (eIdx < edges.size())
    {
      curVert = (int)edges[eIdx][0];
      curMass = edges[eIdx][1];
      curMult = 1;
      eIdx++;
      curMinMass = curMass;
      curMaxMass = curMass;
      while (eIdx < edges.size() and curVert == (int)edges[eIdx][0]
          and edges[eIdx][1] - curMass <= peakTol)
      {
        curMass = (curMult * curMass + edges[eIdx][1]) / (curMult + 1);
        curMult++;
        curMaxMass = edges[eIdx][1];
        eIdx++;
      }
      if (curMult > hmMult or (curMult == hmMult and curMass < hmMass))
      {
        hmVert = curVert;
        hmMult = curMult;
        hmMass = curMass;
        hmMinMass = curMinMass;
        hmMaxMass = curMaxMass;
      }
    }

    //cerr<<"  - scv_findSplitEdge: selected adj. vertex "<<hmVert<<", mass = ["<<hmMinMass<<","<<hmMaxMass<<"], mult = "<<hmMult<<"\n";

    // Build list of spectrum graph edges in the split edge
    scv_splitEdges[cvIdx][direction].clear();
    for (iter = edgesList.begin(); iter != edgesList.end(); iter++)
      if (peakToVertex[iter->specIdx][iter->peaks[direction]] == hmVert
          and iter->edgeMass >= hmMinMass and iter->edgeMass <= hmMaxMass)
        scv_splitEdges[cvIdx][direction].push_back(*iter);
    scv_splitEdges[cvIdx][direction].sort();
  }

  //
  //  scv_useSplitEdge  - splitComposite auxiliary function; uses a split edge to break
  //    a composite vertex into two vertices: a new vertex which is guaranteed not to be composite
  //    and a leftover vertex with the remaining glued spectrum peaks that may still be composite.
  //
  void VertexSet::scv_useSplitEdge(int cvIdx, int direction, vector<TwoValues<
      bool> > &vertexAdjChanged)
  { // cvIdx is index in scv_compVerts
    int vertIdx = scv_compVerts[cvIdx];
    list<SimpleEdge> &edges = scv_splitEdges[cvIdx][direction];

    // Separate vertex into two
    list<TwoValues<int> >::iterator peaksIter =
        vertices[vertIdx].specPeaks.begin();
    list<SimpleEdge>::iterator edgesIter = edges.begin();
    list<TwoValues<int> > newSpecPeaks;
    while (peaksIter != vertices[vertIdx].specPeaks.end() and edgesIter
        != edges.end())
      if ((*peaksIter)[0] == edgesIter->specIdx and (*peaksIter)[1]
          == edgesIter->peaks[1 - direction])
      {
        newSpecPeaks.push_back(*peaksIter);
        peaksIter = vertices[vertIdx].specPeaks.erase(peaksIter);
        edgesIter++;
      }
      else if ((*peaksIter)[0] < edgesIter->specIdx or ((*peaksIter)[0]
          == edgesIter->specIdx and (*peaksIter)[1] < edgesIter->peaks[1
          - direction]))
        peaksIter++;
      else
        edgesIter++;

    //cerr<<"  - scv_useSplitEdge: split edge - "; cerr.flush();
    //for(edgesIter = edges.begin();edgesIter!=edges.end();edgesIter++) cerr<<"("<<edgesIter->specIdx<<","<<edgesIter->peaks[0]<<","<<edgesIter->peaks[1]<<")"; cerr<<endl;
    //cerr<<"  - scv_useSplitEdge: new vertex peaks - "; cerr.flush();
    //for(peaksIter=newSpecPeaks.begin();peaksIter!=newSpecPeaks.end();peaksIter++) cerr<<"("<<(*peaksIter)[0]<<","<<(*peaksIter)[1]<<")"; cerr<<endl;
    //cerr<<"  - scv_useSplitEdge: leftover vertex peaks - "; cerr.flush();
    //for(peaksIter=vertices[vertIdx].specPeaks.begin();peaksIter!=vertices[vertIdx].specPeaks.end();peaksIter++) cerr<<"("<<(*peaksIter)[0]<<","<<(*peaksIter)[1]<<")"; cerr<<endl;

    // Clear refs to this vertex from other vertices compEdges list
    list<int> otherVerts;
    list<int>::iterator vIter;
    for (edgesIter = scv_compEdges[vertIdx].begin(); edgesIter
        != scv_compEdges[vertIdx].end(); edgesIter++)
      if (peakToVertex[edgesIter->specIdx][edgesIter->peaks[0]] == vertIdx)
        otherVerts.push_back(peakToVertex[edgesIter->specIdx][edgesIter->peaks[1]]);
      else
        otherVerts.push_back(peakToVertex[edgesIter->specIdx][edgesIter->peaks[0]]);
    otherVerts.sort();
    vIter = otherVerts.begin();
    while (vIter != otherVerts.end())
    {
      edgesIter = scv_compEdges[*vIter].begin();
      while (edgesIter != scv_compEdges[*vIter].end())
        if (peakToVertex[edgesIter->specIdx][edgesIter->peaks[0]] == vertIdx
            or peakToVertex[edgesIter->specIdx][edgesIter->peaks[1]] == vertIdx)
          edgesIter = scv_compEdges[*vIter].erase(edgesIter);
        else
          edgesIter++;
      int lastVert = *vIter;
      vIter++;
      while (vIter != otherVerts.end() and lastVert == *vIter)
        vIter++;
    }

    // Add new vertex to VertexSet
    int newVertIdx = getFreeVertex();
    numUsedVertices++;
    for (peaksIter = newSpecPeaks.begin(); peaksIter != newSpecPeaks.end(); peaksIter++)
    {
      peakToVertex[(*peaksIter)[0]][(*peaksIter)[1]] = newVertIdx;
      vertices[newVertIdx].addPeak(*peaksIter);
    }

    // Determine whether leftover vertices still define a composite vertex
    int stillComposite = false;
    if (vertices[vertIdx].specPeaks.size() >= 2)
    {
      peaksIter = vertices[vertIdx].specPeaks.begin();
      int prevSpecIdx = (*peaksIter)[0];
      peaksIter++;
      while (peaksIter != vertices[vertIdx].specPeaks.end()
          and not stillComposite)
      {
        stillComposite = stillComposite or (prevSpecIdx == (*peaksIter)[0]);
        prevSpecIdx = (*peaksIter)[0];
        peaksIter++;
      }
    }
    vertices[vertIdx].compositeVertex = stillComposite;

    //cerr<<"  - scv_useSplitEdge: leftover vertex stil composite = "<<stillComposite<<endl;

    // Reprocess all edges connecting to the previous composite vertex
    list<SimpleEdge> leftoverEdges;
    if (not stillComposite)
    {
      int lastVert = scv_compVerts.size() - 1;
      if (cvIdx < lastVert)
      {
        scv_compVerts[cvIdx] = scv_compVerts[lastVert];
        scv_splitEdges[cvIdx][0].clear();
        scv_splitEdges[cvIdx][0].splice(scv_splitEdges[cvIdx][0].begin(),
                                        scv_splitEdges[lastVert][0]);
        scv_splitEdges[cvIdx][1].clear();
        scv_splitEdges[cvIdx][1].splice(scv_splitEdges[cvIdx][1].begin(),
                                        scv_splitEdges[lastVert][1]);
      }
      scv_compVerts.resize(lastVert);
      scv_splitEdges.resize(lastVert);
      scv_predEdges[vertIdx].clear();
      scv_succEdges[vertIdx].clear();
      scv_intnEdges[vertIdx].clear();

      leftoverEdges.splice(leftoverEdges.begin(), scv_compEdges[vertIdx]);
    }
    else
    {
      leftoverEdges.splice(leftoverEdges.begin(), scv_predEdges[vertIdx]);
      leftoverEdges.splice(leftoverEdges.begin(), scv_succEdges[vertIdx]);
      leftoverEdges.splice(leftoverEdges.begin(), scv_intnEdges[vertIdx]);
      leftoverEdges.splice(leftoverEdges.begin(), scv_compEdges[vertIdx]);
      vertexAdjChanged[vertIdx][0] = true;
      vertexAdjChanged[vertIdx][1] = true;
    }
    scv_addEdges(leftoverEdges, vertexAdjChanged);
  }

  //
  //  splitComposite - Uses edges between composite and non-composite vertices to
  //    separate composite vertices into two or more non-composite vertices (similar to RepeatGluer)
  //
  void VertexSet::splitComposite(vector<MSGraph> &spectrumGraphs,
                                 float peakTol,
                                 vector<bool> *usedSpectra)
  {
    vector<TwoValues<bool> > vertexAdjChanged(vertices.size());
    for (unsigned int i = 0; i < vertexAdjChanged.size(); i++)
      vertexAdjChanged[i].set(false, false);
    scv_predEdges.resize(vertices.size());
    scv_succEdges.resize(vertices.size());
    scv_intnEdges.resize(vertices.size());
    scv_compEdges.resize(vertices.size());

    // Build scv_compVerts
    list<int> compVerts;
    for (unsigned int i = 0; i < vertices.size(); i++)
      if (vertices[i].compositeVertex)
        compVerts.push_back(i);
    list<int>::iterator iter = compVerts.begin();
    scv_compVerts.resize(compVerts.size());
    scv_splitEdges.resize(compVerts.size());
    for (unsigned int i = 0; i < scv_compVerts.size(); i++)
    {
      scv_compVerts[i] = (*iter);
      iter++;
    }
    compVerts.clear();

    // Build list<SimpleEdge> with all edges in spectrumGraphs + call addEdges
    list<SimpleEdge> allEdges;
    for (unsigned int s = 0; s < spectrumGraphs.size(); s++)
    { // over all spectrum graphs
      if (usedSpectra and !(*usedSpectra)[s])
        continue;
      for (unsigned int e = 0; e < spectrumGraphs[s].edges.size(); e++) // over all edges
        allEdges.push_back(SimpleEdge(spectrumGraphs[s].eMasses[e],
                                      spectrumGraphs[s].eScores[e],
                                      s,
                                      spectrumGraphs[s].edges[e][0],
                                      spectrumGraphs[s].edges[e][1]));
    }
    //cerr <<"  - Adding "<<allEdges.size()<<" edges to pred/succ/intn/comp Edges\n";
    scv_addEdges(allEdges, vertexAdjChanged);

    // Remove all peaks (in composite vertices) that do _not_ have an edge connecting to it

    int vIdx; // Index (in vertices) of the composite vertex being processed
    while (scv_compVerts.size() > 0)
    {
      int hmVert = -1, hmDir = -1, hmMult = 0; // Highest multiplicity vertex index, direction and split edge multiplicity
      for (unsigned int cvIdx = 0; cvIdx < scv_compVerts.size(); cvIdx++)
      {
        vIdx = scv_compVerts[cvIdx];
        for (unsigned int dir = 0; dir < 2; dir++)
        {
          if (vertexAdjChanged[vIdx][dir])
          {
            scv_findSplitEdge(cvIdx, dir, peakTol);
            vertexAdjChanged[vIdx][dir] = false;
          }
          //cerr<<"  - Vertex "<<vIdx<<" (cvIdx="<<cvIdx<<"), dir "<<dir<<" has "<<scv_splitEdges[cvIdx][dir].size()<<" incident edges\n";
          if ((int)scv_splitEdges[cvIdx][dir].size() > hmMult)
          {
            hmVert = cvIdx;
            hmDir = dir;
            hmMult = scv_splitEdges[cvIdx][dir].size();
          }
        }
      }
      if (hmVert < 0)
      {
        cerr << "Insufficient edges to break composite vertices!\n";
        return;
      }
      //cerr<<"  -> Highest multiplicity split edge is incident on ("<<scv_compVerts[hmVert]<<"/"<<hmVert<<"), dir "<<hmDir<<" with multiplicity "<<hmMult<<endl;
      scv_useSplitEdge(hmVert, hmDir, vertexAdjChanged);
    }

    scv_predEdges.resize(0);
    scv_succEdges.resize(0);
    scv_intnEdges.resize(0);
    scv_compEdges.resize(0);
    scv_compVerts.resize(0);
    scv_splitEdges.resize(0);
  }

  //
  //  consolidatePaths() - transforms paths where all edges and vertices have
  //    multiplicity one into simple edges.
  //
  void VertexSet::consolidatePaths()
  {
    vector<bool> keepVertex(vertices.size());
    vector<int> predVertex(vertices.size());
    list<SimpleEdge>::iterator iter;
    list<vscp_index> verticesToRemove;
    list<vscp_index>::reverse_iterator removeIter;
    vscp_index tmp;

    for (unsigned int i = 0; i < vertices.size(); i++)
    {
      keepVertex[i] = false;
      predVertex[i] = -1;
    }

    // Find predecessors and mark vertices to remove
    int succV;
    for (unsigned int v = 0; v < vertices.size(); v++)
    {
      if (vertices[v].specPeaks.size() != 1 or vertices[v].outMatchEdges.size()
          != 0 or vertices[v].outEdges.size() != 1)
        keepVertex[v] = true;
      for (iter = vertices[v].outEdges.begin(); iter
          != vertices[v].outEdges.end(); iter++)
      {
        succV = peakToVertex[iter->specIdx][iter->peaks[1]];
        if (predVertex[succV] >= 0)
          keepVertex[succV] = true;
        else
          predVertex[succV] = v;
      }
    }

    // Decide which vertices to remove and order them
    for (unsigned int v = 0; v < vertices.size(); v++)
    {
      if (!keepVertex[v] and predVertex[v] >= 0)
      { // predVertex[v]<0 when vertex is a source
        iter = vertices[v].outEdges.begin();
        tmp.specIdx = iter->specIdx;
        tmp.peakIdx = iter->peaks[0];
        tmp.vertexIdx = v;
        verticesToRemove.push_front(tmp);
      }
    }
    verticesToRemove.sort();

    // Remove marked vertices
    for (removeIter = verticesToRemove.rbegin(); removeIter
        != verticesToRemove.rend(); removeIter++)
    {
      int v = removeIter->vertexIdx;
      TwoValues<int> thisPeak = *vertices[v].specPeaks.begin();
      SimpleEdge thisEdge = *vertices[v].outEdges.begin();
      //cerr<<"Removing vertex "<<v<<", specPeak = ("<<thisPeak[0]<<","<<thisPeak[1]<<"), edge = ("<<thisEdge.specIdx<<","<<thisEdge.peaks[0]<<","<<thisEdge.peaks[1]<<")\n";
      vertices[predVertex[v]].replaceEdge(thisPeak[0],
                                          thisPeak[1],
                                          thisEdge.peaks[1],
                                          thisEdge.edgeMass,
                                          thisEdge.edgeScore);
      peakToVertex[thisPeak[0]][thisPeak[1]] = -1;
      vertices[v].reset();
      addFreeVertex(v);
    }
  }

  //
  //  buildGraph - converts VertexSet into a graph structure by eliminating all
  //    references to spectrum peaks and merging edges by their amino acid masses
  //    as defined by jumps (edges are merged if they're closest to the same entry
  //    in jumps within peakTol). vSet_index[i] = index of the i-th vertex of g in this VertexSet.
  //
  void VertexSet::buildGraph(MSGraph &g, AAJumps &jumps, float peakTol, vector<
      int> &vSet_index, vector<SpectrumPeakLabels> &labels, short edgeScoreType)
  {
    // Structures for renaming vertices - new graph only contains vertices
    //   that contain spectrum peaks.
    vector<int> newVertNums(vertices.size());
    int curVert = 0, numEdges = 0;
    for (int i = 0; i < (int)vertices.size(); i++)
    {
      if (vertices[i].size() == 0)
        newVertNums[i] = -1;
      else
      {
        newVertNums[i] = curVert++;
        numEdges += vertices[i].outEdges.size();
      }
    }
    //cout << "  - Number of vertices: " << curVert <<", total number of edges = " << numEdges;
    //if(numEdges==0) { cout <<" (WARNING: graph not built!)\n"; cout.flush(); return; }
    g.vNext.resize(curVert);
    g.vScores.resize(curVert);
    g.vPeakLabels.resize(curVert);
    g.edges.resize(numEdges);
    g.eScores.resize(numEdges);
    g.eMasses.resize(numEdges);
    vSet_index.resize(curVert);
    //	vector<list<MatchEdge> >  curDest(curVert);       // Lists of outgoing edges per target vertex
    vector<list<SimpleEdge> > curDest(curVert); // Lists of outgoing edges per target vertex
    vector<TwoValues<float> > curEdges(jumps.size()); // Used to merge outgoing edges per [source,sink] edge pair
    //   Mean edge mass (pos.0), num. edges (pos.1)
    vector<float> curEdgeScores(jumps.size()); // Keeps track of the merged edge scores for edges in curEdges
    list<float> addtnlJumps; // List of edge masses that don't match anything in jumps
    list<TwoValues<float> > addtnlEdges; // Like curEdges but for addtnlJumps
    list<TwoValues<float> > addtnlEdgeScores; // Like curEdgeScores but for addtnlEdges (original edge score(pos.0)/accumulated edge score(pos.1))

    int gVert, gEdge = 0; // Keeps track of the current vertex and edge entries in g
    int vPivot, ePivot; // Used to iterate over vertices and edges
    curVert = 0; // Keeps track of the current vertex in VertexSet
    list<float>::iterator addtnlIter; // Used to iterate addtnlJumps
    list<TwoValues<float> >::iterator addtnlIter2, addtnlIter2sc; // Used to iterate addtnlEdges / addtnlEdgeScores
    //	list<MatchEdge>::iterator edge;
    list<SimpleEdge>::iterator edge;
    list<int> gSuccList; // List of successor vertices (holds intermediate results before setting g.vNext[gVert]
    bool matchedJumps; // Indicates whether a specific edge mass matched an entry in jumps
    while (curVert < (int)vertices.size())
    {
      // Find vertex to process
      while (curVert < (int)vertices.size() and newVertNums[curVert] == -1)
        curVert++;
      if (curVert == (int)vertices.size())
        break;
      gVert = newVertNums[curVert];
      gSuccList.clear();
      vSet_index[gVert] = curVert;

      // Compute score for gVert - add scores of all matched peaks
      g.vScores[gVert] = 0;
      list<TwoValues<int> >::iterator iter =
          vertices[curVert].specPeaks.begin();
      for (; iter != vertices[curVert].specPeaks.end(); iter++)
      {
        g.vScores[gVert] += (*specSet)[(*iter)[0]][(*iter)[1]][1];
        if (!labels.empty())
          g.vPeakLabels[gVert].merge(labels[(*iter)[0]][(*iter)[1]]);
      }

      // Separate edges per destination vertex (instead of per-source as in vertices)
      for (vPivot = 0; vPivot < (int)curDest.size(); vPivot++)
        curDest[vPivot].clear();
      for (edge = vertices[curVert].outEdges.begin(); edge
          != vertices[curVert].outEdges.end(); edge++)
      {
        // if addParallelPaths was called, need to get destination vertex from edge->vertexIdx (otherwise lookup peak index)
        vPivot = (edge->vertexIdx < 0)
            ? newVertNums[peakToVertex[edge->specIdx][edge->peaks[1]]]
            : newVertNums[edge->vertexIdx];
        curDest[vPivot].push_front(*edge);
      }

      // For each destination vertex, merge edges according to jumps
      for (vPivot = 0; vPivot < (int)curDest.size(); vPivot++)
      {
        if (curDest[vPivot].size() == 0)
          continue;

        // Clear curEdges and addtnlEdges
        addtnlEdges.clear();
        addtnlJumps.clear();
        addtnlEdgeScores.clear();
        for (ePivot = 0; ePivot < (int)curEdges.size(); ePivot++)
        {
          curEdges[ePivot].set(0, 0);
          curEdgeScores[ePivot] = 0;
        }
        // Merge edges
        for (edge = curDest[vPivot].begin(); edge != curDest[vPivot].end(); edge++)
        {
          matchedJumps = false;

          // If possible merge into one of the edge masses in jumps
          for (ePivot = 0; ePivot < (int)curEdges.size(); ePivot++)
            if (abs(edge->edgeMass - jumps[ePivot]) <= peakTol)
            {
              curEdges[ePivot][0] = (edge->edgeMass + curEdges[ePivot][0]
                  * curEdges[ePivot][1]) / (curEdges[ePivot][1] + 1);
              curEdges[ePivot][1]++;
              curEdgeScores[ePivot] += ((edgeScoreType > 0) ? edge->edgeScore
                  : 1);
              matchedJumps = true;
            }

          // Add edge score to all edges with non-standard masses (i.e., not in jumps) within tolerance
          addtnlIter = addtnlJumps.begin();
          addtnlIter2 = addtnlEdges.begin();
          addtnlIter2sc = addtnlEdgeScores.begin();
          TwoValues<float> addtnlEntry(edge->edgeMass, 1.0);
          float addtnlEntryScore = ((edgeScoreType > 0) ? edge->edgeScore : 1); // Cumulative version of this additional entry
          for (; addtnlIter != addtnlJumps.end(); addtnlIter++, addtnlIter2++, addtnlIter2sc++)
            if (abs(edge->edgeMass - (*addtnlIter)) <= peakTol)
            {
              // Add this edgeMass to others within tolerance
              (*addtnlIter2)[0] = (edge->edgeMass + (*addtnlIter2)[0]
                  * (*addtnlIter2)[1]) / ((*addtnlIter2)[1] + 1);
              (*addtnlIter2)[1]++;
              (*addtnlIter2sc)[1]
                  += ((edgeScoreType > 0) ? edge->edgeScore : 1);
              // Update this edgeMass according to others within tolerance
              addtnlEntry[0]
                  = ((*addtnlIter) + addtnlEntry[0] * addtnlEntry[1])
                      / (addtnlEntry[1] + 1);
              addtnlEntry[1]++;
              addtnlEntryScore += ((edgeScoreType > 0) ? (*addtnlIter2sc)[0]
                  : 1);
            }

          if (not matchedJumps)
          {
            addtnlJumps.push_front(edge->edgeMass);
            addtnlEdges.push_front(addtnlEntry);
            addtnlEdgeScores.push_front(TwoValues<float> (edge->edgeScore,
                                                          addtnlEntryScore));
          }
        }

        // Remove repeated edges from addtnlEdges
        addtnlEdges.sort();
        addtnlIter2 = addtnlEdges.begin();
        addtnlIter2sc = addtnlEdgeScores.begin();
        list<TwoValues<float> >::iterator nextAddtnl = addtnlIter2,
            nextAddtnlSc = addtnlIter2sc;
        nextAddtnl++;
        nextAddtnlSc++;
        while (addtnlIter2 != addtnlEdges.end() and nextAddtnl
            != addtnlEdges.end())
        {
          while (nextAddtnl != addtnlEdges.end() and abs((*addtnlIter2)[0]
              - (*nextAddtnl)[0]) <= 0.0001)
          {
            (*addtnlIter2)[1] = max((*addtnlIter2)[1], (*nextAddtnl)[1]);
            (*addtnlIter2sc)[1] = max((*addtnlIter2sc)[1], (*nextAddtnlSc)[1]);
            nextAddtnl = addtnlEdges.erase(nextAddtnl);
            nextAddtnlSc = addtnlEdgeScores.erase(nextAddtnlSc);
          }
          addtnlIter2 = nextAddtnl;
          nextAddtnl++; // Move on to next pair of edges
          addtnlIter2sc = nextAddtnlSc;
          nextAddtnlSc++;
        }

        // Add edges to g
        int numEdges = addtnlEdges.size();
        for (ePivot = 0; ePivot < (int)curEdges.size(); ePivot++)
          if (curEdges[ePivot][1] > 0)
            numEdges++;
        //cerr << "  --- number of edges: " << numEdges << endl;   cerr.flush();
        for (ePivot = 0; ePivot < (int)curEdges.size(); ePivot++)
          if (curEdges[ePivot][1] > 0)
          {
            //cerr<< "   --- edge " << curEdges[ePivot][0] << ", " << curEdges[ePivot][1]<< ", ePivot = "<<ePivot<<", edgeScore = "<<curEdgeScores[ePivot]<<endl;
            if (gEdge == g.edges.size())
            {
              g.edges.resize(2 * g.edges.size());
              g.eMasses.resize(g.edges.size());
              g.eScores.resize(g.edges.size());
            }
            g.edges[gEdge].set(gVert, vPivot);
            g.eMasses[gEdge] = curEdges[ePivot][0];
            //cerr<<"EST curEdges["<<ePivot<<"]:"<<edgeScoreType;
            switch (edgeScoreType)
            {
            case EST_ABVERTEX_SCORES:
              g.eScores[gEdge] = g.vScores[vPivot] + curEdges[ePivot][1];
              break;
            default:
              g.eScores[gEdge] = curEdgeScores[ePivot]; // curEdges[ePivot][1];
            }
            //cerr<<", ok\n";
            gSuccList.push_front(gEdge++);
          }

        addtnlIter2 = addtnlEdges.begin();
        addtnlIter2sc = addtnlEdgeScores.begin(); // Add aditional edges that don't match jumps
        for (; addtnlIter2 != addtnlEdges.end(); addtnlIter2++, addtnlIter2sc++)
        {
          if (gEdge == g.edges.size())
          {
            g.edges.resize(2 * g.edges.size());
            g.eMasses.resize(g.edges.size());
            g.eScores.resize(g.edges.size());
          }
          g.edges[gEdge].set(gVert, vPivot);
          //cerr<< "   --- edge " << gVert << ", " << vPivot<< ", addtnl jump mass = "<<(*addtnlIter2)[0]<<", edgeScore = "<<(*addtnlIter2sc)[1]<<endl;
          g.eMasses[gEdge] = (*addtnlIter2)[0];
          //cerr<<"EST curEdges["<<ePivot<<"]:"<<edgeScoreType;
          switch (edgeScoreType)
          {
          case EST_ABVERTEX_SCORES:
            g.eScores[gEdge] = g.vScores[vPivot] + (*addtnlIter2)[1];
            break;
          default:
            g.eScores[gEdge] = (*addtnlIter2sc)[1]; // (*addtnlIter2)[1];
          }
          //cerr<<", ok\n";
          gSuccList.push_front(gEdge++);
        }
        //cerr<<"   --- gEdge = "<<gEdge<<endl;
        //cerr << " done.\n";   cerr.flush();
      } // end of merging edges
      g.vNext[gVert].resize(gSuccList.size());
      g.vNext[gVert].assign(gSuccList.begin(), gSuccList.end());
      curVert++;
    } // end of vertex processing
    g.edges.resize(gEdge);
    g.eMasses.resize(gEdge);
    g.eScores.resize(gEdge);
  }
  /*
   void VertexSet::addParallelPaths(float peakTol)
   {
   abruijn::AbruijnGraph abGraph;

   map<abruijn::AbruijnNode*, int> nodeToVertIdx;
   map<int, map<int, abruijn::AbruijnNode*> > peakToNode;
   map<int, vector<set<int> > > connectedEdges;
   set<unsigned int> assembledSpecs;

   for (unsigned int i = 0; i < vertices.size(); i++)
   {
   if (vertices[i].specPeaks.size() == 0)
   {
   continue;
   }

   abruijn::AbruijnNode newNode;
   for (list<TwoValues<int> >::iterator pIt = vertices[i].specPeaks.begin(); pIt
   != vertices[i].specPeaks.end(); pIt++)
   {
   int specIdx = (*pIt)[0];
   int peakIdx = (*pIt)[1];
   newNode.addPeak((*specSet)[specIdx], peakIdx);
   assembledSpecs.insert(specIdx);
   }
   newNode.setGreen(false);

   abruijn::AbruijnNode* addedNode = abGraph.addAbruijnNode(newNode);
   nodeToVertIdx[addedNode] = i;

   for (list<TwoValues<int> >::iterator pIt = vertices[i].specPeaks.begin(); pIt
   != vertices[i].specPeaks.end(); pIt++)
   {
   int specIdx = (*pIt)[0];
   int peakIdx = (*pIt)[1];

   if (peakToNode.count(specIdx) == 0)
   {
   map<int, abruijn::AbruijnNode*> tempMap;
   tempMap[peakIdx] = addedNode;
   peakToNode[specIdx] = tempMap;
   }
   else
   {
   peakToNode[specIdx][peakIdx] = addedNode;
   }
   }
   }

   map<string, unsigned int> specIdToIdx;
   for (set<unsigned int>::iterator specIt = assembledSpecs.begin(); specIt
   != assembledSpecs.end(); specIt++)
   {
   unsigned int specIdx = *specIt;
   specIdToIdx[(*specSet)[specIdx].getUniqueID()] = specIdx;
   vector<set<int> > tempMap((*specSet)[specIdx].size());

   if (peakToNode.count(specIdx) == 0)
   {
   continue;
   }

   if (abruijn::AbruijnGraph::DEBUG_EXPANSION)
   {
   DEBUG_MSG("Spectrum " << (*specSet)[specIdx].scan << ":");
   (*specSet)[specIdx].output(cerr);
   cerr << "\n";
   }

   for (unsigned int i = 0; i < (*specSet)[specIdx].size(); i++)
   {
   if (peakToNode[specIdx].count(i) == 0)
   {
   continue;
   }

   peakToVertex[specIdx][i] = nodeToVertIdx[peakToNode[specIdx][i]];

   for (unsigned int j = i + 1; j < (*specSet)[specIdx].size(); j++)
   {
   if (peakToNode[specIdx].count(j) == 0)
   {
   continue;
   }
   if (abruijn::AbruijnGraph::DEBUG_EXPANSION)
   {
   DEBUG_MSG("Calling addAbruijnEdges on " << i << " -> " << j << " (" << (*specSet)[specIdx][i][0] << " -> " << (*specSet)[specIdx][j][0] << ")")
   }
   list<abruijn::AbruijnEdge*> addedEdges;
   abGraph.addAbruijnEdges((*specSet)[specIdx], i, j, &addedEdges);

   if (addedEdges.size() > 0)
   {
   tempMap[i].insert(j);
   }
   }
   }
   connectedEdges[specIdx] = tempMap;
   }

   for (unsigned int i = 0; i < vertices.size(); i++)
   {
   if (vertices[i].specPeaks.size() == 0)
   {
   continue;
   }

   for (list<SimpleEdge>::iterator eIt = vertices[i].outEdges.begin(); eIt
   != vertices[i].outEdges.end(); eIt++)
   {
   int specIdx = eIt->specIdx;
   int peakIdxFrom = eIt->peaks[0];
   int peakIdxTo = eIt->peaks[1];

   if (connectedEdges[specIdx][peakIdxFrom].count(peakIdxTo) == 0)
   {
   abruijn::AbruijnNode* nodeFrom = peakToNode[specIdx][peakIdxFrom];
   abruijn::AbruijnNode* nodeTo = peakToNode[specIdx][peakIdxTo];
   abGraph.addLabelFreeAbEdge(nodeFrom,
   nodeTo,
   eIt->edgeMass,
   eIt->edgeScore);
   }
   }
   }

   DEBUG_MSG("Merging edges ...");
   // merge parallel edges
   abGraph.mergeAllParallelLabelFreeEdgesByMass(peakTol);
   abGraph.mergeAllParallelEdges();

   DEBUG_MSG("Expanding forward paths ...");
   // expand forward parallel paths
   abGraph.expandAllPaths(false);

   DEBUG_MSG("Merging edges ...");
   // merge parallel edges
   abGraph.mergeAllParallelLabelFreeEdgesByMass(peakTol);
   abGraph.mergeAllParallelEdges();

   DEBUG_MSG("Expanding reverse paths ...");
   // expand reverse parallel paths
   abGraph.expandAllPaths(true);

   DEBUG_MSG("Merging edges ...");
   // merge parallel edges
   abGraph.mergeAllParallelLabelFreeEdgesByMass(peakTol);
   abGraph.mergeAllParallelEdges();

   //abGraph.saveGraphviz("postExpand.dot");

   if (abGraph.computeHeaviestPathDAG())
   {
   DEBUG_MSG("Found heaviest path!!");
   for (unsigned int i = 0; i < vertices.size(); i++)
   {
   if (vertices[i].specPeaks.size() == 0)
   {
   continue;
   }
   vertices[i].outEdges.clear();
   }

   bool inBreak = false;
   float totMass = 0;
   float totScore = 0;
   int startVert = -1;
   for (unsigned int i = 0; i < abGraph.m_consensusPath.size(); i++)
   {
   abruijn::AbruijnEdge* abEdge =
   abruijn::AbruijnEdge::castEdgePtr(abGraph.m_consensusPath[i]);
   abruijn::AbruijnNode* abNodeFrom =
   abruijn::AbruijnNode::castNodePtr(abEdge->getFrom());
   abruijn::AbruijnNode* abNodeTo =
   abruijn::AbruijnNode::castNodePtr(abEdge->getTo());

   bool havePeaks = false;
   for (unsigned int i = 0; i < abNodeFrom->size(); i++)
   {
   if ((*abNodeFrom)[i].getIntensity() > -0.01)
   {
   havePeaks = true;
   break;
   }
   }

   if (havePeaks && nodeToVertIdx.count(abNodeFrom) == 0)
   {
   nodeToVertIdx[abNodeFrom] = addAbruijnNode(abNodeFrom, specIdToIdx);
   }

   if (havePeaks)
   {
   startVert = nodeToVertIdx[abNodeFrom];
   }

   havePeaks = false;
   for (unsigned int i = 0; i < abNodeTo->size(); i++)
   {
   if ((*abNodeTo)[i].getIntensity() > -0.01)
   {
   havePeaks = true;
   break;
   }
   }

   if (havePeaks && nodeToVertIdx.count(abNodeTo) == 0)
   {
   nodeToVertIdx[abNodeTo] = addAbruijnNode(abNodeTo, specIdToIdx);
   ;
   }

   if (!havePeaks && !inBreak)
   {
   inBreak = true;
   }
   totScore += abEdge->getWeight();
   totMass += abEdge->getMass();

   if (havePeaks)
   {
   int destVertIdx = nodeToVertIdx[abNodeTo];
   SimpleEdge newEdge(totMass, totScore, -1, -1, -1, destVertIdx);
   vertices[startVert].outEdges.push_back(newEdge);
   totScore = 0;
   totMass = 0;
   }
   if (inBreak && havePeaks)
   {
   inBreak = false;
   }
   }
   return;
   }

   DEBUG_MSG("Removing duplicate edges ...");
   // remove repeated edges between the same nodes with the same jump mass
   abGraph.pruneAllParallelEdgesByMass(peakTol);

   DEBUG_MSG("Removing empty nodes ...");
   // replace 'dummy' nodes with multi-residue jumps
   abGraph.removeEmptyNodes();

   DEBUG_MSG("Removing duplicate edges ...");
   // remove repeated edges between the same nodes with the same jump mass
   abGraph.pruneAllParallelEdgesByMass(peakTol);

   DEBUG_MSG("Adding in expanded paths ...");

   for (NodeSet::iterator nIt = abGraph.beginNodes(); nIt
   != abGraph.endNodes(); nIt++)
   {
   abruijn::AbruijnNode* abNode = abruijn::AbruijnNode::castNodePtr(*nIt);
   if (nodeToVertIdx.count(abNode) > 0)
   {
   continue;
   }

   nodeToVertIdx[abNode] = addAbruijnNode(abNode, specIdToIdx);
   }
   for (unsigned int i = 0; i < vertices.size(); i++)
   {
   if (vertices[i].specPeaks.size() == 0)
   {
   continue;
   }
   vertices[i].outEdges.clear();
   }

   for (EdgeSet::iterator eIt = abGraph.beginEdges(); eIt
   != abGraph.endEdges(); eIt++)
   {
   abruijn::AbruijnNode* abNodeFrom =
   abruijn::AbruijnNode::castNodePtr((*eIt)->getFrom());
   abruijn::AbruijnNode* abNodeTo =
   abruijn::AbruijnNode::castNodePtr((*eIt)->getTo());

   // Only copy edges in expanded paths
   if ((copyNodes.count(abNodeFrom) > 0) && (copyNodes.count(abNodeTo) > 0))
   {
   continue;
   }

   abruijn::AbruijnEdge* abEdge = abruijn::AbruijnEdge::castEdgePtr(*eIt);

   if (nodeToVertIdx.count(abNodeTo) == 0 || nodeToVertIdx[abNodeTo] < 0)
   {
   ERROR_MSG("Failed to locate vertex for " << abNodeTo->getGraphvizLabel());
   abort();
   }
   if (nodeToVertIdx.count(abNodeFrom) == 0 || nodeToVertIdx[abNodeFrom] < 0)
   {
   ERROR_MSG("Failed to locate vertex for " << abNodeFrom->getGraphvizLabel());
   abort();
   }
   int destVertIdx = nodeToVertIdx[abNodeTo];
   int origVertIdx = nodeToVertIdx[abNodeFrom];

   SimpleEdge newEdge(abEdge->getMass(),
   abEdge->getWeight(),
   -1,
   -1,
   -1,
   destVertIdx);

   vertices[origVertIdx].outEdges.push_back(newEdge);
   }
   }*/

  int VertexSet::addAbruijnNode(abruijn::AbruijnNode* abNode, map<string,
      unsigned int>& specIdToIdx)
  {
    int vertIdx = -1;
    for (unsigned int i = 0; i < abNode->size(); i++)
    {
      if ((*abNode)[i].getIntensity() < -0.1)
      {
        continue;
      }
      int specIdx = specIdToIdx[(*abNode)[i].getSpecID()];
      int peakIdx = (*specSet)[specIdx].findClosest((*abNode)[i].getMass());

      if (vertIdx < 0)
      {
        vertIdx = addVertex(specIdx, peakIdx);
      }
      else
      {
        vertices[vertIdx].addPeak(TwoValues<int> (specIdx, peakIdx));
      }
    }

    if (vertIdx < 0)
    {
      ERROR_MSG("Failed to locate a real peak in: " << abNode->toString());
      abort();
    }
    return vertIdx;
  }

  //
  // enumerateHeaviestPath - Uses the heaviest path in g (as marked in g.edgeInPath) to
  //   generate sets of selected spectrum peaks from the current vertex set / specSet
  //   On output, enumSpecs contains the selected peaks and enumSpecsIdx the indices of the
  //     spectra these peaks were selected from.
  void VertexSet::enumerateHeaviestPath(MSGraph &g,
                                        vector<int> &vSet_index,
                                        SpecSet &enumSpecs,
                                        vector<int> &enumSpecsIdx)
  {

    // Convert heaviest path to list of matched vertices
    vector<bool> matchedVerts(g.vNext.size());
    for (unsigned int i = 0; i < matchedVerts.size(); i++)
      matchedVerts[i] = false;
    for (unsigned int i = 0; i < g.edges.size(); i++)
      if (g.edgeInPath[i])
      {
        matchedVerts[g.edges[i][0]] = true;
        matchedVerts[g.edges[i][1]] = true;
      }

    list<int> vSet_verts; // List of matched vertices in the vertex set
    for (unsigned int i = 0; i < matchedVerts.size(); i++)
      if (matchedVerts[i])
        vSet_verts.push_front(vSet_index[i]);

    // Find the indices of the spectra whose peaks were matched, fill in enumSpecsIdx
    list<int> matchedSpecsIdx;
    vector<bool> matchedSpecs(specSet->size());
    for (unsigned int i = 0; i < matchedSpecs.size(); i++)
      matchedSpecs[i] = false;
    vector<vector<bool> > matchedPeaks(specSet->size());
    list<TwoValues<int> >::iterator iter;
    list<int>::iterator vertexIter;
    for (vertexIter = vSet_verts.begin(); vertexIter != vSet_verts.end(); vertexIter++)
      for (iter = vertices[*vertexIter].specPeaks.begin(); iter
          != vertices[*vertexIter].specPeaks.end(); iter++)
      {
        if (!matchedSpecs[(*iter)[0]])
        {
          matchedSpecs[(*iter)[0]] = true;
          matchedSpecsIdx.push_front((*iter)[0]);
          matchedPeaks[(*iter)[0]].resize((*specSet)[(*iter)[0]].size());
        }
        matchedPeaks[(*iter)[0]][(*iter)[1]] = true;
      }

    // Convert list of matched vertices to list of matched peaks per spectrum
    enumSpecs.resize(matchedSpecsIdx.size());
    enumSpecsIdx.resize(matchedSpecsIdx.size());
    list<int>::iterator specIter = matchedSpecsIdx.begin();
    int esIdx = 0;
    for (; specIter != matchedSpecsIdx.end(); specIter++)
    {
      int numMatchedPeaks = 0;
      for (unsigned int i = 0; i < matchedPeaks[*specIter].size(); i++)
        if (matchedPeaks[*specIter][i])
          numMatchedPeaks++;
      enumSpecs[esIdx].copyNP((*specSet)[*specIter]);
      enumSpecs[esIdx].resize(numMatchedPeaks);
      numMatchedPeaks = 0;
      for (unsigned int i = 0; i < matchedPeaks[*specIter].size(); i++)
        if (matchedPeaks[*specIter][i])
        {
          enumSpecs[esIdx][numMatchedPeaks] = (*specSet)[*specIter][i];
          numMatchedPeaks++;
        }
      enumSpecsIdx[esIdx++] = *specIter;
    }
  }

  void VertexSet::getMatchedPeaks(vector<int> &matchedVertices, vector<list<
      TwoValues<int> > > &putHere)
  {
    unsigned int vIdx = 0, phIdx = 0;

    if (matchedVertices.size() > 0)
    {
      putHere.resize(matchedVertices.size());
      for (vIdx = 0; vIdx < matchedVertices.size(); vIdx++)
      {
        putHere[vIdx].clear();
        putHere[vIdx].assign(vertices[matchedVertices[vIdx]].specPeaks.begin(),
                             vertices[matchedVertices[vIdx]].specPeaks.end());
      }
    }
    else
    {
      putHere.resize(vertices.size());
      phIdx = 0;
      for (vIdx = 0; vIdx < vertices.size(); vIdx++)
        if (not vertices[vIdx].specPeaks.empty())
        {
          putHere[phIdx].clear();
          putHere[phIdx].assign(vertices[vIdx].specPeaks.begin(),
                                vertices[vIdx].specPeaks.end());
          phIdx++;
        }
      putHere.resize(phIdx);
    }
  }

  void VertexSet::outputGraph(vector<SpectrumPeakLabels> &labels)
  {
    PeakLabel vLabel;
    string vLabelText("<empty>");

    for (unsigned int v = 0; v < vertices.size(); v++)
    {
      if (vertices[v].size() == 0)
        continue;
      if (labels.size() > 0)
      {
        vLabel.set(0, 0, 0, 0, 0);
        list<TwoValues<int> >::iterator iterP = vertices[v].specPeaks.begin();
        for (; iterP != vertices[v].specPeaks.end(); iterP++)
          vLabel.merge(labels[(*iterP)[0]][(*iterP)[1]]);
        if (vLabel.isEmpty())
          vLabelText.assign("<empty>");
        else
          vLabelText = vLabel.asString();
      }

      cout << "Vertex " << v << ", label " << vLabelText << ":\n";
      cout << " --- vertices: ";
      cout.flush();
      list<TwoValues<int> >::iterator iterP = vertices[v].specPeaks.begin();
      for (; iterP != vertices[v].specPeaks.end(); iterP++)
        cout << "(" << (*iterP)[0] << "," << (*iterP)[1] << ","
            << (*specSet)[(*iterP)[0]][(*iterP)[1]][0] << ") ";
      cout << "\n --- match edges: ";
      cout.flush();
      list<MatchEdge>::iterator iterE = vertices[v].outMatchEdges.begin();
      for (; iterE != vertices[v].outMatchEdges.end(); iterE++)
        cout << "[(" << (*iterE).from.specPeak1[0] << ","
            << (*iterE).from.specPeak1[1] << "," << (*iterE).from.specPeak2[0]
            << "," << (*iterE).from.specPeak2[1] << ")->("
            << (*iterE).to.specPeak1[0] << "," << (*iterE).to.specPeak1[1]
            << "," << (*iterE).to.specPeak2[0] << ","
            << (*iterE).to.specPeak2[1] << ")] ";
      cout << "\n --- simple edges: ";
      cout.flush();
      list<SimpleEdge>::iterator iterE2 = vertices[v].outEdges.begin();
      for (; iterE2 != vertices[v].outEdges.end(); iterE2++)
        cout << "[" << (*iterE2).specIdx << "," << (*iterE2).peaks[0] << "->"
            << (*iterE2).peaks[1] << "," << (*iterE2).edgeMass << ","
            << (*iterE2).edgeScore << "] ";
      cout << endl;
    }
  }

  // Output an A-Bruijn multiple alignment in graphviz format
  void VertexSet::output_graphviz_ma(char *filename,
                                     vector<int> &matchedVertices)
  {
    ofstream out(filename, ios::binary);
    if (!out)
    {
      cerr << "Error opening " << filename << "!\n";
      return;
    }

    /*	for(unsigned int vIdx=0; vIdx<matchedVertices.size(); vIdx++) {
     out << "Vertex " << matchedVertices[vIdx] << ": ";
     for(list<TwoValues<int> >::iterator iter=vertices[vIdx].specPeaks.begin(); iter!=vertices[vIdx].specPeaks.end(); iter++)
     out<<"("<<(*iter)[0]<<","<<(*iter)[1]<<","<<(*specSet)[(*iter)[0]][(*iter)[1]][0]<<"),";
     }*/

    list<int> specIndices; // Unique indices of all spectra with peaks in matchedVertices
    for (unsigned int vIdx = 0; vIdx < matchedVertices.size(); vIdx++)
      for (list<TwoValues<int> >::iterator iter =
          vertices[matchedVertices[vIdx]].specPeaks.begin(); iter
          != vertices[matchedVertices[vIdx]].specPeaks.end(); iter++)
        specIndices.push_back((*iter)[0] + 1);
    specIndices.sort();
    specIndices.unique();

    out << "digraph G {\n";
    out << "  ranksep=.75; size = \"7.5,7.5\";\n";
    out << "  { node [shape=plaintext, fontsize=16];\n";
    list<int>::iterator indicesIter = specIndices.begin();
    out << "  \"" << *indicesIter;
    indicesIter++;
    for (; indicesIter != specIndices.end(); indicesIter++)
      out << "\" -> \"" << *indicesIter;
    out << "\";\n";
    out << "  }\n";

    list<TwoValues<int> >::iterator iterPrev;
    bool firstVertex;
    vector<int> lastUsedPeak(specSet->size());
    for (unsigned int specIdx = 0; specIdx < lastUsedPeak.size(); specIdx++)
      lastUsedPeak[specIdx] = -1;
    for (unsigned int vIdx = 0; vIdx < matchedVertices.size(); vIdx++)
    {
      firstVertex = true;
      for (list<TwoValues<int> >::iterator iter =
          vertices[matchedVertices[vIdx]].specPeaks.begin(); iter
          != vertices[matchedVertices[vIdx]].specPeaks.end(); iter++)
      {
        //			out<<"("<<(*iter)[0]<<","<<(*iter)[1]<<","<<(*specSet)[(*iter)[0]][(*iter)[1]][0]<<"),";
        out << "  \"" << (*iter)[0] + 1 << "-" << (*iter)[1]
            << "\" [shape=ellipse, label=\""
            << (*specSet)[(*iter)[0]][(*iter)[1]][0] << "\"];\n";
        out << "  { rank=same; \"" << (*iter)[0] + 1 << "\"; \"" << (*iter)[0]
            + 1 << "-" << (*iter)[1] << "\" }\n";

        if (lastUsedPeak[(*iter)[0]] >= 0)
          //				out<<"  \""<<(*iter)[0]+1<<"-"<<lastUsedPeak[(*iter)[0]]<<"\" -> \""<<(*iter)[0]+1<<"-"<<(*iter)[1]<<"\" [label=\""<<(*specSet)[(*iter)[0]][(*iter)[1]][0]-(*specSet)[(*iter)[0]][lastUsedPeak[(*iter)[0]]][0]<<"\", minlen=0, constraint=false];\n";
          //				out<<"  \""<<(*iter)[0]+1<<"-"<<lastUsedPeak[(*iter)[0]]<<"\" -> \""<<(*iter)[0]+1<<"-"<<(*iter)[1]<<"\";\n";
          out << "  \"" << (*iter)[0] + 1 << "-" << lastUsedPeak[(*iter)[0]]
              << "-" << (*iter)[1] << "\" [shape=diamond, label=\""
              << (*specSet)[(*iter)[0]][(*iter)[1]][0]
                  - (*specSet)[(*iter)[0]][lastUsedPeak[(*iter)[0]]][0]
              << "\"];\n" << "  { rank=same; \"" << (*iter)[0] + 1 << "\"; \""
              << (*iter)[0] + 1 << "-" << lastUsedPeak[(*iter)[0]] << "-"
              << (*iter)[1] << "\" }\n" << "  \"" << (*iter)[0] + 1 << "-"
              << lastUsedPeak[(*iter)[0]] << "\" -> \"" << (*iter)[0] + 1
              << "-" << lastUsedPeak[(*iter)[0]] << "-" << (*iter)[1]
              << "\";\n" << "  \"" << (*iter)[0] + 1 << "-"
              << lastUsedPeak[(*iter)[0]] << "-" << (*iter)[1] << "\" -> \""
              << (*iter)[0] + 1 << "-" << (*iter)[1] << "\";\n";
        lastUsedPeak[(*iter)[0]] = (*iter)[1];

        if (firstVertex)
          firstVertex = false;
        else
          out << "  \"" << (*iterPrev)[0] + 1 << "-" << (*iterPrev)[1]
              << "\" -> \"" << (*iter)[0] + 1 << "-" << (*iter)[1]
              << "\" [arrowhead=none, minlen=0, constraint=false];\n\n";
        iterPrev = iter;
      }
    }
    out << "}\n";
    out.close();
  }

  void Save_abinfo_v1_0(const char *filename, SpecSet &specSet, vector<
      list<int> > &cSpectra, vector<bool> &specFlipped, vector<vector<list<
      TwoValues<int> > > > &abVertices);
  void Save_abinfo_v1_1(const char *filename, SpecSet &specSet, vector<
      list<int> > &cSpectra, vector<bool> &specFlipped, vector<vector<list<
      TwoValues<int> > > > &abVertices);

  void Save_abinfo(const char *filename,
                   SpecSet &specSet,
                   vector<list<int> > &cSpectra,
                   vector<bool> &specFlipped,
                   vector<vector<list<TwoValues<int> > > > &abVertices)
  {
    Save_abinfo_v1_0(filename, specSet, cSpectra, specFlipped, abVertices);
  }


  void getAssembledShifts(SpecSet& inputContigSpectra,
                          abinfo_t& inputContigAbinfo,
                          vector<vector<sps::tuple<unsigned int, float, bool> > >& outputAssembledShifts)
  {
    outputAssembledShifts.clear();
    outputAssembledShifts.resize(inputContigSpectra.size());

    if (inputContigSpectra.size() == 0)
    {
      return;
    }

    for (unsigned int i = 0; i < inputContigSpectra.size(); i++)
    {
      if (inputContigSpectra[i].size() == 0 || inputContigAbinfo.count(i) == 0)
      {
        continue;
      }

      int scanIndex = inputContigSpectra[i].scan - 1;
      if (scanIndex < 0) {
        continue;  // Should never happen.. but let's be safe
      }
      vector<std::pair<vector<int> , vector<double> > >& abVerts =
          inputContigAbinfo[scanIndex].second;

      vector<int>& specsRev = inputContigAbinfo[scanIndex].first.second;
      vector<int>& specsIdx = inputContigAbinfo[scanIndex].first.first;
      map<int, bool> specMapRev;

      for (unsigned int j = 0; j < specsIdx.size(); j++)
      {
        specMapRev[specsIdx[j]] = specsRev[j] == 1;
      }

      outputAssembledShifts[i].resize(specsIdx.size());
      unsigned int idxUse = 0;

      set<int> addedSpecs;

      for (unsigned int j = 0; j < abVerts.size() && idxUse
          < outputAssembledShifts[i].size(); j++)
      {
        float abShift = inputContigSpectra[i][j][0];

        vector<int>& abSpecs = abVerts[j].first;
        vector<double> abPeakMasses = abVerts[j].second;

        for (unsigned int k = 0; k < abSpecs.size() && idxUse
            < outputAssembledShifts[i].size(); k++)
        {
          if (addedSpecs.count(abSpecs[k]) > 0)
          {
            continue;
          }
          addedSpecs.insert(abSpecs[k]);
          outputAssembledShifts[i][idxUse].m0 = abSpecs[k];
          outputAssembledShifts[i][idxUse].m1 = abShift - abPeakMasses[k];
          outputAssembledShifts[i][idxUse++].m2 = specMapRev[abSpecs[k]];
        }
      }
    }

    return;
  }

  void Merge_abinfo_v1_0(SpecSet& specSet,
                         vector<list<int> >& cSpectra,
                         vector<bool>& specFlipped,
                         vector<vector<list<TwoValues<int> > > >& abVertices,
                         abinfo_t& merged_abinfo)
  {
    merged_abinfo.clear();
    pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> , vector<
        double> > > > comp_info;

    unsigned int cIdx, vIdx;
    int tempIdx, tempSize, flipped;
    list<TwoValues<int> >::iterator ltiIter;

    for (cIdx = 0; cIdx < cSpectra.size(); cIdx++)
    {
      set<int> usedSpecs;
      merged_abinfo[cIdx] = comp_info;
      merged_abinfo[cIdx].second.resize(abVertices[cIdx].size());
      for (vIdx = 0; vIdx < abVertices[cIdx].size(); vIdx++)
      {
        tempSize = abVertices[cIdx][vIdx].size();
        merged_abinfo[cIdx].second[vIdx].first.resize(tempSize);
        merged_abinfo[cIdx].second[vIdx].second.resize(tempSize);
        tempIdx = 0;
        for (ltiIter = abVertices[cIdx][vIdx].begin(); ltiIter
            != abVertices[cIdx][vIdx].end(); ltiIter++)
        {
          merged_abinfo[cIdx].second[vIdx].first[tempIdx] = (*ltiIter)[0];
          usedSpecs.insert((*ltiIter)[0]);
          merged_abinfo[cIdx].second[vIdx].second[tempIdx]
              = specSet[(*ltiIter)[0]][(*ltiIter)[1]][0];
          tempIdx++;
        }
      }
      merged_abinfo[cIdx].first.first.resize(usedSpecs.size());
      merged_abinfo[cIdx].first.second.resize(usedSpecs.size());
      tempIdx = 0;
      for (set<int>::const_iterator idxIt = usedSpecs.begin(); idxIt
          != usedSpecs.end(); idxIt++)
      {
        merged_abinfo[cIdx].first.first[tempIdx] = *idxIt;
        flipped = (specFlipped[*idxIt]) ? 1 : 0;
        merged_abinfo[cIdx].first.second[tempIdx] = flipped;
        tempIdx++;
      }
    }
  }

  void Append_abinfo_v1_0(list<abinfo_t>& inabinfo,
                          list<int>& specCount,
                          abinfo_t& outabinfo)
  {
    outabinfo.clear();
    list<abinfo_t>::iterator abinfo_it = inabinfo.begin();
    list<int>::iterator countIt = specCount.begin();
    int num_spec_in = 0, size_in = 0;
    abinfo_t::iterator compit;
    abinfo_t in_specs;

    while (abinfo_it != inabinfo.end())
    {
      in_specs = *abinfo_it;
      for (compit = in_specs.begin(); compit != in_specs.end(); compit++)
      {
        pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> ,
            vector<double> > > > second = compit->second;
        for (int j = 0; j < second.first.first.size(); j++)
        {
          second.first.first[j] += num_spec_in;
        }
        for (int j = 0; j < second.second.size(); j++)
        {
          pair<vector<int> , vector<double> > vert = second.second[j];
          for (int k = 0; k < vert.first.size(); k++)
            second.second[j].first[k] += num_spec_in;
        }
        outabinfo[compit->first + size_in] = second;
      }
      num_spec_in += *countIt;
      size_in += in_specs.size();
      abinfo_it++;
      countIt++;
    }
  }

  /**
   * Combines abinfo when contigs have been merged to meta-contigs
   * @param childAbinfo abinfo of original contigs
   * @param contig_rev whether original contigs have been reversed
   * @param star_spectra spectra assembled into original contigs
   * @param contigsMerged whether each original contig was merged
   * @param contigs original contigs
   * @param parentAbinfo abinfo of parent contigs (no original contigs
   *   should be marked as reversed here)
   * @param outabinfo abinfo output mapping meta-contigs and their abruijn
   *   vertices to star spectrum indices
   * @return
   */
  void Combine_abinfo_v1_0(abinfo_t& childAbinfo,
                           vector<bool>& contig_rev,
                           SpecSet& star_spectra,
                           vector<bool>& contigsMerged,
                           SpecSet& contigs,
                           abinfo_t& parentAbinfo,
                           abinfo_t& outabinfo)
  {
    outabinfo.clear();
    //int debugIdx = 8;
    abinfo_t::iterator compit;

    vector<bool> stars_reversed_child(star_spectra.size(), false);

    star_spectra.addZPMpeaks(0.1, 0, true);

    for (compit = childAbinfo.begin(); compit != childAbinfo.end(); compit++)
    {
      if (compit->second.second.size() == 0)
        continue;

      map<int, bool> childReversed;
      for (int j = 0; j < compit->second.first.first.size(); j++)
      {
        childReversed[compit->second.first.first[j]]
            = compit->second.first.second[j] == 1;
      }

      for (int j = 0; j < compit->second.second.size(); j++)
      {
        pair<vector<int> , vector<double> >* vert = &compit->second.second[j];
        for (int k = 0; k < vert->first.size(); k++)
        {
          int specIdx = (*vert).first[k];
          float mass = (*vert).second[k];
          int peakIdx = star_spectra[specIdx].findClosest(mass);

          if (childReversed[specIdx]
              || !isEqual(mass, star_spectra[specIdx][peakIdx][0], 0.01))
          {
            stars_reversed_child[specIdx] = true;
          }
        }
      }
    }

    vector<bool> stars_reversed_parent(stars_reversed_child);

    set<int> starIdxIncl;

    for (compit = parentAbinfo.begin(); compit != parentAbinfo.end(); compit++)
    {

      if (compit->second.second.size() == 0)
        continue;

      int masterIdx = compit->first;

      pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> , vector<
          double> > > > second = compit->second;

      vector<pair<vector<int> , vector<double> > > verts;
      starIdxIncl.clear();
      /*    if (masterIdx == debugIdx) {
       cout << "Meta-Component " << debugIdx << ":\n";
       for (int j = 0; j < second.first.first.size(); j++) {
       cout << "contig " << second.first.first[j] << ":\n";
       contigs[second.first.first[j]].output(cout);
       cout << "\n";
       }
       }*/
      for (int j = 0; j < second.second.size(); j++)
      {
        vector<int> specIdxs;
        vector<double> peaksMasses;
        //if (masterIdx == debugIdx) {cout << "Vertex " << j << ": ";}
        for (int k = 0; k < second.second[j].first.size(); k++)
        {
          int contigIdx = second.second[j].first[k];
          float contigPeakMass = (float)second.second[j].second[k];
          Spectrum contig = contigs[contigIdx];
          if (contig.size() == 0)
            continue;
          bool rev = contig_rev[contigIdx] == 1;
          int contigPeakIdx;

          if (rev)
          {
            contig.reverse(0.0 - AAJumps::massH2O, 0);
            contigPeakIdx = contig.findClosest(contigPeakMass);
            contigPeakIdx = contig.size() - 1 - contigPeakIdx;
          }
          else
          {
            contigPeakIdx = contig.findClosest(contigPeakMass);
          }
          /*
           if (!isEqual(peakMass, use_contig[contigPeakIdx][0], 0.01)) {
           cerr << "Expected mass does not match theoretical for contig "
           << contigIdx << " (" << rev << ") : " << peakMass << " != "
           << use_contig[contigPeakIdx][0] << "\n";
           use_contig.output(cerr);
           cerr << "\n";
           }
           */
          /*       if (masterIdx == debugIdx) {
           cout << "(contig=" << contigIdx << " peak=" << contigPeakMass << " peakIdx=" << contigPeakIdx << " rev=" << rev << ":: ";
           }*/
          for (int l = 0; l
              < childAbinfo[(unsigned int)contigIdx].second[contigPeakIdx].first.size(); l++)
          {
            int
                specIdx =
                    childAbinfo[(unsigned int)contigIdx].second[contigPeakIdx].first[l];
            Spectrum star = star_spectra[specIdx];
            float
                peakMass =
                    childAbinfo[(unsigned int)contigIdx].second[contigPeakIdx].second[l];
            /*
             if (masterIdx == debugIdx) {
             cout << "specIdx=" << specIdx << " peakMass=" << peakMass << ", ";
             }*/

            stars_reversed_parent[specIdx] = stars_reversed_child[specIdx]
                != rev;
            if (rev)
            {
              peakMass = star.parentMass - AAJumps::massHion - peakMass;
            }

            starIdxIncl.insert(specIdx);

            specIdxs.push_back(specIdx);
            peaksMasses.push_back(peakMass);
          }
          //if (masterIdx == debugIdx) {cout << ")\n";}
        }
        //if (masterIdx == debugIdx) {cout << "\n\n";}
        pair<vector<int> , vector<double> > peaks;
        peaks.first = specIdxs;
        peaks.second = peaksMasses;
        verts.push_back(peaks);
      }

      vector<int> idxs(starIdxIncl.size());
      vector<int> revs(starIdxIncl.size());
      int locIdx = 0;

      for (set<int>::iterator starIt = starIdxIncl.begin(); starIt
          != starIdxIncl.end(); starIt++)
      {
        /*if (masterIdx == debugIdx) {
         cout << "\nstar " << *starIt << ":\n";
         star_spectra[*starIt].output(cout);
         }*/
        idxs[locIdx] = *starIt;
        revs[locIdx] = (stars_reversed_parent[*starIt]) ? 1 : 0;
        locIdx++;
      }

      for (int j = 0; j < second.first.first.size(); j++)
        if (second.first.second[j] == 1)
          cerr << "ATTENTION: masab reversed contig " << second.first.first[j]
              << "\n";

      pair<vector<int> , vector<int> > firstIdx;
      firstIdx.first = idxs;
      firstIdx.second = revs;
      second.first = firstIdx;
      second.second = verts;
      outabinfo[compit->first] = second;
    }

    for (int i = 0; i < contigsMerged.size(); i++)
    {
      if (!contigsMerged[i])
        outabinfo[i] = childAbinfo[i];
    }
  }

  // *****************************************************************************
  //   Save_abinfo_v1_0 - saves information about the components in binary format:
  //
  //     nNumComponents [int]
  //     per component:
  //       num_specs [short]
  //         per spec: specIndex [int, 0-based], specFlipped [short 0/1]
  //       num_ABruijn_vertices [short]
  //         per vertex:
  //           num_peaks_in_vertex [short]
  //           per peak: specIndex [int, 0-based], peakMass [float]
  //
  //     numUsedSpectra [int]
  //     numUsedSpectra-by-1 unsigned int array of all used spectrum indices
  //     numUsedSpectra-by-2 unsigned short array of [component index,specFlipped]
  //     numComponents [int]
  //     numComponents-by-1 unsigned short array of number of ABruijn vertices per component
  //     numABVertices-by-1 unsigned short array of number of spectrum peaks per ABruijn vertex
  //     totNumSpecPeaks-by-1 unsigned int array of spectrum index per peak
  //     totNumSpecPeaks-by-1 float array of peak masses per spectrum peak
  //
  // *****************************************************************************
  void Save_abinfo_v1_0(const char *filename, SpecSet &specSet, vector<
      list<int> > &cSpectra, vector<bool> &specFlipped, vector<vector<list<
      TwoValues<int> > > > &abVertices)
  {
    FILE *fp = fopen(filename, "wb");
    if (fp == 0)
    {
      cerr << "Error opening " << filename << "!!\n";
      return;
    }
    unsigned int cIdx, vIdx;
    list<TwoValues<int> >::iterator ltiIter;

    // Prepare and write the spectrum info part
    unsigned int specCount = 0;
    for (cIdx = 0; cIdx < cSpectra.size(); cIdx++)
      specCount += cSpectra[cIdx].size();
    unsigned int *specIndices = (unsigned int *)malloc(specCount
        * sizeof(unsigned int)); // Vector of indices of all used spectra
    unsigned short *specInfo = (unsigned short *)malloc(2 * specCount
        * sizeof(unsigned short)); // Sequence of [component index, specFlipped]
    unsigned int indIdx = 0, infoIdx = 0;
    for (cIdx = 0; cIdx < cSpectra.size(); cIdx++)
    {
      for (list<int>::iterator liIter = cSpectra[cIdx].begin(); liIter
          != cSpectra[cIdx].end(); liIter++)
      {
        specIndices[indIdx++] = (unsigned int)*liIter;
        specInfo[infoIdx++] = (unsigned short)cIdx;
        specInfo[infoIdx++] = (unsigned short)specFlipped[*liIter];
      }
    }

    fwrite(&specCount, sizeof(unsigned int), 1, fp);
    fwrite(specIndices, sizeof(unsigned int), specCount, fp);
    free(specIndices);
    fwrite(specInfo, sizeof(unsigned short), 2 * specCount, fp);
    free(specInfo);

    // Prepare and write the ABruijn vertices info part
    unsigned int numComponents = cSpectra.size(), numVertices = 0, numPeaks = 0;
    for (cIdx = 0; cIdx < numComponents; cIdx++)
    {
      numVertices += abVertices[cIdx].size();
      for (vIdx = 0; vIdx < abVertices[cIdx].size(); vIdx++)
        numPeaks += abVertices[cIdx][vIdx].size();
    }
    unsigned short *vertsPerComp = (unsigned short *)malloc(numComponents
        * sizeof(unsigned short));
    unsigned short *peaksPerVert = (unsigned short *)malloc(numVertices
        * sizeof(unsigned short));
    unsigned int *vertsSpecIdx = (unsigned int *)malloc(numPeaks
        * sizeof(unsigned int));
    //	unsigned short *vertsPeakIdx = (unsigned short *) malloc(numPeaks*sizeof(unsigned short));
    float *vertsPeakMass = (float *)malloc(numPeaks * sizeof(float));
    unsigned int gVidx = 0, gPidx = 0; // Global vertex/vertexPeak indices
    for (cIdx = 0; cIdx < numComponents; cIdx++)
    {
      vertsPerComp[cIdx] = abVertices[cIdx].size();
      for (vIdx = 0; vIdx < abVertices[cIdx].size(); vIdx++)
      {
        peaksPerVert[gVidx++] = abVertices[cIdx][vIdx].size();
        for (ltiIter = abVertices[cIdx][vIdx].begin(); ltiIter
            != abVertices[cIdx][vIdx].end(); ltiIter++)
        {
          vertsSpecIdx[gPidx] = (unsigned int)(*ltiIter)[0];
          vertsPeakMass[gPidx] = specSet[(*ltiIter)[0]][(*ltiIter)[1]][0];
          gPidx++;
        }
        //				{ vertsSpecIdx[gPidx]=(unsigned int)(*ltiIter)[0];    vertsPeakIdx[gPidx]=(unsigned short)(*ltiIter)[1];    gPidx++; }
      }
    }
    fwrite(&numComponents, sizeof(unsigned int), 1, fp);
    fwrite(vertsPerComp, sizeof(unsigned short), numComponents, fp);
    free(vertsPerComp);
    fwrite(peaksPerVert, sizeof(unsigned short), numVertices, fp);
    free(peaksPerVert);
    fwrite(vertsSpecIdx, sizeof(unsigned int), numPeaks, fp);
    free(vertsSpecIdx);
    //	fwrite(vertsPeakIdx, sizeof(unsigned short), numPeaks, fp);       free(vertsPeakIdx);
    fwrite(vertsPeakMass, sizeof(float), numPeaks, fp);
    free(vertsPeakMass);

    /*
     unsigned int iValue = cSpectra.size();
     unsigned short sValue;
     fwrite(&iValue, sizeof(unsigned int), 1, fp);
     for(cIdx=0; cIdx<cSpectra.size(); cIdx++) {
     // How many/which spectra in the component
     sValue = (unsigned short) cSpectra[cIdx].size();   fwrite(&sValue, sizeof(unsigned short), 1, fp);
     for(list<int>::iterator liIter=cSpectra[cIdx].begin(); liIter!=cSpectra[cIdx].end(); liIter++) {
     iValue = (unsigned int) *liIter;                fwrite(&iValue, sizeof(unsigned int), 1, fp);
     sValue = (unsigned short)specFlipped[iValue];   fwrite(&sValue, sizeof(unsigned short), 1, fp);
     }

     // ABruijn vertices in the component
     sValue = (unsigned short) abVertices[cIdx].size();   fwrite(&sValue, sizeof(unsigned short), 1, fp);
     for(vIdx=0; vIdx<abVertices[cIdx].size(); vIdx++) {
     // Output peaks per ABruijn vertex
     sValue = (unsigned short) abVertices[cIdx][vIdx].size();   fwrite(&sValue, sizeof(unsigned short), 1, fp);
     for(ltiIter=abVertices[cIdx][vIdx].begin(); ltiIter!=abVertices[cIdx][vIdx].end(); ltiIter++) {
     iValue = (unsigned int)  (*liIter)[0];   fwrite(&iValue, sizeof(unsigned int), 1, fp);
     sValue = (unsigned short)(*liIter)[1];   fwrite(&sValue, sizeof(unsigned short), 1, fp);
     }
     }
     }
     */
    fclose(fp);
  }

  // *****************************************************************************
  //   Save_abinfo_v1_1 - saves information about the components in binary format:
  //
  //     nNumComponents [int]
  //     per component:
  //       num_specs [short]
  //         per spec: specIndex [int, 0-based], specFlipped [short 0/1]
  //       num_ABruijn_vertices [short]
  //         per vertex:
  //           num_peaks_in_vertex [short]
  //           per peak: specIndex [int, 0-based], peakIndex [short, 0-based]
  //
  //    0 (zero) [int]
  //    Major version number [2 byte unsigned short]
  //    Minor version number [2 byte unsigned short]
  //    numUsedSpectra [unsigned int]
  //    numUsedSpectra-by-2 unsigned int array of [used spectrum index,component index]
  //    numUsedSpectra-by-1 1-byte array of [specFlipped = 0/1]
  //    numComponents [unsigned int]
  //    numComponents-by-1 unsigned short array of number of ABruijn vertices per component
  //    numABVertices-by-1 unsigned short array of number of spectrum peaks per ABruijn vertex
  //    totNumSpecPeaks-by-1 unsigned int array of spectrum index per peak (per ABruijn vertex)
  //    totNumSpecPeaks-by-1 float array of peak masses per spectrum peak
  //
  // *****************************************************************************
  void Save_abinfo_v1_1(const char *filename, SpecSet &specSet, vector<
      list<int> > &cSpectra, vector<bool> &specFlipped, vector<vector<list<
      TwoValues<int> > > > &abVertices)
  {
    FILE *fp = fopen(filename, "wb");
    if (fp == 0)
    {
      cerr << "Error opening " << filename << "!!\n";
      return;
    }
    unsigned int cIdx, vIdx;
    uint32_t foo32;
    uint16_t foo16;
    list<TwoValues<int> >::iterator ltiIter;

    // Write version information
    foo32 = 0;
    fwrite(&foo32, sizeof(uint32_t), 1, fp); // For backwards compatibility (numUsedSpectra=0)
    foo16 = 1;
    fwrite(&foo16, sizeof(uint16_t), 1, fp); // Major version number
    foo16 = 1;
    fwrite(&foo16, sizeof(uint16_t), 1, fp); // Minor version number

    // Prepare and write the spectrum info part
    uint32_t specCount = 0;
    for (cIdx = 0; cIdx < cSpectra.size(); cIdx++)
      specCount += cSpectra[cIdx].size();
    uint32_t *specIndices =
        (uint32_t *)malloc(2 * specCount * sizeof(uint32_t)); // Vector of spectrum-index/component index
    char *specFlip = (char *)malloc(specCount * sizeof(char)); // Flipped/not-flipped indicators
    unsigned int indIdx = 0, flipIdx = 0;
    for (cIdx = 0; cIdx < cSpectra.size(); cIdx++)
    {
      for (list<int>::iterator liIter = cSpectra[cIdx].begin(); liIter
          != cSpectra[cIdx].end(); liIter++)
      {
        specIndices[indIdx++] = (uint32_t)*liIter;
        specIndices[indIdx++] = (uint32_t)cIdx;
        specFlip[flipIdx++] = (char)specFlipped[*liIter];
      }
    }
    fwrite(&specCount, sizeof(uint32_t), 1, fp);
    fwrite(specIndices, sizeof(uint32_t), 2 * specCount, fp);
    free(specIndices);
    fwrite(specFlip, sizeof(char), specCount, fp);
    free(specFlip);

    // Prepare and write the ABruijn vertices info part
    uint32_t numComponents = cSpectra.size(), numVertices = 0, numPeaks = 0;
    for (cIdx = 0; cIdx < numComponents; cIdx++)
    {
      numVertices += abVertices[cIdx].size();
      for (vIdx = 0; vIdx < abVertices[cIdx].size(); vIdx++)
        numPeaks += abVertices[cIdx][vIdx].size();
    }
    uint16_t *vertsPerComp = (uint16_t *)malloc(numComponents
        * sizeof(uint16_t));
    uint16_t *peaksPerVert = (uint16_t *)malloc(numVertices * sizeof(uint16_t));
    uint32_t *vertsSpecIdx = (uint32_t *)malloc(numPeaks * sizeof(uint32_t));
    //	unsigned short *vertsPeakIdx = (unsigned short *) malloc(numPeaks*sizeof(unsigned short));
    float *vertsPeakMass = (float *)malloc(numPeaks * sizeof(float));
    unsigned int gVidx = 0, gPidx = 0; // Global vertex/vertexPeak indices
    for (cIdx = 0; cIdx < numComponents; cIdx++)
    {
      vertsPerComp[cIdx] = abVertices[cIdx].size();
      for (vIdx = 0; vIdx < abVertices[cIdx].size(); vIdx++)
      {
        peaksPerVert[gVidx++] = abVertices[cIdx][vIdx].size();
        for (ltiIter = abVertices[cIdx][vIdx].begin(); ltiIter
            != abVertices[cIdx][vIdx].end(); ltiIter++)
        {
          vertsSpecIdx[gPidx] = (uint32_t)(*ltiIter)[0];
          vertsPeakMass[gPidx] = specSet[(*ltiIter)[0]][(*ltiIter)[1]][0];
          gPidx++;
        }
      }
    }
    fwrite(&numComponents, sizeof(uint32_t), 1, fp);
    fwrite(vertsPerComp, sizeof(uint16_t), numComponents, fp);
    free(vertsPerComp);
    fwrite(peaksPerVert, sizeof(uint16_t), numVertices, fp);
    free(peaksPerVert);
    fwrite(vertsSpecIdx, sizeof(uint32_t), numPeaks, fp);
    free(vertsSpecIdx);
    fwrite(vertsPeakMass, sizeof(float), numPeaks, fp);
    free(vertsPeakMass);

    fclose(fp);
  }

  /**
   * Physically merges 2 abinfo's into a new one.
   * ASSUMPTIONS:
   * 1.  The two abinfos were two interpretations of the same set of spectra, (1-N), but now we want one of the abinfos to be on spectra 1-N
   * and the other abinfo to be on spectra N+1 - 2N (spectrum 1 = spectrum N+1)
   */
  bool merge_abinfo(abinfo_t & info1,
                    abinfo_t & info2,
                    const char * filename,
                    unsigned int totalSpecCount)
  {

    FILE *fp = fopen(filename, "wb");
    if (fp == 0 || !fp)
    {
      cerr << "Error opening " << filename << "!!\n";
      return false;
    }

    abinfo_t::iterator compit;
    unsigned numComponents1 = 0;
    unsigned int specCount = 0, numComponents = 0, numVertices = 0, numPeaks =
        0;
    vector<std::pair<vector<int> , vector<double> > > peaksPerVertex1;
    ;
    vector<std::pair<vector<int> , vector<double> > > peaksPerVertex2;

    //Count the info from the first abinfo
    for (compit = info1.begin(); compit != info1.end(); compit++)
    {
      specCount += (compit->second).first.first.size();
      numComponents++;
      numVertices += (compit->second).second.size();
      peaksPerVertex1 = (compit->second).second;
      for (int i = 0; i < peaksPerVertex1.size(); i++)
      {
        numPeaks += peaksPerVertex1[i].second.size();
      }
    }

    numComponents1 = numComponents;

    //Count the info from the second abinfo
    for (compit = info2.begin(); compit != info2.end(); compit++)
    {
      specCount += (compit->second).first.first.size();
      numComponents++;
      numVertices += (compit->second).second.size();
      peaksPerVertex2 = (compit->second).second;
      for (int i = 0; i < peaksPerVertex2.size(); i++)
      {
        numPeaks += peaksPerVertex2[i].second.size();
      }
    }

    unsigned int *specIndices = (unsigned int *)malloc(specCount
        * sizeof(unsigned int)); // Vector of indices of all used spectra
    unsigned short *specInfo = (unsigned short *)malloc(2 * specCount
        * sizeof(unsigned short)); // Sequence of [component index, specFlipped]
    unsigned short *vertsPerComp = (unsigned short *)malloc(numComponents
        * sizeof(unsigned short));
    unsigned short *peaksPerVert = (unsigned short *)malloc(numVertices
        * sizeof(unsigned short));
    unsigned int *vertsSpecIdx = (unsigned int *)malloc(numPeaks
        * sizeof(unsigned int));
    float *vertsPeakMass = (float *)malloc(numPeaks * sizeof(float));

    unsigned int indIdx = 0, infoIdx = 0, gVidx = 0, gPidx = 0;
    for (compit = info1.begin(); compit != info1.end(); compit++)
    {

      //specIndexes for the first info do not need to be updated
      vector<int> specIdx = (compit->second).first.first;
      vector<int> revIdx = (compit->second).first.second;

      //contig indices for the first info do not need to be updated
      vertsPerComp[compit->first] = (compit->second).second.size();

      for (int i = 0; i < (compit->second).second.size(); i++)
      {
        peaksPerVert[gVidx++] = (compit->second).second[i].second.size();

        for (int j = 0; j < (compit->second).second[i].second.size(); j++)
        {
          vertsSpecIdx[gPidx]
              = (unsigned int)(compit->second).second[i].first[j];
          vertsPeakMass[gPidx] = (compit->second).second[i].second[j];
          gPidx++;
        }
      }
      for (int i = 0; i < specIdx.size(); i++)
      {
        specIndices[indIdx++] = (unsigned int)specIdx[i];
        specInfo[infoIdx++] = (unsigned short)compit->first;
        specInfo[infoIdx++] = (unsigned short)revIdx[i];
      }
    }

    for (compit = info2.begin(); compit != info2.end(); compit++)
    {

      vector<int> specIdx = (compit->second).first.first;
      vector<int> revIdx = (compit->second).first.second;

      //Component index needs to be updated by the number of abinfo1 components
      vertsPerComp[compit->first + numComponents1]
          = (compit->second).second.size();
      for (int i = 0; i < (compit->second).second.size(); i++)
      {
        peaksPerVert[gVidx++] = (compit->second).second[i].second.size();

        //Spectrum index needs to be shifted by N
        for (int j = 0; j < (compit->second).second[i].second.size(); j++)
        {
          vertsSpecIdx[gPidx]
              = (unsigned int)(compit->second).second[i].first[j]
                  + totalSpecCount;
          vertsPeakMass[gPidx] = (compit->second).second[i].second[j];
          gPidx++;
        }
      }
      for (int i = 0; i < specIdx.size(); i++)
      {
        specIndices[indIdx++] = (unsigned int)specIdx[i] + totalSpecCount;
        specInfo[infoIdx++] = (unsigned short)compit->first + numComponents1;
        specInfo[infoIdx++] = (unsigned short)revIdx[i];
      }
    }
    fwrite(&specCount, sizeof(unsigned int), 1, fp);
    fwrite(specIndices, sizeof(unsigned int), specCount, fp);
    free(specIndices);
    fwrite(specInfo, sizeof(unsigned short), 2 * specCount, fp);
    free(specInfo);
    fwrite(&numComponents, sizeof(unsigned int), 1, fp);
    fwrite(vertsPerComp, sizeof(unsigned short), numComponents, fp);
    free(vertsPerComp);
    fwrite(peaksPerVert, sizeof(unsigned short), numVertices, fp);
    free(peaksPerVert);
    fwrite(vertsSpecIdx, sizeof(unsigned int), numPeaks, fp);
    free(vertsSpecIdx);
    fwrite(vertsPeakMass, sizeof(float), numPeaks, fp);
    free(vertsPeakMass);
    fclose(fp);
    return true;

  }

  void Copy_abinfo(abinfo_t& from, abinfo_t& to)
  {
    to.clear();
    for (abinfo_t::iterator abIt = from.begin(); abIt != from.end(); abIt++)
    {
      unsigned int idx = abIt->first;
      to[idx] = from[idx];
    }
  }

  int Load_abinfo_v1_0(const char *filename, abinfo_t & abinfo)
  {
    ifstream in(filename, ios::binary);

    if (!in)
      return 0;

    // read data
    uint32_t num_used_spectra;
    in.read(reinterpret_cast<char *> (&num_used_spectra),
            sizeof(num_used_spectra));

    uint32_t (* spec_indices) = new uint32_t[num_used_spectra];
    in.read(reinterpret_cast<char *> (spec_indices),
            sizeof(uint32_t[num_used_spectra]));

    uint16_t (* spec_info)[2] = new uint16_t[num_used_spectra][2];
    in.read(reinterpret_cast<char *> (spec_info),
            sizeof(uint16_t[num_used_spectra][2]));

    uint32_t num_components;
    in.read(reinterpret_cast<char *> (&num_components), sizeof(num_components));

    uint16_t (* verts_per_comp) = new uint16_t[num_components];
    in.read(reinterpret_cast<char *> (verts_per_comp),
            sizeof(uint16_t[num_components]));
    int num_verts = std::accumulate(verts_per_comp, verts_per_comp
        + num_components, 0);

    uint16_t (* peaks_per_vert) = new uint16_t[num_verts];
    in.read(reinterpret_cast<char *> (peaks_per_vert),
            sizeof(uint16_t[num_verts]));
    int num_peaks = std::accumulate(peaks_per_vert,
                                    peaks_per_vert + num_verts,
                                    0);

    uint32_t (* verts_spec_idx) = new uint32_t[num_peaks];
    in.read(reinterpret_cast<char *> (verts_spec_idx),
            sizeof(uint32_t[num_peaks]));

    float (* verts_peak_mass) = new float[num_peaks];
    in.read(reinterpret_cast<char *> (verts_peak_mass),
            sizeof(float[num_peaks]));

    // separate data into more convenient format
    for (int c = 0, cur_p = 0, cur_v = 0; c < num_components; c++)
    {
      pair<int, int> idx;
      for (idx.first = 0; idx.first < num_used_spectra; idx.first++)
        if (spec_info[idx.first][0] == c)
          break;
      for (idx.second = idx.first + 1; idx.second < num_used_spectra; idx.second++)
        if (spec_info[idx.second][0] != c)
          break;

      abinfo[c].first.first.reserve(idx.second - idx.first);
      abinfo[c].first.second.reserve(idx.second - idx.first);
      for (int i = idx.first; i < idx.second; i++)
      {
        abinfo[c].first.first.push_back(spec_indices[i]);
        abinfo[c].first.second.push_back(spec_info[i][1]);
      }

      abinfo[c].second.resize(verts_per_comp[c]);

      for (int v = 0; v < verts_per_comp[c]; v++)
      {
        pair<int, int> peak_range = make_pair(cur_p, cur_p
            + peaks_per_vert[cur_v]);
        cur_p += peaks_per_vert[cur_v];
        cur_v += 1;
        abinfo[c].second[v].first.assign(&verts_spec_idx[peak_range.first],
                                         &verts_spec_idx[peak_range.second]);
        abinfo[c].second[v].second.assign(&verts_peak_mass[peak_range.first],
                                          &verts_peak_mass[peak_range.second]);
      }
    }

    delete[] spec_indices;
    delete[] spec_info;
    delete[] verts_per_comp;
    delete[] peaks_per_vert;
    delete[] verts_spec_idx;
    delete[] verts_peak_mass;

    return 1;
  }

  int Load_abinfo_v1_1(const char *filename, abinfo_t & abinfo)
  {
    ifstream in(filename, ios::binary);

    if (!in)
      return 0;

    uint32_t type;
    in.read(reinterpret_cast<char *> (&type), sizeof(type));

    uint16_t version[2];
    in.read(reinterpret_cast<char *> (&version), sizeof(version));

    // read data
    uint32_t num_used_spectra;
    in.read(reinterpret_cast<char *> (&num_used_spectra),
            sizeof(num_used_spectra));

    uint32_t (* spec_indices)[2] = new uint32_t[num_used_spectra][2];
    in.read(reinterpret_cast<char *> (spec_indices),
            sizeof(uint32_t[num_used_spectra][2]));

    uint8_t (* spec_info) = new uint8_t[num_used_spectra];
    in.read(reinterpret_cast<char *> (spec_info),
            sizeof(uint8_t[num_used_spectra]));

    uint32_t num_components;
    in.read(reinterpret_cast<char *> (&num_components), sizeof(num_components));

    uint16_t (* verts_per_comp) = new uint16_t[num_components];
    in.read(reinterpret_cast<char *> (verts_per_comp),
            sizeof(uint16_t[num_components]));
    int num_verts = std::accumulate(verts_per_comp, verts_per_comp
        + num_components, 0);

    uint16_t (* peaks_per_vert) = new uint16_t[num_verts];
    in.read(reinterpret_cast<char *> (peaks_per_vert),
            sizeof(uint16_t[num_verts]));
    int num_peaks = std::accumulate(peaks_per_vert,
                                    peaks_per_vert + num_verts,
                                    0);

    uint32_t (* verts_spec_idx) = new uint32_t[num_peaks];
    in.read(reinterpret_cast<char *> (verts_spec_idx),
            sizeof(uint32_t[num_peaks]));

    float (* verts_peak_mass) = new float[num_peaks];
    in.read(reinterpret_cast<char *> (verts_peak_mass),
            sizeof(float[num_peaks]));

    // separate data into more convenient format
    for (int c = 0, cur_p = 0, cur_v = 0; c < num_components; c++)
    {
      pair<int, int> idx;
      for (idx.first = 0; idx.first < num_used_spectra; idx.first++)
        if (spec_indices[idx.first][1] == c)
          break;
      for (idx.second = idx.first + 1; idx.second < num_used_spectra; idx.second++)
        if (spec_indices[idx.second][1] != c)
          break;

      abinfo[c].first.first.reserve(idx.second - idx.first);
      abinfo[c].first.second.reserve(idx.second - idx.first);
      for (int i = idx.first; i < idx.second; i++)
      {
        abinfo[c].first.first.push_back(spec_indices[i][0]);
        abinfo[c].first.second.push_back(spec_info[i]);
      }

      abinfo[c].second.resize(verts_per_comp[c]);

      for (int v = 0; v < verts_per_comp[c]; v++)
      {
        pair<int, int> peak_range = make_pair(cur_p, cur_p
            + peaks_per_vert[cur_v]);
        cur_p += peaks_per_vert[cur_v];
        cur_v += 1;
        abinfo[c].second[v].first.assign(&verts_spec_idx[peak_range.first],
                                         &verts_spec_idx[peak_range.second]);
        abinfo[c].second[v].second.assign(&verts_peak_mass[peak_range.first],
                                          &verts_peak_mass[peak_range.second]);
      }
    }

    delete[] spec_indices;
    delete[] spec_info;
    delete[] verts_per_comp;
    delete[] peaks_per_vert;
    delete[] verts_spec_idx;
    delete[] verts_peak_mass;

    return 1;
  }

  bool Save_abinfo_v1_0(const char* filename, abinfo_t& abinfo)
  {
    FILE *fp = fopen(filename, "wb");
    if (fp == 0 || !fp)
    {
      cerr << "Error opening " << filename << "!!\n";
      return false;
    }
    abinfo_t::iterator compit;

    unsigned int specCount = 0, numComponents = 0, numVertices = 0, numPeaks =
        0;
    for (compit = abinfo.begin(); compit != abinfo.end(); compit++)
    {
      specCount += (compit->second).first.first.size();
      numComponents = max(numComponents, compit->first + 1);
      numVertices += (compit->second).second.size();
      vector<std::pair<vector<int> , vector<double> > > peaksPerVertex =
          (compit->second).second;
      for (int i = 0; i < peaksPerVertex.size(); i++)
      {
        numPeaks += peaksPerVertex[i].second.size();
      }
    }

    unsigned int *specIndices = (unsigned int *)malloc(specCount
        * sizeof(unsigned int)); // Vector of indices of all used spectra
    unsigned short *specInfo = (unsigned short *)malloc(2 * specCount
        * sizeof(unsigned short)); // Sequence of [component index, specFlipped]
    unsigned short *vertsPerComp = (unsigned short *)malloc(numComponents
        * sizeof(unsigned short));
    unsigned short *peaksPerVert = (unsigned short *)malloc(numVertices
        * sizeof(unsigned short));
    unsigned int *vertsSpecIdx = (unsigned int *)malloc(numPeaks
        * sizeof(unsigned int));
    float *vertsPeakMass = (float *)malloc(numPeaks * sizeof(float));

    for (unsigned int i = 0; i < numComponents; i++)
    {
      vertsPerComp[i] = 0;
    }

    unsigned int indIdx = 0, infoIdx = 0, gVidx = 0, gPidx = 0;
    for (compit = abinfo.begin(); compit != abinfo.end(); compit++)
    {
      vector<int> specIdx = (compit->second).first.first;
      vector<int> revIdx = (compit->second).first.second;
      vertsPerComp[compit->first] = (compit->second).second.size();
      for (int i = 0; i < (compit->second).second.size(); i++)
      {
        peaksPerVert[gVidx++] = (compit->second).second[i].second.size();
        for (int j = 0; j < (compit->second).second[i].second.size(); j++)
        {
          vertsSpecIdx[gPidx]
              = (unsigned int)(compit->second).second[i].first[j];
          vertsPeakMass[gPidx] = (compit->second).second[i].second[j];
          gPidx++;
        }
      }
      for (int i = 0; i < specIdx.size(); i++)
      {
        specIndices[indIdx++] = (unsigned int)specIdx[i];
        specInfo[infoIdx++] = (unsigned short)compit->first;
        specInfo[infoIdx++] = (unsigned short)revIdx[i];
      }
    }

    fwrite(&specCount, sizeof(unsigned int), 1, fp);
    fwrite(specIndices, sizeof(unsigned int), specCount, fp);
    free(specIndices);
    fwrite(specInfo, sizeof(unsigned short), 2 * specCount, fp);
    free(specInfo);
    fwrite(&numComponents, sizeof(unsigned int), 1, fp);
    fwrite(vertsPerComp, sizeof(unsigned short), numComponents, fp);
    free(vertsPerComp);
    fwrite(peaksPerVert, sizeof(unsigned short), numVertices, fp);
    free(peaksPerVert);
    fwrite(vertsSpecIdx, sizeof(unsigned int), numPeaks, fp);
    free(vertsSpecIdx);
    fwrite(vertsPeakMass, sizeof(float), numPeaks, fp);
    free(vertsPeakMass);
    fclose(fp);
    return true;
  }

  int Load_abinfo(const char *filename, abinfo_t & abinfo)
  {
    ifstream in(filename, ios::binary);

    if (!in)
      return 0;

    uint32_t type;
    in.read(reinterpret_cast<char *> (&type), sizeof(type));

    in.close();

    if (type == 0)
      return Load_abinfo_v1_1(filename, abinfo);
    else
      return Load_abinfo_v1_0(filename, abinfo);
  }

  //////////////////////////////////////////////////////////////////////////////
  int compareAbInfo(abinfo_t & abinfo1, abinfo_t & abinfo2, abDiff_t &diffData)
  {
    // initialize assuming they are equal.
    diffData.different = false;

    // get size of both abruijn graphs (by number of contigs)
    int ab1Size = abinfo1.size();
    int ab2Size = abinfo2.size();

    // check if number of contigs is the same
    if (ab1Size != ab2Size)
    {
      diffData.size1 = ab1Size;
      diffData.size2 = ab2Size;
      diffData.different = true;
    }

    // iterate thru all contigs
    abinfo_t::iterator it1, it2;

    // cycle thru
    for (it1 = abinfo1.begin(), it2 = abinfo2.begin(); it1 != abinfo1.end()
        && it2 != abinfo2.end(); it1++, it2++)
    {

      // structure to hold contig differences
      abDiffContig_t diffContig;

      // check if contig Indices are the same
      diffContig.id1 = it1->first;
      diffContig.id2 = it2->first;

      // if they are not, store the difference and move to the next
      if (diffContig.id1 != diffContig.id2)
      {
        // store data for this contig
        diffData.contigs.push_back(diffContig);
        diffData.different = true;
        continue;
      }

      // compare reversal flags

      // get data
      std::pair<vector<int> , vector<int> > &rev1 = it1->second.first;
      std::pair<vector<int> , vector<int> > &rev2 = it2->second.first;
      vector<int> &revIdx1 = rev1.first;
      vector<int> &revFlg1 = rev1.second;
      vector<int> &revIdx2 = rev2.first;
      vector<int> &revFlg2 = rev2.second;
      // get sizes
      int rev1IdxSize = revIdx1.size();
      int rev1FlgSize = revFlg1.size();
      int rev2IdxSize = revIdx2.size();
      int rev2FlgSize = revFlg2.size();
      // check if sizes are the same. If they are not, leave list empty
      if (rev1IdxSize == rev1FlgSize && rev2IdxSize == rev2FlgSize
          && rev1IdxSize == rev2IdxSize)
        // check element individually
        for (int i = 0; i < rev1IdxSize; i++)
          if (revIdx1[i] == revIdx2[i])
            if (revFlg1[i] != revFlg2[i])
            {
              diffContig.reversedList.push_back(revIdx1[i]);
              diffData.different = true;
            }

      // compare nodes

      // get the nodes
      vector<std::pair<vector<int> , vector<double> > > &nodes1 =
          it1->second.second;
      vector<std::pair<vector<int> , vector<double> > > &nodes2 =
          it2->second.second;

      // get the nodes size
      diffContig.nodesList.size1 = nodes1.size();
      diffContig.nodesList.size2 = nodes2.size();

      // check if nodes lists have the same size
      if (diffContig.nodesList.size1 != diffContig.nodesList.size2)
      {
        diffData.different = true;
        continue;
      }

      // cycle thru the nodes
      for (int i = 0; i < diffContig.nodesList.size1; i++)
      {

        abDiffNode_t node;

        // set node index
        node.idx = i;

        // get the nodes ids and data
        vector<int> &starIdx1 = nodes1[i].first;
        vector<double> &starData1 = nodes1[i].second;
        vector<int> &starIdx2 = nodes2[i].first;
        vector<double> &starData2 = nodes2[i].second;

        // get nodes sizes
        node.size1i = starIdx1.size();
        node.size1d = starData1.size();
        node.size2i = starIdx2.size();
        node.size2d = starData2.size();

        // if nodes have different sizes on ids...
        if (node.size1i != node.size2i || node.size1i != node.size1d
            || node.size2i != node.size2d)
        {
          diffContig.nodesList.nodeData.push_back(node);
          diffData.different = true;
          continue;
        }

        // cycle thru each node's idx and data
        for (int j = 0; j < node.size1i; j++)
        {

          abDiffPeak_t peak;
          peak.starId1 = starIdx1[j];
          peak.starId2 = starIdx2[j];

          // check star ids in the node
          if (starIdx1[j] != starIdx2[j])
          {
            // store peak info
            node.peaks.push_back(peak);
            diffData.different = true;
          }
          else

          // check data in the node
          if (starData1[j] != starData2[j])
          {

            diffData.different = true;
            peak.val1 = starData1[j];
            peak.val2 = starData2[j];
            // store peak info
            node.peaks.push_back(peak);
          }

        }
        // store the node list
        diffContig.nodesList.nodeData.push_back(node);
      }
      // store data for this contig
      diffData.contigs.push_back(diffContig);
    }
    return diffData.different;
  }

  //////////////////////////////////////////////////////////////////////////////
  int dumpAbInfo(const char * filename, abinfo_t & abinfo)
  {
    FILE *fp = fopen(filename, "wb");
    if (fp == 0 || !fp)
    {
      cerr << "Error opening " << filename << "!!\n";
      return false;
    }

    abinfo_t::iterator compit;

    unsigned int specCount = 0, numComponents = 0, numVertices = 0, numPeaks =
        0;
    vector<std::pair<vector<int> , vector<double> > > peaksPerVertex1;
    ;

    cout << "Dumping abinfo to " << filename << endl;

    //Count the info from the first abinfo
    for (compit = abinfo.begin(); compit != abinfo.end(); compit++)
    {
      specCount += (compit->second).first.first.size();
      numComponents++;
      numVertices += (compit->second).second.size();
      peaksPerVertex1 = (compit->second).second;
      for (int i = 0; i < peaksPerVertex1.size(); i++)
      {
        numPeaks += peaksPerVertex1[i].second.size();
      }
    }

    unsigned int *specIndices = (unsigned int *)malloc(specCount
        * sizeof(unsigned int)); // Vector of indices of all used spectra
    unsigned short *specInfo = (unsigned short *)malloc(2 * specCount
        * sizeof(unsigned short)); // Sequence of [component index, specFlipped]
    unsigned short *vertsPerComp = (unsigned short *)malloc(numComponents
        * sizeof(unsigned short));
    unsigned short *peaksPerVert = (unsigned short *)malloc(numVertices
        * sizeof(unsigned short));
    unsigned int *vertsSpecIdx = (unsigned int *)malloc(numPeaks
        * sizeof(unsigned int));
    float *vertsPeakMass = (float *)malloc(numPeaks * sizeof(float));

    unsigned int indIdx = 0, infoIdx = 0, gVidx = 0, gPidx = 0;
    for (compit = abinfo.begin(); compit != abinfo.end(); compit++)
    {

      //specIndexes for the first info do not need to be updated
      vector<int> specIdx = (compit->second).first.first;
      vector<int> revIdx = (compit->second).first.second;

      //contig indices for the first info do not need to be updated
      vertsPerComp[compit->first] = (compit->second).second.size();

      for (int i = 0; i < (compit->second).second.size(); i++)
      {
        peaksPerVert[gVidx++] = (compit->second).second[i].second.size();

        for (int j = 0; j < (compit->second).second[i].second.size(); j++)
        {
          vertsSpecIdx[gPidx]
              = (unsigned int)(compit->second).second[i].first[j];
          vertsPeakMass[gPidx] = (compit->second).second[i].second[j];
          gPidx++;
        }
      }
      for (int i = 0; i < specIdx.size(); i++)
      {
        specIndices[indIdx++] = (unsigned int)specIdx[i];
        specInfo[infoIdx++] = (unsigned short)compit->first;
        specInfo[infoIdx++] = (unsigned short)revIdx[i];
      }
    }
    string line = "Spec count: ";
    line += intToString(specCount);
    line += "\n";

    fwrite(line.c_str(), sizeof(char), line.size(), fp);
    int k = 0;
    for (int i = 0; i < specCount; ++i)
    {
      line = "spectrum ";
      line += intToString(specIndices[i]);
      line += " ";
      line += intToString(specInfo[k++]);
      line += " ";
      line += intToString(specInfo[k++]);
      line += "\n";

      fwrite(line.c_str(), sizeof(char), line.size(), fp);
    }
    line = "Num components: ";
    line += intToString(numComponents);
    line += "\n";
    fwrite(line.c_str(), sizeof(char), line.size(), fp);
    for (k = 0; k < numComponents; k++)
    {
      line = "component ";
      line += intToString(k);
      line += " has ";
      line += intToString(vertsPerComp[k]);
      line += " vertices\n";
      fwrite(line.c_str(), sizeof(char), line.size(), fp);
    }

    for (k = 0; k < numVertices; ++k)
    {
      line = "vertex ";
      line += intToString(k);
      line += " has ";
      line += intToString(peaksPerVert[k]);
      line += " peaks\n";
      fwrite(line.c_str(), sizeof(char), line.size(), fp);
    }

    for (k = 0; k < numPeaks; ++k)
    {
      line = "peak ";
      line += intToString(k);
      line += " is from spectrum ";
      line += intToString(vertsSpecIdx[k]);
      line += " mass ";
      line += intToString(vertsPeakMass[k]);
      line += "\n";
      fwrite(line.c_str(), sizeof(char), line.size(), fp);
    }
    fclose(fp);
  }
}
