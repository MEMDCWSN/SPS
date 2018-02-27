#include "Logger.h"
#include "Filters.h"

namespace specnets
{
  /**
   * TODO: add description
   *
   *@param edge1
   *@param edge2
   *@param edge3
   *@param pmTol
   *@return
   */
  inline bool validateTriangle(SpectrumPair &edge1,
                               SpectrumPair &edge2,
                               SpectrumPair &edge3,
                               float pmTol)
  {

    /* Considers 3 possibilities for c1: edge1, edge2, edge3:
     *          c1
     *   O---------------->O
     *    \------>O------>/
     *       c2      c3
     *  Note that one of the 3 possibilities _has_ to correspond to c1.
     */
    SpectrumPair *short1, *short2; // Pointers to the 2 short edges (c2,c3 above)
    float longShift, // Value of the long shift (c1 above)
        curShift;

    // Find longShift, short1, short2
    longShift = abs(edge1.shift1) > abs(edge1.shift2) ? edge1.shift1
        : edge1.shift2;
    short1 = &edge2;
    short2 = &edge3;
    curShift = abs(edge2.shift1) > abs(edge2.shift2) ? edge2.shift1
        : edge2.shift2;
    if (abs(curShift) > abs(longShift))
    {
      longShift = curShift;
      short1 = &edge1;
      short2 = &edge3;
    }
    curShift = abs(edge3.shift1) > abs(edge3.shift2) ? edge3.shift1
        : edge3.shift2;
    if (abs(curShift) > abs(longShift))
    {
      longShift = curShift;
      short1 = &edge1;
      short2 = &edge2;
    }

    // Test c1 = c2+c3
    if (abs(longShift - short1->shift1 - short2->shift1) <= pmTol
        or abs(longShift - short1->shift1 - short2->shift2) <= pmTol
        or abs(longShift - short1->shift2 - short2->shift1) <= pmTol
        or abs(longShift - short1->shift2 - short2->shift2) <= pmTol)
    {
      return true;
    }

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
  unsigned int filterTriangles(SpectrumPairSet &aligns,
                               unsigned int idxStart,
                               unsigned int idxEnd,
                               float pmTol,
                               vector<unsigned int> &selectedIdx)
  {
    vector<bool> valid(aligns.size()); // Indicates whether an edge is part of some triangle
    vector<TwoValues<unsigned int> > adj; // Start/end indices (in aligns) for each vertex's adjacency list
    unsigned int maxVertex = 0, // Largest vertex index in aligns
        idxMain, // Index of the edge being processed
        idxPair; // Index of the tentative pair for idxMain

    // Initialize adjacency lists
    for (idxMain = 0; idxMain < aligns.size(); idxMain++)
    {
      valid[idxMain] = false;
      if ((unsigned int)aligns[idxMain].spec1 > maxVertex)
      {
        maxVertex = aligns[idxMain].spec1;
      }
      if ((unsigned int)aligns[idxMain].spec2 > maxVertex)
      {
        maxVertex = aligns[idxMain].spec2;
      }
    }
    adj.resize(maxVertex + 1);
    for (unsigned int vertexIdx = 0; vertexIdx <= maxVertex; vertexIdx++)
    {
      adj[vertexIdx].set(1, 0); // Adjacency list is empty if start>end
    }
    for (idxMain = 0; idxMain < aligns.size(); idxMain++)
    {
      if (idxMain == 0 or aligns[idxMain - 1].spec1 != aligns[idxMain].spec1)
      {
        adj[aligns[idxMain].spec1][0] = idxMain;
      }
      adj[aligns[idxMain].spec1][1] = idxMain;
    }

    //cerr << "First edge ("<<aligns[0].spec1<<","<<aligns[0].spec2<<")\n";

    // Find triangles - for each edge (v1,v2), find a pairing edge (v1,v3) from the same vertex
    //  and look for a connector edge between the two destination vertices (v2,v3).
    //  pmTol *= 2;  // *2 = Once for c2 and once for c3
    TwoValues<unsigned int> searchBounds; // lower/upper limits for the binary search for a connector edge
    unsigned int idxMiddle; // index used in the binary search - middle index between searchBounds[0]/searchBounds[1]
    unsigned int targetVertex; // vertex to find in the adjacency list
    idxEnd = min(idxEnd, (unsigned int)aligns.size());
    idxEnd = max((unsigned int)0, idxEnd);
    idxStart = max(idxStart, (unsigned int)0);

    //cerr<<"Looking for triangles with an edge between indices "<<idxStart<<" and "<<idxEnd<<"\n"; cerr.flush();
    for (idxMain = idxStart; idxMain < idxEnd; idxMain++)
    {
      for (idxPair = adj[aligns[idxMain].spec1][0]; idxPair
          < adj[aligns[idxMain].spec1][1]; idxPair++)
      {
        //cerr << "Edge indices ("<<idxMain<<","<<idxPair<<"), edges ("<<aligns[idxMain].spec1<<"->"<<aligns[idxMain].spec2<<"), ("<<aligns[idxPair].spec1<<"->"<<aligns[idxPair].spec2<<")\n";
        if (idxPair == idxMain)
        {
          continue;
        }
        if (idxPair < idxMain)
        {
          searchBounds = adj[aligns[idxPair].spec2];
          targetVertex = aligns[idxMain].spec2;
        }
        else
        {
          searchBounds = adj[aligns[idxMain].spec2];
          targetVertex = aligns[idxPair].spec2;
        }

        //cerr << "Looking for targetVertex "<<aligns[idxPair].spec2<<" between indices "<<searchBounds[0]<<":"<<searchBounds[1]<<endl;
        // Binary search for targetVertex in the adjacency defined by searchBounds
        while (searchBounds[0] <= searchBounds[1])
        {
          //cerr << " --- searchBounds = ("<<searchBounds[0]<<","<<searchBounds[1]<<")\n";
          if ((unsigned int)aligns[searchBounds[1]].spec2 == targetVertex)
          {
            searchBounds[0] = searchBounds[1];
          }
          if ((unsigned int)aligns[searchBounds[0]].spec2 == targetVertex
              or searchBounds[1] - searchBounds[0] < 2)
          {
            break;
          }
          idxMiddle = max(searchBounds[0] + 1,
                          (unsigned int)round((searchBounds[0]
                              + searchBounds[1]) / 2.0));
          if ((unsigned int)aligns[idxMiddle].spec2 <= targetVertex)
          {
            searchBounds[0] = idxMiddle;
          }
          else
          {
            searchBounds[1] = idxMiddle;
          }
        }
        if ((unsigned int)aligns[searchBounds[0]].spec2 != targetVertex)
        {
          continue; // Connector edge not found
        }

        //cerr<<"Found targetVertex at "<<searchBounds[0]<<", edge ("<<aligns[searchBounds[0]].spec1<<","<<aligns[searchBounds[0]].spec2<<")\n";

        if (validateTriangle(aligns[idxMain],
                             aligns[idxPair],
                             aligns[searchBounds[0]],
                             pmTol))
        {
          valid[idxMain] = true;
          valid[idxPair] = true;
          valid[searchBounds[0]] = true;
          //cerr<<"Validated triangle: idxMain="<<idxMain<<", idxPair="<<idxPair<<", searchBounds[0]="<<searchBounds[0]<<"\n"; cerr.flush();
        }
      }

      if ((idxMain % 1000) == 0)
      {
        DEBUG_MSG("Done processing edge " << idxMain << " (" << ((double)idxMain
            - idxStart) / ((double)idxEnd - idxStart) << "% completed)");
      }
    }

    // Remove edges that do not participate in a triangle
    unsigned int idxKept = 0; // Index to the last kept alignment
    for (idxMain = 0; idxMain < aligns.size(); idxMain++)
    {
      if (valid[idxMain])
      {
        aligns[idxKept++] = aligns[idxMain];
      }
    }
    aligns.resize(idxKept);
    selectedIdx.resize(idxKept);
    idxKept = 0;
    for (idxMain = 0; idxMain < valid.size(); idxMain++)
    {
      if (valid[idxMain])
      {
        selectedIdx[idxKept++] = idxMain;
      }
    }
    return idxKept; // Number of retained alignments
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
  unsigned int filterAligns(SpectrumPairSet &aligns,
                            vector<unsigned int> &idxKept,
                            vector<TwoValues<float> > &pvalues,
                            vector<TwoValues<float> > &means,
                            vector<float> &variances,
                            vector<TwoValues<float> > &ratios,
                            float minPValue,
                            float minRatio,
                            float pmTol,
                            bool filterTrigs)
  {
    // Compute pvalues and filter by ratios and p-values
    unsigned int idxPair, idxLast = 0;
    idxKept.resize(aligns.size());
    pvalues.resize(aligns.size());
    for (idxPair = 0; idxPair < aligns.size(); idxPair++)
    {
      pvalues[idxPair][0] = 1
          - Utils::gaussiancdf(aligns[idxPair].score1,
                               means[aligns[idxPair].spec1][0],
                               variances[aligns[idxPair].spec1]);
      pvalues[idxPair][1] = 1
          - Utils::gaussiancdf(aligns[idxPair].score2,
                               means[aligns[idxPair].spec2][0],
                               variances[aligns[idxPair].spec2]);

      if (pvalues[idxPair][0] <= minPValue && pvalues[idxPair][1] <= minPValue
          && (ratios.size() == 0 || (ratios[idxPair][0] >= minRatio
              && ratios[idxPair][1] >= minRatio)))
      {
        if (idxLast < idxPair)
        {
          aligns[idxLast] = aligns[idxPair];
        }
        idxKept[idxLast] = idxPair;
        idxLast++;
      }
    }

    DEBUG_VAR(idxLast);
    aligns.resize(idxLast);
    idxKept.resize(idxLast);

    vector<unsigned int> idxKeptT;
    if (filterTrigs)
    {
      filterTriangles(aligns, 0, aligns.size(), pmTol, idxKeptT);
      idxLast = 0;
      for (unsigned int idxPair = 0; idxPair < idxKeptT.size(); idxPair++)
      {
        if (idxLast < idxPair)
        {
          aligns[idxLast] = aligns[idxPair];
          idxKept[idxLast] = idxKeptT[idxPair];
        }
        idxLast++;
      }
      aligns.resize(idxLast);
      idxKept.resize(idxLast);
    }
    return idxKept.size();
  }

  void filterMatchRatioASP(SpecSet & specSet,
                           SpectrumPairSet & results,
                           short ratioType,
                           vector<TwoValues<float> > & ratios)
  {
    float szOverlap, shift, ratio1, ratio2, totalScore1, totalScore2;
    int spec1, spec2, i, j;

    ratios.resize(results.size());
    for (i = 0; i < results.size(); i++)
    {
      spec1 = results[i].spec1;
      spec2 = results[i].spec2;
      shift = results[i].shift1;
      szOverlap = min(specSet[spec1].parentMass, specSet[spec2].parentMass);

      totalScore1 = 0;
      for (j = 0; j <= specSet[spec1].size() && specSet[spec1][j][0]
          <= szOverlap - 56; j++)
      {
        if (specSet[spec1][j][0] > 56)
        {
          totalScore1 += specSet[spec1][j][1];
        }
      }
      while (j <= specSet[spec1].size() && specSet[spec1][j][0]
          <= specSet[spec1].parentMass - szOverlap + 56)
      {
        j++;
      }
      for (; j <= specSet[spec1].size() && specSet[spec1][j][0]
          <= specSet[spec1].parentMass - 56; j++)
      {
        totalScore1 += specSet[spec1][j][1];
      }
      totalScore2 = 0;
      for (j = 0; j <= specSet[spec2].size() && specSet[spec2][j][0]
          <= szOverlap - 56; j++)
      {
        if (specSet[spec2][j][0] > 56)
        {
          totalScore2 += specSet[spec2][j][1];
        }
      }
      while (j <= specSet[spec2].size() && specSet[spec2][j][0]
          <= specSet[spec2].parentMass - szOverlap + 56)
      {
        j++;
      }
      for (; j <= specSet[spec2].size() && specSet[spec2][j][0]
          <= specSet[spec2].parentMass - 56; j++)
      {
        totalScore2 += specSet[spec2][j][1];
      }

      //cerr << "score1 = " << results[i].score1 << ", score2 = " << results[i].score2 << endl;
      //cerr << "totalScore1 = " << totalScore1 << ", totalScore2 = " << totalScore2 << endl;

      switch (ratioType)
      {
      case 6:
        if (totalScore1 - results[i].score1 > 1)
        {
          ratio1 = results[i].score1 / (totalScore1 - results[i].score1);
        }
        else
        {
          ratio1 = 1;
        }
        if (totalScore2 - results[i].score2 > 1)
        {
          ratio2 = results[i].score2 / (totalScore2 - results[i].score2);
        }
        else
        {
          ratio2 = 1;
        }
        break;
      case 7:
        ratio1 = results[i].score1 / totalScore1;
        ratio2 = results[i].score2 / totalScore2;
        break;
      default:
        ratio1 = 0;
        ratio2 = 0;
      }

      if (ratio1 < 0)
      {
        ratios[i][0] = 0;
      }
      else
      {
        ratios[i][0] = ratio1;
      }
      if (ratio2 < 0)
      {
        ratios[i][1] = 0;
      }
      else
      {
        ratios[i][1] = ratio2;
      }
      //cerr << " --- ratio 1 = " << ratios[i][0] << ", ratio 2 = " << ratios[i][1] << endl;
    }
  }

  void computeScoreMeans(SpecSet &specSet, SpectrumPairSet &results, vector<
      TwoValues<float> > &scoreMeans)
  {
    int i, specIdx, n, numSpecs = specSet.size();

    scoreMeans.resize(numSpecs);
    for (i = 0; i < numSpecs; i++)
    {
      scoreMeans[i][0] = 0;
      scoreMeans[i][1] = 0;
    }

    for (i = 0; i < results.size(); i++)
    {
      specIdx = results[i].spec1;
      n = (int)scoreMeans[specIdx][1];
      scoreMeans[specIdx][0] = scoreMeans[specIdx][0] * n / (n + 1)
          + results[i].score1 / (n + 1);
      scoreMeans[specIdx][1] = n + 1;

      specIdx = results[i].spec2;
      n = (int)scoreMeans[specIdx][1];
      scoreMeans[specIdx][0] = scoreMeans[specIdx][0] * n / (n + 1)
          + results[i].score2 / (n + 1);
      scoreMeans[specIdx][1] = n + 1;
    }
  }

} //namespace specnets
