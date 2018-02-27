// Header Includes
#include "ExecAssembly.h"
#include "ExecParallelAssembly.h"

// Module Includes
#include "Logger.h"
//#include "FileUtils.h"

// SpecNets Includes

using namespace std;
using namespace specnets;

namespace specnets
{
  void SplitASPPA(SpectrumPairSet &pairs,
                  float pmTol,
                  float maxModMass,
                  SpectrumPairSet &pairsASP,
                  SpectrumPairSet &pairsPA)
  {
    pairsASP.resize(pairs.size());
    pairsPA.resize(pairs.size());
    unsigned int idxPair, idxASP = 0, idxPA = 0;
    for (idxPair = 0; idxPair < pairs.size(); idxPair++)
      if ((fabs(pairs[idxPair].shift1) <= pmTol and fabs(pairs[idxPair].shift2)
          <= maxModMass) or (fabs(pairs[idxPair].shift1) <= maxModMass
          and fabs(pairs[idxPair].shift2) <= pmTol))
        pairsASP[idxASP++] = pairs[idxPair];
      else
        pairsPA[idxPA++] = pairs[idxPair];
    pairsASP.resize(idxASP);
    pairsPA.resize(idxPA);
  }

  // -------------------------------------------------------------------------

  ExecAssembly::ExecAssembly(void) :
    m_spectra(0x0), m_spectrumPairs(0x0), ownInput(true), m_contigShifts(0x0),
        m_outputAbruijn(0x0), ownOutput(true)
  {
    m_name = "ExecAssembly";
    m_type = "ExecAssembly";
    m_labels.resize(0);
  }

  // -------------------------------------------------------------------------

  ExecAssembly::ExecAssembly(const ParameterList & inputParams) :
    ExecBase(inputParams), m_spectra(0x0), m_spectrumPairs(0x0), ownInput(true),
        m_contigShifts(0x0), m_outputAbruijn(0x0), ownOutput(true)
  {
    m_name = "ExecAssembly";
    m_type = "ExecAssembly";
    m_labels.resize(0);
  }

  // -------------------------------------------------------------------------

  ExecAssembly::ExecAssembly(const ParameterList & inputParams,
                             SpecSet * spectra,
                             SpectrumPairSet * spectrumPairs,
                             Clusters * outputContigShifts,
                             abinfo_t * outputAbruijn) :
    ExecBase(inputParams), m_spectra(spectra), m_spectrumPairs(spectrumPairs),
        ownInput(false), m_contigShifts(outputContigShifts), ownOutput(false),
        m_outputAbruijn(outputAbruijn)
  {
    m_name = "ExecAssembly";
    m_type = "ExecAssembly";
    m_labels.resize(0);
  }

  // -------------------------------------------------------------------------

  ExecAssembly::~ExecAssembly(void)
  {
    if (ownInput)
    {
      if (m_spectra)
        delete m_spectra;
      if (m_spectrumPairs)
        delete m_spectrumPairs;
    }
    if (ownOutput)
    {
      if (m_contigShifts)
        delete m_contigShifts;
      if (m_outputAbruijn)
        delete m_outputAbruijn;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecAssembly::clone(const ParameterList & inputParams) const
  {
    return new ExecAssembly(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecAssembly::invoke(void)
  {
    DEBUG_TRACE;
    if (m_params.getValueBool("DEBUG_PARAMS")) {
      stringstream aux;
      m_params.print(aux);
      DEBUG_MSG(aux.str());
    }

    float penalty_ptm = (float)m_params.getValueDouble("PENALTY_PTM");
    float penalty_sameVert =
        (float)m_params.getValueDouble("PENALTY_SAME_VERTEX");

    int graphType = m_params.exists("GRAPH_TYPE")
        ? m_params.getValueInt("GRAPH_TYPE") : 2;
    int maxAAjump = m_params.exists("MAX_AA_JUMP")
        ? m_params.getValueInt("MAX_AA_JUMP") : 0;
    float maxModMass = m_params.exists("MAX_MOD_MASS")
        ? (float)m_params.getValueDouble("MAX_MOD_MASS") : 100.0;
    float peakTol = m_params.exists("TOLERANCE_PEAK")
        ? (float)m_params.getValueDouble("TOLERANCE_PEAK") : 0.5;
    float pmTol = m_params.exists("TOLERANCE_PM")
        ? (float)m_params.getValueDouble("TOLERANCE_PM") : 1.0;
    short edgeScoreType = (short)(m_params.exists("EDGE_SCORE_TYPE")
        ? m_params.getValueInt("EDGE_SCORE_TYPE") : 1);
    unsigned int minMatchedPeaks = m_params.exists("MIN_MATCHED_PEAKS")
        ? (int)m_params.getValueInt("MIN_MATCHED_PEAKS") : 1;
    unsigned int minEdgesToComponent =
        m_params.exists("MIN_EDGES_TO_COMPONENT")
            ? (int)m_params.getValueInt("MIN_EDGES_TO_COMPONENT") : 1;
    unsigned int pathMinSpecs =
        (unsigned int)(m_params.exists("PATH_MIN_SPECS")
            ? m_params.getValueInt("PATH_MIN_SPECS") : 1);
    short pathMinPeaks = (short)(m_params.exists("PATH_MIN_PEAKS")
        ? m_params.getValueInt("PATH_MIN_PEAKS") : 1);
    int specType = m_params.exists("SPEC_TYPE_MSMS")
        ? ((int)m_params.getValueInt("SPEC_TYPE_MSMS") ? 1 : 0) : 0;
    float ionOffset = specType ? AAJumps::massHion : 0;

    bool noSequencing = m_params.getValueBool("NO_SEQUENCING", false);
    bool addEndpoints = m_params.getValueBool("ADD_ENDPOINTS", true);
    bool parallelPaths = m_params.getValueBool("PARALLEL_PATHS", false);
    const string wholeABFN = m_params.getValue("OUTPUT_COMPLETE_ABRUIJN");

    bool ignoreReversals = m_params.getValueBool("IGNORE_REVERSALS", 0);

    if (!m_spectra or m_spectra->size() == 0)
    {
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

    if (!m_spectrumPairs or m_spectrumPairs->size() == 0)
    {
      ERROR_MSG("ERROR: empty set of input spectral pairs");
      return false;
    }

    if (addEndpoints)
    {
      // Make sure every spectrum has zero/parentMass-19 nodes
      if (!m_labels.empty())
      {
        for (unsigned int i = 0; i < m_spectra->size(); i++)
          (*m_spectra)[i].addZPMpeaks(peakTol, ionOffset, true, &m_labels[i]);
      }
      else
      {
        for (unsigned int i = 0; i < m_spectra->size(); i++)
          (*m_spectra)[i].addZPMpeaks(peakTol, ionOffset, true);
      }
    }

    // Separate aligns into components
    SpectrumPairSet tmpAlignsPA, tmpAlignsASP; // Quick-fix variables used below
    SetMerger components(m_spectra->size());
    components.createSets(m_spectra->size(), 2, tmpAlignsASP, *m_spectrumPairs);
    components.splitAligns(tmpAlignsASP, *m_spectrumPairs);

    DEBUG_MSG("Got "<<components.size()<<" component(s).");

    vector<MSGraph> spectrumGraphs(m_spectra->size()); // Used to build an MSGraph for each spectrum
    AAJumps jumps2(2);

    // Split aligns per component and process each component
    char sBuf[1024];
    AAJumps jumps(1);
    unsigned int numElemsInSets = components.numElemsInSets(); // Maximum possible number of components (all singletons)
    vector<vector<float> > cStats(numElemsInSets);
    for (unsigned int cIdx = 0; cIdx < cStats.size(); cIdx++)
    {
      cStats[cIdx].resize(9);
      for (unsigned int i = 0; i < 9; i++)
        cStats[cIdx][i] = 0;
    } // Absolute value 9 depends on MSGraph::info_heaviestPath
    m_contigShifts->resize(numElemsInSets); // Keep the de-novo reconstructed heaviestPath sequences as spectra in a Cluster variable

    vector<bool> specFlipped(m_spectra->size());
    for (unsigned int i = 0; i < specFlipped.size(); i++)
    {
      specFlipped[i] = false;
      (*m_spectra)[i].setPeakTolerance(peakTol);
    }

    // Keeps track of which spectrum peaks were matched (dim.3) in each ABruijn vertex (dim.2)
    //   from the de novo sequence (heaviest path) in each component (dim.1)
    vector<vector<list<TwoValues<int> > > > abVertices(numElemsInSets);
    for (unsigned int i = 0; i < abVertices.size(); i++)
      abVertices[i].resize(0);

    // Keeps track of which spectrum peaks were matched (dim.3) for all
    //   ABruijn vertices (dim.2) in each component (dim.1)
    vector<vector<list<TwoValues<int> > > > abVerticesAll(numElemsInSets);
    for (unsigned int i = 0; i < abVertices.size(); i++)
      abVerticesAll[i].resize(0);

    vector<vector<short> > abCounts(numElemsInSets); // Records info on number of vertices and edges per ABruijn graph
    for (unsigned int i = 0; i < abCounts.size(); i++)
    {
      abCounts[i].resize(2);
      abCounts[i][0] = 0;
      abCounts[i][1] = 0;
    }

    //	for(unsigned int cIdx=20; cIdx<21; cIdx++) {
    unsigned int prevNumSpecs;
    for (unsigned int cIdx = 0; cIdx < components.size(); cIdx++)
    {

      DEBUG_MSG("Processing component "<<cIdx);

      prevNumSpecs = 0;
      while (prevNumSpecs != components.sets[cIdx].size())
      { // Iterates ABruijn/sequencing/split until all
        //   remaining spectra match the best ABruijn path
        prevNumSpecs = components.sets[cIdx].size();

        stringstream logMsg;
        logMsg << "  - Spectrum indices: ";
        for (list<int>::iterator iter = components.sets[cIdx].begin(); iter
            != components.sets[cIdx].end(); iter++)
          logMsg << *iter << " ";
        logMsg << endl;
        logMsg << "  - Component defined by "
            << components.cAlignsPA[cIdx].size() << " pairs...\n";
        DEBUG_MSG(logMsg.str());

        vector<vector<TwoValues<int> > > matches;

        //
        // Choose consensus orientations for pairwise alignments - a spectrum can only be
        //   used as-is or reversed, not both. Also determines set of matched peaks for
        //   the consensus orientations.
        //
        vector<float> modPos;
        //		SplitPairs3(*m_spectra, cAlignsASP[cIdx], components.cAlignsPA[cIdx], peakTol, maxAAjump, penalty_sameVert, penalty_ptm, matches, matchesPA, specFlipped, modPos, false, &m_labels);
        SplitPairs(*m_spectra,
                   components.cAlignsPA[cIdx],
                   peakTol,
                   pmTol,
                   maxAAjump,
                   maxModMass,
                   penalty_sameVert,
                   penalty_ptm,
                   matches,
                   specFlipped,
                   modPos,
                   minMatchedPeaks,
                   minEdgesToComponent,
                   false,
                   ignoreReversals,
                   NULL,
                   &m_labels);

        /*if(cIdx==62)
         for(unsigned int i=0; i<components.cAlignsPA[cIdx].size(); i++) {
         stringstream logMsg; logMsg<<"*** Pair ("<<components.cAlignsPA[cIdx][i].spec1<<","<<components.cAlignsPA[cIdx][i].spec2<<"), shifts "<<components.cAlignsPA[cIdx][i].shift1<<"/"<<components.cAlignsPA[cIdx][i].shift2<<": ";
         for(unsigned int j=0; j<matches[i].size(); j++) logMsg<<"("<<matches[i][j][0]<<","<<matches[i][j][1]<<")";
         DEBUG_MSG(logMsg.str());
         }*/

        SplitASPPA(components.cAlignsPA[cIdx],
                   pmTol,
                   maxModMass,
                   tmpAlignsASP,
                   tmpAlignsPA);

        /*			vector<bool> present(m_spectra->size()); for(unsigned int i=0; i<present.size();i++) present[i]=false;
         for(unsigned int i=0;i<components.cAlignsPA[cIdx].size();i++)
         if(matches[i].size()>0) {
         if(!present[components.cAlignsPA[cIdx][i].spec1]) { cSpectra[cIdx].push_back(components.cAlignsPA[cIdx][i].spec1); present[components.cAlignsPA[cIdx][i].spec1]=true; }
         if(!present[components.cAlignsPA[cIdx][i].spec2]) { cSpectra[cIdx].push_back(components.cAlignsPA[cIdx][i].spec2); present[components.cAlignsPA[cIdx][i].spec2]=true; }
         }
         */
        //
        // Maximize endpoint scores for the matched spectra
        //
        /*
         int debugIdx = 12;

         if (cIdx == debugIdx)
         {
         DEBUG_VAR(cIdx);
         }
         */
        list<int>::iterator sicIter;
        if (addEndpoints)
        {
          for (sicIter = components.sets[cIdx].begin(); sicIter
              != components.sets[cIdx].end(); sicIter++)
          {

            (*m_spectra)[*sicIter].maximizeZPMpeaks(peakTol, 0, true);
          }
        }

#ifdef DBG_MASAB
        // Add labels to the split pairs graph for graphviz output
        MSGraph g; g.build(tmpAlignsASP); g.add(tmpAlignsPA);
        g.vLabels.resize(m_spectra->size());
        if(m_spectra->size()==m_spectra->size())
        for(unsigned int i=0; i<m_spectra->size(); i++)
        {
          if(specFlipped[i]) sprintf(sBuf,"v%d_R",i);
          else sprintf(sBuf,"v%d",i);
          g.vLabels[i] = string((const char *)sBuf);
        }
        else for(unsigned int i=0; i<m_spectra->size(); i++)
        {
          sprintf(sBuf,"v%d",i); g.vLabels[2*i] = string((const char *)sBuf);
          sprintf(sBuf,"v%d_R",i); g.vLabels[2*i+1] = string((const char *)sBuf);
        }
        sprintf(sBuf,"split_pairs_graph_%d.txt",cIdx); g.output_graphviz(sBuf);
#endif

        //
        // Build spectrum graphs for the spectra in this component
        //
        if (graphType > 0)
          for (unsigned int i = 0; i < components.cAlignsPA[cIdx].size(); i++)
          {
            if (spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].numVerts()
                == 0)
            {
              if (graphType == 1)
                spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].ConnectConsecutive((*m_spectra)[components.cAlignsPA[cIdx][i].spec1]);
              else
                spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].ConnectJumps((*m_spectra)[components.cAlignsPA[cIdx][i].spec1],
                                                                                 jumps2,
                                                                                 peakTol);
              //sprintf(sBuf,"spectrum_graph_%d.txt",components.cAlignsPA[cIdx][i].spec1);
              //spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].output_graphviz(sBuf);
            }
            if (spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].numVerts()
                == 0)
            {
              if (graphType == 1)
                spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].ConnectConsecutive((*m_spectra)[components.cAlignsPA[cIdx][i].spec2]);
              else
                spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].ConnectJumps((*m_spectra)[components.cAlignsPA[cIdx][i].spec2],
                                                                                 jumps2,
                                                                                 peakTol);
              //sprintf(sBuf,"spectrum_graph_%d.txt",components.cAlignsPA[cIdx][i].spec2);
              //spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].output_graphviz(sBuf);
            }
          }

        //
        // Build A-Bruijn graph
        //
        VertexSet vSet(*m_spectra, 2048);
        vector<bool> usedSpectra(m_spectra->size());
        for (unsigned int i = 0; i < m_spectra->size(); i++)
          usedSpectra[i] = false;

        for (unsigned int i = 0; i < matches.size(); i++)
        {

          if (matches[i].size() > 0)
          {
            float
                specShift =
                    (*m_spectra)[components.cAlignsPA[cIdx][i].spec2][matches[i][0][1]][0]
                        - (*m_spectra)[components.cAlignsPA[cIdx][i].spec1][matches[i][0][0]][0];

            if (MZRange::EqualWithinRange(specShift, 0, peakTol * 2)
                && addEndpoints && matches[i][0][0] != 0 && matches[i][0][1]
                != 0)
            {
              DEBUG_MSG("Adding B end-point glues for (" << components.cAlignsPA[cIdx][i].spec1 << "," << components.cAlignsPA[cIdx][i].spec2 << ")");
              vector<TwoValues<int> > tempMatches;
              matches[i].swap(tempMatches);
              matches[i].push_back(TwoValues<int> (0, 0));
              matches[i].insert(matches[i].end(),
                                tempMatches.begin(),
                                tempMatches.end());
            }

            if (graphType == 0) // adds edges for all consecutive matched spectrum peaks
            {

              vSet.addGlues(components.cAlignsPA[cIdx][i].spec1,
                            components.cAlignsPA[cIdx][i].spec2,
                            matches[i]);
            }
            else
            {
              // add edges later
              vSet.addGlues(components.cAlignsPA[cIdx][i].spec1,
                            components.cAlignsPA[cIdx][i].spec2,
                            matches[i],
                            &spectrumGraphs);
            }
            usedSpectra[components.cAlignsPA[cIdx][i].spec1] = true;
            usedSpectra[components.cAlignsPA[cIdx][i].spec2] = true;
          }
        }
        list<int> includedSpecs;
        includedSpecs.clear();
        for (unsigned int i = 0; i < m_spectra->size(); i++)
          if (usedSpectra[i])
            includedSpecs.push_back(i);

        int specCount = 0;
        for (unsigned int i = 0; i < usedSpectra.size(); i++)
          if (usedSpectra[i])
            specCount++;
        if (specCount == 0)
        {
          DEBUG_MSG("  - component is empty, ABruijn graph not built");
          continue;
        }

        DEBUG_MSG("  - ABruijn graph built on "<<specCount<<" spectra.");

        if (specCount > 10000)
        {
          DEBUG_MSG("Spectrum count per component exceeds maximum of 7500 spectra! Skipping component...");
          continue;
        }

        //
        // Find/count/list/split composite vertices
        //
        int compositeVertexCount = 0;
        list<int> compositeSet;
        for (unsigned int i = 0; i < vSet.vertices.size(); i++)
          if (vSet.vertices[i].size() > 0 and vSet.vertices[i].compositeVertex)
          {
            compositeVertexCount++;
            compositeSet.push_front(i);
            stringstream aux;
            aux << ">>>>> composite vertex " << i << ": ";
            for (list<TwoValues<int> >::iterator iter =
                vSet.vertices[i].specPeaks.begin(); iter
                != vSet.vertices[i].specPeaks.end(); iter++)
              aux << "(" << (*iter)[0] << "," << (*iter)[1] << ")";
            DEBUG_MSG(aux.str());

          }
        DEBUG_MSG("  - Abruijn graph contains "<<compositeVertexCount<<" composite vertices: ");
        stringstream aux;
        list<int>::iterator iter = compositeSet.begin();
        for (; iter != compositeSet.end(); iter++)
          aux << (*iter) << " ";
        if (compositeVertexCount > 0)
        {
          aux << "-> splitting...";
          vSet.splitComposite(spectrumGraphs, peakTol, &usedSpectra);
          aux << "done.";
        }
        DEBUG_MSG(aux.str());

        //
        // Add spectrum graph edges to ABruijn graph (may connect non-aligned peaks, e.g. if
        //   missing in all but one of the aligned spectra)
        //
        if (graphType == 1)
        {
          vSet.addEdges(spectrumGraphs, &usedSpectra);
          vSet.consolidatePaths();
        }
        if (graphType == 2)
        {
          vSet.addEdges(spectrumGraphs, &usedSpectra);
        }
        // Changed 2010/09/14			if(not addEndpoints) vSet.addEndpointEdges(tmpAlignsASP,matches,modPos,jumps2,peakTol);
        //				if(addEndpoints) vSet.addEndpointEdges(tmpAlignsASP,matches,modPos,jumps2,peakTol);
        if (addEndpoints)
          vSet.addEndpointEdges(components.cAlignsPA[cIdx],
                                matches,
                                modPos,
                                jumps2,
                                maxModMass,
                                peakTol,
                                pmTol);

        //
        //  ABruijn path-finishing procedures should go here.
        //

#ifdef DBG_MASAB
        vSet.outputGraph(m_labels);
#endif

        //
        // Create ABruijn MSGraph from VertexSet for b-ion and y-ion endpoints; choose option with highest-scoring path
        //
        VertexSet copy(vSet), *vSetP;
        MSGraph *abg, abgB, abgY, path;
        vector<int> vSet_indexB, vSet_indexY, *vSet_indexP; // Correspondences between vertex indices in the ABruijn graph and simplified graph
        vector<int> pathVertsIdxB, pathVertsIdxY, *pathVertsIdx; // Indices of the vertices in the heaviest path
        Spectrum specWithBep, specWithYep;
        float scoreWithBep, scoreWithYep;

        if (addEndpoints)
          vSet.removeEndpoints(false, peakTol); // Remove y-mass endpoints added above
        //
        //  ABruijn endpoint-finishing procedures should go here.
        //

        vSet.buildGraph(abgB,
                        jumps,
                        peakTol,
                        vSet_indexB,
                        m_labels,
                        edgeScoreType);
        //char filename[2048];   sprintf(filename,"component_%d.txt",cIdx+1);
        //abgB.output_graphviz(filename);
        if (noSequencing)
        {
          pathVertsIdxB.resize(0);
          if (!wholeABFN.empty())
            vSet.getMatchedPeaks(pathVertsIdxB, abVerticesAll[cIdx]);
          continue;
        }

        scoreWithBep = abgB.heaviestPath(path,
                                         false,
                                         &specWithBep,
                                         &pathVertsIdxB);

        if (addEndpoints)
          copy.removeEndpoints(true, peakTol); // Remove b-mass endpoints added above
        //
        //  ABruijn endpoint-finishing procedures should go here.
        //

        copy.buildGraph(abgY,
                        jumps,
                        peakTol,
                        vSet_indexY,
                        m_labels,
                        edgeScoreType);
        //sprintf(filename,"component_%d_y.txt",cIdx+1);
        //abgY.output_graphviz(filename);

        scoreWithYep = abgY.heaviestPath(path,
                                         false,
                                         &specWithYep,
                                         &pathVertsIdxY);

        //			scoreWithBep = max(scoreWithBep,scoreWithYep)+1;  // Force B endpoints


        if (scoreWithBep >= scoreWithYep - 0.00001)
        {
          abg = &abgB;
          m_contigShifts->consensus[cIdx] = specWithBep;
          vSet_indexP = &vSet_indexB;
          vSetP = &vSet;
          pathVertsIdx = &pathVertsIdxB;
          DEBUG_MSG("  - Selected B endpoints.");
        }
        else
        {
          abg = &abgY;
          m_contigShifts->consensus[cIdx] = specWithYep;
          vSet_indexP = &vSet_indexY;
          vSetP = &copy;
          pathVertsIdx = &pathVertsIdxY;
          DEBUG_MSG("  - Selected Y endpoints.");
        }

        for (unsigned int vIdx = 0; vIdx < pathVertsIdx->size(); vIdx++)
          (*pathVertsIdx)[vIdx] = (*vSet_indexP)[(*pathVertsIdx)[vIdx]]; // Convert simplified graph vertex indices to ABruijn vertex indices.

#ifdef DBG_MASAB
        //		sprintf(sBuf,"graph_ma_%d.txt",cIdx);
        //		vSetP->output_graphviz_ma(sBuf, *pathVertsIdx);
#endif
        vSetP->getMatchedPeaks(*pathVertsIdx, abVertices[cIdx]);
        m_contigShifts->set_endpoints(cIdx,
                                      *m_spectra,
                                      abVertices[cIdx],
                                      peakTol,
                                      pmTol);

        if (!wholeABFN.empty())
        {
          pathVertsIdx->resize(0);
          vSetP->getMatchedPeaks(*pathVertsIdx, abVerticesAll[cIdx]);
        }

        abg->info_heaviestPath(cStats[cIdx]);
        abCounts[cIdx][0] = abg->numVertices();
        abCounts[cIdx][1] = abg->numEdges();
        //		cout<<"  - Heaviest path stats: ["; for(unsigned int i=0;i<cStats[cIdx].size();i++) {cout<<cStats[cIdx][i]; if(i<cStats[cIdx].size()-1) cout<<", "; } cout<<"]\n";

#ifdef DBG_MASAB
        sprintf(sBuf,"graph_%d.txt",cIdx); abg->output_graphviz(sBuf);
#endif

        //
        // Find spectra with not enough ABruijn vertices on the heaviest path (if any)
        //   and remove them from the current component. Unused spectra may define new (leftover)
        //   components if connected by at least one edge (ie, at least two spectra in a
        //   connected component).
        //
        list<int> usedSpecs;
        usedSpecs.clear();
        list<int> unusedSpecs;
        unusedSpecs.clear();
        vector<short> numPathPeaks(m_spectra->size());
        for (unsigned int i = 0; i < m_spectra->size(); i++)
          numPathPeaks[i] = 0;
        for (unsigned int i = 0; i < abVertices[cIdx].size(); i++)
        {
          //cerr<<"AB vertex "<<i<<": ";
          list<TwoValues<int> >::iterator vNext; // Used to remove remaining composite vertices
          for (list<TwoValues<int> >::iterator vIter =
              abVertices[cIdx][i].begin(); vIter != abVertices[cIdx][i].end();)
          {
            //cerr<<"("<<(*vIter)[0]<<","<<(*vIter)[1]<<")";
            vNext = vIter;
            if (vNext != abVertices[cIdx][i].end())
              vNext++;
            if (vNext != abVertices[cIdx][i].end() and (*vIter)[0]
                == (*vNext)[0])
            {
              DEBUG_MSG("  - ERROR: inconsistent contig vertex containing ("<<(*vIter)[0]<<","<<(*vIter)[1]<<") and ("<<(*vNext)[0]<<","<<(*vNext)[1]<<")");
              vIter = abVertices[cIdx][i].erase(vIter); // Remove inconsistent spectrum/peaks
              vIter = abVertices[cIdx][i].erase(vIter);
            }
            else
            {
              if (++numPathPeaks[(*vIter)[0]] >= pathMinPeaks)
                usedSpecs.push_back((*vIter)[0]);
              vIter++;
            }
          }
          //cerr<<endl;
        }
        usedSpecs.sort();
        usedSpecs.unique();

        // Prune ABruijn vertices: remove peaks from spectra without enough matches to the consensus path
        unsigned int maxPeaksPerVertex = 0;
        for (unsigned int i = 0; i < abVertices[cIdx].size(); i++)
        {
          for (list<TwoValues<int> >::iterator vIter =
              abVertices[cIdx][i].begin(); vIter != abVertices[cIdx][i].end();)
            if (numPathPeaks[(*vIter)[0]] >= pathMinPeaks)
              vIter++;
            else
            {
              unusedSpecs.push_back((*vIter)[0]);
              vIter = abVertices[cIdx][i].erase(vIter);
            }
          if (abVertices[cIdx][i].size() > maxPeaksPerVertex)
            maxPeaksPerVertex = abVertices[cIdx][i].size();
          if (abVertices[cIdx][i].size() > usedSpecs.size())
          {
            DEBUG_MSG("  - ERROR: inconsistent contig vertex: ");
            for (list<TwoValues<int> >::iterator vIter =
                abVertices[cIdx][i].begin(); vIter != abVertices[cIdx][i].end(); vIter++)
            {
              DEBUG_MSG("("<<(*vIter)[0]<<","<<(*vIter)[1]<<")");
            }
          }
        }
        unusedSpecs.sort();
        unusedSpecs.unique();

        // Output complete set of vertices
        if (!wholeABFN.empty())
        {
          for (unsigned int i = 0; i < abVerticesAll[cIdx].size(); i++)
            for (list<TwoValues<int> >::iterator vIter =
                abVerticesAll[cIdx][i].begin(); vIter
                != abVerticesAll[cIdx][i].end();)
              if (numPathPeaks[(*vIter)[0]] >= pathMinPeaks)
                vIter++;
              else
                vIter = abVerticesAll[cIdx][i].erase(vIter);
        }

        if (usedSpecs.size() >= pathMinSpecs)
        {
          if (usedSpecs.size() < components.sets[cIdx].size() - 1)
          { // Whenever there are at least 2 unused spectra
            stringstream logMsg;
            logMsg << "  - Keeping " << usedSpecs.size()
                << " spectra; number of components: " << components.size()
                << " -> ";
            components.splitSet(cIdx, usedSpecs);
            components.splitSet(cIdx, usedSpecs);
            logMsg << components.size();
            DEBUG_MSG(logMsg.str());
          }
          else if (usedSpecs.size() == components.sets[cIdx].size() - 1)
            components.removeElements(cIdx, unusedSpecs);
          if (maxPeaksPerVertex > usedSpecs.size())
          {
            DEBUG_MSG("  - ERROR: inconsistent contig containing a vertex with "<<maxPeaksPerVertex<<" spectrum peaks but only "<<usedSpecs.size()<<" used spectra (contig deleted)");
            m_contigShifts->consensus[cIdx].resize(0);
            abVertices[cIdx].resize(0);
          }
        }
        else
        {
          if (includedSpecs.size() > 0)
          { // Whenever there is at least 1 included spectrum
            stringstream logMsg;
            logMsg << "  - Keeping " << components.sets[cIdx].size()
                - includedSpecs.size() << " spectra; number of components: "
                << components.size() << " -> ";
            components.splitSet(cIdx, includedSpecs);
            logMsg << components.size();
            DEBUG_MSG(logMsg.str());
          }
          DEBUG_MSG("  - Only "<<usedSpecs.size()<<" spectra with at least "<<pathMinPeaks<<" peaks in the consensus path - component deleted.");
          m_contigShifts->consensus[cIdx].resize(0); // Remove de novo sequences for poor contigs (not enough spectra with enough matched peaks)
        }

        if (components.cAlignsPA[cIdx].size() == 0
            and components.cAlignsASP[cIdx].size() == 0)
        { // No pairs left
          prevNumSpecs = components.sets[cIdx].size();
          DEBUG_MSG("  --> No pairs left for second iteration - keeping ABruijn path");
        }
      } // while (prevNumSpecs!=components.sets[cIdx].size())
    }
    // Resize down to the final number of resulting connected components
    abCounts.resize(components.size());
    cStats.resize(components.size());
    m_contigShifts->resize(components.size());
    abVertices.resize(components.size());
    abVerticesAll.resize(components.size());

    Merge_abinfo_v1_0(*m_spectra,
                      components.sets,
                      specFlipped,
                      abVertices,
                      *m_outputAbruijn);

    // Output debug/development details
    //    components.saveas_binListArray("components.bla");
    //    Save_binArray("component_stats.bna",cStats);
    //    Save_binArray("abCounts.bin", abCounts);
    if (m_params.exists("OUTPUT_COMPLETE_ABRUIJN"))
    {
      string aux = m_params.getValue("OUTPUT_COMPLETE_ABRUIJN");
      if (aux.size())
        Save_abinfo(aux.c_str(),
                    *m_spectra,
                    components.sets,
                    specFlipped,
                    abVerticesAll);
    }

    if (parallelPaths)
    {
      ParameterList pParams(m_params);
      pParams.setValue("ENFORCE_B_ENDPTS", "1");
      pParams.setValue("REVERSE_STARS", "1");

      Clusters newShifts;
      newShifts = *m_contigShifts;
      abinfo_t newAbruijn;

      DEBUG_MSG("Invoking ParallelAssembly ...");
      ExecParallelAssembly pAssembly(pParams,
                                     m_outputAbruijn,
                                     m_spectra,
                                     &newShifts,
                                     &newAbruijn);

      if (!pAssembly.invoke())
      {
        ERROR_MSG("Failed to invoke ParallelAssembly!!!");
        abort();
      }

      m_outputAbruijn->operator =(newAbruijn);
      m_contigShifts->operator =(newShifts);
    }

    // Need to assign scan numbers so that PSMs will associate properly
    for (int i = 0; i < m_contigShifts->consensus.size(); i++)
    {
      m_contigShifts->consensus[i].scan = i + 1;
    }

    DEBUG_VAR(m_contigShifts->size());

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecAssembly::loadInputData(void)
  {
    DEBUG_TRACE;
    if (ownInput)
    {
      if (!m_spectra) {
        m_spectra = new SpecSet;
      }
      if (!m_spectrumPairs) {
        m_spectrumPairs = new SpectrumPairSet;
      }
    }
    DEBUG_VAR(m_spectra);
    m_spectra->resize(0);
    DEBUG_VAR(m_spectrumPairs);
    m_spectrumPairs->resize(0);
    DEBUG_VAR(&m_labels);
    m_labels.resize(0);

    if (ownOutput)
    {
      if (!m_contigShifts)
        m_contigShifts = new Clusters;
      if (!m_outputAbruijn)
        m_outputAbruijn = new abinfo_t;
    }
    m_contigShifts->resize(0);

    if (m_params.exists("AMINO_ACID_MASSES"))
    {
      AAJumps tmpJumps(-1);
      tmpJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true); // Set global defaults for amino acid masses
    }

    if (m_spectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str())
        <= 0 or m_spectra->size() == 0)
    {
      ERROR_MSG("Error reading input spectra from "<< m_params.getValue("INPUT_SPECS_PKLBIN"));
      return false;
    }

    if (!m_spectrumPairs->loadFromBinaryFile(m_params.getValue("INPUT_ALIGNS")))
    {
      ERROR_MSG("Error reading spectral pairs from "<<m_params.getValue("INPUT_ALIGNS"));
      return false;
    }

    if (m_params.exists("INPUT_LABELS"))
    {
      if (LoadLabels(*m_spectra,
                     m_params.getValue("INPUT_LABELS").c_str(),
                     m_labels) <= 0)
      {
        ERROR_MSG("Error reading labels file "<<m_params.getValue("INPUT_LABELS"));
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecAssembly::saveOutputData(void)
  {
    string dataDir;
    if (m_params.exists("OUTPUT_SPECTRA_PATH")) {
      dataDir = m_params.getValue("OUTPUT_SPECTRA_PATH");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      dataDir = dataDir + "/";
    }
    DEBUG_VAR(dataDir);
    if (m_params.exists("OUTPUT_CLUSTERS")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_CLUSTERS");
      DEBUG_MSG("Saving clusters to: [" << fileName << "]");
      m_contigShifts->Save(fileName.c_str());
    }

    if (m_params.exists("OUTPUT_SPECS")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_SPECS");
      DEBUG_MSG("Saving spectra to: [" << fileName << "]");
      m_contigShifts->consensus.savePklBin(fileName.c_str());
    }

    if (m_params.exists("OUTPUT_ABINFO")) {
      string fileName = dataDir + m_params.getValue("OUTPUT_ABINFO");
      DEBUG_MSG("Saving abinfo to: [" << fileName << "]");
      Save_abinfo_v1_0(fileName.c_str(), *m_outputAbruijn);
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecAssembly::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecAssembly::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecAssembly::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecAssembly::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecAssembly::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("PENALTY_PTM");
    VALIDATE_PARAM_EXIST("PENALTY_SAME_VERTEX");
    //    VALIDATE_PARAM_EXIST("GRAPH_TYPE");
    //    VALIDATE_PARAM_EXIST("MAX_AA_JUMP");
    //    VALIDATE_PARAM_EXIST("MAX_MOD_MASS");
    //    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    //    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    //    VALIDATE_PARAM_EXIST("EDGE_SCORE_TYPE");
    //    VALIDATE_PARAM_EXIST("MIN_MATCHED_PEAKS");
    //    VALIDATE_PARAM_EXIST("MIN_EDGES_TO_COMPONENT");
    //    VALIDATE_PARAM_EXIST("PATH_MIN_SPECS");
    //    VALIDATE_PARAM_EXIST("PATH_MIN_PEAKS");
    //    VALIDATE_PARAM_EXIST("SPEC_TYPE_MSMS");
    //    VALIDATE_PARAM_EXIST("NO_SEQUENCING");
    //    VALIDATE_PARAM_EXIST("ADD_ENDPOINTS");
    //    VALIDATE_PARAM_EXIST("OUTPUT_COMPLETE_ABRUIJN");

    m_isValid = true;
    return true;
  }

// -------------------------------------------------------------------------


}
