/*
 * CombineContigs.cpp
 *
 *  Created on: Jun 11, 2012
 *      Author: aguthals
 */

#include "CombineContigs.h"

#include "ExecFramework/ExecAssembly.h"

using namespace std;

namespace specnets
{

  const char* DEBUG_OUTPUT = "debug_combinations.csv";
  const float minContigOverlap = 400.0;
  const float minSecondaryRatio = 0;
  const float minParentMassSep = 10.0;

  struct SortEdges : public std::binary_function<TwoValues<float> , TwoValues<
      float> , bool>
  {
    bool operator()(TwoValues<float> left, TwoValues<float> right) const
    {
      return left[1] > right[1];
    }
    ;
  };

  int FindMatchPeaksAll22(Spectrum &spec1,
                         Spectrum &spec2,
                         float shift,
                         float tolerance,
                         vector<int> &idx1,
                         vector<int> &idx2)
  {
    int i, j; // Iterators over the peaks indices
    int low = 0, high = 0; // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1

    idx1.resize(0);
    idx2.resize(0);
    for (i = 0; i < (int)spec1.size(); i++)
    {
      while (low > 0 && (spec1[i][0] - tolerance - 0.000001) < (spec2[low][0]
          + shift))
        low--;
      while ((low < (int)spec2.size()) && (spec1[i][0] - tolerance - 0.000001)
          > (spec2[low][0] + shift))
        low++;
      while ((high < (int)spec2.size()) && (spec1[i][0] + tolerance + 0.000001)
          >= (spec2[high][0] + shift))
        high++; // high is index of first unreachable peak
      for (j = low; j < high; j++)
      {
        idx1.push_back(i);
        idx2.push_back(j);
      }
    }
    return idx1.size();
  }

  void mergeCommon(Spectrum& candSpec,
                   Spectrum &withSpectrum,
                   Spectrum *toSpec,
                   float shift,
                   float peakTol,
                   short mergeType)
  {
    vector<int> idx1, idx2;
    unsigned int pivot;
    FindMatchPeaksAll22(candSpec, withSpectrum, shift, peakTol, idx1, idx2);
    Spectrum tmpToSpec;
    tmpToSpec.resize(idx1.size());
    tmpToSpec.copyNP(candSpec); // Temporary spectrum, in case toSpec==0

    float s, r1, r2;
    for (pivot = 0; pivot < idx1.size(); pivot++)
    {
      s = candSpec[idx1[pivot]][1] + withSpectrum[idx2[pivot]][1];
      r1 = candSpec[idx1[pivot]][1] / s;
      r2 = withSpectrum[idx2[pivot]][1] / s;
      tmpToSpec[pivot][0] = r1 * candSpec[idx1[pivot]][0] + r2
          * (withSpectrum[idx2[pivot]][0] + shift);

      switch (mergeType)
      {
      case 0:
        tmpToSpec[pivot][1] = candSpec[idx1[pivot]][1]
            + withSpectrum[idx2[pivot]][1];
        break;
      case 1:
        tmpToSpec[pivot][1] = max(candSpec[idx1[pivot]][1],
                                  withSpectrum[idx2[pivot]][1]);
        break;
      case 2:
        tmpToSpec[pivot][1] = min(candSpec[idx1[pivot]][1],
                                  withSpectrum[idx2[pivot]][1]);
        break;
      case 3:
        tmpToSpec[pivot][1] = candSpec[idx1[pivot]][1];
        break;
      case 4:
        tmpToSpec[pivot][1] = withSpectrum[idx2[pivot]][1];
        break;
      }
    }
    tmpToSpec.sortPeaks();
    (*toSpec) = tmpToSpec;
  }

  TwoValues<int> countConsecutiveMP(Spectrum& contig1,
                                    Spectrum& contig2,
                                    float shiftFor,
                                    float shiftRev,
                                    float pmTol,
                                    float minOverlapArea,
                                    bool debug)
  {
    vector<int> idx1_f, idx2_f, idx1_r, idx2_r;
    Spectrum rev2;
    contig2.reverse(0.0 - AAJumps::massH2O, &rev2);
    FindMatchPeaksAll22(contig1, contig2, shiftFor, pmTol, idx1_f, idx2_f);
    FindMatchPeaksAll22(contig1, rev2, shiftRev, pmTol, idx1_r, idx2_r);

    vector<bool> c1_f(contig1.size()), c1_r(contig1.size()),
        c2_f(contig2.size()), c2_r(contig2.size());

    for (int i = 0; i < idx1_f.size(); i++)
    {
      c1_f[idx1_f[i]] = true;
    }
    for (int i = 0; i < idx1_r.size(); i++)
    {
      c1_r[idx1_r[i]] = true;
    }
    for (int i = 0; i < idx2_f.size(); i++)
    {
      c2_f[idx2_f[i]] = true;
    }
    for (int i = 0; i < idx2_r.size(); i++)
    {
      c2_r[idx2_r[i]] = true;
    }

    int cons_peaks = 0, max_cons_f1 = 0, max_cons_f2 = 0, max_cons_r1 = 0,
        max_cons_r2 = 0;
    for (int i = 0; i < c1_f.size(); i++)
    {
      if (c1_f[i])
      {
        cons_peaks++;
        if (cons_peaks > max_cons_f1)
        {
          max_cons_f1 = cons_peaks;
        }
      }
      else
      {
        cons_peaks = 0;
      }
    }
    cons_peaks = 0;
    for (int i = 0; i < c1_r.size(); i++)
    {
      if (c1_r[i])
      {
        cons_peaks++;
        if (cons_peaks > max_cons_r1)
        {
          max_cons_r1 = cons_peaks;
        }
      }
      else
      {
        cons_peaks = 0;
      }
    }
    cons_peaks = 0;
    for (int i = 0; i < c2_f.size(); i++)
    {
      if (c2_f[i])
      {
        cons_peaks++;
        if (cons_peaks > max_cons_f2)
        {
          max_cons_f2 = cons_peaks;
        }
      }
      else
      {
        cons_peaks = 0;
      }
    }
    cons_peaks = 0;
    for (int i = 0; i < c2_r.size(); i++)
    {
      if (c2_r[i])
      {
        cons_peaks++;
        if (cons_peaks > max_cons_r2)
        {
          max_cons_r2 = cons_peaks;
        }
      }
      else
      {
        cons_peaks = 0;
      }
    }
    max_cons_f1 = min(max_cons_f1, max_cons_f2);
    max_cons_r1 = min(max_cons_r1, max_cons_r2);
    return TwoValues<int> (max_cons_f1, max_cons_r2);
  }

  TwoValues<float> getContigOverlapScores(Spectrum& contig1,
                                          Spectrum& contig2,
                                          float shiftFor,
                                          float shiftRev,
                                          float pmTol,
                                          float minOverlapArea,
                                          bool debug)
  {
    Spectrum contig = contig1;
    Spectrum f_contig_o = contig2;

    Spectrum r_contig_o;
    contig2.reverse(0.0 - AAJumps::massH2O, &r_contig_o);

    Spectrum f_contig = f_contig_o;
    Spectrum r_contig = r_contig_o;
    unsigned int idxUse = 0;
    for (unsigned int i = 0; i < f_contig_o.size(); i++)
    {
      float mass = f_contig_o[i][0] + shiftFor;
      if (mass > 0.0 - pmTol && mass < contig.parentMass + pmTol)
      {
        f_contig[idxUse] = f_contig_o[i];
        idxUse++;
      }
    }
    f_contig.resize(idxUse);
    idxUse = 0;
    for (unsigned int i = 0; i < r_contig_o.size(); i++)
    {
      float mass = r_contig_o[i][0] + shiftRev;
      if (mass > 0.0 - pmTol && mass < contig.parentMass + pmTol)
      {
        r_contig[idxUse] = r_contig_o[i];
        idxUse++;
      }
    }
    r_contig.resize(idxUse);
    idxUse = 0;

    if (debug)
    {
      cout << "Forward Shift: " << shiftFor << "\nReverse Shift: " << shiftRev
          << "\nSpectrum:\n";
      contig.output(cout);
      cout << "\nForward Contig2:\n";
      for (unsigned int i = 0; i < f_contig.size(); i++)
      {
        cout << f_contig[i][0] << " + " << shiftFor << " = " << f_contig[i][0]
            + shiftFor << "\n";
      }
      cout << "\nReverse Contig2:\n";
      for (unsigned int i = 0; i < r_contig.size(); i++)
      {
        cout << r_contig[i][0] << " + " << shiftRev << " = " << r_contig[i][0]
            + shiftRev << "\n";
      }
    }

    Spectrum f_merged;
    Spectrum f_merged2;
    Spectrum r_merged;
    Spectrum r_merged2;

    if (f_contig.size() > 0)
    {
      //contig.mergeClosestPeaks();
      mergeCommon(contig, f_contig, &f_merged, shiftFor, pmTol, 3);
      mergeCommon(contig, f_contig, &f_merged2, shiftFor, pmTol, 4);
    }
    if (r_contig.size() > 0)
    {
      mergeCommon(contig, r_contig, &r_merged, shiftRev, pmTol, 3);
      mergeCommon(contig, r_contig, &r_merged2, shiftRev, pmTol, 4);
    }

    float f_2_min = 0.0 - shiftFor - pmTol;
    float f_2_max = (*contig1.back())[0] - shiftFor + pmTol;
    float r_2_min = 0.0 - shiftRev - pmTol;
    float r_2_max = (*contig1.back())[0] - shiftRev + pmTol;
    float f_1_min = shiftFor - pmTol;
    float f_1_max = shiftFor + (*contig2.back())[0] + pmTol;
    float r_1_min = shiftRev - pmTol;
    float r_1_max = shiftRev + (*contig2.back())[0] + pmTol;

    if (debug)
    {
      cout << "pmTol = " << pmTol << endl;
      cout << "f_1_min = " << f_1_min << ", f_1_max = " << f_1_max << endl;
      cout << "f_2_min = " << f_2_min << ", f_2_max = " << f_2_max << endl;
      cout << "r_1_min = " << r_1_min << ", r_1_max = " << r_1_max << endl;
      cout << "r_2_min = " << r_2_min << ", r_2_max = " << r_2_max << endl;

      cout << "\nMerged Forward 1: \n";
      f_merged.output(cout);
      cout << "\nMerged Forward 2: \n";
      f_merged2.output(cout);
      cout << "\nMerged Reverse 1: \n";
      r_merged.output(cout);
      cout << "\nMerged Reverse 2: \n";
      r_merged2.output(cout);
      /*
       cout << "\nContig1  parent mass: " << contig.parentMass << "\nContig2 parent mass: " << f_contig.parentMass << " = " << r_contig.parentMass << "\n";
       cout << "Forward 1 min - max: " << f_1_min << " " << f_1_max << "\n";
       cout << "Reverse 1 min - max: " << r_1_min << " " << r_1_max << "\n";
       cout << "Forward 2 min - max: " << f_2_min << " " << f_2_max << "\n";
       cout << "Reverse 2 min - max: " << r_2_min << " " << r_2_max << "\n\n\n";
       */
    }

    float f_intensity = 0, f_intensity2 = 0;
    for (unsigned int i = 0; i < contig.size(); i++)
    {
      if (contig[i][0] > f_1_min && contig[i][0] < f_1_max)
      {
        f_intensity += contig[i][1];
      }
    }
    for (unsigned int i = 0; i < f_contig.size(); i++)
    {
      if (f_contig[i][0] > f_2_min && f_contig[i][0] < f_2_max)
      {
        f_intensity2 += f_contig[i][1];
      }
    }

    float r_intensity = 0, r_intensity2 = 0;
    for (unsigned int i = 0; i < contig.size(); i++)
    {
      if (contig[i][0] > r_1_min && contig[i][0] < r_1_max)
      {
        r_intensity += contig[i][1];
      }
    }
    for (unsigned int i = 0; i < r_contig.size(); i++)
    {
      if (r_contig[i][0] > r_2_min && r_contig[i][0] < r_2_max)
      {
        r_intensity2 += r_contig[i][1];
      }
    }

    float m_f_intensity = 0;
    for (unsigned int i = 0; i < f_merged.size(); i++)
      m_f_intensity += f_merged[i][1];
    float m_r_intensity = 0;
    for (unsigned int i = 0; i < r_merged.size(); i++)
      m_r_intensity += r_merged[i][1];
    float m_f_intensity2 = 0;
    for (unsigned int i = 0; i < f_merged2.size(); i++)
      m_f_intensity2 += f_merged2[i][1];
    float m_r_intensity2 = 0;
    for (unsigned int i = 0; i < r_merged2.size(); i++)
      m_r_intensity2 += r_merged2[i][1];

    float minF = min((m_f_intensity / f_intensity), (m_f_intensity2
        / f_intensity2));
    float minR = min((m_r_intensity / r_intensity), (m_r_intensity2
        / r_intensity2));

    if (debug)
    {
      cout << "m_f_intensity = " << m_f_intensity << ", m_f_intensity2 = "
          << m_f_intensity2 << endl;
      cout << "f_intensity = " << f_intensity << ", f_intensity2 = "
          << f_intensity2 << endl;
      cout << "m_r_intensity = " << m_r_intensity << ", m_r_intensity2 = "
          << m_r_intensity2 << endl;
      cout << "r_intensity = " << r_intensity << ", r_intensity2 = "
          << r_intensity2 << endl;
      cout << "minF = " << minF << ", minR = " << minR << endl;
    }

    if (min(f_1_max, contig.parentMass) - max((float)0, f_1_min)
        < minOverlapArea)
    {
      minF = -1.0;
    }
    else if (f_merged.size() == 0)
    {
      minF = 0;
    }

    if (min(r_1_max, contig.parentMass) - max((float)0, r_1_min)
        < minOverlapArea)
    {
      minR = -1.0;
    }
    else if (r_merged.size() == 0)
    {
      minR = 0;
    }

    return TwoValues<float> (minF, minR);
  }

  TwoValues<int> getContigOverlapPeaks(Spectrum& contig1,
                                       Spectrum& contig2,
                                       float shiftFor,
                                       float shiftRev,
                                       float pmTol,
                                       float minOverlapArea,
                                       bool debug)
  {
    vector<int> idx1_f, idx2_f, idx1_r, idx2_r;
    Spectrum rev2;
    contig2.reverse(0.0 - AAJumps::massH2O, &rev2);
    FindMatchPeaksAll22(contig1, contig2, shiftFor, pmTol, idx1_f, idx2_f);
    FindMatchPeaksAll22(contig1, rev2, shiftRev, pmTol, idx1_r, idx2_r);

    vector<bool> c1_f(contig1.size()), c1_r(contig1.size()),
        c2_f(contig2.size()), c2_r(contig2.size());

    for (int i = 0; i < idx1_f.size(); i++)
    {
      c1_f[idx1_f[i]] = true;
    }
    for (int i = 0; i < idx1_r.size(); i++)
    {
      c1_r[idx1_r[i]] = true;
    }
    for (int i = 0; i < idx2_f.size(); i++)
    {
      c2_f[idx2_f[i]] = true;
    }
    for (int i = 0; i < idx2_r.size(); i++)
    {
      c2_r[idx2_r[i]] = true;
    }

    int max_f1 = 0, max_f2 = 0, max_r1 = 0, max_r2 = 0;
    for (int i = 0; i < c1_f.size(); i++)
    {
      if (c1_f[i])
      {
        max_f1++;
      }
    }
    for (int i = 0; i < c1_r.size(); i++)
    {
      if (c1_r[i])
      {
        max_r1++;
      }
    }
    for (int i = 0; i < c2_f.size(); i++)
    {
      if (c2_f[i])
      {
        max_f2++;
      }
    }
    for (int i = 0; i < c2_r.size(); i++)
    {
      if (c2_r[i])
      {
        max_r2++;
      }
    }
    max_f1 = min(max_f1, max_f2);
    max_r1 = min(max_r1, max_r2);
    return TwoValues<int> (max_f1, max_r1);
  }

  void copy(CombineContigsParams& from, CombineContigsParams& to)
  {
    to.contigs = from.contigs;
    to.star_spectra = from.star_spectra;
    //to.consensus = from.consensus;
    to.in_abinfo = from.in_abinfo;
    //to.out_abinfo = from.out_abinfo;
    //to.out_reversed = from.out_reversed;
    to.contig_alignments = from.contig_alignments;

    to.peak_tol = from.peak_tol;
    to.parent_mass_tol = from.parent_mass_tol;
    to.contig_overlap_score = from.contig_overlap_score;
    to.min_matched_peaks_contig_align = from.min_matched_peaks_contig_align;
    to.resolution = from.resolution;

    //to.matepairs = from.matepairs;
    // to.matepair_alignments = from.matepair_alignments;

    // to.min_matepair_charge = from.min_matepair_charge;
    //to.min_ratio_matepair_align = from.min_ratio_matepair_align;
    // to.min_ratio_matepair_contig_align = from.min_ratio_matepair_contig_align;
    //to.min_matched_peaks_matepair_align = from.min_matched_peaks_matepair_align;
    //to.min_num_matepairs = from.min_num_matepairs;
    //to.min_edges_to_component = from.min_edges_to_component;
    to.min_component_size = from.min_component_size;

    to.contig_prot_idx = from.contig_prot_idx;
    to.contig_peak_idx = from.contig_peak_idx;
    to.proteins = from.proteins;
    to.protein_idx_ignore = from.protein_idx_ignore;
  }

  /*
   struct CombineContigsParams {
   //required contig-contig merging parameters
   SpecSet* contigs;
   SpecSet consensus;
   abinfo_t* in_abinfo;
   abinfo_t out_abinfo;
   vector<int> out_reversed;
   list<SpectrumPair>* contig_alignments;

   float peak_tol;
   float parent_mass_tol;
   float resolution;
   float contig_overlap_score;
   int min_matched_peaks_contig_align;

   //required contig-matepair merging parameters
   SpecSet* matepairs;
   list<SpectrumPair>* matepair_alignments;

   short min_matepair_charge;
   float min_ratio_matepair_align;
   float min_ratio_matepair_contig_align;
   int min_matched_peaks_matepair_align;
   int min_num_matepairs;
   int min_edges_to_component;
   int min_component_size;

   //required matchma parameters
   vector<vector<int>* contig_prot_idx;
   SpecSet* contig_peak_idx;
   DB_fasta* proteins;
   set<int>* protein_idx_ignore;
   };
   */

  CombineContigs::CombineContigs()
  {
  }
  CombineContigs::~CombineContigs()
  {
  }

  CombineContigs::CombineContigs(CombineContigsParams* params)
  {
    construct(params);
  }

  void CombineContigs::construct(CombineContigsParams* params)
  {
    merging_params = params;

    if (params->contig_prot_idx != NULL && params->contig_peak_idx != NULL
        && params->proteins != NULL && params->protein_idx_ignore != NULL)
    {
      overlaps = *(params->contig_peak_idx);
      prot_match = *(params->contig_prot_idx);
      fasta = *(params->proteins);
      idxNotUsed = *(params->protein_idx_ignore);
      haveRes = true;
    }
    else
    {
      haveRes = false;
    }

    contigs = *(merging_params->contigs);
    contig_abinfo = *(params->in_abinfo);

    if (params->contig_alignments != NULL)
      contig_alignments = *(params->contig_alignments);

    haveMatePairs = false;

    for (unsigned int i = 0; i < contigs.size(); i++)
    {
      if (contigs[i].size() == 0)
        continue;
      contigs[i].parentMass = (*contigs[i].back())[0] + AAJumps::massMH;
      //contigs[i].annotation.resize(0);
      //contigs[i].ionTypes.resize(0);
    }
    peakTol = params->peak_tol;
    parentMassTol = params->parent_mass_tol;
    intPeakTol = floatToInt(peakTol / merging_params->resolution);
    intParentMassTol = floatToInt(parentMassTol / merging_params->resolution);
    init();
  }

  // mergeType = 0 (combine w/ contig shifts), 1 (combine w/ connector shifts)
  void CombineContigs::combineEdges(short mergeType)
  {

    init();

    //short precCharge = merging_params->min_matepair_charge;
    //float minRatioConn = merging_params->min_ratio_matepair_align;
    //float minRatioConnContig = merging_params->min_ratio_matepair_contig_align;
   // int minNumMatchedPeaksConn =
    //    merging_params->min_matched_peaks_matepair_align;
    //int minNumConnectors = merging_params->min_num_matepairs;
    float minCombScore = merging_params->contig_overlap_score;
    int minComp = merging_params->min_component_size;
    int minContigMP = merging_params->min_matched_peaks_contig_align;

    float contigLen = 0, numContigs = 0, maxContigL = 0, contigA = 0;
    int maxI = 0, maxL = 0, maxLI = 0;
    for (int i = 0; i < contigs.size(); i++)
    {
      if (contigs[i].size() == 0)
        continue;
      contigLen += contigs[i].parentMass;
      numContigs += 1.0;
      if (contigs[i].parentMass > maxContigL)
      {
        maxContigL = contigs[i].parentMass;
        maxI = i;
      }
      if (haveRes && prot_match[i][0] >= 0)
      {
        float cover = overlaps[i][overlaps[i].size() - 1][1]
            - overlaps[i][0][1] + 1.0;
        contigA += cover;
        if ((int)cover > maxL)
        {
          maxL = (int)cover;
          maxLI = i;
        }
      }
    }
    printf("\nMaximum contig length (Da) : %.1f (index %d)\n", maxContigL, maxI);
    printf("Average contig length (Da) : %.1f (%.0f contigs found of possible %d)\n",
           contigLen / numContigs,
           numContigs,
           contigs.size());
    if (haveRes)
    {
      printf("Maximum contig AA length : %d (index %d)\n", maxL, maxLI);
      printf("Average contig AA length : %.1f\n", contigA / numContigs);
    }
    root_alignments.clear();
    root_alignments.resize(contigs.size());
    list<TwoValues<int> > comp_alignments;
    TwoValues<int> node_pair;

    graph.clear();
    map<int, map<int, int> > graphBefMerg;
    map<int, int>::iterator mapIt;
    edges.resize(0);
    vector<vector<float> > edgesBefMerg;
    list<TwoValues<float> > scoredEdges;
    list<TwoValues<float> >::iterator scoredEdgeIt;
    components.clear();
    components.resize(contigs.size());
    vector<bool> contigsMerged(contigs.size());
    for (int i = 0; i < contigsMerged.size(); i++)
      contigsMerged[i] = false;
    vector<set<int> > componentsBefMerg;
    set<int> erasedEdges;
    set<int>::iterator nodeIt;
    set<int>::iterator nodeIt2;
    TwoValues<float> res, res1, res2;
    vector<float> edgeTo(6);
    vector<float> edgeFrom(6);
    reversed.clear();
    reversed.resize(contigs.size());
    for (int i = 0; i < reversed.size(); i++)
      reversed[i] = false;
    vector<bool> reversedContigMerg;
    vector<bool> reversedBefMerg;
    SpecSet myContigsContigMerg;
    SpecSet myContigsBefMerg;
    for (int i = 0; i < contigs.size(); i++)
      components[i].insert(i);

    int maxContigEdgeIdx = -1;

    if (mergeType == 0)
    {
      getCondensedContigAlignmentGraph(graph,
                                       edges,
                                       scoredEdges,
                                       minCombScore,
                                       minContigMP);
    }
    vector<vector<float> > originalEdges = edges;

    if (mergeType == 1)
    {
      scoredEdge[0] = -1.0;
      scoredEdge[1] = -1.0;
      scoredEdges.push_back(scoredEdge);
    }

    cout << "Initialized with " << edges.size() << " edges\n";

    graphBefMerg = graph;
    edgesBefMerg = edges;
    componentsBefMerg = components;
    reversedBefMerg = reversed;
    myContigsBefMerg = oriented;
    //outputGraph(cout, graph, edges);

    vector<int> edgeRef(edges.size());
    for (int i = 0; i < edgeRef.size(); i++)
      edgeRef[i] = -1;

    float FFShift, FRShift, RFShift, RRShift, oldFFShift, oldFRShift,
        oldRFShift, oldRRShift, conFFShift, conFRShift, conRFShift, conRRShift,
        edgeScore;
    float goodEdges = 0, totalEdges = 0, edgesAdded = 0, goodEdgesAdded = 0,
        zero = 0;

    FILE* output = fopen(DEBUG_OUTPUT, "wb");

    fprintf(output, "Adding Contig Edges\n");
    cout << "\nadding contig edges...\n";

    if (haveRes)
    {
      fprintf(output,
              "Protein Idx1%sProtein Idx2%sContig Idx1%sContig Idx2%sExperimental Shift FF%sExperimental Shift RR%sExpected Shift FF%sExpected Shift RR%sEqual?%sRemoved?%sContig2 Reversed?%sContig1 Reversed Global?%sContig2 Reversed Global?%sForward Score%sReverse Score%sComponent Overlap (Da)\n",
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP);
    }
    else
    {
      fprintf(output,
              "Contig Idx1%sContig Idx2%sExperimental Shift FF%sExperimental Shift RR%sRemoved?%sContig2 Reversed?%sForward Score%sReverse Score%sComponent Overlap (Da)\n",
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP);
    }

    bool checkOverlap = true;
    bool mergingMatePairs = false;
    int countmatep = 0, countcont = 0, countmerg = 0;

    scoredEdgeIt = scoredEdges.begin();
    while (scoredEdgeIt != scoredEdges.end())
    {
      int edgeIdx = floatToInt((*scoredEdgeIt)[0]);
      int origIdx = edgeIdx;
      edgeScore = (*scoredEdgeIt)[1];
      int idx1, idx2;

      if (edgeIdx >= 0)
      {
        while (edgeRef[edgeIdx] != -1)
          edgeIdx = edgeRef[edgeIdx];
        idx1 = floatToInt(edges[edgeIdx][0]);
        idx2 = floatToInt(edges[edgeIdx][1]);

        if (erasedEdges.count(edgeIdx) > 0 || components[idx1].count(idx2) > 0
            || components[idx2].count(idx1) > 0 || graph.count(idx1) == 0
            || graph[idx1].count(idx2) == 0)
        {
          scoredEdgeIt++;
          continue;
        }

      }
      else
      {
        ERROR_MSG("No longer support meta-pair assembly");
        abort();
      }

      oldFFShift = edges[edgeIdx][2];
      oldFRShift = edges[edgeIdx][3];
      oldRFShift = edges[edgeIdx][4];
      oldRRShift = edges[edgeIdx][5];

      /**
       * Test merge here
       */

      /**
       *
       */

      float fShift = edges[edgeIdx][2], rShift = edges[edgeIdx][3];

      bool debug = false;//(idx1 == 1029 && idx2 == 137);

      if (debug)
      {
        int origIdx1 = floatToInt(originalEdges[edgeIdx][0]);
        int origIdx2 = floatToInt(originalEdges[edgeIdx][1]);
        float origFFShift = originalEdges[edgeIdx][2];
        float origFRShift = originalEdges[edgeIdx][3];
        cout << "\nMerging " << idx1 << " (" << origIdx1 << ") and " << idx2
            << " (" << origIdx2 << "), shift = " << origFFShift << ", "
            << origFRShift << "\n";
        cout << "Contig " << origIdx1 << " (" << reversed[origIdx1] << "):\n";
        contigs[origIdx1].output(cout);
        /*
         cout << "\nContig " << origIdx1 << " (reversed):\n";
         Spectrum tempContig1 = contigs[origIdx1];
         tempContig1.reverse(0.0 - AAJumps::massH2O, 0);
         tempContig1.output(cout);
         */
        cout << "\nContig " << origIdx2 << " (" << reversed[origIdx2] << "):\n";
        contigs[origIdx2].output(cout);
        cout << "\nContig " << origIdx2 << " (reversed):\n";
        Spectrum tempContig = contigs[origIdx2];
        tempContig.reverse(0.0 - AAJumps::massH2O, 0);
        tempContig.output(cout);
      }

      res = getConsensusOverlap(idx1,
                                idx2,
                                fShift,
                                rShift,
                                graph,
                                edges,
                                components,
                                root_alignments,
                                false);

      float avgF = res[0], avgR = res[1];

      bool rev = avgR > avgF;
      //rev = (prot_match[idx1][2] != prot_match[idx2][2]) && (reversed[idx1] == reversed[idx2]);
      bool remove = max(avgF, avgR) <= 0 || (checkOverlap && (max(avgF, avgR)
          < minCombScore));
      //remove = false;

      const char* revStr = (rev) ? "1" : "0";

      float shift_use = (rev) ? rShift : fShift;
      float overlap_area = getComponentOverlapArea(idx1,
                                                   idx2,
                                                   shift_use,
                                                   rev,
                                                   graph,
                                                   edges,
                                                   components);

      if (!remove)
      {
        //if (rev)
        //  reverseNode(idx2, graph, edges, components, reversed);
        remove = !tryMergeContigs(idx1,
                                  idx2,
                                  edgeIdx,
                                  graph,
                                  edges,
                                  components,
                                  root_alignments,
                                  rev,
                                  debug);
        //if (rev)
        //  reverseNode(idx2, graph, edges, components, reversed);
      }

      const char* remStr = (remove) ? "1" : "0";

      conFFShift = (rev) ? oldFRShift : oldFFShift;
      conFRShift = (rev) ? oldFFShift : oldFRShift;
      conRFShift = (rev) ? oldRRShift : oldRFShift;
      conRRShift = (rev) ? oldRFShift : oldRRShift;
      bool correct;
      if (haveRes && prot_match[idx1][0] >= 0 && prot_match[idx2][0] >= 0)
      {
        totalEdges += 1.0;
        bool reverse1 = prot_match[idx1][2] == 1;
        bool reverse2 = prot_match[idx2][2] == 1;
        float shift = conFFShift;
        float shiftR = conRRShift;

        bool sameProt = prot_match[idx1][0] == prot_match[idx2][0];
        float protShift1 = getContigShift(idx1, reverse1);
        float protShift2 = getContigShift(idx2, reverse2);
        float protFFShift = protShift2 - protShift1;
        float protRRShift = (protShift1 + contigs[idx1].parentMass)
            - (protShift2 + contigs[idx2].parentMass);
        float temF = protFFShift, temR = protRRShift;
        protFFShift = (reverse1 == reversed[idx1]) ? temF : temR;
        protRRShift = (reverse1 == reversed[idx1]) ? temR : temF;
        correct = isEqual(shift, protFFShift, 2.0) && prot_match[idx1][0]
            == prot_match[idx2][0];
        if (!remove && !sameProt)
        {
          cout
              << "INCORRECTLY MERGING CONTIGS FROM DIFFERENT PROTEINS [contigs "
              << idx1 << " (prot " << prot_match[idx1][0] << ") and " << idx2
              << " (prot " << prot_match[idx2][0] << "), avgOverlap = "
              << max(avgF, avgR) << "]\n";
        }
        else if (!remove && !correct)
        {
          cout << "CONTIG-CONTIG SHIFT DIFFERS FROM EXPECTED BY "
              << abs(protFFShift - shift) << " [contigs " << idx1 << " and "
              << idx2 << " (prot " << prot_match[idx1][0] << "), avgOverlap = "
              << max(avgF, avgR) << "]\n";
          //cout << fShift << " or " << rShift << ": chose " << shift << "\n";
          //cout << prot_match[idx1][2] << ", " << reversed[idx1] << ": " << prot_match[idx2][2] << ", " << reversed[idx2] << "\n";
        }
        protFFShift = (sameProt) ? protFFShift : 0 / zero;
        protRRShift = (sameProt) ? protRRShift : 0 / zero;
        const char* corrStr = (sameProt && correct) ? "1" : "0";
        fprintf(output,
                "%d%s%d%s%d%s%d%s%.2f%s%.2f%s%.2f%s%.2f%s%s%s%s%s%s%s%d%s%d%s%.2f%s%.2f%s%.2f\n",
                prot_match[idx1][0],
                CSV_SEP,
                prot_match[idx2][0],
                CSV_SEP,
                idx1,
                CSV_SEP,
                idx2,
                CSV_SEP,
                shift,
                CSV_SEP,
                shiftR,
                CSV_SEP,
                protFFShift,
                CSV_SEP,
                protRRShift,
                CSV_SEP,
                corrStr,
                CSV_SEP,
                remStr,
                CSV_SEP,
                revStr,
                CSV_SEP,
                prot_match[idx1][2],
                CSV_SEP,
                prot_match[idx2][2],
                CSV_SEP,
                avgF,
                CSV_SEP,
                avgR,
                CSV_SEP,
                overlap_area);
      }
      else if (haveRes)
      {
        fprintf(output,
                "%d%s%d%s%d%s%d%s%.2f%s%.2f%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%.2f%s%.2f%s%.2f\n",
                prot_match[idx1][0],
                CSV_SEP,
                prot_match[idx2][0],
                CSV_SEP,
                idx1,
                CSV_SEP,
                idx2,
                CSV_SEP,
                conFFShift,
                CSV_SEP,
                conRRShift,
                CSV_SEP,
                "",
                CSV_SEP,
                "",
                CSV_SEP,
                "",
                CSV_SEP,
                remStr,
                CSV_SEP,
                revStr,
                CSV_SEP,
                "",
                CSV_SEP,
                "",
                CSV_SEP,
                avgF,
                CSV_SEP,
                avgR,
                CSV_SEP,
                overlap_area);
      }
      else
      {
        const char* remStr = (remove) ? "1" : "0";
        const char* revStr = (rev) ? "1" : "0";
        float shift = conFFShift;
        float shiftR = conRRShift;
        fprintf(output,
                "%d%s%d%s%.2f%s%.2f%s%s%s%s%s%.2f%s%.2f%s%.2f\n",
                idx1,
                CSV_SEP,
                idx2,
                CSV_SEP,
                shift,
                CSV_SEP,
                shiftR,
                CSV_SEP,
                remStr,
                CSV_SEP,
                revStr,
                CSV_SEP,
                avgF,
                CSV_SEP,
                avgR,
                CSV_SEP,
                overlap_area);
      }
      if (haveRes && correct)
        goodEdges += 1.0;
      if (!remove)
        edgesAdded += 1.0;
      if (haveRes && correct && !remove)
        goodEdgesAdded += 1.0;

      if (remove)
      {
        erasedEdges.insert(edgeIdx);
        if (graph[idx1][idx2] == edgeIdx)
        {
          erasedEdges.insert(graph[idx2][idx1]);
          graph[idx1].erase(idx2);
          graph[idx2].erase(idx1);
        }
        scoredEdgeIt++;
        if (debug)
          break;
        continue;
      }

      if (mergingMatePairs)
      {
        countmatep += edgeNum[edgeIdx];
        countmerg += 1;
      }
      else
      {
        countcont += 1;
      }

      edgeTo[0] = (float)idx1;
      edgeTo[1] = (float)idx2;
      edgeTo[2] = conFFShift;
      edgeTo[3] = conFRShift;
      edgeTo[4] = conRFShift;
      edgeTo[5] = conRRShift;
      edgeFrom[0] = (float)idx2;
      edgeFrom[1] = (float)idx1;
      edgeFrom[2] = 0.0 - conFFShift;
      edgeFrom[3] = 0.0 - conRFShift;
      edgeFrom[4] = 0.0 - conFRShift;
      edgeFrom[5] = 0.0 - conRRShift;
      edges[edgeIdx] = edgeTo;
      edges[graph[idx2][idx1]] = edgeFrom;

      erasedEdges.insert(graph[idx2][idx1]);
      graph[idx2].erase(idx1);

      if (rev)
        reverseNode(idx2, graph, edges, components, reversed);

      int most_matched_peaks = 0;
      TwoValues<int> best_node_pair;
      best_node_pair[0] = idx1;
      best_node_pair[1] = idx2;

      for (nodeIt = components[idx1].begin(); nodeIt != components[idx1].end(); nodeIt++)
      {
        float shift1 = (*nodeIt == idx1) ? 0 : edges[graph[idx1][*nodeIt]][2];
        for (nodeIt2 = components[idx2].begin(); nodeIt2
            != components[idx2].end(); nodeIt2++)
        {
          if (*nodeIt2 == *nodeIt)
            continue;
          float shift2 = (*nodeIt2 == idx2) ? 0
              : edges[graph[idx2][*nodeIt2]][2];
          float locFShift = shift_use - shift1 + shift2;
          TwoValues<int> mp = getContigOverlapPeaks(oriented[*nodeIt],
                                                    oriented[*nodeIt2],
                                                    locFShift,
                                                    0,
                                                    parentMassTol,
                                                    minContigOverlap,
                                                    debug);
          if (mp[0] > most_matched_peaks)
          {
            best_node_pair[0] = *nodeIt;
            best_node_pair[1] = *nodeIt2;
            most_matched_peaks = mp[0];
          }
        }
      }
      if (most_matched_peaks == 0)
      {
        cout << "\nCOULD NOT FIND OVERLAPPING CONTIGS BETWEEN COMONENTS "
            << idx1 << " AND " << idx2 << "\n";
      }
      if (debug)
      {
        cout << "\nbest node pair: " << best_node_pair[0] << " and "
            << best_node_pair[1] << "(" << most_matched_peaks
            << " matching peaks)\n";
      }
      //list<TwoValues<int> > idx2_aligns = root_alignments[idx2];
      //idx2_aligns.push_back(best_node_pair);
      root_alignments[idx1].push_back(best_node_pair);
      for (list<TwoValues<int> >::iterator alIt = root_alignments[idx2].begin(); alIt
          != root_alignments[idx2].end(); alIt++)
      {
        root_alignments[idx1].push_back(*alIt);
      }

      /*if (idx1 == 3247 || idx1 == 1029 || idx1 == 137 || idx1 == 1557 || idx1 == 1293) {
       cout << "idx1 = " << idx1 << "\n";
       cout << "idx2 = " << idx2 << "\n";
       cout << "components[idx1].size() = " << components[idx1].size() << "\n";
       cout << "components[idx2].size() = " << components[idx2].size() << "\n";
       cout << "root_alignments[idx1].size() = " << root_alignments[idx1].size() << "\n";
       cout << "root_alignments[idx2].size() = " << root_alignments[idx2].size() << "\n";
       }*/

      components[idx1].insert(idx2);

      for (nodeIt = components[idx2].begin(); nodeIt != components[idx2].end(); nodeIt++)
      {
        int node2 = *nodeIt;

        if (node2 == idx2)
          continue;

        int edge2 = graph[idx2][node2];

        FFShift = conFFShift + edges[edge2][2];
        FRShift = conFFShift + edges[edge2][3];
        RFShift = conRRShift + edges[edge2][4];
        RRShift = conRRShift + edges[edge2][5];

        edgeTo[0] = (float)idx1;
        edgeTo[1] = (float)node2;
        edgeTo[2] = FFShift;
        edgeTo[3] = FRShift;
        edgeTo[4] = RFShift;
        edgeTo[5] = RRShift;

        components[idx1].insert(node2);
        edges[edge2] = edgeTo;
        graph[idx1][node2] = edge2;
      }

      Spectrum consensusIdx1;
      pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> , vector<
          double> > > > abIdx1;
      TwoValues<float> extraShift1 = mergeContigs(idx1,
                                                  graph,
                                                  edges,
                                                  components,
                                                  root_alignments,
                                                  consensusIdx1,
                                                  abIdx1,
                                                  debug);

      consensus[idx1] = consensusIdx1;
      extraShifts[idx1] = extraShift1;
      parent_abinfo[idx1] = abIdx1;
      contigsMerged[idx1] = true;
      contigsMerged[idx2] = true;

      map<int, map<int, int> > graphSep = graph;
      vector<vector<float> > edgesSep = edges;

      for (mapIt = graph[idx2].begin(); mapIt != graph[idx2].end(); mapIt++)
      {
        int node2 = mapIt->first;
        int edge2 = mapIt->second;
        if (components[idx1].count(node2) > 0)
          continue;

        FFShift = conFFShift + edges[edge2][2];
        FRShift = conFFShift + edges[edge2][3];
        RFShift = conRRShift + edges[edge2][4];
        RRShift = conRRShift + edges[edge2][5];

        edgeTo[0] = (float)idx1;
        edgeTo[1] = (float)node2;
        edgeTo[2] = FFShift;
        edgeTo[3] = FRShift;
        edgeTo[4] = RFShift;
        edgeTo[5] = RRShift;
        edgeFrom[0] = (float)node2;
        edgeFrom[1] = (float)idx1;
        edgeFrom[2] = 0.0 - FFShift;
        edgeFrom[3] = 0.0 - RFShift;
        edgeFrom[4] = 0.0 - FRShift;
        edgeFrom[5] = 0.0 - RRShift;

        if (graph[idx1].count(node2) > 0)
        {
          bool edge1OK = graph[idx1][node2] <= maxContigEdgeIdx;
          bool edge2OK = graph[idx2][node2] <= maxContigEdgeIdx;

          float prevShiftF = edges[graph[idx1][node2]][2];
          float prevShiftR = edges[graph[idx1][node2]][3];
          float prevShiftRF = edges[graph[idx1][node2]][4];
          float prevShiftRR = edges[graph[idx1][node2]][5];
          if (!isEqual(prevShiftF, FFShift, 2.0) && !isEqual(prevShiftR,
                                                             FRShift,
                                                             2.0)
              && !isEqual(prevShiftRF, RFShift, 2.0) && !isEqual(prevShiftRR,
                                                                 RRShift,
                                                                 2.0))
          {
            bool deb = false;//(idx1 == 253 && idx2 == 523);
            //if (deb) {cout << "idx1 = " << idx1 << ", idx2 = " << idx2 << ", node2 = " << node2 << "\n"; cout.flush();}
            res1 = getConsensusOverlap(idx1,
                                       node2,
                                       prevShiftF,
                                       prevShiftR,
                                       graphSep,
                                       edgesSep,
                                       components,
                                       root_alignments,
                                       deb);
            res2 = getConsensusOverlap(idx1,
                                       node2,
                                       FFShift,
                                       FRShift,
                                       graphSep,
                                       edgesSep,
                                       components,
                                       root_alignments,
                                       deb);
            float score1 = max(res1[0], res1[1]);
            float score2 = max(res2[0], res2[1]);

            //cout << "Conflicting edges between " << node2 << " and component (" << idx1 << "," << idx2 << ")\n"; cout << idx1 << " -- " << node2 << ": FFShift = " << prevShiftF << ", FRShift = " << prevShiftR << ", RFShift = " << prevShiftRF << ", RRShift = " << prevShiftRR << "\n" << idx1 << " -- " << idx2 << " -- " << node2 << ": FFShift = " << FFShift << ", FRShift = " << FRShift << ", RFShift = " << RFShift << ", RRShift = " << RRShift << "\n";
            /*
             bool r1 = prot_match[idx1][2] == 1;
             bool r2 = prot_match[node2][2] == 1;
             float ps1 = getContigShift(idx1, r1);
             float ps2f = getContigShift(node2, r2);
             cout << "correct shift = " << ps2f - ps1 << "\n";
             */
            if ((edge1OK && !edge2OK) || (score1 >= score2))
            {
              //cout << "Removing edge from " << idx2 << " to " << node2 << "\n";
              erasedEdges.insert(graph[idx2][node2]);
              erasedEdges.insert(graph[node2][idx2]);
              edgeRef[graph[idx2][node2]] = graph[idx1][node2];
              edgeRef[graph[node2][idx2]] = graph[node2][idx1];
              graph[node2].erase(idx2);
            }
            else
            {//score2 > minSecondaryRatio &&
              //cout << "Removing edge from " << idx1 << " to " << node2 << "\n";
              erasedEdges.insert(graph[idx1][node2]);
              erasedEdges.insert(graph[node2][idx1]);
              edgeRef[graph[idx1][node2]] = graph[idx2][node2];
              edgeRef[graph[node2][idx1]] = graph[node2][idx2];
              //if (components[idx2].count(node2) > 0) idx1sToAdd.insert(node2);//components[idx1].insert(node2);
              int edgeToNode = graph[idx2][node2];
              int edgeFromNode = graph[node2][idx2];
              edges[edgeToNode] = edgeTo;
              edges[edgeFromNode] = edgeFrom;
              graph[node2].erase(idx2);
              graph[node2][idx1] = edgeFromNode;
              graph[idx1][node2] = edgeToNode;
            }/* else {
             //cout << "Removing edge from " << idx2 << " to " << node2 << "\n";
             //cout << "Removing edge from " << idx1 << " to " << node2 << "\n";
             erasedEdges.insert(graph[idx2][node2]);
             erasedEdges.insert(graph[node2][idx2]);
             graph[node2].erase(idx2);
             erasedEdges.insert(graph[idx1][node2]);
             erasedEdges.insert(graph[node2][idx1]);
             graph[idx1].erase(node2);
             graph[node2].erase(idx1);
             }*/
          }
          else
          {
            //cout << "Agreeable edges between " << node2 << " and component (" << idx1 << "," << idx2 << ")\n";
            /*if (components[idx2].count(node2) > 0) {
             //cout << idx2 << " already in component.\n";
             edges[graph[idx2][node2]] = edgeTo;
             idx1sToAdd.insert(node2);//components[idx1].insert(node2);
             }
             else if (components[idx1].count(node2) > 0) {
             //cout << idx1 << " already in component.\n";
             erasedEdges.insert(graph[idx2][node2]);
             erasedEdges.insert(graph[node2][idx2]);
             graph[node2].erase(idx2);
             }*/
            if (edge2OK && !edge1OK)
            {
              //cout << "choosing edge " << edge2 << " 2( " << maxContigEdgeIdx << " ) ( " << edgeTo[2] << " )\n";
              erasedEdges.insert(graph[idx1][node2]);
              erasedEdges.insert(graph[node2][idx1]);
              edgeRef[graph[idx1][node2]] = graph[idx2][node2];
              edgeRef[graph[node2][idx1]] = graph[node2][idx2];
              //if (components[idx2].count(node2) > 0) idx1sToAdd.insert(node2);//components[idx1].insert(node2);
              int edgeToNode = graph[idx2][node2];
              int edgeFromNode = graph[node2][idx2];
              edges[edgeToNode] = edgeTo;
              edges[edgeFromNode] = edgeFrom;
              graph[node2].erase(idx2);
              graph[node2][idx1] = edgeFromNode;
              graph[idx1][node2] = edgeToNode;
            }
            else
            {
              //cout << "choosing edge " << graph[idx1][node2] << " 1( " << maxContigEdgeIdx << " ) ( " << edges[graph[idx1][node2]][2] << " )\n";
              erasedEdges.insert(graph[idx2][node2]);
              erasedEdges.insert(graph[node2][idx2]);
              edgeRef[graph[idx2][node2]] = graph[idx1][node2];
              edgeRef[graph[node2][idx2]] = graph[node2][idx1];
              graph[node2].erase(idx2);
            }
          }

        }
        else
        {
          edges[edge2] = edgeTo;
          graph[idx1][node2] = edge2;

          int node2Edge = graph[node2][idx2];
          edges[node2Edge] = edgeFrom;
          graph[node2].erase(idx2);
          graph[node2][idx1] = node2Edge;
        }
      }
      //cout << "erasing\n";
      //cout.flush();
      graph.erase(idx2);
      components[idx2].clear();
      consensus[idx2].resize(0);
      parent_abinfo.erase(idx2);

      //outputGraph(cout, graph, edges, &erasedEdges, 0);
      //break;
      //if (idx1 == 90 && idx2 == 153) {cout << "\n" << idx1 << " and " << idx2 << " after: " << edgeIdx << "\n\n"; outputGraph(cout, graph, edges, &erasedEdges, &components);}
      if (scoredEdgeIt != scoredEdges.end())
        scoredEdgeIt++;
      scoredEdges.erase(scoredEdges.begin(), scoredEdgeIt);

      recomputeEdges(idx1,
                     graph,
                     edges,
                     components,
                     erasedEdges,
                     edgeRef,
                     root_alignments,
                     scoredEdges);

      scoredEdgeIt = scoredEdges.begin();

      if (debug)
        break;
    }

    set<int> countNodes;
    vector<int> compSize(components.size());
    for (int i = 0; i < components.size(); i++)
    {
      compSize[i] = components[i].size();
      if (components[i].size() == 0)
      {
        consensus[i].resize(0);
        continue;
      }
      if (components[i].size() < 1)
        continue;
      set<int> seen;
      for (nodeIt = components[i].begin(); nodeIt != components[i].end(); nodeIt++)
      {
        seen. insert(*nodeIt);
      }
      //c  ompSize[i] = seen.size();
      if (mergeType == 1 && seen.size() < minComp)
      {
        components[i].clear();
        continue;
      }
      countNodes.insert(i);
      for (nodeIt = components[i].begin(); nodeIt != components[i].end(); nodeIt++)
      {
        countNodes.insert(*nodeIt);
      }
    }
    //printf("\nCombined %s %.2f good edges (%.0f/%.0f)\n", "%", (100.0*goodEdges)/edgesAdded, goodEdges, edgesAdded);
    outputComponents(output,
                     graph,
                     edges,
                     components,
                     oriented,
                     reversed,
                     false,
                     false,
                     0,
                     true);
    fclose(output);

    if (mergeType == 1)
      printf("\nStatistics for merged contigs after contig/mate-pair shifts used:\n");
    else
      printf("\nStatistics for merged contigs after contig/contig shifts used:\n");

    outputComponents(0,
                     graph,
                     edges,
                     components,
                     oriented,
                     reversed,
                     true,
                     false,
                     &countNodes);

    printf("\nStatistics for contigs before any shifts used:\n");
    outputComponents(0,
                     graphBefMerg,
                     edgesBefMerg,
                     componentsBefMerg,
                     myContigsBefMerg,
                     reversedBefMerg,
                     true,
                     false,
                     &countNodes,
                     true);

    printf("\nUsed %d contig/contig alignments and %d contig/mate-pair alignments (%d contigs merged with mate-pairs)\n\n",
           countcont,
           countmatep,
           countmerg);

    cout << "merging_params->min_component_size: "
        << merging_params->min_component_size << "\n";
    for (int i = 0; i < consensus.size(); i++)
    {
      if (compSize[i] < merging_params->min_component_size)
      {
        consensus[i].resize(0);
      }
    }

    merging_params->consensus.resize(consensus.size());
    for (int i = 0; i < consensus.size(); i++)
    {
      merging_params->consensus[i] = consensus[i];
    }

    SpecSet stars = *(merging_params->star_spectra);
    stars.addZPMpeaks(peakTol, 0, false);

    /*
     unsigned int idxcheck = 11;
     cout << "\n" <<idxcheck << ":\n";
     for (int i = 0; i < parent_abinfo[idxcheck].second.size(); i++)
     {
     cout << " -- " << i << " -> ";
     for (int j = 0; j < parent_abinfo[idxcheck].second[i].first.size(); j++)
     {
     cout << "(" << parent_abinfo[idxcheck].second[i].first[j] << ", " << reversed[parent_abinfo[idxcheck].second[i].first[j]] << ", " << parent_abinfo[idxcheck].second[i].second[j] << "); ";
     }
     cout << endl;
     }
     */

    Combine_abinfo_v1_0(contig_abinfo,
                        reversed,
                        stars,
                        contigsMerged,
                        contigs,
                        parent_abinfo,
                        consensus_abinfo);

    merging_params->out_reversed.resize(reversed.size());

    for (int i = 0; i < reversed.size(); i++)
      merging_params->out_reversed[i] = (reversed[i]) ? 1 : 0;

    merging_params->out_abinfo.clear();
    for (unsigned int i = 0; i < consensus.size(); i++)
    {
      merging_params->out_abinfo[i] = consensus_abinfo[i];
    }

  }

  bool CombineContigs::saveComponents(const char* outcomponents)
  {
    FILE* compsOut = fopen(outcomponents, "wb");
    if (compsOut == NULL)
    {
      cerr << "ERROR: Cannot write to file " << outcomponents << endl;
      return false;
    }

    int count = 0;
    fprintf(compsOut,
            "Component Idx%sContig Idx2%sContig Root Idx%sContig2 Shift\n",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
    for (int i = 0; i < components.size(); i++)
    {
      if (components[i].size() < 2)
        continue;
      set<int> seen;
      map<int, set<int> > seenNodes;
      set<int> snode;
      for (set<int>::iterator nodeIt = components[i].begin(); nodeIt
          != components[i].end(); nodeIt++)
      {
        if (seen.count(*nodeIt) > 0)
          continue;
        seen.insert(*nodeIt);
        float FFShift = (i == *nodeIt) ? 0.0 : edges[graph[i][*nodeIt]][2];
        fprintf(compsOut,
                "%d%s%d%s%d%s%.1f\n",
                count,
                CSV_SEP,
                *nodeIt,
                CSV_SEP,
                i,
                CSV_SEP,
                FFShift);
      }
      count++;
    }
    fclose(compsOut);
    return true;
  }

  void CombineContigs::init()
  {
    merging_params->consensus = contigs;
    merging_params->out_abinfo = contig_abinfo;
    (merging_params->out_reversed).resize(contigs.size());
    for (unsigned int i = 0; i < (merging_params->out_reversed).size(); i++)
    {
      (merging_params->out_reversed)[i] = 0;
    }
    consensus = contigs;
    oriented = contigs;
    consensus_abinfo = contig_abinfo;
    extraShifts.resize(contigs.size());
    for (unsigned int i = 0; i < extraShifts.size(); i++)
    {
      extraShifts[i][0] = 0;
      extraShifts[i][1] = 0;
    }
    numCorrectPairs = 0;
  }

  void CombineContigs::outputGraph(ostream& output,
                                   map<int, map<int, int> >& _graph,
                                   vector<vector<float> >& _edges,
                                   set<int>* erasedEdges,
                                   vector<set<int> >* _components)
  {
    map<int, int>::iterator mapIt;
    set<int>::iterator nodeIt;
    output << "Nodes (edge):\n";
    for (int i = 0; i < contigs.size(); i++)
    {
      if (_graph.count(i) == 0)
        continue;
      output << i << " -> ";
      for (mapIt = _graph[i].begin(); mapIt != _graph[i].end(); mapIt++)
      {
        output << mapIt->first << "(" << mapIt->second << "), ";
      }
      output << "\n";
    }
    output << "\nEdges:\n";
    for (int i = 0; i < _edges.size(); i++)
    {
      if (erasedEdges != 0 && (*erasedEdges).count(i) > 0)
        continue;
      output << i << ": ";
      for (int j = 0; j < _edges[i].size(); j++)
      {
        output << _edges[i][j] << ", ";
      }
      output << "\n";
    }
    if (_components != 0)
    {
      output << "\nComponents:\n";
      for (int i = 0; i < (*_components).size(); i++)
      {
        if ((*_components)[i].size() < 2)
          continue;
        output << i << ": ";
        for (nodeIt = (*_components)[i].begin(); nodeIt
            != (*_components)[i].end(); nodeIt++)
        {
          output << *nodeIt;
          if (*nodeIt == i)
            cout << "(0), ";
          else
            cout << "(" << _edges[_graph[i][*nodeIt]][2] << ","
                << _edges[_graph[i][*nodeIt]][5] << "), ";
        }
        output << "\n";
      }
    }
  }

  float CombineContigs::getComponentOverlapArea(int idx1,
                                                int idx2,
                                                float shift_use,
                                                bool rev,
                                                map<int, map<int, int> >& _graph,
                                                vector<vector<float> >& _edges,
                                                vector<set<int> >& _components)
  {

    float idx1_left_bound = 0;
    float idx1_right_bound = oriented[idx1].parentMass;

    float idx2_left_bound = shift_use;
    float idx2_right_bound = shift_use + oriented[idx2].parentMass;

    for (set<int>::iterator nodeIt = _components[idx1].begin(); nodeIt
        != _components[idx1].end(); nodeIt++)
    {
      if (*nodeIt == idx1)
        continue;
      float ffshift = _edges[_graph[idx1][*nodeIt]][2];

      idx1_left_bound = min(idx1_left_bound, ffshift);
      idx1_right_bound = max(idx1_right_bound, ffshift
          + oriented[*nodeIt].parentMass);
    }
    for (set<int>::iterator nodeIt = _components[idx2].begin(); nodeIt
        != _components[idx2].end(); nodeIt++)
    {
      if (*nodeIt == idx2)
        continue;
      float ffshift = (rev) ? _edges[_graph[idx2][*nodeIt]][5]
          : _edges[_graph[idx2][*nodeIt]][2];

      idx2_left_bound = min(idx2_left_bound, ffshift + shift_use);
      idx2_right_bound = max(idx2_right_bound, ffshift + shift_use
          + oriented[*nodeIt].parentMass);
    }

    float left_bound = max(idx1_left_bound, idx2_left_bound);
    float right_bound = min(idx1_right_bound, idx2_right_bound);

    return right_bound - left_bound;
  }

  void CombineContigs::outputComponents(FILE* output,
                                        map<int, map<int, int> >& _graph,
                                        vector<vector<float> >& _edges,
                                        vector<set<int> >& _components,
                                        SpecSet& orientedContigs,
                                        vector<bool>& _reversed,
                                        bool outputStats,
                                        bool printComps,
                                        set<int>* countNodes,
                                        bool justContigs)
  {

    bool checkPossible = false;
    if (outputStats && haveRes && numCorrectPairs == 0)
    {
      checkPossible = true;
      int possible = 0;
      map<int, set<int> > seen;
      set<int> sit;
      for (int i = 0; i < contigs.size(); i++)
      {
        for (int j = 0; j < contigs.size(); j++)
        {
          int idx1 = i, idx2 = j;
          if (idx1 == idx2)
            continue;
          if (idxNotUsed.count(prot_match[idx1][0]) > 0
              || idxNotUsed.count(prot_match[idx2][0]) > 0)
            continue;
          if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0
              || prot_match[idx1][0] != prot_match[idx2][0])
            continue;
          if (seen.count(i) > 0 && seen[i].count(j) > 0)
            continue;
          if (seen.count(j) > 0 && seen[j].count(i) > 0)
            continue;
          if (seen.count(i) > 0)
            seen[i].insert(j);
          else
          {
            sit.clear();
            sit.insert(j);
            seen[i] = sit;
          }
          bool reverse1 = prot_match[idx1][2] == 1;
          bool reverse2 = prot_match[idx2][2] == 1;
          float protShift1 = getContigShift(idx1, reverse1);
          float protShift2 = getContigShift(idx2, reverse2);
          float protFFShift = protShift2 - protShift1;
          //float protRRShift = (protShift1+contigs[idx1].parentMass) - (protShift2+contigs[idx2].parentMass);
          //protFFShift = (reverse1) ? protRRShift : protFFShift;
          //if (protFFShift < 0 - contigs[idx2].parentMass || protFFShift > contigs[idx1].parentMass) continue;
          Spectrum contig1 = contigs[idx1];
          if (reverse1)
          {
            contig1.reverse(0.0 - AAJumps::massH2O, 0);
          }
          Spectrum contig2 = contigs[idx2];
          if (reverse2)
          {
            contig2.reverse(0.0 - AAJumps::massH2O, 0);
          }
          TwoValues<int> mp = getContigOverlapPeaks(contig1,
                                                    contig2,
                                                    protFFShift,
                                                    contigs[idx1].parentMass,
                                                    0.5,
                                                    0,
                                                    false);

          if (mp[0] < 6)
            continue;

          if (correctPairs.count(i) > 0)
            correctPairs[i].insert(j);
          else
          {
            sit.clear();
            sit.insert(j);
            correctPairs[i] = sit;
          }

          possible++;
        }
      }
      numCorrectPairs = possible;
    }
    bool addTo = seenBef.size() == 0;
    float totalLen = 0, numLen = 0, maxLen = 0, totalCover = 0, maxCover = 0;

    set<int>::iterator nodeIt;
    set<int>::iterator nodeIt2;
    if (haveRes)
    {
      if (output != 0)
        fprintf(output,
                "\nComponent Idx%sProtein Idx%sContig Root Idx%sContig Idx2%sExperimental FF Shift%sExperimental RR Shift%sExpected FF Shift%sExpected RR Shift%sHas Good Edge%sComponent Length (Da)%sComponent Length (AA)\n",
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP);
    }
    else
    {
      if (output != 0)
        fprintf(output,
                "\nComponent Idx%sContig Root Idx%sContig Idx2%sExperimental FF Shift%sExperimental RR Shift%sComponent Length (Da)\n",
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP);
    }
    map<int, int> componentRef;
    for (int i = 0; i < contigs.size(); i++)
    {
      componentRef[i] = -1;
    }
    map<int, set<int> > seenNodes;
    set<int> snode;
    float goodEdges = 0, edgesAdded = 0;
    int numCompsWSize = 0;
    int compIdx = 0;
    int minCompSize = 5;
    int good_aligned = 0;
    if (haveRes && printComps)
    {
      printf("\nComponent Index - Length (AA) - Protein Index - Protein Header\n");
    }
    for (int i = 0; i < _components.size(); i++)
    {
      //if (_components[i].size() == 0) continue;
      if (countNodes != 0 && (*countNodes).count(i) == 0)
        continue;

      if (_components[i].size() <= 1)
        continue;

      if (_components[i].size() >= minCompSize)
      {
        numCompsWSize++;
      }

      float len = 0, AAlen = 0, coverAA = 0;
      vector<int> protids(fasta.size());
      for (int k = 0; k < protids.size(); k++)
        protids[k] = 0;

      if (haveRes && prot_match[i][0] >= 0
          && idxNotUsed.count(prot_match[i][0]) == 0)
      {
        coverAA = abs(overlaps[i][overlaps[i].size() - 1][1]
            - overlaps[i][0][1]) + 1;
        AAlen = max(AAlen, coverAA);
        len = contigs[i].parentMass;
      }
      if (!haveRes)
      {
        len = contigs[i].parentMass;
      }
      for (nodeIt = _components[i].begin(); nodeIt != _components[i].end(); nodeIt++)
      {
        componentRef[*nodeIt] = compIdx;
        float FFShift = (i == *nodeIt) ? 0.0 : _edges[_graph[i][*nodeIt]][2];
        float RRShift = (i == *nodeIt) ? 0.0 : _edges[_graph[i][*nodeIt]][5];
        if (haveRes && prot_match[*nodeIt][0] >= 0)
        {
          protids[prot_match[*nodeIt][0]] += 1;
        }
        if (haveRes)
        {
          for (nodeIt2 = _components[i].begin(); nodeIt2
              != _components[i].end(); nodeIt2++)
          {
            int idx1 = *nodeIt, idx2 = *nodeIt2;
            if (idx1 == idx2)
              continue;
            if (idxNotUsed.count(prot_match[idx1][0]) > 0
                || idxNotUsed.count(prot_match[idx2][0]) > 0)
            {
              continue;
            }
            if ((seenNodes.count(idx1) > 0 && seenNodes[idx1].count(idx2) > 0)
                || (seenNodes.count(idx2) > 0 && seenNodes[idx2].count(idx1)
                    > 0))
            {
              continue;
            }
            if (seenNodes.count(idx1) == 0)
            {
              snode.clear();
              snode.insert(idx2);
              seenNodes[idx1] = snode;
            }
            else
            {
              seenNodes[idx1].insert(idx2);
            }

            float shift = abs(_edges[_graph[i][idx1]][2]
                - _edges[_graph[i][idx2]][2]);
            if (_edges[_graph[i][idx1]][2] + contigs[idx1].parentMass
                > _edges[_graph[i][idx2]][2] + contigs[idx2].parentMass)
            {
              shift += contigs[idx1].parentMass;
            }
            else
            {
              shift += contigs[idx2].parentMass;
            }
            len = max(len, shift);

            if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0)
              continue;

            float FFShift2 = (i == idx2) ? 0.0 : _edges[_graph[i][idx2]][2];
            float RRShift2 = (i == idx2) ? 0.0 : _edges[_graph[i][idx2]][5];

            coverAA = overlaps[idx2][overlaps[idx2].size() - 1][1]
                - overlaps[idx2][0][1] + 1.0;
            coverAA = max(coverAA, overlaps[idx1][overlaps[idx1].size() - 1][1]
                - overlaps[idx1][0][1] + (float)1.0);

            if (prot_match[idx1][0] == prot_match[idx2][0]
                && (validShiftMod(FFShift2 - FFShift,
                                  idx1,
                                  idx2,
                                  _reversed[idx2]) || validShiftMod(RRShift2
                    - RRShift, idx1, idx2, !_reversed[idx2])))
            {
              goodEdges += 1.0;

              float maxPos = max(overlaps[idx2][overlaps[idx2].size() - 1][1],
                                 overlaps[idx1][overlaps[idx1].size() - 1][1]);
              float minPos = min(overlaps[idx2][0][1], overlaps[idx1][0][1]);
              coverAA = max(coverAA, maxPos - minPos + (float)1.0);
              TwoValues<int> mp =
                  getContigOverlapPeaks(orientedContigs[idx1],
                                        orientedContigs[idx2],
                                        FFShift2 - FFShift,
                                        contigs[idx1].parentMass,
                                        1.0,
                                        0,
                                        false);
              if (mp[0] >= 6 && ((correctPairs.count(idx1) > 0
                  && correctPairs[idx1].count(idx2) > 0)
                  || (correctPairs.count(idx2) > 0
                      && correctPairs[idx2].count(idx1) > 0)))
              {
                if (addTo)
                {
                  if (seenBef.count(idx1) > 0)
                    seenBef[idx1].insert(idx2);
                  else
                  {
                    snode.clear();
                    snode.insert(idx2);
                    seenBef[idx1] = snode;
                  }
                }
                else
                {
                  if ((seenBef.count(idx1) == 0 || seenBef[idx1].count(idx2)
                      == 0) && (seenBef.count(idx2) == 0
                      || seenBef[idx2].count(idx1) == 0))
                  {
                    cout << "\nFOUND IT: " << idx1 << " " << idx2 << "\n";
                  }
                }

                good_aligned++;
              }
            }

            AAlen = max(AAlen, coverAA);
            edgesAdded += 1.0;
          }
        }
        else
        {
          for (nodeIt2 = _components[i].begin(); nodeIt2
              != _components[i].end(); nodeIt2++)
          {
            int idx1 = *nodeIt, idx2 = *nodeIt2;
            if (idx1 == idx2)
              continue;
            float shift = abs(_edges[_graph[i][idx1]][2]
                - _edges[_graph[i][idx2]][2]);
            if (_edges[_graph[i][idx1]][2] + contigs[idx1].parentMass
                > _edges[_graph[i][idx2]][2] + contigs[idx2].parentMass)
            {
              shift += contigs[idx1].parentMass;
            }
            else
            {
              shift += contigs[idx2].parentMass;
            }
            len = max(len, shift);
          }
        }
      }
      if (_components[i].size() >= minCompSize)
      {
        totalLen += len;
        numLen += 1.0;
        maxLen = max(maxLen, len);
      }

      if (haveRes)
      {
        int bestprot = -1, maxCount = 0;
        for (int k = 0; k < protids.size(); k++)
        {
          if (protids[k] > maxCount)
          {
            bestprot = k;
            maxCount = protids[k];
          }
        }
        if (_components[i].size() >= minCompSize)
        {
          maxCover = max(maxCover, AAlen);
          totalCover += AAlen;
        }

        const char* protident = "none";
        if (bestprot >= 0)
          protident = fasta.getID(bestprot);
        if (printComps)
          printf("%d - %.0f - %d - %s\n", compIdx, AAlen, bestprot, protident);
      }
      for (nodeIt = _components[i].begin(); nodeIt != _components[i].end(); nodeIt++)
      {
        float FFShift = (i == *nodeIt) ? 0.0 : _edges[_graph[i][*nodeIt]][2];
        float RRShift = (i == *nodeIt) ? 0.0 : _edges[_graph[i][*nodeIt]][5];
        if (haveRes)
        {
          string protFFStr = "";
          string protRRStr = "";
          const char* eqStr = "";
          if (prot_match[*nodeIt][0] >= 0 && prot_match[i][0] >= 0)
          {
            bool reverse1 = prot_match[i][2] == 1;
            bool reverse2 = prot_match[*nodeIt][2] == 1;

            bool sameProt = prot_match[i][0] == prot_match[*nodeIt][0];
            float protShift1 = getContigShift(i, reverse1);
            float protShift2 = getContigShift(*nodeIt, reverse2);
            float protFFShift = protShift2 - protShift1;
            float protRRShift = (protShift1 + contigs[i].parentMass)
                - (protShift2 + contigs[*nodeIt].parentMass);
            float temF = protFFShift, temR = protRRShift;
            protFFShift = (reverse1 == _reversed[i]) ? temF : temR;
            protRRShift = (reverse1 == _reversed[*nodeIt]) ? temR : temF;
            if (sameProt)
            {
              protFFStr = parseFloat(protFFShift, 3);
              protRRStr = parseFloat(protRRShift, 3);
              eqStr = (isEqual(protFFShift, FFShift, 58.0)) ? (char*)"1"
                  : (char*)"0";
            }
            else
            {
              protFFStr = "nan";
              protRRStr = "nan";
              eqStr = "0";
            }
          }
          if (output != 0)
            fprintf(output,
                    "%d%s%d%s%d%s%d%s%.1f%s%.1f%s%s%s%s%s%s%s%.1f%s%.0f\n",
                    compIdx,
                    CSV_SEP,
                    prot_match[*nodeIt][0],
                    CSV_SEP,
                    i,
                    CSV_SEP,
                    *nodeIt,
                    CSV_SEP,
                    FFShift,
                    CSV_SEP,
                    RRShift,
                    CSV_SEP,
                    protFFStr.c_str(),
                    CSV_SEP,
                    protRRStr.c_str(),
                    CSV_SEP,
                    eqStr,
                    CSV_SEP,
                    len,
                    CSV_SEP,
                    AAlen);
        }
        else
        {
          if (output != 0)
            fprintf(output,
                    "%d%s%d%s%d%s%.1f%s%.1f%s%.1f\n",
                    compIdx,
                    CSV_SEP,
                    i,
                    CSV_SEP,
                    *nodeIt,
                    CSV_SEP,
                    FFShift,
                    CSV_SEP,
                    RRShift,
                    CSV_SEP,
                    len);
        }
      }
      compIdx++;
    }

    if (addTo)
    {
      for (map<int, set<int> >::iterator possIt = correctPairs.begin(); possIt
          != correctPairs.end(); possIt++)
      {
        int idx1 = possIt->first;
        set<int> idx2s = possIt->second;
        for (set<int>::iterator secIt = idx2s.begin(); secIt != idx2s.end(); secIt++)
        {
          int idx2 = *secIt;
          if (componentRef[idx1] != componentRef[idx2] && componentRef[idx1]
              >= 0 && componentRef[idx2] >= 0)
          {// && (seenNodes.count(idx1) == 0 || seenNodes[idx1].count(idx2) == 0) && (seenNodes.count(idx2) == 0 || seenNodes[idx2].count(idx1) == 0)) {
            //cout << idx1 << "-" << componentRef[idx1] << ", " << idx2 << "-" << componentRef[idx2] << ": Missed true positive contig/contig connection w/ >= 6 matching peaks between " << idx1 << " (protein " << prot_match[idx1][0] << ") and " << idx2 << " (protein " << prot_match[idx2][0] << ")\n";
          }
        }
      }
    }
    /*
     if (checkPossible) {
     for (map<int, set<int> >::iterator possIt = possiblePairs.begin(); possIt != possiblePairs.end(); possIt ++) {
     int idx1 = possIt->first;
     set<int> idx2s = possIt->second;
     for (set<int>::iterator secIt = idx2s.begin(); secIt != idx2s.end(); secIt ++) {
     int idx2 = *secIt;
     if (componentRef[idx1] != componentRef[idx2] && componentRef[idx1] >= 0 && componentRef[idx2] >= 0 && (seenNodes.count(idx1) == 0 || seenNodes[idx1].count(idx2) == 0) && (seenNodes.count(idx2) == 0 || seenNodes[idx2].count(idx1) == 0)) {
     cout << idx1 << "-" << componentRef[idx1] << ", " << idx2 << "-" << componentRef[idx2] << ": True positive contig/contig connection w/ < 4 matching peaks between " << idx1 << " (protein " << prot_match[idx1][0] << ") and " << idx2 << " (protein " << prot_match[idx2][0] << ")\n";
     }
     }
     }
     }
     */
    if (outputStats)
    {
      printf("\nFound %d components (%d with %d contigs or more)\n",
             compIdx,
             numCompsWSize,
             minCompSize);
      if (haveRes)
      {
        printf("Combined %s %.2f good edges (%.0f/%.0f)\n", "%", (100.0
            * goodEdges) / edgesAdded, goodEdges, edgesAdded);
      }
      printf("Maximum component length (Da) : %.1f\n", maxLen);
      printf("Average component length (Da) : %.1f\n", totalLen / numLen);
      if (haveRes)
      {
        printf("Maximum component length (AA) : %.0f\n", maxCover);
        printf("Average component length (AA) : %.1f\n", totalCover / numLen);
        printf("Found %s %.2f (%d/%d) of possible contig pairs with 6 matching peaks\n",
               "%",
               (100.0 * (float)good_aligned) / ((float)numCorrectPairs),
               good_aligned,
               numCorrectPairs);
      }
    }
  }

  void CombineContigs::recomputeEdges(int idx,
                                      map<int, map<int, int> >& _graph,
                                      vector<vector<float> >& _edges,
                                      vector<set<int> >& _components,
                                      set<int>& erasedEdges,
                                      vector<int>& edgeRef,
                                      vector<list<TwoValues<int> > >& _root_alignments,
                                      list<TwoValues<float> >& scoredEdges)
  {
    map<int, float> newEdgeScores;
    map<int, int>::iterator mapIt;
    list<TwoValues<float> >::iterator scoredEdgeIt;
    float newEdgeScore;
    set<int> seen;
    int oldSize = edgeRef.size();

    for (mapIt = _graph[idx].begin(); mapIt != _graph[idx].end(); mapIt++)
    {
      int node2 = mapIt->first;
      int edgeToNode = mapIt->second;
      seen.insert(node2);
      if (_components[idx].count(node2) > 0 || node2 == idx)
        continue;
      int edgeFromNode = _graph[node2][idx];
      float node2FShift = _edges[edgeToNode][2];
      float node2RShift = _edges[edgeToNode][3];
      TwoValues<float> res = getConsensusOverlap(idx,
                                                 node2,
                                                 node2FShift,
                                                 node2RShift,
                                                 _graph,
                                                 _edges,
                                                 _components,
                                                 _root_alignments,
                                                 false);
      newEdgeScore = max(res[0], res[1]);
      newEdgeScores[edgeToNode] = newEdgeScore;
      newEdgeScores[edgeFromNode] = newEdgeScore;
    }

    for (scoredEdgeIt = scoredEdges.begin(); scoredEdgeIt != scoredEdges.end(); scoredEdgeIt++)
    {
      int nextEdgeIdx = floatToInt((*scoredEdgeIt)[0]);
      if (nextEdgeIdx == -1)
        continue;
      while (edgeRef[nextEdgeIdx] != -1)
        nextEdgeIdx = edgeRef[nextEdgeIdx];
      if (erasedEdges.count(nextEdgeIdx) > 0)
        continue;
      if (newEdgeScores.count(nextEdgeIdx) > 0)
        (*scoredEdgeIt)[1] = newEdgeScores[nextEdgeIdx];
    }

    scoredEdges.sort(SortEdges());
  }

  TwoValues<float> CombineContigs::getOverlapScore(int idx1,
                                                   int idx2,
                                                   float fShift,
                                                   float rShift,
                                                   map<int, map<int, int> >& _graph,
                                                   vector<vector<float> >& _edges,
                                                   vector<set<int> >& _components)
  {
    float totalScoreF = 0, numScoreF = 0, totalScoreR = 0, numScoreR = 0;
    set<int>::iterator nodeIt, nodeIt2;
    TwoValues<float> res;
    TwoValues<int> mp;
    for (nodeIt = _components[idx1].begin(); nodeIt != _components[idx1].end(); nodeIt++)
    {
      int node1 = *nodeIt;
      float shift1 = 0;
      if (node1 != idx1)
        shift1 = _edges[_graph[idx1][node1]][2];

      for (nodeIt2 = _components[idx2].begin(); nodeIt2
          != _components[idx2].end(); nodeIt2++)
      {
        int node2 = *nodeIt2;
        if (node1 == node2)
        {
          continue;
        }
        float shift2 = 0, shiftr2 = 0;
        if (node2 != idx2)
        {
          shift2 = _edges[_graph[idx2][node2]][2];
          shiftr2 = _edges[_graph[idx2][node2]][5];
        }

        float locFShift = fShift - shift1 + shift2;
        float locRShift = rShift - shift1 + shiftr2;

        bool debug = false;//(idx1 == 35 and idx2 == 560);
        res = getContigOverlapScores(oriented[node1],
                                     oriented[node2],
                                     locFShift,
                                     locRShift,
                                     parentMassTol,
                                     minContigOverlap,
                                     debug);
        mp = getContigOverlapPeaks(oriented[node1],
                                   oriented[node2],
                                   locFShift,
                                   locRShift,
                                   parentMassTol,
                                   minContigOverlap,
                                   debug);
        //if (idx1 == 35 and idx      2 == 560)
        {
          cout << res[0] << " AND " << res[1] << ", " << mp[0] << " AND "
              << mp[1] << "\n";
        }
        /*
         if (idx1 == 51 and idx2 == 409) {
         float protShift1 = getContigShift(node1, prot_match[node1][2] == 1);
         float protShift2 = getContigShift(node2, prot_match[node2][2] == 1);
         float protFFShift = protShift2 - protShift1;
         float protRRShift = (protShift1+contigs[node1].parentMass) - (protShift2+contigs[node2].parentMass);
         cout << "overlapping " << node1 << " and " << node2 << ": " << res[0] << " " << res[1] << ":\n";
         cout << "shifts are " << locFShift << " and " << locRShift << ": one should match " << protFFShift << " or " << protRRShift << "\n\n";
         }
         */
        if (res[0] > -0.01)
        {
          numScoreF += 1.0;
          totalScoreF += res[0];
        }
        if (res[1] > -0.01)
        {
          numScoreR += 1.0;
          totalScoreR += res[1];
        }
      }
    }
    if (numScoreF < 0.01)
    {
      totalScoreF = 1.0;
      numScoreF = -1.0;
    }
    if (numScoreR < 0.01)
    {
      totalScoreR = 1.0;
      numScoreR = -1.0;
    }

    return TwoValues<float> (totalScoreF / numScoreF, totalScoreR / numScoreR);
  }

  TwoValues<float> CombineContigs::getConsensusShifts(int idx1,
                                                      int idx2,
                                                      float fShift,
                                                      float rShift,
                                                      map<int, map<int, int> >& _graph,
                                                      vector<vector<float> >& _edges,
                                                      vector<set<int> >& _components,
                                                      bool debug)
  {

    float totalScoreF, totalScoreR;
    float leftEdgeF = 0, rightEdgeF = 0, rightEdgeR = 0;
    for (set<int>::iterator nodeIt = _components[idx1].begin(); nodeIt
        != _components[idx1].end(); nodeIt++)
    {
      if (*nodeIt == idx1)
        continue;
      leftEdgeF = min(leftEdgeF, _edges[_graph[idx1][*nodeIt]][2]);
    }
    for (set<int>::iterator nodeIt = _components[idx2].begin(); nodeIt
        != _components[idx2].end(); nodeIt++)
    {
      if (*nodeIt == idx2)
        continue;
      rightEdgeF = min(rightEdgeF, _edges[_graph[idx2][*nodeIt]][2]);
      rightEdgeR = min(rightEdgeR, _edges[_graph[idx2][*nodeIt]][5]);
    }
    if (debug)
    {
      cout << "Root shift between " << idx1 << " and " << idx2 << " = "
          << fShift << ", " << rShift << "\n";
      cout << "consensusFShift = " << fShift << " - " << leftEdgeF << " + "
          << rightEdgeF << " = " << fShift - leftEdgeF + rightEdgeF << "\n";
      cout << "consensusRShift = " << rShift << " - " << leftEdgeF << " + "
          << rightEdgeR << " = " << rShift - leftEdgeF + rightEdgeR << "\n";
    }
    float consensusFShift = fShift - leftEdgeF + rightEdgeF;
    float consensusRShift = rShift - leftEdgeF + rightEdgeR;
    return TwoValues<float> (consensusFShift, consensusRShift);
  }

  TwoValues<float> CombineContigs::getConsensusOverlap(int idx1,
                                                       int idx2,
                                                       float fShift,
                                                       float rShift,
                                                       map<int, map<int, int> >& _graph,
                                                       vector<vector<float> >& _edges,
                                                       vector<set<int> >& _components,
                                                       vector<list<TwoValues<
                                                           int> > >& _root_alignments,
                                                       bool debug)
  {

    TwoValues<float> consensusShifts = getConsensusShifts(idx1,
                                                          idx2,
                                                          fShift,
                                                          rShift,
                                                          _graph,
                                                          _edges,
                                                          _components,
                                                          debug);

    Spectrum spec1 = consensus[idx1], spec2 = consensus[idx2];
    TwoValues<float> extraShift1 = extraShifts[idx1], extraShift2 =
        extraShifts[idx2];

    float consensusFShift = consensusShifts[0] - extraShift1[0]
        + extraShift2[0];
    float consensusRShift = consensusShifts[1] - extraShift1[0]
        + extraShift2[1];
    /*
     if (isEqual(spec1.parentMass - AAJumps::massHion - spec1.peakList.front()[0], AAJumps::massH2O, 0.01)) {
     for (int i = 0; i < spec1.size(); i++) {
     spec1[i][0] -= AAJumps::massH2O;
     }
     }
     */
    if (debug)
    {
      cout << "\nconsensusFShift = " << consensusFShift << " = "
          << consensusShifts[0] << " - " << extraShift1[0] << " + "
          << extraShift2[0] << "\n";
      cout << "consensusRShift = " << consensusRShift << " = "
          << consensusShifts[1] << " - " << extraShift1[0] << " + "
          << extraShift2[1] << "\n";
      cout << "Consensus " << idx1 << ":\n";
      spec1.output(cout);
      cout << "\nConsensus " << idx2 << ":\n";
      spec2.output(cout);
      cout << "\nConsensus " << idx2 << " (reversed):\n";
      Spectrum tempContig = spec2;
      tempContig.reverse(0.0 - AAJumps::massH2O, 0);
      tempContig.output(cout);
    }
    TwoValues<float> res = getContigOverlapScores(spec1,
                                                  spec2,
                                                  consensusFShift,
                                                  consensusRShift,
                                                  parentMassTol,
                                                  minContigOverlap,
                                                  debug);
    TwoValues<int> mp = getContigOverlapPeaks(spec1,
                                              spec2,
                                              consensusFShift,
                                              consensusRShift,
                                              parentMassTol,
                                              minContigOverlap,
                                              false);

    if (debug)
    {
      cout << "Forward score = " << res[0] << "(" << mp[0]
          << " matching peaks)\n";
      cout << "Reverse score = " << res[1] << "(" << mp[1]
          << " matching peaks)\n";
    }
    //if (res[0] > -0.01) res[0] = 0;
    //if (res[1] > -0.01) res[1] = 0;

    if ((idx1 == 67 && idx2 == 712) || (idx1 == 712 && idx2 == 67))
    {
      cout << "new score between " << idx1 << " and " << idx2 << ": " << res[0]
          * ((float)mp[0]) << ", " << res[1] * ((float)mp[1]) << endl;
    }

    return TwoValues<float> (res[0] * ((float)mp[0]), res[1] * ((float)mp[1]));
  }

  bool CombineContigs::tryMergeContigs(int idx1,
                                       int idx2,
                                       int edgeIdx,
                                       map<int, map<int, int> >& _graph,
                                       vector<vector<float> >& _edges,
                                       vector<set<int> >& _components,
                                       vector<list<TwoValues<int> > >& _root_alignments,
                                       bool rev,
                                       bool debug)
  {
    float conFFShift = (rev) ? _edges[edgeIdx][3] : _edges[edgeIdx][2];

    int most_matched_peaks = 0;
    TwoValues<int> best_node_pair;
    best_node_pair[0] = idx1;
    best_node_pair[1] = idx2;
    set<int>::iterator nodeIt, nodeIt2;
    for (nodeIt = _components[idx1].begin(); nodeIt != _components[idx1].end(); nodeIt++)
    {
      if (debug)
      {
        cout << "idx1 = " << idx1 << ", nodeIt = " << *nodeIt << "\n";
      }
      float shift1 = (*nodeIt == idx1) ? 0 : _edges[_graph[idx1][*nodeIt]][2];
      for (nodeIt2 = _components[idx2].begin(); nodeIt2
          != _components[idx2].end(); nodeIt2++)
      {
        if (*nodeIt2 == *nodeIt)
          continue;
        float shift_use = (rev) ? _edges[_graph[idx2][*nodeIt2]][5]
            : _edges[_graph[idx2][*nodeIt2]][2];
        float shift2 = (*nodeIt2 == idx2) ? 0 : shift_use;
        float locFShift = conFFShift - shift1 + shift2;
        Spectrum spec2 = oriented[*nodeIt2];
        if (rev)
          spec2.reverse(0.0 - AAJumps::massH2O);
        TwoValues<int> mp = getContigOverlapPeaks(oriented[*nodeIt],
                                                  spec2,
                                                  locFShift,
                                                  0,
                                                  parentMassTol,
                                                  minContigOverlap,
                                                  debug);
        if (mp[0] > most_matched_peaks)
        {
          best_node_pair[0] = *nodeIt;
          best_node_pair[1] = *nodeIt2;
          most_matched_peaks = mp[0];
        }
      }
    }
    if (most_matched_peaks == 0)
    {
      cout << "\nCOULD NOT FIND OVERLAPPING CONTIGS BETWEEN COMONENTS " << idx1
          << " AND " << idx2 << "\n";
    }
    if (debug)
    {
      cout << "\nbest node pair: " << best_node_pair[0] << " and "
          << best_node_pair[1] << "(" << most_matched_peaks
          << " matching peaks)\n";
    }
    list<TwoValues<int> > idx2_aligns(_root_alignments[idx1]);
    for (list<TwoValues<int> >::iterator alIt = _root_alignments[idx2].begin(); alIt
        != _root_alignments[idx2].end(); alIt++)
    {
      idx2_aligns.push_back(*alIt);
    }
    idx2_aligns.push_back(best_node_pair);

    if (debug)
    {
      cout << "_root_alignments[idx1].size() = "
          << _root_alignments[idx1].size() << "\n";
      cout << "_components[idx1].size() = " << _components[idx1].size() << "\n";
      cout << "_root_alignments[idx2].size() = "
          << _root_alignments[idx2].size() << "\n";
      cout << "_components[idx2].size() = " << _components[idx2].size() << "\n";
      cout << "Merged " << idx2_aligns.size() << " total alignments for "
          << _components[idx1].size() + _components[idx2].size()
          << " spectra\n";
    }
    //list<TwoValues<int> > idx2_aligns = root_alignments[idx2];
    //    idx2_aligns.push_back(best_node_pair);
    //    root_alignments[idx1].merge(idx2_aligns);
    //    components[idx1].insert(idx2);

    SpectrumPair next;
    SpectrumPairSet outAlign;

    SpecSet inspec;
    abinfo_t inabinfo;

    float minLeftEdgeF = 0;
    float maxRightEdge = (*oriented[idx1].back())[0];
    map<int, TwoValues<float> > contigIndent;

    SpecSet inputSpecs = oriented;

    //debug = idx == 0;

    if (debug)
      cout << "\nInput spectra to masab for " << idx1 << ":\n";
    for (nodeIt = _components[idx1].begin(); nodeIt != _components[idx1].end(); nodeIt++)
    {
      //node_to_index[*nodeIt] = outs.specs.size();
      //index_to_node[outs.specs.size()] = *nodeIt;
      //outs.specs.push_back(oriented[*nodeIt]);
      float FFShift = (*nodeIt == idx1) ? 0 : _edges[_graph[idx1][*nodeIt]][2];
      minLeftEdgeF = min(minLeftEdgeF, FFShift);
      float righEdge = FFShift + (*oriented[*nodeIt].back())[0];
      maxRightEdge = max(maxRightEdge, righEdge);
      if (debug)
      {
        cout << "\nContig " << *nodeIt << "(shift=" << FFShift << "):\n";
        oriented[*nodeIt].output(cout);
      }
    }
    for (nodeIt = _components[idx2].begin(); nodeIt != _components[idx2].end(); nodeIt++)
    {
      //node_to_index[*nodeIt] = outs.specs.size();
      //index_to_node[outs.specs.size()] = *nodeIt;
      //outs.specs.push_back(oriented[*nodeIt]);
      if (rev)
      {
        inputSpecs[*nodeIt].reverse(0.0 - AAJumps::massH2O);
      }
      float shift_use = (rev) ? _edges[_graph[idx2][*nodeIt]][5]
          : _edges[_graph[idx2][*nodeIt]][2];
      float FFShift = (*nodeIt == idx2) ? conFFShift : conFFShift + shift_use;
      minLeftEdgeF = min(minLeftEdgeF, FFShift);
      float righEdge = FFShift + (*oriented[*nodeIt].back())[0];
      maxRightEdge = max(maxRightEdge, righEdge);
      if (debug)
      {
        cout << "\nContig " << *nodeIt << "(shift=" << FFShift << "):\n";
        oriented[*nodeIt].output(cout);
      }
    }

    for (set<int>::iterator nodeIt = _components[idx1].begin(); nodeIt
        != _components[idx1].end(); nodeIt++)
    {
      //int locIdx = node_to_index[*nodeIt];
      float FFShift = (*nodeIt == idx1) ? 0 : _edges[_graph[idx1][*nodeIt]][2];
      contigIndent[*nodeIt][0] = FFShift - minLeftEdgeF;
    }

    for (set<int>::iterator nodeIt = _components[idx2].begin(); nodeIt
        != _components[idx2].end(); nodeIt++)
    {
      //int locIdx = node_to_index[*nodeIt];
      float shift_use = (rev) ? _edges[_graph[idx2][*nodeIt]][5]
          : _edges[_graph[idx2][*nodeIt]][2];
      float FFShift = (*nodeIt == idx2) ? conFFShift : conFFShift + shift_use;
      contigIndent[*nodeIt][0] = FFShift - minLeftEdgeF;
    }

    maxRightEdge -= minLeftEdgeF;
    if (debug)
      cout << "Found " << idx2_aligns.size() << " alignments\n";

    for (list<TwoValues<int> >::iterator node_pair_it =

    idx2_aligns.begin(); node_pair_it != idx2_aligns.end(); node_pair_it++)
    {
      int locIdx1 = (*node_pair_it)[0], locIdx2 = (*node_pair_it)[1];
      //int locIdx1 = node_to_index[idx1], locIdx2 = node_to_index[idx2];

      float FFShift;
      if (_components[idx1].count(locIdx1) > 0)
      {
        FFShift = (locIdx1 == idx1) ? 0 : _edges[_graph[idx1][locIdx1]][2];
      }
      else
      {
        float shift_use = (rev) ? _edges[_graph[idx2][locIdx1]][5]
            : _edges[_graph[idx2][locIdx1]][2];
        FFShift = (locIdx1 == idx2) ? conFFShift : conFFShift + shift_use;
      }
      float FFShift2;
      if (_components[idx1].count(locIdx2) > 0)
      {
        FFShift2 = (locIdx2 == idx1) ? 0 : _edges[_graph[idx1][locIdx2]][2];
      }
      else
      {
        float shift_use = (rev) ? _edges[_graph[idx2][locIdx2]][5]
            : _edges[_graph[idx2][locIdx2]][2];
        FFShift2 = (locIdx2 == idx2) ? conFFShift : conFFShift + shift_use;
      }

      float shift = FFShift2 - FFShift;

      next.spec1 = locIdx1;
      next.spec2 = locIdx2;
      next.score1 = 1.0;
      next.score2 = 1.0;
      next.shift1 = shift;
      next.shift2 = oriented[locIdx1].parentMass + 100.0;
      outAlign.push_back(next);
      if (debug)
      {
        cout << "Shift: " << locIdx1 << ", " << locIdx2 << " = " << shift
            << "\n";
      }
    }

    ParameterList inputParams;
    inputParams.addIfDoesntExist("PENALTY_PTM", "-200");
    inputParams.addIfDoesntExist("PENALTY_SAME_VERTEX", "-1000000");
    inputParams.addIfDoesntExist("GRAPH_TYPE", "2");
    inputParams.addIfDoesntExist("MAX_AA_JUMP", "2");
    inputParams.addIfDoesntExist("MAX_MOD_MASS", "100.0");
    inputParams.addIfDoesntExist("TOLERANCE_PEAK",
                                 parseFloat(merging_params->parent_mass_tol, 5));
    inputParams.addIfDoesntExist("TOLERANCE_PM",
                                 parseFloat(merging_params->parent_mass_tol, 5));
    inputParams.addIfDoesntExist("MIN_MATCHED_PEAKS",
                                 parseInt(merging_params->min_matched_peaks_contig_align).c_str());
    inputParams.addIfDoesntExist("MIN_EDGES_TO_COMPONENT", "0");
    inputParams.addIfDoesntExist("PATH_MIN_SPECS", "2");
    inputParams.addIfDoesntExist("PATH_MIN_PEAKS",
                                 parseInt(merging_params->min_matched_peaks_contig_align).c_str());
    inputParams.addIfDoesntExist("SPEC_TYPE_MSMS", "0");
    inputParams.addIfDoesntExist("NO_SEQUENCING", "0");
    inputParams.addIfDoesntExist("ADD_ENDPOINTS", "0");
    inputParams.addIfDoesntExist("OUTPUT_COMPLETE_ABRUIJN", "");
    inputParams.addIfDoesntExist("EDGE_SCORE_TYPE", "1");
    inputParams.addIfDoesntExist("IGNORE_REVERSALS", "1");
    inputParams.addIfDoesntExist("PARALLEL_PATHS", "0");

    //DEBUG_VAR(minNumMatchedPeaks);

    Clusters outputClusters;
    abinfo_t outputAbinfo;

    ExecAssembly assemblyObj(inputParams,
                             &inputSpecs,
                             &outAlign,
                             &outputClusters,
                             &outputAbinfo);

    //Logger& currentLogger = Logger::getDefaultLogger();
    //Logger::setDefaultLogger(Logger::getLogger(1));

    //Logger::setDefaultLogger(Logger::getLogger(2));

    assemblyObj.invoke();

    inspec = outputClusters.consensus;
    inabinfo = outputAbinfo;

    if (inspec.size() == 0 || inspec[0].size() == 0)
    {
      return false;
      /*cout << "\nMASAB RETURNED 0 COMPONENTS, MERGING FOR COMPONENT " << idx
       << " FAILED!\n";
       inspec[0].output(cout);
       putSpec = oriented[idx];
       return TwoValues<float> (0, 0);*/
    }
    else if (inspec.size() > 1)
    {
      cout << "\nCAUGHT MASAB SPLIT for (idx1=" << idx1 << ",idx2=" << idx2
          << ")!!!\n\n";
      return false;
    }

    return true;
  }

  TwoValues<float> CombineContigs::mergeContigs(int idx,
                                                map<int, map<int, int> >& _graph,
                                                vector<vector<float> >& _edges,
                                                vector<set<int> >& _components,
                                                vector<list<TwoValues<int> > >& _root_alignments,
                                                Spectrum & putSpec,
                                                pair<pair<vector<int> , vector<
                                                    int> > , vector<pair<
                                                    vector<int> ,
                                                    vector<double> > > >& putAb,
                                                bool debug)
  {

    SpectrumPair next;
    SpectrumPairSet outAlign;

    SpecSet inspec;
    abinfo_t inabinfo;

    float minLeftEdgeF = 0;
    float maxRightEdge = (*oriented[idx].back())[0];
    map<int, TwoValues<float> > contigIndent;

    SpecSet inputSpecs = oriented;

    //debug = idx == 0;

    if (debug)
      cout << "\nInput spectra to masab for " << idx << ":\n";
    for (set<int>::iterator nodeIt = _components[idx].begin(); nodeIt
        != _components[idx].end(); nodeIt++)
    {
      //node_to_index[*nodeIt] = outs.specs.size();
      //index_to_node[outs.specs.size()] = *nodeIt;
      //outs.specs.push_back(oriented[*nodeIt]);
      float FFShift = (*nodeIt == idx) ? 0 : _edges[_graph[idx][*nodeIt]][2];
      minLeftEdgeF = min(minLeftEdgeF, FFShift);
      float righEdge = FFShift + (*oriented[*nodeIt].back())[0];
      maxRightEdge = max(maxRightEdge, righEdge);
      if (debug)
      {
        cout << "\nContig " << *nodeIt << "(shift=" << FFShift << "):\n";
        oriented[*nodeIt].output(cout);
      }
    }

    for (set<int>::iterator nodeIt = _components[idx].begin(); nodeIt
        != _components[idx].end(); nodeIt++)
    {
      //int locIdx = node_to_index[*nodeIt];
      float FFShift = (*nodeIt == idx) ? 0 : _edges[_graph[idx][*nodeIt]][2];
      contigIndent[*nodeIt][0] = FFShift - minLeftEdgeF;
    }

    maxRightEdge -= minLeftEdgeF;
    if (debug)
      cout << "Found " << _root_alignments[idx].size() << " alignments\n";
    for (list<TwoValues<int> >::iterator node_pair_it =
        _root_alignments[idx].begin(); node_pair_it
        != _root_alignments[idx].end(); node_pair_it++)
    {
      int idx1 = (*node_pair_it)[0], idx2 = (*node_pair_it)[1];
      //int locIdx1 = node_to_index[idx1], locIdx2 = node_to_index[idx2];

      float FFShift = (idx1 == idx) ? 0 : _edges[_graph[idx][idx1]][2];
      float FFShift2 = (idx2 == idx) ? 0 : _edges[_graph[idx][idx2]][2];

      float shift = FFShift2 - FFShift;

      next.spec1 = idx1;
      next.spec2 = idx2;
      next.score1 = 1.0;
      next.score2 = 1.0;
      next.shift1 = shift;
      next.shift2 = oriented[idx1].parentMass + 100.0;
      outAlign.push_back(next);
      if (debug)
      {
        cout << "Shift: " << idx1 << ", " << idx2 << " = " << shift << "\n";
      }
    }

    ParameterList inputParams;
    inputParams.addIfDoesntExist("PENALTY_PTM", "-200");
    inputParams.addIfDoesntExist("PENALTY_SAME_VERTEX", "-1000000");
    inputParams.addIfDoesntExist("GRAPH_TYPE", "2");
    inputParams.addIfDoesntExist("MAX_AA_JUMP", "2");
    inputParams.addIfDoesntExist("MAX_MOD_MASS", "100.0");
    inputParams.addIfDoesntExist("TOLERANCE_PEAK",
                                 parseFloat(merging_params->parent_mass_tol, 5));
    inputParams.addIfDoesntExist("TOLERANCE_PM",
                                 parseFloat(merging_params->parent_mass_tol, 5));
    inputParams.addIfDoesntExist("MIN_MATCHED_PEAKS",
                                 parseInt(merging_params->min_matched_peaks_contig_align).c_str());
    inputParams.addIfDoesntExist("MIN_EDGES_TO_COMPONENT", "0");
    inputParams.addIfDoesntExist("PATH_MIN_SPECS", "2");
    inputParams.addIfDoesntExist("PATH_MIN_PEAKS",
                                 parseInt(merging_params->min_matched_peaks_contig_align).c_str());
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
                             &outAlign,
                             &outputClusters,
                             &outputAbinfo);

    //Logger& currentLogger = Logger::getDefaultLogger();
    //Logger::setDefaultLogger(Logger::getLogger(1));

    //Logger::setDefaultLogger(Logger::getLogger(2));

    assemblyObj.invoke();

    inspec = outputClusters.consensus;
    inabinfo = outputAbinfo;

    if (inspec.size() == 0 || inspec[0].size() == 0)
    {
      cout << "\nMASAB RETURNED 0 COMPONENTS, MERGING FOR COMPONENT " << idx
          << " FAILED!\n";
      inspec[0].output(cout);
      putSpec = oriented[idx];
      return TwoValues<float> (0, 0);
    }
    else if (inspec.size() > 1)
    {
      cout << "\nMASAB SPLIT COMPONENT " << idx << " INTO " << inspec.size()
          << " CONTIGS!\n";
    }
    putSpec = inspec[0];
    putAb = inabinfo[0];

    putSpec.parentMass = (*putSpec.back())[0] + AAJumps::massMH;

    /*
     typedef std::map<
     unsigned,  // contig index
     std::pair<
     std::pair< vector<int>, vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
     vector< std::pair< // ABruijn vertices
     vector<int>, vector<double> > // Spectrum index, peak mass
     >
     >
     > abinfo_t;
     */
    pair<vector<int> , vector<double> > first_mass = inabinfo[0].second.front();
    int locIdxF = first_mass.first.front();
    float gapF = first_mass.second.front() + contigIndent[locIdxF][0];

    //Spectrum revPut = putSpec;
    //revPut.reverse(0.0 - AAJumps::massH2O, 0)
    float gapR = maxRightEdge - (*putSpec.back())[0] - gapF;// - ((putSpec.parentMass - AAJumps::massHion - putSpec.peakList.back()[0]) - putSpec.peakList.front()[0]);

    /*
     pair< vector<int>, vector<double> > last_mass = inabinfo[0].second.back();
     int locIdxR = last_mass.first.front();
     float gapR = putSpec.parentMass - putSpec.peakList.back()[0] - AAJumps::massHion;// - last_mass.second.front();
     */
    if (debug)
    {
      cout << "\ngapF = " << gapF << " = " << first_mass.second.front()
          << " + " << contigIndent[locIdxF][0] << "\n";
      cout << "gapR = " << gapR << " = " << maxRightEdge << " - "
          << (*putSpec.back())[0] << " - " << gapF << "\n";
      cout << "First mass from " << locIdxF << "\n";
      cout << "minLeftEdgeF " << minLeftEdgeF << "\n";
    }

    return TwoValues<float> (gapF, gapR);
  }

  void CombineContigs::reverseNode(int idx,
                                   map<int, map<int, int> >& _graph,
                                   vector<vector<float> >& _edges,
                                   vector<set<int> >& _components,
                                   vector<bool>& _reversed)
  {
    map<int, int>::iterator mapIt;
    float FFShift, FRShift, RFShift, RRShift;

    consensus[idx].reverse(0.0 - AAJumps::massH2O, 0);
    float tempExtr1 = extraShifts[idx][0], tempExtr2 = extraShifts[idx][1];
    extraShifts[idx][0] = tempExtr2;
    extraShifts[idx][1] = tempExtr1;
    oriented[idx].reverse(0.0 - AAJumps::massH2O, 0);
    _reversed[idx] = !_reversed[idx];

    for (mapIt = _graph[idx].begin(); mapIt != _graph[idx].end(); mapIt++)
    {
      int idx2 = mapIt->first;
      int edge2 = mapIt->second;

      if (idx2 == idx)
        continue;

      FFShift = _edges[edge2][2];
      FRShift = _edges[edge2][3];
      RFShift = _edges[edge2][4];
      RRShift = _edges[edge2][5];
      if (_components[idx].count(idx2) > 0)
      {
        _edges[edge2][2] = RRShift;
        _edges[edge2][3] = RFShift;
        _edges[edge2][4] = FRShift;
        _edges[edge2][5] = FFShift;
        oriented[idx2].reverse(0.0 - AAJumps::massH2O, 0);
        consensus[idx2].reverse(0.0 - AAJumps::massH2O, 0);
        _reversed[idx2] = !_reversed[idx2];
      }
      else
      {
        _edges[edge2][2] = RFShift;
        _edges[edge2][3] = RRShift;
        _edges[edge2][4] = FFShift;
        _edges[edge2][5] = FRShift;
        int edgeFrom = _graph[idx2][idx];
        _edges[edgeFrom][2] = 0.0 - RFShift;
        _edges[edgeFrom][3] = 0.0 - FFShift;
        _edges[edgeFrom][4] = 0.0 - RRShift;
        _edges[edgeFrom][5] = 0.0 - FRShift;
      }
    }
  }
  /*
   void CombineContigs::mergeContigs(Spectrum& spec1, Spectrum& spec2, float consensusShift, Spectrum& putSpec) {
   Results_PA next;
   vector<Results_PA> outAlign;
   SpecSet outs;
   SpecSet inspec;
   outs.specs.push_back(spec1);
   outs.specs.push_back(spec2);
   next.spec1 = 0;
   next.spec2 = 1;
   next.score1 = 1.0;
   next.score2 = 1.0;
   next.shift1 = consensusShift;
   next.shift2 = spec1.parentMass + 100.0;
   outAlign.push_back(next);
   chdir("/home/aguthals/tem p_masab");
   Save_results_bin("sps_merge_aligns.bin", (unsigned int)outAlign.size(), outAlign.begin());
   outs.SaveSpecSet_pklbin("sps_orient.pklbin");
   system("/home/aguthals/cpplib/masab assembly.params");
   //execl( ".", "/home/aguthals/cpplib/masab assembly.params");
   inspec.LoadSpecSet_pklbin("sps_seqs_merged.pklbin");
   chdir("/home/aguthals/merge_cid_hcd_karl_4mp");
   putSpec = inspec[0];
   }
   */
  void CombineContigs::getCondensedContigAlignmentGraph(map<int, map<int, int> >& _graph,
                                                        vector<vector<float> >& _edges,
                                                        list<TwoValues<float> >& scoredEdges,
                                                        float minCombScore,
                                                        int minMatchedPeaks)
  {
    list<SpectrumPair>::iterator pivot_res;
    _graph.clear();
    scoredEdges.clear();
    int edgeIdxUse = _edges.size();
    _edges.resize((contigs.size() * contigs.size()) - contigs.size());
    edge.resize(6);
    TwoValues<int> mp;
    TwoValues<float> res;

    int total = 0;
    int used = 0;

    float FFShift, FRShift, RFShift, RRShift, shift, mRatio, score,
        reversedShift;
    for (pivot_res = contig_alignments.begin(); pivot_res
        != contig_alignments.end(); pivot_res++)
    {
      int idx1 = (*pivot_res).spec1, idx2 = (*pivot_res).spec2;

      if (haveRes && (idxNotUsed.count(prot_match[idx1][0]) > 0
          || idxNotUsed.count(prot_match[idx2][0]) > 0))
        continue;
      shift = (*pivot_res).shift1, mRatio = min((*pivot_res).score1,
                                                (*pivot_res).score2);
      bool reverse2 = isEqual((*pivot_res).shift2, 1.0, 0.0001);

      total++;
      if (mRatio < minCombScore)
        continue;
      //bool correct = validShift(shift, idx1, idx2, reverse2);
      //if (! correct || prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0) continue;

      score = mRatio;
      reversedShift = contigs[idx1].parentMass - (shift
          + contigs[idx2].parentMass);

      used++;

      FFShift = (reverse2) ? contigs[idx1].parentMass : shift;
      FRShift = (reverse2) ? shift : contigs[idx1].parentMass;
      RFShift = (reverse2) ? reversedShift : contigs[idx1].parentMass;
      RRShift = (reverse2) ? contigs[idx1].parentMass : reversedShift;

      edge[0] = (float)idx1;
      edge[1] = (float)idx2;
      edge[2] = FFShift;
      edge[3] = FRShift;
      edge[4] = RFShift;
      edge[5] = RRShift;
      scoredEdge[0] = (float)edgeIdxUse;
      scoredEdge[1] = score;
      if (_graph.count(idx1) == 0)
      {
        edgeMap.clear();
        edgeMap[idx2] = edgeIdxUse;
        _graph[idx1] = edgeMap;
      }
      else
        _graph[idx1][idx2] = edgeIdxUse;

      _edges[edgeIdxUse] = edge;
      scoredEdges.push_back(scoredEdge);

      if ((idx1 == 67 && idx2 == 712) || (idx1 == 712 && idx2 == 67))
      {
        cout << "adding edge " << idx1 << " to " << idx2 << ", edgId x == "
            << edgeIdxUse << ", score = " << scoredEdge[1] << endl;
      }

      edgeIdxUse++;

      FFShift = (reverse2) ? contigs[idx2].parentMass : 0.0 - shift;
      FRShift = (reverse2) ? 0.0 - reversedShift : contigs[idx2].parentMass;
      RFShift = (reverse2) ? 0.0 - shift : contigs[idx2].parentMass;
      RRShift = (reverse2) ? contigs[idx2].parentMass : 0.0 - reversedShift;

      edge[0] = (float)idx2;
      edge[1] = (float)idx1;
      edge[2] = FFShift;
      edge[3] = FRShift;
      edge[4] = RFShift;
      edge[5] = RRShift;
      if (_graph.count(idx2) == 0)
      {
        edgeMap.clear();
        edgeMap[idx1] = edgeIdxUse;
        _graph[idx2] = edgeMap;
      }
      else
        _graph[idx2][idx1] = edgeIdxUse;

      _edges[edgeIdxUse] = edge;
      edgeIdxUse++;

    }
    _edges.resize(edgeIdxUse);
    scoredEdges.sort(SortEdges());
    cout << "\nUsed " << used << " of " << total << " contig-contig shifts\n";
  }
  /*
   graph: idx1 --> idx2 --> edgeIdx (undirected)
   edges: edgeIdx --> (idx1, idx2, FFShift, FRShift, RFShift, RRShift)
   scoredEdges: list of (edgeIdx, score) sorted in by decreasing score
   */

  bool CombineContigs::validShift(float FFShift,
                                  int idx1,
                                  int idx2,
                                  bool isReversed)
  {
    bool reverse1 = prot_match[idx1][2] == 1;
    bool reverse2 = prot_match[idx2][2] == 1;
    float protShift1 = getContigShift(idx1, reverse1);
    float protShift2 = getContigShift(idx2, reverse2);
    float protFFShift = protShift2 - protShift1;
    float protRRShift = (protShift1 + contigs[idx1].parentMass) - (protShift2
        + contigs[idx2].parentMass);
    protFFShift = (reverse2 == isReversed) ? protFFShift : protRRShift;
    bool correct = isEqual(FFShift, protFFShift, AAJumps::massHion
        + merging_params->parent_mass_tol) && prot_match[idx1][0]
        == prot_match[idx2][0];
    return correct;
  }

  bool CombineContigs::validShiftMod(float FFShift,
                                     int idx1,
                                     int idx2,
                                     bool isReversed)
  {
    bool reverse1 = prot_match[idx1][2] == 1;
    bool reverse2 = prot_match[idx2][2] == 1;
    float protShift1 = getContigShift(idx1, reverse1);
    float protShift2 = getContigShift(idx2, reverse2);
    float protFFShift = protShift2 - protShift1;
    float protRRShift = (protShift1 + contigs[idx1].parentMass) - (protShift2
        + contigs[idx2].parentMass);
    protFFShift = (reverse2 == isReversed) ? protFFShift : protRRShift;
    bool correct = isEqual(FFShift, protFFShift, 58.0) && prot_match[idx1][0]
        == prot_match[idx2][0];
    return correct;
  }

  float CombineContigs::getContigShift(int index, bool reverse)
  {
    if (!haveRes)
      return 0;
    unsigned int pepIndex = prot_match[index][0];
    char* peptide = fasta.getSequence(pepIndex);
    Spectrum masses = fasta.getMassesSpec(pepIndex);
    float pmTol = parentMassTol;
    Spectrum overlap = overlaps[index];
    Spectrum contig_spec = contigs[index];
    if (reverse)
      contig_spec.reverse(0.0 - AAJumps::massH2O, 0);

    int end = (int)(overlap[overlap.size() - 1][1] + 0.01);
    list<int>::iterator start_other_pivot;
    map<int, float> shiftScore;
    map<int, float>::iterator shiftScoreIt;
    float best_shift1;
    float shift;
    float intensity;
    char sequence[strlen(peptide)];
    strcpy(sequence, peptide);
    int intShift;
    int intPmTol = intParentMassTol;

    for (int i = 0; i < overlap.size(); i++)
    {
      int cIdx = floatToInt(overlap[i][0]);
      int pIdx = floatToInt(overlap[i][1]);
      shift = masses[pIdx][0] - contig_spec[cIdx][0];
      //if (index == 98) cout << cIdx << " : " << contig_spec[cIdx][0] << ", " << pIdx << " : " << masses[pIdx][0] << " = " << shift << "\n";
      intShift = floatToInt(shift / merging_params->resolution);
      intensity = contig_spec[cIdx][1];

      for (int idx = intShift - intPmTol; idx <= intShift + intPmTol; idx++)
      {
        if (shiftScore.count(idx) == 0)
          shiftScore[idx] = intensity - (0.0001 * (float)abs(intShift - idx));
        else
          shiftScore[idx] += intensity - (0.0001 * (float)abs(intShift - idx));
      }
    }
    float maxScore = 0.0;
    for (shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt++)
    {
      if (shiftScoreIt->second > maxScore)
      {
        maxScore = shiftScoreIt->second;
        best_shift1 = ((float)shiftScoreIt->first) * merging_params->resolution;
      }
    }
    return best_shift1;
  }

  void CombineContigs::getContigDistances(map<int, map<int, float> >& contigContigShifts)
  {
    if (!haveRes)
      return;
    contigContigShifts.clear();
    map<int, float> contigShifts;
    map<int, float> computedShifts;
    for (int i = 0; i < contigs.size(); i++)
    {
      int pepIndex = prot_match[i][0];

      float shift1;
      if (computedShifts.count(i) == 0)
      {
        shift1 = getContigShift(i, prot_match[i][2] == 1);
        computedShifts[i] = shift1;
      }
      else
        shift1 = computedShifts[i];

      contigShifts.clear();
      contigContigShifts[i] = contigShifts;

      for (int j = 0; j < contigs.size(); j++)
      {
        if (i == j || pepIndex != prot_match[j][0])
          continue;
        float shift2;
        if (computedShifts.count(j) == 0)
        {
          shift2 = getContigShift(j, prot_match[j][2] == 1);
          computedShifts[j] = shift2;
        }
        else
          shift2 = computedShifts[j];

        contigContigShifts[i][j] = shift2 - shift1;
      }
    }
  }

}
