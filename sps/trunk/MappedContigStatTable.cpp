/*
 * MappedContigStatTable.cpp
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

#include "MappedContigStatTable.h"

namespace specnets
{

  MappedContigStatTable::MappedContigStatTable() :
    OutputTable(), mapped_sps_proj(0x0)
  {

  }

  MappedContigStatTable::MappedContigStatTable(MappedSpecnets* _mapped_sps_proj) :
    OutputTable(), mapped_sps_proj(_mapped_sps_proj)
  {
  }

  /**
   * Prepares output table with all necessary statistics
   * @return
   */
  void MappedContigStatTable::prepareTable(int index)
  {
    //cout << "preparing " << index << "\n"; cout.flush();
    MappedContig* contig = &(*(mapped_sps_proj->contigs))[index];
    if (contig->length == 0)
    {
      values.resize(0);
      return;
    }
    //cout << "preparing 1 " << index << "\n"; cout.flush();
    pair<string, bool> dflt("", false);
    int numRows = ((contig->length - 1) * 12) + 10 + 35;
    vector<pair<string, bool> > dflt_vec(2, dflt);
    values.resize(numRows, dflt_vec);
    //cout << "preparing 2 " << index << "\n"; cout.flush();
    setValue(0, 0, "Contig Index", true);
    setValue(0, 1, parseInt(index));
    values[1].resize(0);

    setValue(2, 0, "# of Vertices", true);
    setValue(2, 1, parseInt(contig->length));

    setValue(3, 0, "# of Assembled Spectra", true);
    setValue(3, 1, parseInt(contig->numSpecs));

    setValue(4, 0, "# of Assembled Peptides", true);
    setValue(4, 1, parseInt(contig->numPeptides));

    setValue(5, 0, "# of Assembled Annotated Spectra", true);
    setValue(5, 1, parseInt(contig->getNumAnnotSpec()));

    setValue(6, 0, "# of Assembled Mapped Spectra", true);
    setValue(6, 1, parseInt(contig->getNumMappedSpec()));
    values[7].resize(0);

    setValue(8, 0, "% Annotated Vertices from b or y", true);
    setValue(8, 1, parseFloat(contig->getPercBYVerts(), 1));

    setValue(9, 0, "% Annotated Vertices from b", true);
    setValue(9, 1, parseFloat(contig->getPercBVerts(), 1));

    setValue(10, 0, "% Annotated Vertices from y", true);
    setValue(10, 1, parseFloat(contig->getPercYVerts(), 1));
    values[11].resize(0);

    setValue(12, 0, "Chimeric", true);
    string val = (contig->chimeric) ? "1" : "0";
    setValue(12, 1, val);
    values[13].resize(0);

    setValue(14, 0, "% Annotated Vertices Correct", true);
    setValue(14, 1, parseFloat(contig->getPercVerts(3), 1));

    setValue(15, 0, "% Annotated Vertices Chimeric", true);
    setValue(15, 1, parseFloat(contig->getPercVerts(2), 1));

    setValue(16, 0, "% Annotated Vertices Incorrect", true);
    setValue(16, 1, parseFloat(contig->getPercVerts(1), 1));

    setValue(17, 0, "% Un-annotated Vertices", true);
    setValue(17, 1, parseFloat(contig->getPercVerts(0), 1));
    values[18].resize(0);

    pair<float, float> res = contig->getPercCallsAcc(2);
    setValue(19, 0, "# Annotated Gaps Correct", true);
    setValue(19, 1, parseFloat(res.first, 0));

    res = contig->getPercCallsAcc(1);
    setValue(20, 0, "# Annotated Gaps Incorrect", true);
    setValue(20, 1, parseFloat(res.first, 0));

    res = contig->getPercCallsAcc(0);
    setValue(21, 0, "# Un-annotated Gaps", true);
    setValue(21, 1, parseFloat(res.first, 0));
    values[22].resize(0);

    res = contig->getPercCallsLen(0);
    setValue(23, 0, "# Single Jump Calls", true);
    setValue(23, 1, parseFloat(res.first, 0));

    res = contig->getPercCallsAcc(1);
    setValue(24, 0, "# Consec Gap Calls", true);
    setValue(24, 1, parseFloat(res.first, 0));

    res = contig->getPercCallsAcc(2);
    setValue(25, 0, "# Non-consec Gap Calls", true);
    setValue(25, 1, parseFloat(res.first, 0));
    values[26].resize(0);

    setValue(27, 0, "Inspect Protein Index", true);
    setValue(27, 1, parseInt(contig->starProtIdx));

    setValue(28, 0, "Matchma Protein Index", true);
    setValue(28, 1, parseInt(contig->vertProtIdx));

    setValue(29, 0, "Matchma Reversed", true);
    val = (contig->reversed) ? "1" : "0";
    setValue(29, 1, val);

    setValue(30, 0, "Matchma AA Coverage", true);
    setValue(30, 1, parseInt(contig->lastResidue - contig->firstResidue + 1));

    setValue(31, 0, "Matchma Beg Residue Index", true);
    setValue(31, 1, parseInt(contig->firstResidue));

    setValue(32, 0, "Matchma End Residue Index", true);
    setValue(32, 1, parseInt(contig->lastResidue));
    values[33].resize(0);

    values[34].resize(7, dflt);
    setValue(34, 0, "Abruijn Vertex/Gap", true);
    setValue(34, 1, "Mass", true);
    setValue(34, 2, "Residue Idx", true);
    setValue(34, 3, "Intensity", true);
    setValue(34, 4, "Label", true);
    setValue(34, 5, "Annotation", true);
    setValue(34, 6, "Assembled Ions", true);

    string sep("-----");
    string list_sep = ";";

    for (int i = 0; i < contig->length; i++)
    {
      int tr = (i * 12) + 35;

      for (int j = 0; j < 10; j++)
      {
        values[tr + j].resize(contig->abruijnVerts[i].starPeaks.size() + 7,
                              dflt);
      }

      setValue(tr, 0, parseInt(i));
      setValue(tr, 1, parseFloat(contig->abruijnVerts[i].getMass(), 1));

      val = (contig->abruijnVerts[i].residueIdx >= 0)
          ? parseInt(contig->abruijnVerts[i].residueIdx) : "";
      setValue(tr, 2, val);
      setValue(tr, 3, parseFloat(contig->abruijnVerts[i].getIntensity(), 1));
      setValue(tr, 4, MappedVertex::Labels[contig->abruijnVerts[i].label]);
      setValue(tr, 5, contig->abruijnVerts[i].annotation);

      setValue(tr, 6, "Spectrum Idx", true);
      setValue(tr + 1, 6, "Spectrum Rev", true);
      setValue(tr + 2, 6, "Spectrum Ident", true);
      setValue(tr + 3, 6, "Peak Idx", true);
      setValue(tr + 4, 6, "Peak Mass", true);
      setValue(tr + 5, 6, "Peak Endpt", true);
      setValue(tr + 6, 6, "Protein Idx", true);
      setValue(tr + 7, 6, "Residue Idx", true);
      setValue(tr + 8, 6, "Annotation", true);

      for (int j = 0; j < contig->abruijnVerts[i].starPeaks.size(); j++)
      {
        int tc = j + 7;
        MappedPeak* pk = contig->abruijnVerts[i].starPeaks[j];
        MappedSpectrum* spec = &((*mapped_sps_proj->spectra)[pk->specIdx]);

        setValue(tr, tc, parseInt(pk->specIdx));
        val = (spec->reversed) ? "1" : "0";
        setValue(tr + 1, tc, val);
        val = (spec->identified) ? "1" : "0";
        setValue(tr + 2, tc, val);
        setValue(tr + 3, tc, parseInt(pk->peakIdx));
        setValue(tr + 4, tc, parseFloat(pk->getMass(), 1));
        val = (pk->endpt) ? "1" : "0";
        setValue(tr + 5, tc, val);

        string prot_str = "";
        string res_str = "";
        if (spec->mapped)
        {

          bool startL = true;
          for (map<int, list<int> >::iterator protIt =
              spec->residueIdxs.begin(); protIt != spec->residueIdxs.end(); protIt++)
          {
            for (list<int>::iterator resIt = protIt->second.begin(); resIt
                != protIt->second.end(); resIt++)
            {
              if (!startL)
              {
                prot_str.append(list_sep);
                if (pk->mapped)
                {
                  res_str.append(list_sep);
                }
              }
              startL = false;
              prot_str.append(parseInt(protIt->first));
              if (pk->mapped)
              {
                if (pk->BpeptideIdx >= 0)
                {
                  res_str.append(parseInt((*resIt) + pk->BpeptideIdx));
                }
                if (pk->YpeptideIdx >= 0)
                {
                  if (pk->BpeptideIdx >= 0)
                  {
                    res_str.append(list_sep);
                  }
                  res_str.append(parseInt((*resIt) + pk->YpeptideIdx));
                }
              }
            }
          }
        }
        setValue(tr + 6, tc, prot_str);
        setValue(tr + 7, tc, res_str);

        setValue(tr + 8, tc, pk->annotation);
        setValue(tr + 9, tc, sep, true);
      }
      for (int j = 0; j < 7; j++)
      {
        setValue(tr + 9, j, sep, true);
      }

      if (i < contig->length - 1)
      {
        for (int j = 10; j < 12; j++)
        {
          values[tr + j].resize(5, dflt);
        }
        string gap = "";
        gap.append(parseInt(i));
        gap.append("-");
        gap.append(parseInt(i + 1));

        setValue(tr + 10, 0, gap);
        setValue(tr + 10, 1, parseFloat(contig->abruijnGaps[i].getMass(), 1));

        setValue(tr + 10, 3, parseFloat(contig->abruijnGaps[i].getIntensity(),
                                        1));

        setValue(tr + 10, 4, MappedGap::Labels[contig->abruijnGaps[i].label]);

        for (int j = 0; j < 5; j++)
        {
          setValue(tr + 11, j, sep, true);
        }
      }
    }

  }
}
