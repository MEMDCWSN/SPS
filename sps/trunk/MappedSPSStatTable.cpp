/*
 * MappedSpecnetsStatTable.cpp
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

#include "MappedSPSStatTable.h"

namespace specnets
{

  MappedSPSStatTable::MappedSPSStatTable() :
    OutputTable(), mapped_sps_proj(0x0)
  {
  }

  MappedSPSStatTable::MappedSPSStatTable(MappedSpecnets* _mapped_sps_proj) :
    OutputTable(), mapped_sps_proj(_mapped_sps_proj)
  {
  }

  /**
   * Prepares output table with all necessary statistics
   * @return
   */
  void MappedSPSStatTable::prepareTable()
  {
    pair<string, bool> deflt("", false);
    int numRows = 36;

    values.resize(numRows);
    int row = 0;

    values[row].resize(2, deflt);
    values[row][0].first = "# Contigs";
    values[row][0].second = true;
    values[row++][1].first = parseInt(mapped_sps_proj->getNumContigs());

    values[row].resize(2, deflt);
    values[row][0].first = "# Star Spectra";
    values[row][0].second = true;
    values[row++][1].first = parseInt(mapped_sps_proj->getNumSpectra());

    values[row].resize(2, deflt);
    values[row][0].first = "# ID Star Spectra";
    values[row][0].second = true;
    values[row++][1].first = parseInt(mapped_sps_proj->getNumSpecIdent());

    values[row].resize(2, deflt);
    values[row][0].first = "# Modified ID Star Spectra";
    values[row][0].second = true;
    values[row++][1].first = parseInt(mapped_sps_proj->getNumSpecModified());

    values[row++].resize(0);

    for (int i = 4; i < numRows; i++)
    {
      values[i].resize(mapped_sps_proj->proteins->size() + 2, deflt);
      //cout << "resizing " << i << " to " << mapped_sps_proj->proteins->size() + 2 << ", actual size is " << values[i].size() << "\n"; cout.flush();
    }

    values[row][0].first = "Protein Index";
    values[row][0].second = true;
    values[row][1].first = "all";
    values[row++][1].second = true;

    values[row][0].first = "Identified or Mapped Contigs";
    values[row++][0].second = true;

    values[row][0].first = "Identified Contigs";
    values[row++][0].second = true;

    values[row][0].first = "Mapped Contigs";
    values[row++][0].second = true;

    values[row][0].first = "Identified and Mapped Contigs";
    values[row++][0].second = true;

    values[row][0].first = "Identified Spectra";
    values[row++][0].second = true;

    values[row][0].first = "Assembled Spectra";
    values[row++][0].second = true;

    values[row][0].first = "Assembled Identified Spectra (%)";
    values[row++][0].second = true;

    values[row][0].first = "Assembled Identified Modified Spectra (%)";
    values[row++][0].second = true;

    values[row][0].first = "Correctly Mapped Verts (%)";
    values[row++][0].second = true;

    values[row][0].first = "Spectrum Coverage (%)";
    values[row++][0].second = true;

    values[row][0].first = "Sequencing Coverage (%)";
    values[row++][0].second = true;

    values[row][0].first = "Sequencing Coverage Redundancy";
    values[row++][0].second = true;

    values[row][0].first = "Spectra Assembled Per Contig";
    values[row++][0].second = true;

    values[row][0].first = "Peptides Assembled Per Contig";
    values[row++][0].second = true;

    values[row][0].first = "Average Length (Da)";
    values[row++][0].second = true;

    values[row][0].first = "Average Length (AA)";
    values[row++][0].second = true;

    values[row][0].first = "Longest Sequence (AA)";
    values[row++][0].second = true;

    values[row][0].first = "Longest Contig";
    values[row++][0].second = true;

    values[row][0].first = "Annotated Vertices Correct (%)";
    values[row++][0].second = true;

    values[row][0].first = "Annotated Vertices Chimeric (%)";
    values[row++][0].second = true;

    values[row][0].first = "Annotated Vertices Incorrect (%)";
    values[row++][0].second = true;

    values[row][0].first = "Un-annotated Vertices (%)";
    values[row++][0].second = true;

    values[row][0].first = "B Vertices (%)";
    values[row++][0].second = true;

    values[row][0].first = "Y Vertices (%)";
    values[row++][0].second = true;

    values[row][0].first = "Calls Correct (%)";
    values[row++][0].second = true;

    values[row][0].first = "Calls Incorrect (%)";
    values[row++][0].second = true;

    values[row][0].first = "Un-annotated Calls (%)";
    values[row++][0].second = true;

    values[row][0].first = "Single Jump Calls (%)";
    values[row++][0].second = true;

    values[row][0].first = "Consec Gap Calls (%)";
    values[row++][0].second = true;

    values[row][0].first = "Non-consec Gap Calls (%)";
    values[row++][0].second = true;

    int numProts = mapped_sps_proj->proteins->size();
    //cout << "proteins size=" << numProts << "\n"; cout.flush();

    for (int prot_idx = -1; prot_idx < numProts; prot_idx++)
    {

      int tp = prot_idx + 2;
      //cout << "on " << prot_idx << "\n";
      if (prot_idx >= 0 && (*mapped_sps_proj->proteins)[prot_idx].length() == 0)
      {
        //        cout << "skipping " << prot_idx << "\n";
        continue;
      }

      row = 5;

      if (prot_idx >= 0)
      {
        values[row++][tp].first = parseInt(prot_idx);
        //        cout << "setting " << 4 << ", " << tp << " to " << parseInt(prot_idx) << ", actually got " << values[4][tp].first << "\n";
      }
      else
      {
        row++;
      }

      values[row++][tp].first
          = parseInt(mapped_sps_proj->getNumContigsMapped(prot_idx));
      values[row++][tp].first
          = parseInt(mapped_sps_proj->getNumContigsSpecMapped(prot_idx));
      values[row++][tp].first
          = parseInt(mapped_sps_proj->getNumContigsVertMapped(prot_idx));
      values[row++][tp].first
          = parseInt(mapped_sps_proj->getNumContigsVertSpecMapped(prot_idx));
      values[row++][tp].first
          = parseInt(mapped_sps_proj->getNumSpecMapped(prot_idx));
      values[row++][tp].first
          = parseInt(mapped_sps_proj->getNumAssembledSpec(prot_idx));
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercAssemSpecIdent(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercAssemIdentSpecMod(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getMatchmaAccuracy(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercSpecCov(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercSeqCov(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getCovRedundancy(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getSpecPerContig(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPepPerContig(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getDaLengthPerContig(prot_idx), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getAALengthPerContig(prot_idx), 1);
      pair<int, int> aaRes = mapped_sps_proj->getLongestAAContig(prot_idx);
      values[row++][tp].first = parseInt(aaRes.first);
      values[row++][tp].first = parseInt(aaRes.second);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercVerts(prot_idx, 3), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercVerts(prot_idx, 2), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercVerts(prot_idx, 1), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercVerts(prot_idx, 0), 1);
      pair<float, float> byRes = mapped_sps_proj->getPercBYVerts(prot_idx);
      values[row++][tp].first = parseFloat(byRes.first, 1);
      values[row++][tp].first = parseFloat(byRes.second, 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercCallsAcc(prot_idx, 2), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercCallsAcc(prot_idx, 1), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercCallsAcc(prot_idx, 0), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercCallsLen(prot_idx, 0), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercCallsLen(prot_idx, 1), 1);
      values[row++][tp].first
          = parseFloat(mapped_sps_proj->getPercCallsLen(prot_idx, 2), 1);
    }
  }
}
