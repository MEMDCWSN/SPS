/*
 * MappedContigSetTable.cpp
 *
 *  Created on: Dec 7, 2012
 *      Author: aguthals
 */

#include "MappedContigSetTable.h"

namespace specnets
{

  MappedContigSetTable::MappedContigSetTable() :
    OutputTable(), mapped_sps_proj(0x0)
  {

  }

  MappedContigSetTable::MappedContigSetTable(MappedSpecnets* _mapped_sps_proj) :
    OutputTable(), mapped_sps_proj(_mapped_sps_proj)
  {
  }

  /**
   * Prepares output table with all necessary statistics
   * @return
   */
  void MappedContigSetTable::prepareTable()
  {
    const unsigned int numRows = mapped_sps_proj->contigs->size();
    const unsigned int numCols = 18;

    pair<string, bool> dflt("", false);
    pair<string, bool> dflt_header("", true);
    vector<pair<string, bool> > dflt_vec(numCols, dflt);

    values.resize(numRows + 1, dflt_vec);

    int row = 0, col = 0;
    pair<float, float> res;

    values[row].assign(numCols, dflt_header);
    values[row][col++].first = "Contig Index";
    values[row][col++].first = "DBSearch Protein";
    values[row][col++].first = "Alignment Protein";
    values[row][col++].first = "Sequence Length (AA)";
    values[row][col++].first = "Contig Parent Mass (Da)";
    values[row][col++].first = "Sequencing Accuracy (%)";
    values[row][col++].first = "Un-annotated Calls (%)";
    values[row][col++].first = "# Single Jump Calls";
    values[row][col++].first = "# Consec Gap Calls";
    values[row][col++].first = "# Non-consec Gap Calls";
    values[row][col++].first = "Reversed";
    values[row][col++].first = "# Vertices";
    values[row][col++].first = "# Spectra";
    values[row][col++].first = "# Peptides";
    values[row][col++].first = "# Annotated Spectra";
    values[row][col++].first = "# Mapped Spectra";
    values[row][col++].first = "Begin Align Idx";
    values[row][col++].first = "End Align Idx";

    for (unsigned int i = 0; i < mapped_sps_proj->contigs->size(); i++)
    {
      row++;
      col = 0;
      MappedContig* contig = &(*(mapped_sps_proj->contigs))[i];

      setValue(row, col++, parseInt(i));
      setValue(row, col++, parseInt(contig->starProtIdx));
      setValue(row, col++, parseInt(contig->vertProtIdx));
      setValue(row, col++, parseInt(contig->lastResidue - contig->firstResidue
          + 1));
      setValue(row, col++, parseFloat(contig->parentMass, 1));
      res = contig->getPercCallsAcc(2);
      float perc = (res.second > 0.1) ? res.first * 100.0 / res.second : 0;
      setValue(row, col++, parseFloat(perc, 1));
      res = contig->getPercCallsAcc(0);
      perc = (res.second > 0.1) ? res.first * 100.0 / res.second : 0;
      setValue(row, col++, parseFloat(perc, 1));
      res = contig->getPercCallsLen(0);
      setValue(row, col++, parseFloat(res.first, 0));
      res = contig->getPercCallsLen(1);
      setValue(row, col++, parseFloat(res.first, 0));
      res = contig->getPercCallsLen(2);
      setValue(row, col++, parseFloat(res.first, 0));
      string val = (contig->reversed) ? "1" : "0";
      setValue(row, col++, val);
      setValue(row, col++, parseInt(contig->length));
      setValue(row, col++, parseInt(contig->numSpecs));
      setValue(row, col++, parseInt(contig->numPeptides));
      setValue(row, col++, parseInt(contig->getNumAnnotSpec()));
      setValue(row, col++, parseInt(contig->getNumMappedSpec()));
      setValue(row, col++, parseInt(contig->firstResidue));
      setValue(row, col++, parseInt(contig->lastResidue));
    }

  }

}
