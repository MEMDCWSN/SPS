////////////////////////////////////////////////////////////////////////////////
#include "ReportTableProtein.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportTableProtein::ReportTableProtein(const string &projectPath, const string &tableFilename, int columnFilter)
: ReportTableBase(projectPath, tableFilename)
{
  // set the filter column
  setFilterColumn(columnFilter);
  // define view (cluster list)
  defineView();
  // set sort column
  setIdColumn(0);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableProtein::defineView(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;

  // colTypes[0] -> (CTstring) Protein name in fasta file
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Protein";
  auxS->text         = "<1>";
  auxS->link         = "protein.<0>.0.html";
  m_colTypes.push_back(auxS);

  // colTypes[1] -> (CTstring) Protein description
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Description";
  auxS->text         = "<2>";
  m_colTypes.push_back(auxS);

  // colTypes[2] -> (CTstring) contigs
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Contigs";
  auxS->text         = "<3>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) spectra
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Spectra";
  auxS->text         = "<4>";
  m_colTypes.push_back(auxS);

  // colTypes[4] -> (CTstring) Amino acids
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Amino acids";
  auxS->text         = "<5>";
  m_colTypes.push_back(auxS);

  // colTypes[5] -> (CTstring) coverage %
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Coverage (%)";
  auxS->text         = "<6>";
  m_colTypes.push_back(auxS);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableProtein::defineView2(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;

  // colTypes[0] -> (CTstring) Protein name in fasta file
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Protein";
  auxS->text         = "<1>";
  auxS->link         = "protein.<0>.0.html";
  m_colTypes.push_back(auxS);
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
