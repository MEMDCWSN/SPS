////////////////////////////////////////////////////////////////////////////////
#include "ReportTableHeader.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportTableHeader::ReportTableHeader(const string &projectPath, const string &tableFilename)
: ReportTableBase(projectPath, tableFilename)
{
  m_direction = 1;
  // define view (cluster list)
  defineView();
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableHeader::defineView(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString          *auxS;
  ReportColumnTypeStringMultiple  *auxM;
  ReportColumnTypeImageOnDemand   *auxI;
  ReportColumnTypeBox             *auxB;

  // colTypes[0] -> (CTstring) Job name
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Job";
  auxS->text         = "<0>";
  m_colTypes.push_back(auxS);

  // colTypes[1] -> (CTstring) User
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "User";
  auxS->text         = "<1>";
  m_colTypes.push_back(auxS);

  // colTypes[2] -> (CTstring) Status
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Status";
  auxS->text         = "<2>";
  m_colTypes.push_back(auxS);

  // colTypes[3] -> (CTstring) Elapsed
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Elapsed";
  auxS->text         = "<SPSTIME>";
  m_colTypes.push_back(auxS);

  // colTypes[4] -> (CTstring) Log
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Log";
  auxS->link         = "reportlog.txt";
  auxS->text         = "<4>";
  m_colTypes.push_back(auxS);



  // colTypes[5] -> (CTstring) Data
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel  = "Data";

  // colTypes[5] -> (CTstring) Data
  auxS = new ReportColumnTypeString();
  auxS->text         = "Group by Contig";
  auxS->link         = "<5>.0.html";
  auxB->sequences.push_back(auxS);

  // colTypes[5] -> (CTstring) Data
  auxS = new ReportColumnTypeString();
  auxS->text         = "Group by Protein";
  auxS->link         = "<6>.html";
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);



  // colTypes[6] -> (CTstring) Cluster Data
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel  = "Cluster Data";

  // colTypes[5] -> (CTstring) Data
//  auxS = new ReportColumnTypeString();
//  auxS->text         = "All Clusters (txt)";
//  auxS->link         = "<7>";
//  auxB->sequences.push_back(auxS);

  // colTypes[5] -> (CTstring) Data
  auxM = new ReportColumnTypeStringMultiple();
  auxM->linkPrefix   = "spectra.";
  auxM->linkSuffix   = ".0.html";
  auxM->text         = "<8>";
  auxM->link         = "<9>";
  auxB->sequences.push_back(auxM);

  m_colTypes.push_back(auxB);
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
