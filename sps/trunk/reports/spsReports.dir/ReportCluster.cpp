///////////////////////////////////////////////////////////////////////////////
#include "ReportCluster.h"
///////////////////////////////////////////////////////////////////////////////
#include "ReportTableClusterConsensus.h"
#include "ReportTableInputSpectra.h"
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
///////////////////////////////////////////////////////////////////////////////
ReportCluster::ReportCluster(const string &projectPath, const string &clusterTableFilename, const string &inputSpectraTableFilename)
{
  // The cluster consensus spectra on top of page. Set the filter column to cluster index column
  ReportTableBase *c = new ReportTableClusterConsensus(projectPath, clusterTableFilename, TABLE_CLUSTER_FILTER_COL_CLUSTER);
  // Define the specific view (spectra on top)
  c->defineView2();
  // no borders and no header
  c->noBorders(); c->noHeaders();
  // Add the table to the report
  addTable(c);

  // The input spectra list. By default, it's filtered at cluster consensus column
  ReportTableBase *c2 = new ReportTableInputSpectra(projectPath, inputSpectraTableFilename);
  // Define the specific view (spectra on top)
  c2->defineView2();
  // Add the table to the report
  addTable(c2);
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
