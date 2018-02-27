///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_PROTEIN_COVERAGE_H__
#define __REPORT_TABLE_PROTEIN_COVERAGE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>

#include "ReportTableBase.h"
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
#define TABLE_PROTEIN_COVERAGE_FILTER_COL_PROTEIN  0

////////////////////////////////////////////////////////////////////////////////
// Generated table, per column:
////////////////////////////////////////
//
// cells[0][0] -> text   --> protein index (1-based)
// cells[1][0] -> text   --> protein name
// cells[2][0] -> text   --> Protein length
// cells[3][0] -> text   --> Protein sequence
// cells[4][0] -> text   --> csps contig data
// cells[5][0] -> text   --> sps contig data
///////////////////////////////////////////////////////////////////////////////
 /*! \brief Report table class for protein coverage

   Defines a protein coverage table for reports

   */

class ReportTableProteinCoverage : public ReportTableBase {

 protected:


 public:

  // Constructors and destructor

  //! \name CONSTRUCTORS
  //@{

  /*! \brief Contructor. Takes a project path and the table file name as parameters
   */
  ReportTableProteinCoverage(const string &projectPath, const string &tableFilename, int columnFilter = TABLE_PROTEIN_COVERAGE_FILTER_COL_PROTEIN);

  /*! \brief Default destructor
   */
  virtual ~ReportTableProteinCoverage() {};

};
///////////////////////////////////////////////////////////////////////////////
};
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
